% run_sLASER_makebasisset
%
% Contributors:
%
% Jamie Near, McGill University, 2015
% Georg Oeltzschner, Johns Hopkins University School of Medicine, 2019
% Muhammad G Saleh, Johns Hopkins University School of Medicine, 2019
% Dana Goerzen and Jamie Near, McGill University, 2021
% Niklaus Zölch, Universität Zürich, 2024
% Jessica Archibald, Weill Cornell Medicine, 2024
% Mark Mikkelsen, Weill Cornell Medicine, 2026
%
% DESCRIPTION & MODIFICATIONS:
%
% This script was modified by Niklaus Zölch & Jessica Archibald to include:
%
%   a. A 0 ppm reference peak, useful depending on the fitting software
%   b. A loop to run through the selected metabolites
%   c. Saving .raw, .png, and .mat files
%   d. Output a .pdf and .basis file using modified functions from Osprey
%
% USAGE:
%
% This script simulates a semi-LASER experiment with fully shaped
% refocusing pulses. Coherence order filtering is employed to only simulate
% desired signals. This results in a 4x speed-up compared to phase cycling
% (see deprecated run_simSemiLASERShaped_fast_phCyc.m). Furthermore,
% simulations are run at various locations in space to account for the
% within-voxel spatial variation of the metabolite signal. Summation across
% spatial positions is performed. The MATLAB parallel computing toolbox
% (parfor loop) was used to accelerate the simulations. Acceleration is
% currently performed in the direction of the slice selective pulse along
% the x-direction, but this can be changed. Up to a factor of 12
% acceleration can be achieved using this approach. To achieve faster
% perfomance compared to the original 'run_simSemiLASER_shaped.m' function,
% this code uses the method described by ﻿Zhang et al. (2017)
% doi:10.1002/mp.12375. Some additional acceleration is currently performed
% using parfor loops in both x and y directions. To enable the use of the
% MATLAB parallel computing toolbox, initialize the multiple worker nodes
% using "matlabpool size X" where "X" is the number of available processing
% nodes.  If the parallel processing toolbox is not available, then replace
% the "parfor" loop with a "for" loop.

% INPUTS:
%
% To run this script, there is technically only one input argument:
% spinSys           = spin system to simulate
% However, the user should also edit the following parameters as
% desired before running the function:
% refocWaveform     = name of refocusing pulse waveform.
% refTp             = duration of refocusing pulses[ms]
% Bfield            = Magnetic field strength in [T]
% Npts              = number of spectral points
% sw                = spectral width [Hz]
% Bfield            = magnetic field strength [Tesla]
% lw                = linewidth of the output spectrum [Hz]
% thkX              = slice thickness of x refocusing pulse [cm]
% thkY              = slice thickness of y refocusing pulse [cm]
% fovX              = full simulation FOV in the x direction [cm]
% fovY              = full simulation FOV in the y direction [cm]
% nX                = number of spatial grid points to simulate in x-direction
% nY                = number of spatial grid points to simulate in y-direction
% taus              = vector of pulse sequence timings  [ms]
%
% OUTPUTS:
%
% out               = Simulation results, summed over all space.

clear;
clc;
close all;



% ************ INPUT PARAMETERS **********************************

% Define the variable Basis_name at the beginning of your script
basis_name='lcm_gamma_new.basis'; %keep "_gamma_"
main_dir=fileparts(mfilename("fullpath"));
addpath(genpath(main_dir));
ToolboxCheck;
output_folder=fullfile(main_dir,'my_basis'); % or select a folder somewhere else e.g. '~/Desktop/makebasisset_output'
save_result=true;
complete_run=true; % 
show_plots=false;
vendor='Philips';
sequence='sLASER';
refocWaveform='standardized_goia.txt'; %name of refocusing pulse waveform
flip_angle=180;
refTp=4.4496; %duration of refocusing pulses[ms]
Npts=4096; %number of spectral points
sw=4000; %spectral width [Hz]
lw=2; %linewidth of the output spectrum [Hz]
Bfield=3; %Magnetic field strength in [T]
thkX=2.4; %slice thickness of x refocusing pulse [cm]
thkY=2.2; %slice thickness of y refocusing pulse [cm]
fovX=3; %size of the full simulation Field of View in the x-direction [cm]
fovY=3; %size of the full simulation Field of View in the y-direction [cm]
nX=64; %Number of grid points to simulate in the x-direction
nY=64; %Number of grid points to simulate in the y-direction
x=linspace(-fovX/2,fovX/2,nX); %X positions to simulate [cm]
y=linspace(-fovY/2,fovY/2,nY);
te=32;%timing of the pulse sequence [ms]
centreFreq=2.02; %Centre frequency of MR spectrum [ppm]
B1max=[22]; %B1max for refocusing pulses; if empty, B1max is calculated automatically

fovX=-x(1)+x(end);
fovY=-y(1)+y(end);

% spin systems
spinSysList={'PE', 'Asc', 'Scyllo','Glu','Cr','NAA','NAAG','PCr','GSH','Gly','Glc','GPC',...
    'PCh','Ala','Asp','GABA', 'Gln', 'Ins', 'Lac', 'Tau'};

% shift
shift_in_ppm=(4.65-centreFreq);

% ************ END OF INPUT PARAMETERS BY USER **********************************

%%JA edit: confirmation popup with current-run and final basis-set contents
current_run_summary=strjoin(spinSysList,', ');
if numel(current_run_summary) > 260
    current_run_summary=[current_run_summary(1:260) ' ...'];
end

existing_basis_folder=fullfile(output_folder,'matfiles_post');
existing_basis_mets={};
if exist(existing_basis_folder,'dir') && ~complete_run
    existing_basis_files=dir(fullfile(existing_basis_folder,'*.mat'));
    existing_basis_mets=cell(size(existing_basis_files));
    for existing_idx=1:numel(existing_basis_files)
        [~,existing_basis_mets{existing_idx},~]=fileparts(existing_basis_files(existing_idx).name);
    end
end

if complete_run
    final_basis_mets=spinSysList(:);
else
    final_basis_mets=unique([existing_basis_mets(:); spinSysList(:)],'stable');
end

final_basis_summary=strjoin(final_basis_mets',', ');
if numel(final_basis_summary) > 260
    final_basis_summary=[final_basis_summary(1:260) ' ...'];
end

confirmation_lines={
    'Are you sure you want to simulate a basis set with:'
    ' '
    ['Basis file: ' basis_name]
    ['Output folder: ' output_folder]
    ['Vendor / Sequence: ' vendor ' / ' sequence]
    ['Refocusing waveform: ' refocWaveform]
    ['Flip angle: ' num2str(flip_angle) ' deg']
    ['Refocusing duration: ' num2str(refTp) ' ms']
    ['B-field: ' num2str(Bfield) ' T']
    ['TE: ' num2str(te) ' ms']
    ['Npts / SW / LW: ' num2str(Npts) ' / ' num2str(sw) ' Hz / ' num2str(lw) ' Hz']
    ['Slice thickness X/Y: ' num2str(thkX) ' / ' num2str(thkY) ' cm']
    ['FOV X/Y: ' num2str(fovX) ' / ' num2str(fovY) ' cm']
    ['Grid points X/Y: ' num2str(nX) ' / ' num2str(nY)]
    ['Centre frequency: ' num2str(centreFreq) ' ppm']
    ['B1max: ' mat2str(B1max)]
    ' '
    ['You are now simulating: ' current_run_summary]
    ['Your final basis set will contain: ' final_basis_summary]
    };

confirmation_text=sprintf('%s\n',confirmation_lines{:});

popup_handle=dialog( ...
    'Name','Confirm Basis Set Simulation', ...
    'Position',[200 120 760 520], ...
    'Color',[0.97 0.97 0.99], ...
    'WindowStyle','modal');
setappdata(popup_handle,'popup_choice','No');

uicontrol( ...
    'Parent',popup_handle, ...
    'Style','text', ...
    'String','Basis Set Simulation Check', ...
    'Position',[25 475 320 26], ...
    'HorizontalAlignment','left', ...
    'FontSize',16, ...
    'FontWeight','bold', ...
    'BackgroundColor',[0.97 0.97 0.99], ...
    'ForegroundColor',[0.12 0.18 0.32]);

uicontrol( ...
    'Parent',popup_handle, ...
    'Style','edit', ...
    'Max',2, ...
    'Min',0, ...
    'Enable','inactive', ...
    'String',confirmation_text, ...
    'Position',[25 85 710 380], ...
    'HorizontalAlignment','left', ...
    'FontSize',12, ...
    'BackgroundColor',[1 1 1]);

uicontrol( ...
    'Parent',popup_handle, ...
    'Style','pushbutton', ...
    'String','No', ...
    'Position',[520 22 90 38], ...
    'FontSize',12, ...
    'Callback',@(src,evt) local_set_popup_choice(popup_handle,'No'));

uicontrol( ...
    'Parent',popup_handle, ...
    'Style','pushbutton', ...
    'String','Yes', ...
    'Position',[625 22 90 38], ...
    'FontSize',12, ...
    'FontWeight','bold', ...
    'BackgroundColor',[0.23 0.56 0.34], ...
    'ForegroundColor',[1 1 1], ...
    'Callback',@(src,evt) local_set_popup_choice(popup_handle,'Yes'));

set(popup_handle,'CloseRequestFcn',@(src,evt) local_set_popup_choice(popup_handle,'No'));
uiwait(popup_handle);

popup_choice=getappdata(popup_handle,'popup_choice');
if ishandle(popup_handle)
    delete(popup_handle);
end

if ~strcmp(popup_choice,'Yes')
    fprintf('\nSimulation cancelled by user before launch.\n\n');
    return;
end


if show_plots
    vis_flag='on'; %#ok<*UNRCH>
else
    vis_flag='off';
end


% if all should be rerun then remove the output folder 
if exist(output_folder,'dir') && complete_run
    rmdir(output_folder,'s');
end

% folders for saving
save_out_mat      = fullfile(output_folder,'matfiles_pre');
save_figure       = fullfile(output_folder,'figures');
save_raw          = fullfile(output_folder,'raw');
save_out_mat_end  = fullfile(output_folder,'matfiles_post');

folders = {save_out_mat, save_figure, save_raw, save_out_mat_end};

% create folders if needed
for k = 1:numel(folders)
    if ~exist(folders{k},'dir')
        mkdir(folders{k});
    end
end
%--------------------------------------------------------------------------
%Load RF waveform
%--------------------------------------------------------------------------
rfPulse=io_loadRFwaveform(refocWaveform,'ref',0,B1max);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
sysRef.J=0;
sysRef.shifts=0;
sysRef.scaleFactor=1;
sysRef.name='Ref_0ppm';
sysRef.centreFreq=centreFreq;
ref=run_mysLASERShaped_fast(rfPulse,refTp,Npts,sw,lw,Bfield,thkX,thkY,x,y,te,sysRef,flip_angle);

tau1=15; %fake timing
tau2=13; %fake timing
refjustforppmrange=sim_press(Npts,sw,Bfield,lw,sysRef,tau1,tau2);
%-------------------------------------------------------------------------

%------------------------------------------------
% Shift
%------------------------------------------------
freqShift_hz=shift_in_ppm*(Bfield*42.577478); % in Hz
%-------------------------------------------------------------------------
% Add shift here for the ref
%-------------------------------------------------------------------------
ref.fids=ref.fids.*exp(-(1i*2*pi*freqShift_hz).*ref.t).';
%--------------------------------------------------------------------------
% Additional Metabolites
[sysETH,sysAcetate,sysAcac,sysSucc,sysGlyc,sysVal,sysAceton,sysbHBHM]=define_spin_systems;

%-------------------------------------------------------------------------
%Load spin systems (for the rest)
load(fullfile(main_dir,'my_mets','my_spinSystem.mat'));
%-------------------------------------------------------------------------

for met_nr=1:size(spinSysList,2)

    spinSys=spinSysList{met_nr}; %spin system to simulate
    sys=eval(['sys' spinSys]);
    % Schreibe die einfach im ersten rein
    sys(1).centreFreq=centreFreq;

    %-------------------------------------------------------------------------
    % Simulation
    %-------------------------------------------------------------------------
    out=run_mysLASERShaped_fast(rfPulse,refTp,Npts,sw,lw,Bfield,thkX,thkY,x,y,te,sys,flip_angle);

    %add w1 max
    out.w1max=rfPulse.w1max;

    % Save before the shift -
    save([save_out_mat,filesep,spinSys],'out');

    %-------------------------------------------------------------------------
    % Add shift here for every simulated metabolite
    %-------------------------------------------------------------------------
    out.fids=out.fids.*exp(-(1i*2*pi*freqShift_hz).*ref.t).';

    %-------------------------------------------------------------------------
    % add TMS ref
    %-------------------------------------------------------------------------
    out=op_addScans(out,ref);

    h=figure('Visible',vis_flag);
    clf(h);
    plot(refjustforppmrange.ppm,real(ifftshift(ifft(out.fids))),'b');
    set(gca,'xdir','reverse');
    colormap;set(gcf,'color','w');
    xlim([-1 5]);
    xlabel('ppm');
    title(['figure with ref',spinSys]);
    print(h,'-dpng','-r300',[save_figure,filesep,spinSys]);

    out.name=spinSys;
    out.centreFreq=centreFreq; % This is needed for the check within fit_LCMmakeBasis.


    if save_result
        RF=io_writelcmraw(out,[save_raw, filesep, spinSys '.raw'],spinSys);
    end

    % Saving after shift
    save([save_out_mat_end,filesep,spinSys],'out');

end

fprintf('\nRunning fit_makeLCMBasis...\n\n');
close(101)
BASIS=fit_makeLCMBasis(save_out_mat_end, false, [output_folder, filesep, basis_name], vendor, sequence, vis_flag);

rmpath(genpath(main_dir));
fprintf('\nDone! Output saved in ''%s''\n\n',output_folder);

function local_set_popup_choice(popup_handle,choice_value)
if ishandle(popup_handle)
    setappdata(popup_handle,'popup_choice',choice_value);
    uiresume(popup_handle);
end
end
