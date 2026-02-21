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

ToolboxCheck;

% ************ INPUT PARAMETERS **********************************

% Define the variable Basis_name at the beginning of your script
basis_name='lcm_gamma_new.basis'; %keep "_gamma_"
maindir=fileparts(mfilename("fullpath"));
addpath(genpath(maindir));
folder_to_save='~/Desktop/makebasisset_output';
save_result=true;
vendor='Philips';
sequence='sLASER';
refocWaveform='standardized_goia.txt'; %name of refocusing pulse waveform
flip_angle=180;
refTp=4.4496; %duration of refocusing pulses[ms]I've I've
Npts=4096; %number of spectral points
sw=4000; %spectral width [Hz]
lw=2; %linewidth of the output spectrum [Hz]
Bfield=3; %Magnetic field strength in [T]
thkX=2.4; %slice thickness of x refocusing pulse [cm]
thkY=2.2; %slice thickness of y refocusing pulse [cm]
fovX=3; %size of the full simulation Field of View in the x-direction [cm]
fovY=3; %size of the full simulation Field of View in the y-direction [cm]
nX=40; %Number of grid points to simulate in the x-direction
nY=40; %Number of grid points to simulate in the y-direction
x=linspace(-fovX/2,fovX/2,nX); %X positions to simulate [cm]
y=linspace(-fovY/2,fovY/2,nY);
te=32;%timing of the pulse sequence [ms]
centreFreq=2.02; %Centre frequency of MR spectrum [ppm]

fovX=-x(1)+x(end);
fovY=-y(1)+y(end);

% spin systems
spinSysList={'PE', 'Asc', 'Scyllo','Glu','Cr','NAA','NAAG','PCr','GSH','Gly','Glc','GPC',...
    'PCh','Ala','Asp','GABA', 'Gln', 'Ins', 'Lac', 'Tau'};

% shift
shift_in_ppm=(4.65-centreFreq);

% ************ END OF INPUT PARAMETERS BY USER **********************************

%--------------------------------------------------------------------------
%Load RF waveform
%--------------------------------------------------------------------------
rfPulse=io_loadRFwaveform(refocWaveform,'ref',0);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
sysRef.J=0;
sysRef.shifts=0;
sysRef.scaleFactor=1;
sysRef.name='Ref_0ppm';
sysRef.centreFreq=centreFreq;
ref = run_mysLASERShaped_fast(rfPulse,refTp,Npts,sw,lw,Bfield,thkX,thkY,x,y,te,sysRef,flip_angle);

tau1=15; %fake timing
tau2=13; %fake timing
refjustforppmrange=sim_press(Npts,sw,Bfield,lw,sysRef,tau1,tau2);
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% shift
ppm_range=ref.ppm(1)-ref.ppm(end);
ppm_per_point=ppm_range/size(ref.ppm,2);
shift_in_points=round(shift_in_ppm/ppm_per_point);
ref.fids=ref.fids.*exp(-1i*2*pi*shift_in_points*(0:1:(size(ref.fids,1)-1)).'/(size(ref.fids,1)));
%--------------------------------------------------------------------------
% Additional Metabolites
[sysETH,sysAcetate,sysAcac,sysSucc,sysGlyc,sysVal,sysAceton,sysbHBHM]=define_spin_systems;

%-------------------------------------------------------------------------
%Load spin systems (for the rest)
load(fullfile(maindir, 'my_mets', 'my_spinSystem.mat'));
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
    save_out_mat=fullfile(folder_to_save,'matfiles_pre');
    if ~exist(save_out_mat,'dir')
        mkdir(save_out_mat);
    end
    save([save_out_mat,filesep,spinSys],'out');

    %-------------------------------------------------------------------------
    % Add shift here and later for every simulated
    %
    %-------------------------------------------------------------------------
    out.fids=out.fids.*exp(-1i*2*pi*shift_in_points*(0:1:(size(out.fids,1)-1)).'/(size(out.fids,1)));

    %-------------------------------------------------------------------------
    % add TMS ref
    %-------------------------------------------------------------------------
    out=op_addScans(out,ref);

    save_figure=fullfile(folder_to_save,'figures');
    if ~exist(save_figure,'dir')
        mkdir(save_figure);
    end

    % figure
    figure;
    plot(refjustforppmrange.ppm,real(ifftshift(ifft(out.fids))),'b');
    set(gca,'xdir','reverse');
    colormap;set(gcf,'color','w');
    xlim([-1 5]);
    xlabel('ppm');
    title(['figure with ref',spinSys]);
    print('-dpng','-r300',[save_figure,filesep,spinSys]);

    out.name=spinSys;
    out.centreFreq=centreFreq; % This is needed for the check within fit_LCMmakeBasis.
    save_raw=fullfile(folder_to_save,'raw');

    if save_result
        if ~exist(save_raw,'dir')
            mkdir(save_raw);
        end
        RF=io_writelcmraw(out,[save_raw, filesep, spinSys '.raw'],spinSys);
    end

    % Saving after shift
    save_out_mat_end=fullfile(folder_to_save,'matfiles_post');
    if ~exist(save_out_mat_end,'dir')
        mkdir(save_out_mat_end);
    end
    save([save_out_mat_end,filesep,spinSys],'out');

end

disp('Running fit_makeLCMBasis...');

BASIS=fit_makeLCMBasis(save_out_mat_end, false, [folder_to_save, filesep, basis_name], vendor, sequence);

rmpath(genpath(maindir));
disp('Done');
