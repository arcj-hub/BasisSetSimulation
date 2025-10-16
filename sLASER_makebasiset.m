% run_simSemiLASERShaped_fast.m
% Jamie Near, McGill University 2015.
% Dana Goerzen and Jamie Near, McGill University 2021.
% Fast version by Muhammad G Saleh (Johns Hopkins University School of Medicine, 2019)

clear all;
close all;
%
curfolder=fileparts(mfilename('fullpath'));

%--------------------------------------------------------------------------
% DESCRIPTION & MODIFICATIONS:
%--------------------------------------------------------------------------

% This script was modified by Niklaus Zoelch & Jessica Archibald 2020-2025 to include:

%     a.	A 0 reference peak
%     b.	A loop to run through all the metabolites
%     c.	saving .raw, .png, and .mat files
%     d.	output a .pdf and .basis file -> using modified functions from Osprey- Georg Oeltzschner, Johns Hopkins University 2019.

% USAGE:

% This script simulates a semi-LASER  experiment with fully shaped refocusing
% pulses. Coherence order filtering is employed to only simulate desired signals
% This results in a 4x speed up compared to phase cycling
% (see deprecated run_simSemiLASERShaped_fast_phCyc.m)
% Furthermore, simulations are run at various locations in space to account for the
% within-voxel spatial variation of the metabolite signal.  Summation
% across spatial positions is performed. The MATLAB parallel computing toolbox
% (parfor loop) was used to accelerate the simulations.  Acceleration
% is currently performed in the direction of the slice selective pulse along
% the x-direction, but this can be changed.  Up to a factor of 12 acceleration
% can be achieved using this approach. To achieve
% faster perfomance compared to the original 'run_simSemiLASER_shaped.m' function,
% this code uses the method described by Yan Zhang et al. Med Phys 2017;44(8):
% 4169-78.  Some additional acceleration is currently performed using parfor
% loops in both x and y directions.  To enable the use of the MATLAB
% parallel computing toolbox, initialize the multiple worked nodes using
% "matlabpool size X" where "X" is the number of available processing
% nodes.  If the parallel processing toolbox is not available, then replace
% the "parfor" loop with a "for" loop.
%
%
% INPUTS:

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
% out               = Simulation results, summed over all space.

ToolboxCheck

% ************INPUT PARAMETERS**********************************
%
% Define the variable Basis_name at the beginning of your script
basis_name = 'lcm_gamma_new.basis' ; %keep "_gamma_" 

% Path to the FID-A Folder
%pathtofida='W:/SharedProgramme/FID-A-git/FID-A';
pathtofida='/Users/jessicaarchibald/BasisREMY/FID-A-Fork';

% Where everything is saved
folder_to_save=[curfolder,'/my_basis/'];

%--------------------------------------------------------------------------
% Metabolite Information
%--------------------------------------------------------------------------
path_to_spinsystem=[curfolder,'/my_mets/my_spinSystem'];
%--------------------------------------------------------------------------

save_result=true;
%
Waveform=[curfolder,'/my_pulse/standardized_goia.txt'];
B1max= 22;% in mikro Tesla // set to [] to be prompted with figure
flip_angle=180;
refTp=4.5008; %duration of refocusing pulses[ms]
Npts=2048; %number of spectral points
sw=2000; %spectral width [Hz]
lw=2; %linewidth of the output spectrum [Hz]
Bfield=3; %Magnetic field strength in [T]
thkX=2.; %slice thickness of x refocusing pulse [cm]
thkY=2.; %slice thickness of y refocusing pulse [cm]
%
fovX=3; %size of the full simulation Field of View in the x-direction [cm]
fovY=3; %size of the full simulation Field of View in the y-direction [cm]
%
nX=64; %Number of grid points to simulate in the x-direction
nY=64; %Number of grid points to simulate in the y-direction
%
% full voxel
x=linspace(-fovX/2,fovX/2,nX); %X positions to simulate [cm]
y=linspace(-fovY/2,fovY/2,nY); %Y positions to simulate [cm]
%
te=30;%timing of the pulse sequence [ms]
centreFreq=2.67; %Centre frequency of MR spectrum [ppm]
%
fovX=-x(1)+x(end);
fovY=-y(1)+y(end);
% select the metabolites you want to simulate
% the basis set is created with all metabolites in the matfiles_post

spinSysList={'PE','Asc','Scyllo','Glu','Cr','NAA','NAAG','PCr','GSH','Gly','Glc','GPC',...
    'PCh','Ala','Asp','Gln', 'Ins', 'Lac'};

% shift
shift_in_ppm=(4.65-centreFreq);

% ************ END OF INPUT PARAMETERS BY USER **********************************
% addpath
addpath(genpath(pathtofida),'-begin');

% Add the scripts adapted from fidA
% Will check here first to find the function
addpath(fullfile(curfolder, 'dependencies'), '-begin');

%--------------------------------------------------------------------------
%Load RF waveform
%--------------------------------------------------------------------------
%Niklaus : set inv / macht keinen unterschied oder?
rfPulse=io_loadRFwaveform(Waveform,'inv',0,B1max);
w1max=rfPulse.w1max; % to save the input value
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
sysRef.J=0;
sysRef.shifts=0;
sysRef.scaleFactor=1;
sysRef.name='Ref_0ppm';
%
sysRef.centreFreq=centreFreq;
%
[ ref] = run_mysLASERShaped_fast(rfPulse,refTp,Npts,sw,lw,Bfield,thkX,thkY,x,y,te,sysRef,flip_angle);
%
tau1=15; tau2=13;%fake timing
refjustforppmrange=sim_press(Npts,sw,Bfield,lw,sysRef,tau1,tau2);
%-------------------------------------------------------------------------
% shift for the ref
%-------------------------------------------------------------------------
% https://ch.mathworks.com/matlabcentral/newsreader/view_thread/243061

freqShift_hz=shift_in_ppm*(Bfield*42.577478); % in Hz
%--------------------------------------------------------------------------------------
ref.fids=ref.fids.*exp(-(1i*2*pi*freqShift_hz).*ref.t).';

%-------------------------------------------------------------------------
%Load spin systems
%-------------------------------------------------------------------------
load(path_to_spinsystem)
%-------------------------------------------------------------------------

for met_nr=1:size(spinSysList,2)
    %
    spinSys=spinSysList{met_nr}; %spin system to simulate
    sys=eval(['sys' spinSys]);
    % Schreibe die einfach im ersten rein
    sys(1).centreFreq=centreFreq;

    %-------------------------------------------------------------------------
    % Simulation
    %-------------------------------------------------------------------------
    [ out] = run_mysLASERShaped_fast(rfPulse,refTp,Npts,sw,lw,Bfield,thkX,thkY,x,y,te,sys,flip_angle);

    % Save before the shift -
    save_out_mat=[folder_to_save,'matfiles_pre'];
    if (exist(save_out_mat,'dir')==0)
                 mkdir(save_out_mat);
    end
    save([save_out_mat,'/',spinSys],'out')

        %-------------------------------------------------------------------------
    % Add shift here and later for every simulated
    %-------------------------------------------------------------------------
    out.fids=out.fids.*exp(-(1i*2*pi*freqShift_hz).*ref.t).';

    %-------------------------------------------------------------------------
    % add tms ref
    %-------------------------------------------------------------------------
    out=op_addScans(out,ref);


    save_figure=[folder_to_save,'figures'];
    if (exist(save_figure,'dir')==0)
                 mkdir(save_figure);
    end
    % figure
    figure;plot(refjustforppmrange.ppm,real(ifftshift(ifft(out.fids))),'b')
    set(gca,'xdir','reverse')
    colormap;set(gcf,'color','w');
    xlim([-1 5])
    xlabel('ppm');
    title(['met with ref ',spinSys])
    print('-dpng','-r300',[save_figure,'\',spinSys])
    %---------------------------------------------------------------------
    % add inforation that then can be stored in BASIS
    %---------------------------------------------------------------------
    out.name=io_sysname({sys.name});
    out.centerFreq=centreFreq; % This is needed for the check within fit_LCMmakeBasis
    out.w1max=w1max;

    save_raw=[folder_to_save,'raw'];
    if save_result
        if (exist(save_raw,'dir')==0)
                 mkdir(save_raw);
        end

        RF=io_writelcmraw(out,[save_raw,'/',spinSys,'.RAW'],spinSys);

    end

    % Saving after shift
    save_out_mat_end=[folder_to_save,'matfiles_post'];
    if (exist(save_out_mat_end,'dir')==0)
                 mkdir(save_out_mat_end);
    end
    save([save_out_mat_end,'/',spinSys],'out')


end

disp('Running fit_makeLCMBasis...');

%
BASIS=fit_makeLCMBasis(save_out_mat_end, false, [folder_to_save,'/', basis_name],'Philips','sLASER');
%


% %Vizualize your created basis set
% figure;plot(BASIS.ppm,real(BASIS.specs));legend(BASIS.name)
% set(gca,'xdir','reverse','XGrid','on')

disp('Done');
