clear;
close all;

%-------------------------------------------------------------------------
% Edit here
%--------------------------------------------------------------------------

%-------------------------------------------------------------------------
%path to the base
orig_spinSystem='W:/SharedProgramme/FID-A-git/FID-A/simulationTools/metabolites/spinSystems';
% name and path for the new spinSystem
% spinSystem ist the collection of all 
new_spinSystem_name='my_spinSystem';
% stores output in the same folder as the script runs
scriptPath = fileparts(mfilename('fullpath'));
goal_path=[scriptPath,'/my_mets'];
if (exist(goal_path,'dir')==0)
            mkdir(goal_path);
end


%load the original spinsystems
load(orig_spinSystem)
%--------------------------------------------------------------------------
%sys GABA
%--------------------------------------------------------------------------
% source of values : https://onlinelibrary.wiley.com/doi/epdf/10.1002/nbm.3336
% order [2,2',3,3',4,4']
sysGABA_govind.name='GABA';
sysGABA_govind.scaleFactor=1;
sysGABA_govind.shifts=[2.2840;2.2840;1.889;1.889;3.0128;3.0128];
%        2,   2',      3,       3',    4,      4' 
% 2     [0, -10.744,   7.775,   6.173, 0,      0; ...       
% 2'     0,   0,       7.432,   7.933, 0,      0; ...       
% 3      0,   0,       0,     -13.121, 5.372, 10.578; ...  
% 3'     0,   0,       0,       0,     7.127,  6.982; ...  
% 4      0,   0,       0,       0,     0,    -12.021; ...  
% 4'     0,   0,       0,       0,     0,      0];          

sysGABA_govind.J=[0, -10.744,   7.775, 6.173, 0,      0; ...       
                0,   0,       7.432, 7.933, 0,      0; ...       
                0,   0,       0,   -13.121, 5.372, 10.578; ...   
                0,   0,       0,     0,     7.127,  6.982; ...  
                0,   0,       0,     0,     0,    -12.021; ...   
                0,   0,       0,     0,     0,      0];         


% save the sys if you want to
save([goal_path,'/GABA_govind'],'sysGABA_govind') 
%--------------------------------------------------------------------------
% end added Met
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%sys Tau
%--------------------------------------------------------------------------
% nach neustem Govin
% source of values : https://onlinelibrary.wiley.com/doi/epdf/10.1002/nbm.3336
% order [1,1',2,2']
sysTau_govind.name='Tau';
sysTau_govind.scaleFactor=1;
sysTau_govind.shifts=[3.4206,3.4206,3.2459,3.2459];
%        1,   1',      2,       2'
% 1     [0, -12.438,   6.742,   6.464;  ...       
% 1'     0,   0,       6.403,   6.792;  ...       
% 2      0,   0,       0,     -12.930;  ...  
% 2'     0,   0,       0,       0];      ...  

sysTau_govind.J=[0, -12.438,   6.742,   6.464;  ...       
              0,   0,        6.403,   6.792;  ...       
              0,   0,        0,     -12.930;  ...  
              0,   0,        0,       0];      ...        

% save the sys if you want to
save([goal_path,'/Tau_govind'],'sysTau_govind') 
%--------------------------------------------------------------------------
% end added Met
%--------------------------------------------------------------------------



