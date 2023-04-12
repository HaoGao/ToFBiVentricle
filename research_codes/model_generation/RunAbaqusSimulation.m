%%% this is trying to run a simulation and post-processing the results 
clear all; close all; clc;

% LVWM_config;
LVWM_config;

%%now need to choose phase to segment
list_phase = {'early_diastole', 'end_diastole', 'end_systole'};
[idx_phase, tf] = listdlg('ListString', list_phase);
cd(resultDir);
if ~exist(list_phase{idx_phase},'dir')
    mkdir(list_phase{idx_phase});
    cd(list_phase{idx_phase});
    phase_resultDir = pwd();    
else
    cd(list_phase{idx_phase});
    phase_resultDir = pwd();
end
cd(workingDir);
phase_selected = list_phase{idx_phase};

cd(phase_resultDir);
cd(solidworksDir);
solidworksDir = pwd();
cd(gmeshDir);
gmeshDir = pwd();
cd(workingDir);


cd(phase_resultDir);
cd(abaqusSimulationDir);
abaqusSimulationDir = pwd();
cd(workingDir);

opt_file_config = 'optimization_config';
cd(resultDir);
run(opt_file_config);
cd(workingDir);



optimize_opt.abaqusSimulationDir = abaqusSimulationDir;
optimize_opt.meshDir = gmeshDir;
optimize_opt.segDir = phase_resultDir;
optimize_opt.workingDir = workingDir;
optimize_opt.pythonOriginalFilesDir = pythonOriginalFilesDir;
optimize_opt.pythonfilename = pythonfilename;
optimize_opt.lineNoForODBName = lineNoForODBName;
optimize_opt.lineNoForDisName = lineNoForDisName;
optimize_opt.abaqus_inputfile = abaqus_inputfile;
optimize_opt.abaqus_dis_out_filename = abaqus_dis_out_filename;

cd(optimize_opt.meshDir);
load abaqusInput;
cd(workingDir);
optimize_opt.abaqusInput = abaqusInput;

DislacementExtractByPython(optimize_opt.abaqusSimulationDir,...
                           optimize_opt.abaqus_inputfile,...
                           optimize_opt.pythonOriginalFilesDir,...
                           optimize_opt.pythonfilename, ...
                           optimize_opt.lineNoForODBName, ...
                           optimize_opt.lineNoForDisName,...
                           optimize_opt.abaqus_dis_out_filename)

ResVol = biVenCavityCalculation(optimize_opt);






