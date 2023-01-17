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

%%% abaqus Simulation dir, so the result file will be saved in this folder
cd(phase_resultDir);
cd(abaqusSimulationDir);
abaqusSimulationDir = pwd();
cd(workingDir);

cd(abaqusSimulationDir);

%%% load mesh files
cd(gmeshDir);
load abaqusInput;
fibreDir = load('fibreDir.txt');
sheetDir = load('sheetDir.txt');
cirDir =   load('circumDir.txt');
radDir =   load('radDir.txt');
load LV_RV_assignment;
load LVMeshSegDivisions; %% elRegionsFull nodeElRegionFull saved the strain locations according to each slice
cd(workingDir);

%% pull out settings for optimization 
cd(resultDir);
run('optimization_config');
cd(resultDir);
% if exist('BiVentricleVolume', 'var')
%     save BiVentricleVolume BiVentricleVolume;
% end
load BiVentricleVolume;
cd(workingDir);

%%find out the node sets
nodeSets = abaqusInput.nodeSets;
for i = 1 : size(nodeSets,1)
    if strcmp( nodeSets(i).str_node_set, 'NODE_LV_ENDO')
        LVendoNodes = nodeSets(i).nodelist;
    end
    if strcmp( nodeSets(i).str_node_set, 'NODE_EPI' )
        BiVenEpiNodes = nodeSets(i).nodelist;
    end
end
%% now need to find out the nodes for RV_SEPTUM_ENDO
node_RV_SEPTUM_ENDO = extract_node_from_surface(abaqusInput, 'SURF_RV_SEPTUM_ENDO');
BiVen_EPI_RVSeptumENDO = unique([BiVenEpiNodes; node_RV_SEPTUM_ENDO]);
node_assign = LV_RV_assignment.node_assign;
LV_EPI_RVSeptumENDO = [];
for i = 1 : length(BiVen_EPI_RVSeptumENDO)
   nodeId =  BiVen_EPI_RVSeptumENDO(i);
   if node_assign(nodeId) == 1 
       LV_EPI_RVSeptumENDO = [LV_EPI_RVSeptumENDO; nodeId]; %% only keep the LV 
   end
end

optimize_opt.abaqus_dis_out_filename = abaqus_dis_out_filename; 
optimize_opt.abaqusSimulationDir = abaqusSimulationDir;
optimize_opt.abaqusInput = abaqusInput;
optimize_opt.cirDir = cirDir;
optimize_opt.radDir = radDir;
optimize_opt.fibreDir = fibreDir;
optimize_opt.sheetDir = sheetDir;
optimize_opt.elRegionsFull = elRegionsFull;
optimize_opt.nodeElRegionFull = nodeElRegionFull;
optimize_opt.LV_RV_assignment = LV_RV_assignment; 
optimize_opt.cpu_number = cpu_number;
optimize_opt.subroutine_file = subroutine_file;
optimize_opt.abaqus_inputfile = abaqus_inputfile;
optimize_opt.LVendoNodes = LVendoNodes;
optimize_opt.LV_EPI_RVSeptumENDO = LV_EPI_RVSeptumENDO;
optimize_opt.abaqus_inputfile_original = abaqus_inputfile_original;
optimize_opt.materialParam_startLine_LV = materialParam_startLine_LV;
optimize_opt.materialParam_startLine_RV = materialParam_startLine_RV;
optimize_opt.pressureParam_startLine = pressureParam_startLine;
optimize_opt.logfile_name = logfile_name;


%% now a test run 
A = 0.2244;
B=1.6215;
Af=2.4267;
Bf=1.8268;
An=0.5562;
Bn=0.7746;
Afs=0.3905;
Bfs=1.695;
Ca_RV=2.0;
press = 8.0;

%% regenerate updated abaqus input file
abaqusInputFileUpdate_MatModel(optimize_opt, A, B, Af, Bf, An, Bn, Afs, Bfs, ...
                                        Ca_RV, press);
[status, result] = RunAbaqusJobFromMatlab(abaqusSimulationDir, abaqus_inputfile, ...
                      cpu_number, subroutine_file);
SuccessB = AbaqusRunningSuccessOrNot(abaqusSimulationDir, abaqus_inputfile);

if (SuccessB == 1)
    %%% now extract the displacement fields
    DislacementExtractByPython(abaqusSimulationDir,abaqus_inputfile,...
                            pythonOriginalFilesDir,pythonfilename, ...
                            lineNoForODBName, lineNoForDisName,...
                            abaqus_dis_out_filename)

    %% volume calculation
    ResVol = biVenCavityCalculation(optimize_opt);

    %% strain calculation
    using_fsn = 0; %% whether to project the strain along myofibre direction or not
    fiberStrain_crl = strainCalculationFromNodalDisplacements(optimize_opt, using_fsn);

    strain_segs= segRegionsStrainSummarization(optimize_opt, fiberStrain_crl, LVendoNodes);
    strain_segs_LVEpi= segRegionsStrainSummarization(optimize_opt, fiberStrain_crl, LV_EPI_RVSeptumENDO);
    
    %% output to the log file 
    cd(abaqusSimulationDir);
    fid_log = fopen(optimize_opt.logfile_name,'a');
    cd(workingDir);
    
    timestr =  datestr(clock());
    fprintf(fid_log, 'running a forward simulation\n');
    fprintf(fid_log, 'Step running on %s\n', timestr);
    fprintf(fid_log, 'abaqus running success: %d\n', SuccessB);
    fprintf(fid_log, 'parameters used for LV: %f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f, \t%f mmHg\n', A, B, Af, Bf, An, Bn, Afs, Bfs, press);
    fprintf(fid_log, 'parameters used for RV: %f\n',Ca_RV);
%     fprintf(fid_log, 'LV volume: %f(target: %f)\n', strainComparison.LVVolumeAba, strainComparison.LVVolumeMRI);
%     %fprintf(fid_log, 'LV strain difference squared: %f\n', sum(feval_total.^2));
%     fprintf(fid_log, 'LV strain difference squared: %f using %d segments (ave: %f vs %f) \n', sum(feval_total.^2)^0.5, length(strainComparison.strainMRITotal), ...
%                                    mean(strainComparison.strainMRITotal),mean(strainComparison.strainAbaTotal) );
   fclose(fid_log);


    
    figure; hold on;
    plot(strain_segs(7:30), 'r-'); hold on;
    plot(strain_segs_LVEpi(7:30), 'b-');
    plot(  (strain_segs(7:30) + strain_segs_LVEpi(7:30))/2, 'b--');
    legend('endo', 'epi', 'average');
    
    %% update the BiVentriculeVolume data
    BiVentricleVolume.RV_endo_correction =  ResVol.RV_vol_ori - BiVentricleVolume.RV_early_diastole;
    cd(resultDir)
    save BiVentricleVolume BiVentricleVolume;
    cd(workingDir);
end

