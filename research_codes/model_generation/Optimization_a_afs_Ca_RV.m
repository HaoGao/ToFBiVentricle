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


cd(resultDir);
load BiVentricleVolume;
cd(workingDir);

%%% read the strain data from images 
cd(resultDir);
cd(bSplineDir);
bSplineDir = pwd();
cd(workingDir);

cd(bSplineDir);
fid_strainMRI = fopen(straininvivoMRI_filename,'r');
cd(workingDir);
%%%read strain from MRI measurement
strainData_struct = readStrainAfterBsplineRecovery(fid_strainMRI);
fclose(fid_strainMRI);
strainData = [];
strainDataFlag = [];
for i = 1 : size(strainData_struct,2)
    strainData = [strainData; strainData_struct(1,i).segStrain];
    strainDataFlag = [strainDataFlag; strainData_struct(1,i).segStrainB];
end



global optimize_opt;
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
optimize_opt.pythonOriginalFilesDir = pythonOriginalFilesDir;
optimize_opt.pythonfilename=pythonfilename;
optimize_opt.lineNoForODBName=lineNoForODBName;
optimize_opt.lineNoForDisName=lineNoForDisName;
optimize_opt.BiVentricleVolume = BiVentricleVolume;
optimize_opt.strainData = strainData;
optimize_opt.strainDataFlag = strainDataFlag;
optimize_opt.LVEDP_High = 20;
%optimize_opt.mpara = mpara;

%% from the previous step, update the parameters

cd(abaqusSimulationDir);
load step_Af_Bf_opt;
cd(workingDir);
mpara.A_opt = A_opt;
mpara.B_opt = B_opt;
mpara.Af_opt = Af_opt;
mpara.Bf_opt = Bf_opt;
mpara.An_opt = An_opt;
mpara.Bn_opt = Bn_opt;
mpara.Afs_opt = Afs_opt;
mpara.Bfs_opt = Bfs_opt;
mpara.Ca_RV_opt = Ca_RV_opt;
optimize_opt.mpara = mpara;


options_Opt = optimset('Algorithm', 'sqp', 'TolFun', 1e-6, ...
                           'TolX',0.000001,'Diffminchange',1e-4,'Diffmaxchange',1e-3, 'MaxIter', 100);
amin = 0.1;
amax = 5;
bmin = 0.1;
bmax = 5;

Lb = [amin bmin amin];
Ub = [amax bmax amax];
x0 = [1.0  1.0 1.0];
[x,fval,exitflag,output] = fmincon(@obj_a_b_Ca_RV,x0,[],[],[],[],Lb,Ub,[],options_Opt);

cd(abaqusSimulationDir);
A_opt = optimize_opt.mpara.A_opt*x(1);
B_opt = optimize_opt.mpara.B_opt*x(2);
Af_opt = optimize_opt.mpara.Af_opt;
Bf_opt = optimize_opt.mpara.Bf_opt;
An_opt = optimize_opt.mpara.An_opt;
Bn_opt = optimize_opt.mpara.Bn_opt;
Afs_opt = optimize_opt.mpara.Afs_opt*x(1);
Bfs_opt = optimize_opt.mpara.Bfs_opt*x(2);
Ca_RV_opt = optimize_opt.mpara.Ca_RV_opt*x(3);
save step_a_b_CaRV_opt A_opt B_opt Af_opt Bf_opt An_opt Bn_opt Afs_opt Bfs_opt Ca_RV_opt;
cd(workingDir);





