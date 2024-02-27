clear all; close all; clc;
% 
% % LVWM_config;
  LVWM_config;
% % 
% % %%now need to choose phase to segment
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
cd(meshDir);
meshDir = pwd();
cd(workingDir)

%%% load mesh files
cd(meshDir);
load abaqusInput;
fibreDir = load('fibreDir.txt');
sheetDir = load('sheetDir.txt');
cirDir =   load('circumDir.txt');
radDir =   load('radDir.txt');
load LV_RV_assignment;
%load LVMeshSegDivisions; %% elRegionsFull nodeElRegionFull saved the strain locations according to each slice
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
% node_RV_SEPTUM_ENDO = extract_node_from_surface(abaqusInput, 'SURF_RV_SEPTUM_ENDO');
% BiVen_EPI_RVSeptumENDO = unique([BiVenEpiNodes; node_RV_SEPTUM_ENDO]);
% node_assign = LV_RV_assignment.node_assign;
% LV_EPI_RVSeptumENDO = [];
% for i = 1 : length(BiVen_EPI_RVSeptumENDO)
%    nodeId =  BiVen_EPI_RVSeptumENDO(i);
%    if node_assign(nodeId) == 1 
%        LV_EPI_RVSeptumENDO = [LV_EPI_RVSeptumENDO; nodeId]; %% only keep the LV 
%    end
% end


cd(resultDir);
load BiVentricleVolume;
cd(workingDir);


global optimize_opt;
optimize_opt.abaqus_dis_out_filename = abaqus_dis_out_filename; 
optimize_opt.abaqusSimulationDir = abaqusSimulationDir;
optimize_opt.abaqusInput = abaqusInput;
optimize_opt.cirDir = cirDir;
optimize_opt.radDir = radDir;
optimize_opt.fibreDir = fibreDir;
optimize_opt.sheetDir = sheetDir;
%optimize_opt.elRegionsFull = elRegionsFull;
%optimize_opt.nodeElRegionFull = nodeElRegionFull;
optimize_opt.LV_RV_assignment = LV_RV_assignment; 
optimize_opt.cpu_number = cpu_number;
optimize_opt.abaqus_command = abaqus_command;
optimize_opt.subroutine_file = subroutine_file;
optimize_opt.abaqus_inputfile = abaqus_inputfile;
optimize_opt.LVendoNodes = LVendoNodes;
%optimize_opt.LV_EPI_RVSeptumENDO = LV_EPI_RVSeptumENDO;
optimize_opt.abaqus_inputfile_original = abaqus_inputfile_original;
optimize_opt.materialParam_startLine_LV = materialParam_startLine_LV;
optimize_opt.materialParam_startLine_RV = materialParam_startLine_RV;
optimize_opt.pressureParam_startLine = pressureParam_startLine;
optimize_opt.logfile_name = logfile_name;
optimize_opt.pythonOriginalFilesDir = templateFilesDir;
optimize_opt.pythonfilename=pythonfilename;
optimize_opt.lineNoForODBName=lineNoForODBName;
optimize_opt.lineNoForDisName=lineNoForDisName;
optimize_opt.BiVentricleVolume = BiVentricleVolume;
%optimize_opt.strainData = strainData;
%optimize_opt.strainDataFlag = strainDataFlag;
optimize_opt.LVEDP_High = 20;
%optimize_opt.mpara = mpara;
mpara.press0=BiVentricleVolume.lvfp;

%% from the previous step, update the parameters

cd(abaqusSimulationDir);
load step_Ca_Cb_opt_refine;
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


cd(abaqusSimulationDir);
fid_log = fopen(logfile_name,'a');
cd(workingDir);
fprintf(fid_log, '\n \n beginning of step for af bf\n \n');

options_Opt = optimset('Algorithm', 'sqp', 'TolFun', 1e-6, ...
                           'TolX',0.000001,'Diffminchange',1e-4,'Diffmaxchange',1e-3, 'MaxIter', 50);
% options_Opt = optimset('Algorithm', 'sqp', 'TolFun', 1e-2, ...
%                            'TolX',0.01,'Diffminchange',1e-2,'Diffmaxchange',1e-2, 'MaxIter', 100);

afmin = 0.1;
afmax = 5;
bfmin = 0.1;
bfmax = 5;

Lb = [afmin bfmin ];
Ub = [afmax bfmax ];
x0 = [1.0  1.0];
[x,fval,exitflag,output] = fmincon(@obj_afbf,x0,[],[],[],[],Lb,Ub,[],options_Opt);

cd(abaqusSimulationDir);
A_opt = optimize_opt.mpara.A_opt;
B_opt = optimize_opt.mpara.B_opt;
Af_opt = optimize_opt.mpara.Af_opt*x(1);
Bf_opt = optimize_opt.mpara.Bf_opt*x(2);
An_opt = optimize_opt.mpara.An_opt;
Bn_opt = optimize_opt.mpara.Bn_opt;
Afs_opt = optimize_opt.mpara.Afs_opt;
Bfs_opt = optimize_opt.mpara.Bfs_opt;
Ca_RV_opt = optimize_opt.mpara.Ca_RV_opt;
save step_Af_Bf_opt A_opt B_opt Af_opt Bf_opt An_opt Bn_opt Afs_opt Bfs_opt Ca_RV_opt;
cd(workingDir);

