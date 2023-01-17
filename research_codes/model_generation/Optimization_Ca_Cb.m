clear all; close all; clc;

% LVWM_config;
LVWM_config;


%%need to figure out the end-diastolic folder and load the mesh
cd(resultDir);
cd('end_diastole');
cd(solidworksDir);
cd(gmeshDir);
end_diastolic_gmeshDir = pwd();
cd(workingDir);
cd(end_diastolic_gmeshDir);
load abaqusInput.mat;
abaqusInput_end_diastole = abaqusInput;
clear abaqusInput;
cd(workingDir);

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
optimize_opt.mpara = mpara;
optimize_opt.weighting_varifold = 0.01;


%% the varifold distance 
disp('calculate the vari_dis between early-diastole to end-diastole, be patient...\n');
[vari_dis, var_dis_LV_endo] = varifold_distance_abaqus(abaqusInput, abaqusInput_end_diastole, []);
optimize_opt.vari_dis_max = vari_dis;


%% now a test run 
A0 = optimize_opt.mpara.A0;
B0=optimize_opt.mpara.B0;
Af0=optimize_opt.mpara.Af0;
Bf0=optimize_opt.mpara.Bf0;
An0=optimize_opt.mpara.An0;
Bn0=optimize_opt.mpara.Bn0;
Afs0=optimize_opt.mpara.Afs0;
Bfs0=optimize_opt.mpara.Bfs0;
Ca_RV0=optimize_opt.mpara.Ca_RV0;
press0 = optimize_opt.mpara.press0;

lb_Ca = 0.1; ub_Ca = 2.0;
lb_Cb = 0.1; ub_Cb = 1.0;

Ca_seq = lb_Ca : 0.2 : ub_Ca;
Cb_seq = lb_Cb : 0.1 : ub_Cb;

simIndex = 0;
for caIndex = 1 : length(Ca_seq)
    for cbIndex = 1 : length(Cb_seq)
        Ca = Ca_seq(caIndex);
        Cb = Cb_seq(cbIndex);
        A = A0*Ca;
        Af = Af0*Ca;
        An = An0*Ca;
        Afs = Afs0*Ca;
        B = B0*Cb;
        Bf = Bf0*Cb;
        Bn = Bn0*Cb;
        Bfs = Bfs0*Cb;
        Ca_RV = Ca_RV0*1.0;
        press = press0*1.0;
        
        cd(optimize_opt.abaqusSimulationDir);
        fid_log = fopen(optimize_opt.logfile_name, 'a');
        cd(workingDir);
        fprintf(fid_log, 'Ca: %f\t Cb: %f\n', Ca, Cb);
        fclose(fid_log);

        ResSimT = BiVenPassiveForwardSimulationMatPres(optimize_opt, ...
                                            A, B, Af, Bf, An, Bn, Afs, Bfs, ...
                                           Ca_RV, press);

       simIndex = simIndex + 1;
       
       
       %%output the strain data
       cd(optimize_opt.abaqusSimulationDir);
       fid_log = fopen(optimize_opt.logfile_name, 'a');
       cd(workingDir);
       
       fprintf(fid_log, 'abaqus strain data at endo \n');
       NStrains = length(strainData);
       diff_strain_Endo = [];
       for i = 1 : length(NStrains)
           if ~isnan(ResSimT.strain_segs(i))
            fprintf(fid_log, '%f\t', ResSimT.strain_segs(i) );
            diff_strain_Endo = [diff_strain_Endo ResSimT.strain_segs(i) - strainData(i)];
           end
       end
       fprintf(fid_log, '\n');	
       fprintf(fid_log, 'abaqus strain data at epi \n');
       diff_strain_Epi = [];
       for i = 1 : length(NStrains)
           if ~isnan(ResSimT.strain_segs_LVEpi(i))
            fprintf(fid_log, '%f\t', ResSimT.strain_segs_LVEpi(i) );
             diff_strain_Epi = [diff_strain_Epi ResSimT.strain_segs_LVEpi(i) - strainData(i)];
           end
       end
       fprintf(fid_log, '\n');	
       fprintf(fid_log, 'strain difference at endo: %f( %f)\n', mean(diff_strain_Endo), std(diff_strain_Endo));
       fprintf(fid_log, 'strain difference at epi: %f( %f)\n', mean(diff_strain_Epi), std(diff_strain_Epi));
       
       %%output volume 
       fprintf(fid_log, 'LV_Volume: %f mL (original: %f mL), target: %f \n', ...
           ResSimT.ResVol.LV_vol_update, ResSimT.ResVol.LV_vol_ori, BiVentricleVolume.LV_end_diastole);
       fprintf(fid_log, 'RV_Volume: %f mL (original: %f mL), target: %f\n', ...
           ResSimT.ResVol.RV_vol_update, ResSimT.ResVol.RV_vol_ori, BiVentricleVolume.RV_end_diastole);
       
       %%output the difference with measured strain data and volume
       
       %% the varifold distance
       %% here we can try to use the varifold_distance to define the distance 
        cd(abaqusSimulationDir);
        dis = load(optimize_opt.abaqus_dis_out_filename); 
        cd(workingDir);
        [vari_dis, vari_dis_LV_endo] = varifold_distance_abaqus(abaqusInput, abaqusInput_end_diastole, dis);
        fprintf(fid_log, 'vari_dis: %f (max: %f), ratio: \n', ...
            vari_dis.dis, optimize_opt.vari_dis_max.dis, (vari_dis.dis)/(optimize_opt.vari_dis_max.dis) );

       
        
       fprintf(fid_log, '\n');
       fprintf(fid_log, '\n\n');
       fclose(fid_log);
       
%         figure; hold on;
%         plot(strain_segs(7:30), 'r-'); hold on;
%         plot(strain_segs_LVEpi(7:30), 'b-');
%         plot(  (strain_segs(7:30) + strain_segs_LVEpi(7:30))/2, 'b--');


        ResSim(simIndex,1).result = ResSimT;
        ResSim(simIndex,1).vari_dis = vari_dis;
        
    end
end

%% save the results
 cd(optimize_opt.abaqusSimulationDir);
 save ResSim_CaCb ResSim;
 cd(workingDir);
 
 cd(optimize_opt.abaqusSimulationDir);
 load ResSim_CaCb;
 cd(workingDir);
 
 
 %% now figure out the best Ca and Cb for potential second step optimization
simIndex = 0;
fitting_error = [];
Ca_grid = []; 
Cb_grid = [];
mse_vol = []; mse_varifold=[]; mse_total=[];
for caIndex = 1 : length(Ca_seq)
    for cbIndex = 1 : length(Cb_seq)
        simIndex = simIndex + 1;
        CaPara(simIndex,1) = Ca_seq(caIndex);
        CbPara(simIndex,1) = Cb_seq(cbIndex);
        Ca_grid(caIndex,cbIndex) = Ca_seq(caIndex);
        Cb_grid(caIndex,cbIndex) = Cb_seq(cbIndex);
        
        LV_vol_ori(simIndex,1) = ResSim(simIndex,1).result.ResVol.LV_vol_ori;
        LV_vol_update(simIndex,1) = ResSim(simIndex,1).result.ResVol.LV_vol_update;
        
        RV_vol_ori(simIndex,1) = ResSim(simIndex,1).result.ResVol.RV_vol_ori;
        RV_vol_update(simIndex,1) = ResSim(simIndex,1).result.ResVol.RV_vol_update;
    
        
       varifoldDis(simIndex,1) = ResSim(simIndex,1).vari_dis.dis/(optimize_opt.vari_dis_max.dis);
       
       %% here we will consider the left side and total varifold distance    
       fitting_error(simIndex,1) =  (1 - LV_vol_update(simIndex,1)/(BiVentricleVolume.LV_end_diastole) )^2  + ...
                             0.01*varifoldDis(simIndex,1)^2; 
       fitting_error(simIndex,2) =  (1 - LV_vol_update(simIndex,1)/(BiVentricleVolume.LV_end_diastole) )^2;
       fitting_error(simIndex,3) =  varifoldDis(simIndex,1)^2;
       
       mse_vol(caIndex, cbIndex) = fitting_error(simIndex,2);
       mse_varifold(caIndex, cbIndex) = fitting_error(simIndex,3);
       mse_total(caIndex, cbIndex) = fitting_error(simIndex,1);
       
       
    end
end

[fitting_error_sort, ind] = sort(fitting_error(:,1));
fitting_error_sort(:,2) = CaPara(ind);
fitting_error_sort(:,3) = CbPara(ind);
%% find out the minimal error with largest Ca and minimal Cb
Ca_optimal = CaPara(ind(1));
Cb_optimal = CbPara(ind(1));
%% there will be other better way to do it
% for i = 1 : length(fitting_error_sort)
%     
% end
 
 cd(optimize_opt.abaqusSimulationDir);
 save ResSim_CaCb ResSim Ca_optimal Cb_optimal fitting_error;
 cd(workingDir);
 
 
 
 
 %% here generate contour plot
%  [Ca_grid_t, Cb_grid_t] = meshgrid(Ca_seq,Cb_seq);
 hcontur = figure; hold on;
 subplot(221)
 contour(Ca_grid, (Cb_grid), log(mse_total) , 'Fill', 'on'); hold on;
 xlabel('C_a');
 ylabel('C_b');
 title('log(mse\_total)');
 %% plot the minimum 
[indi, indj] = find_minimum_index_2D(mse_total);
 plot( Ca_grid(indi,indj), Cb_grid(indi,indj), '*', 'MarkerSize', 6); 
 ResSimT = BiVenPassiveForwardSimulationMatPres(optimize_opt, ...
                     A0*Ca_grid(indi,indj), B0*Cb_grid(indi,indj),...
                     Af0*Ca_grid(indi,indj), Bf0*Cb_grid(indi,indj), ...
                     An0*Ca_grid(indi,indj), Bn0*Cb_grid(indi,indj), ...
                     Afs0*Ca_grid(indi,indj), Bfs0*Cb_grid(indi,indj), ...
                     Ca_RV0*1.0, press0*1.0);
  cd(abaqusSimulationDir);
  dis = load(optimize_opt.abaqus_dis_out_filename); 
  cd(workingDir);
  [vari_dis, vari_dis_LV_endo] = varifold_distance_abaqus(abaqusInput, abaqusInput_end_diastole, dis);
  cd(abaqusSimulationDir);
  mkdir('step1_optim/total_min');
  movefile('*.vtk',  'step1_optim/total_min');
  cd(workingDir);

 
 subplot(223)
 contour(Ca_grid, (Cb_grid), log(mse_vol) , 'Fill', 'on');hold on;
 xlabel('C_a');
 ylabel('C_b');
 title('log(mse\_vol)');
 %% plot the minimum 
 [indi, indj] = find_minimum_index_2D(mse_vol);
 plot( Ca_grid(indi,indj), Cb_grid(indi,indj), '*', 'MarkerSize', 6);
 
 ResSimT = BiVenPassiveForwardSimulationMatPres(optimize_opt, ...
                     A0*Ca_grid(indi,indj), B0*Cb_grid(indi,indj),...
                     Af0*Ca_grid(indi,indj), Bf0*Cb_grid(indi,indj), ...
                     An0*Ca_grid(indi,indj), Bn0*Cb_grid(indi,indj), ...
                     Afs0*Ca_grid(indi,indj), Bfs0*Cb_grid(indi,indj), ...
                     Ca_RV0*1.0, press0*1.0);
  cd(abaqusSimulationDir);
  dis = load(optimize_opt.abaqus_dis_out_filename); 
  cd(workingDir);
  [vari_dis, vari_dis_LV_endo] = varifold_distance_abaqus(abaqusInput, abaqusInput_end_diastole, dis);
  cd(abaqusSimulationDir);
  mkdir('step1_optim/vol_min');
  movefile('*.vtk',  'step1_optim/vol_min');
 cd(workingDir);

 subplot(224)
 contour(Ca_grid, (Cb_grid), log(0.01*mse_varifold) , 'Fill', 'on');hold on;
 xlabel('C_a');
 ylabel('C_b');
 title('log(mse\_varifold)');
 %% plot the minimum 
 [indi, indj] = find_minimum_index_2D(mse_varifold);
 plot( Ca_grid(indi,indj), Cb_grid(indi,indj), '*', 'MarkerSize', 6);
 
 ResSimT = BiVenPassiveForwardSimulationMatPres(optimize_opt, ...
                     A0*Ca_grid(indi,indj), B0*Cb_grid(indi,indj),...
                     Af0*Ca_grid(indi,indj), Bf0*Cb_grid(indi,indj), ...
                     An0*Ca_grid(indi,indj), Bn0*Cb_grid(indi,indj), ...
                     Afs0*Ca_grid(indi,indj), Bfs0*Cb_grid(indi,indj), ...
                     Ca_RV0*1.0, press0*1.0);
  cd(abaqusSimulationDir);
  dis = load(optimize_opt.abaqus_dis_out_filename); 
  cd(workingDir);
  [vari_dis, vari_dis_LV_endo] = varifold_distance_abaqus(abaqusInput, abaqusInput_end_diastole, dis);
  cd(abaqusSimulationDir);
  mkdir('step1_optim/varifold_min');
  movefile('*.vtk',  'step1_optim/varifold_min');
  cd(workingDir);
 
 
 
 
 
 
 
 
 
 
 
 



