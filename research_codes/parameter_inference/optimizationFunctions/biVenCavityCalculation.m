function ResVol = biVenCavityCalculation(optimize_opt)

%%% calculate the volume cavity for both LV and RV
abaqusInput = optimize_opt.abaqusInput;
nodeSets = abaqusInput.nodeSets;
nodesT = abaqusInput.nodes;
nodes = nodesT(:,2:4);
clear nodesT;
for i = 1 : size(nodeSets, 1)
    if strcmp(nodeSets(i).str_node_set, 'NODE_LV_ENDO') 
        endoNodes= nodeSets(i).nodelist;
        [LV_vol_ori,LV_vol_update] = LVCavityVolumeCalculation(nodes, endoNodes,... 
            optimize_opt.abaqus_dis_out_filename,... 
            optimize_opt.abaqusSimulationDir);
        %%convert into mL
        LV_vol_ori = LV_vol_ori/1000;
        LV_vol_update = LV_vol_update/1000;
    end
    
    if strcmp(nodeSets(i).str_node_set, 'NODE_EPI')
       epiNodes = nodeSets(i).nodelist;
       [BiEPI_vol_ori,BiEPI_vol_update] = LVCavityVolumeCalculation(nodes, epiNodes,... 
            optimize_opt.abaqus_dis_out_filename,... 
            optimize_opt.abaqusSimulationDir);
        %%convert into mL
        BiEPI_vol_ori = BiEPI_vol_ori/1000;
        BiEPI_vol_update = BiEPI_vol_update/1000;
    end
     
     %%% now we will calculate the RV endo + LV cavity + septum   
    if strcmp(nodeSets(i).str_node_set, 'NODE_RV_ENDO')
       RVendoNodes = nodeSets(i).nodelist;
       [RV_vol_ori,RV_vol_update] = LVCavityVolumeCalculation(nodes, RVendoNodes,... 
            optimize_opt.abaqus_dis_out_filename,... 
            optimize_opt.abaqusSimulationDir);
        %%convert into mL
        RV_vol_ori = RV_vol_ori/1000;
        RV_vol_update = RV_vol_update/1000;
    end
    
    
end


% %%% this is to calcualte the LV cavity + septum
% node_RV_SEPTUM_ENDO = extract_node_from_surface(abaqusInput, 'SURF_RV_SEPTUM_ENDO');
% LV_RVSEPTUM_endo = [endoNodes; node_RV_SEPTUM_ENDO];
% [LV_Septum_vol_ori,LV_Septum_vol_update] = LVCavityVolumeCalculation(nodes, LV_RVSEPTUM_endo,... 
%             optimize_opt.abaqus_dis_out_filename,... 
%             optimize_opt.abaqusSimulationDir);
%         %%convert into mL
% LV_Septum_vol_ori = LV_Septum_vol_ori/1000;
% LV_Septum_vol_update = LV_Septum_vol_update/1000;
%         
% LV_RV_endo = [RVendoNodes; endoNodes];
% [RV_LV_vol_ori,RV_LV_vol_update] = LVCavityVolumeCalculation(nodes, LV_RV_endo,... 
%             optimize_opt.abaqus_dis_out_filename,... 
%             optimize_opt.abaqusSimulationDir);
%     %%convert into mL
% RV_LV_vol_ori = RV_LV_vol_ori/1000;
% RV_LV_vol_update = RV_LV_vol_update/1000;     
%  
% 
% RV_cavity_vol_ori = RV_LV_vol_ori - LV_Septum_vol_ori;
% RV_cavity_vol_update = RV_LV_vol_update - LV_Septum_vol_update;

%% now we will need to calculate the volume of all elements, a new function is needed here
elemsT = abaqusInput.elems;
elems = elemsT(:, 2:5);
clear elemsT;
[wall_vol_ori, wall_vol_update] = wallVolumeCalculation(nodes, elems,...
    optimize_opt.abaqus_dis_out_filename,... 
            optimize_opt.abaqusSimulationDir);
wall_vol_ori = wall_vol_ori/1000;
wall_vol_update = wall_vol_update/1000;



ResVol.LV_vol_ori = LV_vol_ori;
ResVol.LV_vol_update = LV_vol_update;
ResVol.EPI_vol_ori = BiEPI_vol_ori;
ResVol.RV_vol_ori = RV_vol_ori;
ResVol.RV_vol_update = RV_vol_update;
ResVol.wall_vol_ori = wall_vol_ori;
ResVol.wall_vol_update = wall_vol_update;

