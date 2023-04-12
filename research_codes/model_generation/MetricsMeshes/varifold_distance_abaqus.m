function [vari_dis, vari_dis_selected] = varifold_distance_abaqus(abaqusInput, abaqusInput_end_diastole, dis)
sigma = 10; %% kernel width

global optimize_opt;

surface_nodes_ed = abaqusInput_end_diastole.surface_nodes; 
surface_elements_ed = abaqusInput_end_diastole.surface_elements;

surface_nodes_ori = abaqusInput.surface_nodes;
surface_elements_ori = abaqusInput.surface_elements;
surface_LV_endo_idx = abaqusInput.surface_LV_endo_idx; 


%%need to move surface_nodes_ed to be aligned at z=0 for the basal plane
z_max = max(surface_nodes_ed(:,3));
surface_nodes_ed(:,3) = surface_nodes_ed(:,3) - z_max;

target.nodes = surface_nodes_ed;
target.elems = surface_elements_ed;


surface_nodes_map = abaqusInput.surface_nodes_map;
source.nodes = surface_nodes_ori;
source.elems = surface_elements_ori;
%%update the displacements to the original mesh 
dxdydz = [];
if size(dis, 1)>1
    disT(:,1:3) = dis(:, 2:4); 
    dis = disT;
    nodeT = optimize_opt.abaqusInput.nodes;
    node(:, 1) = nodeT(:,2);
    node(:, 2) = nodeT(:,3);
    node(:, 3) = nodeT(:,4);
    
    %% now regenerate the new surface nodes
    for i = 1 : size(surface_nodes_map,1)
         i_local_surf = surface_nodes_map(i,1);
         i_global = surface_nodes_map(i,2);
         if i_local_surf >1.0e-6 && i_global>1.0e-6
            dxdydz(i_local_surf,1:3) = dis(i_global,1:3);
         end
    end
    
    source.nodes = surface_nodes_ori + dxdydz;
end


abaqusDir = optimize_opt.abaqusSimulationDir;
write_vtk_trigular_surface('ori_viven.vtk', abaqusDir, surface_nodes_ori, surface_elements_ori, []);
write_vtk_trigular_surface('ori_viven_deform.vtk', abaqusDir, source.nodes, surface_elements_ori, []);
write_vtk_trigular_surface('end_diastolic_viven.vtk', abaqusDir, surface_nodes_ed, surface_elements_ed, []);

vari_dis = varifold_distance(source, target, sigma);
vari_dis_selected = varifold_distance_selected(source, target, sigma,surface_LV_endo_idx);
