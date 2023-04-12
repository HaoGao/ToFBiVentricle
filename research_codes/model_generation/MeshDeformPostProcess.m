clear all;
close all;
clc;

% LVWM_config;
LVWM_config;

%%now need to choose phase to segment
cd(resultDir);
cd(deformetricaDir);
deformetricaDir = pwd();
cd(workingDir);

cd(resultDir);
cd(deformetrica_libMeshDir);
deformetrica_libMeshDir = pwd();
cd(workingDir);

%%now need to choose phase to segment
list_phase = {'early_diastole', 'end_diastole', 'end_systole'};
[idx, tf] = listdlg('ListString', list_phase);
cd(resultDir);
if ~exist(list_phase{idx},'dir')
    mkdir(list_phase{idx});
    cd(list_phase{idx});
    phase_resultDir = pwd();    
else
    cd(list_phase{idx});
    phase_resultDir = pwd();
end
cd(workingDir);
phase_selected = list_phase{idx};

cd(phase_resultDir);
cd(solidworksDir);
solidworksDir = pwd();
cd(gmeshDir);
gmeshDir = pwd();
cd(workingDir);

cd(gmeshDir); 
load abaqusInput;
cd(workingDir);

for i = 0 : polyFileName_N-1
    str_msg = sprintf('read from %s\n', deformetricaDir);
    disp(str_msg);
    polyFileName = sprintf('%s%d.vtk', polyFileName_base, i);
    polyMesh_T = polyDataReader(polyFileName, deformetricaDir);
    polyMeshData(i+1,1).polyMesh = polyMesh_T;
end

cd(deformetricaDir);
save polyMeshData polyMeshData;
cd(workingDir);

%%here we can add a check whether the surface mesh is still the same as the
%%one in abaqusInput
surface_nodes = abaqusInput.surface_nodes;
surface_elems = abaqusInput.surface_elements;

poly_nodes = polyMeshData(1).polyMesh.nodes;
poly_elems = polyMeshData(1).polyMesh.elems;

%%check the node number 
node_diff = [max(poly_nodes(:) - surface_nodes(:)), ...
             min(poly_nodes(:) - surface_nodes(:))];
elem_diff_1 = [max( surface_elems(:,1) - poly_elems(:,2) ), ...
               min( surface_elems(:,1) - poly_elems(:,2) ) ];
elem_diff_2 = [max( surface_elems(:,2) - poly_elems(:,3) ), ...
               min( surface_elems(:,2) - poly_elems(:,3) ) ]; 
elem_diff_3 = [max( surface_elems(:,3) - poly_elems(:,4) ), ...
               min( surface_elems(:,3) - poly_elems(:,4) ) ];

%%to figure out the displacement field
node_dxdydz = zeros(size(poly_nodes));
for meshIndex = 1 : size(polyMeshData,1)
    if meshIndex == 1
        polyMeshData(meshIndex).polyMesh.dxdydz = node_dxdydz;
        poly_nodes_ref = polyMeshData(meshIndex).polyMesh.nodes;
    else
        poly_nodes_cur = polyMeshData(meshIndex).polyMesh.nodes;
        node_dxdydz = poly_nodes_cur - poly_nodes_ref;
        polyMeshData(meshIndex).polyMesh.dxdydz = node_dxdydz;
    end
end

%% now will need to transfer the displacement into the volumetric mesh,
%% which will be interpolated later by LibMesh
polyMesh_last = polyMeshData(end,1).polyMesh;
nodes = abaqusInput.nodes;
elems = abaqusInput.elems;
surface_nodes_map = abaqusInput.surface_nodes_map; 

dxdydz_surf = polyMesh_last.dxdydz;
dxdydz = zeros([size(nodes,1), 4]);
for i = 1 : size(surface_nodes_map,1)
    i_local_surf = surface_nodes_map(i,1);
    i_global = surface_nodes_map(i,2);
    if i_local_surf >1.0e-6 && i_global>1.0e-6
        dxdydz_surf_T = dxdydz_surf(i_local_surf,:);
        dxdydz(i_global,1:3) = dxdydz_surf_T;
        dxdydz(i_global,4) = 1; % 1 means in the surface
    end
end

point_val = dxdydz(:,4);
write_vtk_tet_volume('deformetric_disp_bool.vtk', deformetrica_libMeshDir, ...
    nodes(:,2:4), elems(:,2:5), [], point_val, []);
write_vtk_tet_volume('deformetric_disp_dz.vtk', deformetrica_libMeshDir, ...
    nodes(:,2:4), elems(:,2:5), [], dxdydz(:,3), []);
%%need to output the displacement in the surface for the volumetric mesh 
cd(deformetrica_libMeshDir);
fid = fopen('deformetric_dis.txt', 'w');
cd(workingDir);
fprintf(fid, '%d\t 0 0 0\n',size(dxdydz,1) );
for i = 1 : size(dxdydz,1)
    %if dxdydz(i,4) ==1 
       fprintf(fid,'%d\t %f\t %f\t %f\t %d\n',i, ...
                    dxdydz(i,1), dxdydz(i,2),dxdydz(i,3),...
                    dxdydz(i,4));
    %end
end
fclose(fid);

cd(deformetrica_libMeshDir);
save polyMeshData polyMeshData dxdydz;
cd(workingDir);

%% how to check the output displacement is right, just plot the surface
% h3D = figure(); hold on;
% set(h3D, 'Visible', 'off'); 
% for i = 1 : size(nodes,1)
%     if dxdydz(i,4) ==1 
%        plot3(nodes(i,2), nodes(i,3), nodes(i,4), 'Marker', '.', 'MarkerSize', 3); 
%     end
% end
% set(h3D, 'Visible', 'on');
% figure(h3D);
% axial equal;

%% another sanity check there
% j = 0;
% for i = 1 : size(nodes,1)
%      if dxdydz(i,4) ==1 
%          j = j+1;
%          node_surfaceT(j,1) =  nodes(i,1);
%      end
% end
% 
% node_surface_T2 = abaqusInput.surface_elems_unsorted;
% node_surface_T2 = unique(node_surface_T2(:));
% node_diff_surf  = [ max(node_surfaceT - node_surface_T2), ...
%                     min(node_surfaceT - node_surface_T2)];
