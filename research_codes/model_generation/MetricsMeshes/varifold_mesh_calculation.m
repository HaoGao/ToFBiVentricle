clear all; close all; clc
%% calculate teh varifolds distance for two poly surface mesh
load polyMeshData;
sigma = 10; %% kernel width

MaxIndex = size(polyMeshData,1);
for meshIndex = 1 : size(polyMeshData,1)

    nodes_1 = polyMeshData(MaxIndex).polyMesh.nodes;
    elems_1 = polyMeshData(MaxIndex).polyMesh.elems;

    nodes_2 = polyMeshData(meshIndex).polyMesh.nodes;
    elems_2 = polyMeshData(meshIndex).polyMesh.elems;

    elems_1 = elems_1(:, 2:4);
    elems_2 = elems_2(:, 2:4);

    source.nodes = nodes_1;
    source.elems = elems_1;

    target.nodes = nodes_2;
    target.elems = elems_2;

    

    S_dis = varifold_distance(source, target, sigma);
    
    S_dis_seq(meshIndex) = S_dis.dis;
    
    meshIndex

end
 

