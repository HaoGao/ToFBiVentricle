function [vol_ori, vol_updated] = wallVolumeCalculation(nodes, elems,...
    abaqus_dis_out_filename, abaqusDir)

workingDir = pwd();
cd(abaqusDir);
displacement = load(abaqus_dis_out_filename);
cd(workingDir);

nodes_x = nodes(:,1);
nodes_y = nodes(:,2);
nodes_z = nodes(:,3);

dis_x = displacement(:,2);
dis_y = displacement(:,3);
dis_z = displacement(:,4);

nodes_x_update = nodes_x + dis_x;
nodes_y_update = nodes_y + dis_y;
nodes_z_update = nodes_z + dis_z;

nodes_update(:,1) = nodes_x_update;
nodes_update(:,2) = nodes_y_update;
nodes_update(:,3) = nodes_z_update;


vol_ori = 0;
vol_updated = 0;

for el = 1 : size(elems, 1)
    node_list = elems(el, 1:4);
    a = nodes(node_list(1), 1:3);
    b = nodes(node_list(2), 1:3);
    c = nodes(node_list(3), 1:3);
    d = nodes(node_list(4), 1:3);
    v = tetra_vol(a, b, c, d);
    vol_ori = vol_ori + v;
    
    
    a = nodes_update(node_list(1), 1:3);
    b = nodes_update(node_list(2), 1:3);
    c = nodes_update(node_list(3), 1:3);
    d = nodes_update(node_list(4), 1:3);
    v_up = tetra_vol(a, b, c, d);
    vol_updated = vol_updated + v_up;
    
    
end



function v = tetra_vol(a, b, c, d)

a = a - d;
b = b - d;
c = c - d;
bc = cross(b, c);
v = abs( dot(a, bc) )/6.0;

%%ref https://en.wikipedia.org/wiki/Tetrahedron