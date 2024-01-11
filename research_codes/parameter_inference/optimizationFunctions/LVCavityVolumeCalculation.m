function [vol_ori,vol_update] = LVCavityVolumeCalculation(nodes, endoNodes,... 
    abaqus_dis_out_filename, abaqusDir)

workingDir = pwd();
cd(abaqusDir);
displacement = load(abaqus_dis_out_filename);
% re order for explicit file
idx=[];
[~,idx] = sort(displacement(:,1)); % sort just the first column
displacement = displacement(idx,:);
cd(workingDir);



nodes_x = nodes(endoNodes',1);
nodes_y = nodes(endoNodes',2);
nodes_z = nodes(endoNodes',3);

dis_x = displacement(endoNodes',2);
dis_y = displacement(endoNodes',3);
dis_z = displacement(endoNodes',4);

nodes_x_update = nodes_x + dis_x;
nodes_y_update = nodes_y + dis_y;
nodes_z_update = nodes_z + dis_z;

%%%%% modification by DB 02/07/2023


[~, vol_ori] = boundary([nodes_x nodes_y nodes_z],0.9);

[~, vol_update] = boundary([nodes_x_update nodes_y_update nodes_z_update],0.9);



% [~, vol_ori] = convexHull(DT);
% [~, vol_update] = convexHull(DT_update);


% DT =  delaunayTriangulation(nodes_x,nodes_y,nodes_z);
% DT_update =  delaunayTriangulation(nodes_x_update,nodes_y_update,nodes_z_update);
% % figure
% % tetramesh(DT,'FaceColor','r','FaceAlpha',0.3);
% 
% % Points = DT.Points;
% % ConnectivityList = DT.ConnectivityList;
% % 
% % NoOfElems = size(ConnectivityList,1);
% 
% % vol = 0.0;
% % for eleIndex = 1 : NoOfElems
% %     a = Points(ConnectivityList(eleIndex,1),1:3);
% %     b = Points(ConnectivityList(eleIndex,2),1:3);
% %     c = Points(ConnectivityList(eleIndex,3),1:3);
% %     d = Points(ConnectivityList(eleIndex,4),1:3);
% %     
% %     vmax = [a-b; b-c; c-d];
% %     vol = vol + 1/6.0*abs(det(vmax));
% % end
% 
% [~, vol_ori] = convexHull(DT);
% [~, vol_update] = convexHull(DT_update);