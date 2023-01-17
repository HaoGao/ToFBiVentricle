function VTKSurfaceOutput()

clear all; close all; clc;

workingDir = pwd();
mesh_dir = 'C:\Users\sharp\Desktop\work\LVReconstructionAI\Results\HV01\earlyDiastole';
cd(mesh_dir);
load abaqusInputData;
cd(workingDir);


node_xyz = abaqusInputData.node(:, 1:3);
elem = abaqusInputData.elem;
endo_faces = abaqusInputData.endofaces;
epifaces = abaqusInputData.epifaces;
basefaces = abaqusInputData.basefaces;

surface_faces = [basefaces; endo_faces; epifaces];
tri_faces = triangle_surface_extract(surface_faces, elem);

filename = 'HV01.vtk';
resultDir = mesh_dir;
write_vtk_trigular_surface(filename, resultDir, node_xyz, tri_faces);



function tri_faces = triangle_surface_extract(surface_faces, elem)

Nfaces = size(surface_faces, 1);
tri_faces = zeros([2*Nfaces, 3]);

for i = 1 : Nfaces
    faceT = surface_faces(i,:);
    el = faceT(1); faceID = faceT(2);
    
    ndlist = elem(el,1:8);
    
    nlist = [];
    if faceID == 1
        nlist= ndlist([1, 2, 3, 4]);
    elseif faceID == 2
        nlist = ndlist([5 8 7 6]);
    elseif faceID == 3
        nlist = ndlist([1 5 6 2]);
    elseif faceID == 4
        nlist = ndlist([2 6 7 3]);
    elseif faceID == 5
        nlist = ndlist([3 7 8 4]);
    elseif faceID == 6
        nlist = ndlist([4 8 5 1]);
    end
    
    tri_faces(2*(i-1)+1, :) = nlist([1 2 4]);
    tri_faces(2*(i-1)+2, :) = nlist([2 3 4]);
end



function write_vtk_trigular_surface(filename, resultDir, vtx, elem)

% filename = 'biven_surface.vtk';
% resultDir = pwd();
% val = zeros([size(elem,1), 1]);

workingDir = pwd();
num_vtx=size(vtx,1);					
num_elem=size(elem,1);					
					
cd(resultDir);					
fid=fopen(filename,'wt');	
cd(workingDir);
if fid==-1					
error('Cannot open vtk file');					
end % if FID					
					
frewind(fid);					
					
fprintf(fid,'# vtk DataFile Version 3.0\n');					
fprintf(fid,'%s %s exported from MATLAB\n',date,filename);					
fprintf(fid,'ASCII\n');

% Here we write the verticies					
					
fprintf(fid,'DATASET POLYDATA\n');					
fprintf(fid,'POINTS %i float\n',num_vtx);					
					
for i=1:num_vtx					
fprintf(fid,'%f %f %f\n',vtx(i,1),vtx(i,2),vtx(i,3));					
end % for					
					
%fprintf(fid,'\n');	

fprintf(fid,'POLYGONS %i %i\n',num_elem,4*num_elem);					
					
for i=1:num_elem					
fprintf(fid,'%i %i %i %i \n',3,elem(i,1)-1,elem(i,2)-1,elem(i,3)-1);					
end % for

% fprintf(fid,'\n');					
					
%fprintf(fid,'CELL_TYPES %i\n',num_elem);					
					
%for i=1:num_elem					
%fprintf(fid,'5\n');					
%end % for	


%fprintf(fid,'CELL_DATA %i\n',num_elem);					
%fprintf(fid,'SCALARS MatlabExportedScalars double\n');					
%fprintf(fid,'LOOKUP_TABLE default\n');					
					
%for i=1:num_elem					
%fprintf(fid,'%d\n',val(i));					
%end % for

fclose(fid);


