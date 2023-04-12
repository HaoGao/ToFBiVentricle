function write_vtk_trigular_surface(filename, resultDir, vtx, elem, val)

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




