function write_vtk_tet_volume(filename, resultDir, vtx, elem, cell_val, point_val, field_val)

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
					
fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');					
fprintf(fid,'POINTS %i double\n',num_vtx);					
					
for i=1:num_vtx					
fprintf(fid,'%d %d %d\n',vtx(i,1),vtx(i,2),vtx(i,3));					
end % for					
					
fprintf(fid,'\n');	

fprintf(fid,'CELLS %i %i\n',num_elem,5*num_elem); % per ogni tetraedro 4 dati + 1 che specifi	ca numpi	nt=4	i	n	totale 5 vedi pag 357				
					
for i=1:num_elem					
fprintf(fid,'%i %i %i %i %i\n',4,elem(i,1)-1,elem(i,2)-1,elem(i,3)-1,elem(i,4)-1);					
end % for	

fprintf(fid,'\n');					
					
fprintf(fid,'CELL_TYPES %i\n',num_elem);					
					
for i=1:num_elem					
fprintf(fid,'10\n');					
end % for

% if ~isempty(cell_val)
%     
% end

if ~isempty(cell_val) 
        fprintf(fid,'CELL_DATA %i\n',num_elem);					
        fprintf(fid,'SCALARS cell_val double\n');					
        fprintf(fid,'LOOKUP_TABLE default\n');					
        for j=1:num_elem					
            fprintf(fid,'%f\n',cell_val(j));					
        end % for
end

if ~isempty(point_val)
    fprintf(fid,'POINT_DATA %i\n',num_vtx);					
    fprintf(fid,'SCALARS point_val double\n');					
    fprintf(fid,'LOOKUP_TABLE default\n');					

    for i=1:num_vtx					
    fprintf(fid,'%f\n',point_val(i));					
    end % for
end


%%output field value
if ~isempty(field_val)
    fprintf(fid, 'FIELD FieldData %d\n', size(field_val,1));
    for i = 1 : size(field_val,1)
       field_name = field_val(i,1).field_name;
       field_data = field_val(i,1).field_data;
        fprintf(fid,'%s 1 %i double\n',field_name, num_elem);									
        for j=1:num_elem					
            fprintf(fid,'%d\n',field_data(j));					
        end % for
    end
end

fclose(fid);




