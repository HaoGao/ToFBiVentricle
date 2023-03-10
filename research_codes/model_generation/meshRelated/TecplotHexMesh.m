function TecplotHexMesh(Ver, Mesh,segRegions,fid)

if isempty(segRegions)
    segRegions = zeros([size(Ver,1), 1]);
end


fprintf(fid, 'TITLE = "Heart MESH maped with nodal value" \n');
fprintf(fid, 'VARIABLES = "x", "y", "z", "u", "v", "w" \n');
fprintf(fid, 'ZONE T="Fitted LV Mesh", N=%d, E=%d, F=FEPOINT, ET=BRICK\n', size(Ver, 1), size(Mesh,1));


for i = 1 : size(Ver,1)
    fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%f\n',Ver(i,2),Ver(i,3),Ver(i,4),segRegions(i),0,0);
end

for i = 1 : size(Mesh,1)
     fprintf(fid, '%d   %d    %d    %d    %d    %d     %d    %d\n', Mesh(i,1),Mesh(i,2), ...
         Mesh(i,3),Mesh(i,4),Mesh(i,5),Mesh(i,6),Mesh(i,7),Mesh(i,8));
end
