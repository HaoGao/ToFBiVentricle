function TecplotMeshGenWithFileName(Ver, Mesh,Dis,resultDir,filename)

workingDir = pwd();
cd(resultDir);
fileName = sprintf('TecplotTet%s.dat',filename);
% fid = fopen('TecplotTetHelixAngle.dat','w');
fid = fopen(fileName,'w');
cd(workingDir)
fprintf(fid, 'TITLE = "Heart MESH" \n');
fprintf(fid, 'VARIABLES = "x", "y", "z", "u", "v", "w" \n');
fprintf(fid, 'ZONE T="TotalMesh", N = %d, E=%d, F=FEPOINT, ET=TETRAHEDRON \n', size(Ver, 1), size(Mesh,1));


for i = 1 : size(Ver,1)
    fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%f\n',Ver(i,2),Ver(i,3),Ver(i,4),Dis(i,1),Dis(i,2),Dis(i,3));
end

for i = 1 : size(Mesh,1)
    fprintf(fid,'%d\t%d\t%d\t%d\n', Mesh(i,2),Mesh(i,3),Mesh(i,4),Mesh(i,5));
end



fclose(fid);