function WriteAbaqusInput_Tet4(abaqusInput, fileName, fileDir)

workingDir = pwd();
Ver=abaqusInput.nodes;
TMesh=abaqusInput.elems;

[~, Ver, TMesh] = nodeMapAdjust(Ver, TMesh);

%%%generate abaqus input file with adjust node number
cd(fileDir);
fid = fopen(fileName,'w');
cd(workingDir);
fprintf(fid, '*************************************   \n');
fprintf(fid, '*HEADING\n');
fprintf(fid,'*************************************   \n');
fprintf(fid,'*Part, name=PART-1\n');
fprintf(fid,'*NODE\n');
%%%node output
for i = 1 : size(Ver,1)
    fprintf(fid, '%d,\t%f,\t%f,\t%f\n', Ver(i,1), Ver(i,2), Ver(i,3), Ver(i,4));
end
%%%element output
fprintf(fid,'*Element, type=C3D4\n');
for i = 1 : size(TMesh,1)
    fprintf(fid,'%d,\t%d,\t%d,\t%d,\t%d\n', TMesh(i,1),TMesh(i,2),TMesh(i,3),TMesh(i,4),TMesh(i,5));
end
fprintf(fid,'*End Part\n');
fprintf(fid,'**\n');
fprintf(fid,'**\n');
fprintf(fid,'** ASSEMBLY\n');
fprintf(fid,'*Assembly, name=Assembly\n');
fprintf(fid,'**\n');
fprintf(fid,'*Instance, name=PART-1-1, part=PART-1\n');
fprintf(fid,'*End Instance\n');
fprintf(fid,'**\n');
fprintf(fid,'*Nset, nset=ALL, instance=PART-1-1, generate\n');
fprintf(fid,'1, %d, 1\n',size(Ver,1));
fprintf(fid,'*Elset, elset=SOLID_BODY, instance=PART-1-1, generate\n');
fprintf(fid,'1, %d, 1\n', size(TMesh,1));
fprintf(fid,'*End Assembly\n');


fclose(fid);