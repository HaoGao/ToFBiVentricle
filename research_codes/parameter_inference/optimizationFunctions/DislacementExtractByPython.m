function DislacementExtractByPython(abaqusSimulationDir,abaqus_inputfile,...
    pythonOriginalFilesDir,pythonfile, lineNoForODBName, lineNoForDisName,...
    abaqus_dis_out_filename, abaqus_command)

workingDir = pwd();

%%%abaqus odb file name
[~, abaqusInputfileName, ~]= fileparts(abaqus_inputfile);
abaqus_odb_name = sprintf('%s.odb',abaqusInputfileName);

%%%need to change some lines in pythonfile
cd(pythonOriginalFilesDir);
fid_python = fopen(pythonfile,'r');
cd(abaqusSimulationDir);
fid_python_updated = fopen('abaqus_dis_up.py','w');
cd(workingDir);

lineIndex = 1;
tline = fgetl(fid_python);
fprintf(fid_python_updated,'%s\r\n',tline);

while ~feof(fid_python)
    lineIndex = lineIndex + 1;
    tline = fgetl(fid_python);
    if lineIndex == lineNoForODBName
        tline = sprintf('myOdb = openOdb(path=''%s'')',abaqus_odb_name);
    elseif lineIndex == lineNoForDisName
        tline = sprintf('outfilename=''%s''',abaqus_dis_out_filename );
    end
    fprintf(fid_python_updated,'%s\r\n',tline);
end
fclose(fid_python_updated);
fclose(fid_python);
%%%now run the python file to get the displacement
cd(abaqusSimulationDir);
command = sprintf('%s python abaqus_dis_up.py -odb, %s',abaqus_command, abaqus_odb_name);
dos(command);


% eval(command);




cd(workingDir);