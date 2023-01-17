function SuccessB = AbaqusRunningSuccessOrNot(abaqusDir, abaqus_input_main_filename)
%%%this function needs to read Abaqus sta file and find out the last time
%%%step, if there is only one stimulation step, it should work. 
%%if fail, the last sentence is  THE ANALYSIS HAS NOT BEEN COMPLETED
%%if success, the last sentence is  THE ANALYSIS HAS COMPLETED SUCCESSFULLY

SuccessB = 0;

workingDir = pwd();
cd(abaqusDir);
fileName = sprintf('%s.sta',abaqus_input_main_filename);
fid = fopen(fileName, 'r');


tline = fgetl(fid);
while ~feof(fid) 
       tline = fgetl(fid);
       matchStr = regexp(tline, 'SUCCESSFULLY', 'match');
       if ~isempty(matchStr)
           SuccessB = 1;
       end 
end
fclose(fid);

cd(workingDir);
    