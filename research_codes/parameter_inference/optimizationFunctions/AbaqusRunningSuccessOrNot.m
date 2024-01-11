function SuccessB = AbaqusRunningSuccessOrNot(abaqusDir, abaqus_input_main_filename)
%%%this function needs to read Abaqus sta file and find out the last time
%%%step, if there is only one stimulation step, it should work. 
%%if fail, the last sentence is  THE ANALYSIS HAS NOT BEEN COMPLETED
%%if success, the last sentence is  THE ANALYSIS HAS COMPLETED SUCCESSFULLY

SuccessB = 0;SuccessB1 = 0;SuccessB2 = 1;

workingDir = pwd();
[~, filename, ~] = fileparts(abaqus_input_main_filename);
cd(abaqusDir);
fileName = sprintf('%s.sta',filename);
fid = fopen(fileName, 'r');


tline = fgetl(fid);
while ~feof(fid) 
       tline = fgetl(fid);
       matchStr = regexp(tline, 'SUCCESSFULLY', 'match');
       if ~isempty(matchStr)
           SuccessB1 = 1;
           break
       end 
end
fclose(fid);

%%%  modification from DB

fileName = sprintf('%s.dat',filename);
fid = fopen(fileName, 'r');


tline = fgetl(fid);
while ~feof(fid) 
       tline = fgetl(fid);
       matchStr = regexp(tline, 'NaN', 'match');
       if ~isempty(matchStr)
           SuccessB2 = 0;
           break
       end 
end
fclose(fid);

if SuccessB1==1 && SuccessB2==1
    SuccessB=1;
end

cd(workingDir);
    