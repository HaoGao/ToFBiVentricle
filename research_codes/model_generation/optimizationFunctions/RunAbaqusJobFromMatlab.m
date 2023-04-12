function [status, result] = RunAbaqusJobFromMatlab(abaqus_dir, abaqus_inputfile, ...
                      cpu_number, subroutine_file)

disp('abaqus running');

workingDir = pwd();
cd(abaqus_dir);


abaqus_inputfile = sprintf('%s.inp',abaqus_inputfile);

%%%now running the abaqus 
%%%first we delete all old results
delete('*.dat');
delete('*.com');
delete('*.lck');
delete('*.odb');
delete('*.par');
delete('*.msg');
delete('*.pes');
delete('*.pmg');
delete('*.prt');
delete('*.sim');
delete('*.sta');
delete('*.mdl');

command = sprintf('abaqus job=%s user=%s cpus=%d interactive',...
                  abaqus_inputfile, subroutine_file, cpu_number);

[status,result] = dos(command)
%% without ; so that it can post out the message in matlab


cd(workingDir);