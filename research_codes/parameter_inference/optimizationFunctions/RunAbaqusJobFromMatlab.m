function [status, result] = RunAbaqusJobFromMatlab(abaqus_command, abaqus_dir, abaqus_inputfile, ...
                      cpu_number, subroutine_file)

disp('abaqus running');

workingDir = pwd();
cd(abaqus_dir);


%abaqus_inputfile = sprintf('%s.inp',abaqus_inputfile);

%%%now running the abaqus 
%%%first we delete all old results
if ~isempty(dir('*.dat'))
    delete('*.dat');
end

if ~isempty(dir('*.com'))
    delete('*.com');
end

if ~isempty(dir('*.lck'))
    delete('*.lck');
end

if ~isempty(dir('*.odb'))
    delete('*.odb');
end

if ~isempty(dir('*.par'))
    delete('*.par');
end

if ~isempty(dir('*.msg'))
    delete('*.msg');
end

if ~isempty(dir('*.pes'))
    delete('*.pes');
end

if ~isempty(dir('*.pmg'))
    delete('*.pmg');
end

if ~isempty(dir('*.prt'))
    delete('*.prt');
end

if ~isempty(dir('*.sim'))
    delete('*.sim');
end

if ~isempty(dir('*.sta'))
    delete('*.sta');
end

if ~isempty(dir('*.mdl'))
    delete('*.mdl');
end

command = sprintf('%s job=%s user=%s cpus=%d interactive',...
                  abaqus_command, abaqus_inputfile, subroutine_file, cpu_number);

[status,result] = dos(command)
%% without ; so that it can post out the message in matlab


cd(workingDir);