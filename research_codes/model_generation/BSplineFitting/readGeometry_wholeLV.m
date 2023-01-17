function [alpha0, w1, w2, umax, umin] = readGeometry_wholeLV(resultDir, prolateParametersFileName)

% fid = fopen('ProlateParameters.dat','r');
workingDir = pwd();
cd(resultDir);
fid = fopen(prolateParametersFileName,'r');
cd(workingDir);

tline = fgetl(fid);
alpha0 = sscanf(tline, '%f');

tline = fgetl(fid);
w1 = sscanf(tline, '%f');

tline = fgetl(fid);
w2 = sscanf(tline, '%f');

tline = fgetl(fid);
umax = sscanf(tline, '%f');

tline = fgetl(fid);
umin = sscanf(tline, '%f');
 
% w0 = w1;
umax = umax*pi/180;
umin = umin*pi/180;

sline = sprintf('alpha0: \tw1:\t w2: \tumax: \t umin:');
disp(sline);

sline = sprintf('%f \t%f   \t%f  \t%f \t%f', alpha0, w1, w2, umax, umin);
disp(sline);


fclose(fid);