
%%%whole lv manual segmentation configuration
if ispc
    %path(path, '..\model_generation\segmentation');
    %path(path, '.\BSplineFitting');
    %path(path, '.\meshRelated');
    %path(path, '..\model_generation\FEMfunctions');
    path(path, '.\optimizationFunctions');
    path(path, '.\MetricsMeshes');
end

if ismac
    %path(path, './segmentation');
    %path(path, './BSplineFitting');
    path(path, './meshRelated');
    path(path, './FEMfunctions');
    path(path, './optimizationFunctions');
end

workingDir = pwd();



%%load the patient config file
%%resultDirRoot = 'D:\HaoGao\PhDs\PhD_Alan\Results';
resultDirRoot = 'D:\BHFToF\workspace\Results';
cd(resultDirRoot);

if ispc 
    [FileName,PathName,~] = uigetfile('..\..\..\Results\*.m');
elseif ismac || isunix 
    [FileName,PathName,~] = uigetfile('../../../Results/*.m');
end
% [FileName, PathName] = uigetfile( ...
%        {'*.m'}, ...
%         'Pick a file');
projectConfig_dir = PathName;
projectConfig_name = FileName(1:end-2);
cd(workingDir);

cd(projectConfig_dir);
run(projectConfig_name);
cd(workingDir);



