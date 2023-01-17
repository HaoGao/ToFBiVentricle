%%it will try to figure out alpha0 for each case, and saved in the
%%result_root dir, then it will be avaialbe for generating meshes using
%%same alpha0
clear all; close all; clc;

%%%find the proper directory 
LVWM_config;

%%%% load parameter for one study
cd(resultDir);
if exist('SetParameters.m', 'file')
    SetParameters;
    cd(workingDir);
else
    cd(workingDir);
    SetParametersDefault;
end


%%%only for cubic spline
if(kS ~=4)
    disp('The program is ONLY for Cubic B-Spline with kS=4');
    disp('Please change kS in setparameters back to kS=4');
    stop;
end

%%%pre_processing, from segmentated boundary data
readGuidePoints_preProcessing_GlobalAlpha0(outterGuidePointsFileName, innerGuidePointsFileName,...
                outterGuide4FittingFileName,innerGuide4FittingFileName, ...
                prolateParametersFileName, resultDir, ...
                resultDirRoot, patientConfigs.name);
            
cd(workingDir);
            
      