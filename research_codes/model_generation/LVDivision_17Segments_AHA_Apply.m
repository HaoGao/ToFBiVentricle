%%%this is used to divide LV into different regions
%%%warning there must be six short axis slice availabe for post-processing,
%%%if it is not, then it might be problem

clear all;  close all; clc;

% LVWM_config;
LVWM_config;

%%now need to choose phase to segment
% cd(resultDir);
% cd(deformetricaDir);
% deformetricaDir = pwd();
% cd(workingDir);

% cd(resultDir);
% cd(deformetrica_libMeshDir);
% deformetrica_libMeshDir = pwd();
% cd(workingDir);

%%now need to choose phase to segment
list_phase = {'early_diastole', 'end_diastole', 'end_systole'};
[idx, tf] = listdlg('ListString', list_phase);
cd(resultDir);
if ~exist(list_phase{idx},'dir')
    mkdir(list_phase{idx});
    cd(list_phase{idx});
    phase_resultDir = pwd();    
else
    cd(list_phase{idx});
    phase_resultDir = pwd();
end
cd(workingDir);
phase_selected = list_phase{idx};

cd(phase_resultDir);
cd(solidworksDir);
solidworksDir = pwd();
cd(gmeshDir);
gmeshDir = pwd();
cd(workingDir);

cd(gmeshDir); 
load abaqusInput;
load LV_RV_assignment;
cd(workingDir);

cd(resultDir);
load DivisionConfig;
load imDesired; %%load in the images to make sure there is enough data available
cd(workingDir);


%%%now need to extract the right data 
totalSXSliceLocation = size(SXSliceSorted, 2);
if totalSXSliceLocation-patientConfigs(1,1).basalSliceIndex+1 < 6
    errordlg('<6 slices are loaded','Image Number Error');
    pause;
end

%%%this part needs to be checked every time to ensure things are right

SliceThickness = SXSliceSorted(patientConfigs(1,1).basalSliceIndex).SXSlice(1).imInfo.SliceThickness;
Slice_location1 = SXSliceSorted(patientConfigs(1,1).basalSliceIndex).SXSlice(1).imInfo.SliceLocation;
Slice_location2 = SXSliceSorted(patientConfigs(1,1).basalSliceIndex+1).SXSlice(1).imInfo.SliceLocation;
SASliceDistance = round(abs(Slice_location2-Slice_location1));
SliceGap = SASliceDistance - SliceThickness;

%%%the apex region will be > 
usuableSXSlice = totalSXSliceLocation - patientConfigs(1,1).basalSliceIndex+1;
SASliceBase = patientConfigs(1,1).basalSliceIndex;
SASlicePositionApex = patientConfigs(1,1).SASlicePositionApex;
%%%%%%%%%%%%


 LVMeshDivisionAHA(resultDir,gmeshDir, phase_resultDir, MidConfig, ApexConfig, SASliceBase, SASlicePositionApex, ...
                   SASliceDistance,SliceThickness,usuableSXSlice, abaqusInput, LV_RV_assignment);

% LVMeshDivision_17Segments_AHA(resultDir, gmeshDir, ...
%                               MidConfig, ApexConfig, SASlicePositionBase, SASlicePositionMiddle, SASlicePositionApex, ...
%                               SASliceDistance,usuableSXSlice, abaqusInput);

%%mid slice
%% 1: inferior septum
%% 2: anterior septum
%% 3: anterior
%% 4: anterior laterial
%% 5: inferior laterial
%% 6: inferior

%%apical slice
%% 1: septum
%% 2: anterior
%% 3: laterial
%% 4: inferior

%%apex: 7

