%%%LVOT segmentation
%%% purpose: This will segment Aorta and MV, tricupsid valve and others,
%%% basically 4 valves will be segmented, thus in the configure file, needs
%%% to specify valve slices for those

clear all;
close all;
clc;

% LVWM_config;
LVWM_config;
fresh_segB = 1;
segB = 1;

cd(resultDir);
load imDesired;
cd(workingDir); 
%% need to reconsider to segment long-axis image again within LVOTSliceSorted
ValveSliceLocation = size(ValveSliceSorted,2);
LASliceLocation = size(LVOTSliceSorted,2);
% append LVSlices into ValveSlices
for i = 1 : LASliceLocation
    ValveSliceSorted(1,i+ValveSliceLocation).ValveSlice = LVOTSliceSorted(1,i).LXSlice;
    ValveSliceSorted(1,i+ValveSliceLocation).TimeEndOfSystole = LVOTSliceSorted(1,i).TimeEndOfSystole;
    ValveSliceSorted(1,i+ValveSliceLocation).TimeEarlyOfDiastole = LVOTSliceSorted(1,i).TimeEarlyOfDiastole;
    ValveSliceSorted(1,i+ValveSliceLocation).TimeEndOfDiastole = LVOTSliceSorted(1,i).TimeEndOfDiastole;
end

%%now need to choose phase to segment
list_phase = {'early_diastole', 'end_diastole', 'end_systole'};
[idx, tf] = listdlg('ListString', list_phase);

totalValveSliceLocation = size(ValveSliceSorted,2);
sampleN = patientConfigs.sampleN;


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

cd(phase_resultDir);
if ~exist('DataSegValve.mat', 'file')
    fresh_segB = 1;
end

if fresh_segB == 1
   %%%initialize 
    for imIndex = 1 : totalValveSliceLocation
        data.rect = [];
        data.endo_c = [];
        data.epi_c = [];
        DataSegValve(imIndex)=data;
    end
    clear data;
    
    cd(phase_resultDir);
    save DataSegValve DataSegValve;
    cd(workingDir);
    
end


%%to keep it simple, will not sample the BC, but segment from beginning to
%%end, and output to endo and epi, respectively


if segB == 1
    %list_img = cell([totalSXSliceLocation,1]);
%     list_img = 1:totalSXSliceLocation;
    for i = 1: totalValveSliceLocation
        list_img{i} = sprintf('%d', i);
    end
    list_img{totalValveSliceLocation+1} = 'quite_loop';
    imIndex = 0; %% intital imIndex to be 1, and then will be updated

    quite_loop = 0; %%control whether continue to segment or not
    while ~quite_loop
       [imIndex, ~] = listdlg('ListString', list_img,'SelectionMode', 'single',...
           'InitialValue', imIndex+1);
       if strcmp( list_img{imIndex}, 'quite_loop' )
           quite_loop = 1;
           disp('quite the segment loop for all images');
       else
            if strcmp(list_phase{idx}, 'early_diastole')
                timeInstanceSelected =  ValveSliceSorted(1,imIndex).TimeEarlyOfDiastole;
            elseif  strcmp(list_phase{idx}, 'end_diastole')
                timeInstanceSelected = ValveSliceSorted(1,imIndex).TimeEndOfDiastole;
            elseif strcmp(list_phase{idx}, 'end_systole')
                timeInstanceSelected = ValveSliceSorted(1,imIndex).TimeEndOfSystole;
            else
                disp('could not determine the cardiac phase, quite')
                return;
            end
            
            imData =  ValveSliceSorted(1,imIndex).ValveSlice(timeInstanceSelected).imData;
            imInfo1 = ValveSliceSorted(1,imIndex).ValveSlice(timeInstanceSelected).imInfo;
            imInfo = infoExtract(imInfo1);
            sliceLocationStr = sprintf('%s',imInfo.SliceLocation);
            [imCropData, rect] = imcrop(imData,[]);

            hSA=figure();
            imshow(imCropData,[]);hold on;
            
            cd(phase_resultDir);
            load DataSegValve DataSegValve;
            cd(workingDir);
            
            %%need to ask whether to segment or not
            answer = questdlg('would you like to segment','segment or skip');
            if strcmp(answer, 'Yes') %%then segment 
                endo_c=[];
                epi_c=[];
                endo_sample = [];
                epi_sample = [];
                but = 1;
                n=1;
                while but ==1 
                    [x, y, but]=myginput(1,'crosshair');
                    plot(x,y,'b.');hold on;
                    endo_c(1,n)=x;
                    endo_c(2,n)=y;
                    n=n+1;
                end
                plot(endo_c(1,:),endo_c(2,:),'b');
                
                
                n=1;but=1;
                while but==1
                    [x, y, but]=myginput(1, 'crosshair');
                    plot(x,y,'r.');hold on;
                    epi_c(1,n)=x;
                    epi_c(2,n)=y;
                    n=n+1;
                end
                
                endo_sample(1,:) = endo_c(1,:) +rect(1);
                endo_sample(2,:) = endo_c(2,:) +rect(2);
                epi_sample(1,:) = epi_c(1,:) +rect(1);
                epi_sample(2,:) = epi_c(2,:) +rect(2);
            
            
                h1=figure();
                imshow(imData,[]); hold on;
                title(sliceLocationStr);
                if ~isempty(endo_sample)
                    plot(endo_sample(1,:),endo_sample(2,:),'b');
                end
                if ~isempty(epi_sample)
                     plot(epi_sample(1,:),epi_sample(2,:),'r');
                end
                
                if ~isempty(endo_sample)
                    endo_cReal = TransformCurvesFromImToRealSpace(endo_sample,imInfo);
                else
                    endo_cReal = [];
                end
                if ~isempty(epi_sample)
                    epi_cReal = TransformCurvesFromImToRealSpace(epi_sample,imInfo);
                else
                    epi_cReal = [];
                end

 
                DataSegValve(imIndex).endo_c = endo_sample;
                DataSegValve(imIndex).epi_c = epi_sample;
                DataSegValve(imIndex).endo_cReal = endo_cReal;
                DataSegValve(imIndex).epi_cReal = epi_cReal;
               
                pause;
                close all;
                
                cd(phase_resultDir);
                save DataSegValve DataSegValve;
                save DataSegValveOri DataSegValve;
                cd(workingDir); 
    

            end %%segment
               
           
       end %% if strcmp 
    end  %% while
else
    disp('choose not to segment')
    cd(phase_resultDir);
    load DataSegValve DataSegValve;
    cd(workingDir);
end %% segB
    
 
%%%to show the boundaries in 3D with long axis views
% imFileName = LVOTSliceSorted(3).Time(timeInstanceSelected).name;
% imFileName = SXSliceSorted(4).Time(timeInstanceSelected).name;
% imFileName = sprintf('%s/%s',dicomDir,imFileName);
% imData = MRIMapToReal(imFileName);

imData = ValveSliceSorted(1,1).ValveSlice(timeInstanceSelected).imData;
imInfo1 = ValveSliceSorted(1,1).ValveSlice(timeInstanceSelected).imInfo;
imInfo = infoExtract(imInfo1);
imData = MRIMapToRealWithImageAndHeadData(imData, imInfo);


h3D = figure(); hold on;
DicomeRealDisplay(h3D, imData);

for imIndex = 1 : totalValveSliceLocation
% for imIndex = 3:3
    endo_c = DataSegValve(imIndex).endo_cReal;
    epi_c = DataSegValve(imIndex).epi_cReal;
    %%%%plot 3D curves
    if ~isempty(endo_c) 
     plot3(endo_c(1,:),endo_c(2,:), endo_c(3,:),'LineStyle', '-', 'Color', 'b', 'LineWidth',2);
    end
    if ~isempty(epi_c)
        plot3(epi_c(1,:),epi_c(2,:), epi_c(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
    end
end


%%%this is for apex info










