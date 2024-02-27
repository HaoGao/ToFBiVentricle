%%%SA segmentation
clear all;
close all;
clc;

% LVWM_config;
LVWM_config;
fresh_segB = 0;
segB = 1;

cd(resultDir);
load imDesired;
cd(workingDir); 

%%now need to choose phase to segment
list_phase = {'early_diastole', 'end_diastole', 'end_systole'};
[idx, tf] = listdlg('ListString', list_phase);

totalLASliceLocation = size(LVOTSliceSorted,2);
sampleN = patientConfigs.sampleN;

cd(resultDir);
if ~exist(list_phase{idx},'dir')
    mkdir(list_phase{idx});
    cd(list_phase{idx});
    phase_resultDir = pwd();     
end
cd(list_phase{idx});
phase_resultDir = pwd();
cd(workingDir);

cd(phase_resultDir);
if ~exist('DataSegLA.mat', 'file')
    fresh_segB = 1;
end

if fresh_segB == 1
   %%%initialize 
    for imIndex = 1 : totalLASliceLocation
        data.rect = [];
        data.endo_lv = [];
        data.endo_rv = [];
        data.epi_c = [];
        DataSegLA(imIndex)=data;
    end
    clear data;
    
    cd(phase_resultDir);
    save DataSegLA DataSegLA;
    cd(workingDir);
    
end


%%to keep it simple, will not sample the BC, but segment from beginning to
%%end, and output to endo and epi, respectively


if segB == 1
    %list_img = cell([totalSXSliceLocation,1]);
%     list_img = 1:totalSXSliceLocation;
    for i = 1: totalLASliceLocation
        list_img{i} = sprintf('%d', i);
    end
    list_img{totalLASliceLocation+1} = 'quit_loop';
    imIndex = 0;
    
    quit_loop = 0; %%control whether continue to segment or not
    while ~quit_loop
       [imIndex, ~] = listdlg('ListString', list_img,'SelectionMode', 'single',...
           'InitialValue', imIndex+1);
       if strcmp( list_img{imIndex}, 'quit_loop' )
           quit_loop = 1;
           disp('quit the segment loop for all images');
       else
            if strcmp(list_phase{idx}, 'early_diastole')
                timeInstanceSelected =  LVOTSliceSorted(1,imIndex).TimeEarlyOfDiastole;
            elseif  strcmp(list_phase{idx}, 'end_diastole')
                timeInstanceSelected = LVOTSliceSorted(1,imIndex).TimeEndOfDiastole;
            elseif strcmp(list_phase{idx}, 'end_systole')
                timeInstanceSelected = LVOTSliceSorted(1,imIndex).TimeEndOfSystole;
            else
                disp('could not determine the cardiac phase, quit')
                return;
            end
            
            imData =  LVOTSliceSorted(1,imIndex).LXSlice(timeInstanceSelected).imData;
            imInfo1 = LVOTSliceSorted(1,imIndex).LXSlice(timeInstanceSelected).imInfo;
            imInfo = infoExtract(imInfo1);
            sliceLocationStr = sprintf('%s',imInfo.SliceLocation);
            %[imCropData, rect] = imcrop(imData,[]);
            [imCropData, rect]= imcrop( uint8(double(imData)./(max(max(double(imData)))) *255) );

            hLA=figure();
            imshow(imCropData,[]);hold on;
            
            cd(phase_resultDir);
            load DataSegLA DataSegLA;
            cd(workingDir);
            
            %%need to ask whether to segment or not
            answer = questdlg('would you like to segment','segment or skip');
            if strcmp(answer, 'Yes') %%then segment 
                endo_lv=[];
                endo_rv=[];
                epi_c=[];
                endo_sample_lv = [];
                endo_sample_rv = [];
                epi_sample = [];
                
                disp('.............................................');
                disp('predefined boundaries are plotted in red');
                disp('define points in the boundary by click points')
                disp('press stop, then click the last point');
                disp('press replot to generate the curve');
                disp('whenever change the point position, using replot to regenerate the curve');
                disp('double press space key to return from boudnary definition');
                 
                choice = questdlg('LV endo: Would you like to continue?', ...
                        'Dessert Menu', 'Yes','No','Yes');  
                if strcmp(choice, 'Yes')
                    [endo_lv, hpts_endo_lv]=defineBoundaryByimpoint_2023(imCropData, [], 0, 'LV endo');
                    endo_sample_lv(1,:) = endo_lv(1,:) +rect(1);
                    endo_sample_lv(2,:) = endo_lv(2,:) +rect(2);
                end

                choice = questdlg('RV endo: Would you like to continue?', ...
                    'Dessert Menu', 'Yes','No','Yes');
                if strcmp(choice, 'Yes')
                    [endo_rv, hpts_endo_rv]=defineBoundaryByimpoint_2023(imCropData, endo_lv, 0, 'RV endo');
                    endo_sample_rv(1,:) = endo_rv(1,:) +rect(1);
                    endo_sample_rv(2,:) = endo_rv(2,:) +rect(2);
                end

                choice = questdlg('EPI: Would you like to continue?', ...
                    'Dessert Menu', 'Yes','No','Yes');
                if strcmp(choice, 'Yes')
                    [epi_c, hpts_epi_c]=defineBoundaryByimpoint_2023(imCropData, [endo_lv, endo_rv], 0, 'epi bc');
                    epi_sample(1,:) = epi_c(1,:) +rect(1);
                    epi_sample(2,:) = epi_c(2,:) +rect(2);
                end
            
            
                h1=figure();
                imshow(imCropData,[]); hold on;
                title(sliceLocationStr);
                if ~isempty(endo_sample_lv)
                    bc_interp = samplingBCWithoutIm_interp(endo_sample_lv, 0);
                    plot(bc_interp(1,:)-rect(1),bc_interp(2,:)-rect(2),'b-');
                end
                if ~isempty(endo_sample_rv)
                    bc_interp = samplingBCWithoutIm_interp(endo_sample_rv, 0);
                    plot(bc_interp(1,:)-rect(1),bc_interp(2,:)-rect(2), 'b-');
                end
                if ~isempty(epi_sample)
                    bc_interp = samplingBCWithoutIm_interp(epi_sample, 0);
                     plot(bc_interp(1,:)-rect(1),bc_interp(2,:)-rect(2), 'r-');
                end
                
                %save the figure with boundaries
                cd(phase_resultDir);
                figure_name = sprintf('LA-%d.png',imIndex);
                figure(h1);ax=gca;
                exportgraphics(ax,figure_name,'Resolution',300) 
                cd(workingDir);
                
                if ~isempty(endo_sample_lv)
                    endo_lvReal = TransformCurvesFromImToRealSpace(endo_sample_lv,imInfo);
                else
                    endo_lvReal = [];
                end
                if ~isempty(endo_sample_rv)
                    endo_rvReal = TransformCurvesFromImToRealSpace(endo_sample_rv,imInfo);
                else
                    endo_rvReal = [];
                end
                if ~isempty(epi_sample)
                    epi_cReal = TransformCurvesFromImToRealSpace(epi_sample,imInfo);
                else
                    epi_cReal = [];
                end

 
                DataSegLA(imIndex).endo_lv = endo_sample_lv;
                DataSegLA(imIndex).endo_rv = endo_sample_rv;
                DataSegLA(imIndex).epi_c = epi_sample;
                DataSegLA(imIndex).endo_lvReal = endo_lvReal;
                DataSegLA(imIndex).endo_rvReal = endo_rvReal;
                DataSegLA(imIndex).epi_cReal = epi_cReal;
               
                pause;
                close all;
                
                cd(phase_resultDir);
                save DataSegLA DataSegLA;
                save DataSegLAOri DataSegLA;
                cd(workingDir); 
    

            end %%segment
               
           
       end %% if strcmp 
    end  %% while
else
    disp('choose not to segment')
    cd(phase_resultDir);
    load DataSegLA DataSegLA;
    cd(workingDir);
end %% segB
    
 
%%%to show the boundaries in 3D with long axis views
% imFileName = LVOTSliceSorted(3).Time(timeInstanceSelected).name;
% imFileName = SXSliceSorted(4).Time(timeInstanceSelected).name;
% imFileName = sprintf('%s/%s',dicomDir,imFileName);
% imData = MRIMapToReal(imFileName);

sliceIndex = patientConfigs(1,1).basalSliceIndex+2; % the middle ventrilce
imData = SXSliceSorted(1,sliceIndex).SXSlice(timeInstanceSelected).imData;
imInfo1 = SXSliceSorted(1,sliceIndex).SXSlice(timeInstanceSelected).imInfo;
imInfo = infoExtract(imInfo1);
imData = MRIMapToRealWithImageAndHeadData(imData, imInfo);


h3D = figure(); hold on;
DicomeRealDisplay(h3D, imData);

for imIndex = 1 : totalLASliceLocation
% for imIndex = 3:3
    endo_lv = DataSegLA(imIndex).endo_lvReal;
    endo_rv = DataSegLA(imIndex).endo_rvReal;
    epi_c  = DataSegLA(imIndex).epi_cReal;
    %%%%plot 3D curves
    if ~isempty(endo_lv) 
     plot3(endo_lv(1,:),endo_lv(2,:), endo_lv(3,:),'LineStyle', '-', 'Color', 'b', 'LineWidth',2);
    end
    if ~isempty(endo_rv) 
     plot3(endo_rv(1,:),endo_rv(2,:), endo_rv(3,:),'LineStyle', '-', 'Color', 'y', 'LineWidth',2);
    end
    if ~isempty(epi_c)
        plot3(epi_c(1,:),epi_c(2,:), epi_c(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
    end
end

imgTecplotOutput(imData,'LAslice3Tec.dat',phase_resultDir);
%%%this is for apex info










