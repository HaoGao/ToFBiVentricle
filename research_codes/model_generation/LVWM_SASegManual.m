%%%SA segmentation
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

%%now need to choose phase to segment
list_phase = {'early_diastole', 'end_diastole', 'end_systole'};
[idx, tf] = listdlg('ListString', list_phase);

totalSXSliceLocation = size(SXSliceSorted,2);
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
if ~exist('DataSegSA.mat', 'file')
    fresh_segB = 1;
end
fresh_segB=0;

if fresh_segB == 1
   %%%initialize 
    for imIndex = 1 : totalSXSliceLocation
        data.rect = [];
        data.endo_lv = [];
        data.endo_rv = [];
        data.epi_c = [];
        DataSegSA(imIndex)=data;
    end
    clear data;
    
    cd(phase_resultDir);
    save DataSegSA DataSegSA;
    cd(workingDir);
    
end


%%to keep it simple, will not sample the BC, but segment from beginning to
%%end, and output to endo and epi, respectively


if segB == 1
    %list_img = cell([totalSXSliceLocation,1]);
%     list_img = 1:totalSXSliceLocation;
    for i = 1: totalSXSliceLocation
        list_img{i} = sprintf('%d', i);
    end
    list_img{totalSXSliceLocation+1} = 'quit_loop';
    imIndex = 0; %% intital imIndex to be 1, and then will be updated

    quit_loop = 0; %%control whether continue to segment or not
    while ~quit_loop
       [imIndex, ~] = listdlg('ListString', list_img,'SelectionMode', 'single',...
           'InitialValue', imIndex+1);
       if strcmp( list_img{imIndex}, 'quit_loop' )
           quit_loop = 1;
           disp('quit the segment loop for all images');
       else
            if strcmp(list_phase{idx}, 'early_diastole')
                timeInstanceSelected =  SXSliceSorted(1,imIndex).TimeEarlyOfDiastole;
            elseif  strcmp(list_phase{idx}, 'end_diastole')
                timeInstanceSelected = SXSliceSorted(1,imIndex).TimeEndOfDiastole;
            elseif strcmp(list_phase{idx}, 'end_systole')
                timeInstanceSelected = SXSliceSorted(1,imIndex).TimeEndOfSystole;
            else
                disp('could not determine the cardiac phase, quit')
                return;
            end
            
            imData =  SXSliceSorted(1,imIndex).SXSlice(timeInstanceSelected).imData;
            imInfo1 = SXSliceSorted(1,imIndex).SXSlice(timeInstanceSelected).imInfo;
            imInfo = infoExtract(imInfo1);
            sliceLocationStr = sprintf('%s',imInfo.SliceLocation);
            %[imCropData, rect] = imcrop(imData);
            [imCropData, rect]= imcrop( uint8(double(imData)./(max(max(double(imData)))) *255) );

            hSA=figure();
            imshow(imCropData,[]);hold on;
            
            cd(phase_resultDir);
            load DataSegSA DataSegSA;
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
                but = 1;
                n=1;
                while but ==1 
                    [x, y, but]=myginput(1,'crosshair');
                    h2=plot(x,y,'b.');hold on;   
                    choice = questdlg('LV: Would you like to continue?', ...
                        'Dessert Menu', 'Yes','No','Delete','Yes');
                    % Handle response
                    switch choice
                        case 'Yes'
                            endo_lv(1,n)=x;
                            endo_lv(2,n)=y;
                            n=n+1;
                        case 'No'
                            endo_lv(1,n)=x;
                            endo_lv(2,n)=y;
                            n=n+1;
                            endo_lv(1,n)=endo_lv(1,1);
                            endo_lv(2,n)=endo_lv(2,1);
                            break
                        case 'Delete'
                            delete(h2)
                    end
                    
                end
                plot(endo_lv(1,:),endo_lv(2,:),'b');

                but = 1;
                n=1;
                while but ==1 
                    [x, y, but]=myginput(1,'crosshair');
                    h2=plot(x,y,'y.');hold on;   
                    choice = questdlg('RV: Would you like to continue?', ...
                        'Dessert Menu', 'Yes','No','Delete','Yes');
                    % Handle response
                    switch choice
                        case 'Yes'
                            endo_rv(1,n)=x;
                            endo_rv(2,n)=y;
                            n=n+1;
                        case 'No'
                            endo_rv(1,n)=x;
                            endo_rv(2,n)=y;
                            n=n+1;
                            endo_rv(1,n)=endo_rv(1,1);
                            endo_rv(2,n)=endo_rv(2,1);
                            break
                        case 'Delete'
                            delete(h2)
                    end
                    
                end
                plot(endo_rv(1,:),endo_rv(2,:),'y');
                
                
                n=1;but=1;
                while but==1
                    [x, y, but]=myginput(1, 'crosshair');
                    h2=plot(x,y,'r.');hold on;
                    choice = questdlg('Epi: Would you like to continue?', ...
                        'Dessert Menu', 'Yes','No','Delete','Yes');
                    % Handle response
                    switch choice
                        case 'Yes'
                            epi_c(1,n)=x;
                            epi_c(2,n)=y;
                            n=n+1;
                        case 'No'
                            epi_c(1,n)=x;
                            epi_c(2,n)=y;
                            n=n+1;
                            epi_c(1,n)=epi_c(1,1);
                            epi_c(2,n)=epi_c(2,1);
                            break
                        case 'Delete'
                            delete(h2)
                    end
                end
                plot(epi_c(1,:),epi_c(2,:),'r');
                
                endo_sample_lv(1,:) = endo_lv(1,:) +rect(1);
                endo_sample_lv(2,:) = endo_lv(2,:) +rect(2);
                endo_sample_rv(1,:) = endo_rv(1,:) +rect(1);
                endo_sample_rv(2,:) = endo_rv(2,:) +rect(2);
                epi_sample(1,:) = epi_c(1,:) +rect(1);
                epi_sample(2,:) = epi_c(2,:) +rect(2);
            
            
                h1=figure();
                imshow(imData,[]); hold on;
                title(sliceLocationStr);
                if ~isempty(endo_sample_lv)
                    plot(endo_sample_lv(1,:),endo_sample_lv(2,:),'b');
                end
                if ~isempty(endo_sample_rv)
                    plot(endo_sample_rv(1,:),endo_sample_rv(2,:),'b');
                end
                if ~isempty(epi_sample)
                     plot(epi_sample(1,:),epi_sample(2,:),'r');
                end
                
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

 
                DataSegSA(imIndex).endo_lv = endo_sample_lv;
                DataSegSA(imIndex).endo_rv = endo_sample_rv;
                DataSegSA(imIndex).epi_c = epi_sample;
                DataSegSA(imIndex).endo_lvReal = endo_lvReal;
                DataSegSA(imIndex).endo_rvReal = endo_rvReal;
                DataSegSA(imIndex).epi_cReal = epi_cReal;
               
                pause;
                close all;
                
                cd(phase_resultDir);
                save DataSegSA DataSegSA;
                save DataSegSAOri DataSegSA;
                cd(workingDir); 
    

            end %%segment
               
           
       end %% if strcmp 
    end  %% while
else
    disp('choose not to segment')
    cd(phase_resultDir);
    load DataSegSA DataSegSA;
    cd(workingDir);
end %% segB
    
 
%%%to show the boundaries in 3D with long axis views
% imFileName = LVOTSliceSorted(3).Time(timeInstanceSelected).name;
% imFileName = SXSliceSorted(4).Time(timeInstanceSelected).name;
% imFileName = sprintf('%s/%s',dicomDir,imFileName);
% imData = MRIMapToReal(imFileName);

imData = SXSliceSorted(1,5).SXSlice(timeInstanceSelected).imData;
imInfo1 = SXSliceSorted(1,5).SXSlice(timeInstanceSelected).imInfo;
imInfo = infoExtract(imInfo1);
imData = MRIMapToRealWithImageAndHeadData(imData, imInfo);


h3D = figure(); hold on;
DicomeRealDisplay(h3D, imData);

for imIndex = 1 : totalSXSliceLocation
% for imIndex = 3:3
    endo_lv = DataSegSA(imIndex).endo_lvReal;
    endo_rv = DataSegSA(imIndex).endo_rvReal;
    epi_c = DataSegSA(imIndex).epi_cReal;
    %%%%plot 3D curves
    if ~isempty(endo_lv) 
     plot3(endo_lv(1,:),endo_lv(2,:), endo_lv(3,:),'LineStyle', '-', 'Color', 'b', 'LineWidth',2);
    end
    if ~isempty(endo_rv) 
     plot3(endo_rv(1,:),endo_rv(2,:), endo_rv(3,:),'LineStyle', '-', 'Color', 'b', 'LineWidth',2);
    end
    if ~isempty(epi_c)
        plot3(epi_c(1,:),epi_c(2,:), epi_c(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
    end
end

imgTecplotOutput(imData,'SAslice3Tec.dat',phase_resultDir);
%%%this is for apex info










