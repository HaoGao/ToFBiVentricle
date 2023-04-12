% LVWM_measuring mitral inflow velocity

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

TimeInstance = 3;

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
if ~exist('DataSegLA_MVInflow.mat', 'file')
    fresh_segB = 1;
end

if fresh_segB == 1
   %%%initialize 
    for imIndex = 1 : TimeInstance
        data.rect = [];
        data.endo_c = [];
        data.epi_c = [];
        data.img = [];
        DataSegLA_MVInflow(imIndex)=data;
    end
    clear data;
    
    cd(phase_resultDir);
    save DataSegLA_MVInflow DataSegLA_MVInflow;
    cd(workingDir);
    
end



if segB == 1
    %list_img = cell([totalSXSliceLocation,1]);
%     list_img = 1:totalSXSliceLocation;
    for i = 1: TimeInstance
        list_img{i} = sprintf('%d', i);
    end
    list_img{TimeInstance+1} = 'quite_loop';
    imIndex = 0;
    
    quite_loop = 0; %%control whether continue to segment or not
    while ~quite_loop
       [imIndex, ~] = listdlg('ListString', list_img,'SelectionMode', 'single',...
           'InitialValue', imIndex+1);
       if strcmp( list_img{imIndex}, 'quite_loop' )
           quite_loop = 1;
           disp('quite the segment loop for all images');
       else
            if strcmp(list_phase{idx}, 'early_diastole')
                timeInstanceSelected =  LVOTSliceSorted(1,2).TimeEarlyOfDiastole;
            elseif  strcmp(list_phase{idx}, 'end_diastole')
                timeInstanceSelected = LVOTSliceSorted(1,2).TimeEndOfDiastole;
            elseif strcmp(list_phase{idx}, 'end_systole')
                timeInstanceSelected = LVOTSliceSorted(1,2).TimeEndOfSystole;
            else
                disp('could not determine the cardiac phase, quite')
                return;
            end
            
            imData =  LVOTSliceSorted(1,2).LXSlice(timeInstanceSelected+imIndex-1).imData;
            imInfo1 = LVOTSliceSorted(1,2).LXSlice(timeInstanceSelected+imIndex-1).imInfo;
            imInfo = infoExtract(imInfo1);
            sliceLocationStr = sprintf('%s',imInfo.SliceLocation);
            [imCropData, rect] = imcrop(imData,[]);

            hLA=figure();
            imshow(imCropData,[]);hold on;
            %%need to show the previous boundary
            if imIndex > 1
                endo_pre = DataSegLA_MVInflow(imIndex-1).endo_c;
                endo_pre(1,:) = endo_pre(1,:) - rect(1);
                endo_pre(2,:) = endo_pre(2,:) - rect(2);
                figure(hLA); hold on;
                plot(endo_pre(1,:), endo_pre(2,:), 'r-.');
            end
            
            
            cd(phase_resultDir);
            load DataSegLA_MVInflow;
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
                
                endo_sample(1,:) = endo_c(1,:) +rect(1);
                endo_sample(2,:) = endo_c(2,:) +rect(2);
                
            
                h1=figure();
                imshow(imData,[]); hold on;
                title(sliceLocationStr);
                if ~isempty(endo_sample)
                    plot(endo_sample(1,:),endo_sample(2,:),'b');
                end
                if ~isempty(epi_sample)
                     plot(epi_sample(1,:),epi_sample(2,:),'r');
                end

                DataSegLA_MVInflow(imIndex).endo_c = endo_sample;
                DataSegLA_MVInflow(imIndex).img = imData;
                DataSegLA_MVInflow(imIndex).imInfo = imInfo1;
                
                BW = roipoly(imData,endo_sample(1,:), endo_sample(2,:));
                roi_area = sum(sum(BW));
                ori_dis_vec = endo_sample(:,1) - endo_sample(:,end);
                ori_dis = (sum(ori_dis_vec.^2) )^0.5;
                
                DataSegLA_MVInflow(imIndex).roi_area = roi_area;
                DataSegLA_MVInflow(imIndex).ori_dis = ori_dis;
                DataSegLA_MVInflow(imIndex).TriggerTime = imInfo1.TriggerTime;
               
                pause;
                close all;
                
                cd(phase_resultDir);
                save DataSegLA_MVInflow DataSegLA_MVInflow;
                cd(workingDir); 
    

            end %%segment
               
           
       end %% if strcmp 
    end  %% while
else
    disp('choose not to segment')
    cd(phase_resultDir);
    load DataSegLA_MVInflow DataSegLA_MVInflow;
    cd(workingDir);
end %% segB


for imIndex = 1 : TimeInstance-1
    roi_area = DataSegLA_MVInflow(imIndex).roi_area;
    roi_area_n = DataSegLA_MVInflow(imIndex+1).roi_area;
    ori_dis = DataSegLA_MVInflow(imIndex).ori_dis;
    
    dt = DataSegLA_MVInflow(imIndex+1).TriggerTime - DataSegLA_MVInflow(imIndex).TriggerTime;
    dt = dt*0.001;
    dx = DataSegLA_MVInflow(imIndex).imInfo.PixelSpacing;
    dx = dx.*0.1;
    
    vel_mv = (roi_area_n - roi_area)/ori_dis/dt*dx(1)
    
    annu_vec_1 = DataSegLA_MVInflow(imIndex).endo_c(:,1) - DataSegLA_MVInflow(imIndex+1).endo_c(:,1);
    annu_dis_1 = (sum(annu_vec_1.^2))^0.5;
    
    annu_vec_2 = DataSegLA_MVInflow(imIndex).endo_c(:,end) - DataSegLA_MVInflow(imIndex+1).endo_c(:,end);
    annu_dis_2 = (sum(annu_vec_2.^2))^0.5;
    
    vel_annu = (annu_dis_2+annu_dis_1)/2*dx(1)/dt
    
end

%now need to track the annulus ring velocity 



