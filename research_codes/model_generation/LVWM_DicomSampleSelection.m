%%%%select the proper dicom images for short axis and long axis
clear all; 
close all;
clc;

LVWM_config;

totalPatientNo = size(patientConfigs,1); %%%currently only cosidering one patient in 
if totalPatientNo > 1
    errordlg('load only one patient only');
end

SASlicePosition = 0;
LVOTSlicePosition = 0;
FourChSlicePosition = 0;
OneChSlicePosition = 0;
LGESASlicePosition = 0;
ValveSlicePosition = 0;
RVOTSlicePosition = 0;

for patientIndex = 1 : 1 
       patientT = patientConfigs(patientIndex,1);
       patientName = patientT.name;
    
       totalStudyNo = size(patientT.studyName,1);
       for studyIndex = 1 : totalStudyNo
          studyName_sel = patientConfigs(patientIndex,1).studyName(studyIndex,1).studyName;
          studyDir_sel = patientConfigs(patientIndex,1).dir(studyIndex,1).studyDir;
          ImgDir_sel = patientConfigs(patientIndex,1).dirMidSA(studyIndex,1).ImgDir;
          spec_sel = patientConfigs(patientIndex,1).SliceSpec(studyIndex,1).spec;
          
          %%%now decide it is SA or LA, or LAG   
          if strcmp(spec_sel, 'SAcine')
             %%%load the image into SASlice 
              SASlicePosition = SASlicePosition + 1;
              
              %if ismac
                  dicomDir = sprintf('%s\\%s',studyDir_sel, ImgDir_sel);
              %end
              SXSliceT = loadAllDicomSAOrLAFromOneDir(dicomDir);
              SXSlice(SASlicePosition).SXSlice= SXSliceT;
              imDesired.SXSlice = SXSlice;
              imDesired.TimeEndOfSystole = patientT.TimeEndOfSystole;
              imDesired.TimeEndOfDiastole = patientT.TimeEndOfDiastole;
              imDesired.TimeEarlyOfDiastole = patientT.TimeEarlyOfDiastole;
          elseif  strcmp(spec_sel, 'LAcine_LVOT')
              LVOTSlicePosition = LVOTSlicePosition + 1;
              dicomDir = sprintf('%s\\%s',studyDir_sel, ImgDir_sel);
              %if ismac
                  dicomDir = sprintf('%s\\%s',studyDir_sel, ImgDir_sel);
              %end
              LVOTSliceT = loadAllDicomSAOrLAFromOneDir(dicomDir);
              LVOTSlice(LVOTSlicePosition).LVOTSlice= LVOTSliceT;
              imDesired.LVOTSlice = LVOTSlice;
              
          elseif  strcmp(spec_sel, 'LAcine_4CH')
              FourChSlicePosition = FourChSlicePosition + 1;
              dicomDir = sprintf('%s\\%s',studyDir_sel, ImgDir_sel);
              %if ismac
                  dicomDir = sprintf('%s\\%s',studyDir_sel, ImgDir_sel);
              %end
              FourCHSliceT = loadAllDicomSAOrLAFromOneDir(dicomDir);
              FourCHSlice(LVOTSlicePosition).FourCHSlice= FourCHSliceT;
              imDesired.FourCHSlice = FourCHSlice;
              
         elseif  strcmp(spec_sel, 'LAcine_1CH')
              OneChSlicePosition = OneChSlicePosition + 1;
              dicomDir = sprintf('%s\\%s',studyDir_sel, ImgDir_sel);
              %if ismac
                  dicomDir = sprintf('%s\\%s',studyDir_sel, ImgDir_sel);
              %end
              
              OneCHSliceT = loadAllDicomSAOrLAFromOneDir(dicomDir);
              OneCHSlice(OneChSlicePosition).OneCHSlice= OneCHSliceT;
              imDesired.OneCHSlice = OneCHSlice;
         
        %%%adding SA series of LGE images
          elseif strcmp(spec_sel, 'LGE_SA')
              LGESASlicePosition = LGESASlicePosition + 1;
              dicomDir = sprintf('%s\\%s',studyDir_sel, ImgDir_sel);
              %if ismac
                  dicomDir = sprintf('%s\\%s',studyDir_sel, ImgDir_sel);
              %end
              
              LGESASliceT = loadAllDicomLGEImage(dicomDir);
%               LGESASliceT(LGESASlicePosition).LGESASlice= LGESASliceT;
              imDesired.LGESASlice = LGESASliceT;
              
              
          elseif strcmp(spec_sel, 'valve')    
              %%%load the image into SASlice 
              ValveSlicePosition = ValveSlicePosition + 1;
              dicomDir = sprintf('%s\\%s',studyDir_sel, ImgDir_sel);
              %if ismac
                  dicomDir = sprintf('%s\\%s',studyDir_sel, ImgDir_sel);
              %end
              ValveSliceT = loadAllDicomSAOrLAFromOneDir(dicomDir);
              ValveSlice(ValveSlicePosition).ValveSlice= ValveSliceT;
              imDesired.ValveSlice = ValveSlice; 
              
          %% that is for RVOT
          elseif strcmp(spec_sel, 'RVOT')    
              %%%load the image into SASlice 
              RVOTSlicePosition = RVOTSlicePosition + 1;
              dicomDir = sprintf('%s\\%s',studyDir_sel, ImgDir_sel);
              %if ismac
                  dicomDir = sprintf('%s\\%s',studyDir_sel, ImgDir_sel);
              %end
              RVOTSliceT = loadAllDicomSAOrLAFromOneDir(dicomDir);
              RVOTSlice(RVOTSlicePosition).RVOTSlice= RVOTSliceT;
              imDesired.RVOTSlice = RVOTSlice; 
          end
          
          
       end
       
%         imData =  SXSliceSorted(1,imIndex).SXSlice(timeInstanceSelected).imData;
%         imInfo1 = SXSliceSorted(1,imIndex).SXSlice(timeInstanceSelected).imInfo;
%        SXSlice(SASlicePosition).SXSlice= SXSliceT;
         for studyIndex = 1 : SASlicePosition
             for timeInstIndex = 1 :  size(SXSlice(studyIndex).SXSlice, 1)   
                 SXSliceSorted(1,studyIndex).SXSlice(1,timeInstIndex).imData = SXSlice(studyIndex).SXSlice(timeInstIndex).imData;
                 SXSliceSorted(1,studyIndex).SXSlice(1,timeInstIndex).imInfo = SXSlice(studyIndex).SXSlice(timeInstIndex).imInfo;
                 SXSliceSorted(1,studyIndex).SXSlice(1,timeInstIndex).imIndex = SXSlice(studyIndex).SXSlice(timeInstIndex).imgIndex;
             end
         end
         
         for studyIndex = 1 : ValveSlicePosition
             for timeInstIndex = 1 : size(ValveSlice(studyIndex).ValveSlice, 1)
                 ValveSliceSorted(1,studyIndex).ValveSlice(1,timeInstIndex).imData = ValveSlice(studyIndex).ValveSlice(timeInstIndex).imData;
                 ValveSliceSorted(1,studyIndex).ValveSlice(1,timeInstIndex).imInfo = ValveSlice(studyIndex).ValveSlice(timeInstIndex).imInfo;
                 ValveSliceSorted(1,studyIndex).ValveSlice(1,timeInstIndex).imIndex = ValveSlice(studyIndex).ValveSlice(timeInstIndex).imgIndex;
             end
         end

         %LVOTSliceSorted
         %%now combine all the long axial image together into
         %%LVOTSliceSorted
         LATotalNo = 0;
         for LAIndex = 1 : LVOTSlicePosition
             LATotalNo = LATotalNo + 1;
              for timeInstIndex = 1 : size(LVOTSlice(LAIndex).LVOTSlice, 1)   
                  if ~isempty(LVOTSlice(LAIndex).LVOTSlice(timeInstIndex).imData)
                    LVOTSliceSorted(1,LATotalNo).LXSlice(timeInstIndex).imData = LVOTSlice(LAIndex).LVOTSlice(timeInstIndex).imData;
                    LVOTSliceSorted(1,LATotalNo).LXSlice(timeInstIndex).imInfo = LVOTSlice(LAIndex).LVOTSlice(timeInstIndex).imInfo;
                    LVOTSliceSorted(1,LATotalNo).LXSlice(timeInstIndex).imIndex = LVOTSlice(LAIndex).LVOTSlice(timeInstIndex).imgIndex;
                  end
              end
         end
         
         for LAIndex = 1 : FourChSlicePosition
             LATotalNo = LATotalNo + 1;
              for timeInstIndex = 1 :  size(FourCHSlice(LAIndex).FourCHSlice, 1)
                  if ~isempty(FourCHSlice(LAIndex).FourCHSlice(timeInstIndex).imData)
                    LVOTSliceSorted(1,LATotalNo).LXSlice(timeInstIndex).imData = FourCHSlice(LAIndex).FourCHSlice(timeInstIndex).imData;
                    LVOTSliceSorted(1,LATotalNo).LXSlice(timeInstIndex).imInfo = FourCHSlice(LAIndex).FourCHSlice(timeInstIndex).imInfo;
                    LVOTSliceSorted(1,LATotalNo).LXSlice(timeInstIndex).imIndex = FourCHSlice(LAIndex).FourCHSlice(timeInstIndex).imgIndex;
                  end
              end
         end
         
         for LAIndex = 1 : OneChSlicePosition
             LATotalNo = LATotalNo + 1;
              for timeInstIndex = 1 : size(OneCHSlice(LAIndex).OneCHSlice, 1)  
                   if ~isempty(OneCHSlice(LAIndex).OneCHSlice(timeInstIndex).imData)
                    LVOTSliceSorted(1,LATotalNo).LXSlice(timeInstIndex).imData = OneCHSlice(LAIndex).OneCHSlice(timeInstIndex).imData;
                    LVOTSliceSorted(1,LATotalNo).LXSlice(timeInstIndex).imInfo = OneCHSlice(LAIndex).OneCHSlice(timeInstIndex).imInfo;
                    LVOTSliceSorted(1,LATotalNo).LXSlice(timeInstIndex).imIndex = OneCHSlice(LAIndex).OneCHSlice(timeInstIndex).imgIndex;
                   end
              end
         end

         
         %% this for RVOT, to save sorted position, nothing done currently
          for studyIndex = 1 : RVOTSlicePosition
             for timeInstIndex = 1 : size(RVOTSlice(studyIndex).RVOTSlice, 1)
                 RVOTSliceSorted(1,studyIndex).RVOTSlice(1,timeInstIndex).imData = RVOTSlice(studyIndex).RVOTSlice(timeInstIndex).imData;
                 RVOTSliceSorted(1,studyIndex).RVOTSlice(1,timeInstIndex).imInfo = RVOTSlice(studyIndex).RVOTSlice(timeInstIndex).imInfo;
                 RVOTSliceSorted(1,studyIndex).RVOTSlice(1,timeInstIndex).imIndex = RVOTSlice(studyIndex).RVOTSlice(timeInstIndex).imgIndex;
             end
         end

       %%%save all the data
           cd(resultDir);
           save imDesired imDesired SXSliceSorted LVOTSliceSorted;
           cd(workingDir);
          
end



%%%now we will need to decide the end-diastole, early-diastole and
%%%end-systole for each cine series, which may be different in the
%%%acquisition time 

%% get the trigger time for determing early-diastole and end-systole from LVOT images
LVOT_index = 1;
FourCH_index = 2;
OneCH_index = 3;

LVOTSlice = LVOTSliceSorted(1,LVOT_index).LXSlice;
trigger_time_early_diastole = 0.0;
trigger_time_end_systole = 0.0;
for i = 1 : size(LVOTSlice,2)
    if LVOTSlice(1,i).imIndex == patientConfigs(patientIndex,1).TimeEarlyOfDiastole 
        trigger_time_early_diastole = LVOTSlice(1,i).imInfo.TriggerTime;
        LVOTSliceSorted(1,LVOT_index).TimeEarlyOfDiastole = i;
    end
    if LVOTSlice(1,i).imIndex == patientConfigs(patientIndex,1).TimeEndOfSystole 
        trigger_time_end_systole = LVOTSlice(1,i).imInfo.TriggerTime;
        LVOTSliceSorted(1,LVOT_index).TimeEndOfSystole = i;
    end
end
 LVOTSliceSorted(1,LVOT_index).TimeEndOfDiastole = 1;

 trigger_times = [];
FourCHSlice = LVOTSliceSorted(1,FourCH_index).LXSlice;
for i = 1 : size(FourCHSlice,2)
    if ~isempty(FourCHSlice(1,i).imData)
        trigger_times(i) = FourCHSlice(1,i).imInfo.TriggerTime;
    else
        trigger_times(i) = 0;
    end  
end
[~, imIndex] = min(abs(trigger_times-trigger_time_early_diastole));
LVOTSliceSorted(1,FourCH_index).TimeEarlyOfDiastole = imIndex;
[~, imIndex] = min(abs(trigger_times-trigger_time_end_systole));
LVOTSliceSorted(1,FourCH_index).TimeEndOfSystole = imIndex;
LVOTSliceSorted(1,FourCH_index).TimeEndOfDiastole = 1;

trigger_times = [];
FourCHSlice = LVOTSliceSorted(1,FourCH_index).LXSlice;
for i = 1 : size(FourCHSlice,2)
    if ~isempty(FourCHSlice(1,i).imData)
        trigger_times(i) = FourCHSlice(1,i).imInfo.TriggerTime;
    else
        trigger_times(i) = 0;
    end  
end
[~, imIndex] = min(abs(trigger_times-trigger_time_early_diastole));
LVOTSliceSorted(1,FourCH_index).TimeEarlyOfDiastole = imIndex;
[~, imIndex] = min(abs(trigger_times-trigger_time_end_systole));
LVOTSliceSorted(1,FourCH_index).TimeEndOfSystole = imIndex;
LVOTSliceSorted(1,FourCH_index).TimeEndOfDiastole = 1;

trigger_times = [];
OneCHSlice = LVOTSliceSorted(1,OneCH_index).LXSlice;
for i = 1 : size(OneCHSlice,2)
    if ~isempty(OneCHSlice(1,i).imData)
        trigger_times(i) = OneCHSlice(1,i).imInfo.TriggerTime;
    else
        trigger_times(i) = 0;
    end  
end
[~, imIndex] = min(abs(trigger_times-trigger_time_early_diastole));
LVOTSliceSorted(1,OneCH_index).TimeEarlyOfDiastole = imIndex;
[~, imIndex] = min(abs(trigger_times-trigger_time_end_systole));
LVOTSliceSorted(1,OneCH_index).TimeEndOfSystole = imIndex;
LVOTSliceSorted(1,OneCH_index).TimeEndOfDiastole = 1;

%%figure out for short axis images 
for SAIndex = 1 : size(SXSliceSorted,2)
    SXSlice = SXSliceSorted(1,SAIndex).SXSlice;
    trigger_times = [];
    for i = 1 : size(SXSlice,2)
        if ~isempty(SXSlice(1,i).imData)
            trigger_times(i) = SXSlice(1,i).imInfo.TriggerTime;
        else
            trigger_times(i) = 0;
        end  
    end
    [~, imIndex] = min(abs(trigger_times-trigger_time_early_diastole));
    SXSliceSorted(1,SAIndex).TimeEarlyOfDiastole = imIndex;
    [~, imIndex] = min(abs(trigger_times-trigger_time_end_systole));
    SXSliceSorted(1,SAIndex).TimeEndOfSystole = imIndex;
    SXSliceSorted(1,SAIndex).TimeEndOfDiastole = 1;
end



%%figure out for short axis valve images
if exist('ValveSliceSorted', 'var')
    for ValveIndex = 1 : size(ValveSliceSorted,2)
        ValveSlice = ValveSliceSorted(1,ValveIndex).ValveSlice;
        trigger_times = [];
        for i = 1 : size(ValveSlice,2)
            if ~isempty(ValveSlice(1,i).imData)
                trigger_times(i) = ValveSlice(1,i).imInfo.TriggerTime;
            else
                trigger_times(i) = 0;
            end  
        end
        [~, imIndex] = min(abs(trigger_times-trigger_time_early_diastole));
        ValveSliceSorted(1,ValveIndex).TimeEarlyOfDiastole = imIndex;
        [~, imIndex] = min(abs(trigger_times-trigger_time_end_systole));
        ValveSliceSorted(1,ValveIndex).TimeEndOfSystole = imIndex;
        ValveSliceSorted(1,ValveIndex).TimeEndOfDiastole = 1;
    end
else
    ValveSliceSorted = [];
end


%%figure out for RVOT images
if exist('RVOTSliceSorted', 'var')
    for RVOTIndex = 1 : size(RVOTSliceSorted,2)
        RVOTSlice = RVOTSliceSorted(1,RVOTIndex).RVOTSlice;
        trigger_times = [];
        for i = 1 : size(RVOTSlice,2)
            if ~isempty(RVOTSlice(1,i).imData)
                trigger_times(i) = RVOTSlice(1,i).imInfo.TriggerTime;
            else
                trigger_times(i) = 0;
            end  
        end
        [~, imIndex] = min(abs(trigger_times-trigger_time_early_diastole));
        RVOTSliceSorted(1,RVOTIndex).TimeEarlyOfDiastole = imIndex;
        [~, imIndex] = min(abs(trigger_times-trigger_time_end_systole));
        RVOTSliceSorted(1,RVOTIndex).TimeEndOfSystole = imIndex;
        RVOTSliceSorted(1,RVOTIndex).TimeEndOfDiastole = 1;
    end
else
    RVOTSliceSorted = [];
end

%%resave again after determing end-systole and early-diastole
cd(resultDir);
save imDesired imDesired SXSliceSorted LVOTSliceSorted ValveSliceSorted RVOTSliceSorted;
cd(workingDir);













