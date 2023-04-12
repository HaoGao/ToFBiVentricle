clear all; clc; close all
% LVWM_config; %%this is a general set up
LVWM_config;

path(path,'./bsplineDeform');
path(path, './demonCode');

%%patientConfigs(patientIndex,1).basalSliceIndex = 4
sliceSelected = 2; %%% the range is from 1 to usuableSXSlice

EndoEpiBCManualSegBool = 1;

%%%setup a result dir for deformable image Registration 
cd(resultDir);
if ~exist(bSplineDir, 'dir')
    mkdir(bSplineDir);
end
cd(bSplineDir);
dirName = sprintf('deformRes_4CH');
if ~exist(dirName,'dir')
    mkdir(dirName);
end
cd(dirName);
ResDefDir = pwd();
cd(workingDir);

cd(resultDir);
load imDesired;
cd(workingDir);

% %%now need to choose phase to segment
% list_phase = {'early_diastole', 'end_diastole', 'end_systole'};
% [idx, tf] = listdlg('ListString', list_phase);

%%%now need to extract the right data 
totalSXSliceLocation = size(SXSliceSorted, 2);
TimeEndOfSystole = LVOTSliceSorted(sliceSelected).TimeEndOfSystole; 
TimeEndOfDiastole = LVOTSliceSorted(sliceSelected).TimeEndOfDiastole; 
TimeEarlyOfDiastole =  LVOTSliceSorted(sliceSelected).TimeEarlyOfDiastole;
TimeInstanceSelected = size(LVOTSliceSorted(sliceSelected).LXSlice, 2);
% TimeInstanceSelected = TimeEarlyOfDiastole;
% TimeInstanceSelected = 35;

InstanceSelectedSeries = TimeEarlyOfDiastole-1 : 1: TimeEarlyOfDiastole+3;
totalInstanceNum = length(InstanceSelectedSeries);

InstanceNo = 0;
for   i = 1 : totalInstanceNum
    imInstance = InstanceSelectedSeries(i);
    imData = LVOTSliceSorted(1,sliceSelected).LXSlice(imInstance).imData;
    imInfo = LVOTSliceSorted(1,sliceSelected).LXSlice(imInstance).imInfo;
    
    if ~isempty(imData)
        InstanceNo = InstanceNo + 1;
        imgList(InstanceNo,1).imData = imData;
        imgList(InstanceNo,1).imInfo = imInfo;
    end
end
    
cd(workingDir);

%%total Img no
imgTotalNo = size(imgList,1);

[I1Crop, rect] = imcrop(imgList(1,1).imData,[]);


 %%%segment the first image
 imNo = 1;
 
if EndoEpiBCManualSegBool == 1
    for imNo = 1 : size(imgList,1)
        
        
        imDataT = imcrop(imgList(imNo,1).imData,rect); 
        hIm = figure();
        imshow(imDataT, []); hold on;
        hold on;
        endo_c=[];
        epi_c=[];
        endo_sample = [];
        epi_sample = [];
        but = 1;
        n=1;
        
        if imNo >1
            endo_sample = Annu(imNo-1,1).annu_1 ;
            epi_sample = Annu(imNo-1,1).annu_2;
            endo_sample(1,:) = endo_sample(1,:) -rect(1);
            endo_sample(2,:) = endo_sample(2,:) -rect(2);
            epi_sample(1,:) = epi_sample(1,:) -rect(1);
            epi_sample(2,:) = epi_sample(2,:) -rect(2);
            plot(endo_sample(1,:) , endo_sample(2,:) , 'r-o', 'MarkerSize', 6); hold on;
            plot(epi_sample(1,:) , epi_sample(2,:) , 'r-o', 'MarkerSize', 6); hold on;
            clear   endo_sample epi_sample;
        end
        
       
%         while but ==1 
                    [x, y, but]=myginput(1,'crosshair');
                    plot(x,y,'r+');hold on;
                    endo_c(1,n)=x;
                    endo_c(2,n)=y;
                    n=n+1;
%         end
        plot(endo_c(1,:),endo_c(2,:),'b');
                
                
        n=1;but=1;
%                 while but==1
                    [x, y, but]=myginput(1, 'crosshair');
                    plot(x,y,'r+');hold on;
                    epi_c(1,n)=x;
                    epi_c(2,n)=y;
                    n=n+1;
%                 end
                
        endo_sample(1,:) = endo_c(1,:) +rect(1);
        endo_sample(2,:) = endo_c(2,:) +rect(2);
        epi_sample(1,:) = epi_c(1,:) +rect(1);
        epi_sample(2,:) = epi_c(2,:) +rect(2);
        
        
        
        
        Annu(imNo,1).annu_1 = endo_sample;
        Annu(imNo,1).annu_2 = epi_sample;
        Annu(imNo,1).imData = imgList(imNo,1).imData;
        Annu(imNo,1).imInfo = imgList(imNo,1).imInfo;
        
        
        
        pause; 
        close all;
        clear epi_c endo_c endo_sample epi_sample;
    end
    cd(ResDefDir);
    save AnnuMotion Annu Annu;
    cd(workingDir);
end

cd(ResDefDir);
load AnnuMotion;
cd(workingDir);

ave_vel = 0;
for imNo = 1 : size(imgList,1)-1
     endoT = Annu(imNo,1).annu_1 ;
     epiT =  Annu(imNo,1).annu_2 ;
     if size(endoT,2) > 1
        shape_1 = [ mean(endoT(:,1)), mean(endoT(:,2))];
     else
         shape_1 = endoT;
     end
     if size(epiT,2) > 1
        shape_2 = [ mean(endoT(:,1)), mean(epiT(:,2))] ;
     else
         shape_2 = epiT;
     end
     
     
     endoT = Annu(imNo+1,1).annu_1 ;
     epiT =  Annu(imNo+1,1).annu_2 ;
      if size(endoT,2) > 1
        shape_1_next = [ mean(endoT(:,1)), mean(endoT(:,2))];
      else
         shape_1_next = endoT;
      end
      if size(epiT,2) > 1
        shape_2_next = [ mean(endoT(:,1)), mean(epiT(:,2))] ;
      else
         shape_2_next = epiT;
      end
     
     dx = Annu(imNo,1).imInfo.PixelSpacing.*0.1;
     dt = Annu(imNo+1,1).imInfo.TriggerTime - Annu(imNo,1).imInfo.TriggerTime;
     dt = dt*0.001;
     
     dx_annu_1 = (sum( (shape_1_next - shape_1).^2 ))^0.5;
     dx_annu_2 = (sum( (shape_2_next - shape_2).^2 ))^0.5;
     
     vel_annu_1 = dx_annu_1*dx(1)/dt;
     vel_annu_2 = dx_annu_2*dx(1)/dt;
     
     Annu(imNo,1).vel_annu_1 = vel_annu_1;
     Annu(imNo,1).vel_annu_2 = vel_annu_2;
     Annu(imNo,1).vel_annu = ( vel_annu_1+ vel_annu_2)/2;
     
     ave_vel = ave_vel + ( vel_annu_1+ vel_annu_2)/2;
end

ave_vel = ave_vel / (size(imgList,1)-1)

% shapex = mean(endoT(:,1));
% shapey = mean(endoT(:,2));
% 
%  h1 = figure; hold on;
%  imshow(imcrop(imgList(imNo,1).imData,imgDeformed(imNo,1).rect),[]);hold on;
%  plot(endoT(:,1),endoT(:,2),'-');
%  plot(epiT(:,1),epiT(:,2),'-');
%  F(imNo) = getframe;
%  
%  
%  for imNo = 1 : imgTotalNo-1
%     Tx = imgDeformed(imNo,1).Tx;
%     Ty = imgDeformed(imNo,1).Ty;
%     
%     endoTN = boundaryTracking(endoT, Tx, Ty);
%     epiTN = boundaryTracking(epiT, Tx, Ty);
%     endoT = endoTN;
%     epiT = epiTN;
%    
%     
%     pause;
% %     pause;
%     imshow(imcrop(imgList(imNo+1,1).imData,imgDeformed(imNo,1).rect),[]);
%     plot(endoT(:,1),endoT(:,2),'-');
%     plot(epiT(:,1),epiT(:,2),'r-');
%     F(imNo+1) = getframe; 
%     
%   
%  end
%  
%  
%  
%  