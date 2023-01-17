clear all;close all;clc;
%% it only sets up the segments at one mid-ventricular slice and one apical slice;
%% the mid slice: basalSliceIndex + 2
%% the apical slice: basalSliceIndex + 5 this may need to update for some subjects

% LVWM_config;
LVWM_config;

cd(resultDir);
load imDesired;
% load rotationConfig; %%this is from segmentation part
cd(workingDir);

timeInstanceSelectedDiastole = 1;
SASlicePositonMiddle = patientConfigs(1,1).basalSliceIndex + 1;
SASlicePositonApex = patientConfigs(1,1).basalSliceIndex + 5;

for positionIndex = 1 : 2
    if positionIndex == 1
        imInfo1 = SXSliceSorted(SASlicePositonMiddle).SXSlice(timeInstanceSelectedDiastole).imInfo;
        imInfo = infoExtract(imInfo1);
        imData = SXSliceSorted(SASlicePositonMiddle).SXSlice(timeInstanceSelectedDiastole).imData;
    else
        imInfo1 = SXSliceSorted(SASlicePositonApex).SXSlice(timeInstanceSelectedDiastole).imInfo;
        imInfo = infoExtract(imInfo1);
        imData = SXSliceSorted(SASlicePositonApex).SXSlice(timeInstanceSelectedDiastole).imData;
    end
    


    h1=figure(); imshow(imData,[]);hold on;pause;
    %%define 6 points or 4 points, starting from inferior septal
    if positionIndex ==1 
        endo_c=[];
        but = 1;
        n=1;
        while but ==1 && n<=6
                [x y but]=ginput(1);
                if n==1 
                    plot(x,y,'b.');hold on;
                elseif n==2
                    plot(x,y,'r+'); hold on;
                elseif n==3
                    plot(x,y,'y*');hold on;
                elseif n==4
                    plot(x,y,'b<'); hold on;
                else
                    plot(x,y,'r.');hold on;
                end
                endo_c(1,n)=x;
                endo_c(2,n)=y;
                n=n+1;
        end

        %%%degree calculation
        theta = [];
        centerPoint = [mean(endo_c(1,:)) mean(endo_c(2,:))];
        for i = 1 : size(endo_c,2)
            p = [endo_c(1,i),endo_c(2,i)];
            theta(i) = degreeCalculationPointBased(p,centerPoint)*180/pi; %%in the range of 0-360
        end

        %%%segment region
        MidConfig.InfSeptTheta = degreeReOrder(theta(2),theta(1));
        MidConfig.AntSeptTheta = degreeReOrder(theta(3),theta(2));
        MidConfig.AntTheta = degreeReOrder(theta(4),theta(3));
        MidConfig.AntLatTheta = degreeReOrder(theta(5),theta(4));
        MidConfig.InfLatTheta = degreeReOrder(theta(6),theta(5));
        MidConfig.InfTheta = degreeReOrder(theta(6),theta(1));
        MidConfig.endo_c = endo_c;
        MidConfig.theta = theta;
    else %%this is for the apex slices
        endo_c=[];
        but = 1;
        n=1;
        while but ==1 && n<=4
                [x y but]=ginput(1);
                if n==1 
                    plot(x,y,'b.');hold on;
                elseif n==2
                    plot(x,y,'r+'); hold on;
                elseif n==3
                    plot(x,y,'y*');hold on;
                elseif n==4
                    plot(x,y,'b<'); hold on;
                else
                    plot(x,y,'r.');hold on;
                end
                endo_c(1,n)=x;
                endo_c(2,n)=y;
                n=n+1;
        end

        %%%degree calculation
        centerPoint = [mean(endo_c(1,:)) mean(endo_c(2,:))];
        theta = [];
        for i = 1 : size(endo_c,2)
            p = [endo_c(1,i),endo_c(2,i)];
            theta(i) = degreeCalculationPointBased(p,centerPoint)*180/pi; %%in the range of 0-360
        end

        %%%segment region
        ApexConfig.SeptTheta = degreeReOrder(theta(2),theta(1));
        ApexConfig.AntTheta = degreeReOrder(theta(3),theta(2));
        ApexConfig.Lat = degreeReOrder(theta(4),theta(3));
        ApexConfig.Inf = degreeReOrder(theta(1),theta(4)); 
        ApexConfig.endo_c = endo_c;
        ApexConfig.theta = theta;
    end
end


cd(resultDir)
save DivisionConfig MidConfig ApexConfig;
cd(workingDir);

