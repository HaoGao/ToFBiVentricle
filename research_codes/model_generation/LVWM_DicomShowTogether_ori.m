%%%%show the dicom in 3D together for one time instance
clear all;
close all;
clc;

LVWM_config;
Bwritten = 1; %%whether to write out solidworks file

cd(resultDir);
load imDesired;
cd(workingDir); 

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

cd(phase_resultDir);
load DataSegLA;
load DataSegSA;
if exist('DataSegValve.mat', 'file')
    load DataSegValve;
end
if exist('DataSegRVOT.mat', 'file')
    load DataSegRVOT;
end
cd(workingDir);



%%still need to rotate and translate, specify a SA slice which will be the
%%z=0 plane, and the center is in the left ventricle centre.


%%figure out which is the basal SA slice
basalSliceIndex = patientConfigs(patientIndex,1).basalSliceIndex;
if strcmp(list_phase{idx}, 'early_diastole')
     timeInstanceSelected =  SXSliceSorted(1,basalSliceIndex).TimeEarlyOfDiastole;
elseif  strcmp(list_phase{idx}, 'end_diastole')
    timeInstanceSelected = SXSliceSorted(1,basalSliceIndex).TimeEndOfDiastole;
elseif strcmp(list_phase{idx}, 'end_systole')
    timeInstanceSelected = SXSliceSorted(1,basalSliceIndex).TimeEndOfSystole;
else
    disp('could not determine the cardiac phase, quite')
    return;
end

imData =  SXSliceSorted(1,basalSliceIndex).SXSlice(timeInstanceSelected).imData;
imInfo1 = SXSliceSorted(1,basalSliceIndex).SXSlice(timeInstanceSelected).imInfo;
imInfo = infoExtract(imInfo1);

%%need to ask whether to re segment or not
answer = questdlg('would you like to re-segment','segment or skip');
if strcmp(answer, 'Yes') %%then segment 
    %[imCropData, rect] = imcrop(imData,[]);
    [imCropData, rect]= imcrop( uint8(double(imData)./(max(max(double(imData)))) *255) );
    hSA=figure();
    imshow(imCropData,[]);hold on;
    endo_c=[];
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
    nei=1:length(endo_c(1,:))+1;
    nnei=1:0.1:length(endo_c(1,:))+1;
    endo_cc(1,:)=spline(nei,[endo_c(1,:) endo_c(1,1)], nnei);
    endo_cc(2,:)=spline(nei,[endo_c(2,:) endo_c(2,1)], nnei);
    
    [endo_sample, ~] = samplingBCWithoutIm(endo_cc, endo_cc,100);
    endo_sample(1,:) = endo_sample(1,:) +rect(1);
    endo_sample(2,:) = endo_sample(2,:) +rect(2);
    
    clf(hSA); hold on; figure(hSA);
    imshow(imData,[]); hold on;
    plot(endo_sample(1,:),endo_sample(2,:),'r');
    
    endo_cReal = TransformCurvesFromImToRealSpace(endo_sample,imInfo);
   
    
    basal_centre.xcoor = mean(endo_cReal(1,:));
    basal_centre.ycoor = mean(endo_cReal(2,:));
    basal_centre.zcoor = mean(endo_cReal(3,:));
    basal_centre.centre_coor = [basal_centre.xcoor; basal_centre.ycoor; ...
                                basal_centre.zcoor];
    basal_centre.endo = endo_sample;
    basal_centre.endo_cReal = endo_cReal;
    basal_centre.imIndex = basalSliceIndex;
    
    cd(phase_resultDir);
    save basal_centre basal_centre ;
    cd(workingDir);
else
    cd(phase_resultDir);
    load basal_centre ;
    cd(workingDir);

end

if patientConfigs(patientIndex,1).Brotation
    %%we can also define a coordinate system, 
    %% x points to the horizotal direction 
    SAXVec = NormalizationVec(basal_centre.endo_cReal(:,1) - basal_centre.centre_coor);
    SAYVec = NormalizationVec(basal_centre.endo_cReal(:,25) - basal_centre.centre_coor);
    LAVec = cross(SAXVec, SAYVec);
    %% find the apex slice to determine the long-axis direction 
    apex_index = basalSliceIndex+3;
    LVApexCenterApex = [mean(DataSegSA(1,apex_index).endo_cReal(1,:)); ...
        mean(DataSegSA(1,apex_index).endo_cReal(2,:));... 
        mean(DataSegSA(1,apex_index).endo_cReal(3,:))];
    long_axis = NormalizationVec(basal_centre.centre_coor -LVApexCenterApex );
    if dot(LAVec, long_axis)< 0
        LAVec = - LAVec;
    end
    SAYVec = cross(LAVec, SAXVec);
    
    
    leftM = [SAXVec(1) SAYVec(1) LAVec(1);
             SAXVec(2) SAYVec(2) LAVec(2);
             SAXVec(3) SAYVec(3) LAVec(3)];
    rightM = [1 0 0;
              0 1 0;
              0 0 1];
    RotationMatrix = rightM/leftM;
    
    rotationConfig.SAXVec = SAXVec;
    rotationConfig.SAYVec = SAYVec;
    rotationConfig.LAVec = LAVec;
    rotationConfig.RotationMatrix = RotationMatrix;
    
    cd(phase_resultDir);
    save rotationConfig rotationConfig ;
    cd(workingDir);
end

close all;

%%plot original 3D boundaries 
h3D_ori = figure(); hold on; axis equal;
for imIndex = 1 : size(DataSegSA, 2)
    endo_c = DataSegSA(imIndex).endo_cReal;
    epi_c = DataSegSA(imIndex).epi_cReal;
  
    if ~isempty(endo_c)
        plot3(endo_c(1,:),endo_c(2,:), endo_c(3,:),'LineStyle', '-', 'Color', ...
            'b', 'LineWidth',2);
        plot3(epi_c(1,:),epi_c(2,:), epi_c(3,:),'LineStyle', '-', 'Color', ...
            'r', 'LineWidth',2);
    end
end
for imIndex = 1 : size(DataSegLA,2)
    endo_c = DataSegLA(imIndex).endo_cReal;
    epi_c = DataSegLA(imIndex).epi_cReal;
  
    if ~isempty(endo_c) 
        plot3(endo_c(1,:),endo_c(2,:), endo_c(3,:),'LineStyle', '-', 'Color', ...
            'b', 'LineWidth',2);
        plot3(epi_c(1,:),epi_c(2,:), epi_c(3,:),'LineStyle', '-', 'Color', ...
            'r', 'LineWidth',2);
    end
end
figure(h3D_ori); xlabel('x'); ylabel('y'); zlabel('z'); title('Original BCs');


if patientConfigs(patientIndex,1).Brotation 
%%rotate BCs

    cd(phase_resultDir);
    load rotationConfig;
    cd(workingDir);
    RotationMatrix = rotationConfig.RotationMatrix;
    
    DataSegSARotated = DataSegSA;
    DataSegLARotated = DataSegLA;
    if exist('DataSegValve', 'var')
        DataSegValveRotated = DataSegValve;
    end
    h3D_rotated = figure(); hold on; axis equal; 
    LVUpperCenter = basal_centre.centre_coor;
    for imIndex = 1 : size(DataSegSA, 2)
        endo_c = DataSegSA(imIndex).endo_cReal;
        epi_c = DataSegSA(imIndex).epi_cReal;
        
        if ~isempty(endo_c) && ~isempty(epi_c)

            endo_cT(1,:) = endo_c(1,:)-LVUpperCenter(1);
            endo_cT(2,:) = endo_c(2,:)-LVUpperCenter(2);
            endo_cT(3,:) = endo_c(3,:)-LVUpperCenter(3);
            endo_c = RotationMatrix*endo_cT;

            epi_cT(1,:) = epi_c(1,:)-LVUpperCenter(1);
            epi_cT(2,:) = epi_c(2,:)-LVUpperCenter(2);
            epi_cT(3,:) = epi_c(3,:)-LVUpperCenter(3);
            epi_c = RotationMatrix*epi_cT; 

            DataSegSARotated(imIndex).endo_cReal = endo_c;
            DataSegSARotated(imIndex).epi_cReal = epi_c;

            figure(h3D_rotated);hold on;
            plot3(endo_c(1,:),endo_c(2,:), endo_c(3,:),'LineStyle', '-', 'Color', ...
                        'b', 'LineWidth',2);
            plot3(epi_c(1,:),epi_c(2,:), epi_c(3,:),'LineStyle', '-', 'Color', ...
                        'r', 'LineWidth',2);

            if Bwritten
                %%%output the bc boundaries
                cd(phase_resultDir);
                filename = sprintf('SABCdata_slice%d.txt',imIndex);
                fidSA = fopen(filename,'w');
                lined = 0;
                if imIndex >= patientConfigs(patientIndex,1).basalSliceIndex
                    lined = 1;
                end
                cd(workingDir);
                SolidWorksInputForClosedSplineGeneration(fidSA,endo_c, 0, lined, 1);
                SolidWorksInputForClosedSplineGeneration(fidSA,epi_c, 0, lined, 0);
                fclose(fidSA);
            end

            clear endo_cT;
            clear epi_cT;
        end

    end


    for imIndex = 1 : size(DataSegLA, 2)
        endo_c = DataSegLA(imIndex).endo_cReal;
        epi_c = DataSegLA(imIndex).epi_cReal;

        endo_cT(1,:) = endo_c(1,:)-LVUpperCenter(1);
        endo_cT(2,:) = endo_c(2,:)-LVUpperCenter(2);
        endo_cT(3,:) = endo_c(3,:)-LVUpperCenter(3);
        endo_c = RotationMatrix*endo_cT;

        epi_cT(1,:) = epi_c(1,:)-LVUpperCenter(1);
        epi_cT(2,:) = epi_c(2,:)-LVUpperCenter(2);
        epi_cT(3,:) = epi_c(3,:)-LVUpperCenter(3);
        epi_c = RotationMatrix*epi_cT; 

        DataSegLARotated(imIndex).endo_cReal = endo_c;
        DataSegLARotated(imIndex).epi_cReal = epi_c;

        figure(h3D_rotated);hold on;
        plot3(endo_c(1,:),endo_c(2,:), endo_c(3,:),'LineStyle', '-', 'Color', ...
                    'b', 'LineWidth',2);
        plot3(epi_c(1,:),epi_c(2,:), epi_c(3,:),'LineStyle', '-', 'Color', ...
                    'r', 'LineWidth',2);

        if Bwritten
            %%%output the bc boundaries
            cd(phase_resultDir);
            filename = sprintf('LABCdata_slice%d.txt',imIndex);
            fidSA = fopen(filename,'w');
            cd(workingDir);
            SolidWorksInputForClosedSplineGeneration(fidSA,endo_c, 0, 1, 1);
            SolidWorksInputForClosedSplineGeneration(fidSA,epi_c, 0, 1, 0);
            fclose(fidSA);
        end

        clear endo_cT;
        clear epi_cT;

    end
    
    %%for vavle 
    if exist('DataSegValve', 'var')
        for imIndex = 1 : size(DataSegValve, 2)
            
            endo_c = DataSegValve(imIndex).endo_cReal;
            epi_c = DataSegValve(imIndex).epi_cReal;
            if ~isempty(endo_c)

                endo_cT(1,:) = endo_c(1,:)-LVUpperCenter(1);
                endo_cT(2,:) = endo_c(2,:)-LVUpperCenter(2);
                endo_cT(3,:) = endo_c(3,:)-LVUpperCenter(3);
                endo_c = RotationMatrix*endo_cT;

                epi_cT(1,:) = epi_c(1,:)-LVUpperCenter(1);
                epi_cT(2,:) = epi_c(2,:)-LVUpperCenter(2);
                epi_cT(3,:) = epi_c(3,:)-LVUpperCenter(3);
                epi_c = RotationMatrix*epi_cT; 

                DataSegValveRotated(imIndex).endo_cReal = endo_c;
                DataSegValveRotated(imIndex).epi_cReal = epi_c;

                figure(h3D_rotated);hold on;
                plot3(endo_c(1,:),endo_c(2,:), endo_c(3,:),'LineStyle', '-', 'Color', ...
                            'b', 'LineWidth',2);
                plot3(epi_c(1,:),epi_c(2,:), epi_c(3,:),'LineStyle', '-', 'Color', ...
                            'r', 'LineWidth',2);

                if Bwritten
                    %%%output the bc boundaries
                    cd(phase_resultDir);
                    filename = sprintf('ValveBCdata_slice%d.txt',imIndex);
                    fidSA = fopen(filename,'w');
                    cd(workingDir);
                    SolidWorksInputForClosedSplineGeneration(fidSA,endo_c, 0, 1, 1);
                    SolidWorksInputForClosedSplineGeneration(fidSA,epi_c, 0, 1, 0);
                    fclose(fidSA);
                end

                clear endo_cT;
                clear epi_cT;
            end

        end
    end
    
%     cd(phase_resultDir);
%     if exist('DataSegValve', 'var')
%         save rotatedBCs DataSegLARotated DataSegSARotated DataSegValveRotated;
%     else
%         save rotatedBCs DataSegLARotated DataSegSARotated;
%     end
%     cd(workingDir);
    
    
     %%for RVOT 
    if exist('DataSegRVOT', 'var')
        for imIndex = 1 : size(DataSegRVOT, 2)
            
            endo_c = DataSegRVOT(imIndex).endo_cReal;
            epi_c = DataSegRVOT(imIndex).epi_cReal;
            if ~isempty(endo_c)

                endo_cT(1,:) = endo_c(1,:)-LVUpperCenter(1);
                endo_cT(2,:) = endo_c(2,:)-LVUpperCenter(2);
                endo_cT(3,:) = endo_c(3,:)-LVUpperCenter(3);
                endo_c = RotationMatrix*endo_cT;

                epi_cT(1,:) = epi_c(1,:)-LVUpperCenter(1);
                epi_cT(2,:) = epi_c(2,:)-LVUpperCenter(2);
                epi_cT(3,:) = epi_c(3,:)-LVUpperCenter(3);
                epi_c = RotationMatrix*epi_cT; 

                DataSegRVOTRotated(imIndex).endo_cReal = endo_c;
                DataSegRVOTRotated(imIndex).epi_cReal = epi_c;

                figure(h3D_rotated);hold on;
                plot3(endo_c(1,:),endo_c(2,:), endo_c(3,:),'LineStyle', '-', 'Color', ...
                            'b', 'LineWidth',2);
                plot3(epi_c(1,:),epi_c(2,:), epi_c(3,:),'LineStyle', '-', 'Color', ...
                            'r', 'LineWidth',2);

                if Bwritten
                    %%%output the bc boundaries
                    cd(phase_resultDir);
                    filename = sprintf('RVOTBCdata_slice%d.txt',imIndex);
                    fidSA = fopen(filename,'w');
                    cd(workingDir);
                    SolidWorksInputForClosedSplineGeneration(fidSA,endo_c, 0, 1, 1);
                    SolidWorksInputForClosedSplineGeneration(fidSA,epi_c, 0, 1, 0);
                    fclose(fidSA);
                end

                clear endo_cT;
                clear epi_cT;
            end

        end
    end
    
    %% save appropriate results
    cd(phase_resultDir);
    if exist('DataSegRVOT', 'var') && exist('DataSegValve', 'var')
        save rotatedBCs DataSegLARotated DataSegSARotated DataSegValveRotated DataSegRVOTRotated;
    elseif exist('DataSegRVOT', 'var') && ~exist('DataSegValve', 'var')
        save rotatedBCs DataSegLARotated DataSegSARotated DataSegRVOTRotated;
    elseif ~exist('DataSegRVOT', 'var') && exist('DataSegValve', 'var')
        save rotatedBCs DataSegLARotated DataSegSARotated DataSegValveRotated;
    else
        save rotatedBCs DataSegLARotated DataSegSARotated;
    end
    cd(workingDir);

end
















