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
    disp('could not determine the cardiac phase, quit')
    return;
end

imData =  SXSliceSorted(1,basalSliceIndex).SXSlice(timeInstanceSelected).imData;
imInfo1 = SXSliceSorted(1,basalSliceIndex).SXSlice(timeInstanceSelected).imInfo;
imInfo = infoExtract(imInfo1);

%%need to ask whether to re segment or not

cd(phase_resultDir);
load basal_centre ;
cd(workingDir);



if patientConfigs(patientIndex,1).Brotation
    %%we can also define a coordinate system, 
    %% x points to the horizotal direction 
    SAXVec = NormalizationVec(basal_centre.endo_cReal(:,1) - basal_centre.centre_coor);
    SAYVec = NormalizationVec(basal_centre.endo_cReal(:,25) - basal_centre.centre_coor);
    LAVec = cross(SAXVec, SAYVec);
    %% find the apex slice to determine the long-axis direction 
    apex_index = basalSliceIndex+3;
    LVApexCenterApex = [mean([DataSegSA(1,apex_index).endo_lvReal(1,:) DataSegSA(1,apex_index).endo_rvReal(1,:)]); ...
        mean([DataSegSA(1,apex_index).endo_lvReal(2,:) DataSegSA(1,apex_index).endo_rvReal(2,:)]);... 
        mean([DataSegSA(1,apex_index).endo_lvReal(3,:) DataSegSA(1,apex_index).endo_rvReal(3,:)])];
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
    endo_lv = DataSegSA(imIndex).endo_lvReal;
    endo_rv = DataSegSA(imIndex).endo_rvReal;
    endo_pa = DataSegSA(imIndex).endo_paReal;
    epi_c = DataSegSA(imIndex).epi_cReal;
    epi_crv = DataSegSA(imIndex).epi_crvReal;
    epi_cpa = DataSegSA(imIndex).epi_cpaReal;
  
    if ~isempty(endo_lv)
        plot3(endo_lv(1,:),endo_lv(2,:), endo_lv(3,:),'LineStyle', '-', 'Color', ...
            'b', 'LineWidth',2);
    end
    if ~isempty(endo_rv)
        plot3(endo_rv(1,:),endo_rv(2,:), endo_rv(3,:),'LineStyle', '-', 'Color', ...
            'g', 'LineWidth',2);
    end
    if ~isempty(endo_pa)
        plot3(endo_pa(1,:),endo_pa(2,:), endo_pa(3,:),'LineStyle', '-', 'Color', ...
            'g', 'LineWidth',2);
    end
    if ~isempty(epi_c)
        plot3(epi_c(1,:),epi_c(2,:), epi_c(3,:),'LineStyle', '-', 'Color', ...
            'r', 'LineWidth',2);
    end
    if ~isempty(epi_crv)
        plot3(epi_crv(1,:),epi_crv(2,:), epi_crv(3,:),'LineStyle', '-', 'Color', ...
            'r', 'LineWidth',2);
    end
    if ~isempty(epi_cpa)
        plot3(epi_cpa(1,:),epi_cpa(2,:), epi_cpa(3,:),'LineStyle', '-', 'Color', ...
            'r', 'LineWidth',2);
    end
    
end
for imIndex = 1 : size(DataSegLA,2)
    endo_lv = DataSegLA(imIndex).endo_lvReal;
    endo_rv = DataSegLA(imIndex).endo_rvReal;
    epi_c = DataSegLA(imIndex).epi_cReal;
  
    if ~isempty(endo_lv) 
        plot3(endo_lv(1,:),endo_lv(2,:), endo_lv(3,:),'LineStyle', '-', 'Color', ...
            'b', 'LineWidth',2);
        plot3(endo_rv(1,:),endo_rv(2,:), endo_rv(3,:),'LineStyle', '-', 'Color', ...
            'y', 'LineWidth',2);
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
        endo_lv = DataSegSA(imIndex).endo_lvReal;
        endo_rv = DataSegSA(imIndex).endo_rvReal;
        endo_pa = DataSegSA(imIndex).endo_paReal;
        epi_c = DataSegSA(imIndex).epi_cReal;
        epi_crv = DataSegSA(imIndex).epi_crvReal;
        epi_cpa = DataSegSA(imIndex).epi_cpaReal;
        
        figure(h3D_rotated);hold on;
        if ~isempty(endo_lv)|| ~isempty(endo_rv) || ~isempty(epi_c) ...
           || ~isempty(endo_pa) || ~isempty(epi_crv) || ~isempty(endo_pa) 
            if ~isempty(endo_lv)
            endo_lvT(1,:) = endo_lv(1,:)-LVUpperCenter(1);
            endo_lvT(2,:) = endo_lv(2,:)-LVUpperCenter(2);
            endo_lvT(3,:) = endo_lv(3,:)-LVUpperCenter(3);
            endo_lv = RotationMatrix*endo_lvT;
            plot3(endo_lv(1,:),endo_lv(2,:), endo_lv(3,:),'LineStyle', '-', 'Color', ...
                        'b', 'LineWidth',2);
            end
            
            if ~isempty(endo_rv)
            endo_rvT(1,:) = endo_rv(1,:)-LVUpperCenter(1);
            endo_rvT(2,:) = endo_rv(2,:)-LVUpperCenter(2);
            endo_rvT(3,:) = endo_rv(3,:)-LVUpperCenter(3);
            endo_rv = RotationMatrix*endo_rvT;
            plot3(endo_rv(1,:),endo_rv(2,:), endo_rv(3,:),'LineStyle', '-', 'Color', ...
                        'g', 'LineWidth',2);
            end
            
            if ~isempty(endo_pa)
            endo_paT(1,:) = endo_pa(1,:)-LVUpperCenter(1);
            endo_paT(2,:) = endo_pa(2,:)-LVUpperCenter(2);
            endo_paT(3,:) = endo_pa(3,:)-LVUpperCenter(3);
            endo_pa = RotationMatrix*endo_paT;
            plot3(endo_pa(1,:),endo_pa(2,:), endo_pa(3,:),'LineStyle', '-', 'Color', ...
                        'g', 'LineWidth',2);
            end
            
            if ~isempty(epi_c)
            epi_cT(1,:) = epi_c(1,:)-LVUpperCenter(1);
            epi_cT(2,:) = epi_c(2,:)-LVUpperCenter(2);
            epi_cT(3,:) = epi_c(3,:)-LVUpperCenter(3);
            epi_c = RotationMatrix*epi_cT; 
            plot3(epi_c(1,:),epi_c(2,:), epi_c(3,:),'LineStyle', '-', 'Color', ...
                        'r', 'LineWidth',2);
            end

            if ~isempty(epi_crv)
            epi_crvT(1,:) = epi_crv(1,:)-LVUpperCenter(1);
            epi_crvT(2,:) = epi_crv(2,:)-LVUpperCenter(2);
            epi_crvT(3,:) = epi_crv(3,:)-LVUpperCenter(3);
            epi_crv = RotationMatrix*epi_crvT; 
            plot3(epi_crv(1,:),epi_crv(2,:), epi_crv(3,:),'LineStyle', '-', 'Color', ...
                        'r', 'LineWidth',2);
            end
            
            if ~isempty(epi_cpa)
            epi_cpaT(1,:) = epi_cpa(1,:)-LVUpperCenter(1);
            epi_cpaT(2,:) = epi_cpa(2,:)-LVUpperCenter(2);
            epi_cpaT(3,:) = epi_cpa(3,:)-LVUpperCenter(3);
            epi_cpa = RotationMatrix*epi_cpaT; 
            plot3(epi_cpa(1,:),epi_cpa(2,:), epi_cpa(3,:),'LineStyle', '-', 'Color', ...
                        'r', 'LineWidth',2);
            end
            
            DataSegSARotated(imIndex).endo_lvReal = endo_lv;
            DataSegSARotated(imIndex).endo_rvReal = endo_rv;
            DataSegSARotated(imIndex).endo_paReal = endo_pa;
            DataSegSARotated(imIndex).epi_cReal = epi_c;
            DataSegSARotated(imIndex).epi_crvReal = epi_crv;
            DataSegSARotated(imIndex).epi_cpaReal = epi_cpa;


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
                SolidWorksInputForClosedSplineGeneration(fidSA,endo_lv, 0, lined, 1);
                SolidWorksInputForClosedSplineGeneration(fidSA,endo_rv, 0, lined, 1);
                SolidWorksInputForClosedSplineGeneration(fidSA,epi_c, 0, lined, 0);
                fclose(fidSA);
            end

            clear endo_lvT;
            clear endo_rvT;
            clear endo_paT;
            clear epi_cT;
            clear epi_crvT;
            clear epi_cpaT;
        end

    end


    for imIndex = 1 : size(DataSegLA, 2)
        endo_lv = DataSegLA(imIndex).endo_lvReal;
        endo_rv = DataSegLA(imIndex).endo_rvReal;
        epi_c = DataSegLA(imIndex).epi_cReal;
        
        if ~isempty(endo_lv)&& ~isempty(endo_rv) && ~isempty(epi_c)
            
        endo_lvT(1,:) = endo_lv(1,:)-LVUpperCenter(1);
        endo_lvT(2,:) = endo_lv(2,:)-LVUpperCenter(2);
        endo_lvT(3,:) = endo_lv(3,:)-LVUpperCenter(3);
        endo_lv = RotationMatrix*endo_lvT;

        endo_rvT(1,:) = endo_rv(1,:)-LVUpperCenter(1);
        endo_rvT(2,:) = endo_rv(2,:)-LVUpperCenter(2);
        endo_rvT(3,:) = endo_rv(3,:)-LVUpperCenter(3);
        endo_rv = RotationMatrix*endo_rvT;
        
        epi_cT(1,:) = epi_c(1,:)-LVUpperCenter(1);
        epi_cT(2,:) = epi_c(2,:)-LVUpperCenter(2);
        epi_cT(3,:) = epi_c(3,:)-LVUpperCenter(3);
        epi_c = RotationMatrix*epi_cT; 

        DataSegLARotated(imIndex).endo_lvReal = endo_lv;
        DataSegLARotated(imIndex).endo_rvReal = endo_rv;
        DataSegLARotated(imIndex).epi_cReal = epi_c;

        figure(h3D_rotated);hold on;
        plot3(endo_lv(1,:),endo_lv(2,:), endo_lv(3,:),'LineStyle', '-', 'Color', ...
                    'b', 'LineWidth',2);
        plot3(endo_rv(1,:),endo_rv(2,:), endo_rv(3,:),'LineStyle', '-', 'Color', ...
                    'y', 'LineWidth',2);
        plot3(epi_c(1,:),epi_c(2,:), epi_c(3,:),'LineStyle', '-', 'Color', ...
                    'r', 'LineWidth',2);

        if Bwritten
            %%%output the bc boundaries
            cd(phase_resultDir);
            filename = sprintf('LABCdata_slice%d.txt',imIndex);
            fidSA = fopen(filename,'w');
            cd(workingDir);
            SolidWorksInputForClosedSplineGeneration(fidSA,endo_lv, 0, 1, 1);
            SolidWorksInputForClosedSplineGeneration(fidSA,endo_rv, 0, 1, 1);
            SolidWorksInputForClosedSplineGeneration(fidSA,epi_c, 0, 1, 0);
            fclose(fidSA);
        end

        clear endo_lvT;
        clear endo_rvT;
        clear epi_cT;
        end

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

%% write out for Geomagica
divd=100;
figure
hold on
endo_lv_data=[];endo_rv_data=[];endo_pa_data=[];
epi_c_data=[];epi_crv_data=[];epi_cpa_data=[];
theta_new = linspace(-pi, pi-0.01, divd);

for imIndex = 1 : size(DataSegSARotated, 2)
            
    endo_lv = DataSegSARotated(imIndex).endo_lvReal;
    endo_rv = DataSegSARotated(imIndex).endo_rvReal;
    endo_pa = DataSegSARotated(imIndex).endo_paReal;
    epi_c = DataSegSARotated(imIndex).epi_cReal;
    epi_crv = DataSegSARotated(imIndex).epi_crvReal;
    epi_cpa = DataSegSARotated(imIndex).epi_cpaReal;
    
    endo_lv_int=[];endo_rv_int=[];endo_pa_int=[];
    epi_c_int=[];epi_crv_int=[];epi_cpa_int=[];
    if ~isempty(endo_lv)
        %%%%% lv
        points=[];theta=[];rho=[];z=[];
        points = fnplt(cscvn(endo_lv));
        %[theta,rho,z] = cart2pol(points(1,:),points(2,:),points(3,:));
        center=mean(points,2);
        for i=1:size(points,2)
            vec=points(1:2,i)-center(1:2);
            theta(i)=atan2(vec(2),vec(1));
            rho(i)=norm(vec);
        end
   
        [theta_uni,m1,n1] = unique(theta,'stable');
        rho_uni=rho(m1);
        
        r_new = interp1(theta_uni, rho_uni, theta_new, 'linear', 'extrap');
  
        for i=1:length(theta_new)
            x(i)=r_new(i)*cos(theta_new(i))+center(1);
            y(i)=r_new(i)*sin(theta_new(i))+center(2);
        end
        
        z_new(1:length(theta_new))=center(3);
        endo_lv_int=[x; y; z_new];
        endo_lv_data=[endo_lv_data; endo_lv_int];
        
                
        plot3(endo_lv_int(1,:),endo_lv_int(2,:),endo_lv_int(3,:),'LineStyle', '-', 'Color', ...
                            'b', 'LineWidth',2)

    end
    if ~isempty(endo_rv)
        %%%%% rv
        points=[];theta=[];rho=[];z=[];
        points = fnplt(cscvn(endo_rv));
        %[theta,rho,z] = cart2pol(points(1,:),points(2,:),points(3,:));
        center=mean(points,2);
        for i=1:size(points,2)
            vec=points(1:2,i)-center(1:2);
            theta(i)=atan2(vec(2),vec(1));
            rho(i)=norm(vec);
        end
   
        [theta_uni,m1,n1] = unique(theta,'stable');
        rho_uni=rho(m1);
        
        r_new = interp1(theta_uni, rho_uni, theta_new, 'linear', 'extrap');
        
        for i=1:length(theta_new)
            x(i)=r_new(i)*cos(theta_new(i))+center(1);
            y(i)=r_new(i)*sin(theta_new(i))+center(2);
        end
        
        z_new(1:length(theta_new))=center(3);
        endo_rv_int=[x; y; z_new];
        endo_rv_data=[endo_rv_data; endo_rv_int];
        
        plot3(endo_rv_int(1,:),endo_rv_int(2,:),endo_rv_int(3,:),'LineStyle', '-', 'Color', ...
                            'g', 'LineWidth',2)
    end
    if ~isempty(endo_pa)
        %%%%% rv
        points=[];theta=[];rho=[];z=[];
        points = fnplt(cscvn(endo_pa));
        %[theta,rho,z] = cart2pol(points(1,:),points(2,:),points(3,:));
        center=mean(points,2);
        for i=1:size(points,2)
            vec=points(1:2,i)-center(1:2);
            theta(i)=atan2(vec(2),vec(1));
            rho(i)=norm(vec);
        end
   
        [theta_uni,m1,n1] = unique(theta,'stable');
        rho_uni=rho(m1);
        
        r_new = interp1(theta_uni, rho_uni, theta_new, 'linear', 'extrap');
        
        for i=1:length(theta_new)
            x(i)=r_new(i)*cos(theta_new(i))+center(1);
            y(i)=r_new(i)*sin(theta_new(i))+center(2);
        end
        
        z_new(1:length(theta_new))=center(3);
        endo_pa_int=[x; y; z_new];
        endo_pa_data=[endo_pa_data; endo_pa_int];
        
        plot3(endo_pa_int(1,:),endo_pa_int(2,:),endo_pa_int(3,:),'LineStyle', '-', 'Color', ...
                            'g', 'LineWidth',2)
    end   
    
    if ~isempty(epi_c)
        
        %%%%%% epi
        points=[];theta=[];rho=[];z=[];
        points = fnplt(cscvn(epi_c));
        %[theta,rho,z] = cart2pol(points(1,:),points(2,:),points(3,:));
        center=mean(points,2);
        for i=1:size(points,2)
            vec=points(1:2,i)-center(1:2);
            theta(i)=atan2(vec(2),vec(1));
            rho(i)=norm(vec);
        end
   
        [theta_uni,m1,n1] = unique(theta,'stable');
        rho_uni=rho(m1);
        
        r_new = interp1(theta_uni, rho_uni, theta_new, 'linear', 'extrap');
        for i=1:length(theta_new)
            x(i)=r_new(i)*cos(theta_new(i))+center(1);
            y(i)=r_new(i)*sin(theta_new(i))+center(2);
        end
        
        z_new(1:length(theta_new))=center(3);
        epi_c_int=[x; y; z_new];
        epi_c_data=[epi_c_data; epi_c_int];
        
        plot3(epi_c_int(1,:),epi_c_int(2,:),epi_c_int(3,:),'LineStyle', '-', 'Color', ...
                            'r', 'LineWidth',2)

    end
    
    if ~isempty(epi_crv)
        
        %%%%%% epi
        points=[];theta=[];rho=[];z=[];
        points = fnplt(cscvn(epi_crv));
        %[theta,rho,z] = cart2pol(points(1,:),points(2,:),points(3,:));
        center=mean(points,2);
        for i=1:size(points,2)
            vec=points(1:2,i)-center(1:2);
            theta(i)=atan2(vec(2),vec(1));
            rho(i)=norm(vec);
        end
   
        [theta_uni,m1,n1] = unique(theta,'stable');
        rho_uni=rho(m1);
        
        r_new = interp1(theta_uni, rho_uni, theta_new, 'linear', 'extrap');
        for i=1:length(theta_new)
            x(i)=r_new(i)*cos(theta_new(i))+center(1);
            y(i)=r_new(i)*sin(theta_new(i))+center(2);
        end
        
        z_new(1:length(theta_new))=center(3);
        epi_crv_int=[x; y; z_new];
        epi_crv_data=[epi_crv_data; epi_crv_int];
        
        plot3(epi_crv_int(1,:),epi_crv_int(2,:),epi_crv_int(3,:),'LineStyle', '-', 'Color', ...
                            'r', 'LineWidth',2)

     end
    
    if ~isempty(epi_cpa)
        
        %%%%%% epi
        points=[];theta=[];rho=[];z=[];
        points = fnplt(cscvn(epi_cpa));
        %[theta,rho,z] = cart2pol(points(1,:),points(2,:),points(3,:));
        center=mean(points,2);
        for i=1:size(points,2)
            vec=points(1:2,i)-center(1:2);
            theta(i)=atan2(vec(2),vec(1));
            rho(i)=norm(vec);
        end
   
        [theta_uni,m1,n1] = unique(theta,'stable');
        rho_uni=rho(m1);
        
        r_new = interp1(theta_uni, rho_uni, theta_new, 'linear', 'extrap');
        for i=1:length(theta_new)
            x(i)=r_new(i)*cos(theta_new(i))+center(1);
            y(i)=r_new(i)*sin(theta_new(i))+center(2);
        end
        
        z_new(1:length(theta_new))=center(3);
        epi_cpa_int=[x; y; z_new];
        epi_cpa_data=[epi_cpa_data; epi_cpa_int];
        
        plot3(epi_cpa_int(1,:),epi_cpa_int(2,:),epi_cpa_int(3,:),'LineStyle', '-', 'Color', ...
                            'r', 'LineWidth',2)

    end
    
        
end  


%%%%%%interpolation along z axis
figure
hold on
%%%% LV
min_z=min(DataSegLARotated(2).endo_lvReal(3,:));
max_z=max(endo_lv_data(3:3:size(endo_lv_data,1),1));  %max(DataSegLARotated(2).endo_lvReal(3,:)); %1.5;%
apex_n=find(DataSegLARotated(2).endo_lvReal(3,:)==min_z);

z_int= linspace(min_z, max_z, divd/2);
LV_new=[];
for i=1:divd
    new_data=[];
    x=[];y=[];z=[];
    x=endo_lv_data(1:3:size(endo_lv_data),i);
    y=endo_lv_data(2:3:size(endo_lv_data),i);
    z=endo_lv_data(3:3:size(endo_lv_data),i);
    x_len=length(x);
    x(x_len+1)=DataSegLARotated(2).endo_lvReal(1,apex_n);
    y(x_len+1)=DataSegLARotated(2).endo_lvReal(2,apex_n);
    z(x_len+1)=DataSegLARotated(2).endo_lvReal(3,apex_n);
    
    xx = spline(z,x,z_int);
    yy = spline(z,y,z_int);
    
    new_data=[xx; yy; z_int];
    LV_new=[LV_new new_data];
    
    plot3(new_data(1,:),new_data(2,:),new_data(3,:),'LineStyle', '-', 'Color', ...
                            'b', 'LineWidth',2)
    
end

%%%% RV
min_z=min(DataSegLARotated(2).endo_rvReal(3,:));
max_z=(endo_lv_data(3,1));%max(DataSegLARotated(2).endo_rvReal(3,:));
closest_1 = interp1(endo_rv_data(:,1),endo_rv_data(:,1),max_z,'nearest');
closest_2=find(endo_rv_data(:,1)==closest_1);
apex_n=find(DataSegLARotated(2).endo_rvReal(3,:)==min_z);
z_int= linspace(min_z, max_z, divd/2);
RV_new=[];

for i=1:divd
    new_data=[];
    x=[];y=[];z=[];
    x=endo_rv_data(closest_2-2:3:size(endo_rv_data),i);
    y=endo_rv_data(closest_2-1:3:size(endo_rv_data),i);
    z=endo_rv_data(closest_2:3:size(endo_rv_data),i);
    x_len=length(x);
    x(x_len+1)=DataSegLARotated(2).endo_rvReal(1,apex_n);
    y(x_len+1)=DataSegLARotated(2).endo_rvReal(2,apex_n);
    z(x_len+1)=DataSegLARotated(2).endo_rvReal(3,apex_n);
    
    xx = spline(z,x,z_int);
    yy = spline(z,y,z_int);
    
    new_data=[xx; yy; z_int];
    RV_new=[RV_new new_data];
    
    plot3(new_data(1,:),new_data(2,:),new_data(3,:),'LineStyle', '-', 'Color', ...
                            'g', 'LineWidth',2)
end

%%%% RV-top
min_z=(endo_lv_data(3,1)); %max(endo_lv_data(3:3:size(endo_lv_data,1),1));
max_z=max(endo_rv_data(3:3:size(endo_rv_data,1),1));%max(DataSegLARotated(2).endo_rvReal(3,:));
z_int=min_z:0.4:max_z; %linspace(min_z, max_z, divd/4);
RV_new_top=[];
for i=1:divd
    new_data=[];
    x=[];y=[];z=[];
    x=endo_rv_data(1:3:closest_2,i);
    y=endo_rv_data(2:3:closest_2,i);
    z=endo_rv_data(3:3:closest_2,i);

    xx = spline(z,x,z_int);
    yy = spline(z,y,z_int);
    
    new_data=[xx; yy; z_int];
    %RV_new=[RV_new new_data];
    RV_new_top=[RV_new_top new_data];
    
   % plot3(new_data(1,:),new_data(2,:),new_data(3,:),'LineStyle', '-', 'Color', ...
    %                        'b', 'LineWidth',2)
end


%%%% PA-top
min_z=(endo_lv_data(3,1)); %max(endo_lv_data(3:3:size(endo_lv_data,1),1));
max_z=max(endo_pa_data(3:3:size(endo_pa_data,1),1));%max(DataSegLARotated(2).endo_rvReal(3,:));
z_int= min_z:0.4:max_z;%linspace(min_z, max_z, divd/4);
PA_new_top=[];
for i=1:divd
    new_data=[];
    x=[];y=[];z=[];
    x=[endo_pa_data(1:3:size(endo_pa_data),i);  endo_rv_data(closest_2-2,i)];
    y=[endo_pa_data(2:3:size(endo_pa_data),i);  endo_rv_data(closest_2-1,i)];
    z=[endo_pa_data(3:3:size(endo_pa_data),i);  endo_rv_data(closest_2,i)];

    xx = spline(z,x,z_int);
    yy = spline(z,y,z_int);
    
    new_data=[xx; yy; z_int];
    %RV_new=[RV_new new_data];
    PA_new_top=[PA_new_top new_data];
    
    %plot3(new_data(1,:),new_data(2,:),new_data(3,:),'LineStyle', '-', 'Color', ...
    %                        'r', 'LineWidth',2)
end

%%%%%%%%%%%% Delete overlapping nodes
RV_TOP=[PA_new_top RV_new_top];

idxRV=[];
[~,idxRV] = sort(RV_TOP(3,:)); % sort just the first column
RV_TOP = RV_TOP(:,idxRV);
z_order=unique(RV_TOP(3,:));

for i=1:length(z_order)
    %RV
    columes1=find(RV_new_top(3,:)==z_order(i));
    xyplane1=RV_new_top(:,columes1);
    %PA
    columes2=find(PA_new_top(3,:)==z_order(i));
    xyplane2=PA_new_top(:,columes2);
    %
    poly1 = polyshape(xyplane1(1:2,:)');
    poly2 = polyshape(xyplane2(1:2,:)');   
    poly3 = intersect(poly1,poly2);
    overlap_nodes=poly3.Vertices;
    all_nodes=[poly1.Vertices; poly2.Vertices;];
    if ~isempty(overlap_nodes) 
        RowIdx = find(ismember(all_nodes, overlap_nodes,'rows'));
        all_nodes(RowIdx,:)=[]; 
    end
    all_nodes(:,3)=z_order(i);
    RV_new=[RV_new all_nodes'];
    
    plot3(all_nodes(:,1),all_nodes(:,2),all_nodes(:,3),'LineStyle', '-', 'Color', ...
                            'g', 'LineWidth',2)
end



%%%% EPI
min_z=min(DataSegLARotated(2).epi_cReal(3,:));
max_z=max(epi_c_data(3:3:size(epi_c_data,1),1));%max(DataSegLARotated(2).epi_cReal(3,:));
apex_n=find(DataSegLARotated(2).epi_cReal(3,:)==min_z);
z_int= linspace(min_z, max_z, divd/2);
EPI_new=[];
for i=1:divd
    new_data=[];
    x=[];y=[];z=[];
    x=epi_c_data(1:3:size(epi_c_data),i);
    y=epi_c_data(2:3:size(epi_c_data),i);
    z=epi_c_data(3:3:size(epi_c_data),i);
    x_len=length(x);
    x(x_len+1)=DataSegLARotated(2).epi_cReal(1,apex_n);
    y(x_len+1)=DataSegLARotated(2).epi_cReal(2,apex_n);
    z(x_len+1)=DataSegLARotated(2).epi_cReal(3,apex_n);
    
    xx = spline(z,x,z_int);
    yy = spline(z,y,z_int);
    
    new_data=[xx; yy; z_int];
    EPI_new=[EPI_new new_data];
    
    plot3(new_data(1,:),new_data(2,:),new_data(3,:),'LineStyle', '-', 'Color', ...
                            'r', 'LineWidth',2)
end

%%%% RV-EPI-top
min_z=max(epi_c_data(3:3:size(epi_c_data,1),1));
max_z=max(epi_crv_data(3:3:size(epi_crv_data,1),1));%max(DataSegLARotated(2).epi_cReal(3,:));
z_int= min_z:0.4:max_z; linspace(min_z, max_z, divd/2);
RV_EPI_top=[];
for i=1:divd
    new_data=[];
    x=[];y=[];z=[];
    x=[epi_crv_data(1:3:size(epi_crv_data),i); epi_c_data(1,i)]; 
    y=[epi_crv_data(2:3:size(epi_crv_data),i); epi_c_data(2,i)];
    z=[epi_crv_data(3:3:size(epi_crv_data),i); epi_c_data(3,i)];
    
    xx = spline(z,x,z_int);
    yy = spline(z,y,z_int);
    
    new_data=[xx; yy; z_int];
    %EPI_new=[EPI_new new_data];
    RV_EPI_top=[RV_EPI_top new_data];
    
    %plot3(new_data(1,:),new_data(2,:),new_data(3,:),'LineStyle', '-', 'Color', ...
     %                       'r', 'LineWidth',2)
end

%%%% PA-EPI-top
min_z=max(epi_c_data(3:3:size(epi_c_data,1),1));
max_z=max(epi_cpa_data(3:3:size(epi_cpa_data,1),1));%max(DataSegLARotated(2).epi_cReal(3,:));
z_int= min_z:0.4:max_z;%linspace(min_z, max_z, divd/2);
PA_EPI_top=[];
for i=1:divd
    new_data=[];
    x=[];y=[];z=[];
    x=[epi_cpa_data(1:3:size(epi_cpa_data),i); epi_c_data(1,i)]; 
    y=[epi_cpa_data(2:3:size(epi_cpa_data),i); epi_c_data(2,i)]; 
    z=[epi_cpa_data(3:3:size(epi_cpa_data),i); epi_c_data(3,i)]; 
    
    xx = spline(z,x,z_int);
    yy = spline(z,y,z_int);
    
    new_data=[xx; yy; z_int];
    %EPI_new=[EPI_new new_data];
    PA_EPI_top=[PA_EPI_top new_data];
    
    %plot3(new_data(1,:),new_data(2,:),new_data(3,:),'LineStyle', '-', 'Color', ...
       %                     'r', 'LineWidth',2)
end


%%%%%%%%%%%% Delete overlapping nodes
EPI_TOP=[PA_EPI_top RV_EPI_top];

idxEPI=[];
[~,idxEPI] = sort(EPI_TOP(3,:)); % sort just the first column
EPI_TOP = EPI_TOP(:,idxEPI);
z_order=unique(EPI_TOP(3,:));

for i=1:length(z_order)
    %RV
    columes1=find(RV_EPI_top(3,:)==z_order(i));
    xyplane1=RV_EPI_top(:,columes1);
    %PA
    columes2=find(PA_EPI_top(3,:)==z_order(i));
    xyplane2=PA_EPI_top(:,columes2);
    %
    poly1 = polyshape(xyplane1(1:2,:)');
    poly2 = polyshape(xyplane2(1:2,:)');   
    poly3 = intersect(poly1,poly2);
    overlap_nodes=poly3.Vertices;
    all_nodes=[poly1.Vertices; poly2.Vertices;];
    if ~isempty(overlap_nodes) 
        RowIdx = find(ismember(all_nodes, overlap_nodes,'rows'));
        all_nodes(RowIdx,:)=[]; 
    end
    all_nodes(:,3)=z_order(i);
    EPI_new=[EPI_new all_nodes'];
    
    plot3(all_nodes(:,1),all_nodes(:,2),all_nodes(:,3),'LineStyle', '-', 'Color', ...
                            'r', 'LineWidth',2)
end


        
if Bwritten
                    %%%output the bc boundaries
   cd(phase_resultDir);
   filename = sprintf('SA_LV.txt');
   fidSA = fopen(filename,'w');
   for i=1:size(LV_new,2)
       fprintf(fidSA,'%14.10f \t%14.10f \t%14.10f\n',LV_new(1,i),LV_new(2,i),LV_new(3,i));
   end
   fclose(fidSA);
   
   filename = sprintf('SA_RV.txt');
   fidSA = fopen(filename,'w');
   for i=1:size(RV_new,2)
       fprintf(fidSA,'%14.10f \t%14.10f \t%14.10f\n',RV_new(1,i),RV_new(2,i),RV_new(3,i));
   end
   fclose(fidSA);
   
   filename = sprintf('SA_EPI.txt');
   fidSA = fopen(filename,'w');
   for i=1:size(EPI_new,2)
       fprintf(fidSA,'%14.10f \t%14.10f \t%14.10f\n',EPI_new(1,i),EPI_new(2,i),EPI_new(3,i));
   end
   fclose(fidSA);
   
end







