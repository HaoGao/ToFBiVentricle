%% SA segmentation
%% updated on 15th May 2023, based on Debao Guan's modification
% 4 choices, all results are saved in DataSegSA data structure
% option 1: biventricle boundaries: LV+RV, saved in endo_lv, endo_rv and epi_c
% option 2: LV+RV+PA, saved in endo_lv, endo_rv, endo_pa, epi_c, epi_crv, epi_cpa
% option 3: RV + RA, saved in endo_rv, epi_crv, endo_pa, endo_rv
% option 4: PA, saved in in endo_pa, endo_rv

clear all; close all; clc;

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
                endo_pa=[];
                epi_biv=[];%epi boundaries for LV and RV or for LV epi
                epi_crv=[];%epi boundaries for RV
                epi_cpa=[];%epi boundaries for PA

                endo_sample_lv = [];
                endo_sample_rv = [];
                endo_sample_pa = [];
                epi_sample = [];
                epi_sample_rv = [];
                epi_sample_pa = [];

                list_phaseS = {'Case 1: LV-RV-EPI', 'Case 2: LV-EPI-RV-EPI-PA-EPI', ...
                    'Case 3: RV-EPI-PA-EPI','Case 4: PA-EPI'};
                [idxs, tfs] = listdlg('PromptString',{'Select a case.', ...
                    'Only one file can be selected at a time.',...
                    'Segmentation in order'}, 'ListString', list_phaseS);

                switch idxs
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    case 1
                        disp('.............................................');
                        disp('predefined boundaries are plotted in red');
                        disp('define points in the boundary by click points')
                        disp('press stop, then click the last point');
                        disp('press replot to generate the curve');
                        disp('whenever change the point position, using replot to regenerate the curve');
                        disp('double press space key to return from boudnary definition');


                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LV ENDO
                        choice = questdlg('LV-ENDO: Would you like to START?', ...
                            'Dessert Menu', 'Yes','No','Yes');
                        if strcmp(choice, 'Yes')
                            [endo_lv, hpts_endo_lv]=defineBoundaryByimpoint_2023(imCropData, [], 1, 'LV endo');

                            endo_sample_lv(1,:) = endo_lv(1,:) +rect(1);
                            endo_sample_lv(2,:) = endo_lv(2,:) +rect(2);
                        end

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  RV ENDO
                        choice = questdlg('RV_ENDO: Would you like to START?', ...
                            'Dessert Menu', 'Yes','No','Yes');
                        if strcmp(choice, 'Yes')
                            [endo_rv, hpts_endo_rv]=defineBoundaryByimpoint_2023(imCropData, endo_lv, 1, 'RV endo');

                            endo_sample_rv(1,:) = endo_rv(1,:) +rect(1);
                            endo_sample_rv(2,:) = endo_rv(2,:) +rect(2);
                        end

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BIV-EPI
                        choice = questdlg('LV_RV_EPI: Would you like to START?', ...
                            'Dessert Menu', 'Yes','No','Yes');
                        if strcmp(choice, 'Yes')
                            [epi_biv, hpts_epi_biv]=defineBoundaryByimpoint_2023(imCropData, [endo_lv endo_rv], 1, 'LV_RV_EPI');

                            epi_sample(1,:) = epi_biv(1,:) +rect(1);
                            epi_sample(2,:) = epi_biv(2,:) +rect(2);
                        end


                        %                         % load the previous boundaries
                        %                         basalSliceIndex = patientConfigs(patientIndex,1).basalSliceIndex;
                        %                         if imIndex > basalSliceIndex % starting from the second slice
                        %                              [endo_lv_p, endo_rv_p, epi_biv_p] = project_previousSliceBC(imIndex-1, DataSegSA, rect);
                        %                         end
                        % 					    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LV ENDO
                        %                         choice = questdlg('LV-ENDO: Would you like to START?', ...
                        %                             'Dessert Menu', 'Yes','No','Yes');
                        %
                        %                         Bool_manual_correction = 0;
                        %                         choice1 = questdlg('Make manual corrections', ...
                        %                             'Dessert Menu', 'Yes','No','Yes');
                        %                         if strcmp(choice1, 'Yes')
                        %                             Bool_manual_correction = 1;
                        %                         end
                        %
                        %                         if Bool_manual_correction
                        %                           [endo_lv, hpts_endo_lv] = defineBoundaryByimpointWithPreviousBoundary_2023(imCropData, [], 1, 'LV endo',endo_lv_p);
                        %                         else
                        %                            [endo_lv, hpts_endo_lv]=defineBoundaryByimpoint_2023(imCropData, [], 1, 'LV endo');
                        %                         end
                        %
                        % 				        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  RV ENDO
                        %                         choice = questdlg('RV_ENDO: Would you like to START?', ...
                        %                             'Dessert Menu', 'Yes','No','Yes');
                        %                         if Bool_manual_correction
                        %                             [endo_rv, hpts_endo_rv]=defineBoundaryByimpointWithPreviousBoundary_2023(imCropData, endo_lv, 1, 'RV endo', endo_rv_p);
                        %                         else
                        %                             [endo_rv, hpts_endo_rv]=defineBoundaryByimpoint_2023(imCropData, endo_lv, 1, 'RV endo');
                        %                         end
                        %
                        % 						%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BIV-EPI
                        %                         choice = questdlg('LV_RV_EPI: Would you like to START?', ...
                        %                             'Dessert Menu', 'Yes','No','Yes');
                        %                         if Bool_manual_correction
                        %                             [epi_biv, hpts_epi_biv]=defineBoundaryByimpointWithPreviousBoundary_2023(imCropData, [endo_lv endo_rv], 1, 'LV_RV_EPI', epi_biv_p);
                        %                         else
                        %                             [epi_biv, hpts_epi_biv]=defineBoundaryByimpoint_2023(imCropData, [endo_lv endo_rv], 1, 'LV_RV_EPI');
                        %                         end
                        %




                        %endo_sample_pa(1,:) = endo_pa(1,:) +rect(1);
                        %endo_sample_pa(2,:) = endo_pa(2,:) +rect(2);

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    case 2 % this case will may not often to see
                        % RV endo, PA endo, RV+PA endo
                        disp('.............................................');
                        disp('predefined boundaries are plotted in red');
                        disp('define points in the boundary by click points')
                        disp('press stop, then click the last point');
                        disp('press replot to generate the curve');
                        disp('whenever change the point position, using replot to regenerate the curve');
                        disp('double press space key to return from boudnary definition');

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LV ENDO
                        choice = questdlg('LV-ENDO: Would you like to START?', ...
                            'Dessert Menu', 'Yes','No','Yes');
                        if strcmp(choice, 'Yes')
                            [endo_lv, hpts_endo_lv]=defineBoundaryByimpoint_2023(imCropData, [], 1, 'LV endo');

                            endo_sample_lv(1,:) = endo_lv(1,:) +rect(1);
                            endo_sample_lv(2,:) = endo_lv(2,:) +rect(2);
                        end

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LV EPI
                        choice = questdlg('LV_EPI: Would you like to START?', ...
                            'Dessert Menu', 'Yes','No','Yes');
                        if strcmp(choice, 'Yes')
                            [epi_biv, hpts_epi_biv]=defineBoundaryByimpoint_2023(imCropData, endo_lv, 1, 'LV EPI');

                            epi_sample(1,:) = epi_biv(1,:) +rect(1);
                            epi_sample(2,:) = epi_biv(2,:) +rect(2);
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RV ENDO
                        choice = questdlg('RV_ENDO: Would you like to START?', ...
                            'Dessert Menu', 'Yes','No','Yes');
                        if strcmp(choice, 'Yes')
                            [endo_rv, hpts_endo_rv]=defineBoundaryByimpoint_2023(imCropData, [endo_lv epi_biv], 1, 'RV endo');

                            endo_sample_rv(1,:) = endo_rv(1,:) +rect(1);
                            endo_sample_rv(2,:) = endo_rv(2,:) +rect(2);
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RV EPI
                        choice = questdlg('RV_EPI: Would you like to START?', ...
                            'Dessert Menu', 'Yes','No','Yes');
                        if strcmp(choice, 'Yes')
                            [epi_rv, hpts_epi_rv]=defineBoundaryByimpoint_2023(imCropData, [endo_lv epi_biv endo_rv], 1, 'RV EPI');

                            epi_sample_rv(1,:) = epi_rv(1,:) +rect(1);
                            epi_sample_rv(2,:) = epi_rv(2,:) +rect(2);
                        end

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PA ENDO
                        choice = questdlg('PA_ENDO: Would you like to START?', ...
                            'Dessert Menu', 'Yes','No','Yes');
                        if strcmp(choice, 'Yes')
                            [endo_pa, hpts_endo_pa]=defineBoundaryByimpoint_2023(imCropData, [endo_lv epi_biv endo_rv epi_rv], 1, 'PA endo');

                            endo_sample_pa(1,:) = endo_pa(1,:) +rect(1);
                            endo_sample_pa(2,:) = endo_pa(2,:) +rect(2);
                        end

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PA EPI
                        choice = questdlg('PA_EPI: Would you like to START?', ...
                            'Dessert Menu', 'Yes','No','Yes');
                        if strcmp(choice, 'Yes')
                            [epi_pa, hpts_epi_pa]=defineBoundaryByimpoint_2023(imCropData, [endo_lv epi_biv endo_rv epi_rv endo_pa], 1, 'PA EPI');

                            epi_sample_pa(1,:) = epi_pa(1,:) +rect(1);
                            epi_sample_pa(2,:) = epi_pa(2,:) +rect(2);
                        end

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    case 3
                        % RV endo, RV epi, PA endo, PA epi
                        disp('.............................................');
                        disp('predefined boundaries are plotted in red');
                        disp('define points in the boundary by click points')
                        disp('press stop, then click the last point');
                        disp('press replot to generate the curve');
                        disp('whenever change the point position, using replot to regenerate the curve');
                        disp('double press space key to return from boudnary definition');

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RV ENDO
                        choice = questdlg('RV_ENDO: Would you like to START?', ...
                            'Dessert Menu', 'Yes','No','Yes');
                        if strcmp(choice, 'Yes')
                            [endo_rv, hpts_endo_rv]=defineBoundaryByimpoint_2023(imCropData, [], 1, 'RV endo');

                            endo_sample_rv(1,:) = endo_rv(1,:) +rect(1);
                            endo_sample_rv(2,:) = endo_rv(2,:) +rect(2);
                        end

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RV EPI
                        choice = questdlg('RV_EPI: Would you like to START?', ...
                            'Dessert Menu', 'Yes','No','Yes');
                        if strcmp(choice, 'Yes')
                            [epi_rv, hpts_epi_rv]=defineBoundaryByimpoint_2023(imCropData, endo_rv, 1, 'RV epi');

                            epi_sample_rv(1,:) = epi_rv(1,:) +rect(1);
                            epi_sample_rv(2,:) = epi_rv(2,:) +rect(2);
                        end

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PA ENDO
                        choice = questdlg('PA_ENDO: Would you like to START?', ...
                            'Dessert Menu', 'Yes','No','Yes');
                        if strcmp(choice, 'Yes')
                            [endo_pa, hpts_endo_pa]=defineBoundaryByimpoint_2023(imCropData, [endo_rv epi_rv], 1, 'PA endo');

                            endo_sample_pa(1,:) = endo_pa(1,:) +rect(1);
                            endo_sample_pa(2,:) = endo_pa(2,:) +rect(2);
                        end

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% PA EPI
                        choice = questdlg('PA_EPI: Would you like to START?', ...
                            'Dessert Menu', 'Yes','No','Yes');
                        if strcmp(choice, 'Yes')
                            [epi_pa, hpts_epi_pa]=defineBoundaryByimpoint_2023(imCropData, endo_pa, 1, 'PA epi');

                            epi_sample_pa(1,:) = epi_pa(1,:) +rect(1);
                            epi_sample_pa(2,:) = epi_pa(2,:) +rect(2);
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    case 4
                        % PA endo; PA epi
                        disp('.............................................');
                        disp('predefined boundaries are plotted in red');
                        disp('define points in the boundary by click points')
                        disp('press stop, then click the last point');
                        disp('press replot to generate the curve');
                        disp('whenever change the point position, using replot to regenerate the curve');
                        disp('double press space key to return from boudnary definition');

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PA ENDO
                        choice = questdlg('PA_ENDO: Would you like to START?', ...
                            'Dessert Menu', 'Yes','No','Yes');
                        if strcmp(choice, 'Yes')
                            [endo_pa, hpts_endo_pa]=defineBoundaryByimpoint_2023(imCropData, [], 1, 'PA endo');

                            endo_sample_pa(1,:) = endo_pa(1,:) +rect(1);
                            endo_sample_pa(2,:) = endo_pa(2,:) +rect(2);
                        end

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PA EPI
                        choice = questdlg('PA_EPI: Would you like to START?', ...
                            'Dessert Menu', 'Yes','No','Yes');
                        if strcmp(choice, 'Yes')
                            [epi_pa, hpts_epi_pa]=defineBoundaryByimpoint_2023(imCropData, endo_pa, 1, 'PA epi');

                            epi_sample_pa(1,:) = epi_pa(1,:) +rect(1);
                            epi_sample_pa(2,:) = epi_pa(2,:) +rect(2);
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    otherwise
                        disp('Wrong Input, Try Again')
                end


                h1=figure();
                imshow(imCropData,[]); hold on;
                title(sliceLocationStr);
                if ~isempty(endo_sample_lv)
                    bc_interp = samplingBCWithoutIm_interp(endo_sample_lv, 1);
                    plot(bc_interp(1,:)-rect(1),bc_interp(2,:)-rect(2),'b-');
                end
                if ~isempty(endo_sample_rv)
                    bc_interp = samplingBCWithoutIm_interp(endo_sample_rv, 1);
                    plot(bc_interp(1,:)-rect(1),bc_interp(2,:)-rect(2),'b-');
                end
                if ~isempty(endo_sample_pa)
                    bc_interp = samplingBCWithoutIm_interp(endo_sample_pa, 1);
                    plot(bc_interp(1,:)-rect(1),bc_interp(2,:)-rect(2),'g-');
                end
                if ~isempty(epi_sample)
                    bc_interp = samplingBCWithoutIm_interp(epi_sample, 1);
                    plot(bc_interp(1,:)-rect(1),bc_interp(2,:)-rect(2),'r-');
                end
                if ~isempty(epi_sample_rv)
                    bc_interp = samplingBCWithoutIm_interp(epi_sample_rv, 1);
                    plot(bc_interp(1,:)-rect(1),bc_interp(2,:)-rect(2),'r-');
                end
                if ~isempty(epi_sample_pa)
                    bc_interp = samplingBCWithoutIm_interp(epi_sample_pa, 1);
                    plot(bc_interp(1,:)-rect(1),bc_interp(2,:)-rect(2),'r-');
                end

                %save the figure with boundaries
                cd(phase_resultDir);
                figure_name = sprintf('SA-%d.png',imIndex);
                print(h1,figure_name, '-dpng', '-r300');
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
                if ~isempty(endo_sample_pa)
                    endo_paReal = TransformCurvesFromImToRealSpace(endo_sample_pa,imInfo);
                else
                    endo_paReal = [];
                end
                if ~isempty(epi_sample)
                    epi_cReal = TransformCurvesFromImToRealSpace(epi_sample,imInfo);
                else
                    epi_cReal = [];
                end
                if ~isempty(epi_sample_rv)
                    epi_crvReal = TransformCurvesFromImToRealSpace(epi_sample_rv,imInfo);
                else
                    epi_crvReal = [];
                end
                if ~isempty(epi_sample_pa)
                    epi_cpaReal = TransformCurvesFromImToRealSpace(epi_sample_pa,imInfo);
                else
                    epi_cpaReal = [];
                end


                DataSegSA(imIndex).endo_lv = endo_sample_lv;
                DataSegSA(imIndex).endo_rv = endo_sample_rv;
                DataSegSA(imIndex).endo_pa = endo_sample_pa;
                DataSegSA(imIndex).epi_c = epi_sample;% epi_c either save the bi-ventricle epi or only the LV epi when epi_crv is not empty
                DataSegSA(imIndex).epi_crv = epi_sample_rv;
                DataSegSA(imIndex).epi_cpa = epi_sample_pa;
                DataSegSA(imIndex).endo_lvReal = endo_lvReal;
                DataSegSA(imIndex).endo_rvReal = endo_rvReal;
                DataSegSA(imIndex).endo_paReal = endo_paReal;
                DataSegSA(imIndex).epi_cReal = epi_cReal;
                DataSegSA(imIndex).epi_crvReal = epi_crvReal;
                DataSegSA(imIndex).epi_cpaReal = epi_cpaReal;


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
    endo_pa = DataSegSA(imIndex).endo_paReal;
    epi_c = DataSegSA(imIndex).epi_cReal;
    epi_crv = DataSegSA(imIndex).epi_crvReal;
    epi_cpa = DataSegSA(imIndex).epi_cpaReal;
    %%%%plot 3D curves
    if ~isempty(endo_lv)
        plot3(endo_lv(1,:),endo_lv(2,:), endo_lv(3,:),'LineStyle', '-', 'Color', 'b', 'LineWidth',2);
    end
    if ~isempty(endo_rv)
        plot3(endo_rv(1,:),endo_rv(2,:), endo_rv(3,:),'LineStyle', '-', 'Color', 'b', 'LineWidth',2);
    end
    if ~isempty(endo_pa)
        plot3(endo_pa(1,:),endo_pa(2,:), endo_pa(3,:),'LineStyle', '-', 'Color', 'g', 'LineWidth',2);
    end
    if ~isempty(epi_c)
        plot3(epi_c(1,:),epi_c(2,:), epi_c(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
    end
    if ~isempty(epi_crv)
        plot3(epi_crv(1,:),epi_crv(2,:), epi_crv(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
    end
    if ~isempty(epi_cpa)
        plot3(epi_cpa(1,:),epi_cpa(2,:), epi_cpa(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
    end
end

imgTecplotOutput(imData,'SAslice3Tec.dat',phase_resultDir);
%%%this is for apex info










