%%this will be used to read in an abaqus input file and convert into
%%different types of mesh files for future use
clear all;
close all;
clc;

% LVWM_config;
LVWM_config;

%%now need to choose phase to segment
list_phase = {'early_diastole', 'end_diastole', 'end_systole'};
[idx_phase, tf] = listdlg('ListString', list_phase);
cd(resultDir);
if ~exist(list_phase{idx_phase},'dir')
    mkdir(list_phase{idx_phase});
    cd(list_phase{idx_phase});
    phase_resultDir = pwd();    
else
    cd(list_phase{idx_phase});
    phase_resultDir = pwd();
end
cd(workingDir);
phase_selected = list_phase{idx_phase};

cd(phase_resultDir);
cd(solidworksDir);
solidworksDir = pwd();
cd(gmeshDir);
gmeshDir = pwd();
cd(workingDir);


list_steps = {'read gmesh file', 'read abaqus file for surface vtk mesh', ...
             'write abaqus input', 'write libmesh', 'seperate LV and RV'};
[idx, ~] = listdlg('ListString', list_steps);

if idx == 1
    %%first step will read the gmesh output file
    abaqusInput_gmesh = ReadAbaqusInput(abaqusInput_gmesh, gmeshDir);

    %% output the file into another new abaqus input file with natural ordering
    WriteAbaqusInput_Tet4(abaqusInput_gmesh, abaqusInput_matlab_gmesh, gmeshDir);
    %% this will generate AbaqusInputAdjustedNodeNumber.inp
end

if idx == 2
    %%now need to use abaqus to define node sets and surface sets, so that can
    %%be used to generate the surface mesh for deformetric
    %% after done with abaqus, and save the input file into gmesh with name AbaqusInputSets.inp
    msgbox('using abaqus to generate AbaqusInputSets.inp');

    %% now finally read the input file into matalb 
    surfaceList=['SURF_BASE', 'SURF_LV_ENDO', 'SURF_RV_ENDO', 'SURF_EPI', 'SURF_RV_SEPTUM_ENDO'];
    abaqusInput = ReadAbaqusInput(abaqusInput_abaqus_sets, gmeshDir, surfaceList);
    cd(gmeshDir);
    save abaqusInput abaqusInput;
    cd(workingDir);

    %% write out the vtk surface mesh file in PolyGon format
    bi_ven_surf_vtk_name = sprintf('bi_ven_surf_%s.vtk', phase_selected);
    bi_ven_vol_vtk_name = sprintf('bi_ven_vol_%s.vtk', phase_selected);
    val_elem = zeros([size(abaqusInput.surface_elements,1), 1]);
    write_vtk_trigular_surface(bi_ven_surf_vtk_name, gmeshDir, abaqusInput.surface_nodes, ...
                                    abaqusInput.surface_elements, val_elem);
    
                            
    val_elem = zeros([size(abaqusInput.elems,1), 1]);
    write_vtk_tet_volume(bi_ven_vol_vtk_name, gmeshDir, abaqusInput.nodes(:,2:4), ...
                                    abaqusInput.elems(:,2:5), val_elem, [], []);
                                
                                
    %% we also need to translate the geometry so that the top plane at z=0
    %% no change for the connectivity
%     zmax = max(abaqusInput.nodes(:,4));
%     surface_nodes_translated = abaqusInput.surface_nodes;
%     surface_nodes_translated(:,3) = surface_nodes_translated(:,3) - zmax;
    
    %% translate seems imposing difficulty on the motion in the basal, it
    %% might be better just keep the same, since in one cardiac cycle, the
    %% most moved part is near base, and the apex is motionless.
    surface_nodes_translated = abaqusInput.surface_nodes;
    
    %% somehow need to maintain the relation of the centre point
    %% because every geo is tralated to its own centre point, so here we will 
    %% need to displace the end-diastolic geometry
    %% Note the original point in the diastolic geo is in the second plane 
    %% used for reconstruction, usually will be similar as in the early diastole
    
    if (idx_phase == 2) %%retract some results from early_diastole
        cd(resultDir);
        cd(list_phase{idx_phase-1});
        phase_resultDir_early_diastole = pwd();
        cd(workingDir);
        cd(phase_resultDir_early_diastole);
        load basal_centre;
        centre_coor_early_diastole = basal_centre.centre_coor;
        clear basal_centre;
        cd(workingDir);
    end
    
    
    if (idx_phase == 2)
        cd(phase_resultDir);
        load basal_centre;
        centre_coor_end_diastole = basal_centre.centre_coor;
        clear basal_centre;
        cd(workingDir);
    end
    
    if (idx_phase == 2)
        trans_vec_end_early = centre_coor_end_diastole - centre_coor_early_diastole;
    else
        trans_vec_end_early = [0, 0, 0];
    end
    surface_nodes_translated(:,1) = surface_nodes_translated(:,1) + trans_vec_end_early(1);
    surface_nodes_translated(:,2) = surface_nodes_translated(:,2) + trans_vec_end_early(2);
    surface_nodes_translated(:,3) = surface_nodes_translated(:,3) + trans_vec_end_early(3);
    
    bi_ven_surf_vtk_name = sprintf('bi_ven_surf_%s_trans.vtk', phase_selected);
    write_vtk_trigular_surface(bi_ven_surf_vtk_name, gmeshDir, surface_nodes_translated, ...
                                    abaqusInput.surface_elements, val_elem);
end


if idx == 3
    %% in this step, it will generate Mech.inp for further simulation
    cd(gmeshDir);
    load abaqusInput;
    cd(workingDir);
    Ver = abaqusInput.nodes;
    TMesh = abaqusInput.elems;
    %write out an abaqus input file, for setting up livingheart project
    %%%generate abaqus input file with adjust node number
    input_file_name = sprintf('Mech_%s.inp', phase_selected);
    cd(gmeshDir)
    fid = fopen(input_file_name,'w');
    cd(workingDir);
    
    fprintf(fid, '*************************************   \n');
    fprintf(fid, '*HEADING\n');
    fprintf(fid,'*************************************   \n');
    %fprintf(fid,'*Part, name=PART-1\n');
    fprintf(fid,'*NODE\n');
    %%%node output
    for i = 1 : size(Ver,1)
        fprintf(fid, '%d,\t%f,\t%f,\t%f\n', Ver(i,1), Ver(i,2), Ver(i,3), Ver(i,4));
    end
    %%%element output
    fprintf(fid,'*Element, type=C3D4\n');
    for i = 1 : size(TMesh,1)
        fprintf(fid,'%d,\t%d,\t%d,\t%d,\t%d\n', TMesh(i,1),TMesh(i,2),TMesh(i,3),...
            TMesh(i,4),TMesh(i,5));
    end
    %fprintf(fid,'*End Part\n');
    %fprintf(fid,'**\n');
    %fprintf(fid,'**\n');
    %fprintf(fid,'** ASSEMBLY\n');
    %fprintf(fid,'*Assembly, name=Assembly\n');
    %fprintf(fid,'**\n');
    %fprintf(fid,'*Instance, name=PART-1-1, part=PART-1\n');
    %fprintf(fid,'*End Instance\n');
    %fprintf(fid,'**\n');
    %fprintf(fid,'*Nset, nset=ALL, instance=PART-1-1, generate\n');
    %fprintf(fid,'1, %d, 1\n',size(Ver,1));
    %fprintf(fid,'*Elset, elset=SOLID_BODY, instance=PART-1-1, generate\n');
    %fprintf(fid,'1, %d, 1\n', size(TMesh,1));
    %fprintf(fid,'*End Assembly\n');


    fclose(fid);
end

if idx == 4
    cd(gmeshDir);
    load abaqusInput;
    cd(workingDir);
    Ver = abaqusInput.nodes;
    TMesh = abaqusInput.elems;
    
    %%% here are some discrepancies for the libmesh code
    % endo_RV is the acutal free wall
    % endo_rv_freewall is the septum
    TMeshLibmesh = TMesh;
    epiID = 4096;     epi_name = 'SURF_EPI'; 
    baseID = 4098;    base_name = 'SURF_BASE';
    endoID_LV = 4097; endo_LV_name = 'SURF_LV_ENDO';
    endoID_RV = 5090; endo_RV_name = 'SURF_RV_ENDO';
    endoID_septum_RV = 5091; endo_septum_RV_name = 'SURF_RV_SEPTUM_ENDO';
    endoID_freewall_RV = endoID_RV;endo_freewall_RV_name = 'SURF_RV_FREEWALL_ENDO';
    
    domainID = 1;
    
    %%%now figure out the faces for libMesh
    surfaceSets = abaqusInput.surfaceSets;
    
    %%now need to figure out the SURF_LV_FREEWALL_ENDO set 
    for surfIndex = 1 : size(surfaceSets,1)
        surf_name = surfaceSets(surfIndex,1).surf_name;
        if strcmp(surf_name, endo_RV_name)
            faceSets_RV_ENDO = surfaceSets(surfIndex,1).faceSets;
        elseif strcmp(surf_name, endo_septum_RV_name)
            faceSets_RV_SEPTUM = surfaceSets(surfIndex,1).faceSets;
        end
    end
    
    for faceIndex = 1 : size(faceSets_RV_ENDO,1)
        str_face_RV = faceSets_RV_ENDO(faceIndex,1).str_face;
        elem_list_RV = faceSets_RV_ENDO(faceIndex,1).elem_list;
        
        elem_list_septum_to_delete = [];
        for findex = 1 : size(faceSets_RV_SEPTUM,1)
            str_face_septum = faceSets_RV_SEPTUM(findex,1).str_face;
            elem_list_septum = faceSets_RV_SEPTUM(findex,1).elem_list;
            
            if strcmp(str_face_RV, str_face_septum)
               elem_list_septum_to_delete = elem_list_septum;
            end
        end
        
        if ~isempty(elem_list_septum_to_delete)
            elem_list_RV_freewall = setdiff(elem_list_RV, elem_list_septum_to_delete);
        else
            elem_list_RV_freewall = elem_list_RV;
        end
        
        if size(elem_list_RV_freewall,1) == 1
            elem_list_RV_freewall = elem_list_RV_freewall';
        end
        
        faceSets_RV_freewall(faceIndex,1).surf_str_name = '_SURF_RV_FREEWALL_ENDO';
        faceSets_RV_freewall(faceIndex,1).str_face = str_face_RV;
        faceSets_RV_freewall(faceIndex,1).elem_list = elem_list_RV_freewall;
        
    end
    
    %%add this to abaqusInput
    surfIndex = size(surfaceSets,1)+1;
    surfaceSets(surfIndex,1).surf_name = 'SURF_RV_FREEWALL_ENDO';
    surfaceSets(surfIndex,1).faceSets = faceSets_RV_freewall;
    
    

    
    for surfIndex = 1 : size(surfaceSets,1)
        surfID = 0;
        surf_name = surfaceSets(surfIndex,1).surf_name;
        %%now decide which surface 
        if strcmp(surf_name, epi_name)
            surfID = epiID;
        elseif strcmp(surf_name, base_name)
            surfID = baseID;
        elseif strcmp(surf_name, endo_LV_name)
            surfID = endoID_LV;
        elseif strcmp(surf_name, endo_RV_name) % we do not output the%whole RV endo surface
            surfID = endoID_RV;
        elseif strcmp(surf_name, endo_freewall_RV_name)
            surfID = endoID_freewall_RV;
        elseif strcmp(surf_name, endo_septum_RV_name)
            surfID = endoID_septum_RV;
        else
            disp('that surface is not defined, please check');
        end
        
        if ~strcmp(surf_name, endo_RV_name)
            faceSets = surfaceSets(surfIndex,1).faceSets;
            for faceIndex = 1 : size(faceSets,1)
                elem_list = faceSets(faceIndex,1).elem_list;
                str_face  = faceSets(faceIndex,1).str_face;

                if strcmp(str_face, 'S1')
                    for i = 1 : length(elem_list)
                       TMeshLibmesh(elem_list(i),6)=1;
                       TMeshLibmesh(elem_list(i),10)=surfID;
                    end

                elseif strcmp(str_face, 'S2')
                    for i = 1 : length(elem_list)
                       TMeshLibmesh(elem_list(i),7)=1;
                       TMeshLibmesh(elem_list(i),11)=surfID;
                    end

                elseif strcmp(str_face, 'S3')
                    for i = 1 : length(elem_list)
                       TMeshLibmesh(elem_list(i),8)=1;
                       TMeshLibmesh(elem_list(i),12)=surfID;
                    end

                elseif strcmp(str_face, 'S4')
                    for i = 1 : length(elem_list)
                       TMeshLibmesh(elem_list(i),9)=1;
                       TMeshLibmesh(elem_list(i),13)=surfID;
                    end

                else
                    disp('the face is not defined properly, please check');
                end

            end
        end
        
        
    end
    
    %%finally out the files for LibMesh
    cd(gmeshDir);
    fid=fopen('heart_real.1.node','w');
    cd(workingDir)
    fprintf(fid, '%d  %d   %d   %d\n', size(Ver,1), 3, 0, 0);
    for i = 1 : size(Ver,1)
        fprintf(fid, '%d\t%f\t%f\t%f\n', Ver(i,1), Ver(i,2), Ver(i,3), Ver(i,4));
    end
    fclose(fid);

    cd(gmeshDir);
    fid=fopen('heart_real.1.ele','w');
    fidLibmesh = fopen('heart_real.1.ele.Libmesh.dat','w');
    cd(workingDir);
    
    fprintf(fid, '%d   %d   %d\n', size(TMesh,1), 4, 0);
    
    fprintf(fidLibmesh,'%d\t4\t0\n',size(TMesh,1));
    for i = 1 : size(TMesh,1)
        fprintf(fid,'%d\t%d\t%d\t%d\t%d\n', TMesh(i,1),TMesh(i,2),TMesh(i,3),...
                                            TMesh(i,4),TMesh(i,5));
        fprintf(fidLibmesh,'%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n', ....
                  TMeshLibmesh(i,1),TMeshLibmesh(i,2),TMeshLibmesh(i,3),TMeshLibmesh(i,4),...
                  TMeshLibmesh(i,5),TMeshLibmesh(i,6),TMeshLibmesh(i,7),TMeshLibmesh(i,8),...
                  TMeshLibmesh(i,9),TMeshLibmesh(i,10),TMeshLibmesh(i,11),TMeshLibmesh(i,12),...
                  TMeshLibmesh(i,13), domainID);
    end
    fclose(fid);
    fclose(fidLibmesh);
    
    
end

%%% here trying to seperate LV and RV by looking at LV endo nodes and RV
%%% fee wall nodes
if idx == 5
    cd(gmeshDir);
    load abaqusInput;
    cd(workingDir);
    nodes = abaqusInput.nodes;
    elems = abaqusInput.elems;
    
    
    %%%now figure out the faces for libMesh
    surfaceSets = abaqusInput.surfaceSets;
    
    epi_name = 'SURF_EPI'; 
    base_name = 'SURF_BASE';
    endo_LV_name = 'SURF_LV_ENDO';
    endo_RV_name = 'SURF_RV_ENDO';
    endo_septum_RV_name = 'SURF_RV_SEPTUM_ENDO';
    endo_freewall_RV_name = 'SURF_RV_FREEWALL_ENDO';
    
    %%now need to figure out the SURF_LV_FREEWALL_ENDO set 
    for surfIndex = 1 : size(surfaceSets,1)
        surf_name = surfaceSets(surfIndex,1).surf_name;
        if strcmp(surf_name, endo_RV_name)
            faceSets_RV_ENDO = surfaceSets(surfIndex,1).faceSets;
        elseif strcmp(surf_name, endo_septum_RV_name)
            faceSets_RV_SEPTUM = surfaceSets(surfIndex,1).faceSets;
        end
    end
    
    for faceIndex = 1 : size(faceSets_RV_ENDO,1)
        str_face_RV = faceSets_RV_ENDO(faceIndex,1).str_face;
        elem_list_RV = faceSets_RV_ENDO(faceIndex,1).elem_list;
        
        elem_list_septum_to_delete = [];
        for findex = 1 : size(faceSets_RV_SEPTUM,1)
            str_face_septum = faceSets_RV_SEPTUM(findex,1).str_face;
            elem_list_septum = faceSets_RV_SEPTUM(findex,1).elem_list;
            
            if strcmp(str_face_RV, str_face_septum)
               elem_list_septum_to_delete = elem_list_septum;
            end
        end
        
        if ~isempty(elem_list_septum_to_delete)
            elem_list_RV_freewall = setdiff(elem_list_RV, elem_list_septum_to_delete);
        else
            elem_list_RV_freewall = elem_list_RV;
        end
        
        if size(elem_list_RV_freewall,1) == 1
            elem_list_RV_freewall = elem_list_RV_freewall';
        end
        
        faceSets_RV_freewall(faceIndex,1).surf_str_name = '_SURF_RV_FREEWALL_ENDO';
        faceSets_RV_freewall(faceIndex,1).str_face = str_face_RV;
        faceSets_RV_freewall(faceIndex,1).elem_list = elem_list_RV_freewall;
        
    end   
    
    RV_freewall_nodes_list_unsorted = [];
    for faceID = 1 : size(faceSets_RV_freewall,1)
           str_face = faceSets_RV_freewall(faceID,1).str_face;
           elem_list = faceSets_RV_freewall(faceID,1).elem_list;

           for i = 1 : length(elem_list)
               elem_id = elem_list(i);
               tet4_node_list = elems(elem_id, 1:5);

               if strcmp(str_face, 'S1')
                   node_list = [tet4_node_list(2), tet4_node_list(3), tet4_node_list(4)];
               elseif strcmp(str_face, 'S2')
                   node_list = [tet4_node_list(2), tet4_node_list(5), tet4_node_list(3)];
               elseif strcmp(str_face, 'S3')
                   node_list = [tet4_node_list(3), tet4_node_list(5), tet4_node_list(4)];
               elseif strcmp(str_face, 'S4')
                   node_list = [tet4_node_list(4), tet4_node_list(5), tet4_node_list(2)];
               else
                   disp('wrong at function extract_surface_mesh()');
               end

               RV_freewall_nodes_list_unsorted = [RV_freewall_nodes_list_unsorted; node_list];
           end %%elem_list      
    end
    RV_freewall_nodes_list = unique(RV_freewall_nodes_list_unsorted(:));
    
    RV_septum_nodes_list_unsorted = [];
    for faceID = 1 : size(faceSets_RV_SEPTUM, 1)
        stf_face  = faceSets_RV_SEPTUM(faceID,1).str_face;
        elem_list = faceSets_RV_SEPTUM(faceID,1).elem_list;
        node_list = [];
        for i = 1 : length(elem_list)
           elem_id = elem_list(i);
           tet4_node_list = elems(elem_id, 1:5);
           if strcmp(str_face, 'S1')
                   node_list = [tet4_node_list(2), tet4_node_list(3), tet4_node_list(4)];
           elseif strcmp(str_face, 'S2')
                   node_list = [tet4_node_list(2), tet4_node_list(5), tet4_node_list(3)];
           elseif strcmp(str_face, 'S3')
                   node_list = [tet4_node_list(3), tet4_node_list(5), tet4_node_list(4)];
           elseif strcmp(str_face, 'S4')
                   node_list = [tet4_node_list(4), tet4_node_list(5), tet4_node_list(2)];
           else
                   disp('wrong at function extract_surface_mesh()');
           end
           RV_septum_nodes_list_unsorted = [RV_septum_nodes_list_unsorted; node_list];
        end
    end
    RV_septum_nodes_list_unsorted = unique(RV_septum_nodes_list_unsorted(:));
    
    
    %% extract other node sets
    nodeSets = abaqusInput.nodeSets;
    for nodeSetsID = 1 : size(nodeSets, 1)
        if strcmp( nodeSets(nodeSetsID,1).str_node_set, 'NODE_LV_ENDO' )
            LV_endo_nodes_list = nodeSets(nodeSetsID,1).nodelist;
        end
        
        if strcmp( nodeSets(nodeSetsID,1).str_node_set, 'NODE_RV_ENDO')
            RV_endo_nodes_list = nodeSets(nodeSetsID,1).nodelist;
        end
        
    end
    
    RV_freewall_nodes_list_checked=intersect(RV_freewall_nodes_list,RV_endo_nodes_list);
    if length(RV_freewall_nodes_list_checked) ~= length(RV_freewall_nodes_list)
        disp('RV free wall nodes may not be right, please check');
    end
    
    %%update LV_endo_nodes_list by including RV_septum_nodes_list_unsorted;
    LV_ENDO_RV_SEPTUM = [LV_endo_nodes_list; RV_septum_nodes_list_unsorted];
    LV_dis = patientConfigs(patientIndex,1).LV_dis;
    RV_dis = patientConfigs(patientIndex,1).LV_dis;
    LV_RV_assignment = determine_LV_RV_by_distance_to_inner_surface(nodes, elems,...
    LV_ENDO_RV_SEPTUM, RV_freewall_nodes_list, LV_dis, RV_dis);
    

    %% here is to check whether the element assignement is right
    node_assign = zeros(size(LV_RV_assignment.node_assign)); 
    for el = 1 : size(elems,1)
        el_list = elems(el, 2:5);
        for nl = 1 : length(el_list)
           node_assign(el_list(nl),1 ) = LV_RV_assignment.elem_assign(el,1);
        end
    end


    cd(gmeshDir)
    fid_lv_rv_assign = fopen('LV_RV_assignment.dat','w');
    cd(workingDir);
    TecplotMeshRegions(nodes, elems,node_assign,fid_lv_rv_assign);
    fclose(fid_lv_rv_assign);
    
    
    cd(gmeshDir);
    save LV_RV_assignment LV_RV_assignment;
    cd(workingDir);
    
end %% end for idx == 5



























