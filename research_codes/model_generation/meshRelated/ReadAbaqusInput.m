function abaqusInput = ReadAbaqusInput(fileName, fileDir, surfaceList)
%% this will mainly read Node, Element, and different sets if there are

% close all; clear all; clc;
% fileName ='AbaqusInputSets.inp';
% if ismac   
%     fileDir = '/Users/haogao/Dropbox/work/BiVentricularReconstruction/Results/HV1/early_diastole/solidworks/gmesh1';
% else
%     fileDir ='..\..\Results\HV1\early_diastole\solidworks\gmesh1'; 
% end

workingDir = pwd();
cd(fileDir);
fid = fopen(fileName, 'r');
cd(workingDir);

%%first need to locate the position of different sections
tline = fgetl(fid);

nodes = [];
elems = [];
nodeSetID = 0;
elemSetID = 0;
surfaceSetID = 0;
while ischar(tline)
    strT = strsplit(tline, ',');
    if size(strT,2)>=1
       str_1 = strT{1};
       %%now we will need to decide how to process accroding to keywords
       
       %%this is node section
       if strcmp(str_1, '*NODE') || strcmp(str_1, '*Node')
           %%let try to read into node matrix
          disp('reading Node now')
          tline = fgetl(fid);
          nodeIndex = 0;
          while ischar(tline)
              str_nodes = strsplit(tline, ',');
              str_nodes_1 = str_nodes{1};
              if isempty(str2double(str_nodes_1)) || isnan(str2double(str_nodes_1))
                  frewind(fid);
                  fseek(fid, file_previous_position, 'bof');
                  disp('reading Node finished');
                  break;
              else
                 nodeIndex = nodeIndex + 1;
                 nodes(nodeIndex, 1: 4) = [str2num(str_nodes{1}),...
                                        str2double(str_nodes{2}), ...
                                        str2double(str_nodes{3}), ...
                                        str2double(str_nodes{4})];
                 %% since ICEM may not give natural ordering, so need to record the actual node ID
              end
              file_previous_position = ftell(fid);
              tline = fgetl(fid);
             
          end %% for while node loop
          
          abaqusInput.nodes = nodes;
       end %% reading Node finished
       
    
       %%%now read the element mainly C3D4
       if  size(strT,2)>=2 && (strcmp(str_1, '*ELEMENT') || strcmp(str_1, '*Element'))...
               && (strcmp( strtrim(strT{2}), 'type=C3D4')|| strcmp(strtrim(strT{2}), 'TYPE=C3D4')   )
           disp('reading element now');
           elem_type = 'C3D4';
           tline = fgetl(fid);
           elem_id = 0;
           clear file_previous_position;
           while ischar(tline)
               str_elems = strsplit(tline, ',');
               str_elems_1 = str_elems{1};
               if isempty(str2double(str_elems_1)) || isnan(str2double(str_elems_1))
                   disp('reading element finished');
                   frewind(fid);
                   fseek(fid, file_previous_position, 'bof');
                   break;
               else
                   %%now read into the element matrix
                   elemT = [];
                   for i = 1 : size(str_elems,2)
                        elemT(i) = int32(str2double(str_elems{i}));
                   end
                   elem_id = elem_id + 1;
                   elems(elem_id,:) = elemT;
               end
               
               file_previous_position = ftell(fid);
               tline = fgetl(fid);
           end %% while elemenet
           abaqusInput.elems = elems;
           abaqusInput.elem_type = elem_type;
       end % end for reading element
       
       
       %% now read the node sets 
       if strcmp(str_1, '*Nset')
           str_node_set = strsplit(strT{2},'=');
           msg_str = sprintf('reading Nset: %s', str_node_set{2});
           disp(msg_str);
           tline=fgetl(fid);
           nodeSetTemp = [];
           while ischar(tline)
                str_nset_lines = strsplit(tline, ',');
                if isempty(str2double(str_nset_lines{1})) || isnan(str2double(str_nset_lines{1}))
                   disp('reading Nset finished');
                   frewind(fid);
                   fseek(fid, file_previous_position, 'bof');
                   break;
                else
                   %%now read into the Nset matrix
                   nodeTemp=sscanf(tline,'%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f');
                   if ~isempty(nodeTemp)
                       nodeSetTemp=[nodeSetTemp; nodeTemp];
                   end
                end
               file_previous_position = ftell(fid);
                tline = fgetl(fid);
           end %% while Nset
           
            nodeSetID = nodeSetID + 1;
            nodeSets(nodeSetID,1).nodelist = nodeSetTemp;
            nodeSets(nodeSetID,1).str_node_set = strtrim(str_node_set{2});
       end %%read the node sets
       
       %%now need to read surface sets
       if strcmp(str_1, '*Elset')
            str_elem_set = strsplit(strT{2},'=');
            msg_str = sprintf('reading Elem set: %s', str_elem_set{2});
            disp(msg_str);
            tline=fgetl(fid);
            elemSetTemp = [];
            while ischar(tline)
                str_elemset_lines = strsplit(tline, ',');
                if isempty(str2double(str_elemset_lines{1})) || isnan(str2double(str_elemset_lines{1}))
                    disp('raeding Elset finished');
                    frewind(fid);
                    fseek(fid, file_previous_position, 'bof');
                    break
                else
                    elemTemp = sscanf(tline,'%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f');
                    if ~isempty(elemTemp)
                        elemSetTemp = [elemSetTemp; elemTemp];
                    end
                end %% if for readiing a line
                file_previous_position = ftell(fid);
                tline = fgetl(fid);
            end%%end reading a list
            
            elemSetID = elemSetID + 1;
            elemSets(elemSetID,1).elemlist = elemSetTemp;
            elemSets(elemSetID,1).str_elem_set = strtrim(str_elem_set{2});
       end%% read Elset 
       
       %%now need to read surface
       if strcmp(str_1, '*Surface')
            str_surf_set = strsplit(strT{3},'=');
            msg_str = sprintf('reading surface set: %s', str_surf_set{2});
            disp(msg_str);
            tline=fgetl(fid);
            
            faceN = 0;
            while ischar(tline) && size(strsplit(tline, ','),2) >= 2
                surf_str_temp = strsplit(tline, ',');
                surf_str_name = strtrim(surf_str_temp{1});
                surf_str_face = strtrim(surf_str_temp{2});
                surf_str_check_name = sprintf('_%s_%s', strtrim(str_surf_set{2}), surf_str_face);
                %if strcmp(surf_str_check_name, surf_str_name(1:length(surf_str_check_name)))
                if contains(surf_str_name, surf_str_check_name)
                    faceN = faceN + 1;
                    faceSet(faceN,1).surf_str_name = surf_str_name;
                    faceSet(faceN,1).str_face = surf_str_face;
                else
                    disp('reading Surface finished');
                    frewind(fid);
                    fseek(fid, file_previous_position,'bof');
                    break;
                end
                file_previous_position = ftell(fid);
                tline = fgetl(fid);
            end %% while read a line
           surfaceSetID = surfaceSetID + 1;
           surfaceSets(surfaceSetID,1).surf_name = strtrim(str_surf_set{2});
           surfaceSets(surfaceSetID,1).faceSets = faceSet;
           clear faceSet;
       end%% read surface 
       
       
    end %% non blank line
    tline = fgetl(fid);
    
end %%end for while
if exist('nodeSets','var')
    abaqusInput.nodeSets = nodeSets;
end
if exist('elemSets', 'var')
    abaqusInput.elemSets = elemSets;
end
if exist('surfaceSets', 'var')
    abaqusInput.surfaceSets = surfaceSets;
end
fclose(fid);


if exist('surfaceSets', 'var')
    %%now will need to link surface name with surface element list
    surfaceSets = abaqusInput.surfaceSets;
    elemSets = abaqusInput.elemSets;

    for surfaceID = 1 : size(surfaceSets,1)
        surf_name = surfaceSets(surfaceID,1).surf_name;
        faceSets = surfaceSets(surfaceID,1).faceSets;
        for faceID = 1 : size(faceSets,1)
           face_to_surface_name = faceSets(faceID,1).surf_str_name;
           face_str_ID = faceSets(faceID,1).str_face;
           %%%now need to find out te index in the elemlist
           for elemListID = 1 : size(elemSets,1)
               str_elem_set = elemSets(elemListID,1).str_elem_set;
               if strcmp(str_elem_set, face_to_surface_name)
                   faceSets(faceID,1).elem_list = elemSets(elemListID,1).elemlist;
                   break;
               end
           end    
        end
        surfaceSets(surfaceID,1).faceSets = faceSets; 
    end

    abaqusInput.surfaceSets = surfaceSets;
    abaqusInput = extract_surface_mesh(abaqusInput, surfaceList);
end

%%%need to set up a surface mesh
function abaqusInput = extract_surface_mesh(abaqusInput, surfaceList)
nodes = abaqusInput.nodes; 
elems = abaqusInput.elems;
surfaceSets = abaqusInput.surfaceSets;


surface_nodes = [];
surface_elems_unsorted = [];
surface_LV_endo_idx = [];

for surfaceID = 1 : size(surfaceSets,1)
    faceSets = surfaceSets(surfaceID,1).faceSets;
    surfaceName = surfaceSets(surfaceID,1).surf_name;
    
     if any(contains(surfaceList, surfaceName))
        for faceID = 1 : size(faceSets,1)
           str_face = faceSets(faceID,1).str_face;
           elem_list = faceSets(faceID,1).elem_list;

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

               surface_elems_unsorted = [surface_elems_unsorted; node_list];
               surface_LV_endo_idx = [surface_LV_endo_idx; 0];
               if strcmp(surfaceName, 'SURF_LV_ENDO')
                   surface_LV_endo_idx(end) = 1;
               end
               
           end %%elem_list      
        end
     end

end
surface_elems_unsorted_Transpose = surface_elems_unsorted';
surface_nodes_unsorted = unique(surface_elems_unsorted_Transpose(:), 'stable');
for i = 1 : length(surface_nodes_unsorted)
    node_id_global = surface_nodes_unsorted(i);
    surface_nodes_map(node_id_global,1) = i;
    surface_nodes_map(node_id_global,2) = node_id_global;
end


surface_nodes = nodes(surface_nodes_unsorted,2:4);
surface_elems = zeros(size(surface_elems_unsorted));
for elem_id = 1 : size(surface_elems, 1)
    nlist = surface_elems_unsorted(elem_id,1:3);
    surface_elements(elem_id,1:3) = [surface_nodes_map(nlist(1),1), ...
        surface_nodes_map(nlist(2),1), surface_nodes_map(nlist(3),1)];
end

abaqusInput.surface_nodes = surface_nodes;
abaqusInput.surface_elements = surface_elements;
abaqusInput.surface_nodes_map = surface_nodes_map;
abaqusInput.surface_elems_unsorted = surface_elems_unsorted;
abaqusInput.surface_LV_endo_idx = surface_LV_endo_idx;









