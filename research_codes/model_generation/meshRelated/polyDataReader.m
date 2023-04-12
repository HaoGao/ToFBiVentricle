function polyMesh = polyDataReader(fileName, fileDir)

% fileName = 'DeterministicAtlas__flow__putamen__subject_subj1__tp_0.vtk';
% fileDir = '/Users/haogao/Dropbox/work/BiVentricularReconstruction/Results/HV1/deformetrica/output/trial3_kernelWidth10';


workingDir = pwd();
cd(fileDir);
fid = fopen(fileName);
cd(workingDir);


msg_to_show = sprintf('reading: %s', fileName);
disp(msg_to_show);

%%first need to locate the position of different sections
tline = fgetl(fid);

nodes = [];
elems = [];

while ischar(tline)
    strT = strsplit(tline);
    if size(strT,2)>=1
       str_1 = strT{1};
       %%now we will need to decide how to process accroding to keywords
       
       %%this is node section
       if strcmp(str_1, 'POINTS') || strcmp(str_1, 'Points')
           %%let try to read into node matrix
          NNodes = str2num(strT{2});
          nodes = zeros([NNodes, 3]);
          msg_to_show = sprintf('I will read %d nodes\n', NNodes);
          disp(msg_to_show);
          tline = fgetl(fid);
          nodeIndex = 0;
          while ischar(tline)
              str_nodes = strsplit(tline);
              str_nodes_1 = str_nodes{1};
              if isempty(str2double(str_nodes_1)) || isnan(str2double(str_nodes_1))
                  frewind(fid);
                  fseek(fid, file_previous_position, 'bof');
                  disp('reading Node finished');
                  break;
              else
                 nodeIndex = nodeIndex + 1;
                 nodes(nodeIndex, 1: 3) = [str2double(str_nodes{1}), ...
                                        str2double(str_nodes{2}), ...
                                        str2double(str_nodes{3})];
                 %% since ICEM may not give natural ordering, so need to record the actual node ID
              end
              file_previous_position = ftell(fid);
              tline = fgetl(fid);
             
          end %% for while node loop
          
          polyMesh.nodes = nodes;
       end %% reading Node finished
       
    
       %%%now read the element mainly C3D4
       if  size(strT,2)>=2 && strcmp(str_1, 'POLYGONS') 
           
           NElems = str2num(strT{2});
           NEntries = str2num(strT{3});
           nPerElem = int8(NEntries / NElems);
           
           msg_to_show = sprintf('reading %d polygons with %d per elem\n', ...
                          NElems, nPerElem);
           disp(msg_to_show);
          
           elems = zeros([NElems, nPerElem]);
           
           tline = fgetl(fid);
           elem_id = 0;
           clear file_previous_position;
           while ischar(tline)
               str_elems = strsplit(tline);
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
                        elemT(i) = int32(str2double(str_elems{i}))+1; %% index in vtk starts from 0
                   end
                   elem_id = elem_id + 1;
                   elemT(1) = elemT(1) - 1; %% this is elem type
                   elems(elem_id,:) = elemT;
               end
               
               file_previous_position = ftell(fid);
               tline = fgetl(fid);
           end %% while elemenet
           polyMesh.elems = elems;
       end % end for reading element
       
    end
    
     tline = fgetl(fid);
    
end


fclose(fid);