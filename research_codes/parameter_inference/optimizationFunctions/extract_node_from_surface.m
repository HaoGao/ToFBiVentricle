%%% extract the node for the RV septum 
function nodelistTotal = extract_node_from_surface(abaqusInput, surfaceName)

elems = abaqusInput.elems;

surfaceSets = abaqusInput.surfaceSets;

for i = 1 : size(surfaceSets, 1)
   if strcmp (surfaceSets(i).surf_name, surfaceName)
      faceSets =  surfaceSets(i).faceSets;
   end
end

nodelistTotal = [];

for faceID = 1 : size(faceSets,1)
    str_face = faceSets(faceID).str_face;
    elem_list = faceSets(faceID).elem_list;
    
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
          nodelistTotal = [nodelistTotal; node_list]; 
    end
    
    
   
    
    
end

nodelistTotal = unique(nodelistTotal(:));