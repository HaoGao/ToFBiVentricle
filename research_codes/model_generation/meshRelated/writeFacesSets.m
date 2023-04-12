function writeFacesSets(fid, surfaceSets, surfaceName)

for setIndex = 1 : size(surfaceSets, 1)
    surf_str_name = surfaceSets(setIndex).surf_str_name;
    faces = surfaceSets(setIndex).elem_list;
    writeS1234(fid, surf_str_name, faces);
end
fprintf(fid, '*Surface, type=ELEMENT, name=%s\n', surfaceName);
for setIndex = 1 : size(surfaceSets,1)
    surf_str_name = surfaceSets(setIndex).surf_str_name;
    str_face = surfaceSets(setIndex).str_face;
    fprintf(fid, '%s, %s\n', surf_str_name, str_face);
end




function writeS1234(fid, surf_str_name, faces)
    fprintf(fid, '*Elset, elset=%s, internal\n', surf_str_name);
    facesLength = length(faces);
          rows = floor(facesLength/16);
        %   facesIndex = 0;
          for i = 1 : rows
             for j = 1 : 15
                 facesIndex = (i-1)*16 + j; 
                 fprintf(fid, '%d,\t', faces(facesIndex));
             end
                facesIndex = (i-1)*16 + 16;
                fprintf(fid, '%d\n', faces(facesIndex));
          end
          %%%the last row
          if facesIndex < facesLength
              if facesLength-facesIndex>1
                  for i = facesIndex+1: facesLength-1
                      fprintf(fid, '%d,\t', faces(i));
                  end
                     i = i+1;
                     fprintf(fid, '%d\n', faces(i));
              else
                  fprintf(fid, '%d\n', faces(facesIndex+1));
              end
          end
