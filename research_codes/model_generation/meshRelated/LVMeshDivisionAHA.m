function LVMeshDivisionAHA(resultDir,gmeshDir, phase_resultDir,...
                   MidConfig, ApexConfig, SASliceBase, SASlicePositionApex, ...
                   SASliceDistance,SliceThickness,usuableSXSlice, ...
                   abaqusInputData, LV_RV_assignment)
%%
%% this only applies to tet elements
disp('this only work with basal plane at z=0 and apex at -z');
workingDir = pwd();

nodeMat = abaqusInputData.nodes;
elemMat = abaqusInputData.elems;

BSAn = 0; %whether we have the roated normal direction of the image plane
cd(gmeshDir);
if exist('SANormalRotated.mat', 'file')
    load SANormalRotated;
    BSAn = 1;
end
cd(workingDir);

%% the first slice location will be from zMax to zMax-5;
%% the next slice location will be from (sliceLocation-1)*10-5 to (sliceLocation-1)*10+5
%% the cener point always will be (0 0)
%% also provides another way for define the regions with less than 10mm
%% thickness in order to compare with b-spline recovered strain

for sliceIndex = 1 : SASlicePositionApex-1
    if sliceIndex == 1
        zUpper = SASliceDistance/2;
        zLower = -SASliceDistance/2;
        zUpperSlice = SliceThickness/2;
        zLowerSlice = -SliceThickness/2;
    else
        zUpper = -(sliceIndex-1)*SASliceDistance + SASliceDistance/2;
        zLower = -(sliceIndex-1)*SASliceDistance - SASliceDistance/2;
        zUpperSlice = -(sliceIndex-1)*SASliceDistance + SliceThickness/2;
        zLowerSlice = -(sliceIndex-1)*SASliceDistance - SliceThickness/2;
    end
    centerPoint = [0 0];
    
    if BSAn == 1
        %%%needs to implement with SANormalRotated
        %segRegions = MiddleSliceDivision(nodeMat,zUpper,zLower,centerPoint, MidConfig);
        %segRegionsSlicesPlane = MiddleSliceDivision(nodeMat,zUpperSlice,zLowerSlice,centerPoint, MidConfig);
    else %%old version, if the basal plane is not in the z plane, then some errors will be when using the division scheme
        segRegions = MiddleSliceDivision(nodeMat,zUpper,zLower,centerPoint, MidConfig, LV_RV_assignment);
        segRegionsSlicesPlane = MiddleSliceDivision(nodeMat,zUpperSlice,zLowerSlice,centerPoint, MidConfig, LV_RV_assignment);
    end
    
    z_slices_bc(sliceIndex,:) = [zUpper, zLower];
    sliceRegions(sliceIndex).segRegions = segRegions;
    sliceRegions(sliceIndex).segRegionsSlicesPlane = segRegionsSlicesPlane;
    %%%this is for regions just in the slice with 5 mm thickness
    
end

%% this is for apex region
for sliceIndex = SASlicePositionApex : usuableSXSlice
    zUpper = -(sliceIndex-1)*SASliceDistance + SASliceDistance/2;
    zLower = -(sliceIndex-1)*SASliceDistance - SASliceDistance/2;
    zUpperSlice = -(sliceIndex-1)*SASliceDistance + SliceThickness/2;
    zLowerSlice = -(sliceIndex-1)*SASliceDistance - SliceThickness/2;
    centerPoint = [0 0];
    segRegions = ApexSliceDivision(nodeMat,zUpper,zLower,centerPoint, ApexConfig, LV_RV_assignment);
    segRegionsSlicesPlane = ApexSliceDivision(nodeMat,zUpperSlice,zLowerSlice,centerPoint, ApexConfig, LV_RV_assignment);
    sliceRegions(sliceIndex).segRegions = segRegions;
    sliceRegions(sliceIndex).segRegionsSlicesPlane = segRegionsSlicesPlane;
    
    z_slices_bc(sliceIndex,:) = [zUpper, zLower];
end

%% for apex point
segRegions = zeros([1 size(nodeMat,1)]);

for nodeIndex = 1 : size(nodeMat,1)
    node_assign_t = LV_RV_assignment.node_assign(nodeIndex);
    if nodeMat(nodeIndex,4)<zLower && node_assign_t == 1
        segRegions(nodeIndex)=7;
        segRegionsSlicesPlane(nodeIndex) = 7;
    elseif node_assign_t == 2
        segRegions(nodeIndex)=18;
        segRegionsSlicesPlane(nodeIndex) = 18;
    end
end
sliceRegions(usuableSXSlice+1).segRegions = segRegions;
sliceRegions(usuableSXSlice+1).segRegionsSlicesPlane = segRegionsSlicesPlane;


%% now combine segRegions
segRegions = zeros([1 size(nodeMat,1)]);
segRegionsSlicesPlane = zeros([1 size(nodeMat,1)]);
for sliceIndex = 1 : usuableSXSlice+ 1
    segRegionsT = sliceRegions(sliceIndex).segRegions;
    segRegionsSlicesPlaneT = sliceRegions(sliceIndex).segRegionsSlicesPlane;
    for nodeIndex = 1 : length(segRegionsT)
        if segRegionsT(nodeIndex)>0
            segRegions(nodeIndex) = segRegionsT(nodeIndex);
        end
        if segRegionsSlicesPlaneT(nodeIndex)>0
            segRegionsSlicesPlane(nodeIndex) = segRegionsSlicesPlaneT(nodeIndex);
        end
        
    end    
end


%% segRegions and segRegionsSlicesPlane are combined resutls together

%% output the regions
cd(gmeshDir)
fid_regions = fopen('sliceSegments.dat','w');
cd(workingDir);
TecplotMeshRegions(nodeMat, elemMat,segRegions,fid_regions);
fclose(fid_regions);


%% need to plot with images to make sure they are right
%% rotate the nodeMatrx to the MRI coordinate system
nodeMatMRI = rotationBackToMRICoordinateSystemt(nodeMat,phase_resultDir);
cd(gmeshDir)
fid_regions = fopen('sliceSegments_MRICoor.dat','w');
cd(workingDir);
TecplotMeshRegions(nodeMatMRI, elemMat,segRegionsSlicesPlane,fid_regions);
fclose(fid_regions);

% %%%now need to get the 3D dicom image
cd(resultDir);
load imDesired;
cd(workingDir);
% %%one middle short axis
% %one short aixs
% imFileName = SXSlice(3).Time(1).name;
% imFileName = sprintf('%s/%s',dicomDir,imFileName);
%imDataSA = MRIMapToReal(imFileName);
TimeInstanceSelected = 1;
imIndex = SASliceBase+2;
SXSliceSorted = imDesired.SXSlice;
imData = SXSliceSorted(1,imIndex).SXSlice(TimeInstanceSelected,1).imData;
imInfo1 = SXSliceSorted(1,imIndex).SXSlice(TimeInstanceSelected,1).imInfo;
imInfo = infoExtract(imInfo1);
imData = MRIMapToRealWithImageAndHeadData(imData, imInfo);

fileName = 'imDataTecFormat.dat';
imgTecplotOutput(imData,fileName,gmeshDir);


%% now need to figure out the seg regions based on the planes
slices_with_6regions = 1:SASlicePositionApex-1;
slice_with_4regions = SASlicePositionApex : usuableSXSlice;
totalElRegion = 6*length(slices_with_6regions) + 4*length(slice_with_4regions);
%% will need add another two, one for apex and one for RV

%%% we first assign each node a slice region
for ni = 1 : size(nodeMat,1)
    slice_index_node = 0;
    z_node = nodeMat(ni, 4);
    for sliceIndex = 1 : usuableSXSlice
       z_up = z_slices_bc(sliceIndex,1);
       z_bo = z_slices_bc(sliceIndex,2);
       if z_node > z_bo && z_node <= z_up
          slice_index_node = sliceIndex; 
          break;
       end
    end
    
    if z_node < z_slices_bc(usuableSXSlice,2)
        slice_index_node = usuableSXSlice + 1;
    end
    slice_No_assigned_Node(ni,1) = slice_index_node;
end

%%
centerPoint = [0, 0];
elRegionsFull = [];
nodeElRegionFull = [];

for eli = 1 : size(elemMat,1)
   el_nodes_ID = elemMat(eli, 2:5);
   sliceIndex_nodes = slice_No_assigned_Node(el_nodes_ID,1);
   node_RV_LV = LV_RV_assignment.node_assign(el_nodes_ID); %% RV to be 2
   node_RV_LV = max(node_RV_LV);
   
   
   sliceIndex_min = min(sliceIndex_nodes);
   sliceIndex_ave = round(mean(sliceIndex_nodes ));
   if sliceIndex_min > 0
          sliceIndex = sliceIndex_min;
   else
          sliceIndex = sliceIndex_ave;
   end
   
   if node_RV_LV == 2
       elRegionsFull(eli, 1) = totalElRegion+2; %RV
       nodeElRegionFull(el_nodes_ID,1) = elRegionsFull(eli,1);
   elseif sliceIndex > usuableSXSlice && node_RV_LV==1
       elRegionsFull(eli, 1) = totalElRegion+1; %% LV Apex
       nodeElRegionFull(el_nodes_ID,1) = elRegionsFull(eli,1);
   elseif node_RV_LV==1 && sliceIndex <= usuableSXSlice%% dealing with elements in LV
      
      node_coords_xyz = nodeMat(el_nodes_ID,2:4);
      el_shapeCentre_xyz = [mean(node_coords_xyz(:,1)), ...
                            mean(node_coords_xyz(:,2)), ...
                            mean(node_coords_xyz(:,3))];
       
      theta_centre = degreeCalculationPointBased([el_shapeCentre_xyz(1) el_shapeCentre_xyz(2)],...
                                          centerPoint)*180/pi; 
      theta_el = [];
      for i = 1 : length(el_nodes_ID)
          nx = nodeMat(el_nodes_ID(i), 2);
          ny = nodeMat(el_nodes_ID(i), 3);
          theta_el(i) = degreeCalculationPointBased([nx, ny],...
                                          centerPoint)*180/pi; 
      end
      theta_el_max = max(theta_el);
      theta_el_min = min(theta_el);
      el_cross_zero = 0;
      if theta_el_max - theta_el_min > 180
          el_cross_zero = 1;
      end
      
      k_base_mid = find(slices_with_6regions==sliceIndex);
      k_apex = find(slice_with_4regions==sliceIndex);
      if ~isempty(k_base_mid)
        regionValue = assignSegAccordingToThetaForMiddleRegion(theta_centre,MidConfig);
        if el_cross_zero
            regionValue_min = assignSegAccordingToThetaForMiddleRegion(theta_el_min,MidConfig);
            regionValue_max = assignSegAccordingToThetaForMiddleRegion(theta_el_max,MidConfig);
            regionValue = min([ regionValue_max, regionValue_min]);
        end
        if regionValue>= 0.0 && regionValue < 1.5
        %            elem_basa_InfSept = [elem_basa_InfSept; eli];
                   elRegionsFull(eli,1) = 1 + (sliceIndex-1)*6;
                   nodeElRegionFull(el_nodes_ID,1) = elRegionsFull(eli,1);
        elseif regionValue >=1.5 && regionValue<2.5 
        %            elem_basa_AntSept = [elem_basa_AntSept; eli];
                   elRegionsFull(eli,1) = 2 + (sliceIndex-1)*6;
                   nodeElRegionFull(el_nodes_ID,1) = elRegionsFull(eli,1);
        elseif regionValue >=2.5 && regionValue < 3.5
        %            elem_basa_Ant = [elem_basa_Ant; eli];
                   elRegionsFull(eli,1) = 3 + (sliceIndex-1)*6;
                   nodeElRegionFull(el_nodes_ID,1) = elRegionsFull(eli,1);
        elseif regionValue >=3.5 && regionValue < 4.5
        %            elem_base_AntLat = [elem_base_AntLat; eli];
                   elRegionsFull(eli,1) = 4+(sliceIndex-1)*6;
                   nodeElRegionFull(el_nodes_ID,1) = elRegionsFull(eli,1);
        elseif regionValue >=4.5 && regionValue <5.5
        %            elem_base_InfLat = [elem_base_InfLat; eli];
                   elRegionsFull(eli,1) = 5+(sliceIndex-1)*6;
                   nodeElRegionFull(el_nodes_ID,1) = elRegionsFull(eli,1);
        elseif regionValue >=5.5 && regionValue <=6.5
        %             elem_base_Inf = [elem_base_Inf; eli];
                    elRegionsFull(eli,1) = 6+(sliceIndex-1)*6;
                    nodeElRegionFull(el_nodes_ID,1) = elRegionsFull(eli,1);
         end  
        
      end %% for k_base_mid   
      
      
      if ~isempty(k_apex)
        regionValue = assignSegAccordingToThetaForApexRegion(theta_centre,ApexConfig);
        if el_cross_zero
            regionValue_min = assignSegAccordingToThetaForApexRegion(theta_el_min,ApexConfig);
            regionValue_max = assignSegAccordingToThetaForApexRegion(theta_el_max,ApexConfig);
            regionValue = min([ regionValue_max, regionValue_min]);
        end
        
        if regionValue < 2
    %            elem_apex_Sept = [elem_apex_Sept; eli];
               elRegionsFull(eli,1) = 1+(length(slices_with_6regions))*6+ (sliceIndex- length(slices_with_6regions)-1)*4;
               nodeElRegionFull(el_nodes_ID,1) = elRegionsFull(eli,1);
        elseif regionValue >=2 && regionValue <3.5
    %            elem_apex_Ant = [elem_apex_Ant; eli];
               elRegionsFull(eli,1) = 2+(length(slices_with_6regions))*6+ (sliceIndex- length(slices_with_6regions)-1)*4;
               nodeElRegionFull(el_nodes_ID,1) = elRegionsFull(eli,1);
        elseif regionValue >=3.5 && regionValue < 5
    %            elem_apex_Lat = [elem_apex_Lat; eli];
               elRegionsFull(eli,1) = 3+(length(slices_with_6regions))*6+ (sliceIndex- length(slices_with_6regions)-1)*4;
               nodeElRegionFull(el_nodes_ID,1) = elRegionsFull(eli,1);
        elseif regionValue >= 5 && regionValue <= 6.5
    %            elem_apex_Inf = [elem_apex_Inf; eli];
               elRegionsFull(eli,1) = 4+(length(slices_with_6regions))*6+ (sliceIndex- length(slices_with_6regions)-1)*4;
               nodeElRegionFull(el_nodes_ID,1) = elRegionsFull(eli,1);
        end
        
        
      end %% k_apex
                                      
                                      
   end
    
end



%% output the regions
cd(gmeshDir)
fid_regions = fopen('BiVen_AHASegments_full.dat','w');
cd(workingDir);
TecplotMeshRegions(nodeMat, elemMat,nodeElRegionFull,fid_regions);
fclose(fid_regions);


cd(gmeshDir)
save LVMeshSegDivisions sliceRegions segRegions segRegionsSlicesPlane elRegionsFull nodeElRegionFull slice_No_assigned_Node;
cd(workingDir);







