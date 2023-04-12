function LVMeshDivision_17Segments_AHA(resultDir,gmeshDir, ...
                      MidConfig, ApexConfig, ...
                      SASlicePositionBase, SASlicePositionMiddle, SASlicePositionApex, ...
                      SASliceDistance,usuableSXSlice, abaqusInputData)
              
%% here assume the basal plane is the z=0 plane, and apex is in the z- axis
%%%this only applies to tet elements
workingDir = pwd();

nodeMat = abaqusInputData.nodes;
elemMat = abaqusInputData.elems;

vecLA = [0 0 -1];

%%%now we can define a z vector which will be used to get the perpendicular
%%%distance from that point to the top plane
%%here we always assume top plane is centered at (0,0,0)
zDisToBase = nodeMat(:,4);

%%%output the distance whether it is right 
%%%output the regions
% cd(gmeshDir)
% fid_distance = fopen('zDistanceCoordingToBasalPlane.dat','w');
% cd(workingDir);
% TecplotMeshRegions(nodeMat, elemMat,zDisToBase,fid_distance);
% fclose(fid_distance);


%%%now we can define regions according to zDisToBase
InPlanceCenterPosition = [0,0];

%%%define the basal region 
zBaseLower = (SASlicePositionBase(1)-1)*SASliceDistance;
zBaseUpper = (SASlicePositionBase(2)-1)*SASliceDistance;
segRegions = MiddleSliceDivision_AHA17Segments(nodeMat,zDisToBase, zBaseUpper,zBaseLower,InPlanceCenterPosition, MidConfig);
AHA17.baseSegRegions = segRegions;


%%%define the middle region 
zMidLower = (SASlicePositionMiddle(1)-1)*SASliceDistance;
zMidUpper = (SASlicePositionMiddle(2)-1)*SASliceDistance;
segRegions = MiddleSliceDivision_AHA17Segments(nodeMat,zDisToBase, zMidUpper,zMidLower,InPlanceCenterPosition, MidConfig);
AHA17.midSegRegions = segRegions;

%%define the apical region 
zApexLower = (SASlicePositionApex(1)-1)*SASliceDistance;
zApexUpper = (SASlicePositionApex(2)-1)*SASliceDistance;
segRegions = ApexSliceDivision_AHA17Segments(nodeMat,zDisToBase, zApexUpper,zApexLower,InPlanceCenterPosition, ApexConfig);
AHA17.apexSegRegions = segRegions;

%%the apecial point 5mm away from apical region
zApecialPointLower = (SASlicePositionApex(2)-0.5)*SASliceDistance;
segRegions = zeros([1 size(nodeMat,1)]);
for nodeIndex = 1 : size(nodeMat,1)
    if zDisToBase(nodeIndex)>zApecialPointLower
        segRegions(nodeIndex)=7;
    end
end
AHA17.apexPointSegRegions = segRegions;

clear segRegions;


%%%now combine segRegions
segRegions = zeros([1 size(nodeMat,1)]);
%%basal
segRegionsT = AHA17.baseSegRegions;
for nodeIndex = 1 : length(segRegionsT)
        if segRegionsT(nodeIndex)>0
            segRegions(nodeIndex) = segRegionsT(nodeIndex);
        end
end
%%mid
segRegionsT = AHA17.midSegRegions;
for nodeIndex = 1 : length(segRegionsT)
        if segRegionsT(nodeIndex)>0
            segRegions(nodeIndex) = segRegionsT(nodeIndex)+6;
        end
end
%%apical
segRegionsT = AHA17.apexSegRegions;
for nodeIndex = 1 : length(segRegionsT)
        if segRegionsT(nodeIndex)>0
            segRegions(nodeIndex) = segRegionsT(nodeIndex)+12;
        end
end
%%apical point
segRegionsT = AHA17.apexPointSegRegions;
for nodeIndex = 1 : length(segRegionsT)
        if segRegionsT(nodeIndex)>0
            segRegions(nodeIndex) = segRegionsT(nodeIndex)+10;
        end
end



cd(resultDir)
save LVMeshAHA17_SegDivisions  segRegions AHA17;
cd(workingDir);

cd(resultDir)
fid_AHA17 = fopen('AHA17Segment.dat','w');
cd(workingDir);
TecplotMeshRegions(nodeMat, elemMat,segRegions,fid_AHA17);
fclose(fid_AHA17);





% nodeMat(:,2) = nodeMat(:,2)*abaqusInputData.scaleTomm;
% nodeMat(:,3) = nodeMat(:,3)*abaqusInputData.scaleTomm;
% nodeMat(:,4) = nodeMat(:,4)*abaqusInputData.scaleTomm;

% %%%get 1 face from baseface list
% baseface1 = abaqusInputData.baseFaceS1(1);
% baseface1_nodeList = elemMat(baseface1',2:4);
% baseface_nodeCoor = nodeMat(baseface1_nodeList',2:4);
% 
% vec1 = baseface_nodeCoor(2,:) - baseface_nodeCoor(1,:);
% vec2 = baseface_nodeCoor(3,:) - baseface_nodeCoor(1,:);
% 
% vec1 = NormalizationVec(vec1);
% vec2 = NormalizationVec(vec2);
% 
% vecLA = cross(vec1, vec2);
% vecLA = NormalizationVec(vecLA);
% 
% vecZ = [0 0 1];
% if dot(vecZ,vecLA) > 0
%     vecLA = -vecLA;
% end
%%%now vecLA is pointing apex from basal plane














