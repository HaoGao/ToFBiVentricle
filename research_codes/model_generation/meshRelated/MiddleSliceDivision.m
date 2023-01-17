function segRegions = MiddleSliceDivision(nodeMat,zUpper,zLower,centerPoint, MidConfig, LV_RV_assignment)

for nodeIndex = 1 : size(nodeMat,1)
    segRegions(nodeIndex) = 0;
    node_assign_t = LV_RV_assignment.node_assign(nodeIndex);
    pT = nodeMat(nodeIndex,2:4);
    if pT(3)>=zLower && pT(3)<=zUpper && node_assign_t == 1
        theta = degreeCalculationPointBased([pT(1) pT(2)],centerPoint)*180/pi; %%%need to change it to be degree
        regionValue = assignSegAccordingToThetaForMiddleRegion(theta,MidConfig);
        segRegions(nodeIndex)=regionValue;
    elseif node_assign_t == 2
        segRegions(nodeIndex)=18; %% according to AHA17, the RV starts from 18
    end
end


