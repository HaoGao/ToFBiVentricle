function strain_segs= segRegionsStrainSummarization(optimize_opt, fiberStrain_cra, endoNodes)
% AHALVMeshDivision, fiberStrain_cra, elem, endoNodes

% optimize_opt.elRegionsFull = elRegionsFull;
% optimize_opt.nodeElRegionFull = nodeElRegionFull;

elRegionsFull = optimize_opt.elRegionsFull; %%means we can use the top 4 slices, usually the first two slices are combined

for i = 1 : max(elRegionsFull(:))   
    strainDataAbaqusSegRegions_segs(i,1).strain = [];
    strain_segs(i,1) = 0; 
end

elemT = optimize_opt.abaqusInput.elems;
elem(:,1) = elemT(:,2);
elem(:,2) = elemT(:,3);
elem(:,3) = elemT(:,4);
elem(:,4) = elemT(:,5);





for elemIndex = 1 : size(elem,1)
    nodeSeq = elem(elemIndex, :);
    %% decide whether this elemment in endo 
    onEndo = 0;
    for nodeIndex = 1 : length(nodeSeq)
        if (~isempty(find(endoNodes == nodeSeq(nodeIndex), 1)) )
            onEndo = 1;
            break;
        end
    end
    
    if  onEndo  %%that will be save for processing
        regionIndex = elRegionsFull(elemIndex);
        strain = strainDataAbaqusSegRegions_segs(regionIndex,1).strain;
        strain = [strain fiberStrain_cra(elemIndex)];
        strainDataAbaqusSegRegions_segs(regionIndex,1).strain = strain;
    end
    
    
end


for i = 1 : max(elRegionsFull)
    strain_segs(i) = mean(strainDataAbaqusSegRegions_segs(i,1).strain);
end
