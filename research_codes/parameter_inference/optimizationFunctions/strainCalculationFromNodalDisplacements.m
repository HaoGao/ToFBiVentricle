function fiberStrain = strainCalculationFromNodalDisplacements(optimize_opt, using_fsn)
%abaqus_dis_out_filename,abaqusDir,node,elem, fiberDir, sheetDir

workingDir = pwd();

abaqusDir = optimize_opt.abaqusSimulationDir;
cd(abaqusDir);
dis = load(optimize_opt.abaqus_dis_out_filename); 
cd(workingDir);

disT(:,1:3) = dis(:, 2:4); 

nodeT = optimize_opt.abaqusInput.nodes;
node(:, 1) = nodeT(:,2);
node(:, 2) = nodeT(:,3);
node(:, 3) = nodeT(:,4);
clear nodeT;

elemT = optimize_opt.abaqusInput.elems;
elem(:,1) = elemT(:,2);
elem(:,2) = elemT(:,3);
elem(:,3) = elemT(:,4);
elem(:,4) = elemT(:,5);

if using_fsn 
    fibreDir = optimize_opt.fibreDir;
    sheetDir = optimize_opt.sheetDir;
else
    fibreDir = optimize_opt.cirDir;
    sheetDir = optimize_opt.radDir;
end

nodeDef = node + disT;
fiberStrain = zeros([size(elem,1) 1]);

I = [1 0 0; 0 1 0 ; 0 0 1];
for elemIndex = 1 : size(elem,1)
    nodeSequence = elem(elemIndex,:);
    node_0 =  node(nodeSequence,:);
    node_1 =  nodeDef(nodeSequence,:);
    
    
    %%%now define the F is based on node_1,
    node_0(:,1) = node_0(:,1) - node_0(1,1);
    node_0(:,2) = node_0(:,2) - node_0(1,2);
    node_0(:,3) = node_0(:,3) - node_0(1,3);
    node_0 = node_0';
    
    node_1(:,1) = node_1(:,1) - node_1(1,1);
    node_1(:,2) = node_1(:,2) - node_1(1,2);
    node_1(:,3) = node_1(:,3) - node_1(1,3);
    node_1 = node_1';
    
%     F = node_0/node_1;
%     E = 1/2*(F'*F-I);
    
    
    f = fibreDir(elemIndex,2:4);
    s = sheetDir(elemIndex,2:4);
    
    f = NormalizationVec(f);
    s = NormalizationVec(s);
    n = cross(f,s);
    
    %% now using the way of calculating the strain along myofibre
    %% direction or any direction 
    %% the calculated F is based on early diastole, but the measured strain will be based on end-diastole
    F = node_1/node_0;%% node_0 is in early-diastole
    f_ED = F*f';
    f_ED_L = ( f_ED(1)^2 + f_ED(2)^2 + f_ED(3)^2 )^0.5;
    fiberStrain(elemIndex) = ( 1/(f_ED_L^2)-1)/2;
    
    %% that is projection
%     R = [f(1) f(2) f(3); 
%          s(1) s(2) s(3);
%          n(1) n(2) n(3)];
%      
%      Efsn = R*E*R';
     
     
%      fiberStrain(elemIndex) = Efsn(1,1);
    
    
end
