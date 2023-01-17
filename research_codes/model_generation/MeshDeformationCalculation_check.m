%%deformation gradient calculation and others
%% need to run libmesh displacement interoplation first

clear all;
close all;
clc;

% LVWM_config;
LVWM_config;

%%now need to choose phase to segment
cd(resultDir);
cd(deformetricaDir);
deformetricaDir = pwd();
cd(workingDir);

cd(resultDir);
cd(deformetrica_libMeshDir);
deformetrica_libMeshDir = pwd();
cd(workingDir);

%%now need to choose phase to segment
list_phase = {'early_diastole', 'end_diastole', 'end_systole'};
[idx, tf] = listdlg('ListString', list_phase);
cd(resultDir);
if ~exist(list_phase{idx},'dir')
    mkdir(list_phase{idx});
    cd(list_phase{idx});
    phase_resultDir = pwd();    
else
    cd(list_phase{idx});
    phase_resultDir = pwd();
end
cd(workingDir);
phase_selected = list_phase{idx};

cd(phase_resultDir);
cd(solidworksDir);
solidworksDir = pwd();
cd(gmeshDir);
gmeshDir = pwd();
cd(workingDir);

cd(gmeshDir); 
load abaqusInput;
cd(workingDir);

%%%now we will need to calculate F at each element centre
cd(deformetrica_libMeshDir);
dx_deformetrica_t = load('dis_dx.txt');
dy_deformetrica_t = load('dis_dy.txt');
dz_deformetrica_t = load('dis_dz.txt');
cd(workingDir);
dx_deformetrica = zeros([size(dx_deformetrica_t,1) 1]);
dy_deformetrica = zeros([size(dy_deformetrica_t,1) 1]);
dz_deformetrica = zeros([size(dz_deformetrica_t,1) 1]);
for i = 1 : size(dx_deformetrica_t,1)
   dx_deformetrica(dx_deformetrica_t(i,1),1) =  dx_deformetrica_t(i,2);
   dy_deformetrica(dy_deformetrica_t(i,1),1) =  dy_deformetrica_t(i,2);
   dz_deformetrica(dz_deformetrica_t(i,1),1) =  dz_deformetrica_t(i,2);
end
clear dx_deformetrica_t dy_deformetrica_t dz_deformetrica_t;

elems = abaqusInput.elems;
nodes = abaqusInput.nodes;

cd(gmeshDir);
fvec = load('fibreDir.txt');
svec = load('sheetDir.txt');
cvec = load('circumDir.txt');
rvec = load('radDir.txt');
lvec = load('logiDir.txt');
cd(workingDir);

cvec_node = zeros([size(nodes,1) 3]);
rvec_node = cvec_node;
lvec_node = cvec_node;
fvec_node = cvec_node;
svec_node = cvec_node;

for i = 1 : size(cvec,1)
    nel =  elems(i, 2:5);
    for j = 1 : length(nel)
       cvec_node(nel(j),1:3) = cvec(i, 2:4); 
       rvec_node(nel(j),1:3) = rvec(i, 2:4); 
       lvec_node(nel(j),1:3) = lvec(i, 2:4); 
       fvec_node(nel(j),1:3) = fvec(i, 2:4); 
       svec_node(nel(j),1:3) = svec(i, 2:4); 
    end
end



nodes_current(:, 2) =  nodes(:, 2) + dx_deformetrica;
nodes_current(:, 3) =  nodes(:, 3) + dy_deformetrica;
nodes_current(:, 4) =  nodes(:, 4) + dz_deformetrica;

field_data = cvec_node;
TecplotMeshGenWithFileName(nodes, elems,field_data,deformetrica_libMeshDir,...
                                  'BiVenEearlyDiatole');
TecplotMeshGenWithFileName(nodes_current, elems,field_data,deformetrica_libMeshDir,...
                                  'BiVenEndDiatole'); 
field_data = rvec_node;
TecplotMeshGenWithFileName(nodes, elems,field_data,deformetrica_libMeshDir,...
                                  'BiVenEearlyDiatole_radial');
field_data = fvec_node;
TecplotMeshGenWithFileName(nodes, elems,field_data,deformetrica_libMeshDir,...
                                  'BiVenEearlyDiatole_fibre');  
                              
for el_index = 1 : size(elems,1)
    node_el = elems(el_index,2:5);
    XCOOR = zeros([length(node_el), 3]);
    du = XCOOR;
    for i = 1 : length(node_el)
       XCOOR(i,1:3) = nodes(node_el(i),2:4);
       du(i,1:3) = [dx_deformetrica(node_el(i)), ...
                    dy_deformetrica(node_el(i)), ...
                    dz_deformetrica(node_el(i))];
    end
    xcoor = XCOOR+du;
    
 
    [gpcoor, gpweight] = gp_tetra_lin;
    gpNumber = length(gpweight);
    Ftotal = zeros(3);
    for atGP = 1 : gpNumber
        [XquadT, shape, dshape, detvol] = shape_tetra_lin(XCOOR,gpcoor(atGP,:));  
        F = defgrad_3d(du, dshape); %% equal to du' * dshape
        Ftotal = Ftotal + F;
    end
    Fave = Ftotal ./ gpNumber;
    
     xcoor_t = xcoor;
     xcoor_t(:,1) = xcoor_t(:,1) - xcoor(1,1);
     xcoor_t(:,2) = xcoor_t(:,2) - xcoor(1,2);
     xcoor_t(:,3) = xcoor_t(:,3) - xcoor(1,3);
     
     XCOOR_t = XCOOR;
     XCOOR_t(:,1) = XCOOR_t(:,1) - XCOOR(1,1);
     XCOOR_t(:,2) = XCOOR_t(:,2) - XCOOR(1,2);
     XCOOR_t(:,3) = XCOOR_t(:,3) - XCOOR(1,3);
% 
     Fleast = (xcoor_t')/(XCOOR_t');

    DefGrad(el_index,1).F = Fave;
    DefGrad(el_index,1).Fleast = Fleast;
end

IMatrix = [1 0 0; 0 1 0; 0 0 1];
F11=zeros([size(elems,1), 1]);
F12=zeros([size(elems,1), 1]);
F13=zeros([size(elems,1), 1]);
F21=zeros([size(elems,1), 1]);
F22=zeros([size(elems,1), 1]);
F23=zeros([size(elems,1), 1]);
F31=zeros([size(elems,1), 1]);
F32=zeros([size(elems,1), 1]);
F33=zeros([size(elems,1), 1]);
Ecc = F11;
Err = F11;
Ell = F11;
Ecc_in = F11;
for el_index = 1 : size(DefGrad,1)
    F = DefGrad(el_index,1).Fleast;
    F11(el_index,1) = F(1,1);
    F12(el_index,1) = F(1,2);
    F13(el_index,1) = F(1,3);
    F21(el_index,1) = F(2,1);
    F22(el_index,1) = F(2,2);
    F23(el_index,1) = F(2,3);
    F31(el_index,1) = F(3,1);
    F32(el_index,1) = F(3,2);
    F33(el_index,1) = F(3,3);
    
    crl = [cvec(el_index,2), cvec(el_index,3) cvec(el_index,4); ...
           rvec(el_index,2), rvec(el_index,3) rvec(el_index,4); ...
           lvec(el_index,2), lvec(el_index,3) lvec(el_index,4)];
    E = 1/2.0*(F'*F- IMatrix);
    E_crl = crl*E*crl';
    Ecc(el_index,1) = E_crl(1,1);
    Err(el_index,1) = E_crl(2,2);
    Ell(el_index,1) = E_crl(3,3); 
    
    c0 = [cvec(el_index,2), cvec(el_index,3), cvec(el_index,4)];
    lc = (c0*(F'*F)*c0')^0.5;
    lc_in = 1/lc;
    Ecc_in(el_index,1) = lc_in -1;
    
end

for i = 1 : size(nodes,1)
 Ecc_node_t(i,1).Ecc_node = [];
end

Ecc_node = zeros([size(nodes,1) 1]);
Err_node = Ecc_node;
Ell_node = Ecc_node;
Ecc_in_node = Err_node;
for i = 1 : size(Ecc,1)
    nel =  elems(i, 2:5);
    for j = 1 : length(nel)
       Ecc_node_seq =  Ecc_node_t(nel(j),1).Ecc_node;
       Ecc_node_seq = [Ecc_node_seq; Ecc(i, 1)];
       Ecc_node(nel(j),1) = Ecc(i, 1); 
       Err_node(nel(j),1) = Err(i, 1); 
       Ell_node(nel(j),1) = Ell(i, 1); 
       Ecc_in_node(nel(j),1) = Ecc_in(i,1);
       
       Ecc_node_t(nel(j),1).Ecc_node = Ecc_node_seq;
    end
end

for i = 1 : size(nodes,1)
 Ecc_node_seq(i,1) = mean(Ecc_node_t(i,1).Ecc_node);
end

field_data = [Ecc_node, Ecc_node_seq, Ell_node];
TecplotMeshGenWithFileName(nodes, elems,field_data,deformetrica_libMeshDir,...
                                  'BiVenEearlyDiatole_strain');  






field_val(1,1).field_name = 'F_00'; field_val(1,1).field_data=F11;
field_val(2,1).field_name = 'F_01'; field_val(2,1).field_data=F12;
field_val(3,1).field_name = 'F_02'; field_val(3,1).field_data=F13;

field_val(4,1).field_name = 'F_10'; field_val(4,1).field_data=F21;
field_val(5,1).field_name = 'F_11'; field_val(5,1).field_data=F22;
field_val(6,1).field_name = 'F_12'; field_val(6,1).field_data=F23;

field_val(7,1).field_name = 'F_20'; field_val(7,1).field_data=F31;
field_val(8,1).field_name = 'F_21'; field_val(8,1).field_data=F32;
field_val(9,1).field_name = 'F_22'; field_val(9,1).field_data=F33;

% write_vtk_tet_volume('DefGrad.vtk', deformetrica_libMeshDir, nodes(:, 2:4),...
%                     elems(:, 2:5), [], [],field_val);
                
                
%%here we will need to transform F into crl system
disp('start F calculation');
cd(gmeshDir);
fvec = load('fibreDir.txt');
svec = load('sheetDir.txt');
cvec = load('circumDir.txt');
rvec = load('radDir.txt');
lvec = load('logiDir.txt');
cd(workingDir);
fvec = fvec(:, 2:4);
svec = svec(:, 2:4);
cvec = cvec(:, 2:4);
rvec = rvec(:, 2:4);
lvec = lvec(:, 2:4);
IMatrix = [1 0 0; 0 1 0; 0 0 1];
for el_index = 1 : size(fvec, 1)
    f0 = fvec(el_index, :);
    s0 = svec(el_index, :);
    n0 = cross(f0, s0);
    f0 = NormalizationVec(f0);
    s0 = NormalizationVec(s0);
    n0 = NormalizationVec(n0);
    
    c0 = cvec(el_index,:);
    r0 = rvec(el_index,:);
    l0 = lvec(el_index,:);
    
    Fave = DefGrad(el_index,1).F;
    fsn = [f0(1) f0(2) f0(3); ...
           s0(1) s0(2) s0(3); ...
           n0(1) n0(2) n0(3)];
    F_fsn =  fsn*Fave*(fsn');
    
    crl = [c0(1) c0(2) c0(3); ...
           r0(1) r0(2) r0(3); ...
           l0(1) l0(2) l0(3)];
    F_crl = crl*Fave*(crl');
    
    E = 1/2.0*(F'*F - IMatrix);
    E_fsn = fsn*E*(fsn');
    E_crl = crl*E*(crl');
    
    
    DefGrad_fsn(el_index,1).E = E_fsn;
    DefGrad_fsn(el_index,1).F = F_fsn;
    
    DefGrad_crl(el_index,1).F = F_crl;
    DefGrad_crl(el_index,1).E = E_crl;
end
disp('done with F calculation');


F11=zeros([size(elems,1), 1]);
F12=zeros([size(elems,1), 1]);
F13=zeros([size(elems,1), 1]);
F21=zeros([size(elems,1), 1]);
F22=zeros([size(elems,1), 1]);
F23=zeros([size(elems,1), 1]);
F31=zeros([size(elems,1), 1]);
F32=zeros([size(elems,1), 1]);
F33=zeros([size(elems,1), 1]);

for el_index = 1 : size(DefGrad,1)
    F = DefGrad_fsn(el_index,1).F;
    F11(el_index,1) = F(1,1);
    F12(el_index,1) = F(1,2);
    F13(el_index,1) = F(1,3);
    F21(el_index,1) = F(2,1);
    F22(el_index,1) = F(2,2);
    F23(el_index,1) = F(2,3);
    F31(el_index,1) = F(3,1);
    F32(el_index,1) = F(3,2);
    F33(el_index,1) = F(3,3);
    
end
field_fsn_val(1,1).field_name = 'F_00'; field_fsn_val(1,1).field_data=F11;
field_fsn_val(2,1).field_name = 'F_01'; field_fsn_val(2,1).field_data=F12;
field_fsn_val(3,1).field_name = 'F_02'; field_fsn_val(3,1).field_data=F13;

field_fsn_val(4,1).field_name = 'F_10'; field_fsn_val(4,1).field_data=F21;
field_fsn_val(5,1).field_name = 'F_11'; field_fsn_val(5,1).field_data=F22;
field_fsn_val(6,1).field_name = 'F_12'; field_fsn_val(6,1).field_data=F23;

field_fsn_val(7,1).field_name = 'F_20'; field_fsn_val(7,1).field_data=F31;
field_fsn_val(8,1).field_name = 'F_21'; field_fsn_val(8,1).field_data=F32;
field_fsn_val(9,1).field_name = 'F_22'; field_fsn_val(9,1).field_data=F33;

% write_vtk_tet_volume('DefGrad_fsn.vtk', deformetrica_libMeshDir, nodes(:, 2:4),...
%                     elems(:, 2:5), [], [],field_fsn_val);

field_data = [F11, F12, F13];
TecplotMeshGenWithFileName(nodes, elems,field_data,deformetrica_libMeshDir,...
                                  'DefGrad_fsn_F11_F12_F13.dat');
                              
field_data = [F21, F22, F23];
TecplotMeshGenWithFileName(nodes, elems,field_data,deformetrica_libMeshDir,...
                                  'DefGrad_fsn_F21_F22_F23.dat');

field_data = [F31, F32, F33];
TecplotMeshGenWithFileName(nodes, elems,field_data,deformetrica_libMeshDir,...
                                  'DefGrad_fsn_F31_F32_F33.dat');

%%output the strain 
lambda_f=zeros([size(elems,1), 1]);
lambda_s=zeros([size(elems,1), 1]);
lambda_n=zeros([size(elems,1), 1]);

for el_index = 1 : size(DefGrad,1)
    F = DefGrad(el_index,1).F;
    C = F'*F;
    f0 = fvec(el_index, :);
    s0 = svec(el_index, :);
    n0 = cross(f0, s0);
    f0 = NormalizationVec(f0);
    s0 = NormalizationVec(s0);
    n0 = NormalizationVec(n0);
    
    
    lambda_f(el_index,1) = 0.5*(f0*(C*f0')-1);
    lambda_s(el_index,1) = 0.5*(s0*(C*s0')-1);
    lambda_n(el_index,1) = 0.5*(n0*(C*n0')-1);
    
end


field_data = [lambda_f, lambda_s, lambda_n];
TecplotMeshGenWithFileName(nodes, elems,field_data,deformetrica_libMeshDir,...
                                  'DefGrad_fsn_lambda.dat');
                             




F11=zeros([size(elems,1), 1]);
F12=zeros([size(elems,1), 1]);
F13=zeros([size(elems,1), 1]);
F21=zeros([size(elems,1), 1]);
F22=zeros([size(elems,1), 1]);
F23=zeros([size(elems,1), 1]);
F31=zeros([size(elems,1), 1]);
F32=zeros([size(elems,1), 1]);
F33=zeros([size(elems,1), 1]);

for el_index = 1 : size(DefGrad,1)
    F = DefGrad_crl(el_index,1).F;
    F11(el_index,1) = F(1,1);
    F12(el_index,1) = F(1,2);
    F13(el_index,1) = F(1,3);
    F21(el_index,1) = F(2,1);
    F22(el_index,1) = F(2,2);
    F23(el_index,1) = F(2,3);
    F31(el_index,1) = F(3,1);
    F32(el_index,1) = F(3,2);
    F33(el_index,1) = F(3,3);
    
end


field_data = [F11, F12, F13];
TecplotMeshGenWithFileName(nodes, elems,field_data,deformetrica_libMeshDir,...
                                  'DefGrad_crl_F11_F12_F13.dat');
                              
field_data = [F21, F22, F23];
TecplotMeshGenWithFileName(nodes, elems,field_data,deformetrica_libMeshDir,...
                                  'DefGrad_crl_F21_F22_F23.dat');

field_data = [F31, F32, F33];
TecplotMeshGenWithFileName(nodes, elems,field_data,deformetrica_libMeshDir,...
                                  'DefGrad_crl_F31_F32_F33.dat');

%%output the strain 
lambda_c=zeros([size(elems,1), 1]);
lambda_r=zeros([size(elems,1), 1]);
lambda_l=zeros([size(elems,1), 1]);

for el_index = 1 : size(DefGrad,1)
    F = DefGrad(el_index,1).F;
    C = F'*F;
    c0 = cvec(el_index,:);
    r0 = rvec(el_index,:);
    l0 = lvec(el_index,:);
    
    
    lambda_c(el_index,1) = 0.5*(c0*(C*c0') -1);
    lambda_r(el_index,1) = 0.5*(r0*(C*r0') -1);
    lambda_l(el_index,1) = 0.5*(l0*(C*l0') -1);
    
end


field_data = [lambda_c, lambda_r, lambda_l];
TecplotMeshGenWithFileName(nodes, elems,field_data,deformetrica_libMeshDir,...
                                  'DefGrad_crl_lambda.dat');












