clear all; clc;
path(path, '../src');
%%test a tetrahedra element 
%% all position vectors are stored as
%% x: the first column
%% y: the second column
%% z: the third column
XCOOR = [0, 0, 0; 
        1, 0, 0; 
        1, 1, 1; 
        0, 0, 1]; %% in reference configuration
du = [0 0 0; 
     1, 0 ,0; 
     0, 1, 0; 
     0, 0, 1]; 

x = XCOOR + du;
NN = 4;
I = [1, 0, 0; 0, 1, 0; 0, 0, 1];


[gpcoor, gpweight] = gp_tetra_lin;

gpNumber = length(gpweight);
for atGP = 1 : gpNumber
    [Xquad, shape, dshape, detvol] = shape_tetra_lin(XCOOR,gpcoor(atGP,:));  
    F = defgrad_3d(du, dshape); %% equal to du' * dshape
    du_x = du(:,1);
    du_y = du(:,2);
    du_z = du(:,3);
    F_x = [0 0 0];
    F_y = [0 0 0];
    F_z = [0 0 0];
    for i = 1 : NN
       dphi = dshape(i,:);
       F_x  = F_x +  du_x(i)*dphi;
       F_y  = F_y +  du_y(i)*dphi;
       F_z  = F_z +  du_z(i)*dphi; 
    end
    FF = [F_x; F_y; F_z] + I;
    
    dv = gpweight(atGP)*detvol;  
    xquad = Xquad + shape*du;%%xquad at current configuration
end


F_average = tetra_defgrad_at_mass_centre(x, XCOOR);


