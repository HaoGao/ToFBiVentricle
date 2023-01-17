clear all; clc;
path(path, '../');
%%test a hex element 
XCOOR = [1, -1, -1;
         1, 1, -1;
         -1, 1, -1;
         -1, -1, -1;
         1, -1, 1;
         1, 1, 1;
         -1, 1, 1;
         -1, -1, 1];%% in reference configuration
du = [0 0 0; 
     0, 0 ,0; 
     0, 0, 0; 
     0, 0, 0; 
     0 0 2; 
     0, 0 ,2; 
     0, 0, 2; 
     0, 0, 2]; 

x = XCOOR + du;

[gpcoor, gpweight] = gp_hex_lin;

gpNumber = length(gpweight);
Vtotal = 0.0;
Xquad = zeros(gpNumber, 3);
xquad = Xquad;
Ftotal = zeros(3);
for atGP = 1 : gpNumber
    [XquadT, shape, dshape, detvol] = shape_hex_lin(XCOOR,gpcoor(atGP,:));  
    F = defgrad_3d(du, dshape); %% equal to du' * dshape
    dv = gpweight(atGP)*detvol;  
    Vtotal = Vtotal + dv;
    xquadT = XquadT + shape*du;%%xquad at current configuration
    Xquad(atGP,1:3) = XquadT;
    xquad(atGP,1:3) = xquadT;
    
    Ftotal = Ftotal + F;
end

Fave = Ftotal./gpNumber;
F_average = hex_defgrad_at_mass_centre(x, XCOOR);


