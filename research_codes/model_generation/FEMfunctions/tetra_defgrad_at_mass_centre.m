function F = tetra_defgrad_at_mass_centre(xcoor, XCOOR)
%% this will return only one F for a tetra element at mass centre
%% xcoor: current position:  Nel * 3
%% XCOOR: reference position:  Nel * 3

du = xcoor - XCOOR;
[gpcoor, gpweight] = gp_tetra_lin;

gpNumber = length(gpweight);
Ftotal = zeros(3);
for atGP = 1 : gpNumber
    [XquadT, shape, dshape, detvol] = shape_tetra_lin(XCOOR,gpcoor(atGP,:));  
    F = defgrad_3d(du, dshape); %% equal to du' * dshape
    Ftotal = Ftotal + F;
end
F = Ftotal ./ gpNumber;
