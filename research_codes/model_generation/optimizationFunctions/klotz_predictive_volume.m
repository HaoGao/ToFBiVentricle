function Vpred = klotz_predictive_volume(Ped, Ved, V0, Ppred)
An = 27.78;
Bn = 2.76;

V30 = V0 + (Ved-V0)/( (Ped/An)^(1/Bn) );

Vpred = V0 + (V30 - V0)*(Ppred/An)^(1/Bn);
