function endo_interp = samplingBCWithoutIm_interp(endo_bc,bclosed)
%%% the sampling will based on the centre point from endo_c

%%using bspline to interolate
endo_c = endo_bc;
ni=[];nni=[];

if bclosed
    ni=1:length(endo_c(1,:))+1;
    nni=1:0.1:length(endo_c(1,:))+1;
    endo_cc(1,:)=spline(ni,[endo_c(1,:) endo_c(1,1)], nni);
    endo_cc(2,:)=spline(ni,[endo_c(2,:) endo_c(2,1)], nni);
else
    ni=1:length(endo_c(1,:));
    nni=1:0.1:length(endo_c(1,:));
    endo_cc(1,:)=spline(ni,endo_c(1,:), nni);
    endo_cc(2,:)=spline(ni,endo_c(2,:), nni);
end

endo_interp = endo_cc; %redefine

