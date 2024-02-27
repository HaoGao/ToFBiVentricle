function S_dis = varifold_distance_selected(source, target, sigma,surface_LV_endo_idx)


[c_a, n_a] = get_centers_and_normals(source);
[c_b, n_b] = get_centers_and_normals(target);

%% need to prepare area_a, area_b, n_alpha, n_beta
%% area_a is the L-2 norm
area_a = ((n_a(:,1).^2 + n_a(:,2).^2 + n_a(:,3).^2).^0.5);
area_b = ((n_b(:,1).^2 + n_b(:,2).^2 + n_b(:,3).^2).^0.5);
nalpha = n_a./area_a;
nbeta = n_b./area_b;


%%only keep surface_LV_endo_idx==1
c_a = c_a(surface_LV_endo_idx==1,:);
c_b = c_b(surface_LV_endo_idx==1,:);
nalpha = nalpha(surface_LV_endo_idx==1,:);
nbeta = nbeta(surface_LV_endo_idx==1,:);

debug_check = 0;
if debug_check 
    area_aT = [];
    nalphaT = [];
    for i = 1 : size(n_a, 1)
       area_aT(i,1) = (n_a(i,1)^2 + n_a(i,2)^2 + n_a(i,3)^2)^0.5;  
       nalphaT(i,:) = n_a(i,:)/area_aT(i,1);
    end
    
    for j = 1 : size(n_b,1)
        area_bT(j,1) = (n_b(j,1)^2 + n_b(j,2)^2 + n_b(j,3)^2)^0.5;
        nbetaT(j,:) = n_b(j,:)/area_bT(j,1);
    end
    
    max_area_aDiff = max(abs (area_aT - area_a));
    max_nalpha_diff = max( abs(nalphaT - nalpha));
    max_area_bDiff = max(abs( area_bT - area_b));
    max_nbeta_diff = max(abs(nbetaT - nbeta) );
    
    if max_area_aDiff > 1.0e-6 || max(max_nalpha_diff)>1.0e-6 ...
            || max_area_bDiff>1.0e-6  || max(max_nbeta_diff) > 1.0e-6
       disp('area_a, area_b, nalpha, nbeta are not computed properly in function varifold_distance'); 
    end
end


S_dis_a_b = varifold_surface_inner_produce(c_a, c_b, nalpha, nbeta, sigma);
S_dis_a_a = varifold_surface_inner_produce(c_a, c_a, nalpha, nalpha, sigma);
S_dis_b_b = varifold_surface_inner_produce(c_b, c_b, nbeta, nbeta, sigma);

S_dis.a_a = S_dis_a_a;
S_dis.b_b = S_dis_b_b;
S_dis.a_b = S_dis_a_b;
S_dis.dis = S_dis_a_a + S_dis_b_b - 2*S_dis_a_b;