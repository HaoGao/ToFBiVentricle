function  S_dis = varifold_surface_inner_produce(c_a, c_b, nalpha, nbeta, sigma)

S_dis = 0;
for n = 1 : size(c_a,1)
    c_a_0 = c_a(n,:);
    n_a_0 = nalpha(n,:);
    dis_c_b = c_b - c_a_0;
    dis_c_b_r2 = dis_c_b(:,1).^2 + dis_c_b(:,2).^2 + dis_c_b(:,3).^2; 
    K_0 = exp(-dis_c_b_r2./ (sigma * sigma));
    inner_normal = nbeta*( n_a_0');
    inner_normal = inner_normal.^2;
    S_dis = S_dis + sum(K_0 .*inner_normal);
end