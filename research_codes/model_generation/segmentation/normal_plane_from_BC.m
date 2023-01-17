function dis_vec =  normal_plane_from_BC(SX_epi_cReal)

PN  = size(SX_epi_cReal,2);
quatorN = int8(PN/4);
p1 = SX_epi_cReal(:,1);
p2 = SX_epi_cReal(:,quatorN);
p3 = SX_epi_cReal(:,2*quatorN);
p4 = SX_epi_cReal(:, 3*quatorN);

v1 = NormalizationVec(p3 - p1);
v2 = NormalizationVec(p4 - p2);

dis_vec = NormalizationVec( cross(v1, v2) );
