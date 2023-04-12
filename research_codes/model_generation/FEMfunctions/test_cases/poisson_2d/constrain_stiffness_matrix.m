function K = constrain_stiffness_matrix(K, Ktemp, elem_dofs)

N = length(elem_dofs);

for i = 1 : N
    for j = 1 : N
      ig = elem_dofs(i);
      jg = elem_dofs(j);
      
      K(ig,jg) = K(ig,jg) + Ktemp(i,j);
        
    end
end


