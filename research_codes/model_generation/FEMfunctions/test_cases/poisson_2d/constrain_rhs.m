function f = constrain_rhs(f, ftemp, elem_dofs)

N = length(elem_dofs);

for i = 1 : N
      ig = elem_dofs(i);
      f(ig) = f(ig) + ftemp(i);
end