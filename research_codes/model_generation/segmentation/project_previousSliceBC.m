function [endo_lv_p, endo_rv_p, epi_biv_p] = project_previousSliceBC(imIndex, DataSegSA, rect)

   endo_lv = DataSegSA(imIndex).endo_lv;
   endo_rv = DataSegSA(imIndex).endo_rv;
   epi_biv = DataSegSA(imIndex).epi_c;
   
   endo_lv_p(1,:) = endo_lv(1,:) - rect(1);
   endo_lv_p(2,:) = endo_lv(2,:) - rect(2);
   
   endo_rv_p(1,:) = endo_rv(1,:) - rect(1);
   endo_rv_p(2,:) = endo_rv(2,:) - rect(2);
   
   epi_biv_p(1,:) = epi_biv(1,:) - rect(1);
   epi_biv_p(2,:) = epi_biv(2,:) - rect(2);