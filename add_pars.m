i_mdl = 0;
for i_nps = 1:N_pars

    b_is = nchoosek(all_pars,i_nps);
   for j_b = 1:size(b_is,1)
           i_mdl = i_mdl+1;
 par_ij{i_mdl} = b_is(j_b,:);
end
end