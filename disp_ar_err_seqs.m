err_mnt = 20;

kappa_err_mnt = 1/(20*degRad)^2;

err_ms(~done) = mod(err_ms(~done)*ar_error_mnt + ...
                             vmrand(0, kappa_err_mnt, [n_not_done 1]) +pi,2*pi) -pi; 