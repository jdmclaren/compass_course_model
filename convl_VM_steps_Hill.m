% determine equivalent angular deviation of sum of von Mises distrs
% this is technically a convolution of pdf's
% for Normal this is additive, here we need to use a recursive approx
% (Hill 1981 ACM Trans Math Software)
% for convoluting initial and 2nd within-flight step, 
% initialize previous-step convoluted kappa as kappa of initial detect errs


%     if strcmp(plot_opt,'VM')   
%         % initial convolution
%         convl_kap_eq_N_ntr = kap_dtc;

 prod_bes = bes_1_dtc.^n_mnts_stp;
 
convl_kap_eq_N_ntr =  inv_Bess_ratio_Hill(prod_bes);

        % try multiplying by n_mnts_stps
        convl_kap_eq_N_ntr = convl_kap_eq_N_ntr*n_mnts_stp;
            % finally, invert sqrt kappa for equivalent multi-step angular
            % deviation
            equiv_N_errs_ntr_Hill = sqrt(1./convl_kap_eq_N_ntr);

            for kk = 1:n_errs

              equiv_R_Hill(kk) = besseli(1,convl_kap_eq_N_ntr(kk),1)/ ...
                  besseli(0,convl_kap_eq_N_ntr(kk),1);

               equiv_s_Hill(kk) = sqrt(-2*log(equiv_R_Hill(kk))); %  sqrt(2*(1-equiv_R_Hill(kk)));

            end


        plot(equiv_s_Hill*180/pi,s_ntr*180/pi,'--','LineWidth',Line_W)

%     end
