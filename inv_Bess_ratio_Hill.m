function inv_ratio = inv_Bess_ratio_Hill(bess_k)

    % terms in fraction of Hill algorithm
a0 = 2001.035224;
a1 = 4317.5526;
a2 = 2326;


%         for istp = 1:n_mnts_stp-1
            % determine equivalent concentration of recursive
            % first work out new bessel ratios
%             for kk = 1:n_errs
% 
%                   bes_convl_ntr(kk) = besseli(1,inv_ratio(kk),1)/ ...
%                       besseli(0,inv_ratio(kk),1);
% 
%             end
            % now use approx in Hill to estimate equiv kappa next 
            % maintenance step
%             bess_k = bes_convl_ntr.*bes_1_dtc;

            % term X in best fraction determined by concentration
            gt_95 = bess_k >= 0.95;
            lt_642 = bess_k <= 0.642;
            btw_642_95 = ~gt_95 & ~lt_642;

             if sum(~lt_642) > 0

                 % use continued fraction version in Hill 1981
                 y_Hill = 2./(1-bess_k(~lt_642));
                X_Hill = [];
                if sum(gt_95) > 0
                     X_Hill(gt_95) = 32./(y_Hill(gt_95) - 131.5 + ...
                         120*bess_k(gt_95));
                end

                if sum(btw_642_95) > 0
                     X_Hill(btw_642_95) = a0 + a1*bess_k(btw_642_95) - ...
                         a2*bess_k(btw_642_95).^2;
                end

                fract_term = y_Hill-5-12./(y_Hill-10-X_Hill);
                inv_ratio(~lt_642) = 0.25*(y_Hill + 1 + 3./fract_term);

             end

             if sum(lt_642) > 0 % use Taylor expansion version in Hill 1981

                 R_lt = bess_k(lt_642);
                 R_14_term = R_lt.^11./(1-R_lt).* ...
                     (4.6494 - 5.0797*R_lt + 5.6076*R_lt.^2 - R_lt.^3)/180;
                 inv_ratio(lt_642) = (2*R_lt - R_lt.^3 - R_lt.^5/6 - ...
                     R_lt.^7/24 + R_lt.^9/360 + R_14_term)./(1-R_lt.^2);
                 % + 53*R_lt.^11/2160 
            end
                 %         inv_ratio = 1./(3*bess_k - 4*bess_k.^2 + bess_k.^3); %

%         end
