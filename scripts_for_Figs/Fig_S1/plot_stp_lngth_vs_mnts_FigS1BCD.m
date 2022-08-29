plot_opt =  'VM'; %  'Norm';  % 
addpath 'D:\Oldenburg_models\generic_comp_mig_model\circ_stats'
addpath('D:\Oldenburg_models\geomagn_orientn_model\brewer')
cmaps = colormap(brewermap([],'YlOrRd')); % 'Dark2')); % colormap(brewermap(9,'Set1')); 

% option to estimate convolution of several sub-steps within steps
% following approx in Hill 1981 (see below)
% This currently doesn't accomodate drift or inher errors...
Hill_approx_opt = 0; %  1; %   

% set magnitudes of angular deviance to test
dev_dtct_err_degs = [5 10 20 30 40 50 60]; %+err_drft
n_errs = numel(dev_dtct_err_degs);
dev_dtct_errs = dev_dtct_err_degs*pi/180;
dev_mnt_errs = dev_dtct_errs;
dev_trns_errs = dev_dtct_errs;

% include non-compass (drift) error 
dev_drft_err_degs = 0; % [0 15]; % 
dev_drft_errs = dev_drft_err_degs*pi/180;

% inheritance error
dev_inh_err_degs = 0; % 2.5 % [0 2.5 5]; % [0 20];  10 
dev_inh_errs = dev_inh_err_degs*pi/180;

% number maintenance "calibrations" per flight step
n_mnts_stps = [1:10 12:2:50]; % 2.^(0:6); % e.g., once per hour nightly`

% Rayl_fact = sqrt(2-pi/2);
    
% number of migrants 
n_inds = 10000;

Line_W = 1;
FgWd = 250;
FgHt =  225; % 200; %

for i_mnts = 1:numel(n_mnts_stps)
    
    n_mnts_stp = n_mnts_stps(i_mnts);
    

    for j_inh = 1:numel(dev_inh_errs)



        dev_inh_err = dev_inh_errs(j_inh);
        % calculate equivalent Von Mises concentration kappa
        kap_inh = 1/dev_inh_err^2;

        if dev_inh_err > 0

            if strcmp(plot_opt,'VM')
                err_inh =  vmrand(0, kap_inh, [n_inds 1]); 
            else            
                err_inh = randn(n_inds, 1)*dev_inh_err;                     
            end

        else
            err_inh = zeros(n_inds, 1);
        end

        for idr = 1:numel(dev_drft_errs)

            dev_drft_err = dev_drft_errs(idr);   
            % calculate equivalent Von Mises concentration kappa
            kap_drft = 1/dev_drft_err^2;

            if dev_drft_err > 0

                if strcmp(plot_opt,'VM')
                    err_drft =  vmrand(0, kap_drft, [n_inds 1]); 
                else            
                    err_drft = randn(n_inds, 1)*dev_drft_err;                     
                end

            else
                err_drft = zeros(n_inds, 1);
            end

            % simulate stepwise error for each error magnitude
            for i_dtc_err = 1:numel(dev_dtct_errs)

                % set basis deviance and Von Mises conentrations
                dev_i = dev_dtct_errs(i_dtc_err);
                kap_i = 1/dev_i^2;

                % determine initial detection and cue transfer errors
                if strcmp(plot_opt,'VM')
                    err_dtcs = vmrand(0, kap_i, [n_inds 1]); 
                    err_trs = vmrand(err_dtcs, kap_i*ones(n_inds,1), [n_inds 1]);
                else
                    err_dtcs = randn(n_inds, 1)*dev_i; %     
                    err_trs =  err_dtcs  + randn(n_inds, 1)*dev_i; %
                end

                % initialize x and y components (first sub-step) 
                % i.e., before cuemaintenance
                % Assume err = 0 is due South; here is -Y direction
                x_ntrs = -sin(err_dtcs+err_inh+err_drft);
                y_ntrs = -cos(err_dtcs+err_inh+err_drft);
                x_trs = -sin(err_trs+err_inh+err_drft);
                y_trs = -cos(err_trs+err_inh+err_drft);

                % initialize in-flight errors
        %         err_m_ntrs = err_dtcs; % zeros(n_inds,1);
        %         err_m_trs = err_trs; % zeros(n_inds,1); % 

                % simulate maintenance errors based on stepwise preferred heaidng
                for im = 1:n_mnts_stp-1

                    % for non-transferred movement, cue maintenance is unbiased 
                    if strcmp(plot_opt,'VM')
                        err_m_ntrs = vmrand(0, kap_i, [n_inds 1]);
                    else
                        err_m_ntrs = randn(n_inds, 1)*dev_i; % err_m_ntrs  + 
                    end
                    % for maintenace with transferred movement, 
                    % preferred direction is biased by initial detection error
                    if strcmp(plot_opt,'VM')
                        err_m_trs = vmrand(err_trs, kap_i*ones(n_inds,1), [n_inds 1]); 
                    else
                        err_m_trs = err_trs + randn(n_inds, 1)*dev_i;
                    end

                    x_ntrs = x_ntrs - sin(err_m_ntrs+err_inh+err_drft);
                    y_ntrs = y_ntrs - cos(err_m_ntrs+err_inh+err_drft);
                    x_trs = x_trs - sin(err_m_trs+err_inh+err_drft);
                    y_trs = y_trs - cos(err_m_trs+err_inh+err_drft);

                end

                % add any non-compass (drift) component
                % assume uniform through flight step
    %             if idr > 1
    %                 x_ntrs = x_ntrs - n_mnts_stp*sin(err_drft);
    %                 y_ntrs = y_ntrs - n_mnts_stp*cos(err_drft);
    %                 x_trs = x_trs - n_mnts_stp*sin(err_drft);
    %                 y_trs = y_trs - n_mnts_stp*cos(err_drft);
    %             end

                % sum stepwsie errors including any non-compass "drift" error
        %         err_stp_ntrs = mod(err_m_ntrs + err_drft +pi,2*pi) - pi;            
        %         err_stp_trs = mod(err_m_trs + err_drft +pi,2*pi) - pi; 

                % determine overall flight angle (angular error)
                err_stp_ntrs = atan2(-x_ntrs,-y_ntrs); % y_ntrs./x_ntrs; % 
                err_stp_trs = atan2(-x_trs,-y_trs); % y_trs./x_trs; % 

                % calculate total circ std dev for nontr and tr cases
                [~, s_ntr(i_dtc_err)] = circ_std(err_stp_ntrs);
                [~, s_tr(i_dtc_err)] = circ_std(err_stp_trs);  

                r_ntr(i_dtc_err) = circ_r(err_stp_ntrs);
                r_tr(i_dtc_err) = circ_r(err_stp_trs);

                k_sim_ntr(i_dtc_err) = circ_kappa(err_stp_ntrs);
                k_sim_tr(i_dtc_err) = circ_kappa(err_stp_trs);

        %         s_ntr(i_dtc_err) = circ_std(x_ntrs./y_ntrs);
        %         s_tr(i_dtc_err) = circ_std(x_ntrs./y_trs);  
        %         [~,s_ntr(i_dtc_err)] = circ_std(err_stp_ntrs);
        %         [~,s_tr(i_dtc_err)] = circ_std(err_stp_trs);  
        %         s_ntr(i_dtc_err) = std(err_stp_ntrs);
        %         s_tr(i_dtc_err) = std(err_stp_trs);  
        
                stp_lngth(i_mnts,i_dtc_err) = ...
                    mean(sqrt((x_ntrs.^2+y_ntrs.^2)))/n_mnts_stp;
                
                stp_err(i_mnts,i_dtc_err) = s_ntr(i_dtc_err)*180/pi;


            end

            % calculate Normal approx to compare with stepwise summations of
            %  Von Mises errors

            equiv_N_errs_nontr_s(i_mnts,:) = sqrt(dev_dtct_errs.^2/n_mnts_stp + ...
                     dev_drft_err^2 + dev_inh_err^2); % Rayl_fact* .*sin_ss_ntr

            if n_mnts_stp > 1

                 equiv_N_errs_tr_s(i_mnts,:) = sqrt(dev_dtct_errs.^2 + dev_trns_errs.^2 ...
                   + dev_drft_err^2 + ...
                  dev_mnt_errs.^2/(n_mnts_stp-1) + dev_inh_err^2); 
    %                (n_mnts_stp-1)*dev_mnt_errs.^2/n_mnts_stp^2 + dev_inh_err^2); % Rayl_fact* .*sin_ss_tr

            else

                 equiv_N_errs_tr_s(i_mnts,:) = sqrt((dev_dtct_errs.^2 + dev_trns_errs.^2) ...
                    + dev_drft_err^2 + dev_inh_err^2); %          .*sin_ss_tr

            end

            % account for sin alf approx using ratio of mod bessel fn 2nd order
            kap_dtc = 1./dev_dtct_errs.^2; % 
            kap_eq_N_ntr = 1./equiv_N_errs_nontr_s(i_mnts,:).^2;
            kap_eq_N_tr = 1./equiv_N_errs_tr_s(i_mnts,:).^2;   

            for kk = 1:n_errs
              bes_1_dtc(kk) = besseli(1,kap_dtc(kk),1)/besseli(0,kap_dtc(kk),1);
              bes_1_ntr(kk) = besseli(1,kap_eq_N_ntr(kk),1)/besseli(0,kap_eq_N_ntr(kk),1);
              bes_2_ntr(kk) = besseli(2,kap_eq_N_ntr(kk),1)/besseli(0,kap_eq_N_ntr(kk),1);
              sin_ss_ntr(kk) = 0.5*(1-bes_2_ntr(kk));
              cos_ss_ntr(kk) = 0.5*(1+bes_2_ntr(kk)-2*bes_1_ntr(kk).^2);
              bes_2_tr(kk) = besseli(2,kap_eq_N_tr(kk),1)/besseli(0,kap_eq_N_tr(kk),1);
              bes_1_tr(kk) = besseli(1,kap_eq_N_tr(kk),1)/besseli(0,kap_eq_N_tr(kk),1);
              sin_ss_tr(kk) = 0.5*(1-bes_2_tr(kk));
              cos_ss_tr(kk) = 0.5*(1+bes_2_tr(kk)-2*bes_1_tr(kk).^2);
              bes_k_sim_ntr(kk) = besseli(1,k_sim_ntr(kk),1)/besseli(0,k_sim_ntr(kk),1);
            end         

        %     equiv_N_errs_tr_s = sqrt(dev_dtct_errs.^2 + dev_trns_errs.^2 + ...
        %         (n_mnts_stp-1)*(dev_mnt_errs.^2 + dev_drft_err^2)); % 


              kap_div_n = kap_dtc*n_mnts_stp;
              for kk = 1:n_errs
                  bes_div_n(kk) = besseli(1,kap_div_n(kk),1)/ ...
                      besseli(0,kap_div_n(kk),1);
              end
              s_div_n = sqrt(-2*log(bes_div_n)); %  sqrt(2*(1-bes_div_n)); % 
    %           plot(s_div_n*180/pi,s_ntr*180/pi,'.-m','LineWidth',Line_W)

              s_dtc = sqrt(-2*log(bes_1_dtc)); %  sqrt(2*(1-bes_1_dtc)); % 
    %           plot(s_dtc*180/pi,s_ntr*180/pi,'.-g','LineWidth',Line_W)

        %     plot(equiv_N_errs_tr_Hill*180/pi,s_ntr*180/pi,'--','LineWidth',Line_W)

%             sc1 = scatter(equiv_N_errs_tr_s*180/pi,s_tr*180/pi,100, ...
%                dev_dtct_err_degs,'>','fill');
% 
%             set(sc1,'MarkerEdgeColor','k') % LineWidth',Line_W)
% 
%             sc2 = scatter(equiv_N_errs_nontr_s*180/pi,s_ntr*180/pi,60, ...
%                 dev_dtct_err_degs,'o','fill');
% 
%             set(sc2,'MarkerEdgeColor','k')

        %     plot(equiv_N_errs_tr_s*180/pi,s_tr*180/pi,'LineWidth',Line_W)   

        %     plot(sin_ss_ntr*equiv_N_errs_nontr_s*180/pi,s_ntr*180/pi,'--','LineWidth',Line_W)
        %     plot(sin_ss_tr*equiv_N_errs_tr_s*180/pi,s_tr*180/pi,'LineWidth',Line_W)

        %     plot(dev_dtct_err_degs,s_ntr./equiv_N_errs_nontr_s,'--','LineWidth',Line_W)
        %     plot(dev_dtct_err_degs,s_tr./equiv_N_errs_tr_s,'LineWidth',Line_W)  
        %     plot(dev_dtct_err_degs,sin_ss_ntr,':','LineWidth',Line_W)
        %     plot(dev_dtct_err_degs,sin_ss_tr,'.-','LineWidth',Line_W)  

        end

    %     title({'simulated vs. expected angular error,';[plot_opt ' distribution']})
        % legend({'no drift, no transfer','no drift, transfer', ...
        %     'drift, no transfer','drift and transfer'},'AutoUpdate','off')
    %     plot([0 90],[0 90],':k')
    %     axis([0 90 0 90])




    end


%     figure('Position',[700 400 FgWd FgHt]);
%     cb = colorbar; % ('North');
%     title(cb,['Error components $ (^o)$'], ...
%         'Interpreter','latex','FontSize',9) % (cb,{'Error components (^o)'})
%     caxis([0 60])
%     colormap(cmaps)
%     axis('off')
%     box('off')
    
end

% now to compare all three measures for given conc R
H3 = figure('Position',[400 400 FgWd FgHt]);
hold
title({'D',' '},'FontSize',9)
ax = gca;
ax.TitleHorizontalAlignment = 'left';

%     plot(n_mnts_stps,bes_1_dtc(1)*ones(size(n_mnts_stps)),'--k')
for i_err = 1:numel(dev_dtct_errs)
    sc3 = scatter(n_mnts_stps,smooth(stp_lngth(:,i_err)),35, ...
        dev_dtct_err_degs(i_err)*ones(size(n_mnts_stps)),'fill');
         set(sc3,'MarkerEdgeColor','k') % Li'MarkerEdgeColor','k')
    
%     plot(n_mnts_stps,bes_1_dtc(i_err)*ones(size(n_mnts_stps)),'--k')
end
    colormap(cmaps)
%     plot([0 115],[0 115],':k')
    axis([1 52 floor(10*bes_1_dtc(i_err))/10 1])
    set(gca,'XTick',[1 10:10:50],'YTick',0:0.1:1)
    caxis([0 60])

    set(gca,'YAxisLocation','Right')
%     colormap(cmaps)
       
    ylabel('Relative flight distance', 'Interpreter','latex','FontSize',11)
    xlabel('Number of sub-steps', 'Interpreter','latex','FontSize',11)
     set(gca,'YAxisLocation','Right')
     
     % now to compare all three measures for given conc R
H4 = figure('Position',[400 400 FgWd FgHt]);
hold
title({'C',' '},'FontSize',9)
ax = gca;
ax.TitleHorizontalAlignment = 'left';

%     plot(n_mnts_stps,bes_1_dtc(1)*ones(size(n_mnts_stps)),'--k')
for i_err = 1:numel(dev_dtct_errs)
    sc3 = scatter(n_mnts_stps,smooth(stp_err(:,i_err)),35, ...
        dev_dtct_err_degs(i_err)*ones(size(n_mnts_stps)),'fill');
         set(sc3,'MarkerEdgeColor','k') % Li'MarkerEdgeColor','k')
    
%     plot(n_mnts_stps,bes_1_dtc(i_err)*ones(size(n_mnts_stps)),'--k')
end
    colormap(cmaps)
%     plot([0 115],[0 115],':k')
    axis([1 52 0 61])
    set(gca,'XTick',[1 10:10:50],'YTick',0 :10:60)
    caxis([0 60])

%     set(gca,'YAxisLocation','Right')
%     colormap(cmaps)
       
    ylabel('Stepwise error $ (^o)$', 'Interpreter','latex','FontSize',11)
    xlabel('Number of sub-steps', 'Interpreter','latex','FontSize',11)
     set(gca,'YAxisLocation','Left')
     
     
     figure('Position',[700 200 FgWd 2*FgHt]);
    cb = colorbar; % ('North');
    title(cb,['Error components $ (^o)$'], ...
        'Interpreter','latex','FontSize',9) % (cb,{'Error components (^o)'})
    caxis([0 60])
    colormap(cmaps)
    axis('off')
    box('off')
       
