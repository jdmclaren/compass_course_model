% initialize brewer map colour and marker options
init_brewer_opts

p_offset = 0;

n_rows = round(n_base_errs/2);

% rep_arr_lats = trans_array(repmat(arr_lats,[n_base_errs 1]));

figure

for i_err = 1:n_base_errs 
    
    idx_err = (i_err-1)*n_spec + (1:n_spec);
    
    fr_br_geo_Lox = fr_br(idx_err).*sqrt(cos(az_lox).^2.*cosd(arr_lats).^2 + ...
        sin(az_lox).^2)';
    
    fr_br_geo_TCSC = fr_br(idx_err).*sqrt(cos(az_gc_mn).^2.*cosd(arr_lats).^2 + ...
        sin(az_gc_mn).^2)';
    
%     sig_alf_1 = rte{1}.X
%     
%     e_g_br_1 = best_coeffs{1}(4)
%     
%     fr_br_geo = fr_br(idx_err)./ ...
%     (sig_alf_1.*Geo_1.^abs(e_g_br)./(N_0./bes_1_br).^e_N_br).^e_sig_br
%     
%     p_1_geo_N = p_fit{1,sort_mods(1,1)}(idx_err).* ...
%         geo_Facts(idx_err).^best_coeffs{1}(4);
%     p_2_geo_N = p_fit{1,sort_mods(2,1)}(idx_err).* ...
%         geo_Facts(idx_err).^best_coeffs{2}(4);
    
    subplot(2,n_rows,i_err)
    hold
    
%     for k_rte = 1:n_rtes
        
        scatter(fr_br_geo_Lox,p_fit{1,sort_mods(1,1)}(idx_err),dia_sz_1,clr_gr,'d','LineWidth',LW_sm_st) % clr_gg
        scatter(fr_br_geo_Lox,rte{1}.Y(idx_err),dia_sz_2,clr_gg,'d','fill')
% 
%         plot(fb_all_fit(:,1),p_GGL_avg_geo,'Color',clr_gg,'LineWidth',LW_sm) % ,'LineStyle','--')

        scatter(fr_br_geo_Lox,p_fit{2,sort_mods(2,1)}(idx_err),dia_sz_1,clr_gr,'d','LineWidth',LW_sm_st) % clr_gg
        scatter(fr_br_geo_Lox,rte{2}.Y(idx_err),dia_sz_2,clr_gg,'d','fill')
        
        scatter(fr_br_geo_TCSC,p_fit{3,sort_mods(3,1)}(idx_err),star_sz_1,clr_gr,'h','LineWidth',LW_sm_st) % clr_gg
        scatter(fr_br_geo_TCSC,rte{3}.Y(idx_err),star_sz_2,clr_tc,'h','fill')
        
%         scatter(fr_br_nstp(:,1),p_fit_GGL_st+p_offset,star_sz_1,clr_gr,'h','LineWidth',LW_sm_st) %,'fill')} %clr_gg
%         scatter(fr_br_nstp(:,1),y_GGL_st,star_sz_2,clr_gg,'h','fill') 

%         plot(fb_all_fit(:,1),p_GGL_st_avg_geo,'Color',clr_gg,'LineStyle','--','LineWidth',LW_sm)

      ylim([ 0. 1.02])
%         axis([0 0.25 0. 1.02])
%         set(gca,'XTick',XTks,'YTick',0:0.25:1)
    
%     end
    
end