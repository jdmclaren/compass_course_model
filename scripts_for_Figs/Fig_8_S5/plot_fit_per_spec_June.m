% initialize brewer map colour and marker options
init_brewer_opts

letters = {'a','b','c','d','e','f','g','h','i','j'}; 
%{'Aa','Ab','Ac','Ad','Ae','Af','Ag','Ah','Ai','Aj'}; % {'A','B','C','D','E','F','G','H','I','J'};
FSz = 9;

switch sort_opt_spec_plot
    
    case 1 % goal breadth
        
        [ff,i_fbs] = sort(fr_br(1:9)); %

    case 2 % number steps (error free case)
        
        [ff,i_fbs] = sort(n_hat_fls); % sort(dep_lats); % sort(goal_rads./goal_d); % 
        
    case 3 % goal breadth * sqrt(N steps) as in erf for Normal equiv
        
        [ff,i_fbs] = sort(fr_br(1:9)'.*sqrt(n_hat_fls)); 
        
    case 4 % combined factor meridian breadth (secant(Lat)), N_steps and stepwise distance
        
        [ff,i_fbs] = sort(day_m_d.*n_hat_fls.*geo_fact_gc(1:n_spec)'); % ln_fact); % ./geo_mod_fact
        
end

p_offset = 0;

n_rows = round(sqrt(n_spec));

all_errs = equiv_err_tr_s(1:n_spec:end)*180/pi;

left = 200;
bottom = 200;
width = 400; % 450; % 
height = 360; % 405; % 
figure('Position',[left bottom width height]); %clr_gry); % [0.925 0.975 1]) % [.95 .9 .8])


for i_sp = 1:n_spec 
    
    idx_sp = (0:n_base_errs-1)*n_spec + i_fbs(i_sp);
    
%     all_errs  = fr_br(idx_sp); % .*sqrt(cos_alf.^2.*cos_lat_Arr_2 + sin_alf_2);
%     sig_alf_1 = rte{1}.X
%     
%     e_g_br_1 = best_coeffs{1}(4)
%     
%     all_errs  = fr_br(idx_sp)./ ...
%     (sig_alf_1.*Geo_1.^abs(e_g_br)./(N_0./bes_1_br).^e_N_br).^e_sig_br
%     
%     p_1_geo_N = p_fit{1,sort_mods(1,1)}(idx_sp).* ...
%         geo_Facts(idx_sp).^best_coeffs{1}(4);
%     p_2_geo_N = p_fit{1,sort_mods(2,1)}(idx_sp).* ...
%         geo_Facts(idx_sp).^best_coeffs{2}(4);
    
    subplot(n_rows,n_rows,i_sp)
    hold
    
%     for k_rte = 1:n_rtes


        
        scatter(all_errs  -1.5,p_fit{1,sort_mods(1,1)}(idx_sp),star_sz_1,clr_gg,'h','LineWidth',LW_sm_st) % clr_gg
        h2 = scatter(all_errs -1.5 ,rte{1}.Y(idx_sp),star_sz_2,clr_gg,'h','fill');
        set(h2,'MarkerFaceAlpha',transp_val); % 1%
        
%         plot(fb_all_fit(:,1),p_GGL_avg_geo,'Color',clr_gg,'LineWidth',LW_sm) % ,'LineStyle','--')

        
        scatter(all_errs,p_fit{3,sort_mods(3,1)}(idx_sp),cir_sz_1,clr_tc,'o','LineWidth',LW_sm_st) % clr_tc clr_gr
        h1 = scatter(all_errs ,rte{3}.Y(idx_sp),cir_sz_2,clr_tc,'o','fill');
        set(h1,'MarkerFaceAlpha',transp_val); % 1; %
        
        
%         scatter(all_errs +1.,p_fit{2,sort_mods(2,1)}(idx_sp),dia_sz_1,clr_gm,'d','LineWidth',LW_sm_st) % clr_gr
        h3 = scatter(all_errs +1.,rte{2}.Y(idx_sp),dia_sz_2,clr_gm,'d','fill');
        set(h3,'MarkerFaceAlpha',transp_val); % 1   

        
        if i_sp == n_rows*n_rows -1 % mod(i_sp,n_rows) == 1 && i_sp/n_rows > n_rows -1
            xlabel('Stepwise error (^o)','FontSize',FSz)
        end
        
       if i_sp == n_rows + 1 %mod(i_sp,n_rows) == 1 && i_sp/n_rows > 1
            ylabel('Arrival probability','FontSize',FSz)
        end
        
%         scatter(fr_br_nstp(:,1),p_fit_GGL_st+p_offset,star_sz_1,clr_gr,'h','LineWidth',LW_sm_st) %,'fill')} %clr_gg
%         scatter(fr_br_nstp(:,1),y_GGL_st,star_sz_2,clr_gg,'h','fill') 

%         plot(fb_all_fit(:,1),p_GGL_st_avg_geo,'Color',clr_gg,'LineStyle','--','LineWidth',LW_sm)

      ylim([ 0. 1.02])
      
      if plot_titles_opt
          
          if ~titl_letters
              
             title(species_list{i_fbs(i_sp)})
      
          else
      
%              title({letters{i_sp};' '},'FontSize',FSz)
             title(letters{i_sp},'FontSize',FSz);


          end
          
          ax = gca;
          ax.TitleHorizontalAlignment = 'left';
                       set(gca,'Units','normalized')
            titleHandle = get( gca ,'Title' );
            pos  = get( titleHandle , 'position' );
            pos1 = pos + [0 0.05 0];
            set( titleHandle , 'position' , pos1 );
         
      end
      
      set(gca,'XTick',0:20:60)
      
      xlim([0 det_err_regr(j_err)+2.5])
%         axis([0 0.25 0. 1.02])
%         set(gca,'XTick',XTks,'YTick',0:0.25:1)
    
%     end
    
end