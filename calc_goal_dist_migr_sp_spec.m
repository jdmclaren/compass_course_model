% if ~exist('species_list')
    
   species_list = ["Nathusius","Kirtlands Warbler","Ring Ouzel","Monarch", ...
     "Hoopoe","BU Rosefinch","Finn Marsh Warbler","GreyCheekThrush","Sib Will Warb S Hem"];
 % ,"RF Bluetail RU"]; %  , ...     "Siberian Willow Warbler"
  
% end
% "Sib Will Warb S Hem" % 'France Marsh Warbler',
 
 R_Earth_km = 6371;
 
 plot_var_test = 3;
 
 clear specs goal_rads  goal_d goal_d_lox r_lox  ...
     day_m_d overall_mig_sp front_brdth  ...
     front_brdth_adj front_brdth_adj_lox  ...
     n_hat_fls dep_lats dLon_goals  ...
     dLat_goals mn_lats ln_fact az_gc_1 az_gc_2 az_lox ...
      n_hat_fls_min n_hat_fls_max az_gc_mn max_n_sim sig_date stop_durs
 
 for i_spec = 1:numel(species_list)
     
     species = species_list{i_spec};
     specs{i_spec} = species;
     get_species_params;
     
     goal_d(i_spec) = distance('gc',dep_lat_sp_degs,dep_lon_sp_degs, ...
         lat_goal_centr,lon_goal_centr,'degrees')*pi/180*R_Earth_km;
     goal_d_lox(i_spec) = distance('rh',dep_lat_sp_degs,dep_lon_sp_degs, ...
         lat_goal_centr,lon_goal_centr,'degrees')*pi/180*R_Earth_km;   

     az_gc_1(i_spec) = azimuth('gc',dep_lat_sp_degs,dep_lon_sp_degs, ...
         lat_goal_centr,lon_goal_centr,'degrees')*pi/180 - pi;
     
     az_gc_2(i_spec) = mod(pi + azimuth('gc',lat_goal_centr,lon_goal_centr, ...
         dep_lat_sp_degs,dep_lon_sp_degs,'degrees')*pi/180,2*pi) - pi;  
     
     az_gc_mn(i_spec) = circ_mean([az_gc_1(i_spec) az_gc_2(i_spec)]');
     
     az_lox(i_spec) = azimuth('rh',dep_lat_sp_degs,dep_lon_sp_degs, ...
         lat_goal_centr,lon_goal_centr,'degrees')*pi/180 - pi;   
     
     r_lox(i_spec) = goal_d_lox(i_spec)/goal_d(i_spec);
     
     day_m_d(i_spec) = (Va_mps*n_hs*3.6); % *n_fl_seq/(n_fl_seq+stop_dur)
     
     overall_mig_sp(i_spec) = (Va_mps*n_hs*3.6)*n_fl_seq/(n_fl_seq+stop_dur);
     
     front_brdth(i_spec) = goal_rad/goal_d(i_spec);   
     
     goal_rads(i_spec) = goal_rad;
     
     n_hat_fls(i_spec) = (goal_d(i_spec)-goal_rad)/day_m_d(i_spec);
     n_hat_fls_min(i_spec) = (n_hat_fls(i_spec));    
     n_hat_fls_max(i_spec) = round((max_t_bio+stop_dur)/(1+stop_dur/n_fl_seq));  
     
     max_n_sim(i_spec) = max_n;
     
     front_brdth_adj(i_spec) = front_brdth(i_spec)*sqrt(n_hat_fls(i_spec));
     front_brdth_adj_lox(i_spec) = front_brdth_adj(i_spec)*r_lox(i_spec);
     
     dep_lats(i_spec) = dep_lat_sp_degs;
     arr_lats(i_spec) = lat_goal_centr;
     
     dLon_goals(i_spec) = abs(dep_lon_sp_degs-lon_goal_centr);
     dLat_goals(i_spec) = abs(dep_lat_sp_degs-lat_goal_centr);
     
     mn_lats(i_spec) = (dep_lat_sp_degs+lat_goal_centr)/2;
     
     n_mnts(i_spec) = n_hs;
     
     max_ns(i_spec) = max_n;
     
     if lat_goal_centr > 0
         
        ln_fact(i_spec) = 180/pi*(log(sec(dep_lat_sp_degs*pi/180) + ...
         tan(dep_lat_sp_degs*pi/180)) - log(sec(lat_goal_centr*pi/180) + ...
         tan(abs(lat_goal_centr*pi/180))))*cos(lat_goal_centr*pi/180)/ ...
         (dep_lat_sp_degs - lat_goal_centr);
     
     else
         
         ln_fact(i_spec) = 180/pi*(log(sec(dep_lat_sp_degs*pi/180) + ...
         tan(dep_lat_sp_degs*pi/180)) + log(sec(lat_goal_centr*pi/180) + ...
         tan(abs(lat_goal_centr*pi/180))))*cos(lat_goal_centr*pi/180)/ ...
         (dep_lat_sp_degs - lat_goal_centr);
     
     end
     
%           front_brdth_adj_lat(i_spec) = front_brdth_adj(i_spec)*r_lox(i_spec);        
        
       sig_date(i_spec) = sqrt(std_day^2 + std_n_stop^2);
       
       stop_durs(i_spec) = stop_dur;

 end
     
 save('dist_speed_brdth_frnts_per_sp','specs', ...
    'goal_rads','goal_d','goal_d_lox','r_lox', ...
     'day_m_d','overall_mig_sp','front_brdth', ...
     'front_brdth_adj','front_brdth_adj_lox', ...
     'n_hat_fls','dep_lats','dLon_goals', ...
     'dLat_goals','mn_lats','ln_fact', ...
     'az_gc_1','az_gc_2','az_gc_mn','az_lox', ...
     'n_hat_fls_min','n_hat_fls_max','max_n_sim', ...
     'sig_date','stop_durs')
     
