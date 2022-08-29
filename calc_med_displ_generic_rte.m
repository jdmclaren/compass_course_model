function med_d_sgn = calc_med_displ_comp_rte(alph_0,species, ...
    or_progs,tc_clock_az,err_dets,ar_error_cal,ar_break_stop_opt, ...
    default_inh_err,plot_var_test,init_geogr_head,incl_decl, ...
    star_night,dalf_close,n_days_close,opt_model,dep_lat_degs,lat_dist, ...
    lon_dist,max_t,n_fl_seq,stop_dur,std_n_stop,n_hs,Va_mps,dep_date,std_day,n_inds)

%     lat_dist = 67.5;
%     max_t = 60; %  
%     n_fl_seq =   5; % 60 % 
%     stop_dur =  5; % 10; % 2; %
%     std_n_stop = 2;
%     n_hs = 8; % (2:2:10);
%     Va_mps = 12.5; %3 % 10 % 
%     alfs_deg_or_prs = [56 90 95 120]; % 90 % -120
%         dep_date = (plot_var_test==3)*[1 9 2000] + (plot_var_test==5)*[15 8 2000];
%     std_day = 14*(plot_var_test == 3); %
%              err_dets = 15;

% choose orientation programs
% 1 = magn loxodrome,
% 2 = linear compass (slope = 1 deg/deg inclin)
% 3 = magnetoclinic (transverse),
% 4 = magnetoclinic parallel
% 5 = vertical field compass (slope = 1 per sin(degree) lat)
% 6 = fixed sun comp
% 7 = time-comp sun comp 
% 8  = magnclinic trasverse projn of hz compt fixed
% 9 = mixed: half way between sun azimuth and magnetic (here geogr) offsets
% 10 = mixed: halfway between tc sun and magn lox
% 11 = TC Sun comp reset each stopover
% 12 = mixed: halfway between reset tc sun and magn lox
% 13 great circle
% 14 = fixed sun compass reset

% divided alpha_0 by 100 in algorithm
% to help derivative based root finder
% alph_0 = alph_0_100;

or_prog_test = [1 3 6 11]; %  7  12 [1 3] %  
        
% or_prog_idx = 1:4 %   4 % [1 2 4]; %    
% or_progs = or_prog_test(or_prog_idx);

full_run = false; %   true % 
plot_ind_tracks_opt =  false; % true %     
plot_spec_map_opt =  false; % true %
plot_close_dist_opt = false; % true %  
plot_opt = false;
species_goal_opt = true;

% vary or plot headings (3) or departure date (5) or inherited / 
% initial heading (6)
% plot_var_test =  3; %  5 %
% plot_var_disps = [6 6 5 5]; % 6 % 3 %  5 % 
% plot_var_disps = plot_var_disps(or_prog_idx);

% 0 = none, 1 = end of or progs per sp, 2 = each ap
% plot_cbr_opt =  2; % 1 % 0 % 

get_species_params

 dLon_goal = lon_goal_centr - dep_lon_sp_degs;
 % offset for plotting from generic 180 deg departure in run_no_topo_any
spec_lon_offset = dep_lon_sp_degs-180;

% choose headings from specific or prog chosen
alf_init =  alph_0; % alfs_deg_or_prs(or_prog_idx);

%% run simulation and exclude extreme (>8000 km) distances
% (typically can occur with high lat magnetoclinic simulations)
run_no_topo_any

% if all(abs(d_cl_sgn)> 6000)
% 
%     med_d_sgn = 6000; % abs();
% 
% else
    
if strcmp(opt_model,'fzero')
    
%     if all(abs(d_cl_sgn)>8000)
% 0

%         med_d_sgn = 8000; % abs();
% 
%     else

        med_d_sgn = median(d_cl_sgn); % abs();   

%     end
    
else % fminbnd
    
        med_d_sgn = median(d_close{ia,i_err}); % (d_close{ia,i_err}));
    
end

% end
