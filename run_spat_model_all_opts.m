clear

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

or_prog_test = [1 3 6 11]; %  7  12 [1 3] %  

or_prog_test_sp = [1 1 3 6 11 7 7 12];
        
or_prog_idx = 2; % 3; % 4; %    sum(~is_nan) 4 % 1:4 %   [1 2 4]; %    
or_progs = or_prog_test(or_prog_idx)

full_run =  false; %  true %  
plot_ind_tracks_opt =  true ;% false %   true ;% ~full_run; %  
plot_spec_map_opt =  ~full_run &&  ~plot_ind_tracks_opt % false; % 
plot_close_dist_opt =  false; % true %  
plot_opt = plot_spec_map_opt || plot_ind_tracks_opt;
species_goal_opt = ~full_run && ~plot_ind_tracks_opt;

% vary or plot headings (3) or departure date (5) or inherited / 
% initial heading (6)
plot_var_test =  3; %  5 %
plot_var_disps = [6 6 5 5]; % 6 % 3 %  5 % 
plot_var_disps = plot_var_disps(or_prog_idx);

% 0 = none, 1 = end of or progs per sp, 2 = each ap
plot_cbr_opt =  2; % 1 % 0 % 

% species_list = {'Kirtlands Warbler', 'Nathusius','Monarch', ...
%      'Marsh Warbler','GreyCheekThrush','Siberian Willow Warbler'}; 
% species_list = {'SEPacific Humpback'} % 'Marsh Warbler'}; % 'Spotted Flycatcher'};
species_list =  {'Siberian Willow Warbler'}; %,{'Alk Wheatear'};  % }; %'SEPacific Humpback','BluefinTuna'
% species_list = {} % 'Kirtlands Warbler'} % ,Monarch'} 'BluefinTuna'} % ,
%      'Siberian Willow Warbler','BluefinTuna','GreyCheekThrush'};
% species_list = {'GreyCheekThrush NF'} % 'Siberian Willow Warbler','Monarch','GreyCheekThrush','Kirtlands Warbler'}  %'Kirtlands Warbler'} % '} %'Siberian Willow Warbler'}; %  
% 'Whinchat' % 'Blackcap' % 'SEPacific Humpback',

% Jan 2021 added autocorr error - here the corr coeff for first-order
%  (nightly) lag in detectn error
ar_error_det =  0; %   0.25; %
% option to break up autocorrelation after stopover
% eg if related to weather or finsihed detour
% (less so for magnetoic anomaly unless somehow accounted for using other
% cues)
ar_break_stop_opt = true; % false; 

% new Jan 2021 tc clock can either involve local or reference
% sun azimuth adjustments (I believe local makes most sense, i.e., TCSC
% bird adjusts its azimuth according to reference clock but local rate of
% change of sun az
tc_clock_az = 'reference'; % 'local'; %   
% init_geogr_head = true; %
sun_inher_opt = false; % true; %  

% inh magnetoclines is 1) of magn param
% i.e. transverse projection of inclination gamma
% and (2) is geographic (defauilt)
magcl_inh_opt = 2; % 1; % 

% For option full run,
% the variable "err_sources" specifies which of three error sources to vary (all in degrees):
% err_sources=1: vary "err_dets" = circular standard error in stepwise headings 
%   (imprecision in compass "detection" or displacement, e.g. by wind)
% err_sources=2: vary "d_alf_init_degs" = circular standard deviation in initial (inherited) heading 
%   (variable inheritability, natal dispersal = where they return to breed)
% err_sources=3 vary "d_alf_init_degs" = variability in initial latitude 
%   (natal dispersal or post natal dispersal = where they wander after fledging) 

% For the other options, (plotting generic or given species trajectories),
% we will focus on testing (varying) err_dets, i.e., imprecision in stepwise compass headings;

% We also can specify "default_inh_err" and "default_det_err" to include 
% (a single value for) error in the initial headings or stepwise error among individuals
% (I don't vary the start latitude when testing detection or inheritance;
% that is better done in my other (spatial / evolutionary) model which
% incorporates natal dispersal and generational inheritance)

% choose which analysis / simulation to run
if full_run % plot_ind_tracks_opt

    ar_err_mnt = 0; % 0.75; % 0.25; % 
    is_trans_comp = true; %  false; %   
    
    ar_err_hr_drft = 0.75;
    ar_error_drft = 0.25;
    
    % nonstop option negates stopover and flight squence in
    % run_no_topo_any
        non_stop =  false; % true; %  
        err_sources =  1; %2 % 3 %
        default_det_err = 5;
        default_inh_err = 0; %
        
        % number ind migrs per angle
        n_inds = 10000; % 10 % 5000; % 
        one_vec = ones(1,n_inds);
        % error values to test
        err_dets = 0:60; % 0:2.5:60; % 0:5:30; %  [0 15 30 45]; % 
        err_trans = is_trans_comp*err_dets; % zeros(size(err_dets)); % 
        err_mnts = err_dets; % zeros(size(err_dets); 20*ones(size_err_dets); % 
        err_drfts =  15*ones(size(err_dets)); %zeros(size(err_dets)); % 
        std_err_inits =  0:0.5:10; % 1 % 2.5 % 0; % [-2.5 2.5]; % 0; %
        d_alf_init_degs = 0:0.5:15;      

        dep_lat_degs =  45 % 65 %  80 % 
        lat_dist =  20 % 65 % 
        dep_date =  [15 9 2000]; %  [31 10 2000]; % [1 8 2000]; % 
        std_day = 5; % 10; % 
        max_n = 40; % 60 % 
        n_fl_seq = 5;
        n_hs = 8; % (2:2:10);
        n_mnts = n_hs; %
        stop_dur = 5; % 10 % 12
        std_n_stop = 2;
        Va_mps = 12.5; %3 % 10 % 
%          if is_sun_comp || or_prog == 13 

        alfs_deg_full_sun = -180:180; %  -145:5:145; % [-60 -30 0 30 60 90 135 165 195 220 250] %[-60 60 100 170] %  -25:2.5:177.5; %-30:5:145 %  60 % 0:5:175 %   -120:30:120 %  
%         else

        alfs_deg_full_mag =  0:90; %0:5:90; %   [30 85 90]; %[0:2.5:87.5 89.9:0.01:89.99]; %[0:10:80 85:89 89.9999] %  0:2.5:90; % [0 60] %

%          end   
         
         run_no_topo_any
         
         main_plots_opts = [0 0 1 0 1];
         plot_scatter_opt = false;
         plot_output_spat_mod

elseif plot_spec_map_opt
    
    err_mnt =  20; % 0; % 
    ar_err_mnt = 0.75; % 0.25; % 
    trans_comp = false; % true; % 

    err_sources = 1;  % 2 % 3 % 
    % number ind migrs per angle
    n_inds = 10000; % 10 % 300 % 
    one_vec = ones(1,n_inds);
    
    % for simulation (pre mapping) init Lon is arbitraray but based on 180
    dep_lon_degs = 180;
    
    % 15 %  [0 5 10 15 25 35 50] %  0:2.5:60; % 10 %  
    std_err_inits = 0; % degs Lat offset 1 % 2.5 % 0; % [-2.5 2.5]; % 0; %
    err_dets =  20; %[0 15 30]  %
    d_alf_init_degs = NaN;
    std_err_inits = NaN;
    default_inh_err = 2.5 % 5 % 0 %    
    
%     plot_titl = false; % 
    plot_title = false;
    save_figs_opt = false;
    plot_colbr = false; % true %
    plot_cb_end = false; % true %
    clr_perf_trx = 'm'; % 'w'; %

    std_day = 5*(plot_var_test == 3); % 14 %
            
    % transparanecy of trajectory markers on map (0 invisible, 1 solid)
    transp_val = 0.5; % 1 % 
    
    plot_var_disp = 3; % display trajs with color coded headings
    map_proj = 'stereo'; % 'Mercator'; 
    % test for tighter schedule
    % i.e. close number days within mean +/- std 
    % eg if std = 10: try 14, std = 5: try 10
    n_days_close = std_day;
    dalf_close = default_inh_err;

        
    for isp = 1:numel(species_list)
        
        species = species_list{isp};
        
        get_species_params
        n_mnts = n_hs; %
        get_species_init_alphs
        alfs_deg_or_prs = alf_inits;
          dLon_goal = lon_goal_centr - dep_lon_degs;
         % offset for plotting from generic 180 deg departure in run_no_topo_any
        spec_lon_offset = dep_lon_degs-180;

        
        for i_prog = 1:numel(or_progs)
            
            or_prog = or_progs(i_prog);
            
            
            or_init_idx = find(or_prog == or_prog_test_sp);
    
            if numel(or_init_idx) > 1

                or_init_idx = 2;

            end
            alf_init = alf_inits(or_init_idx);
       
            % choose headings from specific or prog chosen
            alfs_deg_or_prs  = alfs_deg_or_prs(or_prog_idx(i_prog));

            incl_decl = or_prog == 3 || or_prog == 1;

            run_no_topo_any
            
        end
        
    end
    
else % plot ind trajs for generic species

    % drift errs auto set to zero
%     err_drfts = 0;
%     ar_err_hr_drft = 0;
%     ar_error_drft = 0;
    
    species = 'generic';
    
    err_sources = 1;  % 2 % 3 % 
    n_inds = 10000; % 10 % 
    default_inh_err = 0;
    dalf_close = default_inh_err;
     
    % stepwsie flight hours
    n_hs = 8; % (2:2:10);
    Va_mps = 12.5; %3 % 10 % 
    
    one_vec = ones(1,n_inds);
    % error values to test
    err_dets = 0; % :2.5:60 %  0:5:30; % 0:.5:1.5; %  10 %  0 %10 %
    
    % turn off auto-corr error option if testing perfect tracks
    if numel(err_dets) == 1 && err_dets == 0
        ar_error_det = 0;
        ar_break_stop_opt = false;
    end
    
    std_err_inits = 0; % 0:0.5:10; % 1 % 2.5 % 0; % [-2.5 2.5]; % 0; %
    d_alf_init_degs = 0; % 0:15;      

    dep_lat_degs = 65; %  50 %  80 %
    % dep lon is arbirary (set to 180 since that used independently elsewhere)
    dep_lon_degs = 180;
    
    % specify which lat and long distance desired
    lat_dist =  dep_lat_degs %   35; %60 % 
    lon_dist = 90; %-10; % 10; %   
    
    % standard options for inheritance of geographic vs. sun compass or
    % magnetic headings
    init_geogr_head = true;
    
    % default nightly headings follow stars (geogr course) rather than
    % magnetic loxodrome (gemagnetic compass)
    % Independant of inheritance and detection, either are possible
    star_night = true;
        
    % dates only relevant for sun compass or if we include declination
    incl_decl = false;
    dep_date =  [1 8 2000]; % [1 10 2000]; %  [15 9 2000];
    std_day = 0; % 5; % 
    n_days_close = std_day;
    max_n = 40; % 60 % 
    % stopover only relevant with sun compass (azimuth will change)
    % Number of flights befoe extended stopover
    n_fl_seq = 100; % 5;
    stop_dur = 0; % 5; % 12
    std_n_stop = 0; %  2;
        
    prompt = 'do you want to enter initial guess? 0  = no, 1 = yes';
    input_alph_str = input(prompt,'s');
    input_alph = str2double(input_alph_str);
    
    if input_alph ==1
        
        prompt = 'enter init heading vs. South (degrees)';
        alpha_0_str = input(prompt,'s');
        alpha_0 = str2num(alpha_0_str);
      
    else
        % try an initial rhumbline 
        disp('init heading (degs) vs. South is')
        if or_prog_idx < 3
            alpha_0 = azimuth('rh',dep_lat_degs,dep_lon_degs, ...
            dep_lat_degs-lat_dist,dep_lon_degs-lon_dist) - 180
        elseif or_prog_idx == 4 % || lon_dist > 0
            alpha_0 = azimuth('gc',dep_lat_degs,dep_lon_degs, ...
            dep_lat_degs-lat_dist,dep_lon_degs-lon_dist) - 180
        else
             alpha_0 = (azimuth('rh',dep_lat_degs,dep_lon_degs, ...
            dep_lat_degs-lat_dist,dep_lon_degs-lon_dist) - 180); % /2
        end
    end
    
    
    opt_mod = 'fzero'; % 'fmnbnd'; % 
%     opt_mod = 
    
    fun = @(x) calc_med_displ_generic_rte(x,species,or_progs,tc_clock_az, ...
            err_dets,ar_error_det,ar_break_stop_opt, ...
            default_inh_err,plot_var_test,init_geogr_head,incl_decl, ...
            star_night,dalf_close,n_days_close,opt_mod, ...
            dep_lat_degs,lat_dist,lon_dist,max_n,n_fl_seq,stop_dur, ...
            std_n_stop,n_hs,Va_mps,dep_date,std_day,n_inds);     
    
    [alf_init,fval,exitflag,output] = fzero(fun,alpha_0);
     
    disp(['optimal heading vs. South = ' num2str(alf_init)])
     % now plot trajs
    plot_ind_tracks_opt =  true;
    plot_var_disp = 3; % colour-code headings along trajectories headings 
    % colour 'perfect' tracks
    clr_perf_trx = 'c';
    plot_titl = true;
    plot_colbr = true;
    plot_cb_end = false;
    % see assess_spat_mod_runs for other options
    run_no_topo_any
    
end
