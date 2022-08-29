% clear 

if ~exist('n_inds')

    n_inds = 10000;

end
   
zero_vec = zeros(n_inds,1);
one_vec = ones(n_inds,1);
false_vec = false(n_inds,1);
true_vec = true(n_inds,1);
                    
    
 % full run loops through (many) inherited headings and
% (detection, transfer or inherited) errors in heading
if ~exist('full_run')

    full_run = true;
    disp('running full run')

end
    
if full_run
    clear alf_ps lat_ps lon_ps alfs_succ_arr p_close ...
        arr_succ arr_ts d_close dec_close  is_nan is_viable lq_discr ...
        magcl_parl_param magcl_tran_param med* mn* one_vec_perf ...
        plot_alf_idx_opt std_discr std_err std_err_all std_ll_err ...
        true_vec_perf uq_discr 

end
  
% geogr map indicates displacements according to
% 1 infininite projected plane (e.g., Mouritsen & Mouritsen 2000)
% 2 local projected plane (summed hourly)
% 3 inverse Haversine equation (most exact)
% In practivce 2 is virtually identical and faster

    geogr_map = 2; % 3; %  1; %   
    
% species goal option finds optimal routes between specified
% departure and arrival locations for given species attributes
    if ~exist('species_goal_opt')
        
        species_goal_opt = false;
        
    end
    
    if ~exist('non_stop')
        
        non_stop = false;
        
    end
    
    if ~exist('goal_rad')
        
        goal_rad = 1000;
        
    end
    
       
% err_source is 1) detection and setting of heading (alpha),
% 2) endogenous heading , 3) init Lat (2 & 3 both from natal dispersal 
% or mismatches between inherited headings and geophysical cues
    if ~full_run
        err_sources = 1; %    1 % 2 %   
    end
    
    % option to transfer heading from mag/star to sun comp at daen 
    % so avoiding the change in sun azimuth
    if ~exist('transfer_to_sun_dawn')
        
        transfer_to_sun_dawn = false;
        
    end
    
% delete(gcp('nocreate'))
% parpool('local',4)

% set orientation program(s) to test (if not done):
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
    if ~exist('or_progs')

        or_progs = 1;

    end
    
% default options for sun compass (resetting, polar strategy),
% and also for programs involving sun compass calibration, specify
% option to inherit sun compass (true) or geogr heading (false)

    % default init migr direction is geographic
    if ~exist('init_geogr_head')
        init_geogr_head = true;
    end

    % same for star compass at night
    if ~exist('star_night')
        star_night = true;
    end 

    if ~exist('tc_clock_az')
          tc_clock_az = 'local'; % 'reference'; % 
    end
    
    if ~exist('sun_inher_opt')
        
        sun_inher_opt =  false; % true %  
        
    end
        
    % options for scaling cal error in polar day and night
    % (1) that sun compass is replaced by magn compass (same precision)
    % until the bird arrives in sunset zone
    polar_magn_option = false; %  true %
    reset_clock_after_arctic =  false; % true %
    % or by the min. daily sun elevation angle above or below the horizon, resp.
    polar_err_scale_opt = false; %  ~polar_magn_option; %  true % 
    max_det_err_no_sunset = 60*pi/180;
    max_det_err_no_sun = 90*pi/180;
    % for TCS with sun reset (11, 12) we have an option, tc_reset_only, 
    % where only the time lag is refreshed,
    % i.e., where the sun heading is not held fixed during stopover
    tc_reset_only = false; % true; %  
    % inher option magcl 1) projected heading 2) geogr haeading to imprint
    % proj head
    
    
    if ~exist('is_trans_comp') % need to update 
        
        is_trans_comp = false;
        
    end
    
    if ~exist('magcl_inh_opt')
        
        magcl_inh_opt = 2; % 1; %   
        
    end
    
    % alf zug frac gives amount change for strats 1 & 2
    alf_frac_zug = 0.5; % NaN %
    % zug strat 1 part way to S (counterclockw)
    % 2 part way to N (clockw)
    % 3 head due S
    % 4 is after 10 degs Lat change switch to ccl Zugkn as in 1
    zug_strat = 1; % 2 % 3 %  
     % initial or final sprint to avoid large changes in sun azimuth esp. 
    % close to the poles
    init_sprint =   false; % true; %   ~polar_magn_option  %
    final_sprint = false; % true % 
    final_spr_lat_thr = 20; % # degs final sprint   

% set default parameters if not already defined
    
    % auto regressive error term for det error
    if ~exist('ar_error_det')
        ar_error_det = 0;
        % option to break up autocorrelation after stopover
        % eg weather or perhaps detour related error will stop
        ar_break_stop_opt = false;
    end
 
    % set stepwise error for transferring heading to in flight compass 
    % (if applicable) 
    if ~exist('err_trans')
        
        err_trans = zeros(size(err_dets));
        
    end
    
    % set error for maintaining heading in flight during each step
    if ~exist('err_mnts')
        
        err_mnts = zero_vec;
        
    end
    
    if ~exist('ar_error_mnt')
        ar_error_mnt = 0;
%         default_mnt_err =  0; % 
    end
 
    % ditto for drift error 
    if ~exist('err_drfts')
        
        err_drfts = zero_vec;
        
    end
    
    if ~exist('ar_error_drft')
         ar_error_drft = 0;
%         default_mnt_err =  0; % 
    end
    
    if ~exist('ar_err_hr_drft')
         ar_err_hr_drft = 0;
%         default_mnt_err =  0; % 
    end   
    
    if ~exist('std_day')
        std_day = 5;
    end

    if ~exist('dep_lat_degs')
        dep_lat_degs = 65;
    end
    
    if ~exist('close_plots')
        close_plots = false;
    end

    if full_run
        % start simulations at longitude = 180 deg 
        % This works for both East -ve and West +ve initial headings
        % (vs South) if Long distance less than 180 degs
            dep_lon = pi;
            
    elseif species_goal_opt
        
         dep_lon = dep_lon_sp_degs*pi/180; 
         
    else
        
        dep_lon = pi;
        
    end
    
% convert to km/h and radians
    Va = Va_mps*3.6; % 0; % 
    dep_lat = dep_lat_degs*pi/180;
    dep_lon_degs = dep_lon*180/pi;
    % DelLatZug gives Lat change (rads) when to switch
    % (option currently not used)
    DelLatZug = 0.5*dep_lat; % dep_lat - 20*pi/180 % 
    
% Initialize dates

    if ~exist('dep_date')

        dep_date = [15 9 2000];
        std_day = 2; %  5 % 0

    end

    yr_0 = dep_date(3);
    mnth_0 =  dep_date(2);  %  8 %  7 % 10 %    11 %  7  % 
    day_0 = dep_date(1); %   30 % 13  %          %    
    
% plotting parameters
    or_pr_str = {'Loxodrome','Linear','Magclinic transverse component', ...
        'Magclinic parallel component',...
    'vert slope 1','fixed Sun compass','Time comp Sun compass','Magcl trans Hz', ...
       'mixed sun magn','mixed TCsun mag','TC Sun reset','mixed TC reset mag', ...
       'great circ','fixed reset'};
    plot_clrs = {'b','r','g','m','k','c'};
    plot_symbs = {'o','x','h','>','d','+'};

 % some constants for numeric efficiency 
    pi_19_18 = pi*19/18;
    % and 10 degs beyond lat dest\
    pi_36 = pi/36;
    % and 90 deg thresholds
    pi_2 = pi/2;
    two_pi = 2*pi;
    radDeg = 180/pi;
    degRad = pi/180;

%     try
        
% loop through the orientation programs 
for iprog = 1:numel(or_progs)
    
    or_prog = or_progs(iprog);
    is_sun_comp = ismember(or_prog, [6 7 9 10 11 12 14]); 
    tc_sun_comp = ismember(or_prog, [7 10 11 12]); 
    comb_comp = ismember(or_prog, [9 10 12]); 
    %
    magcl = or_prog == 3 || or_prog == 4 || or_prog == 8; % true % false %;   
    magcl_nr = (or_prog == 3) + 2*(or_prog == 4) + 3*(or_prog == 8); % 1 %  
    
    % stop_sun_rerest = 1 means flight dirn is held
    % during stopover by adjusting the heading
    % possibly including shifts from previous time compensation 
    % i.e. for TC sun comp
    stop_sun_reset = ismember(or_prog, [11 12 14]); %true % false; % 
    
    if ~full_run % choose species and or prog heading
        
        % option to use magnetic compass as primary comp (incl_decl) 
        % and fly at night (~star_night)
        geogr_heads = ~incl_decl; % is_sun_comp || 

    % alpha is heading relative to geographic South
        alfs_deg = alf_init;

    else % full run 
        
        % don't use magnetic declination for 'generic' runs
        
        geogr_heads = true;
        init_geogr_head = true;
        star_night = true;
        incl_decl = false;

        % alpha is heading relative to geographic South
        %  depends on whether sun or mag compass (latter <= 90 degs)
        if is_sun_comp || or_prog == 13
            
            alfs_deg = alfs_deg_full_sun;
            
        else
            
            alfs_deg = alfs_deg_full_mag;
             
        end
        
    end
    
    % ensure great cicrle is nonstop
    non_stop = non_stop || or_prog == 13; % true %   tc_sun_comp  % 
    
%     if ~init_sprint
%         n_fls_step = non_stop + 5*~non_stop  % 3*~non_stop % 
%         n_flts_init = n_fls_step; % 8; %  
%         n_stop_step = ~non_stop*12 %8 %  
%     else

    if non_stop
        
        n_fls_step = Inf;
        
    else
        
        n_fls_step = n_fl_seq*~non_stop; % 3 % 8 %
        
    end
    
    n_flts_init = init_sprint*10; % 5 % 8 n_fls_step; %  
    n_stop_step = ~non_stop*stop_dur; % 12 % Mon % 4 %
        
%     end

    if ~exist('std_n_stop')
        std_n_stop = 2*~full_run + 0*full_run; %  4 %   0; % 
    end
    min_n_stop = max(n_stop_step-6,2)*~non_stop; % 
    max_n_stop = (n_stop_step+6)*~non_stop; % min(n_stop_step+6,14)*~non_stop;
    

    if geogr_map == 1 && or_prog ~= 1 

        disp('error: inf plane only implemented for Lox headings')
        return

    end

    if ~exist('std_err_inits')
        
        std_err_inits = 0;
        
    end
        
    if ~exist('d_alf_init_degs')
        
        d_alf_init_degs = 0;
        
    end
        
   if ~exist('err_dets')
       
        err_dets = 0;
        
   end
    
    
    n_init_errs = numel(std_err_inits);
    n_alf_inits = numel(d_alf_init_degs);
    n_det_errs = numel(err_dets);
    n_tran_errs = numel(err_trans);
    alfs = alfs_deg'*degRad; % (0:20:80)'*degRad;
    n_alf = numel(alfs);

    % for storing 'perfect' orientation tracks
    zero_vec_perf = zeros(n_alf,1);
    one_vec_perf = ones(n_alf,1);
    true_vec_perf = true(n_alf,1);

    plot_alf_idx_opt = 1:numel(alfs_deg); % [1 2] %[1:3]; % [1 4 8] % [1:5] %  [1:4]; % 1 % 

    if ~exist('default_det_err')
        default_det_err =  0; %   5; % 10; %
    end
    if ~exist('default_inh_err')
        default_inh_err = 0; %  5; %  
    end
    
    
    % for fixed sun compass

%     range_fl_hrs = 2; % 0; %  
%     offset_dep_dusk = 0; % 1;

% scale and threshold distance parameters
    R_Earth_km = 6371;
    % threshold for "close" to goal longitude
    R_thrs = [100 250 500 1000]; %  1500];
    % index of dist threshold used to gauge arrival success
    succ_thr = 4;
    n_thrs = numel(R_thrs);
    arr_lat_degs = dep_lat_degs - lat_dist; % 15; %  = 
    arr_lat = arr_lat_degs*pi/180; %    
    % new Dec 2020 thresh min Lat for stopping (end sprint)
    thr_spr_lat = final_sprint*(arr_lat_degs-final_spr_lat_thr)*pi/180;
    if geogr_map == 1

        arr_lat = -R_Earth_km*(dep_lat-arr_lat); % pi/6; % pi/12; % *(magcl && magcl_nr == 2); % 0; %  0 %
        arr_lat_thr = max(R_thrs(succ_thr),goal_rad);

    else

        arr_lat_thr = max(R_thrs(succ_thr),goal_rad)/R_Earth_km; % arr_lat; % 

    end

    if ~isempty(plot_alf_idx_opt)
        
        % num trajs to plot
        n_plot = min(n_inds,150); % 280
        n_tracks = min(n_inds,30); % 280
        % transparency of plotted trajectories
        trnsp_plot = 0.5;
        sz_mrkr = 50; % 80
        sz_mrkr_pfct = 120; % 180

        % y ticks and labels for Mercator projn of latitudes
        ytik_lbls = [0 20 45 65 80 85]; %[0 30 60 66.5 80]; % (-10:10:dep_lat_degs);
        ytiks = log(abs(tan(ytik_lbls*pi/180) + sec(ytik_lbls*pi/180)))*180/pi;
        min_y_pl = log(abs(tan(arr_lat-15*pi/180) + sec(arr_lat-15*pi/180)))*180/pi; %arr_lat-15
        max_rads_pl = min(2.5+dep_lat_degs,85.5)*pi/180; % 81.5*pi/180 % min(5+dep_lat_degs,89)*pi/180;
        max_y_pl = log(abs(tan(max_rads_pl) + sec(max_rads_pl)))*180/pi;
        min_x_pl = -100;
        max_x_pl = 20;

        
        n_del_dys = 10; % for colorbar scale
        
        % 1 = fly nr, 2 = delay in departure vs. 'adaptive', 3 = heading,
        % 4 == sun az
%         plot_var_sel =  5 % 3 % is_sun_comp*3 + ~is_sun_comp*6 % 2 %   1 2 % 3 % 
        plot_RL_GC_opt = true; %false %  
        if ~exist('max_n')
            max_fl_n = 40; % max_n
            max_n = 40; %
            max_t = 40; %
%         else
%             max_fl_n = max_n;
%             max_t = max_n;
        end
        max_d_cl_plot = 1500;
        d_d_cl_plot = 250;
        
        lab_sz = 10;
        titl_sz = 10;
        Ft_sz = 9;
        
    end
    
    if ~exist('max_t')
        
        max_t = 100;
        
    end

    % start at Intl dateline % Greenwich Meridian
    % dep_lon = 0 * (geogr_map~=1)*pi; 

    % calcuLate hourly change Lat if Southbound
    % Lon will have 1/sin(theta) factor (below)


    % circm_Earth = 2*pi*R_Earth_km;

    if geogr_map ~= 1 % distances measured in Radians

        hourly_del_Lat = Va/R_Earth_km;  

    else % for inf. plane approx use Cartesian distances (km)

        hourly_del_Lat = Va;  

    end

    nightly_del_Lat = n_hs*hourly_del_Lat;
%     step_del_Lat = nightly_del_Lat/(1+n_mnts);

    % initial headings are westerly (+ve); define least westerly heading
    min_head =  0; % -pi_2; % 
    % for magcl, optional to set max head (W) 
    max_head = pi_2;

    plot_syms = {'o','p','s','h','d'};
    plot_cols = {'b','m','g','r','c'};
    plot_lines = {':','--','-','-.','-'};
    bar_style = {'dotted', 'dashed','solid', 'dashdot', 'solid'};
    LWe = 1.5; 

    % initialize dates for all runs
    curr_date_0.year = yr_0*one_vec;     
    curr_date_0.month = mnth_0*one_vec;
    if ~full_run && plot_opt && plot_var_test == 5
        % vary date 
        dep_days = day_0 + mod(cumsum(one_vec),ndates+1) -1;  % 
    else
        dep_days = day_0*one_vec + round(randn(size(one_vec))*std_day); %
    end
    curr_date_0.day = dep_days;
    curr_date_0.min = zero_vec;
    curr_date_0.sec = zero_vec;   
    dates_0 = datenum([curr_date_0.year curr_date_0.month floor(curr_date_0.day)]);
    med_date_0 = median(dates_0);  
    doys_0 = day(datetime(datevec(dates_0)),'dayofyear' );
    
    
    %% start loop through std error values 

    for err_source = err_sources

        if full_run
             disp(['err source ' num2str(err_source)])
             tStart = tic;
        end


        if close_plots ==  true %false %   
            close all
            close_plots = false; %  
        end

        n_errs = (err_source == 1)*n_det_errs + (err_source == 2)*n_alf_inits + ...
        (err_source == 3)*n_init_errs;


        for i_deplat = 1:numel(dep_lat_degs) % n_hs)

%             n_hs_i = n_hs(1);
            dep_lat = dep_lat_degs(i_deplat)*degRad;
            leg_str{i_deplat} = ['departure Lat = ' num2str(dep_lat_degs(i_deplat)) '^o'];


            % first simulate 'perfect' orientation for each heading

                clear curr_date

                for i_err = 1:n_errs
                    
                   % Jan 2021
                   % adjust det error for (any) autocorrealtive error (for equivalence with
                   % non auto correlative case, i.e. so that the overall variance is about
                   % the same
                   
                   if err_source == 1
                           err_det_base = err_dets(i_err);
                           err_mnt_base = err_mnts(i_err);
                           err_drft_base = err_drfts(i_err);
                   else
                           err_det_base = default_det_err;
                           err_mnt_base = err_mnts(i_err);
                           err_drft_base = err_drfts(i_err);
                   end
%                    err_det_base = (err_source == 1) + (err_source ~= 1)*default_det_err;
                    err_det = sqrt((1-ar_error_det^2)*err_det_base^2);
                    err_mnt = sqrt((1-ar_error_mnt^2)*err_mnt_base^2);
                    err_drft = sqrt((1-ar_error_drft^2)*err_drft_base^2);
                    err_hr_drft = sqrt((1-ar_err_hr_drft^2)*err_drft_base^2);                    
                    
                    % assume no AR (detect and drift) error after stopover
                    if ar_break_stop_opt
                         err_det_stop = err_det_base;
                         err_drft_stop = err_drft_base;
                    else
                        err_det_stop = err_det;
                        err_drft_stop = err_drft;
                    end
                    
                    if full_run
                        display(['error = ' num2str(i_err)])
                        tic
                    end
                    
                    if err_source == 1

        %                 err_vals = err_dets;
                        % von mises parameter ~ 1/var
                        kappa_err_det = 1/(err_det*degRad)^2;
                        kappa_err_base = 1/(err_det_base*degRad)^2;
                        kappa_err_trans = 1/(err_trans(i_err)*degRad)^2;
                        kappa_err_mnt = 1/(err_mnt*degRad)^2;
                        kappa_err_drft = 1/(err_drft*degRad)^2;
                        kappa_err_hr_drft = 1/(err_hr_drft*degRad)^2;
                        kappa_err_stop = 1/(err_det_stop*degRad)^2;
                        kappa_err_stop_drft = 1/(err_drft_stop*degRad)^2;
                        kappa_err_alf = 1/(default_inh_err*degRad)^2;
                        err_init_deg = 0; % zeros(n_inds,1);
                        
                    elseif err_source == 2 

        %                 err_vals = alf_inits;
                        kappa_err_det = 1/(default_det_err*degRad)^2;
                        kappa_err_trans = 1/(err_trans*degRad)^2;
                        kappa_err_mnt = 1/(err_mnt*degRad)^2;
                        kappa_err_drft = 1/(err_drft*degRad)^2;
                        kappa_err_hr_drft = 1/(err_hr_drft*degRad)^2;
                        kappa_err_stop = kappa_err_det;
                        kappa_err_stop_drft = kappa_err_drft;
                        kappa_err_base = kappa_err_det;
                        kappa_err_alf = 1/(d_alf_init_degs(i_err)*degRad)^2;
                        err_init_deg = 0; % zeros(n_inds,1);

                    else

        %                 err_vals = std_err_inits;
                        kappa_err_det = 1/(default_det_err*degRad)^2;
                        kappa_err_trans = 1/(err_trans*degRad)^2;
                        kappa_err_mnt = 1/(err_mnt*degRad)^2;
                        kappa_err_drft = 1/(err_drft*degRad)^2;
                        kappa_err_hr_drft = 1/(err_hr_drft*degRad)^2;
                        kappa_err_stop = kappa_err_det;
                        kappa_err_stop_drft = kappa_err_drft;
                        kappa_err_base = kappa_err_det;
                        kappa_err_alf = Inf;
                        err_init_deg = std_err_inits(i_err); %

                    end

                    % set initial latitude error
                    % (e.g. from natal dispersion for given 'adapted' heading)
                    err_init1 = err_init_deg*randn([n_inds,1])*degRad*sqrt(2)/2;
                    err_init2 = err_init_deg*randn([n_inds,1])*degRad*sqrt(2)/2;
                    % set initial locations
                    if geogr_map ~= 1

                        lon_es = dep_lon*one_vec + err_init1; %  - (ia-1)*degRad; % zero_vec;
                        % theta same but with errs
                        lat_es = dep_lat*ones(n_inds,1) + err_init2; % *randn(n_inds,1);

                    else % infinite plane approximation (e.g., Mouritsen & Mouritsen 2000)

                         lon_es = zero_vec; %  - (ia-1)*degRad; % zero_vec;
                         % theta same but with errs
                         lat_es = R_Earth_km*err_init.*ones(n_inds,1) + err_init; % *randn(n_inds,1);                                              

                    end
  
                    % set locn-specific initial geomg field if incl declination
                   if ~geogr_heads || ~star_night % && ~init_geogr_head

                        % use median date (won't make much diff)
%                         tic
                        [Bx, By, Bz] = igrf(med_date_0, ...
                            lat_es*180/pi, lon_es*180/pi, 0);
                        decln = atan2(By,Bx);
                        decln_0 = decln;
                        hz_inten_0 = hypot(Bx,By);
                        incln = min(atan(Bz./hz_inten_0),pi/2);
                        incln_0 = incln;
%                         toc
                        
                        
%                         tic
%                         [~,~,D, I, F] = igrfmagm(zero_vec', ...
%                             lat_es'*180/pi, lon_es'*180/pi, 2018*one_vec');
% %                         decln = atan2(By,Bx);
% %                         decln_0 = decln;
% %                         hz_inten_0 = hypot(Bx,By);
% %                         incln = min(atan(Bz./hz_inten_0),pi/2);
% %                         incln_0 = incln;
%                         toc

                   end

                   % set magcl paramter
                   if geogr_heads
                        magcl_tran_param = 2*tan(dep_lat)./sin(alfs);
                        magcl_parl_param = 2*tan(dep_lat)./cos(alfs);
                        magcl_hz_tran_param = cos(dep_lat).*sin(alfs);
                   else
                        magcl_tran_param = tan(incln_0)./sin(alfs);
                        magcl_parl_param = tan(incln_0)./cos(alfs);
                        magcl_hz_tran_param = cos(dep_lat).*sin(alfs);               
                   end
           
                    if full_run
                       disp(['std err degs ' num2str(err_init_deg)])
                    end

%                     tic

                    % now simulate each heading with replicates to assess error
                    
                    mn_discr =  NaN*ones(n_alf,1);
        %                 geomn_discr(i_err,ia) = geomean(abs(discr_pe)); % (:,ia)
                    std_discr =  NaN*ones(n_alf,1);             
                    mean_dur_es = NaN*ones(n_alf,1);
                    arr_succ = zeros(n_alf,1);
                    lq_discr =  NaN*ones(n_alf,1);
                    md_discr =  NaN*ones(n_alf,1);
                    uq_discr =  NaN*ones(n_alf,1);
                    std_ll_err =  NaN*ones(n_alf,1);
                    mn_arr_lon_ps = NaN*ones(n_alf,1);
                    alfs_succ_arr = NaN*ones(n_alf,1);

                    stops_over = false_vec;
                    
                    if err_source == 1       

                        alf_ps(1:n_alf,1:n_inds,1)= repmat(alfs,[1 n_inds]);
                        lat_ps(1:n_alf,1:n_inds,1) = dep_lat_degs*pi/180;

                        alf_ps(1:n_alf,1:n_inds,2:max_n)= NaN;
                        lat_ps(1:n_alf,1:n_inds,2:max_n) = NaN;
                        lon_ps(1:n_alf,1:n_inds,2:max_n) = NaN;
                        lon_ps(1:n_alf,1:n_inds,1) = dep_lon;   

                    end
                    
                    % set error for cases of 24 hour light or dark 
                    if is_sun_comp

                        if ~comb_comp

                            max_det_err_light = max_det_err_no_sunset;
                            max_det_err_dark = max_det_err_no_sun;

                        else

                            max_det_err_light = NaN;
                            max_det_err_dark = NaN;

                        end


                    end


                    for ia = 1:n_alf % _arr par   
                        
                        % (re-)initialize date fields (for sun compass &
                        % magn compass if incl dec;lination)
                        curr_date = curr_date_0;
                        dates = dates_0;
                        med_date = med_date_0;
                        doys = doys_0;

                        if isinf(kappa_err_alf)
                            
                            alf0_s = alfs(ia)*one_vec; 
                            magn_params = magcl_tran_param(ia)*one_vec;

                        else

%                             if comb_comp % combine two compasses 
%                                 % assume here only one used on 1st depart
% 
%                                 alf0_s = alfs(ia) + vmrand(0, kappa_err_alf, [n_inds 1]);
%                                 magn_params = magcl_tran_param(ia)*one_vec;                                                                    
% 
%                             else

                            if or_prog ~= 3 && or_prog ~= 4 && or_prog ~= 8

                                alf0_s = alfs(ia) + vmrand(0, kappa_err_alf, [n_inds 1]); 
                                magn_params = magcl_tran_param(ia)*one_vec;

                            else

                                if or_prog == 3

                                    if magcl_inh_opt == 1 % vary inh projn

                                        magn_arg = min(atan(magcl_tran_param(ia))+ vmrand(0, kappa_err_alf, [n_inds 1]),pi_2);
                                        magn_params = tan(magn_arg); % magcl_tran_param(ia) ;
                                        alf_arg = 2*tan(lat_es(:,1))./magn_params;
                                        ok_magcl = abs(alf_arg) <= 1;
                                        alf0_s(ok_magcl,1) = asin(alf_arg(ok_magcl));
                                        alf0_s(~ok_magcl,1) = pi_2*sign(alfs(ia))*ones(sum(~ok_magcl),1);

                                    else % vary initial heading to imprint projn

                                        alf0_s = alfs(ia) + vmrand(0, kappa_err_alf, [n_inds 1]);
                                        
                                        if ~incl_decl % dipole model
                                            magn_params = 2*tan(dep_lat)./sin(alf0_s);
                                        else
                                             magn_params = tan(incln_0)./sin(alf0_s);                                           
                                        end
    %                                     magn_arg = min(atan(magcl_tran_param(ia)),pi/2);
    %                                     magn_params = tan(magn_arg)*one_vec; % magcl_tran_param(ia) ;
    %                                     alf_arg = 2*tan(lat_es(~done,i_t))./magn_params(~done) ....
    %                                         + vmrand(0, kappa_err_alf, [n_inds 1]);                                    

                                    end             

                                elseif or_prog == 4

                                    magn_arg = min(atan(magcl_tran_param(ia))+ vmrand(0, kappa_err_alf, [n_inds 1]),pi_2);
                                    if ~incl_decl % dipole model
                                        magn_params = tan(magn_arg); % magcl_tran_param(ia) ;
                                        alf_arg = 2*tan(lat_es(~done,1))./magn_params(~done);
                                        ok_magcl = abs(alf_arg) <= 1;
                                        alf0_s(ok_magcl) = acos(alf_arg(ok_magcl));
                                        alf0_s(~ok_magcl) = 0;
                                    else
                                        % to do: use incln and cosine
                                        keyboard
                                        % magn_params = tan(incln_0)./cos(alf0_s);                                           
                                    end
                                    
    %                                 magn_params = magcl_tran_param(ia) + vmrand(0, kappa_err_alf, [n_inds 1]);

                                end

                            end

                        end
                        

                        %% initialize simulations
                        i_t = 1;
                        done = false_vec;
                        n_not_done = n_inds;

                        % initialize closest approach and index
        %                 i_close = one_vec;
        %                 d_close = Inf*one_vec;

                        % set migr duration to infinite (i.e. failed)
                        mig_durs = Inf*one_vec;

                        % flag for already passed half (and switched headings)
                         % not currently used
                        has_passed =  false_vec;                
                         
                         % for sun compass inherited option, adjust
                          % initial headings
                           if sun_inher_opt  || plot_var_test >= 1 % use sun azimuth to determine initial headings

                              % first compute sun az for mean dep date                              
                               date_0_mn = datenum([yr_0 mnth_0 day_0]);
                               doy_0_mn = day(datetime(datevec(date_0_mn)),'dayofyear' );
                               
                           end 
 
                       % initialize det error (use kappa_err_base since not
                       % autocorrelated)
                        if ~isinf(kappa_err_det)
                            if ~comb_comp
                                
                                err_dtcs = vmrand(0, kappa_err_base, [n_inds 1]); % err_det*randn(n_not_done,1);
                                
                                % set for first hour based on err_s (i.e. CALIBR error) 
                                % and stepwise auto-corr

                            else % combined magn & (tc or fixed) sun comp = 2 random errors
                                err_dtcs = (vmrand(0, kappa_err_base, [n_inds 1]) + ...
                                    vmrand(0, kappa_err_base, [n_inds 1]))/2;
                            end
                            
                            
                        else
                            err_dtcs = zeros(n_inds,1);
                             
                        end
                        
                        if is_trans_comp && ~isinf(kappa_err_trans)
                            
                            err_trs = vmrand(0, kappa_err_trans, [n_inds 1]); % err_det*randn(n_not_done,1);
                                                       
                        else
                            
                            err_trs = zeros(n_inds,1);
                            
                        end
                        
                        
                       % initialize drift error (use kappa_err_drft_stop since not
                       % autocorrelated)
                        if ~isinf(kappa_err_drft)
                               
                            err_dfs = vmrand(0,  kappa_err_stop_drft, [n_inds 1]); % err_det*randn(n_not_done,1);
                                                       
                        else
                            err_dfs = zeros(n_inds,1);
                             
                        end  
                        
                                                
%                         if ~isinf(kappa_err_mnt)
%                                 
%                           err_ms = vmrand(0, kappa_err_mnt, [n_inds 1]); 
% 
%                         else
                            
                        % maintenance error in first sub-step (e.g., flight hour) is zero (included
                        % in cal error)
                        err_ms = zeros(n_inds,1);

%                         end
                        
                        % no previous error for transfer to tc sun comp
                        % option
                        err_tc_pr_dec = zero_vec;

%                        dates_0 = dates;

                        if is_sun_comp % fixed or TC sun compass progs possibly combined

                            % sun_offset is due to changing sun azimuth
                            sun_offset = zero_vec;
                            tc_offset = zero_vec;
                            tc_ref_lon = dep_lon*one_vec;
                            tc_ref_lat = dep_lat*one_vec;

%                           [sun_az, err_init] = detc_sun_az(lat_es,doys, ...
%                               err_init,max_det_err_light,max_det_err_dark);

    % sun az can augment initial det error in polar night or day
                            if polar_err_scale_opt % assume use sun compass despite being 
                                % in polar summer or winter

                                [sun_az, err_dtcs, is_sunset]  = calc_sun_az(lat_es(:,i_t),doys, ...
                                  err_dtcs,max_det_err_light,max_det_err_dark); 

%                                 is_sunset = true_vec;
                                
                            else % if polar_magn_option, no_sunset option might indicate switch 
                                 % to magnetic compass (or perhaps polarized light if that is possible but error prone for
                                 % sun farter from horizon)
%                                  try
                                  [sun_az, ~, is_sunset]  = calc_sun_az(lat_es(:,i_t),doys); 
                                  if polar_magn_option && comb_comp && ~isinf(kappa_err_base) && sum(~is_sunset) > 0 
                                      % undo 'halving of error
                                      err_dtcs(~is_sunset) = vmrand(0, kappa_err_det, [sum(~is_sunset) 1]);
                                  end
%                                  catch
%                                      keyboard
%                                  end

                            end

                                sun_az_0_mn = calc_sun_az(dep_lat,doy_0_mn);

                               alf0_s = alf0_s +  sun_inher_opt*(sun_az - sun_az_0_mn);

       %                        toc
                               sun_az_0 = sun_az;
                               sun_az_pr = sun_az_0;

                               % ref azimuth for changing heading to sun az 

                               sun_ref_az = sun_az;
                               
                               is_curr_sun_comp = is_sunset | ~polar_magn_option;
                                is_curr_magn_comp = ~is_sunset & polar_magn_option;

                        else                

                          
                           % set no sunset option to false 
                            is_sunset = true_vec;
                            is_curr_sun_comp = false_vec;
                            is_curr_magn_comp = true_vec;
                            
                        end
                                                
                        % offset by magn_decl if magn compass
                        if ~init_geogr_head % if geomagn heads,
                            % convert to geogr heads
                            alf0_s = mod(pi + ...
                                alf0_s + decln_0,2*pi)-pi;
                        end
                        
                        if ~full_run % plot_opt 
                            
%                             if plot_var_disp == 3
                                
                                all_alphs = 180 + (alf0_s+err_dtcs)*180/pi-180;
                                
%                             elseif plot_var_disp == 4
                              if is_sun_comp   
                                 all_sun_azs = 180 + sun_az_0*180/pi;
                              end
                              
%                             elseif plot_var_disp == 5 % doy
                                
                                all_doys = doys;
                                
%                             elseif plot_var_disp == 6 % residual heading                              
                                 
                                all_errs = err_dtcs*180/pi;
                                
%                              end
                            
                        end

                        % initialize counter for plotting 'perfect' orientation 
                        % with short option
                        days_perf = 0;
                        
                        % first flight is error and update free
                        dalf = zero_vec;

                        
                        %% loop stepwise until migration or time are done
                        while n_not_done > 0 && i_t < max_n

%                             idx_not_done = find(~done);

                            % update location given current heading for
                            % unfinished migrants
                            update_locn_spat_mod
                                             
                            % now increment date and update orientn relevant parameters
                            i_t = i_t + 1;
                            
                            % diagnose stopovers and (next flight step) dates
                             if non_stop

                                    stop_durs_it = zero_vec; 

                            else

                                % first night is still "day zero"
                               stop_durs_it = round(min(max(n_stop_step-1 + randn(size(one_vec))*std_n_stop,min_n_stop),max_n_stop));

                             end                                
                             % add is sunset for polar_magn_option
%                              try
                            stops_over(~done,1) = i_t >=  n_flts_init  &  mod((i_t-1),n_fls_step)==0 & ...
                                lat_es(~done,i_t) >= thr_spr_lat & (~polar_magn_option | is_sunset(~done));
                            curr_date.day(~done) = curr_date.day(~done) + ...
                            ((mod(i_t,1)==0) + stops_over(~done).*stop_durs_it(~done)); %  ...
                             
                            n_stops = sum(~done & stops_over);
                            n_not_stops = sum(~done & ~stops_over);
                           
%                              catch
%                                  keyboard
%                              end
%                             curr_date.hour(~done) = civ_dusk.*one_vec(~done);
                            dates(:,i_t) = datenum([curr_date.year curr_date.month curr_date.day]); %  curr_date.hour ...
%                             curr_date.min curr_date.sec]);

                            doys(~done) = day(datetime(datevec(dates(~done,i_t))),'dayofyear' );

                            % 'perfect' migration follows exact stopover
                            % schedule
                            days_perf(i_t) = days_perf(i_t-1) + 1 + ...
                                (mod((i_t-1),n_fls_step)==0).*(n_stop_step-1);
                            
                            % determine stepwise error in compass headings
                            % Will depend on stopover or not if error is
                            % stepwise autocorrelated
                            if ~isinf(kappa_err_base)
                                % Jan 2021 added autocorrelative term for detection
                                if ~comb_comp
                                    
                                    if ar_error_det == 0 || ~ar_break_stop_opt
                                        err_dtcs(~done) = ...
                                            vmrand(0, kappa_err_det, [n_not_done 1]); % err_det*randn(n_not_done,1);
                                    else
                                        
                                        if n_not_stops ~= 0
                                            err_dtcs(~done & ~stops_over) = mod(err_dtcs(~done & ~stops_over)*ar_error_det + ...
                                                vmrand(0, kappa_err_det, [n_not_stops 1]) +pi,2*pi) -pi; 
                                        end
                                        if n_stops ~=0
                                            
                                            err_dtcs(~done & stops_over) = vmrand(0, kappa_err_stop, [n_stops 1]); 
                                        end
                                            
                                    end
                                    
                                elseif ar_error_det == 0 || ~ar_break_stop_opt % combined magn & (tc or fixed) sun comp = 2 indep random errors
                                    
                                    err_dtcs(~done) = (vmrand(0, kappa_err_det, [n_not_done 1]) + ...
                                        vmrand(0, kappa_err_det, [n_not_done 1]))/2;
                                    
                                else
                                    
                                     if n_not_stops ~= 0                                   
                                        err_dtcs(~done & ~stops_over) = (vmrand(0, kappa_err_det, [sum(~done & ~stops_over) 1]) + ...
                                            vmrand(0, kappa_err_det, [n_not_stops 1]))/2;
                                     end
                                     if n_stops ~= 0
                                         err_dtcs(~done & stops_over) = (vmrand(0, kappa_err_det, [sum(~done & stops_over) 1]) + ...
                                            vmrand(0, kappa_err_stop, [n_stops 1]))/2;  
                                     end
                                     
                                end

                            else
                                err_dtcs(~done) = zeros(n_not_done,1);
                            end        
                            
                            if ~isinf(kappa_err_drft)
                                % Jan 2021 added autocorrelative term for detection
                                if ~comb_comp
                                    
                                    if ar_error_drft == 0 || ~ar_break_stop_opt
                                        err_dfs(~done) = ...
                                            vmrand(0, kappa_err_drft, [n_not_done 1]); % err_det*randn(n_not_done,1);
                                    else
                                        
                                        if n_not_stops ~= 0
                                            err_dfs(~done & ~stops_over) = mod(err_dfs(~done & ~stops_over)*ar_error_drft + ...
                                                vmrand(0, kappa_err_drft, [n_not_stops 1]) +pi,2*pi) -pi; 
                                        end
                                        if n_stops ~=0
                                            err_dfs(~done & stops_over) = vmrand(0, kappa_err_stop_drft, [n_stops 1]); 
                                        end
                                            
                                    end
                                    
                                elseif ar_error_drft == 0 || ~ar_break_stop_opt % combined magn & (tc or fixed) sun comp = 2 indep random errors
                                    
                                    err_dfs(~done) = (vmrand(0, kappa_err_drft, [n_not_done 1]) + ...
                                        vmrand(0, kappa_err_drft, [n_not_done 1]))/2;
                                    
                                else
                                    
                                     if n_not_stops ~= 0                                   
                                        err_dfs(~done & ~stops_over) = (vmrand(0, kappa_err_drft, [sum(~done & ~stops_over) 1]) + ...
                                            vmrand(0, kappa_err_drft, [n_not_stops 1]))/2;
                                     end
                                     if n_stops ~= 0
                                         err_dtcs(~done & stops_over) = (vmrand(0, kappa_err_drft, [sum(~done & stops_over) 1]) + ...
                                            vmrand(0, kappa_err_stop_drft, [n_stops 1]))/2;  
                                     end
                                     
                                end

                            else
                                err_dfs(~done) = zeros(n_not_done,1);
                            end   
                            
                            
                            if is_trans_comp && ~isinf(kappa_err_trans) %% tranfer error not stepwsie autocorrelated
                            
                                err_trs(~done) = vmrand(0, kappa_err_trans, [n_not_done 1]); % err_det*randn(n_not_done,1);

                            else

                                err_trs(~done) = zeros(n_not_done,1);

                            end
                            
                            
                            % reset first step's maintenance error (zero
                            % since included in detection
                              err_ms = zeros(n_inds,1); 

                            
%                             if is_sun_comp || non_stop || plot_var_test == 5 || ~star_night || incl_decl
                                                                
%                             end
                            
                            % if geomagn heads will need decln and incln for 
                            % Lox and mgclheads
                            if ~geogr_heads || ~star_night
                            try
                                % use median date (won't make much diff)
                                med_date = median(dates(~done,i_t));
                                [Bx, By, Bz] = igrf(med_date, ...
                                    lat_es(~done,i_t)*180/pi, ...
                                    lon_es(~done,i_t)*180/pi, 0);
                                decln(~done) = atan2(By,Bx);
                                hz_inten = hypot(Bx,By);
                                incln(~done,1) = min(atan(Bz./hz_inten),pi/2);
                            catch
                                keyboard
                            end
                            end
                            
                            if magcl && magcl_nr == 1
                                
                                if ~incl_decl
                                    dalf_arg = 2*tan(lat_es(:,i_t))./magn_params;
                                else
                                    dalf_arg = tan(incln)./magn_params;
                                end

                                magcl_ok = abs(dalf_arg)<=1;
    %                             clear dalf  
                                dalf(~done & magcl_ok,1) = asin(dalf_arg(~done & magcl_ok))-alf0_s(~done & magcl_ok);
                                dalf(~done & ~magcl_ok,1) = -alf0_s(~done & ~magcl_ok) + pi_2.*sign(alf0_s(~done & ~magcl_ok)); % ; %

                            elseif magcl && magcl_nr == 2

                                if ~incl_decl
                                    dalf_arg = (2*tan(lat_es(:,i_t)))./magn_params;
                                else
                                    dalf_arg = tan(incln)./magn_params;
                                end
        %                         dalf_arg = magcl_parl_param(ia)./(2*tan(lat_es(~done,i_t)));
        %                          dalf_arg = magcl_parl_param(ia)*sqrt(1+4*tan(lat_es(~done,i_t)).^2);
                                magcl_ok = abs(dalf_arg)<=1;   
    %                             clear dalf
                                dalf(~done & magcl_ok,1) = ...
                                    acos(dalf_arg(~done & magcl_ok)).*sign(alf0_s(~done & magcl_ok)) ...
                                    -alf0_s(~done & magcl_ok);
                                dalf(~done & ~magcl_ok,1) = -alf0_s(~done & ~magcl_ok); % +pi/2; %; % -pi/4; %; %        

                            elseif or_prog == 1

                                     dalf = zero_vec;

                            elseif is_sun_comp

                                if polar_err_scale_opt
                                    
                                       [sun_az(~done),err_dtcs(~done)]  = calc_sun_az(lat_es(~done,i_t),doys(~done), ...
                                  err_dtcs(~done),max_det_err_light,max_det_err_dark);
                             
                                else

                                      [sun_az(~done),~,is_sunset(~done)] = calc_sun_az(lat_es(~done,i_t),doys(~done));   
                                      
                                end

                               if i_t > 1 % && stop_sun_reset == 1 % (mod(i_t-1,n_fls_step)==0 && stop_sun_reset == 1)
%                                     if n_stops > 0      
%                                         keyboard
%                                     end

                                   if stop_sun_reset == 1 
                                       % add offset to account for keeping
                                       % preferred direction fixed over 
                                       % extended stopover
                                      sun_offset(stops_over) = ~tc_reset_only*dalf(stops_over,1) ...
                                          + tc_reset_only*(sun_az(stops_over) - sun_ref_az(stops_over)); % sun_offset(~done) + 0; % 
                                      tc_offset(stops_over) = 0;
                                      tc_ref_lon(stops_over) = lon_es(stops_over,i_t);

                                      if strcmp(tc_clock_az,'reference')
                                          % adjust departure latitude for rate
                                          % change sun azimuth (i.e., not done
                                          % locally at each departure
                                          tc_ref_lat(stops_over) = lat_es(stops_over,i_t);
                                      end

    %                                   sun_ref_az(~done) = sun_az(~done); 
                                      sun_ref_az(stops_over) = sun_az(stops_over); 
                                      
                                   end
%                                   hh = 0;

%                                     if tc_sun_comp
%                                          tc_offset(~stops_over) = (lon_es(~stops_over,i_t) -  ...
%                                            tc_ref_lon(~stops_over)).*sin(tc_ref_lat(~stops_over,i_t));  
%                                     end
%                                   %                             catch
%     %                                 keyboard
%     %                             end
%                                elseif i_t > 1 && tc_sun_comp
                                         % TCS comp, add TC shift 

                                     if tc_sun_comp && strcmp(tc_clock_az,'local')
                                         % the time is from the original
                                         % (or ref) site but the rate
                                         % change az inferred locally
                                         % I.e., bird adjusts compass to
                                         % local sun az but at the time (of
                                         % sunset) at reference 
                                         % (original or stopover) site
                                         tc_offset(~done) = (lon_es(~done,i_t) -  ...
                                           tc_ref_lon(~done)).*sin(lat_es(~done,i_t)); %
%                           
                                     elseif tc_sun_comp % 'reference' rate change azimuth of tc clock
                                         % I.e., the bird adjusts its
                                         % compass to both the time and
                                         % rate change sun az at the
                                         % departure (original or stopover) site
                                         
                                         tc_offset(~done) = (lon_es(~done,i_t) -  ...
                                           tc_ref_lon(~done)).*sin(tc_ref_lat(~done,1)); %       
%                                          tc_offset(~done) = (lon_es(~done,i_t) -  ...
%                                            lon_es(~done,i_t-1)).*sin(lat_es(~done,i_t-1)); %  
                                       
                                     end

                               end

                               if ~polar_magn_option || sum(~is_sunset) == 0

                                   % (no change to compass type)
                                     dalf(~done,1) = (~transfer_to_sun_dawn*( ...
                                         sun_az(~done) - sun_ref_az(~done)) + ...
                                         transfer_to_sun_dawn*(err_tc_pr_dec(~done)) + ...
                                         sun_offset(~done) + tc_offset(~done))/(1+comb_comp); %
                                     
%                                                             sun_offset(~done) + tc_offset(~done))/(1+comb_comp);  
                                                        
%                                   (transfer_to_sun_dawn*(dalf(~done,1) + all_errs(:,i_t)) - decln(~done)) + ...

                               else % polar magn opt and there are non-sunset (i.e. magn compass) ind's

                                   is_curr_sun_comp = ~done & is_sunset;
                                   is_curr_magn_comp = ~done & ~is_sunset;

                                    dalf(is_curr_sun_comp,1) = (sun_az(is_curr_sun_comp) - sun_ref_az(is_curr_sun_comp) + ...
                                         sun_offset(is_curr_sun_comp) + tc_offset(is_curr_sun_comp))/(1+comb_comp); 


                                    dalf(is_curr_magn_comp,1) = 0;  

                                 % set up new (first) beginning of sun compass
                                 % past polar summer / winter
                                    if sum(is_curr_sun_comp & mod(sun_ref_az,pi) == 0) > 0
                                        try   
                                         set_ref_comp = is_curr_sun_comp & mod(sun_ref_az,pi) == 0;
                                        dalf(set_ref_comp,1) = ~reset_clock_after_arctic*tc_offset(set_ref_comp)/(1+comb_comp); %  

                                        if reset_clock_after_arctic          
                                            tc_offset(set_ref_comp) = 0;
                                            tc_ref_lon(set_ref_comp) = lon_es(set_ref_comp);
                                        end

                                        % in any case set new post-arctic
                                        % azimuth as reference for offset
                                        % inhereited heading
                                         sun_ref_az(set_ref_comp) = sun_az(set_ref_comp);                                   


                                     catch
                                         keyboard
                                     end

                                    end

                               end

                            else % great circle: subtract delta(Lon)*sin(thseta) from heading
                                 % as long as not being reset every stopover 
                                if i_t > 1 

    %                                 dalf = dalf_tc_sun_nxt + ...
    %                                     (lon_es(~done,i_t) -  ...
    %                                     lon_es(~done,i_t-1)).*sin(lat_es(~done,i_t));
                                    dalf(~done,1) = dalf(~done,1) + ...
                                        (lon_es(~done,i_t) -  ...
                                        lon_es(~done,i_t-1)).*sin(lat_es(~done,i_t));

                                else

                                    dalf = zero_vec;

                                end

                            end
                            
                                                       
                            % if magn compass, offset dep geogr head by magn_decl 
                            if ~geogr_heads % uses magn comp to calibrate
                                dalf(~done,1) = mod(pi + ...
                                    dalf(~done,1) + (decln(~done) ...
                                    -decln_0(~done))/(1+comb_comp),2*pi)-pi; % 
%                                 else % calibrate to magn comp, heading shifts at night with decl
%                                       dalf(~done,1) = mod(pi + ...
%                                         dalf(~done,1) + decln(~done) ...
%                                         -decln_0(~done),2*pi)-pi;                                                                      
                            end

                            if i_err == 1 && err_init_deg == 0 
    %                             try
                                    alf_p = alf0_s(~done) + dalf(~done);
                                    alf_ps(ia,~done,i_t)= alf_p;
                                    lat_ps(ia,~done,i_t) = lat_es(~done,i_t);
                                    lon_ps(ia,~done,i_t) = lon_es(~done,i_t);


    %                             catch
    %                                 keyboard
    %                             end                           
    
                            end
                            
                            if ~full_run % plot_opt
                                
%                                 if plot_var_disp == 3
                                
                                     all_alphs(:,i_t) = ...
                                         mod((dalf + alf0_s + err_dtcs)*180/pi+180,360) -180;
                                     
%                                 elseif plot_var_disp == 4
                                 if is_sun_comp
                                     all_sun_azs(:,i_t) = 180 + sun_az*180/pi;
                                 end
%                                 elseif plot_var_disp == 5 % doys
                                    
                                    all_doys(:,i_t) = doys;
                                    
%                                 elseif plot_var_disp == 6 % doys    
                                    
                                    all_errs(:,i_t) = err_dtcs*180/pi; % abs(err_dtcs)*180/pi;
   
%                                 end
                                
                            end


                            % allow to fall below arr lat by 5 degs to catch
                            % 'closest' landing with stochasticity
                            if geogr_map ~= 1

%                                 if plot_opt
%                                     
% %                                     if plot_ind_tracks_opt
%                                        new_done = (lat_es(~done,i_t) <= arr_lat - arr_lat_thr) | ...
%                                        (lat_es(~done,i_t) <= lat_es(~done,1) - arr_lat_thr) & ...
%                                        (abs(lon_es(~done,i_t)-dep_lon) >= pi_19_18) | ...
%                                        ((lat_es(~done,i_t) <= Inf) & i_t == max_n);
%                                    
% %                                     else % plot maps 
% %                                         new_done = mod(lon_es(~done,i_t),2*pi) >= mod(lon_goal_centr,360)*pi/180 & ...
% %                                             lat_es(~done,i_t) >= lat_goal_centr*pi/180;
% %                                     end
%                                 
%                                 else
                                    
%                                        new_done = (lat_es(~done,i_t) <= arr_lat - arr_lat_thr) | ...
%                                        (lat_es(~done,i_t) <= lat_es(~done,1) - arr_lat_thr) & ...
%                                        (abs(lon_es(~done,i_t)-dep_lon) >= pi_19_18) | ...
%                                        i_t == max_n;
                                       new_done = (lat_es(~done,i_t) <= arr_lat - arr_lat_thr) | ...
                                                                              i_t == max_n;

%                                 new_done = (lat_es(~done,i_t) <= arr_lat - arr_lat_thr) | ...
%                                                           i_t == max_n;
                                                                   
%                                 end
    %                            if any(abs(lon_es(~done,i_t)-dep_lon) >= pi_19_18)
    %                                
    %                                keyboard
    %                                
    %                            end

                            else % inf plane approx

                               new_done = lat_es(~done,i_t) <= arr_lat-arr_lat_thr;

                            end

                            done(~done) = new_done;
                            
%                             if sum(~done) ==0
%                                 keyboard
%                             end

                            % store previous change in heading for T-comp sun
                            % comp (won't be empty until all done)

                            if is_sun_comp
    %                             dalf_tc_sun_nxt = dalf(~new_done);
    %                             try
                                sun_az_pr = sun_az;
                                if incl_decl 
                                  err_tc_pr_dec = err_dtcs + tc_offset + decln;
                                else
                                  err_tc_pr_dec = err_dtcs + tc_offset;  
                                end
                                
    %                             catch
    %                                 keyboard
    %                             end
                            end

                            mig_durs(done) = min(mig_durs(done),i_t);

                             n_not_done = sum(~done);

                        end

                    %% update (and optionally plot) the results for this ia, i_err combo 

                            assess_spat_mod_runs

                    end


                if full_run
                    toc
                end

                end

        end

        if full_run
             toc(tStart)
        end


    end
    
    if ~full_run && plot_spec_map_opt
        
       plot_spat_mod_conts_species % plot_spat_mod_ind_species
        
    end
    
end

%     catch
%         keyboard
%     end

