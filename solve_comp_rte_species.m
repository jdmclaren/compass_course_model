
clear 

 plot_ms_opt = true; 
 
 
 
 % option to save3 output
%  you need to set up5 some directories (I should automate)
% see bottom of this script
% if ~exist('save_opt')
  save_opt = true; % false;  %%        
% end

species_list = ["Kirtlands Warbler", "Nathusius","Monarch", ...
     "Finn Marsh Warbler","GreyCheekThrush","Siberian Willow Warbler", ...
      "France Marsh Warbler","Alk Wheatear","Sib Will Warb S Hem", ...
      "Iql Wheatear","Eleonoras Falcon","Ring Ouzel","Hoopoe","BU Rosefinch","RF Bluetail RU"]; 

disp('Choose one or more species from dialog box');
species_list_prompt = {'Kirtlands Warbler Michigan', 'Nathusius Bat Latvia', ...
    'Monarch Butterfly Quebec', 'Marsh Warbler Finland',  ...
    'Gray Cheeked Thrush Yukon','Willow Warbler East Siberia', 'Marsh Warbler France', ...
    'Wheatear Alaska','Sib Willow W to S Hemisph','Wheatear Iqaluit', ...
    'Eleonoras Falcon','Ring Ouzel','Hoopoe','BU Rosefinch','RF Bluetail RU'};
% choose species
prompt = {'Choose one or more species'};
species_nr = listdlg('PromptString',prompt,'ListString',species_list_prompt);
species_tests =  species_list(species_nr);

disp('Enter one or more compass routes from dialog box');
prompt = {'Choose one or more routes'};
or_prog_list = {'geographic loxodrome', 'geomagnetic loxodrome', ...
    'Magnetoclinic', 'Fixed sun compass', ...
    't-comp sun-comp,clock reset','t-comp sun comp, no reset', ...
    't-comp sun comp, nonstop','mix geogr lox tc-comp', ...
    'mix geomag lox tc-comp','fix sun reset'};
or_prog_idx = listdlg('PromptString',prompt,'ListString',or_prog_list); 

% choose orientation program (or_prog_idx)
% 1 = geogr loxodrome
% 2 = geomag loxodrome
% 3 = magnetoclinic compass
% 4 = fixed sun compass
% 5 = time-comp sun compass resetting clock during stopover
% 6 = time comp sun compass without clock resetting
% or_prog_idx = 5 %[ 1 2 3 4 5 6]; 

% choose map projection:
% Stereographic: great circles are straight lines)
% Mercator: loxodromes (rhumblines) are straight lines

% option for geogr (vs geomagnetic) within-flight (nightly) compass
disp('Choose star or magnetic maintenance comp in Dialogue box');
prompt = {'Star or magnetic compass'; 'during (night) flight?'};
star_list = {'Star compass','Magnetic compass'};
star_idx = listdlg('PromptString',prompt,'ListString',star_list); 
star_night =  star_idx == 1; % true; %    false; %  

% option for temporal variability in initial departure date and stopover durations
disp('Select if temporal varying schedule in Dialogue box');
prompt = {'Include (sp-appropriate)'; 'variation schedule?'};
sched_list = {'false','true','double'};
sched_idx = listdlg('PromptString',prompt,'ListString',sched_list); 
var_sched =  sched_idx -1; % true; %    false; %  

% option for geogr (vs geomagnetic) within-flight (nightly) compass

disp('Confirm if comp transferred to maintenance');
prompt = {'So is maint comp transf'; 'from detect comp for flt?'};
trans_idx = listdlg('PromptString',prompt,'ListString',{'no', 'yes'}); 
is_trans_comp =  trans_idx == 2; % true; %    false; %  

disp('Enter a map projection from dialog box');
prompt = {'Enter a map projection'};
map_pr_list = {'stereogr (gt circle = str line)','Meractor (lox = straight line)'};
map_idx = listdlg('PromptString',prompt,'ListString',map_pr_list); 

map_projs = {'stereo', 'Mercator'}; % 
map_proj = map_projs{map_idx};

 % 'Marsh Warbler' %   'Monarch'; %  'Siberian Willow Warbler'; % 'Kirtlands Warbler' % 'Monarch'; % 'GreyCheekThrush' %    
% 5 % [1 5:6] %5 %  [2 3 5] % [1 4 5] %  %  1:5; % 5; %  1:4 %  1:2; % 1 % 4 %  2; %  2 %  3 %  5 %  'Kirtlands Warbler'; %  'Marsh Warbler'; %; % 
or_prog_test = [1 1 3 6 11 7 7 12 12 14]; %  7]; %  7  12 [1 3] %        
% or_prog_idx = 1:4 %   4 % [1 2 4]; %  

% new Jan 2021 option to transfer head to sun comp at dawn so avoiding 
% the change in sun azimuth (since same at dawn and dusk
transfer_dawns = [0 0 0 0 0 0 0 0 0 0];

non_stops = [0 0 0 0 0 0 1 0 0 0];

% prompt = 'Enter std stepwise calibration error (degs)';
% input_err_str = input(prompt,'s');
% err_det_tests = str2double(input_err_str);

disp('Enter std compass error in dialog box');
prompt = {'Enter std compass detection error (degs)'};
dlgtitle = 'Stepwise error (degs)';
dims = [1 35];
definput = {'20'}; % {'0 10 20 30'}; % 
err_det_cel = inputdlg(prompt,dlgtitle,dims,definput);

err_det_tests = cell2mat(cellfun(@str2num,err_det_cel,'UniformOutput',false));

disp('Enter autoregr lags in detection error in dialog box');
prompt = {'Enter ar lags detect err to test (0-1)'};
dlgtitle = 'ar lag coeff';
dims = [1 35];
definput = {'0.'}; % {'0 0.5 0.9'}; % {'0.25'}; % 
ar_cel =inputdlg(prompt,dlgtitle,dims,definput);

ar_tests = cell2mat(cellfun(@str2num,ar_cel,'UniformOutput',false));

disp('Enter std transfer (calibrn) error in dialog box');
prompt = {'Enter std compass transfer error (degs)'};
dlgtitle = 'Stepwise error (degs)';
dims = [1 35];
definput = {'20'}; % {'0 10 20 30'}; % 
err_tran_cel = inputdlg(prompt,dlgtitle,dims,definput);

err_tran_tests = cell2mat(cellfun(@str2num,err_tran_cel,'UniformOutput',false));

disp('Enter std sub-step (hourly) maintence error');
prompt = {'Enter sub-step (hourly) maintence err (degs)'};
dlgtitle = '(single) Sub-step (maint.) err (degs)';
dims = [1 35];
definput = {'20'}; % {'0 10 20 30'}; % 
err_mnt_cel = inputdlg(prompt,dlgtitle,dims,definput);

err_mnt_tests = cell2mat(cellfun(@str2num,err_mnt_cel,'UniformOutput',false));

% star night is sun day --> name file name celestial for saving
star_night_strs = {'mag','celest'};
star_night_str = star_night_strs{1+star_night};
trans_comp_strs = {'same','transf'};
trans_comp_str = trans_comp_strs{1+is_trans_comp};

disp('Enter autoregr lags in maint. error in dialog box');
prompt = {'Enter ar lags maint. err to test (0-1)'};
dlgtitle = 'ar lag maint. coeff';
dims = [1 35];
definput = {'0.'}; %{'0.5'}; %  
ar_mnt_cel =inputdlg(prompt,dlgtitle,dims,definput);

ar_error_mnt = cell2mat(cellfun(@str2num,ar_mnt_cel,'UniformOutput',false));

disp('Enter std sub-step (hourly) (wind, topog) drift');
prompt = {'Enter std (hourly wind/topo) drift (degs)'};
dlgtitle = 'Sub-step drift err (degs)';
dims = [1 35];
definput =  {'0'}; % {'0 10 20 30'}; % {'20'}; %
err_drft_cel = inputdlg(prompt,dlgtitle,dims,definput);

err_drfts = cell2mat(cellfun(@str2num,err_drft_cel,'UniformOutput',false));

disp('Enter autoregr stepwise lags in drift error in dialog box');
prompt = {'Enter ar lags drift err to test (0-1)'};
dlgtitle = 'ar lag drift coeff';
dims = [1 35];
definput = {'0.25'}; % {'0.25'}; % 
ar_drft_cel =inputdlg(prompt,dlgtitle,dims,definput);

ar_error_drft = cell2mat(cellfun(@str2num,ar_drft_cel,'UniformOutput',false));

disp('Enter autoregr hourly lags in in-flight drift error in dialog box');
prompt = {'Enter hourly ar lags drift err to test 0-1)'};
dlgtitle = 'ar hourly lag drift coeff';
dims = [1 35];
definput = {'0.75'}; % {'0.25'}; % 
ar_h_drft_cel =inputdlg(prompt,dlgtitle,dims,definput);

ar_err_hr_drft = cell2mat(cellfun(@str2num,ar_h_drft_cel,'UniformOutput',false));

disp('Enter std error preferred (inherited) heading in dialog box');
prompt = {'Std inherited heading (degs)'};
dlgtitle = 'std inherited heading (degs)';
dims = [1 35];
definput = {'2.5'}; % {'0'}; % 
std_inh_err_cel =inputdlg(prompt,dlgtitle,dims,definput);

default_inh_err = cell2mat(cellfun(@str2num,std_inh_err_cel,'UniformOutput',false));

% default_inh_err = 2.5; %

% err_det_tests = 20; % 0 %  5 % [15 30]; %[0 15 30]  %

% new Jan 2021 tc clock can either involve local or reference
% sun azimuth adjustments (I believe local makes most sense, i.e., TCSC
% bird adjusts its azimuth according to reference clock but local rate of
% change of sun titlePos(2) = 2e7;

disp('Enter if local or ref azimuth in dialog box');
prompt = {'Local or ref sun azimuth shift?'};
sun_az_list = {'local','reference'};
sun_az_idx = listdlg('PromptString',prompt,'ListString',sun_az_list); 

tc_clock_az = sun_az_list{sun_az_idx};

disp('Enter a plot variable from dialog box');
prompt = {'Enter a plot variable to highlight'};
map_plv_list = {'bias to heading','departure date','current heading', ...
    'current date','compass error','arrival date / time'};
plv_idx = listdlg('PromptString',prompt,'ListString',map_plv_list); 

switch plv_idx
    case 1
        plot_disp_tests = [1 1 1 1 1 1 1 1 1 1]; % 
    case 2
        plot_disp_tests = [5 5 5 5 5 5 5 5 5]; % 
    case 3
         plot_disp_tests = [3 3 3 3 3 3 3 3 3 3]; %   
    case 4
        plot_disp_tests = [4 4 4 4 4 4 4 4 4 4]; %  
    case 5
        plot_disp_tests = [7 7 7 7 7 7 7 7 7 7]; %   
    case 6 
        plot_disp_tests = 8*[1 1 1 1 1 1 1 1 1 1]; %   
end

opt_models = {'fzero','fminbnd'}; %
opt_mod_tests = [1 1 2 1 1 1 1 1 1 1];

options = optimset;  
optnew = optimset(options,'TolX',1e-10);
% option for non-geogr headings (include declinasp_tion) for 
% Lox and magnclinic compass routes

% use mag comp for mag lox and tc comp with mag at night then tranf to sun comp at dawn
incl_decl_tests = [0 1 1 0 0 0 0 0 1 0] == 1; % false %  or_prog_idx <= 2; %  true % 

% regardless, init dir could be 'geographic'
init_geogr_head = true; % false; % 
sun_inher_opt =  false; % true; %  
sun_inh_strs = {'geo','sun'};
sun_inh_str = sun_inh_strs{1+sun_inher_opt};

% Jan 2021 added autocorr error - here the corr coeff for first-order
%  (nightly) lag in detectn error
% ar_error_det =  0.75; %0; % 0.9; % 0.25; %    0.5; %
% option to break up autocorrelation after stopover
% eg if related to weather or finsihed detour
% (less so for magnetoic anomaly unless somehow accounted for using other
% cues)
ar_break_stop_opt = true; % false; %  

plot_optimal_spec_maps = false; %  true; % 
plot_title = false; % true; % 

% small marker size for ease of viewing routes in small (publication) format
small_mrkrs = true;

save_figs_opt = false; % save_opt; %  true; % false; % 

% number ind migrs per angle
n_inds = 10000; % 10 % 300 % 
one_vec = ones(1,n_inds);
 
% set 'other' (non-species-related) options off
full_run = false;
plot_ind_tracks_opt =  false; % true %   

% if ~exist('plot_ms_opt')
%     
%     plot_ms_opt = false; % true; % 
%     
% end
% retain this for post-optimization (defaults to not while optimizing)
plot_ms_opt_save = plot_ms_opt;

%% get species and run specific parameters
plot_close_dist_opt = false; % true; % 
species_goal_opt = true;

% vary or plot headings (3) or departure date (5) or inherited / 
% initial heading (6)
plot_var_test = 3; %  4; %  5 % 
% 4*  6 % 3 %  5 % [6 6 6 5 5]; %
 % 0 = none, 1 = end of or progs per sp, 2 = each ap
plot_cbr_opt = 1; % 0; % 1; %    

std_err_inits = 1; % degs Lat offset 1 % 2.5 % 0; % [-2.5 2.5]; % 0; %

d_alf_init_degs = NaN;
std_err_inits = NaN;
%   5  % 0 %   0.5; % 

plot_titl = false; % 
plot_colbr = false; % true %
plot_cb_end = false; % true %
clr_perf_trx = 'm'; % 'w'; %

if ~exist('std_day')
    std_day = 14*(plot_var_test == 3); % 14 %
end

% test for tighter schedule
% i.e. close number days within mean +/- std 
% eg if std = 10: try 14, std = 5: try 10
n_days_close = std_day;
dalf_close = default_inh_err;

or_pr_strs = {'Geogr_Lox','Geomag_Lox','Magclinic', ...
    'fixed Sun compass','TC Sun reset','TCS no reset', ...
    'TCS nonstop','mix geogr lox TCSC', ...
    'mix geomag lox TCSC','fixed reset'}; %  % Mag Nt TCSC day'}; % ,'TC no reset'};
   
% transparanecy of trajectory markers on map (0 invisible, 1 solid)
transp_val = 0.2; % 1 % 

for l_sp = 1:numel(species_tests)
    
    species = species_tests{l_sp}
    
    get_species_params
    
    if var_sched == 0

        std_n_stop = 0;
        std_day = 0;
        
    elseif var_sched == 2
        
        std_n_stop = 2*std_n_stop;
        std_day = 2*std_day;      

    end
     
    get_species_init_alphs       
    alfs_deg_or_prs = alf_inits;

     %% loop through compass orientation programs

        for i_prog = 1:numel(or_prog_idx)

            incl_decl = incl_decl_tests(or_prog_idx(i_prog));
            transfer_to_sun_dawn = transfer_dawns(or_prog_idx(i_prog));
            

            or_prog_i = or_prog_idx(i_prog);

            or_prog  = or_prog_test(or_prog_i);
            non_stop = non_stops(or_prog_i);
            
            % choose optimization model for program
            % Magcl may be better with fminbnd (keeps bounded, uses abs(distance)
            opt_mod =  'fminbnd'; % opt_models{opt_mod_tests(or_prog_idx(i_prog))}; % 

            %         plot_var_disps = plot_var_disps(or_prog_i);
            % optimize initial heading
            alpha_0 =  alf_inits(or_prog_i); % 60; 

            % simulation model uses plural input or_progs
            or_progs = or_prog_test(or_prog_idx(i_prog));

            for j_err = 1:numel(err_det_tests)

               err_dets = err_det_tests(j_err);        
               err_trans = err_tran_tests(j_err); % (m_tran);
               err_mnts = err_mnt_tests(j_err); % (m_tran);
                
               % % mult by 100 to enourage algorithm to explore larger range
                % alpha_0_100 = alpha_0/100;

               for k_ar = 1:numel(ar_tests)

                   ar_error_det = ar_tests(k_ar);

%                    for m_tran = 1:numel(err_tran_tests)

      
                        fun = @(x) calc_med_displ_comp_rte(x,species, ...
                                or_progs,non_stop,tc_clock_az, ...
                                err_dets,ar_error_det,ar_break_stop_opt, ...
                                err_trans,err_mnts, ar_error_mnt, ...
                               err_drfts,ar_error_drft,ar_err_hr_drft, ...
                                default_inh_err,plot_var_test,init_geogr_head, ...
                                    sun_inher_opt,incl_decl, ...
                                star_night,is_trans_comp,transfer_to_sun_dawn, ...
                               dalf_close, n_days_close,opt_mod,var_sched); 

                        tic
                        % 
                             plot_spec_map_opt = false; % true % 
                             plot_opt = false;
                             disp([species ', ' or_pr_strs{or_prog_idx(i_prog)} ', det err ' num2str(err_dets) ' degs'])

                             if strcmp(opt_mod,'fzero') % && strcmp(species,'Iql Wheatear')
                                  [alph_opt(i_prog,j_err,k_ar),fval,exitflag,output] = fzero(fun,alpha_0,optnew);
                             else % use fminbnd (better for magnetoclinic)
                                 if ~ismember(or_prog,[1 3])
%                                      if or_prog ~= 6
                                         [alph_opt(i_prog,j_err,k_ar),fval,exitflag,output] = fminbnd(fun,max(alpha_0-65,-165), ...
                                         min(alpha_0+65,165),optnew); %  fsolve(fun,alpha_0); %
%                                      else
%                                           [alph_opt(i_prog,j_err,k_ar),fval,exitflag,output] = fminbnd(fun,max(alpha_0-65,-145), ...
%                                          min(alpha_0+65,145),optnew); %                                         
%                                      end
                                 else % geomagn Lox and magnetoclinic restricted to within 90 degrees of (geomagn) South
                                     % adjust bounds for "sensitive" ~E-W species
                                     sens_sp = strcmp(species,'Siberian Willow Warbler') ...
                                         || strcmp(species,'BU Rosefinch') || strcmp(species,'Alk Wheatear');
                                     if alpha_0 > 0
                                         min_alf = (or_prog ==3 && sens_sp)*80;
                                         max_alf = 90 + (or_prog ==3)*5;
                                     else
                                         min_alf = -90 - (or_prog ==3 && sens_sp)*5;
                                         max_alf = -(or_prog ==3 && sens_sp)*80;                            
                                     end
                                     
                                         [alph_opt(i_prog,j_err,k_ar),fval,exitflag,output] = fminbnd(fun,min_alf, ...
                                         max_alf,optnew); %  fsolve(fun,alpha_0); %    

                                 end
                             end

                        toc

                        %%
                        % now run again to plot maps and save output
        %                 if plot_optimal_spec_maps

                            plot_spec_map_opt = plot_optimal_spec_maps; % false; % 
                            plot_opt =  plot_spec_map_opt; % true; %false; %
                            plot_ms_opt = plot_ms_opt_save;
                            % update guess with optimum
                            alf_init = alph_opt(i_prog,j_err,k_ar);
                            dLon_goal = lon_goal_centr - dep_lon_sp_degs;
                            % offset for plotting from generic 180 deg departure in run_no_topo_any
                            spec_lon_offset = 0; % dep_lon_sp_degs-180;
            %                 or_progs = or_prog;
                            plot_var_disp = plot_disp_tests(or_prog_idx(i_prog));
                            run_no_topo_any
                            med_err_test(i_prog,j_err,k_ar) = med_d_close_all;
                            p_gt_500(i_prog,j_err,k_ar) = p_gt_500_km;

                             p_goal(i_prog,j_err,k_ar) = p_within_goal;
                             p_goal_tight_inh(i_prog,j_err,k_ar) = p_goal_tight_inher;
                             med_err_tight_inh(i_prog,j_err,k_ar) = med_d_tight_inher;

                             if is_sun_comp
                                 p_goal_tight_sch(i_prog,j_err,k_ar) = p_goal_tight_sched;
                                 med_err_tight_sch(i_prog,j_err,k_ar) = med_d_tight_sched;
                                 p_goal_tight_all(i_prog,j_err,k_ar) = p_goal_tight_both;
                                 med_err_all(i_prog,j_err,k_ar) = med_d_tight_both;
                             else
                                 p_goal_tight_sch(i_prog,j_err,k_ar) = NaN;
                                 med_err_tight_sch(i_prog,j_err,k_ar) = NaN;
                                 p_goal_tight_all(i_prog,j_err,k_ar) = NaN;
                                 med_err_all(i_prog,j_err,k_ar) = NaN;
                             end

                             plot_spec_map_opt = false; % true % 
                             plot_opt = false;
                             % Get a list of all variables
                             allvars = whos;

                            % Identify the variables that ARE NOT graphics handles. This uses a regular
                            % expression on the class of each variable to check if it's a graphics object
                            tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));

                            if save_opt

                                 dir_nm = ['output_species_runs/' num2str(err_dets) ...
                                      '_' num2str(err_trans) '_' ...
                                      num2str(err_mnts) '_'  num2str(err_drfts) ...
                                     '_deg_errs/' species '/'];
                                [status, msg, msgID] = mkdir(dir_nm);
                                
                                if var_sched == 0 && or_prog >= 6
                                    
                                    schd_str = '_cnst_Schdl';
                                    
                                elseif var_sched == 2 && or_prog >= 6
                                    
                                    schd_str = '_dbl_Schdl';
                                    
                                else
                                    
                                    schd_str = '';
                                    
                                end
                                    

                                % Pass these variable names to save
            %                     save(['output_species_runs/' num2str(err_dets) '_deg_err/' species '_' ...
            %                          or_pr_strs{or_prog_idx(i_prog)}], allvars(tosave).name)               
                                save([dir_nm ...
                                     or_pr_strs{or_prog_idx(i_prog)} '_ars_' num2str(100*ar_error_drft) ...
                                     '_' num2str(100*ar_err_hr_drft) '_' star_night_str '_flt_' trans_comp_str ... 
                                     '_' tc_clock_az '_az_' sun_inh_str '_inh' schd_str '_ind_err'  ...
                                     num2str(10*default_inh_err)], allvars(tosave).name)      

                            end
                        
%                        end

                  end

               end

            end

    %     end
    
end
