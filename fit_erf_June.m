
clear

% option to plot wth error on x axis
plot_err_opt = false; % true; % 

% option for order in per species plot
% 1 = breadth 2 = n steps 3 = breadth*sqrt(n steps) as per erf function 
% in Normal equivalent
% 4 =  combined factor meridian breadth (secant(Lat)), N_steps 
% & stepwise distance

sort_opt_spec_plot = 4; % 2; % 3; % 1; %

plot_titles_opt = true; % false; % 
titl_letters = true;

mdl_perf = @(b,X) prod_erf_June_sched(b,X);

% e_N_Lat = bs(1);
% e_N_Lon = bs(2);
% e_N_br = bs(3);

% e_g_Lat = bs(4);
% e_g_Lon = bs(5);
% e_g_br = bs(6);

% e_sig_Lat = bs(7);
% e_sig_Lon = bs(8);
% e_sig_br = bs(9);

% % slopes of geo factor n-exp's
% e_bn_Lat = bs(10);
% e_bn_Lon = bs(11);
% e_bn_br = bs(12);

% e_frbr = bs(13);

all_par_nm_s = {'n1','n2','g1','g2','sig1','sig2'};

% disp('Enter test params in dialog box');
% dlgtitle = 'test pars (expons)';
% test_par_prompt = {'n_Lat','n_br','g_Lat','g_br',...
%     'sig_Lat','sig_br','dnds_Lat','dnds_br', ...
%     'frbr','const sun az','exp sched','exp d step'}; %, ...
% % test_par_prompt = {'n_Lat','n_Lon','n_br','g_Lat','g_Lon','g_br',...
% %     'sig_Lat','sig_Lon','sig_br','dnds_Lat','dnds_Lon','dnds_br', ...
% %     'frbr','const sun az','exp sched'}; %, ...
% %     'd(g1)/d(sig)','d(g2)/d(sig)'};
% prompt = {'Choose one or more params'};
% test_pars  = listdlg('PromptString',prompt,'ListString',test_par_prompt);
% N_pars = numel(test_par_prompt);

test_pars = [2 4 8 12];
N_pars = 12;
n_test_pars = numel(test_pars);

Exp_e_vals = [0.5 0 0 1 1  1 0  0 1 0 0 0]; % 0 0]; % [0.5 0.5 0.5 0.5 1 1]; %
Exp_e_vals_TCSC = [0.5 0 0 1 1  1 0 0 1 0 0 0]; %

% null_vals = [0.5 0.5 1 1 0 0];
% Exp_e_vals = [0.5 0.5 0.5 0 1 1 1 1 1 0 0 0 1 0 0]; % 0 0]; % [0.5 0.5 0.5 0.5 1 1]; %
% Exp_e_vals_TCSC = [0.5 0.5 0.5 0 0 1 1 1 1 0 0 0 1 pi/180 1]; %

% new: slopes d(sig)/d(g_Lat) and /d(g_br)
% b_g_Lat = Exp_e_vals(7);
% b_g_br = Exp_e_vals(8);

species_list = ["Nathusius","Kirtlands Warbler","Ring Ouzel","Monarch", ...
 "Hoopoe","BU Rosefinch","Finn Marsh Warbler","GreyCheekThrush","Sib Will Warb S Hem"];
% Bluetail RU"]; % clear

% number of regr options
sz_bs = [1]; % 1 2 3 4]; % [0 1 2 3]; % [1 1 1 2 2 2 2 3]; % [3 3 3]; % 
n_opts = numel(sz_bs);

sq_2 = sqrt(2);

load_rte_strs = ["Geogr_Lox","Geogr_Lox","Geomag_Lox", ...
    "Geomag_Lox","TC Sun reset","TC Sun reset"];
load_trans_strs = ["same","transf"];
load_comp_strs = ["mag","celest"];
xtra_str = '_ars_25_75_';

out_rte_strs = ["GGL","GGL_st","GML","GML_st","TCSC_su","TCSC"];
is_trans_rtes = [0 1 0 1 0 1];
      %  
%
plot_opt = true; % false; %  
save_opt =  false; % true; %

qnts = [0.05 0.25 0.5 0.75 0.95];

% load dist_speed_brdth_frnts_per_sp

calc_goal_dist_migr_sp_spec

n_spec = numel(species_list);

disp('Enter folder name in dialog box');
prompt = {'Folder name data'};
dlgtitle = 'Folder name data';
dims = [1 35];
definput = {'output_const_inh/output_det/'}; %{'output_2_half_inh/'}; %  {'Apr_no_std_inh_sch_or_drft/'}; % {'Apr_no_inh_sch_mnt/'}; % {'0 10 20 30'}; % 
dir_nm_regr_cel = inputdlg(prompt,dlgtitle,dims,definput);
dir_nm_regr = char(dir_nm_regr_cel);

disp('Select compass rtes in dialog box');
prompt = {'Select compass routes'};
dlgtitle = 'Compass routes';
comp_rtes_prompt = {'Geogr Lox','Geogr Lox trans', ...
    'Geomag Lox','Geomag Lox trans','Time-comp Sun','Time-comp Sun trans'};
rte_nrs = listdlg('PromptString',prompt,'ListString',comp_rtes_prompt);
n_rtes = numel(rte_nrs);

out_str_regr = load_rte_strs(rte_nrs);
is_trans_regr = is_trans_rtes(rte_nrs);

is_mag_comp_regr = ismember(rte_nrs, [3 4]) & ~is_trans_regr;
load_comp_regr = load_comp_strs(2-is_mag_comp_regr);

% option for temporal variability in initial departure date and stopover durations
disp('Select if temporal varying schedule in Dialogue box');
prompt = {'Include (sp-appropriate)'; 'variation schedule?'};
sched_list = {'constant','normal std','double std'};
sched_idx = listdlg('PromptString',prompt,'ListString',sched_list); 
% var_sched =  sched_idx == 1; % true; %    false; %  

disp('Enter std inheritance error in dialog box');
prompt = {'Inheritance error (degs)'};
dlgtitle = 'Inher error (degs)';
dims = [1 35];
definput = {'2.5'}; % {'0 10 20 30'}; % 
err_inh_cel = inputdlg(prompt,dlgtitle,dims,definput);
inh_err_regr = cell2mat(cellfun(@str2num,err_inh_cel,'UniformOutput',false));

disp('Enter std compass det error in dialog box');
prompt = {'Std det error (degs)'};
dlgtitle = 'Stepwise comp error (degs)';
dims = [1 35];
definput = {'5 10 20 30 40 50 60'}; % {'20'}; % 
err_det_cel = inputdlg(prompt,dlgtitle,dims,definput);
det_err_regr = cell2mat(cellfun(@str2num,err_det_cel,'UniformOutput',false));

disp('Enter std compass trans error in dialog box');
prompt = {'Std trans error (degs)'};
dlgtitle = 'Stepwise comp error (degs)';
dims = [1 35];
definput = {'0 0 0 0 0 0 0'}; % {'20'}; % {'0 10 20 30'}; % 
err_trans_cel = inputdlg(prompt,dlgtitle,dims,definput);
trans_err_regr = cell2mat(cellfun(@str2num,err_trans_cel,'UniformOutput',false));

disp('Enter std compass mnt error in dialog box');
prompt = {'Std mnt error (degs)'};
dlgtitle = 'Stepwise comp error (degs)';
dims = [1 35];
definput = {'0 0 0 0 0 0 0'}; % {'20'}; % {'0 10 20 30'}; % 
err_mnt_cel = inputdlg(prompt,dlgtitle,dims,definput);
mnt_err_regr = cell2mat(cellfun(@str2num,err_mnt_cel,'UniformOutput',false));

n_base_errs = numel(det_err_regr);

% replicate inh error across base error combis
inh_err_s_rep = (pi/180)*repmat(inh_err_regr,[n_base_errs*n_spec 1]);
% inh_err_2s_rep = trans_array(inh_err_2s_rep);

n_fls_min_rep = repmat(n_hat_fls_min,[n_base_errs 1]);
N_fl_mins = trans_array(n_fls_min_rep);

n_fls_max_rep = repmat(n_hat_fls_max,[n_base_errs 1]);
N_fl_maxs = trans_array(n_fls_max_rep);

% n_fls_max_rep = repmat(max_n_sim,[n_base_errs 1]);
% N_fl_maxs = trans_array(n_fls_max_rep);

% N_fl_mins = n_fls_min_rep;
% N_fl_maxs = n_fls_max_rep;   
    
% fr_br_Lat = 180/pi*repmat(goal_rads./(R_Earth_km*dLat_goals),[n_base_errs 1]);
% fr_br_Lon = 180/pi*repmat(goal_rads./(R_Earth_km*dLon_goals.*cosd(arr_lats)),[n_base_errs 1]);
% fr_br_Lat = trans_array(fr_br_Lat);
% fr_br_Lon = trans_array(fr_br_Lon);

fr_br = repmat(goal_rads./goal_d,[n_base_errs 1]);
fr_br = trans_array(fr_br);

disp('Enter std drft error in dialog box');
prompt = {'Std drft error (degs)'};
dlgtitle = 'Stepwise drft error (degs)';
dims = [1 35];
definput = {'0'}; % {'0 10 20 30'}; % 
err_drft_cel = inputdlg(prompt,dlgtitle,dims,definput);
drft_err_regr = cell2mat(cellfun(@str2num,err_drft_cel,'UniformOutput',false));

% err_det = base_err_regr; % 20; %  
% err_trans = base_err_regr; % 20; %    
% err_mnt = base_err_regr; % 20; %
err_drft = drft_err_regr; % 0; %15; % 
ar_det = 0.25*100; % 0; % 
ar_mnt = 0.75*100;

n_reps = n_spec*n_base_errs;

day_m_d_reps = trans_array(repmat(day_m_d,[n_base_errs 1])); % /R_Earth_km; %  
    
% create 2 cols of zeros for parameter (exponent) selection and default
% values
zero_par_cols = zeros(n_reps,2);

% equiv_ang_err_nontr_s = equiv_err_nontr_s(:);

equiv_err_tr_s = sqrt(repmat(det_err_regr.^2,[n_spec 1]) + ...
    repmat(trans_err_regr.^2,[n_spec 1]) + ...
    mnt_err_regr.^2./(n_mnts'-1) + drft_err_regr^2); % 
equiv_err_tr_s = pi/180*equiv_err_tr_s(:);
% equiv_ang_err_tr_s = equiv_err_tr_s(:);

if mnt_err_regr(1) ~= 0 
    equiv_err_nontr_s = sqrt(repmat(det_err_regr.^2, ...
        [n_spec 1])./n_mnts' + drft_err_regr^2); % ./sqrt(n_hat_fls) 
    equiv_err_nontr_s = pi/180*equiv_err_nontr_s(:);
else
    equiv_err_nontr_s = equiv_err_tr_s;
end

% compute n_steps factor based on n flights and whether trans
% n_flt_Fact_nt = 1-equiv_err_nontr_s.^2/2;
% n_flt_Fact_tr = 1-equiv_err_tr_s.^2/2;

eq_kap_ntr = (180/pi./det_err_regr).^2; % repmat( ,[n_spec 1]);
eq_kap_ntrs = repmat(eq_kap_ntr, [n_spec 1])./n_mnts';
rep_kap_rat_ntr = repmat(besseli(1,eq_kap_ntr)./besseli(0,eq_kap_ntr), ...
    [n_spec 1]);
rep_2_kap_rat_ntr = repmat(besseli(2,eq_kap_ntr)./besseli(0,eq_kap_ntr), ...
    [n_spec 1]);

% if all(n_mnts > 1)
    eq_kap_tr = (180/pi).^2./(det_err_regr.^2./ones(1,n_base_errs) + ...
        trans_err_regr.^2./ones(1,n_base_errs) + ...
        (n_mnts'-1)*mnt_err_regr.^2./n_mnts'.^2);
% else
%     eq_kap_tr = (180/pi).^2./(det_err_regr.^2./ones(1,n_base_errs) + ...
%         trans_err_regr.^2./ones(1,n_base_errs));  
% end

rep_kap_rat_tr = besseli(1,eq_kap_tr)./besseli(0,eq_kap_tr);
rep_2_kap_rat_tr = besseli(2,eq_kap_tr)./besseli(0,eq_kap_tr);

repmat_base = repmat(ones(1,n_base_errs),[n_spec 1]);

% assume "route-mean" sun compass angle (avg of initial and final)
az_gc_regr = az_gc_mn'; % az_gc_1; % 

% ...except use final heading for computing breadth of goal
az_gc_br = repmat(az_gc_2',[n_base_errs 1]); %
% for daily dist maybe 1st head most influential
az_gc_dist = repmat(az_gc_1',[n_base_errs 1]); % 
% for lox always same per definition
az_lox_br = repmat(az_lox',[n_base_errs 1]); % 

% calculate stepwise std error for N-S (cos alf) and E-W (sin alf) components
cos_alf_gc_rep = repmat(cos(az_gc_regr),[n_base_errs 1]);
cos_2_alf_gc_rep = repmat(cos(2*az_gc_regr),[n_base_errs 1]);
sig_c1_gc_rep_tr = sqrt(1+cos_2_alf_gc_rep.*rep_2_kap_rat_tr(:) - ...
    2*(cos_alf_gc_rep.*rep_kap_rat_tr(:)).^2)/sq_2;
sig_c1_gc_rep_ntr = sqrt(1+cos_2_alf_gc_rep.*rep_2_kap_rat_ntr(:) - ...
    2*(cos_alf_gc_rep.*rep_kap_rat_ntr(:)).^2)/sq_2;

cos_alf_lox_rep = repmat(cos(az_lox'),[n_base_errs 1]);
cos_2_alf_lox_rep = repmat(cos(2*az_lox'),[n_base_errs 1]);
sig_c1_lox_rep_tr = sqrt(1+cos_2_alf_lox_rep.*rep_2_kap_rat_tr(:) - ...
    2*(cos_alf_lox_rep.*rep_kap_rat_tr(:)).^2)/sq_2;
sig_c1_lox_rep_ntr = sqrt(1+cos_2_alf_lox_rep.*rep_2_kap_rat_ntr(:) - ...
    2*(cos_alf_lox_rep.*rep_kap_rat_ntr(:)).^2)/sq_2;

sin_alf_gc_rep = repmat(sin(az_gc_regr),[n_base_errs 1]);
sig_s1_gc_rep_tr = sqrt(1-cos_2_alf_gc_rep.*rep_2_kap_rat_tr(:) - ...
    2*(sin_alf_gc_rep .*rep_kap_rat_tr(:)).^2)/sq_2;
sig_s1_gc_rep_ntr = sqrt(1-cos_2_alf_gc_rep.*rep_2_kap_rat_ntr(:) - ...
    2*(sin_alf_gc_rep .*rep_kap_rat_ntr(:)).^2)/sq_2;

sin_alf_lox_rep = repmat(sin(az_lox(:)),[n_base_errs 1]);
sig_s1_lox_rep_tr = sqrt(1-cos_2_alf_lox_rep.*rep_2_kap_rat_tr(:) - ...
    2*(sin_alf_lox_rep.*rep_kap_rat_tr(:)).^2)/sq_2;
sig_s1_lox_rep_ntr = sqrt(1-cos_2_alf_lox_rep.*rep_2_kap_rat_ntr(:) - ...
    2*(sin_alf_lox_rep.*rep_kap_rat_ntr(:)).^2)/sq_2;

% calculate geographic effect from secant(Lat) "log" factor and headings
% cos and sin alf
log_Fact = repmat(ln_fact',[n_base_errs 1]);
log_Fact = log_Fact(:);
% weight secant factor by cos arr Lat
c_arrL_rep = trans_array(repmat(cosd(arr_lats),[n_base_errs 1]));

% geogr factor to account for convergence of meridians regarding
% effects of stepwise, scheduling and inheritance errors 
% geo_fact_lox = sqrt((sin_alf_lox_rep.*log_Fact).^2 + ...
%     cos_alf_lox_rep.^2);
% geo_fact_gc = sqrt((sin_alf_gc_rep.*log_Fact).^2 + ...
%     cos_alf_gc_rep.^2);
% geo_fact_lox = sqrt((sin_alf_lox_rep.*log_Fact.*c_arrL_rep).^2 + ...
%     cos_alf_lox_rep.^2);
% geo_fact_gc = sqrt((sin_alf_gc_rep.*log_Fact.*c_arrL_rep).^2 + ...
%     cos_alf_gc_rep.^2);
% 
% geo_fact_lox = sqrt((cos_alf_lox_rep.*log_Fact).^2 + ...
%     sin_alf_lox_rep.^2);
% geo_fact_gc = sqrt((cos_alf_gc_rep.*log_Fact).^2 + ...
%     sin_alf_gc_rep.^2);
geo_fact_lox = sqrt((sin_alf_lox_rep.*log_Fact).^2 + ...
    cos_alf_lox_rep.^2);
geo_fact_gc = sqrt((sin_alf_gc_rep.*log_Fact).^2 + ...
    cos_alf_gc_rep.^2);

% geo_fact_lox = sqrt((cos_alf_lox_rep.^2.*log_Fact).^2 + ...
%     sin_alf_lox_rep.^4);
% geo_fact_gc =  sqrt((cos_alf_gc_rep.^2.*log_Fact).^2 + ...
%     sin_alf_gc_rep.^4);

% geo_fact_lox = log_Fact; % 
% geo_fact_gc = log_Fact; % 
% 
% geo_fact_lox = log_Fact; % 
% geo_fact_gc = log_Fact; % 

% add colours for GM Lox (orange) and TCSC (green)
addpath('brewer')
cmaps = colormap(brewermap(9,'Set1')); % 6,'Dark2'));
% point to coulours for current "solve" orien % tn prog
% map_id = [5 5 7 2 3]; % [1 1 3 6 4];
clr_gm = cmaps(5,:);
clr_gg = cmaps(4,:); % cmaps(9,:); % cmaps(8,:); %
clr_tc = cmaps(3,:);
clr_gr = cmaps(9,:); % 

for k_rte = 1:n_rtes 
    
    rte_k = rte_nrs(k_rte);
    
    TC_rte = rte_k >= 5;
     
    sig_alfs = is_trans_regr(k_rte)*equiv_err_tr_s + ...
        ~is_trans_regr(k_rte)*equiv_err_nontr_s;
    
    if sched_idx == 1 && TC_rte % var_inh_err_2s_rep 
        schd_str = '_cnst_Schdl';
        sig_dates = zeros(n_base_errs*n_spec,1);
    elseif sched_idx == 2 || ~TC_rte 
        schd_str = '';
        sig_dates = repmat(sig_date,[n_base_errs 1]);
        sig_dates = trans_array(sig_dates);
    else
        schd_str = '_dbl_Schdl';
        sig_dates = repmat(2*sig_date,[n_base_errs 1]);
        sig_dates = trans_array(sig_dates);
    end
    
%     sig_cos = is_trans_regr(k_rte)*(~TC_rte*sig_c1_lox_rep_tr + ...
%         TC_rte*sig_c1_gc_rep_tr) + ~is_trans_regr(k_rte)* ...
%         (~TC_rte*sig_c1_lox_rep_ntr + TC_rte*sig_c1_gc_rep_ntr);
%         
%     sig_sin = is_trans_regr(k_rte)*(~TC_rte*sig_s1_lox_rep_tr + ...
%         TC_rte*sig_s1_gc_rep_tr) + ~is_trans_regr(k_rte)* ...
%         (~TC_rte*sig_s1_lox_rep_ntr + TC_rte*sig_s1_gc_rep_ntr);
    
    cos_alfs = ~TC_rte*cos_alf_lox_rep + TC_rte*cos_alf_gc_rep;
    cos_2alfs = ~TC_rte*cos_2_alf_lox_rep + TC_rte*cos_2_alf_gc_rep;   
    sin_alfs = abs(~TC_rte*sin_alf_lox_rep + TC_rte*sin_alf_gc_rep);  
    
%     bes_rat = is_trans_regr(k_rte)*rep_kap_rat_tr + ...
%         ~is_trans_regr(k_rte)*rep_kap_rat_ntr;
%     bes_rat = bes_rat(:);

    eq_kaps = is_trans_regr(k_rte)*eq_kap_tr + ...
            ~is_trans_regr(k_rte)*eq_kap_ntrs;
    
%     p_arrs_Lats = is_trans_regr(k_rte)*(~TC_rte*p_arrs_Lats_tr_lox + ...
%         TC_rte*p_arrs_Lats_tr_lox) + ~is_trans_regr(k_rte)* ...
%         (~TC_rte*p_arrs_Lats_ntr_lox + TC_rte*p_arrs_Lats_ntr_gc);
    
    eq_kap = is_trans_regr(k_rte)*eq_kap_tr(:) + ~is_trans_regr(k_rte)*eq_kap_ntrs(:);
    
    geo_Facts  = TC_rte*geo_fact_lox + ~TC_rte*geo_fact_gc;
    
         % and final heading for goal breadth 
    goal_head = TC_rte*az_gc_br + ~TC_rte*az_lox_br;

    dist_head = abs(TC_rte*az_gc_dist + ~TC_rte*az_lox_br);
%     day_m_d_fact = day_m_d_reps./median(day_m_d_reps); % /R_Earth_km; % 
  
    day_m_d_fact = day_m_d_reps./median(day_m_d_reps); % /R_Earth_km; % 
    
%     day_m_d_fact = day_m_d_reps.*sin_alfs./ ...
%         median(day_m_d_reps.*sin_alfs); % /R_Earth_km; % 

%     day_m_d_fact = day_m_d_reps.*sin(dist_head)./ ...
%         median(day_m_d_reps.*sin(dist_head)); % /R_Ea
    
    % June 2021 - add daily dist for TCSC routes
%     day_mds = TC_rte*day_m_d_reps + ~TC_rte*repmat_base(:);
     day_mds = TC_rte*day_m_d_fact + ~TC_rte*repmat_base(:);
     

     
%     rte{k_rte}.X = [fr_br_Lat fr_br_Lon cos_alfs sin_alfs bes_rat ...
%         sig_cos sig_sin geo_Fact N_fl_mins N_fl_maxs equiv_err_tr_s ...
%         fr_br p_arrs_Lats];  
      
% include inher std and log Fact for compass 'bias' (non-stepwise reduced)
    rte{k_rte}.X = [fr_br N_fl_mins N_fl_maxs sig_alfs ...
          geo_Facts cos_alfs cos_2alfs c_arrL_rep.^2 ...
          inh_err_s_rep sig_dates log_Fact day_mds goal_head zero_par_cols];    
                    
    for j_err = 1:n_base_errs

        err_det = det_err_regr(j_err); % 20; %  
        err_trans = trans_err_regr(j_err); % 20; %    
        err_mnt = mnt_err_regr(j_err); % 20; %

        for isp = 1:n_spec
          
%             is_trans_regr load_trans_strs load_rte_strs
            try
                
                if isp ~= 4 || k_rte ~= 3
                    outpt = load(['output_species_runs/' dir_nm_regr num2str(err_det) '_' ...
                    num2str(err_trans) '_' num2str(err_mnt) '_' num2str(err_drft) ...
                    '_deg_errs/' species_list{isp} '/' load_rte_strs{rte_nrs(k_rte)} ...
                    xtra_str load_comp_regr{k_rte} '_flt_' load_trans_strs{1+is_trans_regr(k_rte)} ...
                    '_local_az_geo_inh' schd_str]); % 
                else % Monarch not transferred
                    outpt = load(['output_species_runs/' dir_nm_regr num2str(err_det) '_' ...
                    num2str(err_trans) '_' num2str(err_mnt) '_' num2str(err_drft) ...
                    '_deg_errs/' species_list{isp} '/' load_rte_strs{rte_nrs(k_rte)} ...
                    xtra_str load_comp_regr{k_rte} '_flt_same' ...
                    '_local_az_geo_inh' schd_str]); %                    
                end
            catch
                keyboard
            end
            rte{k_rte}.out{j_err,isp} = outpt;
            rte{k_rte}.p_arr(j_err,isp) = sum(outpt.dec_close{1}<=goal_rads(isp) & ...
                    outpt.arr_ts{1}' <= n_hat_fls_max(isp))/outpt.n_inds;  % outpt.p_within_goal;
            rte{k_rte}.med_d(j_err,isp) = outpt.med_d_close_all;
            
            arrd = outpt.dec_close{1}<=goal_rads(isp);
            rte{k_rte}.med_N_arr(j_err,isp) = nanmean(outpt.arr_ts{1}(arrd));

            rte{k_rte}.alf_inits(j_err,isp) = outpt.alf_init;        
            rte{k_rte}.cos_alf(j_err,isp) = abs(cosd(outpt.alf_init)); % abs(tand(GGL.alf_init/2)); % abs(cosd(GGL.alf_init)/cosd(GGL.dep_lat_degs));
            rte{k_rte}.sin_alf(j_err,isp) = sqrt(1-rte{k_rte}.cos_alf(j_err,isp)^2);
            
            rte{k_rte}.med_nfl(j_err,isp) = nanmean(outpt.arr_ts{1})-1;
            rte{k_rte}.qnts_nfl{j_err,isp} = quantile(outpt.arr_ts{1},qnts);       
            
            rte{k_rte}.alf_inits(j_err,isp) = outpt.alf_init;        
            rte{k_rte}.cos_alf(j_err,isp) = abs(cosd(outpt.alf_init)); % abs(tand(GGL.alf_init/2)); % abs(cosd(GGL.alf_init)/cosd(GGL.dep_lat_degs));
            rte{k_rte}.sin_alf(j_err,isp) = sqrt(1-rte{k_rte}.cos_alf(j_err,isp)^2);
            
            rte{k_rte}.med_nfl(j_err,isp) = nanmedian(outpt.arr_ts{1});    
            
        end
                    
    end
    
            rte{k_rte}.Y = trans_array(rte{k_rte}.p_arr);
            
            fit_all_combos

end

if plot_opt
    
    
%     geo_mod_fact =  sqrt(cos(goal_head).^2./c_arrL_rep.^2 + sin(goal_head).^2);
    
    if plot_err_opt
         plot_regr_June
    end
    
    plot_fit_per_spec_June
 
% plot gain per species and compass precision (detection error)
n_errs = numel(det_err_regr);
n_e_m1 = n_errs - 1;

for isp = 1:n_spec
    
    all_gain_TCSC(isp,1:n_errs) = 100*(rte{3}.Y(isp:n_spec:isp+n_e_m1*n_spec)./ ...
        rte{1}.Y(isp:n_spec:isp+n_e_m1*n_spec)-1);
    all_gain_TCSC_GM(isp,1:n_errs) = 100*(rte{3}.Y(isp:n_spec:isp+n_e_m1*n_spec)./ ...
        rte{2}.Y(isp:n_spec:isp+n_e_m1*n_spec)-1);
    
    med_gain_TCSC(isp) = median(100*(rte{3}.Y(isp:n_spec:end)./ ...
        rte{1}.Y(isp:n_spec:end)-1));
    
    med_gain_TCSC_GM(isp) = median(100*(rte{3}.Y(isp:n_spec:end)./ ...
        rte{2}.Y(isp:n_spec:end)-1));
    
end

CrSz = 60; % 70; %
DndSz =  135; % 90; %
StrSz = 70;

n_hat_fls(5) = n_hat_fls(5) + 0.5;

   H = figure('Position',[400 200 215 350]); % ,clr_gry); %[0.925 0.975 1]) % [.95 .9 .8])

       hold
    
        h = scatter(n_hat_fls([1:3 5:n_spec]), ...
            geo_fact_gc([1:3 5:n_spec]) ... %sqrt(n_hat_fls) ./geo_mod_fact(1:9)' .*ln_fact
       ,DndSz,med_gain_TCSC_GM([1:3 5:n_spec]),'d','fill'); % a
       h.MarkerEdgeColor = 'k'; 
       
         h2 = scatter(n_hat_fls([1:3 5:n_spec]), ... % +0.5
            geo_fact_gc([1:3 5:n_spec])  ... %sqrt(n_hat_fls) ./geo_mod_fact(1:9)' .*ln_fact +.0125
       ,StrSz,med_gain_TCSC([1:3 5:n_spec]),'h','fill'); % a
       h2.MarkerEdgeColor = 'k'; 
           
%     set(gca,'YTick',sqrt([7 15 25 40]),'YTickLabel',{'7','15','25','40'})
%     xlabel('stepwise distance (km)') % ,'Rotation',50)
    xlabel('Minimum number steps') % ylabel
    ylabel('Spherical geometry factor') % zlabel
    
    set(gca,'YAxisLocation','Right')
  
    cax = caxis;
 
    if cax(1) < 0
                
         colormap(brewermap([],'*RdYlBu')) % YlOrRd'))
         
%           hold

          h4 = scatter(n_hat_fls(4), ... % +0.5
            geo_fact_gc(4)  ... %sqrt(n_hat_fls) ./geo_mod_fact(1:9)' .*ln_fact +.0125
                   ,DndSz,30.83,'d','fill'); % TCSC 5-25 biol vs GML incl 15 drft 29.7314 % p_s in output_biol
           h4.MarkerEdgeColor = 'k'; 
           
          h3 = scatter(n_hat_fls(4), ...
            geo_fact_gc(4) ... %sqrt(n_hat_fls) ./geo_mod_fact(1:9)' .*ln_fact +.0125
                        ,StrSz,46.01,'h','fill'); % TCSC 5-25 biol vs GGL incl 15 drft 9.6042
            h3.MarkerEdgeColor = 'w'; 
                    

            axis([0 45 0.95 1.3])
         set(gca,'XTick',0:15:45,'YTick',1:.1:1.3)
        caxis([-35 35])
        
    else
        
       colormap(brewermap([],'YlOrRd'))        
    
    end
    
     figure('Position',[400 200 135 160]); % ,c
  
     cb = colorbar('North');
     title(cb,{'Median gain (%)'},'FontSize',9) % 
     caxis(cax)
     
     if cax(1) < 0
         
         set(cb,'XTick',-50:10:50,'FontSize',7)
         colormap(brewermap([],'*RdYlBu')) % 'YlOrRd'))
         caxis([-35 35])
         
%          
%          ax = axes('Position',cb.Position);
%          tickangle(ax,0)

     else
         
         set(cb,'XTick',0:15:55)
         colormap(brewermap([],'YlOrRd')) % 'YlOrRd'))
       
     end
     

     %         title('J','FontSize',FSz)
     axis off
     box off
%      tight

figure
scatter(fr_br(1:9).*sqrt(N_fl_mins(1:9)), ...
    rte{2}.Y(19:27),90,clr_gm,'d','fill')
hold
scatter(fr_br(1:9).*sqrt(N_fl_mins(1:9)), ...
    rte{3}.Y(19:27),120,clr_tc,'h','fill')

     
end

% end 

%% save

if save_opt
    
       allvars = whos;

       tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));

       save('regr_output_biol',allvars(tosave).name,'-v7.3')
    
end

