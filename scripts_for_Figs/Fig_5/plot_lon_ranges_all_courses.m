addpath 'D:\Oldenburg_models\generic_comp_mig_model'

addpath('D:\Oldenburg_models\generic_comp_mig_model\brewer')
addpath('D:\Oldenburg_models\geomagn_orientn_model\Violinplot-Matlab-master')

cmaps = colormap(brewermap(7,'Set1')); % 6,'Dark2'));

% point to coulours for current "solve" orien % tn prog
% map_id = [5 5 7 2 3]; % [1 1 3 6 4];
clr_Lox = cmaps(5,:); % use mag for dipole, orange looks better :)
% clr_Lox = cmaps(4,:);
clr_tc = cmaps(3,:);
clr_fx = cmaps(2,:);
clr_mcl = cmaps(7,:);

% option to cloe interim plots (plus all others!)
 close_opt = true; % false; % 

% paranters for plot routine called 

% which threshold distance to use for performance (fraction arrival)
plot_thr = 3; % 3 is 500 km, 4 is 1000 km ; 

plot_Lons_opt = 1 % 2 % 3 %

% minimum threshold for consideration as nonzero: 
% i.e., p(dist) < 1000 > p_thr 
 p_thr = 0.25; %  0; %  
 
load_opt = 1 % 0 %  
 
if load_opt == 1 && ~exist('ps_Lox_65')
    if p_thr == 0.25
       if plot_thr == 3
         load all_courses_mnt_drft_500km_p_gt_025
       elseif plot_thr == 4
         load all_courses_mnt_drft_1000km_p_gt_025
       else % == 1
         load all_courses_mnt_drft_100km_p_gt_025
     end
    elseif p_thr == 0
      load all_courses_mnt_drft_500km_p_gt_0
    end
end

figure('Position',[300 100 220 210]); hold
plot_lon_range(ps_Lox_65,clr_Lox, ':') %'--'); %
plot_lon_range(ps_Fix_65,clr_fx,':') %'--'); %
plot_lon_range(ps_Mgcl_65,clr_mcl,':') % '--'); %
plot_lon_range(ps_TC_65,clr_tc,':') % '--'); %
plot_lon_range(ps_Lox_45,clr_Lox,'-') % '--'); %
plot_lon_range(ps_Fix_45,clr_fx,'-') % '--'); %
plot_lon_range(ps_Mgcl_45,clr_mcl,'-') % '--'); %
plot_lon_range(ps_TC_45,clr_tc,'-') % '--'); %

xlabel('Compass precision (^o)')
ylabel('Maximum longitude range (^o)')
box on
set(gca,'FontSize',9)
title('a','FontSize',10)
% title({'Compass course ranges','with biologically-relevant error'})

% plot distribution of fractional gains for 20 degree error = 81st columns
err_plots = [15 30];
letts = {'d','e'};

% longitudes are stored from -180 to 180 in 0.5 deg steps
lon_vec = -180:0.5:180;

% find indices of stats for simulations with long ranges up to 60 degs
% this assumes size is 721, i.e., per half degree
% so first and last 240 indicies will be 60.5-180
long_Lon = 45 % 60 % 
long_Lon_ch = num2str(long_Lon);
fst_short_Lon = find(lon_vec >= -long_Lon,1,'first');
lst_short_Lon = 721-fst_short_Lon; 
length_Short = lst_short_Lon-fst_short_Lon+1;
length_Short_p1 = length_Short+1;
length_Long_p1 = 721 - length_Short+1;
dLon_S = fst_short_Lon:lst_short_Lon;
dLon_L = [1:(fst_short_Lon-1) (lst_short_Lon+1):721];

for i_err = 1:numel(err_plots)
    
    err_plot = err_plots(i_err);
    idx_pl = 1 + 4*err_plot;
    
    p_TC_err_45 = ps_TC_45(:,idx_pl);
    
    p_TC_45_Short = ps_TC_45(dLon_S,idx_pl);   
    p_TC_45_Long = ps_TC_45(dLon_L,idx_pl);   
    gains_45 = 100*[(1-ps_Lox_45(:,idx_pl)./p_TC_err_45) (1-ps_Mgcl_45(:,idx_pl)./p_TC_err_45) (1-ps_Fix_45(:,idx_pl)./p_TC_err_45)];
    gns_45_lon_S = 100*[(1-ps_Lox_45(dLon_S,idx_pl)./p_TC_45_Short) ...
        (1-ps_Mgcl_45(dLon_S,idx_pl)./p_TC_45_Short) (1-ps_Fix_45(dLon_S,idx_pl)./p_TC_45_Short)];
    gns_45_lon_L = 100*[(1-ps_Lox_45(dLon_L,idx_pl)./p_TC_45_Long) ...
        (1-ps_Fix_45(dLon_L,idx_pl)./p_TC_45_Long)];
    % pad the <= and > 60 arrays with NaNs for easiest use with boxplot
    % (don't have to concatenate and use separate grouping if all same
    % length)
    gns_45_lon_S(length_Short_p1:721,:) = NaN;
    gns_45_lon_L(length_Long_p1:721,:) = NaN;  
        
    p_TC_err_65 = ps_TC_65(:,idx_pl);
    p_TC_65_Short = ps_TC_65(dLon_S,idx_pl);   
    p_TC_65_Long = ps_TC_65(dLon_L,idx_pl);   
    gains_65 = 100*[(1-ps_Lox_65(:,idx_pl)./p_TC_err_65) (1-ps_Mgcl_65(:,idx_pl)./p_TC_err_65) (1-ps_Fix_65(:,idx_pl)./p_TC_err_65)];
    gns_65_lon_S = 100*[(1-ps_Lox_65(dLon_S,idx_pl)./p_TC_65_Short) ...
        (1-ps_Mgcl_65(dLon_S,idx_pl)./p_TC_65_Short) (1-ps_Fix_65(dLon_S,idx_pl)./p_TC_65_Short)];
    gns_65_lon_L = 100*[(1-ps_Lox_65(dLon_L,idx_pl)./p_TC_65_Long) ...
        (1-ps_Fix_65(dLon_L,idx_pl)./p_TC_65_Long)];
    % pad the <= and > 60 arrays with NaNs for easiest use with boxplot
    % (don't have to concatenate and use separate grouping if all same
    % length)
    gns_65_lon_S(length_Short_p1:721,:) = NaN;
    gns_65_lon_L(length_Long_p1:721,:) = NaN;  
    
    perf_all = [ps_Lox_45(:,idx_pl) ps_Lox_65(:,idx_pl) ...
               ps_Mgcl_45(:,idx_pl) ps_Mgcl_65(:,idx_pl)  ...
               ps_Fix_45(:,idx_pl) ps_Fix_65(:,idx_pl) ...
               ps_TC_45(:,idx_pl) ps_TC_65(:,idx_pl)];
    
    gains_all = [gains_45(:,1) gains_65(:,1) gains_45(:,2) ...
        gains_65(:,2) gains_45(:,3) gains_65(:,3)];
    
    gains_all(abs(gains_all) > 1e3) = NaN;
    
    gains_all_Lon_cats = [gns_45_lon_S(:,1) gns_45_lon_L(:,1) ...
        gns_65_lon_S(:,1) gns_65_lon_L(:,1) ...
        gns_45_lon_S(:,2) gns_65_lon_S(:,2)  ...
        gns_45_lon_S(:,3) gns_45_lon_L(:,2) ...
        gns_65_lon_S(:,3) gns_65_lon_L(:,2)];
    
    deg = char(176);
    Delta =  char(916);
    MedD_lab = ['45' deg 'N']; % ['45' deg '-25' deg];
    MedD_lab_short_Lon = ['45' deg 'N, ' Delta 'Lon <' long_Lon_ch deg];
    lab_long_Lon = [Delta 'Lon >' long_Lon_ch deg];
    LonD_lab = ['65' deg 'N']; %  '-0' deg];
    LonD_lab_short_Lon = ['65' deg 'N, ' Delta 'Lon <' long_Lon_ch deg];
%     LonD_lab_long_Lon = ['65' deg 'N, ' Delta 'Lon > 60' deg];
    Labs = {MedD_lab, LonD_lab, MedD_lab, LonD_lab, MedD_lab, LonD_lab};
    Labs_perf = {MedD_lab, LonD_lab, MedD_lab, LonD_lab, MedD_lab, ...
        LonD_lab ,MedD_lab, LonD_lab};
    Labs_Lon_cats = {MedD_lab_short_Lon, lab_long_Lon, ...
        LonD_lab_short_Lon lab_long_Lon, ...
        MedD_lab, LonD_lab, ...
        MedD_lab_short_Lon, lab_long_Lon, ...
        LonD_lab_short_Lon lab_long_Lon};
%     Labs = {'Lox 45','Lox 65','Mgcl 45','Mgcl 65','Fix 45','Fix 65'};
    Cols = [clr_Lox;clr_Lox;clr_mcl;clr_mcl;clr_fx;clr_fx]; % ['r','r','k','k','b','b'];
    Cols_perf = ['r','r','k','k','b','b','g','g'];
%     Cols_Lons60 = ['r','r','r','r','k','k','b','b','b','b'];
    Cols_Lons60 = [clr_Lox;clr_Lox;clr_Lox;clr_Lox; ...
        clr_mcl;clr_mcl;clr_fx;clr_fx;clr_fx;clr_fx];
    figure('Position',[300 100 230 200]); hold
%     h = boxplot(gains_all,'Labels',Labs,'Colors',Cols, ...
%          'Notch','on','Symbol','x');
    if plot_Lons_opt == 1
        h = violinplot(gains_all,Labs,'ViolinColor',Cols, ...
            'Width',0.45,'ShowData',false,'ViolinAlpha',1); % ,'ShowMean',true); % ,'Labels',Labs,'Colors',Cols, ...
        set(gca,'FontSize',9)
        ylim([-75 85]) % ylim([-125 125])
        xlim([0.5 6.5])
%      set(gca,'XTickLabel',{'1','2',['$(\Delta)$ $(\le)$ 30'],'4','5','6','7','8'},'interpreter','latex')
    elseif plot_Lons_opt == 2
        h = boxplot(gains_all_Lon_cats,'Labels',Labs_Lons60,'Colors',Cols_Lons60, ...
         'Notch','on','Symbol','x','Whisker',1);
         ylim([-150 100])
    elseif plot_Lons_opt == 3
           h = boxplot(perf_all,'Labels',Labs_perf,'Colors',Cols_perf, ...
             'Notch','on','Symbol','x','Whisker',1);
        ylim([0 1])
    else
        h = boxplot(perf_all_Lon_cats, ... 'Labels',Labs_perf,'Colors',Cols_perf
             'Notch','on','Symbol','x','Whisker',1);
        ylim([0 1])           
    end
     xtickangle(45)
%      set(h,'LineWidth',0.75)

     if i_err == 2
         
         set(gca,'YAxisLocation','Right')

     end
     box on
    % xlabel('substep error magnitudes (^o)')
    ylabel('Performance gain TCSC (%)')
    
    ttl = title(letts{i_err},'FontSize',12);
    ttl.Units = 'Normalize'; 
    ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
    ttl.HorizontalAlignment = 'left'; 
%     title(letts{i_err},'FontSize',10)
%     title({'Relative gain TCSC',['with ' num2str(err_plot) ' biologically-relevant error']})
    
end

idx_errs = 1 + 4*err_plots;
LStys = {'-',':'}; % '--'};
LW = 1.75;
bw = 0.5; % 
pl_w = 1;
bw_res = pl_w/bw;

letts2 = {'b','c'};

for i_err = 1:numel(err_plots)
    
    figure('Position',[600 100 220 110]); hold
    idx_pl = idx_errs(i_err);
    
    sm_TC = ps_TC_65(:,idx_pl); % smooth(ps_TC_65(:,idx_pl),'rlowess'); % movmean(ps_TC_65(:,idx_pl),bw); %  
    sm_TC(isnan(sm_TC)) = 0;
    lst_z0 = find(sm_TC(1:360) == 0,1,'last');
    sm_TC(1:lst_z0-1) = NaN;
    fst_z0 = find(sm_TC(361:end) == 0,1,'first');
    sm_TC(361+fst_z0:end) = NaN;  
    plot(-180:pl_w/2:180,sm_TC(1:pl_w:end), ...
        'Color',clr_tc,'LineStyle',LStys{2},'LineWidth',LW) 
    
    sm_TC = ps_TC_45(:,idx_pl); % smooth(ps_TC_65(:,idx_pl),'rlowess'); % movmean(ps_TC_65(:,idx_pl),bw); %  
    sm_TC(isnan(sm_TC)) = 0;
    lst_z0 = find(sm_TC(1:360) == 0,1,'last');
    sm_TC(1:lst_z0-1) = NaN;
    fst_z0 = find(sm_TC(361:end) == 0,1,'first');
    sm_TC(361+fst_z0:end) = NaN;  
    plot(-180:pl_w/2:180,sm_TC(1:pl_w:end), ...
        'Color',clr_tc,'LineStyle',LStys{1},'LineWidth',LW) 
    
    sm_Lox = ps_Lox_65(:,idx_pl); % smooth(ps_Lox_65(:,idx_pl),'rlowess'); % movmean(ps_Lox_65(:,idx_pl),bw); % smooth( )
    sm_Lox(isnan(sm_Lox)) = 0;
    lst_z0 = find(sm_Lox(1:360) == 0,1,'last');
    sm_Lox(1:lst_z0-1) = NaN;
    fst_z0 = find(sm_Lox(361:end) == 0,1,'first');
    sm_Lox(361+fst_z0:end) = NaN;  
    plot(-180:pl_w/2:180,sm_Lox(1:pl_w:end), ...
        'Color',clr_Lox,'LineStyle',LStys{2},'LineWidth',LW)
  
    sm_Lox = ps_Lox_45(:,idx_pl); % smooth(ps_Lox_65(:,idx_pl),'rlowess'); % movmean(ps_Lox_65(:,idx_pl),bw); % smooth( )
    sm_Lox(isnan(sm_Lox)) = 0;
    lst_z0 = find(sm_Lox(1:360) == 0,1,'last');
    sm_Lox(1:lst_z0-1) = NaN;
    fst_z0 = find(sm_Lox(361:end) == 0,1,'first');
    sm_Lox(361+fst_z0:end) = NaN;  
    plot(-180:pl_w/2:180,sm_Lox(1:pl_w:end), ...
        'Color',clr_Lox,'LineStyle',LStys{1},'LineWidth',LW)
    
    sm_Mgcl = ps_Mgcl_65(:,idx_pl); % smooth(ps_Mgcl_65(:,idx_pl),'rlowess'); % movmean(ps_Mgcl_65(:,idx_pl),bw); %
    sm_Mgcl(isnan(sm_Mgcl)) = 0;
    lst_z0 = find(sm_Mgcl(1:360) == 0,1,'last');
    sm_Mgcl(1:lst_z0-1) = NaN;
    fst_z0 = find(sm_Mgcl(361:end) == 0,1,'first');
    sm_Mgcl(361+fst_z0:end) = NaN; 
    plot(-180:pl_w/2:180,sm_Mgcl(1:pl_w:end), ...
        'Color',clr_mcl,'LineStyle',LStys{2},'LineWidth',LW) 
  
    sm_Mgcl = ps_Mgcl_45(:,idx_pl); % smooth(ps_Mgcl_65(:,idx_pl),'rlowess'); % movmean(ps_Mgcl_65(:,idx_pl),bw); %
    sm_Mgcl(isnan(sm_Mgcl)) = 0;
    lst_z0 = find(sm_Mgcl(1:360) == 0,1,'last');
    sm_Mgcl(1:lst_z0-1) = NaN;
    fst_z0 = find(sm_Mgcl(361:end) == 0,1,'first');
    sm_Mgcl(361+fst_z0:end) = NaN; 
    plot(-180:pl_w/2:180,sm_Mgcl(1:pl_w:end), ...
        'Color',clr_mcl,'LineStyle',LStys{1},'LineWidth',LW) 
    
    sm_Fix = smooth(ps_Fix_65(:,idx_pl),0.05,'rlowess'); % movmean(ps_Fix_65(:,idx_pl),bw); % smooth( )
    sm_Fix(isnan(sm_Fix)) = 0;
    lst_z0 = find(sm_Fix(1:360) == 0,1,'last');
    sm_Fix(1:lst_z0-1) = NaN;
    fst_z0 = find(sm_Fix(361:end) == 0,1,'first');
    sm_Fix(361+fst_z0:end) = NaN; 
    plot(-180:pl_w/2:180,sm_Fix(1:pl_w:end), ...
        'Color',clr_fx,'LineStyle',LStys{2},'LineWidth',LW)
    
    sm_Fix = smooth(ps_Fix_45(:,idx_pl),0.05,'rlowess'); % movmean(ps_Fix_65(:,idx_pl),bw); % smooth( )
    sm_Fix(isnan(sm_Fix)) = 0;
    lst_z0 = find(sm_Fix(1:360) == 0,1,'last');
    sm_Fix(1:lst_z0-1) = NaN;
    fst_z0 = find(sm_Fix(361:end) == 0,1,'first');
    sm_Fix(361+fst_z0:end) = NaN; 
    plot(-180:pl_w/2:180,sm_Fix(1:pl_w:end), ...
        'Color',clr_fx,'LineStyle',LStys{1},'LineWidth',LW)
    
    xlim([-180 180])
    set(gca,'XTick',-180:45:180)
    xtickangle(45)
    
    set(gca,'YAxisLocation','Right')
    ylim([0 1])
    set(gca,'YTick',0:0.25:1,'YTickLabels',{'0','','0.5','','1'})
    
    set(gca,'FontSize',9)

    ttl = title(letts2{i_err},'FontSize',12);
    ttl.Units = 'Normalize'; 
    ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
    ttl.HorizontalAlignment = 'left'; 
%     title(letts2{i_err},'FontSize',10)
    box on
    
end

figure('Position',[600 100 250 100]); 
ylabel('Performance','FontSize',9)
set(gca,'YAxisLocation','Right')
 
figure('Position',[600 100 250 100]); 
xlabel('Longitude distance (^o)','FontSize',9)