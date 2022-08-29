addpath 'D:\Oldenburg_models\generic_comp_mig_model'

% option to cloe interim plots (plus all others!)
 close_opt = true; % false; % 

% paranters for plot routine called 

% which threshold distance to use for performance (fraction arrival)
dist_thr_2022 = 3: % 1 == 100 km, 2 == 250 km, 3 == 500 km, 4 = 1000 km

% minimum threshold for consideration as nonzero: 
% i.e., p(dist) < 1000 > p_thr 
 p_arr_thr_2022 =  0.25; %  0; % 

all_ps_45=nan(721,241,3); % nan(91,61,3); % 

load('Lox_45N_Err_1_mnt_same_drft_15')
main_plots_opts = [0 0 0 1 0]; %
p_thr = p_arr_thr_2022;
i_dist_p_thr = dist_thr_2022;
plot_output_spat_mod
% znew(isnan(znew)) = 0;
ps_Lox_45 = znew; % squeeze(p_close(:,:,3));

load('Mgcl_45N_Err_1_drft_15')
main_plots_opts = [0 0 0 1 0]; %
p_thr = p_arr_thr_2022;
i_dist_p_thr = dist_thr_2022;
plot_output_spat_mod
% znew(isnan(znew)) = 0;
ps_Mgcl_45 = znew; % squeeze(p_close(:,:,3));

load('Fix_45N_Err_1_trans_drft_15_Sep_15')
main_plots_opts = [0 0 0 1 0]; %
p_thr = p_arr_thr_2022;
i_dist_p_thr = dist_thr_2022;
plot_output_spat_mod
% znew(isnan(znew)) = 0;
ps_Fix_45 = znew; % squeeze(p_close(:,:,3));

load('TCSCr_45N_Err_1_trans_drft_15_Sep_15')
main_plots_opts = [0 0 0 1 0]; %
p_thr = p_arr_thr_2022;
i_dist_p_thr = dist_thr_2022;
plot_output_spat_mod
% znew(isnan(znew)) = 0;
ps_TC_45 = znew; % squeeze(p_close(:,:,3));

all_ps_45(:,:,1)=ps_Lox_45;
all_ps_45(:,:,2)=ps_Mgcl_45;
all_ps_45(:,:,3)=ps_Fix_45;
all_ps_45(:,:,4)=ps_TC_45;

% result
min_ps_45=nanmin(all_ps_45,[],3);
max_ps_45=nanmax(all_ps_45,[],3);
% 
% % ps_Lox_45(isnan(ps_Lox_45) & ~isnan(min_ps_45)) = 0;
% % ps_Mgcl_65(isnan(ps_Mgcl_45) & ~isnan(min_ps_45)) = 0;
% % ps_Fix_45(isnan(ps_Fix_45) & ~isnan(min_ps_45)) = 0;
% % ps_TC_45(isnan(ps_TC_45) & ~isnan(min_ps_45)) = 0;

all_ps_65=nan(721,241,3); % nan(91,61,3); % 

load('Lox_Err_1_mnt_same_drft_15_65N')
main_plots_opts = [0 0 0 1 0]; %
p_thr = p_arr_thr_2022;
i_dist_p_thr = dist_thr_2022;
plot_output_spat_mod
% znew(isnan(znew)) = 0;
ps_Lox_65 = znew; % squeeze(p_close(:,:,3));

load('Mgcl_Err_1_mnt_same_drft_15_65N')
main_plots_opts = [0 0 0 1 0]; %
p_thr = p_arr_thr_2022;
i_dist_p_thr = dist_thr_2022;
plot_output_spat_mod
% znew(isnan(znew)) = 0;
ps_Mgcl_65 = znew; % squeeze(p_close(:,:,3));

load('Fix_Err_1_ar_mnt_trans_drft_15_65N_Sep_15')
main_plots_opts = [0 0 0 1 0]; %
p_thr = p_arr_thr_2022;
i_dist_p_thr = dist_thr_2022;
plot_output_spat_mod
% znew(isnan(znew)) = 0;
ps_Fix_65 = znew; % squeeze(p_close(:,:,3));

load('TCSCr_Err_1_mnt_trans_drft_65N_Sep_15')
main_plots_opts = [0 0 0 1 0]; %
p_thr = p_arr_thr_2022;
i_dist_p_thr = dist_thr_2022;
plot_output_spat_mod
% znew(isnan(znew)) = 0;
ps_TC_65 = znew; % squeeze(p_close(:,:,3));

all_ps_65(:,:,1)=ps_Lox_65;
all_ps_65(:,:,3)=ps_Fix_65;
all_ps_65(:,:,4)=ps_TC_65;
all_ps_65(:,:,2)=ps_Mgcl_65;
% result
min_ps_65=nanmin(all_ps_65,[],3);
max_ps_65=nanmax(all_ps_65,[],3);

% ps_Lox_65(isnan(ps_Lox_65) & ~isnan(min_ps_65)) = 0;
% ps_Mgcl_65(isnan(ps_Mgcl_65) & ~isnan(min_ps_65)) = 0;
% ps_Fix_65(isnan(ps_Fix_65) & ~isnan(min_ps_65)) = 0;
% ps_TC_65(isnan(ps_TC_65) & ~isnan(min_ps_65)) = 0;


% YTks = [0 1.1) 1.25) 1.5) 2)]; % ;0:250:1000; % 
% YTkLbs = {'0', '10', '25', '50', '100'};

FgSzX = 250;
FgSzY = 150;
YTks = 100*(-1:0.25:1);

MaxXtraD = 2;

% colorbar limits in percentages
cmin = -75;  %-75; % -1*100; % -0.5
cmax = 75; % 75; % 1*100; % 0.5

FSz = 12;

if close_opt
    close all
end

figure('Position',[250 350 FgSzX FgSzY]); % ,'Color','k'
% del_p_TC_45 = min_ps';
imagesc((-180:d_alf:180),(0:d_err:max_err_val),max_ps_45','AlphaData',double(~isnan(max_ps_45'))) %1+ min_ps')
set(gca, 'ydir', 'Normal','XTick',-180:90:180);
cb = colorbar;
set(cb,'YTick',0:.25:1); % ,
caxis([0 1])
% caxis([0 1000]) % caxis([0 1001)])
% set(cb,'YTick',0:250:1000) %,YTks,'YTickLabel',YTkLbs)
colormap(brewermap([],'*YlOrRd')) % colormap(flipud(parula))
set(gca,'FontSize',FSz)

figure('Position',[350 350 FgSzX FgSzY]); % ,'Color','k'
% del_p_Lox_45 = ps_Lox_45'; % -max_ps_45';
imagesc((-180:d_alf:180),(0:d_err:max_err_val),100*ps_Lox_45','AlphaData',double(~isnan(ps_Lox_45')))
set(gca, 'ydir', 'Normal','XTick',-180:90:180);
caxis([0 100]) % caxis([0 MaxXtraD])
% cb = colorbar;
% set(cb,'YTick',0:25:100); % ,,'YTickLabel',YTkLbs)
colormap(brewermap([],'*YlOrRd')) %,,'*Blues')) % c
set(gca,'FontSize',FSz)

% figure('Position',[350 350 FgSzX FgSzY]); % ,'Color','k'
% del_p_Lox_45 = ps_Lox_45'./min_ps_45'-1;
% imagesc((-180:d_alf:180),(0:d_err:max_err_val),del_p_Lox_45,'AlphaData',double(~isnan(del_p_Lox_45)))
% set(gca, 'ydir', 'Normal');
% % caxis([cmin cmax]) % caxis([0 MaxXtraD])
% cb = colorbar;
% % set(cb,'YTick',YTks,'YTickLabel',YTkLbs)
% colormap(brewermap([],'*YlOrRd')) % c

figure('Position',[450 350 FgSzX FgSzY]); % ,'Col
del_p_Mgcl_45 = ps_Mgcl_45'-ps_Lox_45';
imagesc((-180:d_alf:180),(0:d_err:max_err_val),100*del_p_Mgcl_45,'AlphaData',double(~isnan(del_p_Mgcl_45)))
set(gca, 'ydir', 'Normal','XTick',-180:90:180);
caxis([cmin cmax]) % caxis([0 MaxXtraD])
% cb = colorbar;
% set(cb,'YTick',YTks); % ,'YTickLabel',YTkLbs)
colormap(brewermap([],'*RdYlBu')) % '*YlOrRd')) % ,colormap(flipud(parula))
set(gca,'FontSize',FSz)

figure('Position',[550 350 FgSzX FgSzY]); % ,'Col
del_p_Fix_45 = ps_Fix_45'-ps_Lox_45';
imagesc((-180:d_alf:180),(0:d_err:max_err_val),100*del_p_Fix_45,'AlphaData',double(~isnan(del_p_Fix_45)))
set(gca, 'ydir', 'Normal','XTick',-180:90:180);
caxis([cmin cmax]) % caxis([0 MaxXtraD])
% cb = colorbar;
set(cb,'YTick',YTks); % ,'YTickLabel',YTkLbs)
colormap(brewermap([],'*RdYlBu')) % '*YlOrRd')) % ,'*RdYlBu')) % colormap(flipud(parula))
set(gca,'FontSize',FSz)

figure('Position',[650 350 FgSzX FgSzY]); % ,'Color','k'
del_p_TC_45 = ps_TC_45'-ps_Lox_45'; % ps_TC_45'-min_ps_45');
imagesc((-180:d_alf:180),(0:d_err:max_err_val),100*del_p_TC_45,'AlphaData',double(~isnan(del_p_TC_45)))
set(gca, 'ydir', 'Normal','XTick',-180:90:180);
caxis([cmin cmax]) % caxis([0 MaxXtraD])
% cb = colorbar;
% set(cb,'YTick',YTks); % ,'YTickLabel',YTkLbs)
colormap(brewermap([],'*RdYlBu')) % '*YlOrRd')) % colormap(flipud(parula))
set(gca,'FontSize',FSz)

%% now for 65N

figure('Position',[250 550 FgSzX FgSzY]); % ,'Color','k'
% del_p_TC_45 = min_ps_45';
imagesc((-180:d_alf:180),(0:d_err:max_err_val),max_ps_65','AlphaData',double(~isnan(max_ps_65'))) %1+ min_ps_65')
set(gca, 'ydir', 'Normal','XTick',-180:90:180);
cb = colorbar;
set(cb,'YTick',0:.25:1); % ,
caxis([0 1]) %caxis([0 1000]) % caxis([0 1001)])
% set(cb,'YTick',0:250:1000) %,YTks,'YTickLabel',YTkLbs)
colormap(brewermap([],'*YlOrRd')) %'*RdYlBu')) %  colormap(flipud(parula))
set(gca,'FontSize',FSz)

figure('Position',[350 550 FgSzX FgSzY]); % ,'Color','k'
% del_p_Lox_65 = ps_Lox_65'; % -max_ps_65';
imagesc((-180:d_alf:180),(0:d_err:max_err_val),100*ps_Lox_65','AlphaData',double(~isnan(ps_Lox_65')))
set(gca, 'ydir', 'Normal','XTick',-180:90:180);
caxis([0 100]) % caxis([0 MaxXtraD])
% cb = colorbar;
% set(cb,'YTick',0:25:100); % ,,'YTickLabel',YTkLbs)
colormap(brewermap([],'*YlOrRd')) %,'*Blues')) % c
set(gca,'FontSize',FSz)

% figure('Position',[350 350 FgSzX FgSzY]); % ,'Color','k'
% del_p_Lox_45 = ps_65_Lox_45'./min_ps_65'-1;
% imagesc((-180:d_alf:180),(0:d_err:max_err_val),del_p_Lox_45,'AlphaData',double(~isnan(del_p_Lox_45)))
% set(gca, 'ydir', 'Normal');
% % caxis([cmin cmax]) % caxis([0 MaxXtraD])
% cb = colorbar;
% % set(cb,'YTick',YTks,'YTickLabel',YTkLbs)
% colormap(brewermap([],'*YlOrRd')) % c

figure('Position',[450 550 FgSzX FgSzY]); % ,'Col
del_p_Mgcl_65 = ps_Mgcl_65'-ps_Lox_65';
imagesc((-180:d_alf:180),(0:d_err:max_err_val),100*del_p_Mgcl_65,'AlphaData',double(~isnan(del_p_Mgcl_65)))
set(gca, 'ydir', 'Normal','XTick',-180:90:180);
caxis([cmin cmax]) % caxis([0 MaxXtraD])
% cb = colorbar;
% set(cb,'YTick',YTks); % ,'YTickLabel',YTkLbs)
colormap(brewermap([],'*RdYlBu')) % '*YlOrRd')) % ,'*RdYlBu')) % colormap(flipud(parula))
set(gca,'FontSize',FSz)

figure('Position',[550 550 FgSzX FgSzY]); % ,'Col
del_p_Fix_65 = ps_Fix_65'-ps_Lox_65';
imagesc((-180:d_alf:180),(0:d_err:max_err_val),100*del_p_Fix_65,'AlphaData',double(~isnan(del_p_Fix_65)))
set(gca, 'ydir', 'Normal','XTick',-180:90:180);
caxis([cmin cmax]) % caxis([0 MaxXtraD])
% cb = colorbar;
% set(cb,'YTick',YTks); % ,'YTickLabel',YTkLbs)
colormap(brewermap([],'*RdYlBu')) % '*YlOrRd')) % ,'*RdYlBu')) % colormap(flipud(parula))
set(gca,'FontSize',FSz)

figure('Position',[650 550 FgSzX FgSzY]); % ,'Color','k'
del_p_TC_65 = ps_TC_65'-ps_Lox_65'; % ps_65_TC_45'-min_ps_65');
imagesc((-180:d_alf:180),(0:d_err:max_err_val),100*del_p_TC_65,'AlphaData',double(~isnan(del_p_TC_65)))
set(gca, 'ydir', 'Normal','XTick',-180:90:180);
caxis([cmin cmax]) % caxis([0 MaxXtraD])
% cb = colorbar;
% set(cb,'YTick',YTks); % ,'YTickLabel',YTkLbs)
colormap(brewermap([],'*RdYlBu')) % '*YlOrRd')) % colormap(flipud(parula))
set(gca,'FontSize',FSz)


%% colorbars

figure('Position',[850 350 FgSzX 2*FgSzY]); % ,'Color','k'
caxis([0 100]) % caxis([0 MaxXtraD])
cb = colorbar;
set(cb,'YTick',0:25:100); % ,,'YTickLabel',YTkLbs)
colormap(brewermap([],'*YlOrRd')) %,'*
set(gca,'FontSize',FSz)
axis('off')

figure('Position',[950 350 FgSzX 2*FgSzY]); % ,'Color','k'
caxis([cmin cmax]) % caxis([0 MaxXtraD])
cb = colorbar;
set(cb,'YTick',YTks); % ,'YTickLabel',YTkLbs)
colormap(brewermap([],'*RdYlBu')) % 
axis('off')
set(gca,'FontSize',FSz)