
addpath('D:\Oldenburg_models\generic_comp_mig_model\brewer')

gg_str = 'Geogr_Lox';
gm_str = 'Geomag_Lox';
tc_str = 'TC Sun reset';
f_str_1_no_drf = '_ars_25_75_'; % '_ars_0_75_';
f_str_1_drf = '_ars_25_75_';
f_str_2 = '_local_az_geo_inh';
f_str_3 = '_local_az_geo_inh_cnst_Schdl';
tr_str = '_flt_transf';
sm_str = '_flt_same';

species = 'Finn Marsh Warbler'; %  'Kirtlands Warbler'; % 'Ring Ouzel'; % 'Nathusius'; % 'GreyCheekThrush'; %'Sib Will Warb S Hem'; %  
err_dets =   [5 10 20 30 40]; % [5 10 15 20 25 30]; %0:5:40; %17.5 
err_mnts =  err_dets; % [0 5 10 15 20 30 40]; %; % err_mnts = err_dets; % 
err_trans_s = err_dets;

min_p = 0.25;
max_p = 1.05;
min_d = 0;
max_d = 500;

err_drfts = [0 20];

cmaps = colormap(brewermap(7,'Set1')); % 6,'Dark2'));
% point to coulours for current "solve" orien % tn prog
% map_id = [5 5 7 2 3]; % [1 1 3 6 4];
clr_gm = cmaps(5,:);
clr_gg = cmaps(4,:);
clr_tc = cmaps(3,:);


% for k_dr = 1:numel(err_drfts)
%     
%     err_d = err_drfts(k_dr);
    
    for i_det = 1:numel(err_dets)
        
        err_det = err_dets(i_det);
        err_trans = err_trans_s(i_det);
         err_mnt = err_mnts(i_det); % j_mnt);
         
        for j_mnt = 1:1
            
           
            
            % determine ratios for drift-free migration
            
            % for no wind no autocorr needed
            dir_ij = [num2str(floor(err_det)) '_' num2str(floor(err_trans))  ...
                '_' num2str(err_mnt) ...
                '_0_deg_errs/' species '/'];
            
            GG_tr = load([dir_ij gg_str f_str_1_no_drf 'celest' tr_str f_str_2]);
            GG_same = load([dir_ij gg_str f_str_1_no_drf 'celest' sm_str f_str_2]);
            
            GM_tr = load([dir_ij gm_str f_str_1_no_drf 'celest' tr_str f_str_2]);
            GM_same = load([dir_ij gm_str f_str_1_no_drf 'mag' sm_str f_str_2]);

            
%             TC_mag = load([dir_ij tc_str f_str_1 'celest' tr_str f_str_2]);
            TC_star = load([dir_ij tc_str f_str_1_no_drf 'celest' tr_str f_str_3]);
            
            p_gg_star(i_det,j_mnt,1) = GG_tr.p_within_goal;
            d_gg_star(i_det,j_mnt,1) = GG_tr.med_d_close_all;
            p_gg_mag(i_det,j_mnt,1) = GG_same.p_within_goal;
            d_gg_mag(i_det,j_mnt,1) = GG_same.med_d_close_all;
            
            p_gm_star(i_det,j_mnt,1) = GM_tr.p_within_goal;
%             d_gm_star(i_det,j_mnt,1) = GM_tr.med_d_close_all;
            p_gm_mag(i_det,j_mnt,1) = GM_same.p_within_goal;
            d_gm_mag(i_det,j_mnt,1) = GM_same.med_d_close_all;
            
%             p_tc_mag(i_det,j_mnt,1) = TC_mag.p_within_goal;
%             d_tc_mag(i_det,j_mnt,1) = TC_mag.med_d_close_all;
            p_tc_star(i_det,j_mnt,1) = TC_star.p_within_goal;
            d_tc_star(i_det,j_mnt,1) = TC_star.med_d_close_all;
            
            rp_ggs(i_det,j_mnt,1) = GG_same.p_within_goal/GG_tr.p_within_goal;
            rp_tc_gg(i_det,j_mnt,1) = TC_star.p_within_goal/GG_tr.p_within_goal;
            
            rd_ggs(i_det,j_mnt,1) = GG_same.med_d_close_all/GG_tr.med_d_close_all;        
            rd_tc_gg(i_det,j_mnt,1) = TC_star.med_d_close_all/GG_tr.med_d_close_all;
            
%             rp_gms(i_det,j_mnt,1) = GM_same.p_within_goal/GM_tr.p_within_goal;
%             rp_tc_gm(i_det,j_mnt,1) = TC_star.p_within_goal/GM_tr.p_within_goal;
            
%             rd_gms(i_det,j_mnt,1) = GM_same.med_d_close_all/GM_tr.med_d_close_all;        
%             rd_tc_gm(i_det,j_mnt,1) = TC_star.med_d_close_all/GM_tr.med_d_close_all;

            
            % now repeat for 20 degree drift cases
            
            dir_ij = [num2str(floor(err_det)) '_' num2str(floor(err_trans)) ...
                '_' num2str(err_mnt) ...
                '_' num2str(err_drfts(2)) '_deg_errs/' species '/'];
            
            GG_tr = load([dir_ij gg_str f_str_1_drf 'celest' tr_str f_str_2]);
            GG_same = load([dir_ij gg_str f_str_1_drf 'celest' sm_str f_str_2]);
            
            GM_tr = load([dir_ij gm_str f_str_1_drf 'celest' tr_str f_str_2]);
            GM_same = load([dir_ij gm_str f_str_1_drf 'mag' sm_str f_str_2]);

            
%             TC_mag = load([dir_ij tc_str f_str_1 'celest' tr_str f_str_2]);
            TC_star = load([dir_ij tc_str f_str_1_drf 'celest' tr_str f_str_3]);
            
            p_gg_star(i_det,j_mnt,2) = GG_tr.p_within_goal;
            d_gg_star(i_det,j_mnt,2) = GG_tr.med_d_close_all;
            p_gg_mag(i_det,j_mnt,2) = GG_same.p_within_goal;
            d_gg_mag(i_det,j_mnt,2) = GG_same.med_d_close_all;
            
            p_gm_star(i_det,j_mnt,2) = GM_tr.p_within_goal;
%             d_gm_star(i_det,j_mnt,2) = GM_tr.med_d_close_all;
            p_gm_mag(i_det,j_mnt,2) = GM_same.p_within_goal;
            d_gm_mag(i_det,j_mnt,2) = GM_same.med_d_close_all;
            
%             p_tc_mag(i_det,j_mnt,2) = TC_mag.p_within_goal;
%             d_tc_mag(i_det,j_mnt,2) = TC_mag.med_d_close_all;
            p_tc_star(i_det,j_mnt,2) = TC_star.p_within_goal;
            d_tc_star(i_det,j_mnt,2) = TC_star.med_d_close_all;
            
            rp_ggs(i_det,j_mnt,2) = GG_same.p_within_goal/GG_tr.p_within_goal;
            rp_tc_gg(i_det,j_mnt,2) = TC_star.p_within_goal/GG_tr.p_within_goal;
            
            rd_ggs(i_det,j_mnt,2) = GG_same.med_d_close_all/GG_tr.med_d_close_all;        
            rd_tc_gg(i_det,j_mnt,2) = TC_star.med_d_close_all/GG_tr.med_d_close_all;
            
            % first one last - simulations ignoring maintenance and drift errors
            
            
        end
        
    end

LW = 1.25;
FS = 10.5; % 9;
FS_txt = 10.5; % 9; % 7.5;

H = figure('Position',[200 200 185 185]); %[200 200 350 165]); %
% subplot(1,3,2)
hold
% plot(err_dets,p_gm_star(:,1,1),'Color',clr_gm,'LineStyle','--','LineWidth',LW*1.15)
plot(err_dets,p_gm_mag(:,1,1),'Color',clr_gm,'LineWidth',LW)

plot(err_dets,p_gg_star(:,1,1),'Color',clr_gg,'LineStyle','--','LineWidth',LW*1.15)
plot(err_dets,p_gg_mag(:,1,1),'Color',clr_gg,'LineWidth',LW)
plot(err_dets,p_tc_star(:,1,1),'Color',clr_tc,'LineStyle','--','LineWidth',LW*1.15)
% set(gca,'YTick',0:100:500)
set(gca,'YTick',0.3:0.2:1)
set(gca,'XTick',0:10:40)
set(gca,'FontSize',FS)

    scatter(err_dets(3)-1,p_tc_star(3,1,1),100,'w','fill') % d
text(err_dets(3)-3,p_tc_star(3,1,1)+.025,'(e)','Color',clr_tc,'LineWidth',1.15, 'FontSize', FS_txt) % p_gg_star(3,1,1)+0.05
% if strcmp(species,'Kirtlands Warbler')
%     scatter(10,p_gg_mag(3,1,1)+0.005,150,'w','fill') % d +20 +0.0375
%     text(10-2.,p_gg_mag(3,1,1)+0.005,'(a)','Color',clr_gg,'LineWidth',1.3) % d +20 +0.0375
% else
      scatter(err_dets(3),p_gg_mag(3,1,1)-.05,100,'w','fill') % d
    text(err_dets(3)-3,p_gg_mag(3,1,1)-.05,'(d)','Color',clr_gg,'LineWidth',1.3, 'FontSize', FS_txt) % d  d_gg_mag(3,1,1)+0.05-25
    
     scatter(err_dets(3),p_gg_star(3,1,1),100,'w','fill') 
    text(err_dets(3)-3,p_gg_star(3,1,1)+.025,'(e)','Color',clr_gg,'LineWidth',1.3, 'FontSize', FS_txt) % 
%        scatter(err_dets(3),p_gg_star(3,1,2),100,'w','fill') % d
%     text(err_dets(3)-3,p_gg_star(3,1,2)+.03,'(e)','Color',clr_gg,'LineWidth',1.3, 'FontSize', FS_txt) % p_gg_mag(3,1,2)+0.035+15+10.
 
%     scatter(err_dets(3),p_gm_mag(3,1,1),100,'w','fill') % d
%     text(err_dets(3)-3,p_gm_mag(3,1,1)+.0225,'(c)','Color',clr_gm,'LineWidth',1.3, 'FontSize', FS_txt) % d  d_gm_mag(3,1,1)+0.05-25

    
% end
    
    set(gca,'YAxisLocation','Right')

axis([0 err_dets(end) min_p max_p]) % min_d max_d])  %   
set(gca,'YTick',0.25:0.25:1)
% title('(d)','FontSize',FS) % '(b)'
xlabel('Compass precision (^o)','FontSize',FS)
ylabel('Performance','FontSize',FS)
% ylabel('Arrival probability','FontSize',FS)
% ylabel('Closest approach (km)','FontSize',FS)
% ylabel('Fraction within goal','FontSize',FS)
% set(gca, 'YAxisLocation', 'right')    

H = figure('Position',[350 350 185 185]); %[200 200 350 165]); %
% subplot(1,3,3)
hold
plot(err_dets,p_gg_star(:,1,2),'Color',clr_gg,'LineStyle','--','LineWidth',LW*1.15)
plot(err_dets,p_gg_mag(:,1,2),'Color',clr_gg,'LineWidth',LW)
% plot(err_dets,p_gm_star(:,1,2),'Color',clr_gm,'LineStyle','--','LineWidth',LW*1.15)
plot(err_dets,p_gm_mag(:,1,2),'Color',clr_gm,'LineWidth',LW)

plot(err_dets,p_tc_star(:,1,2),'Color',clr_tc,'LineStyle','--','LineWidth',LW*1.15)
set(gca,'FontSize',FS)
% plot(err_dets,p_tc_mag(:,1,2),'--r','LineWidth',LW)
 scatter(err_dets(3),p_tc_star(3,1,2),120,'w','fill') % d +
text(err_dets(3)-3,p_tc_star(3,1,2)-.005,'(f)','Color',clr_tc,'LineWidth',1.3, 'FontSize', FS_txt) % p_tc_star(3,1,2)+0.05-32
%  scatter(err_dets(2),p_gg_mag(3,1,2)-.01,120,'w','fill') % d +20 +0.0375

%  scatter(err_dets(3),p_gm_mag(2,1,2),120,'w','fill') % d +20 +0.0375
% text(err_dets(3)-3,p_gm_mag(3,1,2)+.05,'(h)','Color',clr_gm,'LineWidth',1.3, 'FontSize', FS_txt) % p_gm_mag(3,1,2)+0.035+15+10.


  scatter(err_dets(3),p_gg_mag(3,1,2),80,'w','fill') % d +
 text(err_dets(3)-3,p_gg_mag(3,1,2),'(f)','Color',clr_gg,'LineWidth',1.3, 'FontSize', FS_txt) % p_gg_mag(3,1,2)+0.035+15+10.

%   scatter(err_dets(3),p_gg_star(3,1,2),120,'w','fill') % d +
% text(err_dets(3)-3,p_gg_star(3,1,2)-.02,'(g)','Color',clr_gg,'LineWidth',1.3, 'FontSize', FS_txt) % p_gg_mag(3,1,2)+0.035+15+10.
% set(gca,'YTick',0:100:500)

set(gca,'YTick',0.25:0.25:1)
set(gca,'XTick',0:10:40)

% add = gca;
% extra = axes('Position',get(add,'Position'),'Color','none', ...
%     'XTick',0:10:40,'XTickLabel',{' ',' ',' ',' ',' '}, ...
%     'YAxisLocation','right','YTick',100:100:350);
% linkaxes([add extra],'xy');
% axis([0 20 120 360])     
axis([0 err_dets(end) min_p max_p]) %  min_d max_d]) %    
set(gca,'YAxisLocation','Left')
% title('(e)','FontSize',FS) % '(c)'
xlabel('Compass precision (^o)','FontSize',FS)   
ylabel('Performance','FontSize',FS)
% ylabel('Arrival probability','FontSize',FS)
% ylabel('Fraction within goal','FontSize',FS)
% ylabel('Closest approach (km)','FontSize',FS)