addpath('brewer')
cmaps = colormap(brewermap(9,'Set1')); % 6,'Dark2'));
% point to coulours for current "solve" orien % tn prog
% map_id = [5 5 7 2 3]; % [1 1 3 6 4];
clr_gm = cmaps(5,:);
clr_gg = cmaps(4,:); % cmaps(9,:); % cmaps(8,:); %
clr_tc = cmaps(3,:);
clr_gr = cmaps(9,:); % 
clr_fx = cmaps(2,:);
clr_mc = cmaps(7,:); % 

calc_goal_dist_migr_sp_spec
adj_br_N = front_brdth.*sqrt(n_hat_fls); % [0.18 0.26 0.28 0.32 0.34 0.48 0.48 0.54 0.65];

plt_mcl = true; % false;
plt_fix = true;
plot_ggl = true; 

% err opt 1 is 20 deg stepwise (ie all effectively det only)
% otherwise 20/20/20/15 degs (det/tran/mnt/drft)
err_opt = 1; % 2; % 

if err_opt == 1
    p_GGL = [83 68 49 39 96 59 55 79 70]; % [88 78 60 38 99 67 38 81 67]; % 
    p_GML = [86 60 55 34 97 62 60 63 59]; % [91 72 67 31 99 70 31 67 55]; % 
    p_MCL = [77 57 55 29 96 8 58 29 1]; % [87 73 67 34 85 8 34 39 0.4]; % 
    p_FIX = [82 67 50 38 95 63 57 83 63]; % [66 52 38 38 82 49 37 57 39]; % 
    p_TCSC = [91 81 68 55 99 76 81 98 94]; % [71 58 45 42 91 54 42 80 69]; % 
else
    p_GGL = [88 78 60 38 99 67 38 81 67]; % 
    p_GML = [91 72 67 31 99 70 31 67 55]; % 
    p_MCL = [87 73 67 34 85 8 34 39 0.4]; % 
    p_FIX = [66 52 38 38 82 49 37 57 39]; % 
    p_TCSC = [71 58 45 42 91 54 42 80 69]; % 

end

idx = 1:9; %  [1 4 6 8 9]; % 
idx_suppl = [2 4 6 8 9]; % [1 3 5 7]; % 
idx_main = [1 3 5 7]; % [2 4 6 8 9]; % 
LW = 1.25; 

offsts = [-0.0075 0.0075 0 0  0]; % [-0.0075:.00375:0.00375 0];

H = figure('Position',[300 100 240 575]); % ,clr_gry); %[0.925 0.975 1]) % [.95 .9 .8])
% H = figure('Position',[300 100 300 120 ]); %[300 100 400 220 ]); %

hold

scatter(offsts(1)+adj_br_N(idx_main),p_GGL(idx_main)/100,180,'Marker','h', ...
    'MarkerFaceColor',clr_gg,'MarkerEdgeColor','w'); % ,'fill') % ,'LineWidth',LW)
scatter(offsts(2)+adj_br_N(idx_main),p_GML(idx_main)/100,120,'Marker','d', ...
    'MarkerFaceColor',clr_gm,'MarkerEdgeColor','w')
scatter(offsts(3)+adj_br_N(idx_main),p_MCL(idx_main)/120,110,'Marker','>', ...
    'MarkerFaceColor',clr_mc,'MarkerEdgeColor','w')
scatter(offsts(4)+adj_br_N(idx_main),p_FIX(idx_main)/100,120,'Marker','s', ...
    'MarkerFaceColor',clr_fx,'MarkerEdgeColor','w') % ,'LineWidth',LW)
scatter(offsts(5)+adj_br_N(idx_main),p_TCSC(idx_main)/100,95,'Marker','o', ...
    'MarkerFaceColor',clr_tc,'MarkerEdgeColor','w')

scatter(offsts(1)+adj_br_N(idx_suppl),p_GGL(idx_suppl)/100,75, ...
    clr_gg,'h','LineWidth',LW); % ,'fill') % ,'LineWidth',LW)
scatter(offsts(2)+adj_br_N(idx_suppl),p_GML(idx_suppl)/100,60, ...
    clr_gm,'d','LineWidth',LW)
scatter(offsts(3)+adj_br_N(idx_suppl),p_MCL(idx_suppl)/120,50, ...
    clr_mc,'>','LineWidth',LW)
scatter(offsts(4)+adj_br_N(idx_suppl),p_FIX(idx_suppl)/100,65, ...
    clr_fx,'s','LineWidth',LW)
scatter(offsts(5)+adj_br_N(idx_suppl),p_TCSC(idx_suppl)/100,55, ...
    clr_tc,'o','LineWidth',LW)    

fb_adjs = min(adj_br_N):0.01:max(adj_br_N);

kap_br = 1./(20*pi/180).^2;

bes_0_br = besseli(0,kap_br,1);
sq_bes_1_br = sqrt(besseli(1,kap_br,1)/bes_0_br);

plot(fb_adjs,erf(sq_bes_1_br*fb_adjs.*sqrt(0.5)./(20*pi/180)),'--k','LineWidth',1);

% legend({'Geogr Lox','Geomg Lox','Magnetocl','Fixed Sun','TC Sun','planar Norm.'})

% xlabel({'Length-adjusted'; 'goal-breadth'},'FontSize',9)
xlabel({'Length-adjusted goal-breadth, \beta_{adj}'},'FontSize',9)
ylabel('Arrival probability','FontSize',9)
% set(gca,'YTick',0:.25:1) % ,'YScale','Linear')
set(gca,'YAxisLocation','Left','YTick',0:.25:1,'XTick',0.1:.1:0.7) % ,'YScale','Linear')
xlim([0.15 0.7])