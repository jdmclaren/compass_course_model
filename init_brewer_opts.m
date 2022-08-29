addpath('brewer')
cmaps = colormap(brewermap(9,'Set1')); % 6,'Dark2'));

transp_val = 0.8;

% point to coulours for current "solve" orien % tn prog
% map_id = [5 5 7 2 3]; % [1 1 3 6 4];
clr_gm = cmaps(5,:);
clr_gg = cmaps(4,:); % cmaps(9,:); % cmaps(8,:); %
clr_tc = cmaps(3,:);
clr_gr = cmaps(9,:); % 
% for "actual" modelled probs plot smaller icons
LW_sm = 1.;
LW_sm_st = 0.75;
dia_sz_1 = 67.5;
star_sz_1 = 105; % 160; % 120; %
cir_sz_1 = 50; % 80; % 75;
sq_sz_1 = 85; % 100;

dia_sz_2 = 0.4*dia_sz_1; % 50;
star_sz_2 = 0.5*star_sz_1; % 75;
cir_sz_2 = 0.5*cir_sz_1; % 45;
sq_sz_2 = 0.4*sq_sz_1;

dia_sz_3 = 0.7*dia_sz_1; % 50;
star_sz_3 = 0.7*star_sz_1; % 75;
cir_sz_3 = 0.7*cir_sz_1; % 45;
sq_sz_3 = 0.7*sq_sz_1;