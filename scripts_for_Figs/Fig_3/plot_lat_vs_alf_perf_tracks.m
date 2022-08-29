% plot_magn_tran
addpath('D:\Oldenburg_models\geomagn_orientn_model\brewer')
%   colormap(brewermap([],'YlOrRd')); %'Reds')); %'YlOrRd')); % colormap(flipud(hot(512)));
% cmap_Accs = colormap(brewermap(8,'Accent')); % 'Dark2')); %'Set2')); % 6 ,'Set1')); % 6
% % point to coulours for current "solve" orien % tn prog
% % map_id = [5 5 7 2 3]; % [1 1 3 6 4];
% cmap_Set3 = colormap(brewermap(12,'Set3')); % 'Dark2')); %'Set2')); % 6 ,'Set1')); % 6
% cmap_Set1 = colormap(brewermap(2,'Set1')); % 'Dark2')); %'Set2')); % 6 ,'Set1')); % 6
% cmap_Set2 = colormap(brewermap(8,'Set2')); % 'Dark2')); %'Set2')); % 6 ,'Set1')); % 6
% cmap_Prd = colormap(brewermap(8,'Paired')); % 

% cont_clr = cmap_Prd(2,:); %cmap_Accs(5,:); % cmap_Set1(9,:); % cmap_Set3(9,:); % cmap_Accs(8,:); %  cmap_Set2(8,:); %     
% cont_clr = [0.1216    0.4706    0.95]; % 
% cont_clr = [0.35    0.45    0.95]; % 

n_trx = numel(alfs_deg);
spacing = floor(size(alf_ps,1)/43);

cont_clr = [0.1216    0.4706    0.95]; % 

% for ii = 1:n_trx
%     
%     iii = ii + spacing*(ii-1);
    iii =1;
 
    alfs_i = squeeze(alf_ps(iii,1,:))*180/pi;
    ths_i = squeeze(lat_ps(iii,1,:))*180/pi;
    d_lat_goal = abs(ths_i(:) - (dep_lat_degs - lat_dist));
    max_ii = find(d_lat_goal == min(d_lat_goal),1,'first');
       
    plot(-alfs_i(1:max_ii),ths_i(1:max_ii),'Color',cont_clr, ... %'b',... % cont_clr, ... %'c',... % 
       'LineStyle',':','linewidth',3) % 'LineStyle','-','linewidth',2) %     

%     plot(-GC_azs+180,latGC*180/pi,'Color',cont_clr, ... %'b',... % cont_clr, ... %'c',... % 
%           'LineStyle','-','linewidth',2) %  'LineStyle',':','linewidth',3) % 
% end
    
plot(RL_az*ones(numel(latRL)-1,1)-180,latRL(1:end-1)*180/pi,'Color',cont_clr,'linewidth',2) %  2.
plot(RL_az*ones(numel(latRL)-1,1)-180,latRL(1:end-1)*180/pi,'Color',cont_clr,'linewidth',3,'LineStyle',':') 

% plot(GC_azs-180,latGC*180/pi,'-.c','linewidth',2.5) %  2.-.
% plot(RL_az*ones(size(latRL))-180,latRL*180/pi,'-.w','linewidth',2.5) %  2.
% 
% plot(-GC_azs+180,latGC*180/pi,'-.c','linewidth',2.5) %  2.
% plot(-RL_az*ones(size(latRL))+180,latRL*180/pi,'-.w','linewidth',2.5) %  2.


allvars = whos;

% Identify the variables that ARE NOT graphics handles. This uses a regular
% expression on the class of each variable to check if it's a graphics object
tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
