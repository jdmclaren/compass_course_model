orng = [0.9100    0.4100    0.1700];
olive = [0.2 0.55 0.05]; % [0.3 0.45 0];
% orng = [0.9300    0.4500    0.2400];

addpath('D:\Oldenburg_models\generic_comp_mig_model\brewer')

% if ~exist('contour_opt')
    
    contour_opt =  false; % true; %  
    
% end

if ~exist('spec_lon_offset')
    spec_lon_offset = 20;
end

if ~exist('small_mrkrs')
    
    small_mrkrs = false;
    
end

if ~exist('plot_lines_opt')
    plot_lines_opt = false;
end

% if ~exist('plot_ms_opt')
    
    plot_ms_opt =  true; %  false; %
    
       plot_var_disp = 5; % 3
       
% end

mrkr_fact = 0.85*small_mrkrs+3*~small_mrkrs; %0.75*

% color schemes Parula (1) Orange Red Yellow (2) and blue purlple (3)
if ~exist('color_opt')
    color_opt =  3; %  2; %    1; % 
end

MkrSz_base = mrkr_fact*12.55; % 45; %  
min_MkrSz =  mrkr_fact*3.5; % 20; %
max_MkrSz =  mrkr_fact*9; % 7.5; % 50; % 

% size star for departure locn
sz_start = 100;

get_spec_plot_params
map_Lat_range = max_lat-min_lat;
     

less_trx = false; % true; % strcmp(species,'Siberian Willow Warbler') ...
%     || strcmp(species,'Marsh Warbler') ||   ...
%     strcmp(species,'GreyCheekThrush');


% if plot_var_disp ~= 8
   
% else
%     n_lines = 0;
% end

% for determining land Bluefin Tuna
cst_accuracy = 95;

% quantiles to plot for visibility of ranges (if plot_var ~=5, i.e. entire
% date range) 
% median extremes at end to be a bit 'covered up' by middling values
if plot_var_test == 3 || plot_var_test == 4
    quants_plt = [0.25 0.75 0.05 0.95 0.5]; % [0.1 0.25 0.4 0.45 0.55 0.6 0.75 0.9 0.05 0.95 0.5]; % 0.01  0.99 0.5]; % 0.05 0.95 0.5]; % 0.5]; %    0.05 0.95 0 1];
    n_qnts = numel(quants_plt);
else
    quants_plt = [];
    n_qnts = 0;
end


if ~plot_ms_opt
    
    MkrSz = 50;

    n_trax = min(round((~less_trx*100 + less_trx*30)*sqrt(map_Lat_range/45)),25); %30; %
    n_lines = ~strcmp(species,'Iql Wheatear')*n_trax;
    transp_val = 1;
    MarkerEdgeAlpha = 0.25;    
    mult_fact = 1;
     
elseif ~contour_opt
    
    switch species
        
        case 'Iql Wheatear'
            
            MkrSz = 15;
            n_trax = 300;
            transp_val = 1;
            n_lines = 0;
            mult_fact = n_hs;
            MarkerEdgeAlpha = 0.25;
           
        case 'Kirtlands Warbler'
            
            MkrSz = 7.5;
            n_trax = 25;
            transp_val =  1;
            n_lines = n_trax + n_qnts;
            trnsp_lin = 0.15;
%             plot_lines_opt = true;
            mult_fact = n_hs;
            MarkerEdgeAlpha = 0.15;
%             map_proj = 'Mercator';
            plot_var_disp = 6; % 3;
            % not all runs used same goal rad! (dist
%             goal_rad = 300;
            
        otherwise
            
            MkrSz = 7.5; % 20; %
            n_trax = 25; % 5; % 
            transp_val =  1; % 0.3;
            n_lines = n_trax + n_qnts;
            trnsp_lin = 1; % 0.5;
            mult_fact = 4; % 1; % 
            MarkerEdgeAlpha =  0.15; % 0.5; % 0.75;
          
    end
            
%     MkrSz = ~strcmp(species,'Iql Wheatear')*10 + strcmp(species,'Iql Wheatear')*15;
% 
%     n_trax = ~strcmp(species,'Iql Wheatear')*45  + strcmp(species,'Iql Wheatear')*300; %  30; % 150; % 25; %
%      n_lines = ~strcmp(species,'Iql Wheatear')*n_trax;
%     transp_val =  ~strcmp(species,'Iql Wheatear')*1 + strcmp(species,'Iql Wheatear')*0.3; % 0.2; %   
% elseif strcmp(species,'Kirtlands Warbler')
% %     plot_lines_opt = true;
% %     n_lines = n_qnts;
%     n_trax = 100; % 70; %  30; % 150; % 25; %
%     transp_val = (or_prog_i~=5)*0.25 + (or_prog_i==5)*0.25; % 0.05; % 0.2; % 
%     MkrSz = 2.; % 1.5; % 3.;
%     plot_cbr_opt = 0;
%     mult_fact = max(round(4*80/(dep_lat_degs-lat_goal_centr)),8);
% 
%     cmaps = colormap(brewermap(7,'Set1')); % 6,'Dark2'));
%     % point to coulours for current "solve" orien % tn prog
%     map_id = [5 5 7 2 3]; % [1 1 3 6 4];
%     cont_clr = cmaps(map_id(or_prog_i),:);
    
else
    
    n_trax = (or_prog_i~=3 || ~strcmp(species,'Sib Will Warb S Hem'))*30 + (or_prog_i==3 && strcmp(species,'Sib Will Warb S Hem'))*50; %60; %  210; %  100; % 30; % 150; % 25; %
    transp_val =1; % 0.75; % (or_prog_i~=3)*0.8 + (or_prog_i==3)*0.6; % 0.05; % 0.2; % 
     if strcmp(map_proj, 'stereo')
        MkrSz = min(cosd(map_lat_cntr)*(65/FLatLims(2))^0.25,0.65); % min(2*cosd(map_lat_cntr),1); %/1.5; % 3.;
     else
        MkrSz =0.75*cosd(map_lat_cntr)*(max_lat-min_lat)/180; % min(2*cosd(map_lat_cntr),1); %/1.5; % 3.;
     end
    plot_cbr_opt = 0;
    mult_fact = max(round(8*80/(dep_lat_degs-lat_goal_centr)),8);

    cmaps = colormap(brewermap(7,'Set1')); % 6,'Dark2'));
    % point to coulours for current "solve" orien % tn prog
    map_id = [5 5 7 2 3 3 3]; % [1 1 3 6 4];
    cont_clr = cmaps(map_id(or_prog_i),:);    
        
end

% addpath('brewer')

% select appropriate colorbar field to plot
% plot_var_disp = plot_var_disps(iprog);

gry = 0.935; % 0.875;
clr_gry = [gry gry gry]; %
clr_gr_gry = [gry-.01 gry gry-.01]; %

med_gry = 0.91;
clr_med_gry = [med_gry med_gry med_gry]; % [med_gry-.05 med_gry med_gry-.05];

dk_gry = 0.5; % 0.575; % 0.625; % 0.85; % 0.725; % 
clr_dk_gry = [dk_gry dk_gry dk_gry];

clr_gr_dk_gry = [dk_gry-.025 dk_gry dk_gry-.025]; %  [gry+.01 gry-.025 gry-.01]; %


% if iprog == 1
map_sz = ~strcmp(species,'Iql Wheatear')*190 + strcmp(species,'Iql Wheatear')*265; % 200 % 300 %
left = map_sz;
bottom = map_sz;
width = map_sz; % 200;
height = map_sz; % 220;

    
    if color_opt == 2
        H = figure('Position',[left bottom width height],'Color','w'); % ,clr_gry); %[0.925 0.975 1]) % [.95 .9 .8])
    elseif color_opt == 1
        H = figure('Position',[left bottom width height],'Color',[0.975 0.995 0.965]); % blueish [0.875 0.935 0.85])
    else
        H = figure('Position',[left bottom width height],'Color','w'); %clr_gry); % [0.925 0.975 1]) % [.95 .9 .8])
    end

    if strcmp(map_proj, 'stereo')

           ax = axesm ('stereo', 'Frame', 'on', 'Grid', 'off','Origin',[map_lat_cntr map_lon_cntr 0],'FlineWidth',1);
%            ax.PositionConstraint = 'innerposition';
           axis('off')
           setm(gca,'FLatLimit',FLatLims) % ,'FLonLimit',[-30 90])
          
%            setm(gca,'FLatLimit',[35 62.5],'FLonLimit',[60 162])
%                    setm(gca,'FLatLimit',[0 35]) % ,[-18 21.5],'FLonLimit',[-21.5 25])
%             setm(gca,'FLatLimit',[25 50],'FLonLimit',[-90 -80])

            % turn off lat lon clipping (needed for Mercator)
            min_lat = min(-5,min_lat);
            max_lat = 90; %77.5;
            min_lon = -180;
            max_lon = 180;
%             min_lon = -180; % -20; % 2; % + (magn_comp_nr == 4)*40;
%             max_lon = 180; % 25
            
    else     
        
               ax =   worldmap([min_lat max_lat],[min_lon max_lon]); % worldmap([40 85],[-100 20]); %
                setm(ax,'mapprojection',map_proj); % ,'MapLatLimit',[-10 90], ...

    end
    
     tightmap
     landareas = shaperead('landareas.shp','UseGeoCoords',true);
     LW_cst = 0.3;
      grn_cst = [0.85 0.945 0.81]; %      
    if color_opt  ~= 2 && color_opt  ~= 4 
%         grn_cst = [0.75 0.8 0.7]; % [0.85 0.965 0.8];
%         geoshow(landareas, 'FaceColor',grn_cst,'EdgeColor',clr_dk_gry) %
        dkst_gry = 0.4; % 375;
        clr_dkst_gry = [dkst_gry dkst_gry dkst_gry];
        geoshow(landareas, 'FaceColor',clr_med_gry,'EdgeColor',clr_dk_gry,'LineWidth',LW_cst) %
%         geoshow(landareas, 'FaceColor',clr_gr_gry,'EdgeColor',clr_gr_dk_gry ) %[0.65 0.75 0.65]) % [0.5 0.7 0.5]) %  [0.0588    0.5020    0.0706]) % [0 0.5 0.15])  %
    elseif color_opt ~= 2
        geoshow(landareas, 'FaceColor',clr_gry,'EdgeColor',clr_dk_gry,'LineWidth',LW_cst) %[0.65 0.75 0.65]) % [0.5 0.7 0.5]) %  [0.0588    0.5020    0.0706]) % [0 0.5 0.15])  %    
%     elseif strcmp(species,'Iql Wheatear')
% %         grn_cst = [0.85 0.945 0.81]; % 0.915  % [0.85 0.965 0.8];
%          % [gy_val gy_val gy_val]) %[0.5 0.7 0.5]) %  [0.0588    0.5020    0.0706]) % [0 0.5 0.15])  % 
%       geoshow(landAreas,'FaceColor',[1 1 .5],'EdgeColor',[.6 .6 .6],'DisplayType','polygon');
      
    else
               
        geoshow(landareas, 'FaceColor',grn_cst,'EdgeColor',clr_dk_gry,'LineWidth',LW_cst)
        
    end

    mlabel('off'); plabel('off'); gridm('off')

xtra_trax = or_prog == 3 && (strcmp(species,'Siberian Willow Warbler') || strcmp(species,'Wheatear'));


if plot_var_disp == 1 % offset inherited from mean inherited direction
    
        var_plot = (alf0_s*180/pi - alf_init);    

elseif plot_var_disp == 2
    
    var_plot = all_sun_azs;
            
elseif plot_var_disp == 3
    
    if mod_alphs

        all_alphs_pl = mod(180+all_alphs,360); % (:,1)-all_errs

    elseif ~strcmp(species,'SEPacific Humpback')

        all_alphs_pl = 180+all_alphs;

    else 

        all_alphs_pl = all_alphs;

    end
    
    var_plot = all_alphs_pl;

elseif plot_var_disp == 4
        
    var_plot = all_doys;
    var_plot(var_plot < 120) = all_doys(var_plot < 120) +365;
    
elseif plot_var_disp == 5
    
    var_plot = doys_0;
    
elseif plot_var_disp == 6
    
    var_plot = 180 + alf0_s*180/pi;
    
elseif plot_var_disp == 7
    
        var_plot = all_errs;
        
elseif plot_var_disp == 8
    
    % will set to arrival date / time
    var_plot = NaN*one_vec;
    
end

half_lat = (max_lat+min_lat)/2;
% MkrSz = max(min(MkrSz_base*sqrt(45/map_Lat_range)* ... %  45/map_Lat_range 
%     (cosd(45)/cosd(half_lat)),max_MkrSz),min_MkrSz); % *sqrt(n_hs*Va_mps/80);
% if iprog ==1
%     Mkr = 'o';
% else
%     Mkr = 'd';
% end

skip_recs = ~strcmp(species,'Iql Wheatear')*(1+(max_lat-min_lat>60)+2*(n_hs*Va_mps  < 40)) + strcmp(species,'Iql Wheatear')*1;

quant_vars = quantile(var_plot(:),quants_plt);

if ~isnan(quant_vars(1))
    for iq = 1:n_qnts
        abs_difq = abs(var_plot - quant_vars(iq));
        quant_recs(iq) = find(abs_difq == min(abs_difq),1,'first');
    end
end

% 
% n_trax = 0;

edge_dist = map_Lat_range/50; % max(map_Lat_range/25,1.5);

% mult_fact is multiplication factor for number of points for 
% for more continuous contours


if strcmp(species,'Iql Wheatear')
    
     max_hs = Inf; % 96;
     
else
    
    max_hs = Inf;
    
end

% can 'fill in' the plot field (date, initial heading easily
% Lat Lons done within loop per traj
n_recs = size(lon_es,2);
% var_plot_orig = var_plot(1:n_trax,:);
var_plot_pl = NaN*ones(n_inds,mult_fact*n_recs);
if plot_var_disp == 2 || plot_var_disp == 3 || plot_var_disp == 4 || plot_var_disp == 7
    
    for jj = 1:n_recs

        var_plot_pl(:,(jj-1)*mult_fact + (1:mult_fact)) = ...
            repmat(var_plot(:,jj),[1 mult_fact]);    
        
    end
    
elseif plot_var_disp ~= 8 % determine arrival as closest approach below
    
    for ii = 1:n_inds
        
        var_plot_pl(ii,:) = var_plot(ii);
    
    end
    
end

i_cl_dist_pl = i_cl_dist_arr*mult_fact;
i_cl_lat_pl = i_cl_lat_arr*mult_fact;

if color_opt == 2
    colormap(brewermap([],'RdBu')) %,'*YlOrRd')) %,'RdBu')) %,'*YlOrRd')) %'RdYlBu')) %,'PiYG')) % 'YlGnBu')) % 
% elseif color_opt == 4
%     colormap(brewermap([],'*YlOrRd')) %,'*YlOrRd')) 
else % if color_opt == 3
    if ~ismember(plot_var_disp,[4 5])
         colormap(brewermap([],'*RdYlBu')) %,'RdBu')) %,'RdYlBu')) %,'PRGn')) % ,'RdYlBu')) %'PiYG')) %'RdBu')) %  'YlGnBu')) % 
    else % dates start red for hotter?
           colormap(brewermap([],'RdYlBu')) %,'R
    end
    
end

% now can fill in quants for plot_var_disp == 8 (arrival date / time)
% if isnan(quant_vars(1))
%     quant_vars = quantile(var_plot(:),quants_plt);
%     for iq = 1:n_qnts
%         abs_difq = abs(var_plot - quant_vars(iq));
%         quant_recs(iq) = find(abs_difq == min(abs_difq),1,'first');
%     end
% end

 i_plt_rec = -1;
 
for ii = 1:n_trax +n_qnts % :-1:1
    % choose every 2nd migrant (date) for viewability of dates with 
    % plot_var_test = 5
    if ii <= n_trax
        i_plt_rec = i_plt_rec + 2;
    else
        i_plt_rec = quant_recs(ii-n_trax); % n_trax+n_qnts+1-ii);
    end
    
    % reselect of WAY off (avoid magnetoclinic nonsense - it's already
    % clearly inferior without those...)
%     if ~strcmp(species,'Siberian Willow Warbler') && ~strcmp(species,'Wheatear') && ~strcmp(species,'Sib Will Warb S Hem')
%         if median(d_close{ia,i_err}(i_plt_rec)) > 6000 || lat_es(i_plt_rec,end-5) > pi/180
% 
%             too_off = true;
%             while too_off && i_plt_rec < n_inds
%                i_plt_rec = i_plt_rec + 1;
%                too_off = median(d_close{ia,i_err}(i_plt_rec)) > 6000;
%             end
% 
%         end
%     end
       
    if ~strcmp(species,'Iql Wheatear')
          icl_dist = max(i_cl_dist_pl(i_plt_rec),i_cl_lat_pl(i_plt_rec));
    else
        idx_nan = find(isnan(lon_es(i_plt_rec,:)),1,'first');
        if isempty(idx_nan)
            icl_dist = size(lon_es,2)*mult_fact;     
        else
            icl_dist = (idx_nan-1)*mult_fact; 
        end
    end
    
    % avoid map edges
    
   
if dep_lon_sp_degs <= 180
        
        Lon_degs = mod(spec_lon_offset+lon_es(i_plt_rec,:)*180/pi+180,360)-180;
        
    else
        
        Lon_degs = spec_lon_offset+lon_es(i_plt_rec,:)*180/pi;
        
    end
    
    Lat_degs = lat_es(i_plt_rec,:)*180/pi;
    
%         figure
%     plot(Lon_degs,Lat_degs)
        
    % remove recs above 89 degs
    lost_at_pole  = find(lat_es(i_plt_rec,:)>80*pi/180,1,'first');
    if ~isempty(lost_at_pole)
%      Lat_degs = lat_es(i_plt_rec,1:lost_at_pole-1)*180/pi;
     icl_dist = min(icl_dist,lost_at_pole)*mult_fact;
    end
    
    % cut off for Sib Willow Warbler if West of -20
    if strcmp(species,'Sib Will Warb S Hem') || strcmp(species,'Siberian Willow Warbler')
        too_West = find(lon_es(i_plt_rec,:)<-60*pi/180,1,'first');
        if ~isempty(too_West)
            icl_dist = min(icl_dist,too_West)*mult_fact;
        end
    end
    
    try
%     if ~strcmp(species,'Wheatear')
        % fill in the field for more continuous contours
        pathXY = [Lon_degs(1:max(icl_dist/mult_fact,2)); Lat_degs(1:max(icl_dist/mult_fact,2))]';
        finalPathXY = interp2D_even(pathXY,mult_fact);
        Lon_degs = finalPathXY(:,1);
        Lat_degs = finalPathXY(:,2); 
        
        % now transform to -180:180
        Lon_degs = mod(Lon_degs+180,360) - 180;
    catch
        keyboard
    end
        
%     end
    
%     if any(Lon_degs > -50)
%         keyboard
%     end

% boost latitudinal edge distance by Latitude (Mercator projection)
    edge_dist_Lats = edge_dist.*cosd(Lat_degs); % *ones(size(Lat_degs)); % 
    
    off_edge = find((Lon_degs(2:end) < min_lon + edge_dist) | ...
        (Lon_degs(2:end) > max_lon - edge_dist) ...
        | (Lat_degs(2:end) < min_lat + edge_dist_Lats(2:end).*sign(Lat_degs(2:end))) | ...
        (Lat_degs(2:end) > max_lat - edge_dist_Lats(2:end)),1,'first'); % | (Lat_degs(2:end) > max_lat - edge_dist)

    if ~isempty(off_edge) && ~strcmp(species,'Iql Wheatear') && ~strcmp(species,'Alk Wheatear') ...
            && ~strcmp(species,'Sib Will Warb S Hem') && ~strcmp(species,'Kirtlands Warbler')
        icl_dist_pl = min(icl_dist,1+off_edge);

        within_edge_later = find((Lon_degs(icl_dist_pl:icl_dist) > min_lon + ...
            edge_dist) & (Lon_degs(icl_dist_pl:icl_dist) < max_lon - edge_dist) ...
            & (Lat_degs(icl_dist_pl:icl_dist) > min_lat + edge_dist_Lats(icl_dist_pl:icl_dist).*sign(Lat_degs(icl_dist_pl:icl_dist))) ...
            & (Lat_degs(icl_dist_pl:icl_dist) < max_lat - edge_dist_Lats(icl_dist_pl:icl_dist)),1,'first'); % _Lats(icl_dist_pl:icl_dist)
    
        if ~isempty(within_edge_later) 
            
            within_edge_later = icl_dist_pl + within_edge_later -1;
            within_edge_last = find((Lon_degs(within_edge_later:icl_dist) > min_lon + edge_dist) & ...
                (Lon_degs(within_edge_later:icl_dist) < max_lon - edge_dist) ...
            & (Lat_degs(within_edge_later:icl_dist) > min_lat + edge_dist_Lats(within_edge_later:icl_dist)) ...
            & (Lat_degs(within_edge_later:icl_dist) < max_lat - edge_dist_Lats(within_edge_later:icl_dist)),1,'last'); % _Lats(within_edge_later:icl_dist)
%             if  ~isempty(within_edge_last)
                within_edge_last = within_edge_later + within_edge_last -1;
%             end
        
        else
            
%             icl_dist_pl = icl_dist;
            within_edge_last = icl_dist;

        end
        
    else
        
        icl_dist_pl = icl_dist;
        within_edge_last = icl_dist;
        within_edge_later = [];
        off_edge_last = [];
    
    end
    
    
%     if iprog ==1

        if plot_var_disp == 2 || plot_var_disp == 3 || plot_var_disp == 4 ||plot_var_disp == 7

            var_pl_i = var_plot_pl(i_plt_rec,1:end);

        else % departure doy doesn't change

            var_pl_i = var_plot_pl(i_plt_rec)*ones(1,mult_fact*n_recs);

        end
        
%         % convert Lons from model to geographic values
%         Lons_geogr = mod(spec_lon_offset+lon_es(i_plt_rec,1:icl_dist_pl)*180/pi +180,360) - 180;
        
        if strcmp(species,'BluefinTuna') || strcmp(species,'Iql Wheatear')  
            
            over_land = landmask(Lat_degs,Lon_degs,cst_accuracy);
            
            if strcmp(species,'BluefinTuna')
                over_land = over_land & Lon_degs > -125;
            else
                over_land = over_land & (Lon_degs > -12 | Lat_degs < 30);
            end
             if sum(over_land) > 0
                 fst_land = find(over_land,1,'first');
                 if strcmp(species,'Iql Wheatear')
                     icl_dist_pl = fst_land;
%                      icl_dist_pl = 41;
                 else % Tuna
                     icl_dist_pl = min(icl_dist,fst_land-1);           
                 end
%                    var_pl_i(1:icl_dist_pl);
             else 
                 icl_dist_pl =  numel(Lat_degs); % *mult_fact
             end
             
        end
        
        if plot_var_disp == 8
            
              var_pl_i = icl_dist_pl*ones(icl_dist_pl,1); %
              
        end
       
try
    if icl_dist_pl/skip_recs < 4 && ~plot_ms_opt
%         if skip_recs > 1
            skip_recs_pl = 1;
%         end
    else
        skip_recs_pl = skip_recs;
    end
    
%     if any(Lat_degs(1:skip_recs_pl:icl_dist_pl) < min_lat+1)
%         
%         keyboard
%         
%     end


         icl_dist_pl = min(icl_dist_pl,max_hs);
         
       if ~contour_opt

            h = scatterm(Lat_degs(1:skip_recs_pl:icl_dist_pl), ...
                Lon_degs(1:skip_recs_pl:icl_dist_pl), ...
                MkrSz,var_pl_i(1:skip_recs_pl:icl_dist_pl),'fill','o'); % ,'MarkerEdgeColor','k'); % ,'LineWidth',1);
            h.Children.MarkerFaceAlpha = transp_val; % 1; % 
    %         if ~strcmp(species,'Iql Wheatear')
                h.Children.MarkerEdgeColor = 'k'; % clr_dkst_gry; % 
                h.Children.MarkerEdgeAlpha = MarkerEdgeAlpha;
                h.Children.LineWidth = 0.05; %
    %         end

            if ~isempty(within_edge_later)

                if within_edge_last > within_edge_later
                    if numel(within_edge_later:skip_recs:within_edge_last) < 4
                        skip_recs_pl = 1;
                    else
                        skip_recs_pl = skip_recs;
                    end


    %                 if any(Lat_degs(within_edge_later:skip_recs_pl:within_edge_last) < min_lat+1)
    %         
    %                     keyboard
    % 
    %                 end

                    h2 = scatterm(Lat_degs(within_edge_later:skip_recs_pl:within_edge_last), ...
                        Lon_degs(within_edge_later:skip_recs_pl:within_edge_last), ...
                        MkrSz,var_pl_i(within_edge_later:skip_recs_pl:within_edge_last),'fill','o'); % ,'MarkerEdgeColor','k'); % ,'LineWidth',1);
                    h2.Children.MarkerFaceAlpha = transp_val; % 1; % 
                    h2.Children.MarkerEdgeColor = 'k';  
                    h2.Children.MarkerEdgeAlpha = MarkerEdgeAlpha;  
                    h.Children.LineWidth = 0.05; %
                    
                end

            end
            
       else
           
          
               h = scatterm(Lat_degs(1:skip_recs_pl:icl_dist_pl), ...
                Lon_degs(1:skip_recs_pl:icl_dist_pl), ...
                MkrSz,'MarkerEdgeColor',cont_clr,'MarkerFaceColor',cont_clr,'Marker','o'); % 'g' olive orng 'fill','MarkerEdgeColor','k'); % ,'LineWidth',1);
            h.Children.MarkerFaceAlpha = transp_val; % 1; % 
            h.Children.MarkerEdgeAlpha = 0; % transp_val; % 1; % 
    %         if ~strcmp(species,'Iql Wheatear')
%                 h.Children.MarkerEdgeColor = 'k'; % clr_dkst_gry; % 
%                 h.Children.MarkerEdgeAlpha = ~strcmp(species,'Iql Wheatear')*0.75 + strcmp(species,'Iql Wheatear')*0.25; % transp_val; % 0.75; %    
%                 h.Children.LineWidth = 0.05; %
        
       end
        
catch
    keyboard
end
      
        
%     else
%          h = scatterm(rad2deg(lat_es(i_plt_rec,1:icl_dist)), ...
%             spec_lon_offset+rad2deg(lon_es(i_plt_rec,1:icl_dist)), ...
%             MkrSz,all_alphs_pl(i_plt_rec,1:icl_dist),'fill','d'); % ,'LineWidth',1);              
%     end
    
    if plot_lines_opt &&  ii > n_trax + n_qnts - n_lines
%         if color_opt == 1     
% try

             p = plotm(Lat_degs(1:icl_dist_pl),Lon_degs(1:icl_dist_pl), ...
                'Color',clr_dk_gry,'LineStyle','-','LineWidth',1); % '--' 0
              p.Color = [0.3,0.3,0.3,trnsp_lin]; %[orng 1]; % [0 0.5020 0.5020 0.5]; % [0.3,0.3,0.3,0.5];
%             p = plotm(Lat_degs(1:within_edge_last),Lon_degs(1:within_edge_last), ...
%                 'Color',clr_dk_gry,'LineStyle','-','LineWidth',1); % '--' 0.75 1.25 clr_dk_gry 'Color','w','LineStyle',
          
%             p.Color = [dk_gry, dk_gry, dk_gry,max(transp_val,0.8)];
%             p.Color = [0.3,0.3,0.3,transp_val];
% catch
%     keyboard
% end
            %         else
%             plotm(rad2deg(lat_es(i_plt_rec,1:icl_dist)),spec_lon_offset+rad2deg(lon_es(i_plt_rec,1:icl_dist)), ...
%                 '--w','LineWidth',1.25);  
%         end
    end
end

if ~contour_opt
    caxis([min(quant_vars(:)) max(quant_vars(:))])
end
%     if strcmp(species,'Iql Wheatear')
% %         grn_cst = [0.85 0.945 0.81]; % 0.915  % [0.85 0.965 0.8];
%         geoshow(landareas, 'FaceColor',grn_cst,'EdgeColor',clr_dk_gry) %
%     end

load coastlines

plotm(coastlat,coastlon,'Color',clr_dk_gry)
    
% plot arrival area
% if color_opt == 1   
%     scatterm(lat_goal_centr,lon_goal_centr,sz_goal,'ws','fill','MarkerEdgeColor','r','LineWidth',1.5)
%     scatterm(dep_lat_degs,spec_lon_offset+dep_lon_degs,sz_start,'wh','fill','MarkerEdgeColor','r','LineWidth',0.5)
% else

% plot departure area witt white star
   sz_fact_lat = sqrt(45/map_Lat_range);
    scatterm(dep_lat_degs,spec_lon_offset+dep_lon_degs,sz_start*sz_fact_lat,'wo','fill') % ,'MarkerEdgeColor','k','LineWidth',0.75)
   scatterm(dep_lat_degs,spec_lon_offset+dep_lon_degs,sz_start*sz_fact_lat,'kh','fill') % ,'MarkerEdgeColor','k','LineWidth',0.5)
   
% plot goal area circle with withe circle (previously used sz_goal)
%     scatterm(lat_goal_centr,lon_goal_centr,sz_start,'wo','fill','MarkerEdgeColor','k','LineWidth',1.5)    
  if ~strcmp(species,'Iql Wheatear')
    goal_rad_arc = goal_rad/R_Earth_km*180/pi;
   [lat_goal_area,lon_goal_area] = scircle1(lat_goal_centr,lon_goal_centr,goal_rad_arc);% end
    plotm(lat_goal_area,lon_goal_area,'w','LineWidth',2)    %2.5*sz_fact_lat
    plotm(lat_goal_area,lon_goal_area,'k','LineWidth',max(1.*sz_fact_lat,1))   
  end
% title({['closest approach km',  ...
%    ' '  num2str(round(0.1*median(dec_close{1}))*10)],  ....
%    ['E-W bias ' num2str(round(0.1*median(d_cl_sgn))*10) ' km']} )

% if plot_title 
%     
%      title({[or_pr_str{or_prog}], ['closest approach km',  ...
%        ' '  num2str(round(0.1*median(dec_close{1}))*10)],  ...
%        ['p(>500 km) ' num2str(round(100*p_gt_500_km)/100)]} )
% 
% end
% , ...    num2str(round(0.1*mean(dec_close{1}))*10)
if plot_cbr_opt == 2 || (plot_cbr_opt == 1 && iprog == numel(or_progs))

    if plot_cbr_opt == 1
        
        figure('Position',[left bottom width height/1.4])
    
    end
    
    cc = colorbar;
      caxis([min(quant_vars(:)) max(quant_vars(:))]);
    if plot_var_disp == 6
        title(cc,{'inherited', 'heading (^o)'})
%        caxis(180 + [quantile(alf0_s*180/pi,0.) quantile(alf0_s*180/pi,1)])
    elseif plot_var_disp == 4 || plot_var_disp == 5
        set(cc,'YTick',[196 213 227 244 258 274 288 304 319 335 349 366 380 397 411], ...
        'YTickLabel',{'Jul 15','Aug 1','Aug 15','Sep 1','Sep 15','Oct 1','Oct 15', ...
        'Oct 31','Nov 15','Dec 1','Dec 15','Jan 1','Jan 15','Feb 1','Feb 15'})
%         caxis([quantile(doys_0,0) quantile(doys_0,1)])
    elseif plot_var_disp == 1
         title(cc,'offset heading (^o)')
    elseif  plot_var_disp == 2
        title(cc,'sun azimuth (^o)')
    elseif plot_var_disp == 7
        title(cc,'compass error (^o)')
    end
    
    set(gca,'FontSize',9)

    axis('off')
    

    if ismember(plot_var_disp,[4 5])
        if color_opt == 2
             colormap(brewermap([],'RdBu')) %,'*YlOrRd')) %,'RdBu')) %,'*YlOrRd')) %'RdYlBu')) %,'PiYG')) % 'YlGnBu')) % 
        else
            colormap(brewermap([],'RdYlBu'))
        end
    else % if color_opt == 3
        if color_opt == 2
             colormap(brewermap([],'*RdBu')) %,'*YlOrRd')) %,'RdBu')) %,'*YlOrRd')) %'RdYlBu')) %,'PiYG')) % 'YlGnBu')) % 
        else
            colormap(brewermap([],'*RdYlBu'))
        end
    end
  

end
   
if save_figs_opt
    
        dir_nm_Fig = ['output_species_runs/' num2str(err_dets) '_deg_err/Figs/' species '/'];
        [status, msg, msgID] = mkdir(dir_nm_Fig);
        
%     savefig(H,['output_species_runs/' num2str(err_det_base) '_deg_err/Figs/traj_' species '_' ...
%                      or_pr_strs{or_prog_idx(i_prog)}],'compact') 
        savefig(H,[dir_nm_Fig ...
                     or_pr_strs{or_prog_idx(i_prog)} '_ar_' num2str(100*ar_error_det)],'compact');
end

% caxis(clims)

if plot_close_dist_opt
    
    figure
    histogram(dec_close{ia,i_err},0:50:3000) % (d_cl_sgn) % 
    title(['med closest approach ' num2str(median(dec_close{ia,i_err})) ' km'])

%     title(['med closest approach ' num2str(median(d_cl_sgn)) ' km'])

    n_pts_plot = 2500;

    if is_sun_comp
        
        H = figure('Position',[left bottom width*1.33 height]);
        
        within_qnts = dec_close{1} > quantile(dec_close{1},0.01) & ...
             dec_close{1} < quantile(dec_close{1},0.99);
         
         within_qnts(n_pts_plot+1:end) = false;
         
         % determine "inherited" geo heading equivalent for sun comp
         if sun_inher_opt
             
             % all_sun_azs already has 180 added so need to subtract
             % twice!)
             alf_plot_var = 180+alf0_s(within_qnts)*180/pi  - ...
                 all_sun_azs(within_qnts,1); %  +sun_az_0_mn*180/pi;
             
         else % geo/magn inheritance
             
             alf_plot_var = 180+alf0_s(within_qnts)*180/pi;
             
         end
         
         if plot_ms_opt
             
            scatter(alf_plot_var,dec_close{1}(within_qnts), ...
                30,doys_0(within_qnts),'fill', ...
                'MarkerEdgeColor','k','MarkerEdgeAlpha',0.25);
%             'MarkerFaceAlpha',min(transp_val,0.35), ...


         else
             
            scatter(alf_plot_var,dec_close{1}(within_qnts), ...
                60,doys_0(within_qnts),'fill', ...
            'MarkerEdgeColor','k','MarkerFaceAlpha',transp_val,...
            'MarkerEdgeAlpha',max(transp_val,0.5));
         
         end
         
%         title(['med closest approach ' num2str(median(dec_close{ia,i_err})) ' km'])
         ylabel('Closest approach (km)')
         xlabel('Initial heading (^o)')

        %      set(gca, 'YScale', 'log')
         set(gca,'FontSize',9)
        %     xtickangle(45)
        %     ylim([10 max(d_close{1})])
%         if color_opt == 2
%             colormap(brewermap([],'RdBu')) %,'*YlOrRd')) %,'*YlOrRd')) %'RdYlBu')) %,'PiYG')) % 'YlGnBu')) % 
%         else % if color_opt == 3
             colormap(brewermap([],'RdBu')) %,'RdYlBu')) %,'RdYlBu')) %,'PRGn')) % ,'RdYlBu')) %'PiYG')) %'RdBu')) %  'YlGnBu')) % 
%         end
        %         caxis([215 280])
        
%         if plot_cbr_opt == 2 
            
             cc = colorbar;
            title(cc,{'Departure date'},'FontSize',9)
            set(cc,'YTick',[196 213 227 244 258 274 288 304 319 335 349 366 380 397 411], ...
            'YTickLabel',{'Jul 15','Aug 1','Aug 15','Sep 1','Sep 15','Oct 1','Oct 15', ...
            'Oct 31','Nov 15','Dec 1','Dec 15','Jan 1','Jan 15','Feb 1','Feb 15'})            
            if log_cl_dist_opt
                set(gca, 'YScale', 'log') % ColorScale'
                 ylim([quantile((dec_close{1}),0.005) quantile((dec_close{1}),0.995)])
            end
            
%         end
        
        if save_figs_opt
    

%             
%             savefig(H,['output_species_runs/' num2str(err_dets) '_deg_err/Figs/dist_hd_dt_' species '_' ...
%                              or_pr_strs{or_prog_idx(i_prog)}],'compact') 
            savefig(H,[dir_nm_Fig  ...
                             or_pr_strs{or_prog_idx(i_prog)} '_ar_' num2str(100*ar_error_det)],'compact')
                         
        end
        
    else
%         set(gca, 'YScale', 'linear')       
%     end
%     else
%         figure
%         scatter(180+alf0_s*180/pi,dec_close{1},120,all_errs(:,1),'fill', ...
%         'MarkerEdgeColor','k','MarkerFaceAlpha',max(transp_val,0.5),...
%         'MarkerEdgeAlpha',max(transp_val,0.5)); 
%          cc = colorbar;
%         title(cc,{'Initial error (^o)'})
%     if log_cl_dist_opt
%         set(gca, 'YScale', 'log') % ColorScale'
%          ylim([10 max(dec_close{1})])
% %     else
% %         set(gca, 'YScale', 'linear')       
%     end
%      ylabel('Closest approach (km)')
%      xlabel('Initial heading (^o)')
%      
% %      set(gca, 'YScale', 'log')
%      set(gca,'FontSize',14)
% %     xtickangle(45)
% %     ylim([10 max(d_close{1})])
%     colormap(brewermap([],'*RdYlBu')) 

    end

end