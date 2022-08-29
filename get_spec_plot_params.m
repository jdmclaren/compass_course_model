switch species

case 'Iql Wheatear'
    
    min_lat = 21.5;
    max_lat = 71; % 72.5;
    min_lon = -70; % -72.5; % -85; %  + (magn_comp_nr == 4)*40;
    max_lon = 10; % 25;
%     map_proj = 'Mercator'; % 'eqdazim'; %
    
    clims = (plot_var_disp==3)*[150 330] + (plot_var_disp==5)*[min(doys_0) max(doys_0)];
    mod_alphs = true;
    sz_goal = 300;
    dep_lat_pl = dep_lat_degs;
     xmax = 5000;
    log_cl_dist_opt = true;
    % parameters for stereo projections
     map_lat_cntr = 55; % 47.5; %58; %  67 ; % 45;
     map_lon_cntr = -27.5; % -8; % 72.5
    FLatLims = [0 25]; %[0 30]; % [0 20]; % 
%     plot_lines_opt = true;
        
case 'Alk Wheatear'

    min_lat = -25;
    max_lat = 75;
    min_lon = 5; % + (magn_comp_nr == 4)*40;
    max_lon = 200; % 25;
%     map_proj = 'Mercator'; % 'eqdazim'; %
    
    clims = (plot_var_disp==3)*[150 330] + (plot_var_disp==5)*[min(doys_0) max(doys_0)];
    mod_alphs = true;
    sz_goal = 300;
    dep_lat_pl = dep_lat_degs;
     xmax = 5000;
    log_cl_dist_opt = true;
    % parameters for stereo projections
     map_lat_cntr = 52; % 48.5;
     map_lon_cntr = 77.5; % 70; % 72.5
    FLatLims = [0 62.5];
        
case 'Monarch'

    min_lat = 17; % 16;
    max_lat = 46; % 48.5;
    min_lon = -107.; % + (magn_comp_nr == 4)*40;
    max_lon = -77.5; % -67.5; %-62.5; %  25;
%     map_proj = 'Mercator'; % 'eqdazim'; %

    clims = (plot_var_disp==3)*[180 260] + (plot_var_disp==5)*[min(doys_0) max(doys_0)];
    mod_alphs = true;
%     sz_goal = 150;
    dep_lat_pl = dep_lat_degs;
     xmax = 1000;
     log_cl_dist_opt = false;
     map_lat_cntr = 31.5; % 33.5;
     map_lon_cntr = -92; % -88; % 72.5
    FLatLims = [0 17.5]; % [0 22];
    
case 'Blackcap'

    min_lat = 15;
    max_lat = 55;
    min_lon = -5; % + (magn_comp_nr == 4)*40;
    max_lon = 25; % 25;
%     map_proj = 'Mercator'; % 'eqdazim'; %
%     spec_lon_offset = -170;
    dep_lat_pl = dep_lat_degs;
    xmax = 1000;
    log_cl_dist_opt = false;
     map_lat_cntr = 30;
     map_lon_cntr = 10; % 72.5
    FLatLims = [0 40];
    
case 'GreyCheekThrush'  

    min_lat = -10;
    max_lat = 70;
    min_lon = -145; % + (magn_comp_nr == 4)*40;
    max_lon = -45; % 25;
%     map_proj = 'Mercator'; % 'eqdazim'; %

    clims = (plot_var_disp==3)*[90 170] + (plot_var_disp==5)*[min(doys_0) max(doys_0)];
    mod_alphs = true;
%     sz_goal = 1000;
    dep_lat_pl = dep_lat_degs;
    xmax = 3000;
    log_cl_dist_opt = true;
     map_lat_cntr = 32.5;
     map_lon_cntr = -90; % 72.5
    FLatLims = [0 46.75];
    
    case 'Whinchat'  

    min_lat = -10;
    max_lat = 70;
    min_lon = 5; % + (magn_comp_nr == 4)*40;
    max_lon = 80; % 25;
%     map_proj = 'Mercator'; % 'eqdazim'; %
%     spec_lon_offset = -105;          
    dep_lat_pl = dep_lat_degs;
     xmax = 1000;
    log_cl_dist_opt = false;
     map_lat_cntr = 34;
     map_lon_cntr = 20; % 72.5
    FLatLims = [0 40];
    
case 'BluefinTuna'

    min_lat = 10;
    max_lat = 60;
    min_lon = 130; % + (magn_comp_nr == 4)*40;
    max_lon = 260; % 25;
    map_proj = 'Mercator'; % 'eqdazim'; %

    clims = (plot_var_disp==3)*[60 120] + (plot_var_disp==5)*[min(doys_0) max(doys_0)];
    mod_alphs = false;
%     sz_goal = 1000;
    dep_lat_pl = dep_lat_degs;
     xmax = 3000;
    log_cl_dist_opt = false;
     map_lat_cntr = 20;
     map_lon_cntr = -160; % 72.5
    FLatLims = [0 50];
    
case 'Siberian Willow Warbler'

    min_lat = -5;
    max_lat = 82.5; % 77.5; % 95; %
    min_lon = 25; % 30; % -20; % 2; % + (magn_comp_nr == 4)*40;
    max_lon = 175; % 25;
%     map_proj = 'ortho'; % 'Mercator'; % 'stereo'; %  'eqdazim'; %  
    map_lat_cntr = 50;
    map_lon_cntr = 55;
    
    clims = (plot_var_disp==3)*[150 330] + (plot_var_disp==5)*[min(doys_0) max(doys_0)];
    mod_alphs = true;
    sz_goal = 300;
    dep_lat_pl = dep_lat_degs;
     xmax = 5000;
    log_cl_dist_opt = true;
     map_lat_cntr = 48.5;
     map_lon_cntr = 65; % 72.5
    FLatLims = [0 57.5];
    
case 'Sib Will Warb S Hem'
    
    min_lat = -30;
    max_lat = 82.5; % 77.5; % 95; %
    min_lon = 25; % 30; % -20; % 2; % + (magn_comp_nr == 4)*40;
    max_lon = 175; % 25;
%     map_proj = 'ortho'; % 'Mercator'; % 'stereo'; %  'eqdazim'; %  
    
    clims = (plot_var_disp==3)*[150 330] + (plot_var_disp==5)*[min(doys_0) max(doys_0)];
    mod_alphs = true;
    sz_goal = 300;
    dep_lat_pl = dep_lat_degs;
     xmax = 5000;
    log_cl_dist_opt = true;
     map_lat_cntr = 40;
     map_lon_cntr = 60; % 72.5
    FLatLims = [0 67.5];
            
case 'SEPacific Humpback'

    min_lat = -70;
    max_lat = 10;
    min_lon = -120; % + (magn_comp_nr == 4)*40;
    max_lon = -50; % 25;
%     map_proj = 'Mercator'; % 'eqdazim'; %'gnomic'; % 'stereo'; %  'Mercator'; % 
    % flip lats for plotting on geogr map
    lat_es = -lat_es;
    % flip dep lat
    dep_lat_pl = -dep_lat_degs;
    clims = (plot_var_disp==3)*[-20 40] + (plot_var_disp==5)*[min(doys_0) max(doys_0)];
    mod_alphs = false;
    sz_goal = 1000;
     xmax = 3000;
    log_cl_dist_opt = false;
    map_lat_cntr = -30; % 35;
    map_lon_cntr = -115; % -88; % -90;
    FLatLims = [0 35];
    
    case 'Kirtlands Warbler'

    min_lat = 21.5; % 22.; %
    max_lat = 46.5; % 47; % 
    min_lon = -86; % -89; % + (magn_comp_nr == 4)*40;
    max_lon = -71.5; % -71; % 25;
%     map_proj = 'Mercator'; % 'eqdazim'; %'gnomic'; % 'stereo'; %  'Mercator'; % 

    clims = (plot_var_disp==3)*[145 190] + (plot_var_disp==5)*[min(doys_0) max(doys_0)];
    mod_alphs = true;
%     sz_goal = 150;
    dep_lat_pl = dep_lat_degs;
    xmax = 1000;
    log_cl_dist_opt = false;

    map_lat_cntr = 34.25; % 35;
    map_lon_cntr =  -78; %-80; % -88; % -90;
    FLatLims = [0 12.5]; % [0 13.5];
    
    case 'Spotted Flycatcher'

    min_lat = -5;
    max_lat = 70;
    min_lon = 0; % + (magn_comp_nr == 4)*40;
    max_lon = 40; % 25;
%     map_proj = 'Mercator'; % 'eqdazim'; %'gnomic'; % 'stereo'; %  'Mercator'; % 

    clims = (plot_var_disp==3)*[145 190] + (plot_var_disp==5)*[min(doys_0) max(doys_0)];
    mod_alphs = true;
%     sz_goal = 350;
    dep_lat_pl = dep_lat_degs;
    xmax = 2000;
    log_cl_dist_opt = false;  
     map_lat_cntr = 15; % 35;
    map_lon_cntr = 5; % -88; % -90;
    FLatLims = [0 55];
    
case 'Finn Marsh Warbler'

    min_lat = -2.5; % -12.5;
    max_lat = 65;
    min_lon = 5; % + (magn_comp_nr == 4)*40;
    max_lon = 52.5; % 55.; % 25;
%     map_proj = 'Mercator'; % 'eqdazim'; %'gnomic'; % 'stereo'; %  'Mercator'; % 

    clims = (plot_var_disp==3)*[145 190] + (plot_var_disp==5)*[min(doys_0) max(doys_0)];
    mod_alphs = true;
%     sz_goal = 350;
    dep_lat_pl = dep_lat_degs;
    xmax = 2000;
    log_cl_dist_opt = false;  
    map_lat_cntr = 32; % 35;
    map_lon_cntr = 30; % -88; % -90;
    FLatLims = [0 34];
    
case 'France Marsh Warbler'

    min_lat = 45; % -2.5; % -12.5;
    max_lat = 70; % 52.5;
    min_lon = 10; % -4; % + (magn_comp_nr == 4)*40;
    max_lon = 50; % 46; % 25;
%     map_proj = 'Mercator'; % 'eqdazim'; %'gnomic'; % 'stereo'; %  'Mercator'; % 

    clims = (plot_var_disp==3)*[145 190] + (plot_var_disp==5)*[min(doys_0) max(doys_0)];
    mod_alphs = true;
%     sz_goal = 350;
    dep_lat_pl = dep_lat_degs;
    xmax = 2000;
    log_cl_dist_opt = false;  

    map_lat_cntr = 60; % 28.5; % 35;
    map_lon_cntr = 30; % 25; % -88; % -90;
    FLatLims = [0 15]; % [0 35];
    
case 'Nathusius'

    min_lat = 35;
    max_lat = 60;
    min_lon = -8; % + (magn_comp_nr == 4)*40;
    max_lon = 25; % 25;
%     map_proj = 'Mercator'; % 'eqdazim'; %'gnomic'; % 'stereo'; %  'Mercator'; % 

    clims = (plot_var_disp==3)*[195 240] + (plot_var_disp==5)*[min(doys_0) max(doys_0)];
    mod_alphs = true;
%     sz_goal = 500;
    dep_lat_pl = dep_lat_degs;
    xmax = 1000;
    log_cl_dist_opt = false;  
    map_lat_cntr = 47.5; % 35;
    map_lon_cntr = 10; % -88; % -90;
    FLatLims = [0 12.5];
    
case 'GreyCheekThrush NF'  

    min_lat = 5;
    max_lat = 55;
    min_lon = -85; % + (magn_comp_nr == 4)*40;
    max_lon = -45; % 25;
%     map_proj = 'Mercator'; % 'eqdazim'; %'gnomic'; % 'stereo'; %  'Mercator'; % 

    clims = (plot_var_disp==3)*[195 240] + (plot_var_disp==5)*[min(doys_0) max(doys_0)];
    mod_alphs = true;
%     sz_goal = 100;
    dep_lat_pl = dep_lat_degs;
    xmax = 2000;
    log_cl_dist_opt = false;         
    map_lat_cntr = 35; % 35;
    map_lon_cntr = -60; % -88; % -90;
    FLatLims = [0 50];
    
case 'Eleonoras Falcon'

    min_lat = -25;
    max_lat = 45; % 77.5; % 95; %
    min_lon = 5; % 30; % -20; % 2; % + (magn_comp_nr == 4)*40;
    max_lon = 55; % 25;
%     map_proj = 'ortho'; % 'Mercator'; % 'stereo'; %  'eqdazim'; %  

    clims = (plot_var_disp==3)*[150 330] + (plot_var_disp==5)*[min(doys_0) max(doys_0)];
    mod_alphs = true;
    sz_goal = 600;
    dep_lat_pl = dep_lat_degs;
     xmax = 5000;
    log_cl_dist_opt = true;
    map_lat_cntr = 0;
    map_lon_cntr = 25;

    FLatLims = [0 50];
    
 case 'Ring Ouzel'  

    min_lat = 30;
    max_lat = 60;
    min_lon = -10; % + (magn_comp_nr == 4)*40;
    max_lon = 5; % 25;
%     map_proj = 'Mercator'; % 'eqdazim'; %

    clims = (plot_var_disp==3)*[90 170] + (plot_var_disp==5)*[min(doys_0) max(doys_0)];
    mod_alphs = true;
%     sz_goal = 1000;
    dep_lat_pl = dep_lat_degs;
    xmax = 3000;
    log_cl_dist_opt = true;
     map_lat_cntr = 44.75;
     map_lon_cntr = -4; % 72.5
    FLatLims = [0 13.75];
    
case 'Hoopoe'

    min_lat = 10;
    max_lat = 50;
    min_lon = -15; % + (magn_comp_nr == 4)*40;
    max_lon = 10; % 25;
%     map_proj = 'Mercator'; % 'eqdazim'; %'gnomic'; % 'stereo'; %  'Mercator'; % 

    clims = (plot_var_disp==3)*[145 190] + (plot_var_disp==5)*[min(doys_0) max(doys_0)];
    mod_alphs = true;
%     sz_goal = 350;
    dep_lat_pl = dep_lat_degs;
    xmax = 2000;
    log_cl_dist_opt = false;  
     map_lat_cntr = 29.5; % 35;
    map_lon_cntr = 0; % -88; % -90;
    FLatLims = [0 20];
    
case 'Finn Rosefinch'
            
    min_lat = 27;
    max_lat = 62.5;
    min_lon = 20; % + (magn_comp_nr == 4)*40;
    max_lon = 80; % 25;
%     map_proj = 'Mercator'; % 'eqdazim'; %'gnomic'; % 'stereo'; %  'Mercator'; % 

    clims = (plot_var_disp==3)*[145 190] + (plot_var_disp==5)*[min(doys_0) max(doys_0)];
    mod_alphs = true;
%     sz_goal = 350;
    dep_lat_pl = dep_lat_degs;
    xmax = 2000;
    log_cl_dist_opt = false;  
     map_lat_cntr = 40; % 35;
    map_lon_cntr = 55; % -88; % -90;
    FLatLims = [0 45]; 
    
case 'BU Rosefinch'
            
    min_lat = 12;
    max_lat = 55;
    min_lon = 10; % + (magn_comp_nr == 4)*40;
    max_lon = 100; % 25;
%     map_proj = 'Mercator'; % 'eqdazim'; %'gnomic'; % 'stereo'; %  'Mercator'; % 

    clims = (plot_var_disp==3)*[145 190] + (plot_var_disp==5)*[min(doys_0) max(doys_0)];
    mod_alphs = true;
%     sz_goal = 350;
    dep_lat_pl = dep_lat_degs;
    xmax = 2000;
    log_cl_dist_opt = false;  
     map_lat_cntr = 30; % 35;
    map_lon_cntr = 50; % -88; % -90;
    FLatLims = [0 30]; 
    
    case 'RF Bluetail RU'
            
    min_lat = 24;
    max_lat = 55;
    min_lon = 100; % + (magn_comp_nr == 4)*40;
    max_lon = 150; % 25;
%     map_proj = 'Mercator'; % 'eqdazim'; %'gnomic'; % 'stereo'; %  'Mercator'; % 

    clims = (plot_var_disp==3)*[145 190] + (plot_var_disp==5)*[min(doys_0) max(doys_0)];
    mod_alphs = true;
%     sz_goal = 350;
    dep_lat_pl = dep_lat_degs;
    xmax = 2000;
    log_cl_dist_opt = false;  
     map_lat_cntr = 35; % 35;
    map_lon_cntr = 125; % -88; % -90;
    FLatLims = [0 25]; 
    
    
end