switch species

case 'generic' 
    
    % lats lons etc. already given 
    dep_lat_sp_degs = dep_lat_degs;
    % longitude arbirary for generic species
    dep_lon_sp_degs = 180;  % 36; %     -55; %      

    % goal depends on specified lat and lon ranges
    lat_goal_centr = dep_lat_degs - lat_dist;
    lon_goal_centr = 180 - lon_dist; % 170;  %  
    
    % goal radius and ndates also arbitrary for this option
    goal_rad = 1000;
    ndates = 30;    
    
%     lat_dist = 67.5;
%     max_n =  60; %  
%     n_fl_seq =   5; % 60 % 
%     stop_dur =  5; % 10; % 2; %
%     std_n_stop = 2;
%     n_hs = 8; % (2:2:10);
%     Va_mps = 12.5; %3 % 10 % 
%     alfs_deg_or_prs = [56 90 95 120]; % 90 % -120
%         dep_date = (plot_var_test==3)*[1 9 2000] + (plot_var_test==5)*[15 8 2000];
%     std_day = 14; %
%              err_cals = 15;

case 'Iql Wheatear'

        dep_lat_sp_degs =  63.7; % 82; %  
        dep_lon_sp_degs = -68.5; % -71; %  -62.5; %  210;  % 36; %     -55; %       
    %     lat_dist = 67.5;
        n_hs = 4; % 1; % 24; % 
        max_n =   120; % max(120/n_hs,40); %  
        max_t_bio = 60;
        n_fl_seq =  120;  % 
        stop_dur =  0; % 10; % 2; %
        std_n_stop = 0;

        Va_mps = 22.5; % 15; % 12.5; %17; % 3 % 10 % 
    %     alfs_deg_or_prs = [56 90 95 120]; % 90 % -120
        dep_date = [24 9 2009]; % [21 7 2000]; %  (plot_var_test==3)* + (plot_var_test==5)*[15 8 2000];
        std_day = 0; %
    %              err_cals = 15;
        lat_goal_centr = 53; % 42; %   18;
        lon_goal_centr =  -6; % -5; % -8; %  8; % -14; % 170;  %   
        goal_rad = 2500;
        % for plotting on geogr map
    %             spec_lon_offset = -10;
        ndates = 30;

case 'Alk Wheatear'
        
%         alf_inits = [-45 0 -20 -55 -55 -55]; %  
    dep_lat_sp_degs = 64; % 47.5; %  
    dep_lon_sp_degs = 202.5; %  130;  %180; %  110; %  215; % 190; %  215; %  % 36; %     -55; %       
%     lat_dist = 67.5;
    max_n =   160; % 60; %  
    max_t_bio = 160;
    n_fl_seq = 30; % 20; % 25; % 30; %    3; % 60 % 
    stop_dur = 15; % 10; % 20; % 2; % 10; %  1; %  10; %
    std_n_stop = 4; % 1; % 2; %
    n_hs = 8; % (2:2:10);
    Va_mps = 12.5; %3 % 10 % 
%     alfs_deg_or_prs = [56 90 95 120]; % 90 % -120
        dep_date = [20 8 2013]; %[1 9 2000]; % (plot_var_test==3)* + (plot_var_test==5)*[15 8 2000];
    std_day = 7; % 14; %
%              err_cals = 15;
    lat_goal_centr = 6.5; % 0; % 10; % 42; % 47.5; %65.5; % %7; %  4; %   E823 finishes at higher Lat
    lon_goal_centr = 34; % 58; % 110; % 130; %  30; %36; % 170;  %  30; %  
    goal_rad = 1000; % 500; %  250; %250; % 
    % for plotting on geogr map
%             spec_lon_offset = -10;
    ndates = 30;

case 'Monarch'

    dep_lat_sp_degs = 44.5; % 47.5;
    dep_lon_sp_degs = -81; % -130; % -70; % -100; %  -100 %  -59 %
%     lat_dist = 27.5;
    max_n =  180; % 120; %   90; %  
    max_t_bio = 120;
    n_fl_seq = 5; % 8; % 
    stop_dur = 3; % 4; % 10; %  
    std_n_stop = 1; % 2; %
    n_hs = 8; % (2:2:10);
    Va_mps = 3; % 10; %  
%     alfs_deg_or_prs = [42. 70 60 60]; %80]; % [0. 0 10 5]; % [-12. -50 -20 -27.5]; %   [42. 89.9999 60 60]; %
        dep_date = [15 8 2000]; % (plot_var_test==3)* + (plot_var_test==5)*[1 8 2000];
    std_day = 14; % 0; %   7; %
%             err_cals = 15;
    lat_goal_centr = 19+36/60;
    lon_goal_centr = -100-14/60;
    goal_rad = 100;    
    ndates = 30;

case 'Blackcap'

    dep_lat_sp_degs = 50;
%     lat_dist = 30;
    max_n =   120; % 40; % 60 % 
    max_t_bio = 90;
    n_fl_seq = 3;
    n_hs = 8; % (2:2:10);
    stop_dur = 10;
    std_n_stop = 2;
     Va_mps = 10; %3 % 10 % 
%     alfs_deg = -25;
        dep_date = [22 8 2000]; % (plot_var_test==3)* + (plot_var_test==5)*[1 8 2000];
    std_day = 14; %
%             err_cals = 15;
    goal_rad = 500;
     ndates = 60;

   case 'GreyCheekThrush NF' 
       % or some such, only nonstop 'works' along route

    dep_lat_sp_degs = 50;
    dep_lon_sp_degs = -57;
%     lat_dist = 39;
    max_n =   120; % 39; % 60 % 
    max_t_bio = 60;
    n_fl_seq = 3;
    n_hs = 24; % (2:2:10);
    stop_dur = 10;
    std_n_stop = 2;
    Va_mps = 11.5; %3 % 10 % 
%     alfs_deg_or_prs = [15.5 25.5 18.5 20.5]; % -40]; %-62.5] % 
        dep_date = [15 9 2000]; % (plot_var_test==3)* + (plot_var_test==5)*[22 8 2000];
        std_day = 7; %
    lat_goal_centr = 11;
    lon_goal_centr = -71;
    goal_rad = 1000;
%             err_cals = 15;
%             spec_lon_offset = 40;
    ndates = 30;

case 'GreyCheekThrush'  

    dep_lat_sp_degs =  65; %
    dep_lon_sp_degs = -140;
%     lat_dist = 60;
    max_n = 180; % 60; %   120; %  40; % 
    max_t_bio = 120;
    n_fl_seq = 5;
    n_hs = 8; % (2:2:10);
    stop_dur = 5;
    std_n_stop = 2; % 0; % 
    Va_mps = 11.5; %3 % 10 % 
%     alfs_deg_or_prs = [-40 -87.5 -5 -55]; % -40]; %-62.5] % 
        dep_date = [10 9 2000]; % (plot_var_test==3)* + (plot_var_test==5)*[24 8 2000];
    std_day =  7; % 14; %  3; % 
    lat_goal_centr =  0; %40; %
    lon_goal_centr = -70;
    goal_rad = 1000;
%             err_cals = 15;
%             spec_lon_offset = 40;
    ndates = 30;

case 'Whinchat'  

    dep_lat_sp_degs = 65;
%     lat_dist = 65;
    max_n =   120; % 40; % 60 % 
    max_t_bio = 60;
    n_fl_seq = 5;
    n_hs = 8; % (2:2:10);
    stop_dur = 10;
    std_n_stop = 2;
    Va_mps = 10; %3 % 10 % 
%     alfs_deg = 75; % 55 %65 % -50 % -15 % -40 %           
        dep_date = [22 8 2000] ; %  (plot_var_test==3)*+ (plot_var_test==5)*[1 8 2000];
    std_day = 7; %
%             err_cals = 15;
    goal_rad = 1000;
    ndates = 30;

case 'BluefinTuna'

    dep_lat_sp_degs = 37.5; % 45;
    dep_lon_sp_degs = 142.5;
%     lat_dist = 5 % 25;
    max_n =  200; % 60 % 
    max_t_bio = 100;
    n_fl_seq = 5;
    n_hs = 8; % 12; % (2:2:10);
    stop_dur = 10; % 5;
    std_n_stop = 2;
    Va_mps = 25; %3 % 10 % 
%     alfs_deg_or_prs = [-80 -90 -75 -90]; %  
        dep_date = [1 9 2000]; % (plot_var_test==3)* + (plot_var_test==5)*[1 6 2000]; % [1 8 2000] % [1 9 2000] %
    std_day = 14; %
%             err_cals = 15;
    lat_goal_centr = 32.5 % 27.5;
    lon_goal_centr = -125;
    goal_rad = 1000;
%              spec_lon_offset = -40;  
     ndates = 60;

case 'Siberian Willow Warbler'

    dep_lat_sp_degs = 67.5; 
    dep_lon_sp_degs = 170;  % 36; %           
%     lat_dist = 67.5;
    max_n =   120; % 60; % 
    max_t_bio = 120;
    n_fl_seq = 5; %  8; %60 % 
    stop_dur = 2; % 5; % 2; % 5; % 
    std_n_stop = 2;
    n_hs = 8; % (2:2:10);
    Va_mps = 10.5; %3 % 10 % 
%     alfs_deg_or_prs = [56 90 95 120]; % 90 % -120
        dep_date = [1 9 2000]; % (plot_var_test==3)*[1 9 2000] + (plot_var_test==5)*[15 8 2000];
    std_day = 7; %
%              err_cals = 15;
    lat_goal_centr = 7.5;
    lon_goal_centr = 35; % 37.5; % 36; % 170;  %   
    goal_rad = 600; % 750;
    % for plotting on geogr map
%             spec_lon_offset = -10;
    ndates = 30;
    
    case 'Sib Will Warb S Hem'

    dep_lat_sp_degs = 67.5; 
    dep_lon_sp_degs = 170;  % 36; %           
%     lat_dist = 67.5;
    max_n = 180; % 120; %    60; %  80; %
     max_t_bio = 120;
    n_fl_seq =  5; % 2; %  8; % 60 %
    stop_dur =   2; % 
    std_n_stop = 2;
    n_hs = 8; % (2:2:10);
    Va_mps = 10.5; %3 % 10 % 
%     alfs_deg_or_prs = [56 90 95 120]; % 90 % -120
        dep_date = [1 9 2000]; % (plot_var_test==3)*[1 9 2000] + (plot_var_test==5)*[15 8 2000];
    std_day = 7; %
%              err_cals = 15;
    lat_goal_centr = -12.5;
    lon_goal_centr = 32.5; % 36; % 170;  %   
    goal_rad = 1000;
    % for plotting on geogr map
%             spec_lon_offset = -10;
    ndates = 30;

case 'SEPacific Humpback' % need to flip for plot and also change season!

    dep_lat_sp_degs = 65; 
    dep_lon_sp_degs = -70;
%     lat_dist = 65;
    max_n =  120; %  
    max_t_bio = 120;
    n_fl_seq = 5;
    stop_dur = 1;
    std_n_stop = 1;
    n_hs = 24; % (2:2:10);
    Va_mps = 1.25; %3 % 10 % 
%     alfs_deg_or_prs = [13. 31 34 32.5] % 40]; % 50 % 60 %  
        dep_date = [15 9 2000]; % (plot_var_test==3)* + (plot_var_test==5)*[15 8 2000];
    std_day = 14; %
%             err_cals = 15;
    lat_goal_centr = 0;
    lon_goal_centr = -90;        
    goal_rad = 1000;
%             spec_lon_offset = 110;
     ndates = 60;

case 'Kirtlands Warbler'

    dep_lat_sp_degs = 45;
    dep_lon_sp_degs =  -85; % -79.5; % 
%     lat_dist = 20;
    max_n =  120; % 180; % 60; %  
    n_fl_seq = 5; % 6;
    max_t_bio = 30; % 60;
    stop_dur =  5; %15 2;
    std_n_stop = 2; % 1;
    n_hs =  8; % 6; %(2:2:10);
    Va_mps =  10; % 
%     alfs_deg_or_prs =  [-12 -15 -10 -13]; % [-5 -10 -2 -5]; %              
        dep_date = [24 9 2000]; %[6 10 2000]; % (plot_var_test==3)* + (plot_var_test==5)*[15 9 2000];
    std_day = 5; % 7; %
%             err_cals = 15;
    lat_goal_centr = 25; % 24.5;
    lon_goal_centr = -76; % -78;
    goal_rad = 300;
    % for plotting on geogr map
%             spec_lon_offset = 95; % 101; % 
       ndates = 30;

    case 'Spotted Flycatcher'

        dep_lat_sp_degs = 65;
        dep_lon_sp_degs = 20;
%         lat_dist = 60;
        max_n =   120; % 60; % 60 % 
        max_t_bio = 90;
        n_fl_seq = 5;
        n_hs = 8; % (2:2:10);
        stop_dur = 10;
        std_n_stop = 2;
        Va_mps = 11.5; %3 % 10 % 
%         alfs_deg_or_prs = [-10 -20 7.5 -2.5]; % -40]; %-62.5] % 
        dep_date = [22 8 2000]; % (plot_var_test==5)* + (plot_var_test==3)*[1 8 2000];
        std_day = 7; %
        lat_goal_centr = 0;
        lon_goal_centr = 35;
        goal_rad = 1000;
%             err_cals = 15;
%             spec_lon_offset = 40;
        ndates = 30;

        case 'Finn Marsh Warbler'

        dep_lat_sp_degs = 62.5;
        dep_lon_sp_degs = 22.5;
%         lat_dist = 60;
        max_n =  120; % 180; % 60; % 70 % 
        max_t_bio = 90;
        n_fl_seq = 3;
        n_hs = 8; % (2:2:10);
        stop_dur =  3; % 1; %
        std_n_stop =  2; % 0; %
        Va_mps = 11.5; %3 % 10 % 
%         alfs_deg_or_prs = [-12.5 -25 10 -2.5]; % -40]; %-62.5] % 
        dep_date = [15 8 2000]; %(plot_var_test==3)*[1 9 2000] + (plot_var_test==5)*[1 8 2000];
        std_day = 7; %
        lat_goal_centr = 2.5; % -2.5;
        lon_goal_centr = 22.5; % 40;
        goal_rad = 500;
%             err_cals = 15;
%             spec_lon_offset = 40;
          ndates = 30;

        case 'France Marsh Warbler'

        dep_lat_sp_degs = 65; % 49.5;
        dep_lon_sp_degs = 42.5; % -1.5;
%         lat_dist = 60;
        max_n =   120; % 60; % 60 % 
        max_t_bio = 90;
        n_fl_seq = 5;
        n_hs = 8; % (2:2:10);
        stop_dur = 5;
        std_n_stop = 2;
        Va_mps = 11.5; %3 % 10 % 
%         alfs_deg_or_prs = [-12.5 -25 10 -2.5]; % -40]; %-62.5] % 
        dep_date = [15 8 2000]; %(plot_var_test==3)*[1 9 2000] + (plot_var_test==5)*[1 8 2000];
        std_day = 7; %
        lat_goal_centr = 54.5; % 2.5; % -2.5
        lon_goal_centr = 18; % 32.5; % 40 
        goal_rad = 100;
%             err_cals = 15;
%             spec_lon_offset = 40;
          ndates = 30;

    case 'Nathusius'

         dep_lat_sp_degs = 56.2;
        dep_lon_sp_degs = 21;
%         lat_dist = 13.8;
        max_n =  120; % 180; % 60; % 60 % 
        max_t_bio = 90;
        n_fl_seq = 3;
        n_hs = 6; % (2:2:10);
        stop_dur = 5;
        std_n_stop = 2;
        Va_mps = 7.5; %3 % 10 % 
%         alfs_deg_or_prs = [43. 55.5 53 53]; % -40]; %-62.5] % 
        dep_date = [25 8 2000]; % (plot_var_test==3)* + (plot_var_test==5)*[16 8 2000];
        std_day = 14; %
        lat_goal_centr = 43; % 42.4;
        lon_goal_centr = 1; % 2.5;
        goal_rad = 300;
%             err_cals = 15;
%             spec_lon_offset = 40;
          ndates = 20;             
          
       case 'Eleonoras Falcon'

        dep_lat_sp_degs = 38.5;
        dep_lon_sp_degs = 8;
%         lat_dist = 13.8;
        max_n =   120; % 60; % 60 % 
        max_t_bio = 90;
        n_fl_seq = 9;
        n_hs = 11; % (2:2:10);
        stop_dur = 14;
        std_n_stop = 2;
        Va_mps = 20; %3 % 10 % 
%         alfs_deg_or_prs = [43. 55.5 53 53]; % -40]; %-62.5] % 
        dep_date = [28 10 2000]; % (plot_var_test==3)* + (plot_var_test==5)*[16 8 2000];
        std_day = 2; %
        lat_goal_centr = -17; % 42.4;
        lon_goal_centr = 45; % 2.5;
        goal_rad = 600;
%             err_cals = 15;
%             spec_lon_offset = 40;
          ndates = 20;      
          
    case 'Ring Ouzel'  
        
        dep_lat_sp_degs =  57; %
        dep_lon_sp_degs = -3.5;
    %     lat_dist = 60;
        max_n =  120; % 180; % 40; % 60 % 
        max_t_bio = 90;
        n_fl_seq = 10; % 5; % 
        n_hs = 8; % (2:2:10);
        stop_dur =  15; % 5; %
        std_n_stop =  5; % 2; %
        Va_mps = 11.5; %3 % 10 % 
    %     alfs_deg_or_prs = [-40 -87.5 -5 -55]; % -40]; %-62.5] % 
        dep_date = [31 8 2000]; % (plot_var_test==3)* + (plot_var_test==5)*[24 8 2000];
        std_day =  7; % 14; %  3; % 
        lat_goal_centr =  33.5; %40; %
        lon_goal_centr = -4;
        goal_rad = 250;
        
    case 'Hoopoe'

        dep_lat_sp_degs = 46;
        dep_lon_sp_degs = 7.5;
    %     lat_dist = 30;
        max_n =  120; % 180; %40; % 60 % 
        max_t_bio = 60;
        n_fl_seq = 5;
        n_hs = 8; % (2:2:10);
        stop_dur = 5;
        std_n_stop = 2;
         Va_mps = 12; %3 % 10 % 
    %     alfs_deg = -25;
          dep_date = [10 8 2000]; % (plot_var_test==3)* + (plot_var_test==5)*[1 8 2000];
        std_day = 7; %
    %             err_cals = 15;
        goal_rad = 800;
         ndates = 60;
         lon_goal_centr = -5; % -2.5;
        lat_goal_centr = 17.5; % 40;
        goal_rad = 800;
        
case 'Finn Rosefinch'

    dep_lat_sp_degs = 60; 
    dep_lon_sp_degs = 23;  % 36; %           
%     lat_dist = 67.5;
    max_n =  120; % 60; % 
    max_t_bio = 90;
    n_fl_seq = 5; %  8; %60 % 
    stop_dur = 2; % 5; % 2; % 5; % 
    std_n_stop = 5;
    n_hs = 8; % (2:2:10);
    Va_mps = 12.5; %3 % 10 % 
%     alfs_deg_or_prs = [56 90 95 120]; % 90 % -120
        dep_date = [1 8 2000]; % (plot_var_test==3)*[1 9 2000] + (plot_var_test==5)*[15 8 2000];
    std_day = 7; %
%              err_cals = 15;
    lat_goal_centr = 15;
    lon_goal_centr = 77; % 36; % 170;  %   
    goal_rad = 300;
    % for plotting on geogr map
%             spec_lon_offset = -10;
    ndates = 30;
    
case 'BU Rosefinch'

    dep_lat_sp_degs = 42; 
    dep_lon_sp_degs = 23.5;  % 36; %           
%     lat_dist = 67.5;
    max_n =  120; % 180; % 60; % 80; % 
    max_t_bio = 60;
    n_fl_seq = 5; %  8; %60 % 
    stop_dur = 2; % 5; % 2; % 5; % 
    std_n_stop = 5;
    n_hs = 8; % (2:2:10);
    Va_mps = 12.5; %3 % 10 % 
%     alfs_deg_or_prs = [56 90 95 120]; % 90 % -120
        dep_date = [7 8 2000]; % (plot_var_test==3)*[1 9 2000] + (plot_var_test==5)*[15 8 2000];
    std_day = 7; %
%              err_cals = 15;
    lat_goal_centr = 24;
    lon_goal_centr = 75; % 36; % 170;  %   
    goal_rad = 400;
    % for plotting on geogr map
%             spec_lon_offset = -10;
    ndates = 30;
    
case 'RF Bluetail RU'
        
    dep_lat_sp_degs = 50; 
    dep_lon_sp_degs = 145;  % 36; %           
%     lat_dist = 67.5;
    max_n =   120; % 40; % 
    max_t_bio = 90;
    n_fl_seq = 5; %  8; %60 % 
    stop_dur = 2; % 5; % 2; % 5; % 
    std_n_stop = 5;
    n_hs = 8; % (2:2:10);
    Va_mps = 10.5; %3 % 10 % 
%     alfs_deg_or_prs = [56 90 95 120]; % 90 % -120
        dep_date = [15 9 2000]; % (plot_var_test==3)*[1 9 2000] + (plot_var_test==5)*[15 8 2000];
    std_day = 7; %
%              err_cals = 15;
    lat_goal_centr = 28;
    lon_goal_centr = 110; % 36; % 170;  %   
    goal_rad = 500;
    % for plotting on geogr map
%             spec_lon_offset = -10;
    ndates = 30;
    
end

lat_dist = dep_lat_sp_degs - lat_goal_centr;

dep_lat_degs = dep_lat_sp_degs;
