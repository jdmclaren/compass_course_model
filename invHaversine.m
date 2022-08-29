function [lat_bs, lon_bs] = invHaversine(arc_dist,track_dirns,lat_bs,lon_bs)         
% use inverse Haversine eqn to update position
               
%             arc_dist = g_sp*TR; % arc dist travelled % 
            sdtb = sin(arc_dist);
            cdtb = cos(arc_dist);
            c_pd = cos(track_dirns);
            s_pd = sin(track_dirns);
            
            s_lat = sin(lat_bs);
            c_lat =  cos(lat_bs);
            
            lat_bs = asin(s_lat.*cdtb + c_lat.*sdtb.*c_pd);
            
%             Y_tb = s_pd.*sdtb.*c_lat;
%             X_tb= cdtb - s_lat.*sin(lat_bs(idx_locs));
%             dlon_tb = atan2(Y_tb,X_tb) + (X_tb<0)*pi;
%             lon_bs(idx_locs) = ~Stopped.*(mod(lon_bs(idx_locs_pre) ...
%                 + dlon_tb+pi,2*pi) -pi) + Stopped.* ...
%                 Stop_poly(min_Stops,1);

% use cheaper formula for Haversine projection
           as_term =  asin(s_pd.*sdtb./cos(lat_bs));
           
           % for simulations with no topography
           % don't use modulo for longitude
           % as we want to rule out global wraparounds
          lon_bs = lon_bs + as_term;
%            
%            lon_bs = mod(lon_bs + as_term+pi,2*pi) -pi;