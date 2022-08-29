% initiate updated location to stepwise departure location
lat_es(:,i_t+1) = lat_es(:,i_t);
lon_es(:,i_t+1) = lon_es(:,i_t);

current_alf = NaN*one_vec;
pref_alf = current_alf;

% Feb 2021 added 'drift' error which we for simplicity assume remains
% constant during flight
% try
if ~polar_magn_option || sum(~is_sunset) == 0
    
    % add to pref dir if no hourly drift is incorporated
    if err_hr_drft == 0
    
         pref_alf(~done) = alf0_s(~done)+dalf(~done)+err_dfs(~done);
         
    else
        
         pref_alf(~done) = alf0_s(~done)+dalf(~done);

    end

else
                          
    %     if sum(~is_sunset & lat_es < 66*pi/180) > 0
    %         keyboard
    %     end
    pref_alf(is_curr_sun_comp,1) = alf0_s(is_curr_sun_comp)+ ...
        dalf(is_curr_sun_comp) +err_dfs(is_curr_sun_comp);
    % note while uses (default) mag compass it still (if tc_comp)
    % accumulates 'jetlag' for when it switches back (tc_sun_comp)
    pref_alf(is_curr_magn_comp,1) = alf0_s(is_curr_magn_comp) + ...
        err_dfs(is_curr_magn_comp);

end
% catch
%         keyboard
% end

% if compass is transferred between detection and maintenance,
% add det err (err_s) to preferred in-flight headings
if ~is_trans_comp 
    
%     if ~isinf(kappa_err_mnt) % no 'offset' error due to detection
        % Det error is first maintenace error (though could be different std's)
        err_ms(~done,1) = err_dtcs(~done); % (is_curr_sun_comp);
%     else
%         err_ms(~done,1) = zeros(n_not_done,1); 
%     end
        
    % pref_alf remains same (i.e., not 'biased' by first detection)

else % transferred to other comp

   % transfer not maintenance error on departure
   err_ms(~done,1) = zeros(n_not_done,1); 
     
   pref_alf(~done) =  pref_alf(~done) + err_dtcs(~done) + err_trs(~done);
        
end

% initialize the current direction
% add drift if incorporated hourly
if err_hr_drft ~= 0

   current_alf(~done) = pref_alf(~done) + err_ms(~done,1) + err_dfs(~done); %
       
else % err_dfs is already incorporated in preferred dirn
    
    current_alf(~done) = pref_alf(~done) + err_ms(~done,1) ; %
    
end

 if geogr_map == 1

    % inifinite plane - no lat dependence of 
    % change in Lon
        lat_es(~done,i_t+1) = lat_es(~done,i_t+1) - ...
            nightly_del_Lat*cos(current_alf(~done));

        lon_es(~done,i_t+1) = lon_es(~done,i_t+1) - ...
            nightly_del_Lat*sin(current_alf(~done));

elseif geogr_map == 2
    
        % use igrf field to obtain geographic locations from geomagnetic headings
        
        % initial headings already adjusted for decln 
        % (if geomagn heads, dalf was adjusted by decln at departure, if
        % geogr heads with star compass, dalf was unadjusted for departure)
%         try
            d_decln = zeros(sum(~done),1);

            for iH = 1:n_hs
%                             
%                current_alf(~done) = mod(pi + ...
%                    current_alf(~done) + d_decln,2*pi)-pi; % 
               
                dLon = hourly_del_Lat.* ...
                    (sin(pi+current_alf(~done)))./ ...
                    cos(lat_es(~done,i_t+1));
                dLat = hourly_del_Lat.*(cos(pi+current_alf(~done)));

                lon_es(~done,i_t+1) = lon_es(~done,i_t+1) + dLon; % mod( ,2*pi);
                lat_es(~done,i_t+1) = lat_es(~done,i_t+1) + dLat;

                % account for polar crossing
                lat_es(lat_es > pi/2) = pi - lat_es(lat_es > pi/2);
                lat_es(lat_es < -pi/2) = -pi - lat_es(lat_es < -pi/2);

                % account for datelne passing
                % don't shift angles since start at pi and head to 0 or
                % 2*pi depending on heading
%                 lon_es(~done,i_t+1) = shiftAnglesFromMinus180To180(lon_es(~done,i_t+1)*180/pi)*pi/180;
                clear Bx By
                if ~star_night
                    % update declination hourly if magnetic loxodrome
                     % unique(date_jul(~finished));
        %                       dy_no_fin = date_jul(~finished);
                    [Bx, By, ~] = igrf(med_date, ...
                           lat_es(~done,i_t+1)*180/pi, lon_es(~done,i_t+1)*180/pi, 0);

                    decln_new = atan2(By,Bx);
                    d_decln = decln_new-decln(~done);      

                end
                
                % first update (maintenance) error, then adjust heading
                % first sub-step accounted for above
                if ~isinf(kappa_err_mnt)
                    
%                     if ~comb_comp
                        err_ms(~done) = mod(err_ms(~done)*ar_error_mnt + ...
                         vmrand(0, kappa_err_mnt, [n_not_done 1]) +pi,2*pi) -pi; 
%                     else
%                         err_ms(~done) = mod(err_ms(~done)*ar_error_mnt + ...
%                              (vmrand(0, kappa_err_mnt, [n_not_done 1]) + ...
%                          vmrand(0, kappa_err_mnt, [n_not_done 1]))/2 +pi,2*pi) -pi;                        
%                     end
                     
                end
                
                
                if err_hr_drft ~= 0 
                    
%                     if ~isinf(kappa_err_hr_drft)

                        err_dfs(~done) = mod(err_dfs(~done)*ar_err_hr_drft + ...
                         vmrand(0, kappa_err_hr_drft, [n_not_done 1]) +pi,2*pi) -pi; 
                        
%                     else
%                         
%                         err_dfs(~done) =  ...
%                          vmrand(0, kappa_err_hr_drft, [n_not_done 1]);                         
%                         
%                     end

                end
                
                % update current direction, given change in decln and
                % maintenance error
                 current_alf(~done) = mod(pi + pref_alf(~done) + ...
                 err_ms(~done)+ d_decln + err_dfs(~done),2*pi)-pi; % 
                 
            end
            
%         catch
%             keyboard
%         end


else % geogr_map == 3 

    % use inverse Haversine eqn to update position

    for ih = 1:n_hs

%         try
        [lat_es(~done,i_t+1), lon_es(~done,i_t+1)] = ...
            invHaversine(hourly_del_Lat,current_alf(~done)+pi, ...
            lat_es(~done,i_t+1), lon_es(~done,i_t+1));
%         catch
%             keyboard
%         end

        % compute goal distance and update if minimal
%                             dgs = distance('gc',[arr_lat,lon_arr_p(ia)], ...
%                             [lat_es(~done,i_t+1),lon_es(~done,i_t+1)],'radians');
%                             closer = dgs < d_close(~done);                
%                             d_close(~done) = min(d_close(~done),dgs);

    end   

 end

% lon_es = mod(lon_es+pi,2*pi)-pi;
lat_es(done,i_t+1) = NaN; % th_ps(done,i_t);
lon_es(done,i_t+1) = NaN; % ll_ps(done,i_t);
% lat_es(~done,i_t+1) = lat_es(~done,i_t);
% lon_es(~done,i_t+1) =  lon_es(~done,i_t);  
