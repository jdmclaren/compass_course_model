if n_not_done  < n_inds % == 0  %  %
    %% find mean arr Lon, then frac arr then closest approach to that target Lon in "perf" case and also the dist / time covered

        if geogr_map ~= 1

            % find inds within 5 degs of arr lat
    %         finished = find(~isinf(mig_durs) & lat_es(:,i_t) <= arr_lat + pi/36);

            % find nr of steps to closest pt to arr lat 
            % (not quite same as closest dist)
            [~, i_cl_lat_arr] = min(abs(lat_es-arr_lat),[],2);   
            lin_cl = sub2ind(size(lon_es), 1:size(lon_es,1),i_cl_lat_arr');

            % ensure lat within R_thr kms of arrival (e.g., 1000 km)
%             if alfs_deg ~= 0
                finished = lat_es(lin_cl') <= arr_lat + arr_lat_thr & ...
                    abs(lon_es(lin_cl') - dep_lon) <= pi;
%             else
%                 finished = lat_es(lin_cl') <= arr_lat + arr_lat_thr & ...
%                     abs(lon_es(lin_cl') - dep_lon) <= pi;
%             end
            lin_cl = lin_cl(finished);

            arr_ts{ia,i_err}(finished) = i_cl_lat_arr(finished);
            arr_ts{ia,i_err}(~finished) = NaN;        

            % find indices of closest to arr lat
    %         [~, i_cl_lat_arr] = min(abs(lat_es-arr_lat),[],2); % ,'linear'); 


            % compute median lon when closest to arr lat
            % use this as 'target' for this configuration
            % of angles, loctions and std error

    %                         lon_es(:) = pi*shiftAnglesFromMinus180To180 ...
    %                                     (lon_es(:)*radDeg)/180;  

    %         if dep_lon == pi


                     % valid if dep_lon = pi

%                      if dep_lon == pi

                         if species_goal_opt % plot_spec_map_opt 
                             
                             % if not for given species and goal,
                            % use median long at target lat
                             % can use circ mean or median
                             % chose median since it supresses outliers which are
                             % not biologically relevant (advantageous)
                             med_lon(ia,i_err) = lon_goal_centr*pi/180; % mod((lon_goal_centr - dep_lon_sp_degs+180),360)*pi/180-pi; %; %  ( ... - spec_lon_offset) % mod(circ_median(lon_es(lin_cl)'),2*pi);
                             % circ_median(lon_es(lin_cl)'); 

%                              med_lon_disp(ia,i_err) =  med_lon(ia,i_err) - dLon_goal*pi/180;
                         
                         else
  
                             med_lon(ia,i_err) = mod(circ_median(lon_es(lin_cl)'),2*pi);
                             med_lon_disp(ia,i_err) = med_lon(ia,i_err) - dep_lon;
                             
                         end

%                      else
% 
%                          disp('need to update if dep_lon not 180')
%                          keyboard
% 
%                      end

             med_arr_ts(ia,i_err) = nanmedian(arr_ts{ia,i_err});

             % consider deciles of flight steps for closest approach
             ten_it_1 = 10*(i_t-1);
             for i_dec = 1:10
%                 dec_lons(:,i_dec:10:10*i_t) = lon_es;
                dec_lons(:,i_dec:10:ten_it_1) = ...
                    ((11-i_dec)*lon_es(:,1:i_t-1)+(i_dec-1)*lon_es(:,2:i_t))/10;
                dec_lats(:,i_dec:10:ten_it_1) = ...
                    ((11-i_dec)*lat_es(:,1:i_t-1)+(i_dec-1)*lat_es(:,2:i_t))/10;

             end
             
            dec_lons(:,ten_it_1+1) = lon_es(:,i_t);
            dec_lats(:,ten_it_1+1) = lat_es(:,i_t);
             
            % compute closest distance to goal

                dgs = distance('gc',arr_lat,med_lon(ia,i_err), ...
                     lat_es,lon_es,'radians')*R_Earth_km;
%             else
%                 
%             end
            
    %                 closer = dgs < d_close(finished);    

            % store closest distances, median dist 
            % and p(within thresh), eg 500km
            [d_close{ia,i_err},i_cl_dist_arr] = min(dgs,[],2);
            
            % compute closest distance to goal considering deciles
            dg_decs = distance('gc',arr_lat,med_lon(ia,i_err), ...
            dec_lats,dec_lons,'radians')*R_Earth_km;            
            [dec_close{ia,i_err}] = min(dg_decs,[],2);
            
            ind_cl = sub2ind(size(dgs),(1:size(dgs,1))',i_cl_dist_arr);
            
            if ~species_goal_opt
                
                  d_cl_sgn = d_close{ia,i_err}.*sign(lon_es(ind_cl)-med_lon(ia,i_err));
                  
            else
                
                 d_cl_sgn = d_close{ia,i_err}.*sign(lon_es(ind_cl) + spec_lon_offset*pi/180 -lon_goal_centr*pi/180);
                
            end
            
            for ith = 1:n_thrs
                R_thr = R_thrs(ith);
                within_thr = d_close{ia,i_err} < R_thr;
                p_close(ia,i_err,ith) = sum(within_thr)/n_inds;
                med_d_close(ia,i_err,ith) = median(d_close{ia,i_err}(within_thr));
                std_err(ia,i_err,ith) = std(d_close{ia,i_err}(within_thr));
            end

%                 med_d_close_all(ia,i_err) = median(d_close{ia,i_err});
%                 std_err_all(ia,i_err) = std(d_close{ia,i_err});
                med_d_close_all(ia,i_err) = median(dec_close{ia,i_err});
                std_err_all(ia,i_err) = std(dec_close{ia,i_err});
                
                if species_goal_opt
                    
                    p_gt_500_km = sum(dec_close{1}>500)/n_inds;
                    p_within_goal = sum(dec_close{1}<=goal_rad)/n_inds;
                    
                    tight_inh = abs(alf0_s - alfs)*180/pi <= dalf_close;
                    p_goal_tight_inher = sum((dec_close{1}<goal_rad) & ...
                        tight_inh)/sum(tight_inh);
                    med_d_tight_inher = median(dec_close{ia,i_err}(tight_inh));
                 
                    
                    
                    if is_sun_comp
                        
                        tight_sched = abs(doys_0(:) - doy_0_mn) <= n_days_close;
                         p_goal_tight_sched = ...
                             sum((dec_close{1}<goal_rad) & tight_sched)/ ...
                              sum(tight_sched);
                         
                            med_d_tight_sched(ia,i_err) = ...
                                median(dec_close{ia,i_err}(tight_sched));
                            
                            p_goal_tight_both = sum((dec_close{1}<goal_rad) & ...
                                tight_sched & tight_inh)/ ...
                              sum(tight_sched & tight_inh);
                            
                            med_d_tight_both(ia,i_err) =...
                                median(dec_close{ia,i_err}(tight_sched & tight_inh));
                            
                    else
                        
                            p_goal_tight_sched = NaN;
                         
                            med_d_tight_sched(ia,i_err) = NaN;
                            
                            p_goal_tight_all = NaN;
                            
                            med_d_tight_all = NaN;
                    
                    end
                end
                
        else

             % find inds within 500 km of arr lat                       
                finished = find(~isinf(mig_durs) & lat_es(:,i_t) <= arr_lat + arr_lat_thr);

             % find nr of steps to closest pt to arr lat 
            % (not quite same as closest dist)
            [~, arr_ts{ia,i_err}] = min(abs(lat_es-arr_lat),[],2);   

            % find indices of closest to arr lat
            % (here, dist in kms)
            [~, i_cl_lat_arr] = min(abs(lat_es-arr_lat),[],2,'linear'); 

            % compue median lon when closest to arr lat
            % use this as 'target' for this configuration
            % of angkes, loctions and std error
             med_lon(ia,i_err) = nanmedian(lon_es(i_cl_lat_arr(abs(lat_es-arr_lat)<pi/18))); 

             if isnan(med_lon(ia,i_err))
                 keyboard
             end

             med_arr_ts(ia,i_err) = nanmedian(arr_ts{ia,i_err});

             if sim_perf_opt == 1
                closest_perf_idx(ia,i_err) = ...
                    find(abs(lon_arr_p - med_lon(ia,i_err)) == ...
                    min(abs(lon_arr_p-med_lon(ia,i_err))),1,'first');
             end

            % compute closest distance to goal
            dgs = sqrt((arr_lat-lat_es).^2 + (med_lon(ia,i_err)-lon_es).^2);
    %                         dgs = distance('gc',arr_lat,med_lon(ia,i_err), ...
    %                         lat_es,lon_es,'radians')*R_Earth_km;
    %                 closer = dgs < d_close(finished);    

            % store closest distances, median dist 
            % and p(within thresh), eg 500km
            [d_close{ia,i_err},i_cl_dist_arr] = min(dgs,[],2);
            for ith = 1:n_thrs
                R_thr = R_thrs(ith);
                within_thr = d_close{ia,i_err} < R_thr;
                p_close(ia,i_err,ith) = sum(within_thr)/n_inds;
                med_d_close(ia,i_err,ith) = median(d_close{ia,i_err}(within_thr));                       
                std_err(ia,i_err,ith) = std(d_close{ia,i_err}(within_thr));
            end
    %                         within_thr = d_close{ia,i_err} < R_thr;
    %                         p_close(ia,i_err,ith) = sum(within_thr)/n_inds;

                med_d_close_all(ia,i_err) = median(d_close{ia,i_err});
            std_err_all(ia,i_err) = std(d_close{ia,i_err});


        end


        qs_i = quantile(d_close{ia,i_err}(within_thr),[0.25 0.5 0.75]);
        lq_discr(ia) = qs_i(1);
        md_discr(ia) = qs_i(2);
        uq_discr(ia) = qs_i(3);

    %                 mad_discr(i_err,ia) = mad(discr_pe);

    %                 arr_succ(ia) = 100*(1-n_not_done/n_inds);
        arr_succ(ia) = 100*p_close(ia,i_err,succ_thr);

        if plot_ind_tracks_opt % ismember(ia,plot_alf_idx_opt)

            plot_ind_tracks_spat_generic
        
        end
                    
        mean_dur_es(ia) = mean(mig_durs);

%         figure; 
%         histogram(lon_es(:,max(i_cl_lat_arr, ...
%             i_cl_dist_arr))*180/pi-dep_lon_degs)

%         figure; histogram(d_close{ia,i_err},0:250:5000,'Normalization','probability') %
        

%         figure; histogram(d_close{ia,i_err},'Normalization','probability') % ,'FaceColor','none') % (within_thr),0:100:4000)
%         hold
%         histogram(round(md_discr(ia)/100)*100,0:100:4000,'FaceColor','blue'); % d_close{ia,i_err}(~within_thr),0:100:4000)
            
else % fill in with NaN values
    
      mn_discr(ia) =  NaN;
%                 geomn_discr(i_err,ia) = geomean(abs(discr_pe)); % (:,ia)
      std_discr(ia) =  NaN;               
       mean_dur_es(ia) = NaN;
       arr_succ(ia) = 0;
       lq_discr(ia) =  NaN;
        md_discr(ia) =  NaN;
        uq_discr(ia) =  NaN;
       discr_pe = NaN;

       d_close{ia,i_err} = NaN;
%                     within_thr = NaN;
        p_close(ia,i_err,1:n_thrs) = 0;
        med_d_close(ia,i_err,1:n_thrs) = NaN;

        med_lon_disp(ia,i_err) = NaN;
        arr_ts{ia,i_err} = NaN;
       med_lon(ia,i_err) = NaN;
       med_arr_ts(ia,i_err) = NaN;
        std_err(ia,i_err,1:n_thrs) = NaN;

        med_d_close_all(ia,i_err) = NaN;
        std_err_all(ia,i_err) = NaN;
        d_cl_sgn = NaN*one_vec;
    
end

if plot_ind_tracks_opt && plot_cb_end && iprog == numel(or_progs)
    
    figure('Position',[400 100 400 300])
    cc = colorbar('EastOutside');
    
    colormap(flipud(parula))
    
    if plot_var_disp == 1
        
            caxis([1 max_fl_n])
            title(cc,'flight number','FontSize',titl_sz)

    elseif plot_var_disp == 5

           caxis([-n_del_dys n_del_dys]) %  caxis([0 100])
            title(cc,'departure day','FontSize',titl_sz)
            
    else

            caxis([160 300]) % 170 270])
            title(cc,'heading (^o Clock. N)','FontSize',titl_sz)
                        
    end
    
    axis off
    
end
