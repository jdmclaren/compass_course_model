% plot_magn_tran
addpath('D:\Oldenburg_models\geomagn_orientn_model\brewer')

if ~exist('plot_for_ms')
    
    plot_for_ms = true; % false;
    
end

if plot_for_ms
    left = 300;
    bottom = 300;
    width = 260; % 200;
    height = 210; % 220;
    FntSz = 9;
    incl_cb = false;
    incl_titl = false;
else
    FntSz = 16;
    incl_cb = true;
    incl_titl = true;  
end

show_polar_night_day = false; % true; % 

% for d[alf(i)]/d[alf(i-1)] sensitivity
% 1 is mgcl, 9 is general (lox) dLong/dAlph(1 deg off)
% 2, 3 are fixed and tc (all) and 11 is tc components of tc
% 12 = zero field representing orientation sens of Lox
plot_choices = [1 2 3 9] %  1 % 1:3 % 12 % 9 % [1 2 3 9] %9  11 % [1 9] % 11 %  9 %  [1  9 ] %  [1  9 ] % 

n_fls = 1; % 10 % to (artificially) extend the Long. shift delta(llamda)

SE_SW =  2 %  1 %
% sign_alf = (SE_SW ==2) - (SE_SW ==1);

max_lat = 90;
r_max_lat = round(max_lat);
max_alf =  120; % 135; %  89.9;
r_max_alf = round(max_alf);

% nr of angles per degree
resln = 10; % 1; % 

% half is small angle to avoid singularities
half = 0.5;

is_mixed_sun = false % true
mixed_fact = (is_mixed_sun) +1;

cax_lo = -2;
cax_hi = 2;


R_Earth_km = 6317;
Va_ms = 12.5;
fl_hrs =  8; % 5 gives 225 km and 8 gives 360 km per night

% stepwise flight distance
r = fl_hrs*Va_ms*3.6/R_Earth_km; % corrcted previous erroneous multiplication by 180/pi !!
% since using sine cosine we need R in radians 

dates =  [213 274] % 270 % [235 258] %  [213 265 274] %  272] %[213 265 304] % 213  % 265 %   228 % 264 % [345 350 355 360] % 262:2:270  % 288 % 228 % 258  % Sep 15th non-leap year

for idate  = 1:numel(dates)
    
    date = dates(idate);

    for iilat = 1:1+resln*max_lat % 30 % 

        ilat = (iilat-1)/resln;
        
         [sun_az, ~, is_sunset] = calc_sun_az(ilat*pi/180,date);
         
           if is_sunset
                      
                is_suns(iilat) = true;
                    
           else
               
               is_suns(iilat) = false;
                    
           end
            
            for jjalf = 1:(2*resln*max_alf+1)

                jalf = -max_alf + (jjalf-1)/resln;

                ij_lat = ilat-r*cosd(jalf)*180/pi;
                ij_lon = -n_fls*r*sind(jalf)/cosd(ilat)*180/pi;
                if ij_lat > 90
                    ij_lat = 180 - ij_lat;
    %                 ij_lon = -ij_lon;
                end


                tan_alph_sun = tan(calc_sun_az(ij_lat*pi/180,date)); % in rads
                if imag(tan_alph_sun) ~= 0
                    keyboard
                end

                % long error loxodrome 
                dlon_alf_lox(iilat,jjalf) = -r*cosd(jalf)/cosd(ilat);

                % vs prev heading alf (for sun compass, 
                % using first step since iterative)

%                 dalf_alf_mgcl(iilat,jjalf) = r*sind(jalf)^2/ ...
%                      (cosd(ilat)*sqrt( ...
%                      tand(ilat)^2*cosd(ij_lat)^2-sind(jalf)^2*sind(ij_lat)^2));

                dalf_alf_mgcl(iilat,jjalf) = r*sind(jalf)^2/ ...
                     (cosd(ilat)*sind(ilat)*cosd(jalf));
                 
                  if is_sunset

                    dalf_alf_fix(iilat,jjalf) =  r*sind(jalf)* ...
                        tand(ij_lat)/tan_alph_sun/mixed_fact;

                    dalf_alf_tcs(iilat,jjalf) =  dalf_alf_fix(iilat,jjalf) + ...
                        r*(sind(jalf)*ij_lon*cosd(ij_lat)*pi/180  ...  %           
                                   - cosd(jalf)*sind(ij_lat)/cosd(ilat))/mixed_fact;
                                % 
                    dalf_tcs_1(iilat,jjalf) =   ...
                        r*(sind(jalf)*ij_lon*cosd(ij_lat)*pi/180); % 0; %

                    dalf_tcs_2(iilat,jjalf) =  ...         
                                    -r*(cosd(jalf)*sind(ij_lat)/cosd(ilat));
                            
                 else
             
                    dalf_alf_fix(iilat,jjalf) =  NaN;

                    dalf_alf_tcs(iilat,jjalf) =  NaN;
                                % 
                    dalf_tcs_1(iilat,jjalf) =  NaN; % 0; %

                    dalf_tcs_2(iilat,jjalf) =  NaN;                     

                 end

                % vs init heading alf_0

                 tan_gam_2 = 4*tand(ilat)^2/sind(jalf)^2;    

                dalf_alf0_mgcl(iilat,jjalf) = -2*tand(ij_lat)* ...
                     (1+1/tan_gam_2)/cosd(jalf);


                dalf_alf0_sun(iilat,jjalf) =  1-dalf_alf_fix(iilat,jjalf);


                % vs lat perturbn

                dalf_lat_mgcl(iilat,jjalf) = 2/(cosd(ilat)^2* ...
                    sqrt(tan_gam_2 - 4*tand(ilat)^2));

                dalf_lat_fix(iilat,jjalf) = -tand(ilat)/tan_alph_sun/mixed_fact;

                dalf_lat_tcs(iilat,jjalf) = dalf_lat_fix(iilat,jjalf) + ...
                     ij_lon*cosd(ij_lat)/mixed_fact;

                % lon from lon perturbation
                dlon_lon_tcs(iilat,jjalf) = dlon_alf_lox(iilat,jjalf)*sind(ilat);
             

                 
            end
                 
    end

    dalf_alf_mgcl = dalf_alf_mgcl*180/pi;
%     dalf_alf_tcs = dalf_alf_tcs*180/pi;
%     dalf_alf_fix = dalf_alf_fix*180/pi;

    all_lats = 0:1/resln:max_lat;
    if sum(~is_suns) == 0 || show_polar_night_day 
        % always a sunset or we plot polar nt / day

        all_lat_suns =  all_lats;
        Lat_sun_tiks = 0:30:90;
        
    else
        
     % remove high lats from consideration
        last_lat_rec = find(is_suns,1,'last');
        all_lat_suns = all_lats(1:last_lat_rec);
        last_sun_lat = (last_lat_rec-1)/resln;
        
        dlon_lon_tcs = dlon_lon_tcs(1:last_lat_rec,:);
        dalf_lat_tcs = dalf_lat_tcs(1:last_lat_rec,:);
        dalf_lat_fix = dalf_lat_fix(1:last_lat_rec,:);
        dalf_alf0_sun = dalf_alf0_sun(1:last_lat_rec,:);
        dalf_alf_tcs = dalf_alf_tcs(1:last_lat_rec,:);
        dalf_alf_fix = dalf_alf_fix(1:last_lat_rec,:);

        % add last 10s of degrees e.g. 70 or 80 
        Lat_sun_tiks = 0:30:last_sun_lat;
        if Lat_sun_tiks(end) < floor(last_sun_lat/10)*10
           Lat_sun_tiks(end+1) =  floor(last_sun_lat/10)*10;
        end
        
    end
 
    % [xx, yy] = meshgrid(-max_alf:max_alf,1:max_lat); % (-pi/2+pi/180:pi/180:pi/2-pi/180,-89.5:89.5); %

    %% plots vs previous heading alf

    % first plot any non-solar ones for first date only
    if idate == 1
            if ismember(1,plot_choices)
                if plot_for_ms
                    figure('Position',[left bottom width height]); %
                else
                    figure
                end
                % contourf(xx,yy,log10(dalf_alf_mgcl)) % .*sign(dthdg))
                imagesc(-max_alf:1/resln:max_alf,all_lats,log10(1+abs(dalf_alf_mgcl))) % .*sign(dthdg))
                set(gca,'YDir','normal')
                % xlabel('transverse projection inclination (^o)')
                ylabel('Latitude (^o)')
                set(gca,'YTick',0:30:90); 
%                 set(gca,'XTick',-90:30:90,'XTickLabel',{'90^o','60^o','30^o','0^o','-30^o','-60^o','-90^o'})
                set(gca,'XTick',0:30:90,'XTickLabel',{'0^o','30^o','60^o','90^o'})
%                 if incl_cb

            %     colormap(flipud(summer))
                % title({'log sensitivity latitude', 'to transv magncl param'})
                
                
%                     ytklabs = {[0.01 0.1 1 10 100 500]}; % {[0.001 0.01 0.1 1 10 100 1000]}; % {'-
%                     set(cc,'ytick',[-2:2 log10(500)],'yticklabel',ytklabs); % -10.^(-2:0.5:2))
%                 end
%                 if incl_title
                  
%                 end

                colormap(brewermap([],'YlOrRd')); % ,'Reds')); % % colormap(flipud(hot(512)));
%                 cm = cm(80:450,:);
%                 colormap(cm)
                % title({'log sensitivity latitude', 'to transv magncl param'})
%                 caxis([cax_lo log10(500)])
                axis([0 90 0 90])
                % title('(e)')
                set(gca,'FontSize',FntSz)
                hold
                
                caxis([0 log(5)]) % log10(101)
                
                if  ~plot_for_ms
                    
                     xlabel('West                   Heading                    East')

                    
                else
                    
                    xlabel('East or West Heading') % xlabel('West     Heading     East')
%                     xtickangle(45)
                    
                    % new figure for colorbar
                    figure('Position',[left bottom width height*0.85]); %
                   colormap(brewermap([],'YlOrRd')); % ,'Reds')); % 'OrRd')); % 'Reds' 'YlOrRd' '*RdYlBu')); % c
                   caxis([0 log(4)]) % log10(101)
                   axis('off')
                   
                end
                
                     cb_mgcl = colorbar;
                    
                   %     set(cc,'ytick',-2:2,'yticklabel',ytklabs); %
                    ytklabs = {'0','25','100','250'}; % [0 1 2.5 5 10 20 50]}; %  100
                    
                    set(cb_mgcl,'ytick',[0 log(1.25) log(2.) log(3.5) ],'yticklabel',ytklabs); %  log10(101)
                    title(cb_mgcl,{'Sensitivity (%)'})
                
            end
            
            if ismember(4,plot_choices)
                if plot_for_ms
                    figure('Position',[left bottom width height]); %
                else
                    figure
                end
                % contourf(xx,yy,log10(dalf_alf_mgcl)) % .*sign(dthdg))
                imagesc(-max_alf:1/resln:max_alf,all_lats,log10(-dalf_alf0_mgcl)) % .*sign(dthdg))
                set(gca,'YDir','normal')
                % xlabel('transverse projection inclination (^o)')
                xlabel('West                   Heading                    East')
                ylabel('Latitude (^o)')
                set(gca,'YTick',0:30:90); 
                set(gca,'XTick',-90:30:90,'XTickLabel',{'90^oE','60^oSE','30^oSE','0^o','30^oSW','60^oSW','90^oW'})
                cc = colorbar;
                ytklabs = {[0.01 0.1 1 10 100]}; % {[0.001 0.01 0.1 1 10 100 1000]}; % {'-
                set(cc,'ytick',cax_lo:cax_hi,'yticklabel',ytklabs); % -10.^(-2:0.5:2))
                cm = colormap(flipud(hot(512)));
                cm = cm(80:450,:);
                colormap(cm)
                % title({'log sensitivity latitude', 'to transv magncl param'})
                caxis([cax_lo cax_hi])
                % title('(e)')
                set(gca,'FontSize',FntSz)
                axis([-r_max_alf r_max_alf 0 r_max_lat])
                hold
            end
            
            if ismember(6,plot_choices)
                if plot_for_ms
                    figure('Position',[left bottom width height]); %
                else
                    figure
                end
                % contourf(xx,yy,log10(dalf_alf_mgcl)) % .*sign(dthdg))
                imagesc(-max_alf:1/resln:max_alf,all_lat_suns,log10(dalf_lat_mgcl)) % .*sign(dthdg))
                set(gca,'YDir','normal')
                % xlabel('transverse projection inclination (^o)')
                xlabel('West                   Heading                    East')
                ylabel('Latitude (^o)')
                set(gca,'YTick',0:30:90); 
                set(gca,'XTick',-90:30:90,'XTickLabel',{'90^oE','60^oSE','30^oSE','0^o','30^oSW','60^oSW','90^oW'})
                cc = colorbar;
                ytklabs = {[0.01 0.1 1 10 100]}; % {[0.001 0.01 0.1 1 10 100 1000]}; % {'-
                set(cc,'ytick',cax_lo:cax_hi,'yticklabel',ytklabs); % -10.^(-2:0.5:2))
                cm = colormap(flipud(hot(512)));
                cm = cm(80:450,:);
                colormap(cm)
                % title({'log sensitivity latitude', 'to transv magncl param'})
                caxis([cax_lo cax_hi])
                % title('(e)')
                set(gca,'FontSize',FntSz)
                axis([-r_max_alf r_max_alf 0 r_max_lat])
                hold
            end
            
            if ismember(9,plot_choices)
                if plot_for_ms
                    figure('Position',[left bottom width height]); %
                else
                    figure
                end
                   imagesc(-max_alf:1/resln:max_alf,all_lats,dlon_alf_lox) % 
                % contourf(xx,yy,log10(-dalf_alf_tcs)) % .*sign(dthdg))
%                 imagesc(-max_alf:1/resln:max_alf,all_lats,log(1-dlon_alf_lox)) % .*sign(dthdg))
                set(gca,'YDir','normal')
                % xlabel('transverse projection inclination (^o)')
                xlabel('West                   Heading                    East')
                ylabel('Latitude (^o)')
                set(gca,'YTick',0:30:90); 
%                 set(gca,'XTick',-90:30:90,'XTickLabel',{'90^o','60^o','30^o','0^o','-30^o','-60^o','-90^o'})
                set(gca,'XTick',0:30:90,'XTickLabel',{'0^o','30^o','60^o','90^o'}) % 
                
            %     colormap(flipud(summer))
                % title({'log sensitivity latitude', 'to transv magncl param'})

            %         cm = colormap(flipud(hot(512)));
            %     cm = cm(80:350,:);
            %     
                colormap(brewermap([],'YlOrRd')); %'Reds')); % cm = colormap(flipud(hot(128)));
%                 cm = cm(30:99,:);
%                 colormap(cm)
            %     colormap(flipud(autumn))
                % title({'log sensitivity latitude', 'to transv magncl param'})
            %     caxis([0 10])
                % title('(e)')
                set(gca,'FontSize',FntSz)    
                axis([0 90 0 90])
                hold
               caxis([0 log(13)]) % log10(101)
                
                if  ~plot_for_ms
                    
                     xlabel('West                   Heading                    East')

                    
                else
                    
                    xlabel('East or West Heading') %  xlabel('West     Heading     East')
%                     xtickangle(40)
                    
                    % new figure for colorbar
                    figure('Position',[left bottom width height]); %
                    caxis([0 log(13)]) % log10(101)
                    colormap(brewermap([],'YlOrRd')); %'Reds')); %
                    axis('off')
                end
                cb_all = colorbar('Location','North');
%                 cb_all = colorbar;
                title(cb_all,{'Sensitivity (^o / ^o)'}) %{'Sensitivity', 'Longitude ^o / ^o'})
            %        ytklabs = {[-0.01 -0.1 -1 -10 -100]}; 
            %     set(cc,'ytick',-2:2,'yticklabel',ytklabs); %
                ytklabs = {[0 1 2.5 5 10 20 50]}; %  100
                set(cb_all,'ytick',[0 log(2) log(3.5) log(6) log(11) log(21) log(51) ],'yticklabel',ytklabs); %  log10(101)          
                                
            end
            

            if ismember(12,plot_choices)
                if plot_for_ms
                    figure('Position',[left bottom width height]); %
                else
                    figure
                end
                   imagesc(-max_alf:1/resln:max_alf,all_lats,zeros(size(dlon_alf_lox)) + randn(size(dlon_alf_lox))/1000) % 
                % contourf(xx,yy,log10(-dalf_alf_tcs)) % .*sign(dthdg))
                imagesc(-max_alf:1/resln:max_alf,all_lats,log(1-dlon_alf_lox)) % .*sign(dthdg))
                set(gca,'YDir','normal')
                % xlabel('transverse projection inclination (^o)')
                xlabel('West                   Heading                    East')
                ylabel('Latitude (^o)')
                set(gca,'YTick',0:30:90); 
%                 set(gca,'XTick',-90:30:90,'XTickLabel',{'90^o','60^o','30^o','0^o','-30^o','-60^o','-90^o'})
                set(gca,'XTick',0:30:90,'XTickLabel',{'0^o','30^o','60^o','90^o'}) % 
                
            %     colormap(flipud(summer))
                % title({'log sensitivity latitude', 'to transv magncl param'})

            %         cm = colormap(flipud(hot(512)));
            %     cm = cm(80:350,:);
            %     
                colormap(brewermap([],'YlOrRd')); %'Reds')); % cm = colormap(flipud(hot(128)));
%                 cm = cm(30:99,:);
%                 colormap(cm)
            %     colormap(flipud(autumn))
                % title({'log sensitivity latitude', 'to transv magncl param'})
            %     caxis([0 10])
                % title('(e)')
                set(gca,'FontSize',FntSz)    
                axis([0 90 0 90])
                hold
               caxis([0 log(13)]) % log10(101)
                
                if  ~plot_for_ms
                    
                     xlabel('West                   Heading                    East')

                    
                else
                    
                    xlabel('East or West Heading') %  xlabel('West     Heading     East')
%                     xtickangle(40)
                    
                    % new figure for colorbar
                    figure('Position',[left bottom width height]); %
                    caxis([0 log(13)]) % log10(101)
                    colormap(brewermap([],'YlOrRd')); %'Reds')); %
                    axis('off')
                end
                cb_all = colorbar('Location','North');
%                 cb_all = colorbar;
                title(cb_all,{'Sensitivity (^o / ^o)'}) %{'Sensitivity', 'Longitude ^o / ^o'})
            %        ytklabs = {[-0.01 -0.1 -1 -10 -100]}; 
            %     set(cc,'ytick',-2:2,'yticklabel',ytklabs); %
                ytklabs = {[0 1 2.5 5 10 20 50]}; %  100
                set(cb_all,'ytick',[0 log(2) log(3.5) log(6) log(11) log(21) log(51) ],'yticklabel',ytklabs); %  log10(101)          
                                
            end
            
            if ismember(11,plot_choices)

                if plot_for_ms
                    figure('Position',[left bottom width height]); %
                else
                    figure
                end
                % contourf(xx,yy,log10(-dalf_alf_tcs)) % .*sign(dthdg))
                dalf_tcs_both = dalf_tcs_1 + dalf_tcs_2;
                imagesc(-max_alf:1/resln:max_alf,all_lat_suns,sign(dalf_tcs_both).*log(1+abs(dalf_tcs_both))) %dalf_tcs_both) %  .*sign(dthdg))
                set(gca,'YDir','normal')
                set(gca, 'XDir','reverse')
                % xlabel('transverse projection inclination (^o)')
                if  ~plot_for_ms
                     xlabel('West                   Heading                    East')
                else
                    xlabel('West      Heading      East')
                    xtickangle(40)
                    a=gca;
                    a.XRuler.TickLabelGapOffset = -1.5;  
%                     set(gca,'TickDir','in'); % 'out');
                end
                ylabel('Latitude (^o)')
                set(gca,'YTick',Lat_sun_tiks); 
                set(gca,'XTick',-120:30:120,'XTickLabel',{'-120^o','-90^o','-60^o','-30^o','0^o','30^o','60^o','90^o','120^o'})
                ytklabs = {[-10 -5 -1 0 1 5 10]}; 
                if incl_cb
                    cc = colorbar;
                    set(cc,'ytick',[-log(11) -log(6) -log(2) 0 log(2) log(6) log(11)],'yticklabel',ytklabs); %

                    title(cc,{'Sensitivity ^o / ^o'})
                end
%                 set(cc,'ytick',-2:2);

%% ooo lway to combine colormaps but Brewer is better
%                 cm = colormap(flipud(hot(128)));
%                 cm = cm(30:79,:);
%                 cmap = [summer(50); cm]; % flipud(autumn(32))]; % flipud(hot(32))]; % 
%                 colormap(cmap)
                colormap(brewermap([],'*RdYlBu')) %,'RdBu')) %,'RdYlBu')) %,'PRGn')) % ,'RdYlBu')) %'PiYG')) %'RdBu')) %  'YlGnBu')) %
           
                caxis([-log(6)-0.5 log(11)+0.5]) %-10 10]) % 
    
%                 caxis([-log(11) log(11)])
                set(gca,'FontSize',FntSz)
            
%                 set(gca,'FontSize',FntSz)
                axis([-r_max_alf r_max_alf 0 last_sun_lat])
%                 hold

            end


    end
    
    % next plot sensitivities all solar factors for this date

    % use black for polar summer or winter
        cm_bl = [0 0 0];
        cm_cy = [0 1 1];
        
        cm_yl = [1 1 0.25]; % [0.95 0.95 0.6];
                   
        if ismember(2,plot_choices)
                if plot_for_ms
                    figure('Position',[left bottom width height]); %
                else
                    figure
                end
%             if date > 265
                
%                 dalf_alf_fix(isinf(dalf_alf_fix)) = Inf; 
                
%             else
%                 
%                 cmap = [cm_cy; summer(50); cm; cm_cy]; % flipud(autumn(32))]; % flipud(hot(32))]; % 
%                          
%             end
            imagesc(-max_alf:1/resln:max_alf,all_lat_suns,sign(dalf_alf_fix).*log(1+abs(dalf_alf_fix))) % dalf_alf_fix) % .*sign(dthdg))
            set(gca,'YDir','normal')
            set(gca, 'XDir','reverse')
            % contourf(xx,yy,dalf_alf_fix) % .*sign(dthdg))
            % xlabel('transverse projection inclination (^o)')
            if  ~plot_for_ms
                 xlabel('West                   Heading                    East')
            else
                xlabel('West      Heading      East')
                xtickangle(40)
                a=gca;
                a.XRuler.TickLabelGapOffset = -1.5; 
            end
            ylabel('Latitude (^o)')
            set(gca,'YTick',Lat_sun_tiks); 
            set(gca,'XTick',-120:30:120,'XTickLabel',{'-120^o','-90^o','-60^o','-30^o','0^o','30^o','60^o','90^o','120^o'})
             if incl_cb
                cc = colorbar;
                ytklabs = {'25','10','0','10','25'};
                set(cc,'ytick',[-log(1.25) -log(1.1) 0 log(1.1) log(1.25)],'yticklabel',ytklabs); %

                title(cc,{'Sensitivity ^%'})
            end          
 
        %     set(cc,'ytick',[-5 -2 0 2 5]); % -10.^(-2:0.5:2))
%             cm = colormap(flipud(hot(1024)));
%             cm = cm(240:632,:); % cm(30:79,:);
%             if date > 265
%                 
%                 cmap = [cm_bl; summer(400); cm; cm_bl]; % flipud(autumn(32))]; % flipud(hot(32))]; % 
%                 
%             else
%                 
%                 cmap = [cm_cy; summer(400); cm; cm_cy]; % flipud(autumn(32))]; % flipud(hot(32))]; % 
%                          
%             end
            cmap = colormap(brewermap(1024,'*RdYlBu')); %'*RdBu')); %'*RdYlBu')); %,'RdBu')) %,'RdYlBu')) %,'PRGn')) % ,'RdYlBu')) %'PiYG')) %'RdBu')) %  'YlGnBu')) %
            
            
%              if date < 265
%                  cmap = [cm_yl; cm; cm_yl]; %
%              else
%                  cmap = [cm_bl; cm; cm_bl]; %
%              end
            
            colormap(cmap) % flipud(parula))
   
            % title({'log sensitivity latitude', 'to transv magncl param'})
            caxis([-log(1.25) log(1.25)]) %caxis([-log(11) log(11)]) % ([-log(6)-0.5 log(6)+0.5]) %caxis([-5 5]) % 
            set(gca,'FontSize',FntSz)
            axis([-r_max_alf r_max_alf 0 last_sun_lat])
            hold
        end

        if ismember(3,plot_choices)
                if plot_for_ms
                    figure('Position',[left bottom width height]); %
                else
                    figure
                end
            % contourf(xx,yy,log10(-dalf_alf_tcs)) % .*sign(dthdg))
            imagesc(-max_alf:1/resln:max_alf,all_lat_suns,sign(dalf_alf_tcs).*log(1+abs(dalf_alf_tcs))) %dalf_alf_tcs) %  .*sign(dthdg))
            set(gca,'YDir','normal')
            set(gca, 'XDir','reverse')
            % xlabel('transverse projection inclination (^o)')
            if  ~plot_for_ms
                 xlabel('West                   Heading                    East')
            else
                xlabel('West      Heading      East')
                xtickangle(40)
                a=gca;
                a.XRuler.TickLabelGapOffset = -1.5; 
            end    

            ylabel('Latitude (^o)')
            set(gca,'YTick',Lat_sun_tiks); 
            set(gca,'XTick',-120:30:120,'XTickLabel',{'-120^o','-90^o','-60^o','-30^o','0^o','30^o','60^o','90^o','120^o'})
            if incl_cb
                cc = colorbar;
                ytklabs = {'-25','-10','0','10','25'};
                set(cc,'ytick',[-log(1.25) -log(1.1) 0 log(1.1) log(1.25)],'yticklabel',ytklabs); % log(11)

                title(cc,{'Sensitivity (%)'})
            end
        %     set(cc,'ytick',[-5 -2 0 2 5]); % -10.^(-2:0.5:2))
%             cm = colormap(flipud(hot(1024)));
%             cm = cm(240:632,:); % cm(30:79,:);
%             if date > 265
%                 
%                 cmap = [cm_bl; summer(400); cm; cm_bl]; % flipud(autumn(32))]; % flipud(hot(32))]; % 
%                 
%             else
%                 
%                 cmap = [cm_cy; summer(400); cm; cm_cy]; % flipud(autumn(32))]; % flipud(hot(32))]; % 
%                          
%             end
% %             cmap = [summer(60); cm]; % flipud(autumn(32))]; % flipud(hot(32))]; % 

            cmap = colormap(brewermap(1024,'*RdYlBu')); %'*RdYlBu')); %,'RdBu')) %,'RdYlBu')) %,'PRGn')) % ,'RdYlBu')) %'PiYG')) %'RdBu')) %  'YlGnBu')) %
            
            
%              if date < 265
%                  cmap = [cm_yl; cm; cm_yl]; %
%              else
%                  cmap = [cm_bl; cm; cm_bl]; %
%              end
            
            colormap(cmap) % flipud(parula))
            
            % title({'log sensitivity latitude', 'to transv magncl param'})
            caxis([-log(1.25) log(1.25)]) % ([-log(11) log(11)]) % ([-log(11)-0.5 log(11)+0.5]) %caxis([-10 10]) % 
            % title('(e)')
            set(gca,'FontSize',FntSz)
            axis([-r_max_alf r_max_alf 0 last_sun_lat])
            hold

        end

        %% plots vs previous offset initial heading alf_0



        if ismember(5,plot_choices)
                if plot_for_ms
                    figure('Position',[left bottom width height]); %
                else
                    figure
                end
            imagesc(-max_alf:1/resln:max_alf,all_lat_suns,dalf_alf0_sun) % .*sign(dthdg))
            set(gca,'YDir','normal')
            set(gca, 'XDir','reverse')
            % contourf(xx,yy,dalf_alf_fix) % .*sign(dthdg))
            % xlabel('transverse projection inclination (^o)')
            xlabel('West                   Heading                    East')
            ylabel('Latitude (^o)')
            set(gca,'YTick',0:30:90); 
            set(gca,'XTick',-120:30:120,'XTickLabel',{'-120^o','-90^o','-60^o','-30^o','0^o','30^o','60^o','90^o','120^o'})
            cc = colorbar;
            % ytklabs = {[-0.01 -0.1 -1 -10 -100]}; % {'-10^-^2','-10^-^1','-1','-10','-100'};
            % set(cc,'ytick',-2:2,'yticklabel',ytklabs); % -10.^(-2:0.5:2))
            caxis([-2 1.1])
            colormap(flipud(parula))
            % title({'log sensitivity latitude', 'to transv magncl param'})
            % caxis([cax_lo cax_hi])
            % title('(e)')
            set(gca,'FontSize',FntSz)
            axis([-r_max_alf r_max_alf 0 last_sun_lat])
            hold
        end


        %% plots vs perturbation latitude

        if ismember(7,plot_choices)
                if plot_for_ms
                    figure('Position',[left bottom width height]); %
                else
                    figure
                end
            imagesc(-max_alf:1/resln:max_alf,all_lat_suns,log10(abs(dalf_lat_fix))) % .*sign(dthdg))
            set(gca,'YDir','normal')
            set(gca, 'XDir','reverse')
            % contourf(xx,yy,dalf_alf_fix) % .*sign(dthdg))
            % xlabel('transverse projection inclination (^o)')
            xlabel('West                   Heading                    East')
            ylabel('Latitude (^o)')
             set(gca,'YTick',0:30:90); 
             set(gca,'XTick',-120:30:120,'XTickLabel',{'-120^o','-90^o','-60^o','-30^o','0^o','30^o','60^o','90^o','120^o'})
            cc = colorbar;
            ytklabs = {[-0.01 -0.1 -1 -10 -100]}; % {'-10^-^2','-10^-^1','-1','-10','-100'};
            set(cc,'ytick',-2:2,'yticklabel',ytklabs); % -10.^(-2:0.5:2))
            caxis([-2 2])
            colormap(flipud(parula))
            % title({'log sensitivity latitude', 'to transv magncl param'})
            % caxis([cax_lo cax_hi])
            % title('(e)')
            set(gca,'FontSize',FntSz)
            axis([-r_max_alf r_max_alf 0 last_sun_lat])
            hold
        end

        if ismember(8,plot_choices)
                if plot_for_ms
                    figure('Position',[left bottom width height]); %
                else
                    figure
                end
            % contourf(xx,yy,log10(-dalf_alf_tcs)) % .*sign(dthdg))
            imagesc(-max_alf:1/resln:max_alf,all_lats,log10(abs(dalf_lat_tcs))) % .*sign(dthdg))
            set(gca,'YDir','normal')
            set(gca, 'XDir','reverse')
            % xlabel('transverse projection inclination (^o)')
            xlabel('West                   Heading                    East')
            ylabel('Latitude (^o)')
             set(gca,'YTick',0:30:90); 
             set(gca,'XTick',-120:30:120,'XTickLabel',{'-120^o','-90^o','-60^o','-30^o','0^o','30^o','60^o','90^o','120^o'})
            cc = colorbar;
            ytklabs = {[-0.01 -0.1 -1 -10 -100]}; 
            set(cc,'ytick',-2:2,'yticklabel',ytklabs); % -10.^(-2:0.5:2))
            colormap(flipud(parula))
            % title({'log sensitivity latitude', 'to transv magncl param'})
            caxis([-2 2])
            % title('(e)')
            set(gca,'FontSize',FntSz)
            axis([-r_max_alf r_max_alf 0 r_max_lat])
            hold
        end
            
        if ismember(10,plot_choices)
                if plot_for_ms
                    figure('Position',[left bottom width height]); %
                else
                    figure
                end
            % contourf(xx,yy,log10(-dalf_alf_tcs)) % .*sign(dthdg))
            imagesc(-max_alf:1/resln:max_alf,all_lats,dlon_lon_tcs) % .*sign(dthdg))
            set(gca,'YDir','normal')
            set(gca, 'XDir','reverse')
            % xlabel('transverse projection inclination (^o)')
            xlabel('West                   Heading                    East')
            ylabel('Latitude (^o)')
             set(gca,'YTick',0:30:90); 
             set(gca,'XTick',-120:30:120,'XTickLabel',{'-120^o','-90^o','-60^o','-30^o','0^o','30^o','60^o','90^o','120^o'})
            cc = colorbar;
        %     ytklabs = {[-0.01 -0.1 -1 -10 -100]}; 
        %     set(cc,'ytick',[-10 -5 -2 -1 0 1]); % -10.^(-2:0.5:2))
        %     ytklabs = {[-10 -5 -1 0]}; 
        %     set(cc,'ytick',[-log(11) -log(6) -log(2) 0],'yticklabel',ytklabs); % -10.^(-2:0.5:2))
        %     colormap(summer) % flipud(parula))
        %     % title({'log sensitivity latitude', 'to transv magncl param'})
        %     caxis([-log(11) 0])
            caxis([-2 0])
            % title('(e)')
            set(gca,'FontSize',FntSz)
            axis([-r_max_alf r_max_alf 0 r_max_lat])
            hold
        end
        
end

if ismember(2,plot_choices) || ismember(3,plot_choices) || ismember(11,plot_choices)
    
    if plot_for_ms
        figure('Position',[left bottom*0.75 width 2.15*height]); %
    else
        figure
    end
    cc = colorbar; axis('off')
   
    colormap(brewermap([],'*RdYlBu')) 
    set(cc,'ytick',[-log(1.25) -log(1.1)  0 log(1.1) log(1.25)],'yticklabel',[-25 -10 0 10 25]); %
     caxis([-log(1.25) log(1.25)])
    set(cc,'FontSize',FntSz-1)
    title(cc,{'Sensitivity (%)'})
    
end
