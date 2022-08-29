lat_deps = [45 65]*pi/180;
lat_arrs = [25 0]*pi/180;

lin_stys = {'-','--'};

for i_lat = 1:2
    
    ln_fact_gen(i_lat) = (log(sec(lat_deps(i_lat)) + ...
             tan(lat_deps(i_lat))) - log(sec(lat_arrs(i_lat)) + ...
             tan(abs(lat_arrs(i_lat)))))*cos(lat_arrs(i_lat))/ ...
             (lat_deps(i_lat) - lat_arrs(i_lat));

end

alfs = [0:135]*pi/180;
sin_alfs = sin(alfs);
cos_alfs = cos(alfs);
% geo_fact_lox = sqrt((sin_alf_lox_rep.*log_Fact).^2 + ...
%     cos_alf_lox_rep.^2);

H = figure('Position',[300 300 240 240]);
hold

for i_lat = 1:2

    geo_fact_gens(i_lat,:) = sqrt((sin_alfs*ln_fact_gen(i_lat)).^2 + ...
        cos_alfs.^2);

    plot(alfs*180/pi,geo_fact_gens(i_lat,:),'LineWidth',1.5, ...
        'LineStyle',lin_stys{i_lat})

end

xlabel('Heading vs. South (^o)')
ylabel('Geographic factor')

% stick to South headings alf < 90
xlim([0 90])
set(gca,'XTick',0:30:90,'YTick',1:.1:1.35)