
%%  revisions Comm Bio - plot range of gain with compass precision among 
addpath 'brewer'

if ~exist('rte')
    load regr_output_biol
end

% plot gain per species and compass precision (detection error)
n_errs = numel(det_err_regr);
n_e_m1 = n_errs - 1;

for isp = 1:n_spec
    
    all_gain_TCSC(isp,1:n_errs) = 100*(rte{3}.Y(isp:n_spec:isp+n_e_m1*n_spec)./ ...
        rte{1}.Y(isp:n_spec:isp+n_e_m1*n_spec)-1);
    all_gain_TCSC_GM(isp,1:n_errs) = 100*(rte{3}.Y(isp:n_spec:isp+n_e_m1*n_spec)./ ...
        rte{2}.Y(isp:n_spec:isp+n_e_m1*n_spec)-1);
    
    med_gain_TCSC(isp) = median(100*(rte{3}.Y(isp:n_spec:end)./ ...
        rte{1}.Y(isp:n_spec:end)-1));
    
    med_gain_TCSC_GM(isp) = median(100*(rte{3}.Y(isp:n_spec:end)./ ...
        rte{2}.Y(isp:n_spec:end)-1));
    
end

noc_bs = 1:9; %  [2:3 5 6:9];
[gain_fact,i_gns] = sort(log(day_m_d.*n_hat_fls.*geo_fact_gc(1:n_spec)')); 
% gain_fact = log(day_m_d.*n_hat_fls.*geo_fact_gc(1:n_spec)'); % log(ff(noc_bs));

DndSz2 = 210; % 180; % 165;
StrSz2 = 150; % 85;

figure('Position',[300 300 225 340]);
hold

% now add 1st & 2nd species (offset for visability)
offst = 0.5;
% off1 = 2; % 1; % 
% off2 = 3; % 2; %

% hD = scatter(det_err_regr-offst,gain_fact(off1)*ones(1,n_errs),DndSz2, ...
%     all_gain_TCSC_GM(noc_bs(off1),:),'d','fill');
% hD.MarkerEdgeColor = 'k'; 
% 
% hS = scatter(det_err_regr-offst,gain_fact(off1)*ones(1,n_errs),StrSz2, ...
%     all_gain_TCSC(noc_bs(off1),:),'h','fill');
% hS.MarkerEdgeColor =  'k'; % 'w'; %
% 
% hD = scatter(det_err_regr+offst,gain_fact(off2)*ones(1,n_errs),DndSz2, ...
%     all_gain_TCSC_GM(noc_bs(off2),:),'d','fill');
% hD.MarkerEdgeColor = 'k'; 
% 
% hS = scatter(det_err_regr+offst,gain_fact(off2)*ones(1,n_errs),StrSz2, ...
%     all_gain_TCSC(noc_bs(off2),:),'h','fill');
% hS.MarkerEdgeColor =  'k'; % 'w'; %
% 

for i_spec = [9 8 6 5 2] % 9:-1:1 %  3:7
    
    i_sp = i_gns(i_spec);
    
    offs = (i_sp == 1 || i_sp == 3)*offst + ...
        - (i_sp == 2 || i_sp == 4)*offst; 

    hD = scatter(det_err_regr+offs,gain_fact(i_spec)*ones(1,n_errs),DndSz2, ...
        all_gain_TCSC_GM(i_sp,:),'d','fill');
    hD.MarkerEdgeColor = 'k'; 
    
    hS = scatter(det_err_regr+offs,gain_fact(i_spec)*ones(1,n_errs),StrSz2, ...
        all_gain_TCSC(i_sp,:),'h','fill');
    hS.MarkerEdgeColor = 'k'; %'w'; %  
    
end

colormap(brewermap([],'*RdYlBu')) % Y
set(gca,'YTick',gain_fact,'YTickLabel',[]) %,'TickDir','out')
set(gca,'YAxisLocation','right','FontSize',10)
xlabel('Compass precision (^o)','FontSize',10)
ylabel('                  Product gain factors','FontSize',10)
xlim([0 34])
ylim([7.35 9.75])
caxis([-40 40])

xh = get(gca,'xlabel') % handle to the label object
p = get(xh,'position') % get the current position property
p(2) = 7.2;        % double the distance, 
                     % negative values put the label below the axis
set(xh,'position',p)  

%% colorbar
figure('Position',[400 200 215 160]); %175]); % ,c
% rotate tick angle colorbar :-)
cb = colorbar('North');
title(cb,{'Gain TCSC (%)'},'FontSize',9) % 
set(cb,'XTick',-40:20:40,'FontSize',9)
colormap(brewermap([],'*RdYlBu')) % 'YlOrRd'))
caxis([-40 40]) % ([-35 35]) % 
axis off

keepLabels = ~cellfun(@isempty,cb.TickLabels);
cb.TickLabels = [];
ax = axes('Position', cb.Position,...
    'Color', 'none',...
    'YTick', [],...
    'XLim', cb.Limits,...
    'XTick', cb.Ticks(keepLabels),...
    'FontSize', cb.FontSize);

xtickangle(ax,0);




         