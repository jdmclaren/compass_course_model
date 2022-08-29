clear 
d_ang = 50; %  1; % 
max_ang = 10000; % 100; % 
angs_ds = [1:100 100:1000:100000]; % 1:d_ang:max_ang;
angs = angs_ds*pi/180;

rt_12 = 2*sqrt(3);
rt_8 = 2*sqrt(2);
rt_2 = sqrt(2);
rt_pi = sqrt(pi);
pi_m1 = 1/pi;
pi_d3 = pi/3;
sq_pi_fct = sqrt(pi/sqrt(2));

for ia = 1:numel(angs)
    
    xx = vmrand(0,1/angs(ia)^2,[100000 1]);
    [s(ia), s0(ia)] = circ_std(xx);
    
end

last_idx = find(angs_ds>=100,1,'first');
idx = 1:last_idx;

figure
plot(angs_ds(idx),s(idx),'LineWidth',1)
hold
plot(angs_ds(idx),s0(idx),'LineWidth',1)
% plot(angs_ds,1./(1/12+angs.^-rt_8/3+angs.^-rt_2/3),'--g');
 plot(angs_ds(idx),(rt_12*angs(idx).^sq_pi_fct./(angs(idx).^sq_pi_fct+rt_pi)),'--g','LineWidth',1.5)
legend({'ang dev squared','circ std dev squared','best quotient fit'}, ...
    'AutoUpdate','off')

plot(angs_ds(idx),angs(idx),'--k','LineWidth',1)
ylim([0 2.2])

figure
plot(angs_ds(idx),s(idx),'LineWidth',1.)
hold
plot(angs_ds(idx),s0(idx),'LineWidth',1.)
plot(angs_ds(idx),rt_12*(erf(pi_m1*angs(idx).^pi_d3)),'--g','LineWidth',1.5);
legend({'ang dev squared','circ std dev squared','best erf fit'}, ...
    'AutoUpdate','off')
plot(angs_ds(idx),angs(idx),'--k')
ylim([0 1.2*rt_12])

figure
plot(angs_ds(idx),s0(idx)*180/pi./angs_ds(idx)*100-100)
title(' % diff s and s0')

% figure
% hold
% plot(angs_ds(idx), ...
%     besseli(1,1./angs(idx).^2)./besseli(0,1./angs(idx).^2), ...
%     'b','LineWidth',1.)
% eq_sigs = ([5 10 20 30 40 50 60]);
% eq_kaps = 1./(pi/180*eq_sigs).^2;
% scatter(eq_sigs,  ...
%     besseli(1,eq_kaps)./besseli(0,eq_kaps), ...
%     50,'bo','LineWidth',1.)

% approx b1 = 1.5, b2 = sqrt(pi)
mdl_crc_ang = @(b,X) rt_12*(X/pi).^b(1)./((X/pi).^b(1)+b(2));
fitnlm(angs,s0,mdl_crc_ang,[1 1 ])

% approx b1 = 1/pi, b2 = 1
mdl_crc_erf = @(b,X) rt_12*(erf(b(1)*X.^b(2)));
fitnlm(angs,s0,mdl_crc_erf,[1 1 ])

% approx b1 = 1.5, b2 = sqrt(pi)
mdl_crc_ang = @(b,X) rt_12*X.^b(1)./(X.^b(1)+b(2));
fitnlm(angs(idx),s0(idx),mdl_crc_ang,[1 1 ])

% approx b1 = 1/pi, b2 = pi/3
mdl_crc_erf = @(b,X) rt_12*(erf(b(1)*X.^b(2)));
fitnlm(angs(idx),s0(idx),mdl_crc_erf,[1 1 ])