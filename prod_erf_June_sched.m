function p_Arr = prod_erf_June_sched(b,X_plus)

try
% parameters and default values in first two rows
% ndft = sum(X_plus(1,:)==-1);
idx_fit = X_plus(:,end-1) == 1;
idx_dft = X_plus(:,end-1) == -1;
bs(idx_fit) = b;
bs(idx_dft) = X_plus(idx_dft,end);

% predictor fields after first two rows

fr_br =  X_plus(:,1);
N_0 =  X_plus(:,2);
N_max = X_plus(:,3);
sig_alf =  X_plus(:,4);
Geo =  X_plus(:,5);
cos_alf = X_plus(:,6);
cos_2alf = X_plus(:,7);
cos_lat_Arr_2 = X_plus(:,8);
inh_err_s = X_plus(:,9);
sched_err_s = X_plus(:,10);
secant_fact = X_plus(:,11);
% day mig d: if not sun comp = ones
day_mig_d_fact = X_plus(:,12);
% final headings for effective breadth
goal_head = X_plus(:,13);

% sin_alf = sqrt(1-cos_alf.^2);
% sin_alf_2 = 1-cos_alf.^2;

% geo factor breadth (note opposite in polarity to geo fact N_hat - 
% cos and sin terms 'switched')
fr_br_geo = fr_br.*sqrt(cos(goal_head).^2./cos_lat_Arr_2 + sin(goal_head).^2);
% fr_br_geo = fr_br.*sqrt(cos_alf.^2.*cos_lat_Arr_2 + sin_alf_2);
% fr_br_geo = fr_br.*sqrt(sin_alf_2.*cos_lat_Arr_2 + cos_alf.^2);

% exponents for assessing arr at goal Lat and within breadth (br)

e_N_Lat = bs(1);
% e_N_Lon = bs(2);
% here we use offset to 0.5
de_N_br = bs(2);

e_g_Lat = bs(3);
% e_g_Lon = bs(5);
e_g_br = bs(4);

e_sig_Lat = bs(5);
% e_sig_Lon = bs(8);
e_sig_br = bs(6);

% slopes of geo factor n-exp's
e_bn_Lat = bs(7);
% e_bn_Lon = bs(11);
e_bn_br = bs(8);

e_frbr = bs(9);
c_sunaz = bs(10)*pi/180;
e_secnt = bs(11);

% exponent for day mig dist effect on error reduction with N steps
e_m_dist = bs(12);

% slopes of geo factor g-exp's
% e_bg_Lat = bs(7);
% e_bg_br = bs(8);

% quantify non-precision compass variability, i.e., compass bias, 
% from inheritance and variability in schedules
sig_fix_2 = inh_err_s.^2 + (c_sunaz*sched_err_s.*secant_fact.^e_secnt).^2 ...
    + (c_sunaz<0)*1e7;

sig_adj =  sqrt(sig_fix_2 + sig_alf.^2);

% first determine Geo factors using slope and constant
Geo_Lat = Geo.^e_g_Lat; %.*log(1+e_bg_Lat.*sig_alf.^2)
% Geo_Lon= Geo.^e_g_Lon; %.*log(1+e_bg_Lat.*sig_alf.^2)
Geo_br = Geo.^e_g_br; %.*log(1+e_bg_br.*sig_alf.^2)

% estimate kappas for each measure incl possible Geo effects
kap_Lat = 1./(sig_adj.*Geo_Lat).^2;
% kap_Lon = 1./(sig_adj.*Geo_Lon).^2;

% effct_sig = (0.5*sig_alf.*(1-bes_2_br)); % sig_S1
kap_br = 1./(sig_adj.*Geo_br).^2;

bes_0_br = besseli(0,kap_br,1);
bes_1_br = besseli(1,kap_br,1)./bes_0_br;
bes_2_br = besseli(2,kap_br,1)./bes_0_br;

% kap_Lat = 1./(sig_alf).^2;
% kap_br = 1./(sig_alf).^2;

% ... and Bessel fn. ratios including Geo effects
bes_0_Lat = besseli(0,kap_Lat,1);
bes_1_Lat = besseli(1,kap_Lat,1)./bes_0_Lat;
bes_2_Lat = besseli(2,kap_Lat,1)./bes_0_Lat;
% bes_1_Lon = besseli(1,kap_Lon,1)./besseli(0,kap_Lon,1);


% and mean moments for std cos(heading)
% sig_C1 = 0.5*(1 + bes_2_Lat.*cos_2alf - 2*(bes_1_Lat.*cos_alf).^2);
% sig_S1 = 0.5*(1 - bes_2_Lat.*cos_2alf - 2*(bes_1_Lat.*sin_alf).^2);
sig_C1 = 0.5*(1 + bes_2_br.*cos_2alf - 2*(bes_1_br.*cos_alf).^2);
% sig_S1 = 0.5*(1 - bes_2_br.*cos_2alf - 2*(bes_1_br.*sin_alf).^2);

% resolver n-exps
% exp_N_Lat = abs(e_N_Lat).*(exp(-e_bn_Lat./kap_Lat));
% exp_N_br = abs(e_N_br).*(exp(-e_bn_br./kap_br));
exp_N_Lat = e_N_Lat.*(exp(-e_bn_Lat./kap_Lat));
% exp_N_Lon = e_N_Lon.*(exp(-e_bn_Lon./kap_Lon));

% N_hat_arr_Lat = N_0./bes_1_Lat;

% now estimate prob arrival at goal Lat
p_arrs_Lats = 0.5*(1 - erf(sqrt(0.5)* ...
    (N_0./N_max - bes_1_Lat).*cos_alf.* ...
    N_max.^exp_N_Lat./sig_C1.^e_sig_Lat));    

% p_arrs_Lons = 0.5*(1 - erf(sqrt(0.5)* ... % 
%     (N_0./N_max - bes_1_Lon).*sin_alf.* ...
%     N_max.^exp_N_Lon./sig_S1.^e_sig_Lon));    

% p_arrs_Lats = 1; % 
% p_arrs_Lons = 1; % 

% p_arrs_Lats = 0.5 - erf(sqrt(0.5).* ...
%     (N_0./N_hat_arr_Lat - bes_1_Lat).*cos_alf.* ...
%     N_hat_arr_Lat.^exp_N_Lat./sig_C1.^e_sig_Lat)/2; 

% ... and - on estinated arrival at goal - whether within required 
% migratory breadth

% now combine 'sigma' and N_hat with same exponent
% p_br = erf(sqrt(0.5)*fr_br_geo.*(kap_br.*N_0./bes_1_br).^e_N_br); % ./ ...

% expected number steps and step-adjusted ang dev
N_hat_arr = N_0./bes_1_br;
% exp_N_br = (0.5+de_N_br).*(exp(-e_bn_br./kap_br));
exp_N_br = (0.5+de_N_br.*day_mig_d_fact.^e_m_dist).*(exp(-e_bn_br./kap_br));
sig_N_stps = sig_alf./N_hat_arr.^exp_N_br;

p_br = erf(sqrt(0.5)*fr_br_geo.^e_frbr./ ...
    sqrt(Geo_br.^2.*(sig_fix_2 + sig_N_stps.^2)).^e_sig_br); % ./ ...

%     (sig_alf.*Geo_br).^e_sig_br); % 
    
% multiply two marginal probs
% p_Arr =  p_arrs_Lats.*p_arrs_Lons.*p_br; % p_br; % 
p_Arr =   p_arrs_Lats.*p_br; % mean(p_arrs_Lats,p_arrs_Lons).*p_br;

if any(isnan(p_Arr) | isinf(p_Arr))
    p_Arr(isnan(p_Arr)) = 0;
    p_Arr(isinf(p_Arr)) = 1;
end

catch
    
    keyboard
    
end