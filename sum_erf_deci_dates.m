function p_Arr = sum_erf_deci_dates(b,X)

% X(:,1:2) are latitudinal and longitudinal breadth
% X(:,3:4) cos and sin heads
% X(:,5,:) Bessel ratios
% X(:,6:7) are sigma cos and sin heads
% X(:,8) are (log(Lat) ratio) Geo facts for sin(head)
% X(:,9:10) are N_min_fls and N_max_fls

% sub-step factor (to avoid leapfrogging the goal area between steps
substp_frc =  1; % 10; %
frc_substp = 1/substp_frc;

% initialize arrival probs
p_Arr = zeros(size(X,1),1);

one_sq2 = 1/sqrt(2);
N_s =  X(:,9:10);

% first, modulate the deviance factors for Lat and Lon erf functions,
% depending on options

% non-modulated options for Lat and Lon erf functions,
% Will be updated per species/std and step number depending on options
if numel(b) == 2
    dev_frc_Lat_1 = one_sq2*X(:,3)./(X(:,6).^b(1));
    dev_frc_Lon_1 = one_sq2*X(:,4)./(X(:,7).^b(1));    
elseif numel(b) < 4
    dev_frc_Lat_1 = one_sq2*X(:,3)./X(:,6);
    dev_frc_Lon_1 = one_sq2*X(:,4)./X(:,7);
else
    dev_frc_Lat_1 = one_sq2*X(:,3)./(X(:,6).^b(4));
    dev_frc_Lon_1 = one_sq2*X(:,4)./(X(:,7).^b(4));
end    

% other factors not modulated per species / deviance combo
Geo_fact = X(:,8);
Bess_rat = X(:,5);
Brdth_Lat = X(:,1);
Brdth_Lon = X(:,2);

sig_alfs = X(:,11);
betas = X(:,12);

% p_Arr = zeros(size(X,2));

for i_rep = 1:size(X,1)
    
    % initial step number and number of (add'l) steps
    N_0 = N_s(i_rep,1);
    N_substs = substp_frc*(N_s(i_rep,2)-N_0)+1;
    
    % keep track of fractional steps
    frc_Ns = N_0:frc_substp:N_s(i_rep,2);
    
    % and number of full steps
    full_Ns = floor(N_0:frc_substp:N_s(i_rep,2));
    
    G_i = Geo_fact(i_rep);

    iN = 1;
    p_Lon = -1;
    p_Lat = -1;
    
    % prob beyond goal latitude (ignore relatifvely 
    % small chance or 'reverse' steps
    p_not_out = 1;
%     p_not_out(2) = 99;
    p_stpw_arr = [];
    
%     p_Lat_sum = 0;
    
    try 
    
    while iN <= N_substs && p_not_out(max(iN-1,1)) > 0 %  && (1 - p_Arr(i_rep,1)) > 1e-10 && p_Lon*p_Lat > 0
        
        N_i = full_Ns(iN);
%         sqN_i = sqrt(N_i);
        rat_fr_N_i = N_0/frc_Ns(iN);
        
        % first, modulate any variables, depending on options
        switch numel(b)

        case 0 % no modulation beyond sqrt(N) from multistep std err in denominator

            dev_frc_Lat = dev_frc_Lat_1*sqrt(N_i);
            dev_frc_Lon = dev_frc_Lon_1*sqrt(N_i);
            
        case 1 % geo fact to sin(heads)
                                        
            dev_frc_Lat = dev_frc_Lat_1*sqrt(N_i);
            dev_frc_Lon = (dev_frc_Lon_1*sqrt(N_i))*G_i^b(1);

%             dev_frc_Lat = dev_frc_Lat_1*N_i^b(1);
%             dev_frc_Lon = dev_frc_Lon_1*N_i^b(1);

%             dev_frc_Lat = b(1)*dev_frc_Lat_1;
%             dev_frc_Lon = b(1)*dev_frc_Lon_1;

% 
        case 2 % variable exponent for stepwise dependence sigmas 


%         case 3 % geo fact and sigmas
            
            dev_frc_Lat = dev_frc_Lat_1*N_i^b(2);
            dev_frc_Lon = dev_frc_Lon_1*N_i^b(2);
%             dev_frc_Lat = dev_frc_Lat_1*N_i^b(1);
%             dev_frc_Lon = dev_frc_Lon_1*G_i^b(2)*N_i^b(1);
            
        case 3 % geo fact to sin(heads)
                                        
%             dev_frc_Lat = dev_frc_Lat_1*sqrt(N_i);
%             dev_frc_Lon = (dev_frc_Lon_1*sqrt(N_i))*G_i^b(1);

            dev_frc_Lat = dev_frc_Lat_1*N_i^b(1);
            dev_frc_Lon = dev_frc_Lon_1*N_i^b(2)*G_i^b(3);

%             dev_frc_Lat = b(1)*dev_frc_Lat_1*N_i^b(2);
%             dev_frc_Lon = b(1)*dev_frc_Lon_1*N_i^b(2)*G_i^b(3);
% 
        case 4 % variable exponent for stepwise dependence sigmas 

%         case 3 % geo fact and sigmas
            
            dev_frc_Lat = dev_frc_Lat_1*N_i^b(1);
            dev_frc_Lon = dev_frc_Lon_1*G_i^b(3)*N_i^b(2);

        end

        % next, determine marginal arrival probabilities within
        % goal latitude and longitude bands

        erf_Lat_arg_1 = rat_fr_N_i*(1-Brdth_Lat(i_rep))-Bess_rat(i_rep);
        erf_Lat_arg_2 = rat_fr_N_i*(1+Brdth_Lat(i_rep))-Bess_rat(i_rep);
        
        % chance beyond goal radius of arr Lat
        p_not_out(iN) = 0.5*(1+erf(dev_frc_Lat(i_rep)*(erf_Lat_arg_2)));
        
        % chance Lat at least immeditaley below goal radius
        p_Lat = 0.5*(1-erf(dev_frc_Lat(i_rep)*(erf_Lat_arg_1)));
                
        % chance within Lon band
        erf_Lon_arg_1 = rat_fr_N_i*(1-Brdth_Lon(i_rep))-Bess_rat(i_rep);
        erf_Lon_arg_2 = rat_fr_N_i*(1+Brdth_Lon(i_rep))-Bess_rat(i_rep);
        p_Lon = 0.25*(1-erf(dev_frc_Lon(i_rep).*(erf_Lon_arg_1))).* ...
            (1+erf(dev_frc_Lon(i_rep)*(erf_Lon_arg_2)));
        

%         p_Lat_sum = p_Lat_sum +  (1-p_Lat_2_rads)*p_Lat;
%         p_Lat_sum = p_Lat_sum +  (1-p_Lat_sum)*p_Lat;
%         p_Arr_1 = p_Arr(i_rep,1) + (1-p_Arr(i_rep,1))* ...
%                        p_Lat.*p_Lon;

%         p_Arr(i_rep,1) = p_Arr(i_rep,1) + (1-p_Lat_2_rads)* ...
%                        p_Lat.*p_Lon;
                   
%         p_Arr(i_rep,1) = p_Arr(i_rep,1) + (1-p_Lat_sum)* ...
%                        p_Lat.*p_Lon;

        p_stpw_arr(iN) = (1-p_Arr(i_rep,1))* ...
                    prod(p_not_out)*p_Lat.*p_Lon;

        p_Arr(i_rep,1) = p_Arr(i_rep,1) +  p_stpw_arr(iN);
                                                  
       iN = iN + 1;
                   
    end
    
    p_Lon_2 = erf(betas(i_rep)*sqrt(N_i)/sig_alfs(i_rep)/sqrt(2));
    
        catch
    
    keyboard
    

    end
    
    try
    p_med(i_rep) = N_0 - 1+ find(p_stpw_arr == ...
        max(p_stpw_arr),1,'first');
    
    p_med_hat(i_rep) = N_0./(1-sig_alfs(i_rep)^2/2);
    catch
        keyboard
    end
    
end

    if any(isnan(p_Arr) | isinf(p_Arr))
        p_Arr(isnan(p_Arr)) = 0;
        p_Arr(isinf(p_Arr)) = 1;
    end
