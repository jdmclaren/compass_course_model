function  [sun_az, cal_errs, is_sunset] = calc_sun_az(lats,dates,varargin)

del_max = 0.4091; % 23.44 degs in rads
% used for calc severity of polar night to augment cal error

twelve_degs = 0.2094; % use 12 degs as angle for nautical twilight 
% at this angle birds can potentially still determine
% the sun azimuth as a glow

if ~isempty(varargin)
    
    std_cal_err = varargin{1}; % in rads
    max_cal_err_light = varargin{2};
    max_cal_err_dark = varargin{3};

    if ~isnan(max_cal_err_light)
        max_d_err_light = max_cal_err_light - std_cal_err;
        max_d_err_dark = max_cal_err_dark - std_cal_err;
    else % mixed compass - assume uses magnetic compass 
        % (but then double the std error since algorithm later halves it)
        std_cal_err = std_cal_err*2;
        max_d_err_light = std_cal_err*0;
        max_d_err_dark = std_cal_err*0;            
    end
    
else % assume no error
    
    std_cal_err = 0*lats;
    max_d_err_light = 0*lats;
    max_d_err_dark = 0*lats;
    
end

% 360 - acosd(sind(-23.44* ...
%     cosd(360*255/360))/cosd(60))

c1 = 0.39779; % *pi/180;
c2 = 1.914*pi/180;
c3 = 0.98565*pi/180;

sin_decls = -c1*cos(c3*(dates+10) + ...
    c2*sin(c3*(dates-2)));

cos_lats = cos(lats);

is_sunset = abs(sin_decls) < cos_lats;
% no_sunset_s = ~is_sunset;

no_sunset = ~is_sunset & sin_decls >= 0;
no_sun = ~is_sunset & sin_decls < 0;

sun_az(is_sunset,1) = acos( ...
    -sin_decls(is_sunset)./cos_lats(is_sunset));

% error normal when sunset exists
cal_errs(is_sunset,1) = std_cal_err(is_sunset,1);

if sum(~is_sunset) > 0
    
    sun_az(no_sunset,1) = pi; % pi/2; 
    
    % sign of cal errs used below
    sign_errs = sign(std_cal_err);
    
    % pi/2's to explore how using N/S could work in polar summer/winter

    sun_az(no_sun,1) =  0; % pi/2; %
    
    % calculate error magnification relative to max when 
    % sun never reaches nautical twilight (12 degrees below horizon)
    elev_no_sun = asin(cos(lats(no_sun,1)-asin(sin_decls(no_sun,1))));
    if ~isempty(elev_no_sun)
        % -ve sign since elev -ve
           cal_errs(no_sun,1) = std_cal_err(no_sun,1) + ...
            sign_errs(no_sun).*min(elev_no_sun/twelve_degs,1).*max_d_err_dark(no_sun,1);
    end
    % calcultae error magnification when sun never sets relative to direct
    % overhead sun, i.e., elev = 90
    % -ve sign here comes from min angle being at sunset (vs. above 
    % formula polar night)
    elev_no_sunset = -asin(cos(lats(no_sunset,1)+asin(sin_decls(no_sunset,1))));
    if ~isempty(elev_no_sunset)
          cal_errs(no_sunset,1) = std_cal_err(no_sunset,1) + ...
            sign_errs(no_sunset).*min(elev_no_sunset/del_max,1).*max_d_err_light(no_sunset,1);    
    end
end
