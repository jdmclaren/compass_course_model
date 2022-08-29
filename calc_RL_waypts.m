function [latRL,lonRL,n_flts_arr] = calc_RL_waypts(lat_lons,day_fd,max_fls)

% calcs waypts for rhumblines between two locations on sphere 
% all angles in radians
% day flight dist in arc rads

% formulas from Ed Wiliams online "formulary"
% https://edwilliams.org/avform.htm
% accessed Nov 10, 2020

lat1 = lat_lons(1);
lon1 = lat_lons(2);
lat2 = lat_lons(3);
lon2 = lat_lons(4);

% calc rhumbline dist between pts

TOL = 1e-15;
dlon_W=mod(lon2-lon1,2*pi);
dlon_E=mod(lon1-lon2,2*pi);
dphi=log(tan(lat2/2+pi/4)/tan(lat1/2+pi/4));
if (abs(lat2-lat1) < sqrt(TOL))
 q=cos(lat1);
else 
 q= (lat2-lat1)/dphi;
end

if (dlon_W < dlon_E) % Westerly rhumb line is the shortest
  tc=mod(atan2(-dlon_W,dphi),2*pi);
  d= sqrt(q^2*dlon_W^2 + (lat2-lat1)^2);
else
  tc=mod(atan2(dlon_E,dphi),2*pi);
  d= sqrt(q^2*dlon_E^2 + (lat2-lat1)^2);
end

% calc true course between points
tc= mod(atan2(lon1-lon2,log(tan(lat2/2+pi/4)/tan(lat1/2+pi/4))),2*pi);    

% calc fractional distances along sphere for each day
f = 1- (d - (0:max_fls)*day_fd)/d;    

% calc numnber of flights to arrive (or close to it)
n_flts_arr = find(f >=1, 1,'first');
if isempty(n_flts_arr)
    n_flts_arr = max_fls;
else
    n_flts_arr = n_flts_arr-1;   
end
f = f(1:n_flts_arr);

% now calc locs for each day

latRL = lat1+f*d*cos(tc);

for i_d = 1:n_flts_arr
    if (abs(latRL(i_d)-lat1) < sqrt(TOL))
     q=cos(lat1);
    else
     dphi=log(tan(latRL(i_d)/2+pi/4)/tan(lat1/2+pi/4));
     q= (latRL(i_d)-lat1)/dphi;
    end
    dlon =-f(i_d)*d*sin(tc)/q;
    lonRL(i_d)=mod(lon1+dlon,2*pi);
end