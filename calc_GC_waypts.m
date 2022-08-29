function [latGC,lonGC,n_flts_arr] = calc_GC_waypts(lat_lons,day_fd,max_fls)

% calcs waypts for great circle between two locations on sphere
% all angles in radians
% day flight dist in arc rads

% formulas from Ed Wiliams online "formulary"
% https://edwilliams.org/avform.htm
% accessed Nov 10, 2020

lat1 = lat_lons(1);
lon1 = lat_lons(2);
lat2 = lat_lons(3);
lon2 = lat_lons(4);


% distance between (scalar) points lat/lon 1&2
d=2*asin(sqrt((sin((lat1-lat2)/2))^2 + ...
                 cos(lat1)*cos(lat2)*(sin((lon1-lon2)/2))^2));
  
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

A=sin((1-f)*d)/sin(d);
B=sin(f*d)/sin(d);
x = A*cos(lat1)*cos(lon1) +  B*cos(lat2)*cos(lon2);
y = A*cos(lat1)*sin(lon1) +  B*cos(lat2)*sin(lon2);
z = A*sin(lat1) +  B*sin(lat2);

% lat/lon of waypoints
latGC=atan2(z,sqrt(x.^2+y.^2));
lonGC=mod(atan2(y,x),2*pi);
    
    