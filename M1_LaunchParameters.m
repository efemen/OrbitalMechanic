% Launch will take place from Foca Izmir, Turkey. On 20th of October 2024.

H = 500; % km, parking orbit altitude.
mu = 398600;

lat = 38.741383;
lon = 26.740603;


jdt = juliandate([2024, 10, 20, 18, 0, 0]);
sdt = siderealTime(jdt) + lon;

% i = lat; % Launch azimuth = 90 deg.
R_E = 6371 + H; % km
A = 90;
rho = [sind(A), cosd(A), 0];
V_hat = t_xX(lat, sdt) * rho';

R_i = [R_E * sind(90 - lat) * cosd(sdt), R_E * sind(90 - lat) * sind(sdt), R_E * cosd(90 - lat)];

magV = sqrt(mu / norm(R_i));
V_i = (magV * V_hat)';


% clear i R_s H f lon R_c R_E R_phi V R_hat V_hat e rho A H