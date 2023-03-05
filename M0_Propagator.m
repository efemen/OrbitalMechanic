clear; clc; close all;
%% Constants

mu = 398600; % Earth, km^3/s^2 
mu_Sun = 1.327124400e11; % Sun, km^3/s^2
earth_w = 0.004178074; % deg/s
m_Sun = 1.989E30;
m_Earth = 5.97219E24;
R_SOI = 150e6 * (m_Earth/m_Sun) ^ (2/5);

%% Setup Geometry and Plots
[X,Y,Z] = sphere;
R_E = 6371; % Earth Radius, km

X_E = X * R_E;
Y_E = Y * R_E;
Z_E = Z * R_E;
c_RotX = mean(mean(X_E));
c_RotY = mean(mean(Y_E));
c_Rot = [c_RotX c_RotY 0];



figure(1);
earth = surf(X_E,Y_E,-Z_E);
earthMap = imread("world_Map.jpg");



set(earth,'CData', earthMap,'FaceColor','texturemap',"EdgeColor","none")

hold on
colormap white
axis equal
set(gca,'Color','#BEBEBE');
set(gca, 'GridColor', 'white'); 


view(135,40)
run('M1_LaunchParameters.m')   % Get launch parameters
rotate(earth, [0 0 1], siderealTime(jdt))

quiver3(6378,0,0, 4000, 0, 0,"filled","LineWidth",3,"ShowArrowHead","on", "Color","green","MaxHeadSize",10);
text(6378 * 2,0,0,"Vernal Eq. ♈")

%% Earth Ephemeris
[earth_Pos, earth_V] = planetEphemeris(jdt, "SolarSystem","Earth");

earth_V_hat = earth_V / norm(earth_V);
rootX = earth_V_hat * 6378;
tipX = earth_V_hat * 4000;
quiver3(rootX(1),rootX(2), 0, tipX(1), tipX(2), 0,"filled","LineWidth",3,"ShowArrowHead","on", "Color","cyan","MaxHeadSize",10);
text(rootX(1) * 2 ,rootX(2) * 2, 0, "Earth Prograde")

earth_R_hat = -earth_Pos / norm(earth_Pos);
rootX = earth_R_hat * 6378;
tipX = earth_R_hat * 4000;
quiver3(rootX(1),rootX(2), 0, tipX(1), tipX(2), 0,"filled","LineWidth",3,"ShowArrowHead","on", "Color","yellow","MaxHeadSize",10);

text(rootX(1) * 2, rootX(2) * 2, 0, "⊙")

earth_Orbit_R = norm(earth_Pos); % this will be accepted as the circular radius of heliocentric Earth orbit.
earth_Orbit_T = 2 * pi * earth_Orbit_R^(1.5) / sqrt(mu_Sun);
earth_Orbit_w = 2 * pi / earth_Orbit_T;

%% Mars Ephemeris
[mars_Pos, mars_V] = planetEphemeris(jdt, "SolarSystem","Mars");
mars_Orbit_R = norm(mars_Pos);  % this will be accepted as the circular radius of heliocentric Mars orbit.
mars_R_hat = mars_Pos / norm(mars_Pos);

mars_Orbit_T = 2 * pi * mars_Orbit_R^(1.5) / sqrt(mu_Sun);
mars_Orbit_w = 2 * pi / mars_Orbit_T;


%% Maneuvering Parameters
planeThreshold = 20;            % desired plane is accepted as reached if SC is closer than the treshold.
planeChangeStatus = 0;          % plane change status
analytical_dV_planeChange = 2 * magV * sind(lat / 2); % analytical value of the plane change delta-v.

departureStatus = 0;

% dp = dot(mars_Pos(1,1:2), earth_Pos(1,1:2));
% angle = acosd(dp / (earth_Orbit_R * mars_Orbit_R));
SC_R = 6378 + 500;

transfer_e = (mars_Orbit_R - earth_Orbit_R) / (mars_Orbit_R + earth_Orbit_R);
perihelion_velocity = sqrt(mu_Sun * (1 + transfer_e) / earth_Orbit_R);
hyperbolic_excess = perihelion_velocity - earth_Orbit_w * earth_Orbit_R;
injection_velocity = sqrt(hyperbolic_excess^2 + 2 * mu / SC_R);
injection_dv = injection_velocity - sqrt(mu / SC_R);
departure_angle = 180 + acosd(1 / (1 + (SC_R * hyperbolic_excess^2 / mu)));

departureThreshold = 500; % km
departure_coordinates = SC_R * [cosd(departure_angle) -sind(departure_angle); sind(departure_angle) cosd(departure_angle)] * earth_V_hat(1, 1:2)';
departure_coordinates = [departure_coordinates' 0];
departure_coordinates = departure_coordinates / (norm(departure_coordinates) / SC_R);

timestepThreshold = 2e4; % km
timestepStatus = 0;
extendedTimeStep = 300; % s

%% Setup Spacecraft Initial Conditions

dt = 20;          % seconds
T = 0:dt:2e5;     % Time matrix
N = length(T);    % Iteration length


X_SC = zeros(N, 3);
V_SC = X_SC;
A_SC = V_SC;


X_SC(1,:) = R_i;
V_SC(1,:) = V_i;

e = zeros(N,1);
u = e;
ke = u;

a = @(X) -mu * X / norm(X)^3;

for i = 1:N
   A_SC(i,:) = a(X_SC(i,:));
   [X_SC, V_SC] = RK4(a, dt, X_SC,V_SC,i);


   % KE, PE and SE calculated as a metric of accuracy.
   u(i) = -mu / norm(X_SC(i,:));
   ke(i) = 0.5 * norm(V_SC(i,:))^2;
   e(i) = u(i) + ke(i);

   % Plot current position.
   figure(1)
   plot3(X_SC(i,1), X_SC(i, 2), X_SC(i, 3),".","Color","#FF3131");
   pause(0.001)
   rotate(earth, [0 0 1], earth_w*dt, c_Rot)

   % this if block checks whether plane change conditions are met.
   if abs(X_SC(i,3)) < planeThreshold
       disp("Plane change time!")
       X_SC(i,3) = 0;
       X_SC(i + 1, :) = X_SC(i, :);
       [ra, dec] = ECI2raDec(X_SC(i,:));
       V_hat = t_xX(dec, ra) * [1; 0; 0];
       
       V_SC(i + 1, :) = (V_hat * magV)';
       dV_planeChange = V_SC(i + 1,:) - V_SC(i,:);
       V_SC(i,:) =  V_SC(i + 1, :);
       planeThreshold = 0;
       planeChangeStatus = 1;
       disp("Plane change successful!")
       disp("Achieved dV = " + string(norm(dV_planeChange)) + " km/s")
       disp("Analytical value dV = " + string(analytical_dV_planeChange) + " km/s")
       disp("Difference = " + string(abs(norm(dV_planeChange) - analytical_dV_planeChange)) + " km/s");
   end

   if planeChangeStatus == 1
       if abs(norm(X_SC(i,:) - departure_coordinates)) < departureThreshold
           disp("Departure time!")
           X_SC(i,:) = departure_coordinates;
           X_SC(i + 1, :) = X_SC(i, :);
           [ra, dec] = ECI2raDec(X_SC(i,:));
           V_hat = t_xX(dec, ra) * [1; 0; 0];
           V_SC(i + 1, :) = (V_hat * injection_velocity)';
           V_SC(i,:) =  V_SC(i + 1, :);
           departureThreshold = 0;
           departureStatus = 1;
           disp("Departure trajectory achieved!")
           disp("dV = " + injection_dv + " km/s")
       end
      if norm(X_SC(i,:)) > timestepThreshold && timestepStatus == 0
           dt = extendedTimeStep;
           disp("Time step has been increased to " + dt + " seconds.")
           ts_left = length(T) - i;
           T(i:end) = T(i):dt:(T(i) + dt*ts_left);
           timestepStatus = 1;
      end
      if norm(X_SC(i,:)) > R_SOI
          disp("SOI radius reached.")
          soi_timestep = i;
          break
      end
   end

   % Real-time energy plots for monitoring.
   figure(2)
   subplot(3,1,1)
   plot(T(i),e(i),".r");
   xlabel("Time (s)")
   ylabel("SE (km^2/s^2)")
  ylim([-30 5])
   grid on
   hold on
   
   subplot(3,1,2)
   plot(T(i),norm(X_SC(i,:)),'.r')
%    ylim([6700, 7000])
   xlabel("Time (s)")
   ylabel("R (km)")
   grid on
   hold on

   subplot(3,1,3)
   plot(T(i),norm(A_SC(i,:)),'.r')
   xlabel("Time (s)")
   ylabel("Acceleration (km/s^2)")
%    ylim([0,12])
   grid on
   hold on
end

T_SOI = T(soi_timestep);
V_SOI = V_SC(soi_timestep);
X_SOI = X_SC(soi_timestep);
disp("Time elapsed since launch is " + string(T_SOI) + " seconds.")

