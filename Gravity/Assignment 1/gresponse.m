%Purpose:
%   The purpose of this function is to plot the gravity response of a point
%   mass below the surface as a function of your distance from the point
%   mass using the point mass' mass, location, and depth.

%Setup:
%   Coordinate System:
%       Cartesian
%           X
%               Positive    ->  Right
%               Negative    ->  Left 
%           Y
%               Positive    ->  Out of Page
%               Negative    ->  Into Page
%           Z
%               Positive    ->  Down
%               Negative    ->  Up

%Input: 
%   mass    ->  Mass at the point
%                   Units: Kilograms
%   depth   ->  Magnitude of the depth of the point below the x-y plane 
%                   Units: Meters
%   x_mass  ->  Displacement of the mass in x direction (Positive Right)
%                   Units: Meters
%   y_mass  ->  Displacement of the mass in y direction (Positive Out)
%                   Units: Meters
%   spread  ->  Distance of axis extension from the origin
%                   Units: Meters
%   res     ->  Distance between each gravity response measurement
%                   Units: Meters

function gresponse(mass, depth, x_mass, y_mass, spread, res)

% Number of observations along a line of measurement
obs_total = 2 * spread / res;

% Setup arrays to store x and y coordinates and the corresponding gravity
% response
G_z = zeros(obs_total, obs_total);
x_obs = zeros(obs_total, obs_total);
y_obs = zeros(obs_total, obs_total);

% Iterate through every observation point on the observation grid and store
% data representing each measurement
for i = 1:obs_total
    for j = 1:obs_total
        % Calculate the x and y coordinates of each observation point
        x_obs(i, j) = -1.0 * spread + i * res;
        y_obs(i, j) = -1.0 * spread + j * res;
        
        % Calculate the gravity response at a given observation location,
        % (x,y)
        g = gresponse_observation(mass, depth, x_mass, y_mass, x_obs(i, j), y_obs(i, j), 0);
        G_z(i,j) = g(3);
    end;
end;

% Graph the data as a surface contour. Each (x, y, gravity response)
surfc(x_obs, y_obs, G_z, 'EdgeColor', 'none');
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Gravity Response (mGal)');
title('Gravity Response of a Point Mass');

%Input: 
%   mass    ->  Mass at the point
%                   Units: Kilograms
%   depth   ->  Magnitude of the depth of the point below the x-y plane 
%                   Units: Meters
%   x_mass  ->  Displacement of the mass in the x direction
%                   Units: Meters
%   y_mass  ->  Displacement of the mass in the y direction
%                   Units: Meters
%   x_obs   ->  Displacement of the observation point in the x direction
%                   Units: Meters
%   y_obs   ->  Displacement of the observation point in the y direction
%                   Units: Meters
%   z_obs   ->  Displacement of the observation point in the z direction
%                   Units: Meters

function g = gresponse_observation(mass, depth, x_mass, y_mass, x_obs, y_obs, z_obs)

%Constants
%   big_g   ->  The Universal Gravitational Constant
%                   Units: (Meters)^(3) * (Kilograms)^(-1) * (Seconds)^(-2)
%   g_conv  ->  Conversion ratio from (Meter) * (Seconds)^(-2) to milliGals
%                   Units: None

big_g   = 6.67384E-11;
g_conv  = 100000;

coefficient = (g_conv * big_g * abs(mass)) / (((x_obs - x_mass)^2 + (y_obs - y_mass)^2 + (z_obs + abs(depth))^2)^(3/2));

g(1) = coefficient * (x_obs - x_mass);
g(2) = coefficient * (y_obs - y_mass);
g(3) = coefficient * (z_obs + abs(depth));
