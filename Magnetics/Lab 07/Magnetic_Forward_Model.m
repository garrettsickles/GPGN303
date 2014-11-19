%For the problem in lab 07, run the function:
%   Magnetic_Forward_Model(50000, 45, 25, 0, 0, 250, 200, 0.05, -775, 800, -775, 800, 25)

%Setup:
%   Coordinate System:
%       Cartesian
%           Easting
%               Positive    ->  Right
%               Negative    ->  Left 
%           Northing
%               Positive    ->  Out of Page
%               Negative    ->  Into Page
%           Z
%               Positive    ->  Down
%               Negative    ->  Up

%Input: 
%   strength->  Strength of the inducing field
%                   Units: nanoTeslas
%   inc     ->  Inclination of the inducing field 
%                   Units: Degrees
%   dec     ->  Declination of the inducing field 
%                   Units: Degrees
%   N_body  ->  Displacement of the body in northing direction
%                   Units: Meters
%   E_body  ->  Displacement of the body in easting direction
%                   Units: Meters
%   depth   ->  Distance bleow the surface of the center of the sphere
%                   Units: Meters
%   radius  ->  Radius of the subsurficial sphere
%                   Units: Meters
%   kappa   ->  Magnetic susceptibility of the subsurficial sphere
%                   Units: Meters
%   N_min   ->  Minimum northing value of the observation grid
%                   Units: Meters
%   N_max   ->  Maximum northing value of the observation grid
%                   Units: Meters
%   E_min   ->  Minimum easting value of the observation grid
%                   Units: Meters
%   E_max   ->  Maximum easting value of the observation grid
%                   Units: Meters
%   res     ->  Resolution of measurement point on the observation grid
%                   Units: Meters

function Magnetic_Forward_Model(strength, inc, dec, N_body, E_body, depth, radius, kappa, N_min, N_max, E_min, E_max, res)

B_naught(1, 1) = strength * cos(inc * pi / 180) * cos(dec * pi / 180);
B_naught(2, 1) = strength * cos(inc * pi / 180) * sin(dec * pi / 180);
B_naught(3, 1) = strength * sin(inc * pi / 180);
B_naught_hat = transpose(B_naught) ./ strength;

%Establish Northing and Easting ranges that determine the gridded area 
%   over which calculations will be taken.

N = N_min:res:N_max; %m
E = E_min:res:E_max; %m

%Establish the total field anomaly
Total = zeros(length(N),length(E));

%Iterate over all of the observation points and calculate the total-field
%   anomaly at each point
for i=1:length(N)
    for j=1:length(E)
        B_a = (kappa * radius^3 / 3) .* (Magnetic_Tensor_Matrix(N_body, E_body, depth, N(i), E(j), 0) * B_naught);
        Total(i, j) = dot(transpose(B_a), B_naught_hat);
    end
end

%Plot the total-field anomaly as a interpolated 2D plot
surf(E, N, Total,'EdgeColor', 'None', 'facecolor', 'interp');
view(2);
title('Total Field Anomaly');
cb = colorbar('location','eastoutside');
xlabel(cb, 'nT')



%Input:
%   N_body  ->  Displacement of the body in northing direction
%                   Units: Meters
%   E_body  ->  Displacement of the body in easting direction
%                   Units: Meters
%   depth   ->  Distance bleow the surface of the center of the sphere
%                   Units: Meters
%   N_obs   ->  Displacement of the observation point in northing direction
%                   Units: Meters
%   E_obs   ->  Displacement of the observation point in easting direction
%                   Units: Meters
%   z_obs   ->  Distance above the surface of the observation point
%                   Units: Meters
function T = Magnetic_Tensor_Matrix(N_body, E_body, depth, N_obs, E_obs, z_obs)

%Constants
%   r_sqrd  ->  Radius from the observation point squared
%                   Units: (Meters)^2

r_sqrd = ((N_obs - N_body)^2 + (E_obs - E_body)^2 + (z_obs + abs(depth))^2);

T = zeros([3,3]);

%Define the Magnetic Tensor Matrix
T(1, 1) = 3 * (N_obs - N_body)^(2) * r_sqrd^(-5/2) - r_sqrd^(-3/2);
T(1, 2) = 3 * (E_obs - E_body) * (N_obs - N_body) * r_sqrd^(-5/2);
T(1, 3) = -3 * (z_obs + abs(depth)) * (N_obs - N_body) * r_sqrd^(-5/2);
T(2, 1) = T(1, 2);
T(2, 2) = 3 * (E_obs - E_body)^(2) * r_sqrd^(-5/2) - r_sqrd^(-3/2);
T(2, 3) = -3 * (z_obs + abs(depth)) * (E_obs - E_body) * r_sqrd^(-5/2);
T(3, 1) = T(1, 3);
T(3, 2) = T(2, 3);
T(3, 3) = 3 * (z_obs + abs(depth))^(2) * r_sqrd^(-5/2) - r_sqrd^(-3/2);