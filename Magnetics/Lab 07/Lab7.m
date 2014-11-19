%This code calculates the total field magnetic response of a spherical target 
%beneath the ground surface and plots a 3D graph of the anomaly.

%Establish northing and easting observation locations on surface.
Nobs = 0;
Eobs = 0;

%Establish depth, radius and susceptibility of buried body.
depth = 250; %m
radius = 200; %m
sus = .05; %unitless

%Establish the strentgh, inclination, and declination of the Earth's
%inducing field
strength = 50000; %nT
inc = 45; %degrees
dec = 25; %degrees
Bx = (strength .* cos(inc .* pi ./ 180) .* cos(dec .* pi ./ 180));
By = (strength .* cos(inc .* pi ./ 180) .* sin(dec .* pi ./ 180));
Bz = (strength .* sin(inc .* pi ./ 180));

%Establish Northing and Easting ranges that determine the gridded area 
%over which calculations will be taken.
N = -775:25:800; %m
E = -775:25:800; %m

%Establish the total field anomaly
Total = zeros(length(N),length(E));

%Establish initial z value at the ground surface and an
%array of 'zeros' for the magnetic response, T.
z = 0;
Bax = zeros(length(N),length(E));
Bay = zeros(length(N),length(E));
Baz = zeros(length(N),length(E));

%Establish each tensor as an array of 'zeros'
Tnn = zeros(length(N),length(E));
Tee = zeros(length(N),length(E));
Tne = zeros(length(N),length(E));
Ten = zeros(length(N),length(E));
Tez = zeros(length(N),length(E));
Tze = zeros(length(N),length(E));
Tnz = zeros(length(N),length(E));
Tzn = zeros(length(N),length(E));
Tzz = zeros(length(N),length(E));

%Establish a vector of zeros to represent the magnetic dipole moment.
m = zeros(3);

%Permeability of Free Space constant
u = 1.2566E-6; %Wb/(Am)

%The following for loop navigates through each N and E point on the grid,
%calculates the components of the tensors, anomalous field and the total
%field anomaly

for i=1:length(N)
    for j=1:length(E)
        
%Equations for the tensors.
Tee(i,j) = ((2 .* (( N(i) - Nobs )^2))-((z - depth)^2)-((E(j) - Eobs)^2))./ ((( (N(i)) - Nobs ) .^2 + ( (E(j)) - Eobs ) .^2 + ( z - depth ) .^2 ) .^ (5/2));
Tnn(i,j) = ((2 .* (( E(j) - Eobs )^2))-((z - depth)^2)-((N(i) - Nobs)^2))./ ((( (N(i)) - Nobs ) .^2 + ( (E(j)) - Eobs ) .^2 + ( z - depth ) .^2 ) .^ (5/2));
Tne(i,j) = ( 3 .* ( (N(i)) - Nobs ) .* ( (E(j)) - Eobs ) ./ ((( (N(i)) - Nobs ) .^2 + ( (E(j)) - Eobs ) .^2 + ( z - depth ) .^2 ) .^ (5/2)));
Ten(i,j) = Tne(i,j);
Tez(i,j) = (3 .* ( (N(i)) - Nobs ) .* ( z - depth ) ./ ((( (N(i)) - Nobs ) .^2 + ( (E(j)) - Eobs ) .^2 + ( z - depth ) .^2 ) .^ (5/2)));
Tze(i,j) = Tez(i,j);
Tnz(i,j) = (3 .* ( (E(j)) - Eobs ) .* ( z - depth ) ./ ((( (N(i)) - Nobs ) .^2 + ( (E(j)) - Eobs ) .^2 + ( z - depth ) .^2 ) .^ (5/2)));
Tzn(i,j) = Tnz(i,j);
Tzz(i,j) = ((2 .* (( z - depth )^2))-((N(i) - Nobs)^2)-((E(j) - Eobs)^2))./ ((( (N(i)) - Nobs ) .^2 + ( (E(j)) - Eobs ) .^2 + ( z - depth ) .^2 ) .^ (5/2));

%Equation for the components of the anomalous field.
Bax(i,j) = (sus ./ 3) .* (radius^3) .* (Tee(i,j) .* Bx + Ten(i,j) .* By + Tez(i,j) .* Bz);
Bay(i,j) = (sus ./ 3) .* (radius^3) .* (Tne(i,j) .* Bx + Tnn(i,j) .* By + Tnz(i,j) .* Bz);
Baz(i,j) = (sus ./ 3) .* (radius^3) .* (Tze(i,j) .* Bx + Tzn(i,j) .* By + Tzz(i,j) .* Bz);

%Equation for the amplitude of the total field anomaly.
Total(i,j) = Bax(i,j) .* (cos(inc .* pi ./ 180) .* cos(dec .* pi ./ 180)) + Bay(i,j) .* cos(inc .* pi ./ 180) .* sin(dec .* pi ./ 180) + Baz(i,j) .* sin(inc .* pi ./ 180);
    end
end

%A figure of the total field anomaly
figure;
surf(E, N, Total,'EdgeColor', 'None', 'facecolor', 'interp');
view(2);
title('Total Field Anomaly');
cb = colorbar('location','eastoutside');
xlabel(cb, 'nT')

