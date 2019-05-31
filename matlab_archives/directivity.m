% directivity.m
% Directivity is the max power divided by the average power.
% (power is the square of the signal at each angle).
%
% d = directivity(radii);
% Where the radii could for example be firing_rates of neurons

function d = directivity(rho)

rho_squared = rho.*rho;
d = max(rho_squared)/mean(rho_squared);
