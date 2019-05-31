% polar_area.m
%
% Usage:
% area_in_curve = polar_area([phi(1) phi(2) ... phi(n) phi(1)],[r(1) r(2) ... r(n) r(1)]);
% Computes the area under the polar coordinates (phi, r) where phi and r
% are the angles in radians and and corresponding radii (magnitudes) of the points in the plot.
% As you can see in the usage example this function assumes that the
% first and last point in each vector are the same value (closed curve).

function area = polar_area(phi, r)

% see http://www.mathopenref.com/triangleareasas.html for formula derivation

C=diff(phi); % find the angles in all the "triangles"
area = 0; %
for index=1:length(C)
    area = area + (0.5) * r(index) * r(index+1) * sin(C(index));
end
