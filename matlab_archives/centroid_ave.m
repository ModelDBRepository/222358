% centroid_ave.m
% Computes polar average firing rate of a collection of polar plots by
% first aligning their centroids and then calculating the average plus or
% minus ste.
%
% Usage:
% [phi, ave_r, ste_r] = centroid_based_average(phi_cell, r_cell)
% where phi_cell{i}, r_cell{i} are the ith polar plots angles and radii
% Here we assume that an additional point at phi_cell{i}(end)+delta_phi{i} would
% be the same as that at phi_cell{i}(0), and that this additional point is
% not supplied in the function call. This additional point is supplied
% internally for interpolation purposes.
%

function [phi, ave_r, ste_r] = centroid_ave(phi_cell, r_cell)
gen_plots=1;

num_of_plots = length(phi_cell);
% first find the centroids:
centxs=zeros(1, num_of_plots);
centys=zeros(1,num_of_plots); % either phi or r has the right length
cent_angles=zeros(1,num_of_plots);
delta_phi = zeros(1, num_of_plots); % the delta angle or breath increment can be different for each plot

for index=1:num_of_plots
    [centxs(index), centys(index)] = centroid(phi_cell{index}, r_cell{index});
    clear i; % reassigns i to square root of minus 1
    Z=centxs(index)+centys(index)*i;
    clear angle; % makes the angle function available
    phi=angle(Z); % returns angle of Z in -pi to pi
    if phi<0
        phi = phi+2*pi;
        % returns pi to 0 <= phi < 2*phi
    end
    cent_angles(index)=phi;
    
    % for use later to supply an additional point to interpolate to
    many_delta_phi = diff(phi_cell{index});
    delta_phi(index)=many_delta_phi(1); % they are all the same so just pick first one
end

% rotate all the polars so their centroids are aligned at 0

for index=1:num_of_plots
    new_phi_cell{index} = mod(phi_cell{index} - cent_angles(index), 2*pi);
end
if gen_plots
    figure
    for index=1:num_of_plots
        h=polar([new_phi_cell{index} new_phi_cell{index}(1)],[r_cell{index} r_cell{index}(1)],'--');
        h.Color=[index/(2*num_of_plots) .5 (num_of_plots+1-index)/(2*num_of_plots)];
        hold on
    end
end

% find all used angles and then use those to calculate the mean and ste
% interpolated radii
all_angles=[];
for index=1:num_of_plots
    all_angles = [all_angles; new_phi_cell{index}];
end
unique_all_angles = unique(all_angles, 'sorted');

% interpolate the centroids to the new angles

new_r = zeros(num_of_plots, length(unique_all_angles));

for cell_index=1:num_of_plots % the cell index here can refer to an electrophysiological cell
    for angle_index = 1:length(unique_all_angles)
        % in the below extend the interpolated regions on both sides (angles less than 0 and greater than 2pi) to
        % prevent NaN's from showing up in the linear interpolation:
        new_r(cell_index, angle_index)=interp1([phi_cell{cell_index}(1)-delta_phi(cell_index) phi_cell{cell_index} phi_cell{cell_index}(end)+delta_phi(cell_index)], ...
            [r_cell{cell_index}(end) r_cell{cell_index} r_cell{cell_index}(1)], mod(unique_all_angles(angle_index)+cent_angles(cell_index),2*pi));
    end
end

phi=mod(unique_all_angles', 2*pi);
ave_r = mean(new_r);
ste_r = std(new_r)/sqrt(num_of_plots);

if gen_plots
    h1=polar([phi phi(1)], [ave_r ave_r(1)],'r');
    h1.Color=[.5 0 0];
    h2=polar([phi phi(1)], [ave_r+ste_r ave_r(1)+ste_r(1)],'-.r');
    h2.Color=[.5 .5 0];
    h3=polar([phi phi(1)], [ave_r-ste_r ave_r(1)-ste_r(1)],'-.r');
    h3.Color=[0 .5 .5];
end

    