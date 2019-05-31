% centroid_ave.m
% Computes polar average firing rate of a collection of polar plots by
% first aligning their centroids and then calculating the average plus or
% minus ste.
%
% Usage:
% [phi, ave_r, ste_r] = centroid_based_average(phi_cell, r_cell)
% where phi_cell{i}, r_cell{i} are the ith polar plots angles and radii

function [phi, ave_r, ste_r] = centroid_ave(phi_cell, r_cell)

% first find the centroids:
centxs=zeros(1, length(phi_cell));
centys=zeros(1,length(phi_cell)); % either phi or r has the right length
cent_angles=zeros(1,length(phi_cell));

for index=1:length(centxs)
    [centxs(index), centys(index)] = centroid(phi_cell{index}, r_cell{index});
    clear i; % reassigns i to square root of minus 1
    Z=centx+centy*i;
    clear angle; % makes the angle function available
    phi=angle(Z); % returns angle of Z in -pi to pi
    if phi<0
        phi = phi+2*pi;
        % returns pi to 0 <= phi < 2*phi
    end
    cent_angles(index)=phi;
end

% rotate all the polars so their centroids are aligned at 0

for index=1:length(centxs)
    phi_cell{index} = phi_cell{index} - cent_angles(index);
end

% find all used angles and then use those to calculate the mean and ste
% interpolated radii
all_angles=[];
for index=1:length(centxs)
    all_angles = [all_angles; phi_cell{index}];
end
unique_all_angles = unique(all_angles, 'sorted');

% interpolate the centroids to the new angles

new_r = zeros(length(r_cell), length(unique_all_angles));

for cell_index=1:length(r_cell)
    for angle_index = 1:length(unique_all_angles)
        new_r(cell_index, angle_index)=interpolate(unique_all_angles(angle_index), r_cell{cell_index}, phi_cell{cell_index});
    end
end

phi=unique_all_angles;
ave_r = mean(new_r);
ste_r = std(new_r)/sqrt(length(r_cell);
