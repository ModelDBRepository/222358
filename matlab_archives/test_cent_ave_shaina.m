% test_cent_ave.m
%
% Test 1
% Evoked (1, 1, 1, 10)
% Control (2, 2, 10, 2)
% 
% Test 2
% Evoked (3, 3, 4, 3)
% Control (1, 2, 1, 5)

% [phi, ave_r, ste_r] = centroid_based_average(phi_cell, r_cell)

phi_cell{1}=[0 pi/2 pi 3*pi/2];
phi_cell{2}=phi_cell{1};
phi_cell{3}=phi_cell{2};

evoked_cell{1}=[1 1 1 10];
evoked_cell{2}=[1 1 1 5];
evoked_cell{3}=[1 1 1 15];

control_cell{1} = [2 2 2 5];
control_cell{2} = [2 2 2 5];
control_cell{3} = [2 2 2 5];
r_cell=evoked_cell;

%centroid_ave
[phi, ave_r, ste_r] = centroid_ave(phi_cell, evoked_cell);
[phi_c, ave_r_c, ste_r_c] = centroid_ave(phi_cell, control_cell);

% upgrades:
% make the control and evoked so that the partner evoke's are shifted with the
% same angle as it's associated control.
