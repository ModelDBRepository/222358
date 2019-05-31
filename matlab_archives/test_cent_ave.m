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

% first test batch
% evoked_cell{1}=[1 1 1 10];
% evoked_cell{2}=[3 3 4 3];
% 
% control_cell{1} = [2 2 10 2];
% control_cell{2} = [1 2 1 5];

% second test batch
control_cell{1} = [2   2   10   2];
evoked_cell{1}=   [2.5 2.5 4   11];

control_cell{2} = [1 2 1 5];
evoked_cell{2}=   [5 3 1 1];

% third test batch
% control_cell{1} = [2   2   10   2];
% evoked_cell{1}=   [2.5 2.5 4   11];
% 
% for index=2:9
%     control_cell{index} = control_cell{1} + rand(1,4);
%     evoked_cell{index}  = evoked_cell{1}  + rand(1,4);
%     phi_cell{index} = phi_cell{1};
% end

% fourth test batch
control_cell{1} = [2   2   10   2];
evoked_cell{1}=   [2.5 2.5 4   11];

for index=2:9
    control_cell{index} = control_cell{1} + [0 0 3*rand(1,1) 0];
    evoked_cell{index}  = evoked_cell{1}  + [0 0 0 4*rand(1,1)];
    phi_cell{index} = phi_cell{1};
end

% 
% control_cell{2} = [1 2 1 5];
% evoked_cell{2}=   [5 3 1 1];
% 
% r_cell=evoked_cell;

%centroid_ave
% [phi, ave_r, ste_r] = centroid_ave(phi_cell, evoked_cell);
% [phi_c, ave_r_c, ste_r_c] = centroid_ave(phi_cell, control_cell);

% centroid_ave_ctrl
% [phi_ctrl, ave_r_ctrl, ste_r_ctrl, phi, ave_r, ste_r] = ...
%     centroid_ave_ctrl(phi_cell_ctrl, r_cell_ctrl, phi_cell, r_cell);
[phi_ctrl, ave_r_ctrl, ste_r_ctrl, phi, ave_r, ste_r] = ...
    centroid_ave_ctrl(phi_cell, control_cell, phi_cell, evoked_cell);
