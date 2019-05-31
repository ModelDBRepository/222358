% centroid_ave_ctrl.m
% Control's centroid based averaging:
% Computes polar average firing rates of a collection of paired control and evoked
% firing rate polar plots by first aligning the control centroids and then
% calculating the average plus or minus ste for the control and evoked
% data.
% 
% The control polar plots are all aligned with their peak firing rate at
% zero angle.  Their associated evoked polar plots are rotated at their
% paired control polar plots shifted angle. The first control polar plot
% (phi_cell_ctrl{1}, r_cell_ctrl{1}) is assumed
% to be paired with the first evoked (phi_cell{1}, r_cell{1}), the second
% control (phi_cell_ctrl{2}, r_cell_ctrl{2}) with the second evoked
% (phi_cell{2}, r_cell{2}) and so forth.
%
% Usage:
% [phi_ctrl, ave_r_ctrl, ste_r_ctrl, phi, ave_r, ste_r] = ...
%     centroid_ave_ctrl(phi_cell_ctrl, r_cell_ctrl, phi_cell, r_cell)
% where phi_cell_ctrl{i}, r_cell_ctrl{i} are the ith control polar plot
% angles and radii, and phi_cell{i}, r_cell{i} are the ith evoked polar
% plot angles and radii.
%
% The returned angle vector phi_ctrl corresponds to each of the other
% control arguments, the average and standard error radius vectors
% ave_r_ctrl, ste_r_ctrl, and the angle vector phi corresponds to each of
% the evoked average and standard error vectors ave_r, ste_r.
% 

% Here we assume that an additional point at phi_cell{i}(end)+delta_phi{i} would
% be the same as that at phi_cell{i}(0), and that this additional point is
% not supplied in the function call. This additional point is supplied
% (created) internally for interpolation purposes.

function [phi_ctrl, ave_r_ctrl, ste_r_ctrl, phi, ave_r, ste_r] = ...
    centroid_ave_ctrl(phi_cell_ctrl, r_cell_ctrl, phi_cell, r_cell)
% function [phi, ave_r, ste_r] = centroid_ave(phi_cell, r_cell)

% assume that the lengths of all the arguments is the same (the number of
% columns can be different but has to match (seperately) in the first or
% second pairs.

gen_plots=1; % graphing switch: 0 for no graphs

num_of_plots = length(phi_cell_ctrl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% control (ctrl) processing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first find the control centroids shift angle and statistics:
centxs_ctrl=zeros(1, num_of_plots);
centys_ctrl=zeros(1,num_of_plots); % either phi or r has the right length
cent_angles_ctrl=zeros(1,num_of_plots);
delta_phi_ctrl = zeros(1, num_of_plots); % the delta angle or breath 
% increment can be different for each plot.  This is the angle between
% adjacent polar points.

for index=1:num_of_plots
    [centxs_ctrl(index), centys_ctrl(index)] = centroid(phi_cell_ctrl{index}, r_cell_ctrl{index});
    clear i; % reassigns i to square root of minus 1
    Z=centxs_ctrl(index)+centys_ctrl(index)*1i;
    clear angle; % makes the angle function available
    phi=angle(Z); % returns angle of Z in -pi to pi
    if phi<0
        phi = phi+2*pi;
        % returns pi to 0 <= phi < 2*phi
    end
    cent_angles_ctrl(index)=phi;
    
    % for use later to supply an additional point to interpolate to. It is
    % possible that these are different sizes (different numbers of angles
    % per polar plot) in which case would need to figure out whether it
    % would work to choose the smallest or need the maximum.
    many_delta_phi_ctrl = diff(phi_cell_ctrl{index});
    delta_phi_ctrl(index)=many_delta_phi_ctrl(1); % for now just pick the
    % first one until need to figure the more complicated cases of
    % different numbers of angles per polar plot
end

% rotate all the polars so their centroids are aligned at 0

for index=1:num_of_plots
    new_phi_cell_ctrl{index} = mod(phi_cell_ctrl{index} - cent_angles_ctrl(index), 2*pi);
end

% find all used angles and then use those to calculate the mean and ste
% interpolated radii
all_angles_ctrl=[];
for index=1:num_of_plots
    all_angles_ctrl = [all_angles_ctrl; new_phi_cell_ctrl{index}];
end
unique_all_angles_ctrl = unique(all_angles_ctrl, 'sorted');

% interpolate the centroids to the new angles

new_r_ctrl = zeros(num_of_plots, length(unique_all_angles_ctrl));

for cell_index=1:num_of_plots % the cell index here can refer to an electrophysiological cell
    for angle_index = 1:length(unique_all_angles_ctrl)
        % in the below extend the interpolated regions on both sides (angles less than 0 and greater than 2pi) to
        % prevent NaN's from showing up in the linear interpolation:
        new_r_ctrl(cell_index, angle_index)=interp1([phi_cell_ctrl{cell_index}(1)-delta_phi_ctrl(cell_index) phi_cell_ctrl{cell_index} phi_cell_ctrl{cell_index}(end)+delta_phi_ctrl(cell_index)], ...
            [r_cell_ctrl{cell_index}(end) r_cell_ctrl{cell_index} r_cell_ctrl{cell_index}(1)], mod(unique_all_angles_ctrl(angle_index)+cent_angles_ctrl(cell_index),2*pi));
    end
end

phi_ctrl=mod(unique_all_angles_ctrl', 2*pi);
ave_r_ctrl = mean(new_r_ctrl);
ste_r_ctrl = std(new_r_ctrl)/sqrt(num_of_plots);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% evoked processing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first find the evoked centroids shift angle and statistics:
centxs=zeros(1, num_of_plots);
centys=zeros(1,num_of_plots); % either phi or r has the right length
cent_angles=zeros(1,num_of_plots);
delta_phi = zeros(1, num_of_plots); % the delta angle or breath 
% increment can be different for each plot.  This is the angle between
% adjacent polar points.

for index=1:num_of_plots
    [centxs(index), centys(index)] = centroid(phi_cell{index}, r_cell{index});
    clear i; % reassigns i to square root of minus 1
    Z=centxs(index)+centys(index)*1i;
    clear angle; % makes the angle function available
    phi=angle(Z); % returns angle of Z in -pi to pi
    if phi<0
        phi = phi+2*pi;
        % returns pi to 0 <= phi < 2*phi
    end
    cent_angles(index)=phi;
    
    % for use later to supply an additional point to interpolate to. These
    % are done on a polar plot by polar plot basis so it is possible that
    % these delta_phi's are different at different "index" values
    many_delta_phi = diff(phi_cell{index});
    delta_phi(index)=many_delta_phi(1); % use first because presumably uniform delta angles for each polar plot
end

% now rotate the evoked using the same angle as their associated ctrl centroid and find evoked statistics

for index=1:num_of_plots
    new_phi_cell{index} = mod(phi_cell{index} - cent_angles_ctrl(index), 2*pi);
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
            [r_cell{cell_index}(end) r_cell{cell_index} r_cell{cell_index}(1)], mod(unique_all_angles(angle_index)+cent_angles_ctrl(cell_index),2*pi)); % note the
        % shift by cent_angles_ctrl instead of cent_angles so that the
        % evoked is shifted by its associated control centroid.  Also note
        % the shift is in the opposite direction since assigning this way
        % reverse rotates the radii to where they need to be.
    end
end

phi=mod(unique_all_angles', 2*pi);
ave_r = mean(new_r);
ste_r = std(new_r)/sqrt(num_of_plots);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% graphing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% delay graphing until the largest radii point can be found and graphed
% first

% control graphing

% make all controls bluish and all evoked redish
if gen_plots
    figure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Find a scale for the polar plot which is set by the first "polar"
    % executed.  We arbitrarily choose to plot a white (invisible) dot at
    % the maximum in set: {all of the polar plot radii for evoked and control
    % and also ave+ste radius over both the evoked and control}.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    all_radii=[]; % find all big radii so a max can be plotted to scale graph
    for index=1:length(r_cell_ctrl) % r_cell_ctrl and r_cell same length hopefully
        all_radii=[all_radii r_cell_ctrl{index}];
        all_radii=[all_radii r_cell{index}];
    end
    all_radii=[all_radii ave_r_ctrl+ste_r_ctrl];
    all_radii=[all_radii ave_r+ste_r];

    max_radius=max(all_radii);
    polar([0],[max_radius],'w'); % plot a white dot at max_radius to allow room for all
    hold on

    for index=1:num_of_plots
        h=polar([phi_cell_ctrl{index} phi_cell_ctrl{index}(1)],[r_cell_ctrl{index} r_cell_ctrl{index}(1)],'--');
        h.Color=[0 .1 (num_of_plots+1-index)/(2*num_of_plots)];
        h=polar([phi_cell{index} phi_cell{index}(1)],[r_cell{index} r_cell{index}(1)],'--');
%         h.Color=[index/(2*num_of_plots) .75 (num_of_plots+1-index)/(2*num_of_plots)];
        h.Color=[(num_of_plots+1-index)/(2*num_of_plots) .1 0];
    end
    title('(Unrotated) Control and evoked polar plots')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % plot control
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    polar([0],[max_radius],'w'); % plot a white dot at max_radius to allow room for all
    hold on

    for index=1:num_of_plots
        h=polar([new_phi_cell_ctrl{index} new_phi_cell_ctrl{index}(1)],[r_cell_ctrl{index} r_cell_ctrl{index}(1)],'--');
        h.Color=[index/(2*num_of_plots) .5 (num_of_plots+1-index)/(2*num_of_plots)];
    end
    h1=polar([phi_ctrl phi_ctrl(1)], [ave_r_ctrl ave_r_ctrl(1)],'r');
    h1.Color=[.5 0 0];
    h2=polar([phi_ctrl phi_ctrl(1)], [ave_r_ctrl+ste_r_ctrl ave_r_ctrl(1)+ste_r_ctrl(1)],'-.r');
    h2.Color=[.3 0 0];
    h3=polar([phi_ctrl phi_ctrl(1)], [ave_r_ctrl-ste_r_ctrl ave_r_ctrl(1)-ste_r_ctrl(1)],'-.r');
    h3.Color=[.3 0 0];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % plot evoked
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for index=1:num_of_plots
        h=polar([new_phi_cell{index} new_phi_cell{index}(1)],[r_cell{index} r_cell{index}(1)],'--');
%         h.Color=[index/(2*num_of_plots) 0.25 (num_of_plots+1-index)/(2*num_of_plots)];
        h.Color=[0 0.1 (num_of_plots+1-index)/(2*num_of_plots)];
    end
    h1=polar([phi phi(1)], [ave_r ave_r(1)],'r');
    h1.Color=[0 0 .5];
    h2=polar([phi phi(1)], [ave_r+ste_r ave_r(1)+ste_r(1)],'-.r');
    h2.Color=[0 0 .25];
    h3=polar([phi phi(1)], [ave_r-ste_r ave_r(1)-ste_r(1)],'-.r');
    h3.Color=[0 0 .25];
    title('Rotated polar plots and ave plus minus stes')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % just ave's +- ste polar
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    polar([0],[max_radius],'w'); % plot a white dot at max_radius to allow room for all
    hold on
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % plot control
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    h1=polar([phi_ctrl phi_ctrl(1)], [ave_r_ctrl ave_r_ctrl(1)],'r');
    h1.Color=[.5 0 0];
    h2=polar([phi_ctrl phi_ctrl(1)], [ave_r_ctrl+ste_r_ctrl ave_r_ctrl(1)+ste_r_ctrl(1)],'-.r');
    h2.Color=[.25 0 0];
    h3=polar([phi_ctrl phi_ctrl(1)], [ave_r_ctrl-ste_r_ctrl ave_r_ctrl(1)-ste_r_ctrl(1)],'-.r');
    h3.Color=[.25 0 0 ];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % plot evoked
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    h1=polar([phi phi(1)], [ave_r ave_r(1)],'r');
    h1.Color=[0.5 0 0];
    h2=polar([phi phi(1)], [ave_r+ste_r ave_r(1)+ste_r(1)],'-.r');
    h2.Color=[0.25 0 0];
    h3=polar([phi phi(1)], [ave_r-ste_r ave_r(1)-ste_r(1)],'-.r');
    h3.Color=[0.25 0 0];
    title('Ave plus minus ste for control and evoked')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % just all the rotated centroids
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Find a scale for the polar plot which is set by the first "polar"
    % executed.  We arbitrarily choose to plot a white (invisible) dot at
    % the maximum in set: {all of the polar plot radii for evoked and control
    % and also ave+ste radius over both the evoked and control}.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    polar([0],[max_radius],'w'); % plot a white dot at max_radius to allow room for all
    hold on
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % plot control
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for index=1:num_of_plots
        h=polar([new_phi_cell_ctrl{index} new_phi_cell_ctrl{index}(1)],[r_cell_ctrl{index} r_cell_ctrl{index}(1)],'--');
        % h.Color=[index/(2*num_of_plots) .5 (num_of_plots+1-index)/(2*num_of_plots)];
        h.Color=[0 .1 (num_of_plots+1-index)/(2*num_of_plots)];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % plot evoked
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for index=1:num_of_plots
        h=polar([new_phi_cell{index} new_phi_cell{index}(1)],[r_cell{index} r_cell{index}(1)],'--');
        %h.Color=[index/(2*num_of_plots) 0.25 (num_of_plots+1-index)/(2*num_of_plots)];
        h.Color=[index/(2*num_of_plots) 0.1 0];
    end
    title('Rotated polar plots')
    
end
disp('Centroid tables')
disp('Evoked')
fprintf('x\ty\tphi (degrees)\n')
for index=1:length(centxs)
    fprintf([num2str(centxs(index)) '\t' num2str(centys(index)) '\t' num2str(rad2deg(cent_angles(index))) '\n'])
end
disp('Control')
fprintf('x\ty\tphi (degrees)\n')
for index=1:length(centxs)
    fprintf([num2str(centxs_ctrl(index)) '\t' num2str(centys_ctrl(index)) '\t' num2str(rad2deg(cent_angles_ctrl(index))) '\n'])
end

% disp(['max_radius = ' num2str(max_radius)])
