% polar_plot1a.m
% experimental breath sorted instantaneous firing rate
% and control breath sorted instantaneous firing rate polar plots
figure

%title ('Experimental Breath Sorted Instantaneous Firing Rate','fontsize',label_fontsize)
polar([xcenters xcenters(1)], [FR_breath_angles+ste_error_exp FR_breath_angles(1)+ste_error_exp(1)],'y')

hold on
polar([xcenters xcenters(1)], [FR_breath_angles FR_breath_angles(1)],'r')
polar([xcenters xcenters(1)], [FR_breath_angles-ste_error_exp FR_breath_angles(1)-ste_error_exp(1)],'y')
axis equal
% calculate centroid
[centx, centy] = centroid(xcenters, FR_breath_angles);

plot([0 centx],[0 centy],'r','linewidth',4)

% Control

%title('Control Breath Sorted Instantaneous Firing Rate','fontsize',label_fontsize)
polar([xcenters xcenters(1)], [FR_breath_angles_control+ste_error_control FR_breath_angles_control(1)+ste_error_control(1)], 'y')

hold on
% xcenters(1) repeats the first value to close the polar plot
polar([xcenters xcenters(1)],[FR_breath_angles_control FR_breath_angles_control(1)],'k')
polar([xcenters xcenters(1)], [FR_breath_angles_control-ste_error_control FR_breath_angles_control(1)-ste_error_control(1)], 'y')
% ylabel('spikes/(num_of_stims*control_response_period) (Hz)', 'Interpreter', 'none','fontsize',label_fontsize)
% xlabel(['breath angle bins (radians) 2*pi/' num2str(num_of_breath_increments) ' bin widths','fontsize',label_fontsize])
[centx_control, centy_control]=centroid(xcenters, FR_breath_angles_control);
axis equal
plot([0 centx_control],[0 centy_control],'k','linewidth',4)
% centroid(xcenters, FR_breath_angles_control)

clear i; % reassigns i to square root of minus 1
Z=centx+centy*i;
clear angle; % makes the angle function available
phi=angle(Z); % returns angle of Z in -pi to pi
if phi<0
    phi = phi+2*pi;
    % returns pi to 0 <= phi < 2*phi
end
disp(['For response: ********************************'])
disp(['response_centroid_angle = ' num2str(phi)])
disp(['response_centroid_magnitude = ' num2str(abs(Z))])
disp(['response_average_polar_radius = ' num2str(mean(FR_breath_angles))])
disp(['response_ratio_of_centroid_magnitude_to_average_polar_radius = ' num2str(abs(Z)/mean(FR_breath_angles))])
disp(['response_polar_area = ' num2str(polar_area([xcenters xcenters(1)], [FR_breath_angles FR_breath_angles(1)]))])
disp(['response_peak_firing = ' num2str( max(FR_breath_angles)) ])
disp(['response_directivity = ' num2str(directivity(FR_breath_angles))])

% now for control
Z=centx_control+centy_control*i;
phi=angle(Z); % returns angle of Z in -pi to pi
if phi<0
    phi = phi+2*pi;
    % returns pi to 0 <= phi < 2*phi
end
disp(['For control: ********************************'])
disp(['control_centroid_angle = ' num2str(phi)])
disp(['control_centroid_magnitude = ' num2str(abs(Z))])
disp(['control_average_polar_radius = ' num2str(mean(FR_breath_angles_control))])
disp(['control_ratio_of_centroid_magnitude_to_average_polar_radius = ' num2str(abs(Z)/mean(FR_breath_angles_control))])
disp(['control_polar_area = ' num2str(polar_area([xcenters xcenters(1)], [FR_breath_angles_control FR_breath_angles_control(1)]))])
disp(['control_peak_firing = ' num2str( max(FR_breath_angles_control)) ])
disp(['control_directivity = ' num2str(directivity(FR_breath_angles_control)) ])


