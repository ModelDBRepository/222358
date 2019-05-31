% diagnostic_script_5.m
% exercises the diagnostic functions to create an artificial tank
% This fifth version uses the breath dependent control and a breath
% dependent stimulus from version 3 (_3) and adds a variable duration
% (see intended durs) for both the light on and off periods.
% The stimuli on and off durations are first chosen from the 8 values in in
% in intended_durs and then smeared out with gaussian noise
% The light on period has a different firing rate gaussian than the light
% off period


intended_durs = [16    33    66   134   267   534  1067  2134]./1000;

disp(['creating stimuli'])
% time is in seconds

num_of_stimuli = 3560; % 3560; % 1760; % 560; % estimated time of experiment = start_time + number of stimuli * inter_stimulus_interval + extra_time
start_time = 20; % time before stimulations
extra_time = 2; % time after stimulations till end of simulated experiment recordings
stimulus_width = 0.2;
inter_stimulus_interval =  .400; % must be greater than stimulus_width or 1
dt = 0.001;

% breathing variables:
gauss_placement_factor = .5; % .2; % a linear factor.  This factor
on_gauss_placement_factor = .25; % .2; % a linear factor.  This factor
% determines the placement of the center of the gaussians.  .5 has the
% center in the middle of the peaks.  Less than .5 moves the center to the
% previous peak and greater than .5, towards the next peak.
width_percentage =  20;  % if 20 then the width at half the peak height is
on_width_percentage = 20 % equal to 20 percent of the distance between the peaks

% The Gaussian's rise up from these background levels:
% The following two firing rate numbers are the elevation of the plateaus that the two
% gaussian mountains, the breath response, and the (light-caused) phase shifted breath response, rise up from
breathing_off_rate = 0.0; % This represents the amount of firing that is spontaneous and breath independent when dark
light_on_off_rate = 0.0; % This represents the amount of firing that is breath independent when light shines
% (The _off_rate part of the variable names refers to being off the
% gaussian.  Since the gaussian goes on forever it really means off the
% part of the gaussian that is significantly large.)

breathing_peak_rate = 40; % Hz for the peak of the gaussian
light_on_peak_rate = 40; % Hz peak for the light on firing rate gaussian 
%[ S_ON, SOFF ] = stimulus_maker( num_of_stimuli, start_time, stimulus_width, inter_stimulus_interval );
[ S_ON, SOFF ] = stimulus_maker( num_of_stimuli, start_time, intended_durs );

total_time = SOFF(end) + extra_time; % end two seconds after end of SOFF
disp(['This data set will include ' num2str(total_time) ' seconds of simulated data'])
on_rate = 20;
off_rate = 2;

%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %
%
%   make a breathing trace
%
%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %

fs=1000; % simplier than 1024 although 1024 should work; % frequency of sampling rate - real data had 24,100 Hz I believe
time_vec = 0:1/fs:total_time;
% breath_period = .6; %
% choose a random breath period in the range .3 to .9

small_breath_period=.3;
large_breath_period=.9;
breath_period  = small_breath_period+(large_breath_period-small_breath_period)*rand(1,1);

theta = 2*pi*time_vec/breath_period; % for use with a breath dependent stimulus on rate
BRTH_data = sin(theta); % start simple with a uniform breath

stimids=[1:num_of_stimuli];
% let glomon be the stimids and glomoff be -1
glomon = stimids';
save('data/glomon.dat','glomon','-ascii')
glomoff = [-1];
save('data/glomoff.dat','glomoff','-ascii')

%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %
%
%   Do analyze_breathing now so can create a breathing_poisson_rate that
%   will get added to the poisson_rate if a stimulus is turned on at the
%   same time.
%
%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %
clear tdt2mat_data

tdt2mat_data.streams.BRTH.data = BRTH_data;
tdt2mat_data.streams.BRTH.fs = fs; % Hz, make the sampling rate handy
addpath('/Users/morse/Documents/projects/VerhagenLab/20140814dev/ORNstim')
%cd ../ORNstim
tdt2mat_data_index = 1;
analyze_breathing
% cd ../diagnostic_simulator

% create the breathing_poisson_rate with a model where the center of the
% gaussian is at some 0 < gauss_placement_factor < 1 between the peaks
% where near 0 is the center of the gaussian is close to the last breath
% peak and near one makes the gaussian close to the next gauss peak
% and the
% gaussian width is width_percentage of the distance between the breath
% peaks.
% these variables are assigned at the top of the script for clarity.

% initialize then set in while loop before:

breathing_poisson_rate=breathing_off_rate*ones(1,floor(total_time/dt)+1);

% loop over the peaks and create the breathing_poisson_rate
store_gauss_centers = zeros(1, length(p_times)-1);
p_index = 1;
while p_index < length(p_times)
    start_p_time = p_times(p_index);
    start_p_index = p_indicies(p_index);
    end_p_time = p_times(p_index+1);
    end_p_index = p_indicies(p_index+1);
    % first get a handle on the scale of gaussian relative to these
    gauss_width  = (width_percentage/100)*(end_p_time - start_p_time);
    gauss_sigma = gauss_width/2; % sigma is the half width
    
    % take care of the starting index for the gaussian x values first
    % start_gauss and end_gauss will hold the indicies in
    % breathing_poisson_rate that need to be assigned to the gaussians for
    % each breath
    
    if p_index == 1
        start_gauss_index = 1; % start the gaussian at the begining for the first peak
    else
        start_gauss_index = p_indicies(p_index-1); % include the gaussian
        % contribution back through the previous peak
    end
    
    end_gauss_index = p_indicies(p_index+1); % this point is always available because
    % we end the while loop on the peak before the last peak
    
    gauss_center_time = gauss_placement_factor*(end_p_time - start_p_time)+start_p_time;
    store_gauss_centers(p_index) = gauss_center_time;
    t_vec = ([start_gauss_index:end_gauss_index]-1)/fs;
    gauss_values = breathing_peak_rate .* gaussmf(t_vec, [gauss_sigma gauss_center_time]);
    %     if p_index< 5
    %         disp(['p_times(p_index) = ' num2str(p_times(p_index)) ', p_times(p_index+1) = ' num2str(p_times(p_index+1))])
    %         disp(['p_indicies(p_index) = ' num2str(p_indicies(p_index)) ', p_indicies(p_index+1) = ' num2str(p_indicies(p_index+1))])
    %         disp([num2str(p_index) ': start_gauss_index = ' num2str(start_gauss_index) ', end_gauss_index = ' num2str(end_gauss_index)])
    %     end
    for index=start_gauss_index:end_gauss_index
        if index<=length(breathing_poisson_rate)
            breathing_poisson_rate(index)=breathing_poisson_rate(index)+gauss_values(index-start_gauss_index+1);
        else
            breathing_poisson_rate(index)=gauss_values(index-start_gauss_index+1);
        end
    end
    p_index = p_index + 1;
end

%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %
%
%   create gaussian rates that correspond to the light on
%
%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %
%
% create the light_on_poisson_rate with a model where the center of the
% gaussian is at some 0 < on_gauss_placement_factor < 1 between the peaks
% where near 0 is the center of the gaussian is close to the last breath
% peak and near one makes the gaussian close to the next gauss peak
% and the
% gaussian width is width_percentage of the distance between the breath
% peaks.
% these variables are assigned at the top of the script for clarity.

% initialize then set in while loop before:

light_on_poisson_rate=light_on_off_rate*ones(1,floor(total_time/dt)+1);

% loop over the peaks and create the light_on_poisson_rate
store_gauss_centers = zeros(1, length(p_times)-1);
p_index = 1;
while p_index < length(p_times)
    start_p_time = p_times(p_index);
    start_p_index = p_indicies(p_index);
    end_p_time = p_times(p_index+1);
    end_p_index = p_indicies(p_index+1);
    % first get a handle on the scale of gaussian relative to these
    gauss_width  = (width_percentage/100)*(end_p_time - start_p_time);
    gauss_sigma = gauss_width/2; % sigma is the half width
    
    % take care of the starting index for the gaussian x values first
    % start_gauss and end_gauss will hold the indicies in
    % light_on_poisson_rate that need to be assigned to the gaussians for
    % each breath
    
    if p_index == 1
        start_gauss_index = 1; % start the gaussian at the begining for the first peak
    else
        start_gauss_index = p_indicies(p_index-1); % include the gaussian
        % contribution back through the previous peak
    end
    
    end_gauss_index = p_indicies(p_index+1); % this point is always available because
    % we end the while loop on the peak before the last peak
    
    gauss_center_time = on_gauss_placement_factor*(end_p_time - start_p_time)+start_p_time;
    store_gauss_centers(p_index) = gauss_center_time;
    t_vec = ([start_gauss_index:end_gauss_index]-1)/fs;
    gauss_values = light_on_peak_rate .* gaussmf(t_vec, [gauss_sigma gauss_center_time]);
    %     if p_index< 5
    %         disp(['p_times(p_index) = ' num2str(p_times(p_index)) ', p_times(p_index+1) = ' num2str(p_times(p_index+1))])
    %         disp(['p_indicies(p_index) = ' num2str(p_indicies(p_index)) ', p_indicies(p_index+1) = ' num2str(p_indicies(p_index+1))])
    %         disp([num2str(p_index) ': start_gauss_index = ' num2str(start_gauss_index) ', end_gauss_index = ' num2str(end_gauss_index)])
    %     end
    for index=start_gauss_index:end_gauss_index
        if index<=length(light_on_poisson_rate)
            light_on_poisson_rate(index)=light_on_poisson_rate(index)+gauss_values(index-start_gauss_index+1);
        else
            light_on_poisson_rate(index)=gauss_values(index-start_gauss_index+1);
        end
    end
    p_index = p_index + 1;
end


%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %
%
%   combine the gaussian rates by selecting from the light on distribution
%   when the light is on (light_on_poisson_rate) and the light off one
%   (breathing_poisson_rate) when the light is off to make phase_shifted_poisson_rate.
%
%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %

phase_shifted_poisson_rate = breathing_poisson_rate; % start with breathing
% poisson rate and then substitute in the phase shifted poisson rate when the light is on

for stim_index=1:length(S_ON)
    [start_index,cv] = searchclosest_fixed(t, S_ON(stim_index));
    [stop_index,cv] = searchclosest_fixed(t, SOFF(stim_index));
    phase_shifted_poisson_rate(start_index:stop_index) = light_on_poisson_rate(start_index:stop_index); % replace rate with light_on_poisson_rate where appropriate
end
total_poisson_rate = phase_shifted_poisson_rate;
    %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %
%
%   finally create the simulated spike output
%
%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %
% soon the spikes will also depend on the breath cycle but for now it
% doesn't:

% for a breathing dependent on_rate it must be assigned values for each
% breath angle
if 0
    breath_cycle_theta_time = 0:1/fs:(breath_period+1/fs); % adding 1/fs insures
    % that the interpolated values include the end point of the breath cycle at 2pi
    
    breath_cycle_theta = (2*pi/breath_period) * breath_cycle_theta_time;
    
    % on_rate as a function of theta
    on_rate = breath_cycle_theta * 60/ (2*pi); % arbitrary function approaches 60 hz at end of breath cycle
    %
    % pi_over_four_shift = floor(length(on_rate)/4);
    %
    % on_rate = circshift(on_rate, [0 pi_over_four_shift]);
    [ poisson_rate ] = monotonic_breath_dependent_on_rate_poisson_rate(total_time, dt, on_rate, breath_cycle_theta, off_rate, S_ON, SOFF, theta ); %
end % stuff specific to monotonic_breath_dependent_rate

% no longer use the simple on_rate when there is light and off rate when there is not
%[ poisson_rate ] = stimulus_implied_poisson_rate(total_time, dt, on_rate, off_rate, S_ON, SOFF );
%max_index=min([length(breathing_poisson_rate) length(poisson_rate)]);

%total_poisson_rate = breathing_poisson_rate(1:max_index) + poisson_rate(1:max_index);
%total_poisson_rate = breathing_poisson_rate(1:max_index)

spikes = poisson(total_poisson_rate, dt); % the spikes output has a value (1=spike, 0 otherwise) for each time step dt

spike_indicies = find(spikes==1); % to convert to a list of spike_times find the indicies of spikes and since
spike_times = spike_indicies*dt; % each index is at a time step of length dt, that number of indices is proportional to time

%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %
%
% create tank tdt2mat
%
% to feed evoke_half_act_... which has commands like the following
%
% ts =     tdt2mat_data.snips.eNeu.ts;
% sortcode=tdt2mat_data.snips.eNeu.sortcode;
% SOFF  =  tdt2mat_data.epocs.SOFF.onset;
% S_ON  =  tdt2mat_data.epocs.S_ON.onset(1:2:end);
% stimids=tdt2mat_data.epocs.S_ON.data(1:2:end);
%
% and analyze_breathing.m which has commands like
%
% a=tdt2mat_data(tdt2mat_data_index).streams.BRTH.data;
% fs=tdt2mat_data(tdt2mat_data_index).streams.BRTH.fs; % Hz, make the sampling rate handy
%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %

tdt2mat_data.snips.eNeu.ts=spike_times;
tdt2mat_data.snips.eNeu.sortcode=ones(1, length(spike_times)); % arbitrarily assign a sortcode of 1 for each spike time
tdt2mat_data.epocs.SOFF.onset = SOFF;
duplicates_S_ON = kron(S_ON, [1 1]);
tdt2mat_data.epocs.S_ON.onset=duplicates_S_ON;
stimids=1:length(S_ON);  % arbitrarily assign whole numbers as stimids
duplicates_stimids = kron(stimids, [1 1]);
tdt2mat_data.epocs.S_ON.data = duplicates_stimids;

% in this diagnostic simulator the below commands were already excuted
% above:
%tdt2mat_data.streams.BRTH.data = BRTH_data;
%tdt2mat_data.streams.BRTH.fs = fs; % Hz, make the sampling rate handy

dtimestamp = timestamp;
filename=['data/tdt2mat_data_' dtimestamp(2:end) '.mat'];
save(filename, 'tdt2mat_data');
disp(['Wrote tdt2mat_data to ' filename]);
cmd = ['cp ' filename ' ../Static20140616/data'];
system(cmd);
disp(['executed: ' cmd]);
fid=fopen('../Static20140616/data/shared_file_pointers.dat','w'); % the first row is the tank the second the breath_info filenames with paths
fprintf(fid, [filename '\n']);
fprintf(fid, ['data/sort_breath/' filename(6:end) '_stimid_ev_angles.dat\n']);
fclose(fid);

