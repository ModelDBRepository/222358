function [ S_ON, SOFF ] = stimulus_maker( num_of_stimuli, start_time, varargin )
%stimulus_maker returns S_ON and SOFF vectors as appear in tdt tank
% Usage:
% [ S_ON, SOFF ] = stimulus_maker( num_of_stimuli, start_time, stimulus_width, inter_stimulus_interval )
%   At a start_time the first stimulus appears.  Subsequent stimuli of
%   length stimulus_width appear at the inter_stimulus_interval to provide
%   num_of_stimuli stimulations
% or
% % [ S_ON, SOFF ] = stimulus_maker( num_of_stimuli, start_time, intended_durations )
% where intended_durations is something like
% intended_durs = [16    33    66   134   267   534  1067  2134]./1000;

if length(varargin)>1
    stimulus_width=varargin{1};
    inter_stimulus_interval=varargin{2};
else
    intended_durations = varargin{1};
end
S_ON = zeros(1, length(num_of_stimuli));
SOFF = zeros(1, length(num_of_stimuli));
if length(varargin)==2
    if length(inter_stimulus_interval)==1
        for index=1:num_of_stimuli
            S_ON(index) = start_time + (index-1) * inter_stimulus_interval;
            SOFF(index) = start_time + (index-1) * inter_stimulus_interval + stimulus_width;
        end
    else
        if length(inter_stimulus_interval)>1
            index=1;
            S_ON(index) = start_time;
            SOFF(index) = start_time  + stimulus_width;
            while index < num_of_stimuli
                index = index + 1;
                % find next S_ON time with a selection from the variable
                % durations in inter_stimulus_interval
                random_index=floor(rand(1,1)*length(inter_stimulus_interval))+1;
                S_ON(index) = SOFF(index-1) + inter_stimulus_interval(random_index);
                SOFF(index) = S_ON(index) + stimulus_width;
            end
        else
            error(['Error in stimulus_make: length(inter_stimulus_interval) = ' num2str(length(inter_stimulus_interval)) ])
        end
    end
else
    offset = randn(2,num_of_stimuli)./1000; % makes intervals slightly off intended values
    cutoff =0.005;
    [row, col, val]=find(abs(offset)>cutoff); % shift offset's greater than 5 ms back to less than 5 ms
    for ind = 1:length(row)
        offset(row(ind),col(ind)) = rem(offset(row(ind),col(ind)),cutoff);
    end
    intended_index = floor(length(intended_durations)*rand(2,num_of_stimuli))+1;
    current_time = start_time;
    for index=1:num_of_stimuli
        % choose a random dark time before S_ON
        current_time = current_time+intended_durations(intended_index(1,index))+offset(1,index);
        S_ON(index)=current_time;
        current_time = current_time+intended_durations(intended_index(2,index))+offset(2,index);
        SOFF(index) = current_time;
    end
end
