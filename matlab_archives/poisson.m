%poisson.m
% Creates a spike train from a poisson rate.
% Usage:
% spike = poisson(poisson_rate, dt)
% The spike vector elements are one at the times of spikes generated
% from the poisson rate whose vector elements are spaced dt apart.
% The units are seconds.

function [spike]=poisson(poisson_rate, dt)
pspike = dt*poisson_rate;
% create spike_bins = 1 if spike, 0 if no spike in dt bins
spike=pspike > rand(size(poisson_rate)); % derived from web example:
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/85119

end
