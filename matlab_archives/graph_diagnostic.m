% graph_diagnostic.m
% graphs relevant info for ceating gaussians in breath cycle to represent
% firing rate
%
% graph the breath_cycle, its peaks, and the gaussians as they fit in them
%
figure
hold on
plot(time_vec, BRTH_data)
plot(p_times, ones(1, length(p_times)), 'rx')
top=max(total_poisson_rate);
for i=1:length(p_times);
    plot([p_times(i) p_times(i)], [0 top],'r')
end

plot(t(1:length(total_poisson_rate)), total_poisson_rate,'color' ,[1 0 1])
plot(store_gauss_centers, ones(1, length(store_gauss_centers)),'o','color',[0 1 1])
ax = axis;
axis([19 25 -1 breathing_peak_rate+on_rate]);
plot(breathmin_times, .9*ones(1, length(breathmin_times)), 'x', 'color',[.5 .5 0])
