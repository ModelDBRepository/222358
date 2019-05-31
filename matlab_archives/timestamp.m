% timestamp.m
% function returns date and time without space

function t=timestamp
c=clock;
t=datestr(c);
p=strfind(t,' ');
t=[t(1:p-1) '_' t(p+1:end)];
p=strfind(t,':');
for i=1:length(p)
    t=[t(1:p(i)-1) '_' t(p(i)+1:end)];
end

p=strfind(t,'_');
t=['d' int2str(c(1)) '_' int2str(c(2)) '_' int2str(c(3)) t(p(1):end)];
