% startup.m
% sets local defaults
global crlf
%extra_paths='/Users/morse/Documents/projects/HigleyLab/GidonSegev_w_InhibSpines/Matlab';	% replace this path to include your own scripts folder
%extra_paths=['
extra_paths='/home/tmm46/matlab';
extra_paths1 = '/Users/morse/Documents/Matlab/my';
path(path,extra_paths)
%path(path,extra_paths1)

crlf=[char(13) char(10)];
filename=[datestr(clock,30) '.txt'];

% uncomment and modify the below lines if you would like to set your current directory
% each time you start
% directory='C:\Documents and Settings\Thomas Morse\My Documents\MATLAB\';
% cd(directory)

format shortG
format compact
disp('"format shortG" and "format compact" were executed')

%directory='/Users/morse/Documents/projects/'; % HigleyLab/GidonSegev_w_InhibSpines/'
%cd(directory)

disp(['Working directory currently set to ' crlf pwd]);

% If you want to use a diaries subdirectory like I do uncomment below
% otherwise the diaries are placed at top level
%cmd=['diary /Users/morse/Documents/MATLAB/diaries/' filename];
cmd=['diary /home/tmm46/matlab/diaries/' filename];
% cmd=['diary ' filename]; % comment this line out if you uncomment the above line
eval(cmd);

diary on
disp(['Begining of diary ' filename crlf ...
    'Welcome to matlab configured by Tom with extra paths:' crlf extra_paths crlf extra_paths1]) % change Tom to your name!

% more on	% can be left out.  The advantage of this is that you can quit a more with a q if
	% if you accidentally leave a semicolon off a matrix equation and are getting a 
	% print out that would bloat your diary file.
%edit	% can be left out.  I personally like to start editing something 99% of the time
	% so I automatically start up the matlab editor
