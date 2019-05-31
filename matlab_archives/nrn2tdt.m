% nrn2dtd.m
% Converts data files written from NEURON into a tdt2mat_data tank
% and saves it in tdt2mat_data_yyyymmdd.mat
% usage:
% tdt2mat_data = nrn2tdt(foldername)
% sets tdt2mat_data and writes a timestamped tank in foldername
% This version
% tdt2mat_data = nrn2tdt(foldername,'no write tank')
% sets tdt2mat_data and but does not write a tank (supply any second
% argument to not write a tank)

function tdt2mat_data = nrn2tdt(foldername,varargin)
% nrn2tdt(foldername) is passed a folder where there are the following
% files written by NEURON:
% _d_epocs_d_SOFF_d_onset.dat	_d_snips_d_eNeu_d_ts.dat
% _d_epocs_d_S_ON_d_data.dat	_d_streams_d_BRTH_d_data.dat
% _d_epocs_d_S_ON_d_onset.dat	_d_streams_d_BRTH_d_fs.dat
% _d_snips_d_eNeu_d_sortcode.dat
% These files are read in and the _d_ is a guide to map these to the
% structure "." in matlab
% The tdt2mat_data is then written with a time stamp to the same folder
% and returned from the function
here=pwd;
cd(foldername)
load _d_epocs_d_SOFF_d_onset.dat;
load _d_snips_d_eNeu_d_ts.dat;
load _d_epocs_d_S_ON_d_data.dat;
load _d_streams_d_BRTH_d_data.dat;
load _d_epocs_d_S_ON_d_onset.dat;
load _d_streams_d_BRTH_d_fs.dat;
load _d_snips_d_eNeu_d_sortcode.dat;
clear tdt2mat_data
tdt2mat_data.epocs.SOFF.onset = X_d_epocs_d_SOFF_d_onset;
tdt2mat_data.snips.eNeu.ts = X_d_snips_d_eNeu_d_ts;
tdt2mat_data.epocs.S_ON.data = X_d_epocs_d_S_ON_d_data;
tdt2mat_data.streams.BRTH.data = X_d_streams_d_BRTH_d_data;
tdt2mat_data.epocs.S_ON.onset=X_d_epocs_d_S_ON_d_onset;
tdt2mat_data.streams.BRTH.fs = X_d_streams_d_BRTH_d_fs;
tdt2mat_data.snips.eNeu.sortcode = X_d_snips_d_eNeu_d_sortcode;
if nargin==1
    filename=['tdt2mat_' timestamp '.mat'];
    save(filename,'tdt2mat_data');
end
cd(here); % go back to folder came from
end
