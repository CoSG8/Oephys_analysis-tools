% Programmed by Akito Kosugi
% v.1.5.1 07.07.2023

clear all
close all
clc


%% Initialization

screenSize = get(0,'screensize');


%% Set GUI

open_ephys_disp_gui;


%% Data loading

correntFolder = pwd;
[fName,dPath] = uigetfile('*.oebin');
filepath = [dPath,fName];
if PCIdx == 1
    idx= findstr(dPath,'\');
else
    idx= findstr(dPath,'/');
end
idxNum = length(idx);
saveName = [dPath(idx(idxNum-4)+1:idx(idxNum-3)-1) '_' dPath(idx(idxNum-3)+1:idx(idxNum-2)-1) '_' dPath(idx(idxNum-2)+1:idx(idxNum-1)-1) '_' dPath(idx(idxNum-1)+1:idx(idxNum)-1)];

%---check the number of NeuroPixels---%
dataPath = [dPath(1:idx(idxNum)) 'continuous'];
dataNum = length(dir(dataPath))-2;
daqIdx = dataNum;
switch probeIdx
    case 1
        NPIdx_AP = 1;
        NPIdx_LFP = 2;
    case 2
        NPIdx_AP = 3;
        NPIdx_LFP = 4;
    case 3
        NPIdx_AP = 5;
        NPIdx_LFP = 6;
end

d1 = load_open_ephys_binary(filepath,'continuous',NPIdx_AP); % AP
data_AP = d1.Data;
fs_AP = d1.Header.sample_rate;
num_channels = d1.Header.num_channels;
d1_2 = load_open_ephys_binary(filepath,'events',NPIdx_AP); % AP
disp('AP data loading complete');

d2 = load_open_ephys_binary(filepath,'continuous',NPIdx_LFP); % LFP
data_LFP = d2.Data;
fs_LFP = d2.Header.sample_rate;
num_channels = d2.Header.num_channels;
disp('LFP data loading complete');

d3 = load_open_ephys_binary(filepath,'continuous',daqIdx); % DAQ
data_daq =  d3.Data;
fs_daq = d3.Header.sample_rate;
d3_2 = load_open_ephys_binary(filepath,'events',daqIdx); % DAQ
disp('DAQ data loading complete');


%% Pre-processing

% time_window = [0,180]; % [s]
% data_AP = data_AP(:,time_window(1)*fs_AP+1:time_window(2)*fs_AP);
% data_LFP = data_AP(:,time_window(1)*fs_LFP+1:time_window(2)*fs_LFP);
% data_daq = data_daq(:,time_window(1)*fs_daq+1:time_window(2)*fs_daq);

filtData_AP = data_AP;
passband = [15,10000]/(fs_AP/2);
dataTemp = data_AP./(fs_AP/2);
[b,a] = butter(2,passband,'bandpass');
filtData_AP = filtfilt(b,a,dataTemp)*(fs_AP/2);

filtData_LFP = data_LFP;
passband = [1,250]/(fs_LFP/2);
dataTemp = data_LFP./(fs_LFP/2);
[b,a] = butter(2,passband,'bandpass');
filtData_LFP = filtfilt(b,a,dataTemp)*(fs_LFP/2);


%% Analysis and plot

switch uiIndex
    case 1
        open_ephys_disp_resting;
    case 2
        open_ephys_disp_LFP;
    case 3
        open_ephys_disp_spike_task;        
end
