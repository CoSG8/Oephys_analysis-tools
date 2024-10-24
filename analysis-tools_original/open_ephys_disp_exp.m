% Programmed by Akito Kosugi
% v.3.0 07.30.2024

clear all
close all
clc


%% Initialization

screenSize = get(0,'screensize');
NPName_AP = 'AP';
NPName_LFP = 'LFP';
daqName = 'NI-DAQmx';


%% Set GUI

open_ephys_disp_gui;


%% Data loading

%---search .oebin file---%
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

%---file loading---%
list_continuous = list_open_ephys_binary(filepath,'continuous');
for i = 1:length(list_continuous)
    if (contains(char(list_continuous{i}),daqName)) == 1
        daqIdx = i;
    elseif (contains(char(list_continuous{i}),NPName_AP)) == 1
        NPIdx_AP = i;
    elseif (contains(char(list_continuous{i}),NPName_LFP)) == 1
        NPIdx_LFP = i;
    end
end

if uiIndex > 1
    d_AP = load_open_ephys_binary(filepath,'continuous',NPIdx_AP); % AP
    fs_AP = d_AP.Header.sample_rate;
    num_channels = d_AP.Header.num_channels;
    bit_volts_AP = d_AP.Header.channels(1).bit_volts;
    data_AP = d_AP.Data.*bit_volts_AP;
    time_AP = d_AP.Timestamps;
    d_AP_2 = load_open_ephys_binary(filepath,'events',NPIdx_AP); % AP
    disp('AP data loading complete');
    
    d_LFP = load_open_ephys_binary(filepath,'continuous',NPIdx_LFP); % LFP
    fs_LFP = d_LFP.Header.sample_rate;
    num_channels = d_LFP.Header.num_channels;
    bit_volts_LFP = d_LFP.Header.channels(1).bit_volts;
    data_LFP = d_LFP.Data.*bit_volts_LFP;
    time_LFP = d_LFP.Timestamps;
    disp('LFP data loading complete');
end

d_daq = load_open_ephys_binary(filepath,'continuous',daqIdx); % DAQ
bit_volts_daq = d_daq.Header.channels(1).bit_volts;
fs_daq = d_daq.Header.sample_rate;
data_daq =  d_daq.Data.*bit_volts_daq;
time_daq = d_daq.Timestamps;

list_event = list_open_ephys_binary(filepath,'events');
for i = 1:length(list_event)
    if (contains(char(list_event{i}),daqName)) == 1
        daqIdx = i;
    elseif (contains(char(list_event{i}),'MessageCenter')) == 1
        eventIdx = i;
    end        
end
d_daq_2 = load_open_ephys_binary(filepath,'events',daqIdx); % DAQ
disp('DAQ data loading complete');
d_event = load_open_ephys_binary(filepath,'events',eventIdx); % Wordinput
disp('Event data loading complete');


%% Analysis and plot

switch uiIndex
    case 1
        open_ephys_disp_trig;        
    case 2    
        open_ephys_disp_resting;
    case 3
        open_ephys_disp_LFP;
end
