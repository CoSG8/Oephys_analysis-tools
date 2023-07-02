% Neuropixel analysis for lever pulling task
% Programmed by Akito Kosugi
% v.1.1 07.02.2023

clc

%% Initialization

time_window = [-1000,500]; % [ms]
data_lever = data_daq(leverCh,:);


%% Epoching

%---synchronization---%
syncTrigIdx_d1 = find(d1_2.Data == 1);
syncTrigIdx_d3 = find(d3_2.Data == 1);
syncTrigTime_d1 = d1_2.Timestamps(syncTrigIdx_d1);
syncTrigTime_d3 = d3_2.Timestamps(syncTrigIdx_d3);
recStartTime_d1 = d1.Timestamps(1);
recStartTime_d3 = d3.Timestamps(1);

for n = 1:length(syncTrigIdx_d1)
    syncTrigDiff_all(n) = syncTrigTime_d3(n)-syncTrigTime_d1(n);
end

recStartDiff = (recStartTime_d3-syncTrigTime_d3)-(recStartTime_d1-syncTrigTime_d1);

%--- devided by trigger---%
trigIdxTemp = find(d3_2.Data == trigCh);
trigTime = d3_2.Timestamps(trigIdxTemp);
trigIdx = round((trigTime-recStartTime_d3)*fs_AP);
num_trig = length(trigIdx);

%---find the first synchronization timing before stimulation---%
n = 1;
while -1
    trigTime = recStartTime_d3+trigIdx(n)/fs_daq;      
    idx = find(syncTrigTime_d3-trigTime <0);
    if isempty(idx) == 0
        break
    else
        n = n+1;
    end
end
firstTrig = n;

%---AP---%
k = 0;
for n = firstTrig:num_maxTrig+firstTrig-1
    k = k+1;
    % find the nearest synchronization timing
    trigTime = recStartTime_d3+trigIdx(n)/fs_daq;      
    idxTemp = find(syncTrigTime_d3-trigTime <0);
    idx = max(idxTemp);
    delay = recStartDiff-syncTrigDiff_all(idx)-n*(1/fs_AP); %[s]
    temp = trigIdx(n)+floor(delay*fs_AP)+fs_AP*time_window(1)/1000+1:trigIdx(n)+floor(delay*fs_AP)+fs_AP*time_window(2)/1000;
    data_AP_epoch_trig(:,:,k) = filtData_AP(:,temp);
end
data_AP_epoch_trig_mean = squeeze(mean(data_AP_epoch_trig,3));

%---Lever---%
k = 0;
for n = firstTrig:num_maxTrig+firstTrig-1
    k = k+1;
    temp = trigIdx(n)+fs_daq*time_window(1)/1000+1:trigIdx(n)+fs_daq*time_window(2)/1000;
    lever_epoch_trig(:,k) = data_lever(temp);
end
lever_epoch_trig_mean = squeeze(mean(lever_epoch_trig,2));


%% Spike detection

passband = [600,3000]/(fs_AP/2);
[b,a] = butter(2,passband);
baseTime_ms = 100;
baseIdx = 1:baseTime_ms*fs_AP/1000;
spikeThCoef = 5;
intervalTh = 0.8; % [ms]
analysisTime = time_window(1)+1000/fs_AP:1000/fs_AP:time_window(2);

for ch = spikeDetectionCh
    X = data_AP_epoch_trig_mean(ch,:);
    spikeAllData.allSpikeTime = [];
    for n = 1:num_maxTrig
        %---Detrend---%
        Y = squeeze(data_AP_epoch_trig(ch,:,n));
        B = X\Y;
        e = Y-X*B;
        e = Y-X;
        detrendData_AP_epoch_trig(ch,:,n) = e;

        %---Filter---%
        f = e./(fs_AP/2);
        ff = filtfilt(b,a,f)*(fs_AP/2);
        filtData_AP_epoch_trig(ch,:,n) = ff;
    
        %---Spike detection---%
        base = ff(baseIdx);
        sigma_n = median(abs(base)/0.6745);
        spikeTh = spikeThCoef*sigma_n;

%         [pks,locs] = findpeaks(abs(ff),'minPeakHeight',spikeTh,'minPeakDistance',intervalTh*fs_AP/1000);
        [pks,locs] = findpeaks(-ff,'minPeakHeight',spikeTh,'minPeakDistance',intervalTh*fs_AP/1000);
        if isempty(pks)
            spikeAllData.spikeData(n).spikeNum = 0;
            spikeAllData.spikeData(n).spikeTh = spikeTh;
            spikeAllData.spikeData(n).locs = 0;
            spikeAllData.spikeData(n).pks = 0;
        else
            spikeAllData.spikeData(n).spikeNum = length(pks);
            spikeAllData.spikeData(n).spikeTh = spikeTh;
            spikeAllData.spikeData(n).locs = locs;
            spikeAllData.spikeData(n).pks = pks;
            spikeAllData.allSpikeTime = horzcat(spikeAllData.allSpikeTime,analysisTime(locs));
        end
        clear B e f ff base pks locs
    end
end


%% Plot

%---Spike---%
plotTime = time_window(1)+1000/fs_AP:1000/fs_AP:time_window(2);
ch = spikeDetectionCh;
% for trial = 1:num_maxTrig
%     f = figure('position',[screenSize(1)+screenSize(3)*1/10,screenSize(2)+screenSize(4)*1/10,screenSize(3)*1/4,screenSize(4)*2/3]);
%     set(f,'name',['stim ' num2str(trial) ', spike ' num2str(spikeAllData(ch).spikeData(trial).spikeNum)])
%     subplot(4,1,1)
%     hold on
%     plot(plotTime,data_AP_epoch_art_mean(ch,:),'r','linewidth',1);
%     plot(plotTime,squeeze(data_AP_epoch_art(ch,:,trial)),'k');
%     plot([0,0],[-3000,1000],'g');
%     xlim([-10,20])
%     ylim([-2500,800])
% %     ylim([-800,800])
%     title([num2str(ch) ' ch, stim ' num2str(trial) ', raw data']);
%     xlabel('Time [ms]');
%     ylabel('Amplitude [uV]');
%     subplot(4,1,2)
%     hold on
%     plot(plotTime,squeeze(detrendData_AP_epoch_art(ch,:,trial)),'k');
%     plot([0,0],[-1000,1000],'g');
%     xlim([-10,20])
%     ylim([-600,600])
%     title([num2str(ch) ' ch, stim ' num2str(trial) ', detrend']);
%     xlabel('Time [ms]');
%     ylabel('Amplitude [uV]');
%     subplot(4,1,3)
%     hold on
%     plot(plotTime,squeeze(filtData_AP_epoch_art(ch,:,trial)),'k');
%     plot([0,0],[-1000,1000],'g');
%     xlim([-10,20])
%     ylim([-600,600])
%     title([num2str(ch) ' ch, stim ' num2str(trial) ' filter']);
%     xlabel('Time [ms]');
%     ylabel('Amplitude [uV]');
%     subplot(4,1,4)
%     hold on
%     plot(plotTime,squeeze(filtData_AP_epoch_art(ch,:,trial)),'k');
%     plot([-100,100],[-spikeAllData(ch).spikeData(trial).spikeTh,-spikeAllData(ch).spikeData(trial).spikeTh],'r');
%     plot([-100,100],[spikeAllData(ch).spikeData(trial).spikeTh,spikeAllData(ch).spikeData(trial).spikeTh],'r');
%     for n = 1:spikeAllData(ch).spikeData(trial).spikeNum
%         idx = spikeAllData(ch).spikeData(trial).locs(n);
%         scatter(plotTime(idx),filtData_AP_epoch_art(ch,idx,trial),20,'g','filled');
%     end
%     plot([0,0],[-1000,1000],'g');
%     xlim([-10,20]);
%     ylim([-600,600])
%     title([num2str(ch) ' ch, stim ' num2str(trial) ', spike']);
%     xlabel('Time [ms]');
%     ylabel('Amplitude [uV]');
% %     saveas(gcf,[saveName '_spike_ch' num2str(ch) '_stim' num2str(trial) '.fig']);
% %     saveas(gcf,[saveName '_spike_ch' num2str(ch) '_stim' num2str(trial) '.bmp']);
%     pause(5)
%     close(gcf);
% end

%---Spike all---%
f = figure('position',screenSize);
set(f,'name',['ch ' num2str(ch)])
for trial = 1:num_maxTrig
    subplot(8,7,trial)
    hold on
    plot(plotTime,squeeze(filtData_AP_epoch_trig(ch,:,trial)),'k');
    plot([time_window(1),time_window(2)],[-spikeAllData.spikeData(trial).spikeTh,-spikeAllData.spikeData(trial).spikeTh],'r');
    plot([time_window(1),time_window(2)],[spikeAllData.spikeData(trial).spikeTh,spikeAllData.spikeData(trial).spikeTh],'r');
    for n = 1:spikeAllData.spikeData(trial).spikeNum
        idx = spikeAllData.spikeData(trial).locs(n);
        scatter(plotTime(idx),filtData_AP_epoch_trig(ch,idx,trial),20,'g','filled');
    end
    plot([0,0],[-1000,1000],'g');
    xlim([-1000,500]);
    ylim([-600,600])
    title([num2str(ch) ' ch, trig ' num2str(trial) ', spike ' num2str(spikeAllData.spikeData(trial).spikeNum)]);
end
xlabel('Time [ms]');
ylabel('Amplitude [uV]');
saveas(gcf,[saveName '_spike_ch' num2str(ch) '_all.fig']);
saveas(gcf,[saveName '_spike_ch' num2str(ch) '_all.bmp']);
pause(3)
close(gcf);


%---Rasrer plot---%
f = figure('position',[screenSize(1)+screenSize(3)*1/10,screenSize(2)+screenSize(4)*1/10,screenSize(3)*1/4,screenSize(4)*2/3]);
set(f,'name',['ch ' num2str(ch)])
subplot(5,1,1)
hold on
for trial = 1:num_maxTrig
    plot(plotTime,squeeze(data_AP_epoch_trig(ch,:,trial)),'k');
    for n = 1:spikeAllData.spikeData(trial).spikeNum
        idx = spikeAllData.spikeData(trial).locs(n);
        scatter(plotTime(idx),data_AP_epoch_trig(ch,idx,trial),20,'g','filled');
    end
end
plot(plotTime,data_AP_epoch_trig_mean(ch,:),'r','linewidth',1);
plot([0,0],[-2500,2500],'g');
xlim([-1000,500]);
ylim([-2000,2000])
% ylim([-800,800])
title([num2str(ch) ' ch, raw data']);
xlabel('Time [ms]');
ylabel('Amplitude [uV]');
subplot(5,1,2)
hold on
for trial = 1:num_maxTrig
    plot(plotTime,squeeze(filtData_AP_epoch_trig(ch,:,trial)),'k');
    for n = 1:spikeAllData.spikeData(trial).spikeNum
        idx = spikeAllData.spikeData(trial).locs(n);
        scatter(plotTime(idx),filtData_AP_epoch_trig(ch,idx,trial),20,'g','filled');
    end
end
plot([0,0],[-1000,1000],'g');
xlim([-1000,500]);
ylim([-600,600])
title([num2str(ch) ' ch, filter']);
xlabel('Time [ms]');
ylabel('Amplitude [uV]');
subplot(5,1,3)
hold on
plot([0,0],[0,num_maxTrig+5],'g');
for trial = 1:num_maxTrig
    for n = 1:spikeAllData.spikeData(trial).spikeNum
        idx = spikeAllData.spikeData(trial).locs;
        scatter(plotTime(idx),trial,20,'k','filled');
    end
end
xlim([-1000,500]);
ylim([-0,num_maxTrig+1]);
title([num2str(ch) ' ch, raster plot']);
xlabel('Time [ms]');
ylabel('Stim No.');
subplot(5,1,4)
hold on
plot([0,0],[0,num_maxTrig+5],'g');
g = histogram(spikeAllData.allSpikeTime,'BinWidth',0.5);
set(g,'facecolor','w');
title([num2str(ch) ' ch, PETH']);
xlim([-1000,500]);
% ylim([0,25]);
ylim([0,10]);
xlabel('Time [ms]');
ylabel('Spike counts');
subplot(5,1,5)
hold on
plot(plotTime,lever_epoch_trig,'k');
title('Lever');
xlim([-1000,500]);
% ylim([0,25]);
ylim([0,12000]);
xlabel('Time [ms]');
ylabel('Lever position');
saveas(gcf,[saveName '_spike_ch' num2str(ch) '_raster.fig']);
saveas(gcf,[saveName '_spike_ch' num2str(ch) '_raster.bmp']);
pause(3)
close(gcf);
