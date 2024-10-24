% Neuropixel analysis for peripheral nerve stimulation
% Programmed by Akito Kosugi
% v.3.0 07.30.2024

clc

%% Initialization

time_window = [-100,100]; % [ms]
cdpPeakTh = 5;
cdpGain = 0.2;

data_cdp = data_daq(cdpCh,:);

checkCh = 225;
artifactTh = 200;


%% Check matlab version

ver_temp = version('-release');
ver = str2num(ver_temp(1:4));


%% Pre-processing

filtData_AP = data_AP;
passband = [15,3000]/(fs_AP/2);
dataTemp = filtData_AP./(fs_AP/2);
[b,a] = butter(2,passband,'bandpass');
filtData_AP = filtfilt(b,a,dataTemp)*(fs_AP/2);

filtData_LFP = data_LFP;
passband = [1,250]/(fs_LFP/2);
dataTemp = filtData_AP./(fs_LFP/2);
[b,a] = butter(2,passband,'bandpass');
filtData_LFP = filtfilt(b,a,dataTemp)*(fs_LFP/2);

filtData_cdp = dataTemp;
passband = [15,10000]/(fs_daq/2);
dataTemp = data_cdp./(fs_AP/2);
[b,a] = butter(2,passband,'bandpass');
filtData_cdp = filtfilt(b,a,dataTemp)*(fs_daq/2);


%% Epoching

%---Find the synchronization timing---%
syncTrigIdx_AP = find(d_AP_2.Data == syncTrigCh );
syncTrigIdx_daq = find(d_daq_2.Data == syncTrigCh );
syncTrigTime_AP = d_AP_2.Timestamps(syncTrigIdx_AP);
syncTrigTime_daq = d_daq_2.Timestamps(syncTrigIdx_daq);

%---Find the stimulus timing---%
stimTrigIdx = find(d_daq_2.Data == stimTrigCh);
stimTrigTime = d_daq_2.Timestamps(stimTrigIdx);
stimTrigNum = length(stimTrigIdx);

%---Epoching by trigger---%
k = 0;
for n = 2:stimTrigNum-1
    [val,idx] = min(abs(time_daq -stimTrigTime(n)));
    trigIdx_daq = idx;
    
    %---Find the nearest syncTrig---%
    [val,idx] = min(abs(syncTrigTime_daq-stimTrigTime(n)));
    useSyncTrigTime_daq = syncTrigTime_daq(idx);
    addTime = stimTrigTime(n)-useSyncTrigTime_daq;
    [val,idx] = min(abs(syncTrigTime_AP-useSyncTrigTime_daq));
    useSyncTrigTime_AP = syncTrigTime_AP(idx);
    trigTime_AP = useSyncTrigTime_AP+addTime;

    %---Check stim timing---%
    [val,idx] = min(abs(time_AP -trigTime_AP));
    trigIdx_AP = idx;
    checkWindow = ([-2,-0.05]); % [ms]
    if max(abs(filtData_AP(checkCh,trigIdx_AP+checkWindow(1)/1000*fs_AP:trigIdx_AP+checkWindow(2)/1000*fs_AP)))<artifactTh
        checkWindow = ([0.25,1]); % [ms]
        if max(abs(filtData_AP(checkCh,trigIdx_AP+checkWindow(1)/1000*fs_AP:trigIdx_AP+checkWindow(2)/1000*fs_AP)))<artifactTh
            %---Epoching AP by stimTrig---%
            k = k+1;
            data_AP_epoch_trig_temp(:,:,k) = filtData_AP(:,trigIdx_AP+time_window(1)/1000*fs_AP:(trigIdx_AP+time_window(2)/1000*fs_AP)-1,:);
            %---Epoching CDP by stimTrig---%
            [val,idx] = min(abs(time_daq -(stimTrigTime(n))));
            trigIdx_daq = idx;
            data_cdp_epoch_trig_temp(:,k) = filtData_cdp(trigIdx_daq+time_window(1)/1000*fs_daq:(trigIdx_daq+time_window(2)/1000*fs_daq)-1);
        end
    end
end

data_AP_epoch_trig = data_AP_epoch_trig_temp(:,:,1:num_maxTrig);
data_cdp_epoch_trig = data_cdp_epoch_trig_temp(:,1:num_maxTrig);

data_AP_epoch_trig_mean = squeeze(mean(data_AP_epoch_trig(:,:,1:num_maxTrig),3));
cdp_epoch_trig_mean = mean(data_cdp_epoch_trig,2);


%% Spike detection

if spikeDetectionCh > 0
    passband = [600,3000]/(fs_AP/2);
    [b,a] = butter(2,passband);
    baseIdx = 1:(0-time_window(1)-2)*fs_AP/1000;
    spikeThCoef = 5;
    intervalTh = 1; % [ms]
    analysisTime = time_window(1)+1000/fs_AP:1000/fs_AP:time_window(2);

    h= waitbar(0,'Data processing...');
    for ch = 1:num_channels
        X = data_AP_epoch_trig_mean(ch,:);
        spikeAllData(ch).allSpikeTime = [];
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
            ff((0-time_window(1)-1.1)*fs_AP/1000:(0-time_window(1)+1.1)*fs_AP/1000) = 0;

            %---Spike detection---%
            base = ff(baseIdx);
            sigma_n = median(abs(base)/0.6745);
            spikeTh = spikeThCoef*sigma_n;
            if ver < 2022
    %         [pks,locs] = findpeaks(abs(ff),'minPeakHeight',spikeTh,'minPeakDistance',intervalTh*fs_AP/1000);
                [pks,locs] = findpeaks(-ff,'minPeakHeight',spikeTh,'minPeakDistance',intervalTh*fs_AP/1000);
            else
                loc_temp = findpeaks(-ff,spikeTh);
                locs = loc_temp.loc;
                pks = -ff(locs);
            end
            if isempty(pks)
                spikeAllData(ch).spikeData(n).spikeNum = 0;
                spikeAllData(ch).spikeData(n).spikeTh = spikeTh;
                spikeAllData(ch).spikeData(n).locs = 0;
                spikeAllData(ch).spikeData(n).pks = 0;
            else
                spikeAllData(ch).spikeData(n).spikeNum = length(pks);
                spikeAllData(ch).spikeData(n).spikeTh = spikeTh;
                spikeAllData(ch).spikeData(n).locs = locs;
                spikeAllData(ch).spikeData(n).pks = pks;
                spikeAllData(ch).allSpikeTime = horzcat(spikeAllData(ch).allSpikeTime,analysisTime(locs));
            end
            clear B e f ff base pks locs
        end
        waitbar(ch/num_channels,h);
    end
    close(h)
end


%% Signal processsing

%---Amplitude analysis---%
time_window_analysis = [1.5,3]; % [ms]
analysisIdx = (time_window_analysis(1)-time_window(1))*fs_AP/1000+1:(time_window_analysis(2)-time_window(1))*fs_AP/1000;
for ch = 1:num_channels
    amp(ch) = mean(data_AP_epoch_trig_mean(ch,analysisIdx),2);
end

%---CDP analysis---%
time_window = [-100,100]; % [ms]
time_window_analysis = [0.3,1]; % [ms]

analysisIdx = (time_window_analysis(1)-time_window(1))*fs_daq/1000+1:(time_window_analysis(2)-time_window(1))*fs_daq/1000;
if ver < 2022
    [pks,locs] = findpeaks(-cdp_epoch_trig_mean(analysisIdx),'minPeakHeight',cdpPeakTh);
else
    loc_temp = findpeaks(-cdp_epoch_trig_mean(analysisIdx),cdpPeakTh);
    locs = loc_temp.loc;
end
cdpOnset = time_window_analysis(1)+(locs)*1000/fs_daq;


%% Plot

%---Spinal LFP---%
time_window = [-100,100]; % [ms]
plotTime = time_window(1)+1000/fs_AP:1000/fs_AP:time_window(2);
figure
hold on
plot(plotTime,squeeze(data_AP_epoch_trig(checkCh ,:,:)));
xlim([-1,5]);
ylim([-250,250]);
title([num2str(checkCh ) ' ch']);
xlabel('Time [ms]');
ylabel('Amplitude [uV]');
set(gca,'fontsize',15);

pause(1)

for i = 1:columnNum
    figure('position',screenSize)
    if spikeDetectionCh > 0    
        subplot(1,3,1)
    else
        subplot(1,2,1)
    end
    hold on
    iter = 0;
    y_label{1} = 'CDP';
    for ch = i:4:num_channels
        iter = iter+1;
        plot(plotTime,data_AP_epoch_trig_mean(ch,:).*gain+500*iter,'k','linewidth',1);
        y_label{iter+1} = [num2str(ch) ' ch'];
    end
    plot(plotTime,cdp_epoch_trig_mean.*cdpGain-500,'k','linewidth',1);
    plot([0,0],[-1000,500*(iter+4)],'g');
    %     plot([cdpOnset,cdpOnset],[-3000,2500*(iter+4)],'b');
    set(gca,'ytick',[-500,500:500:500*iter]);
    set(gca,'yticklabel',y_label);
    ylim([-1000,500*(iter+1)])
    xlim([-2,10])
    title('Spinal LFP');
    xlabel('Time from stimulus onset [ms]');
    clear y_label

    if spikeDetectionCh > 0    
        subplot(1,3,2)
    else
        subplot(1,2,2)
    end    
    hold on
    plotCh = i:4:num_channels;
    imagesc(plotTime,1:1:iter,data_AP_epoch_trig_mean(plotCh,:));
    axis xy
    colormap('jet')
    caxis([round(-1000/gain),round(1000/gain)]);
    xlim([-2,10])
    ylim([-1,iter+1])
    clear y_label
    iter = 0;
    for ch = plotCh
        iter = iter+1;
        y_label{iter} = [num2str(ch) ' ch'];
    end
    set(gca,'ytick',[1:1:iter]);
    set(gca,'yticklabel',y_label);
    xlabel('Time from stimulus onset [ms]');
    title('Spinal LFP, amplitude');
    clear y_label

    if spikeDetectionCh > 0
        subplot(1,3,3)
        hold on
        iter = 0;
        y_label{1} = 'CDP';
        for ch = plotCh
            iter = iter+1;
            plot(plotTime,squeeze(filtData_AP_epoch_trig(ch,:,:)).*3+500*iter,'k','linewidth',1);
            y_label{iter+1} = [num2str(ch) ' ch'];
        end
        plot(plotTime,cdp_epoch_trig_mean.*cdpGain-500,'k','linewidth',1);
        plot([0,0],[-1000,500*(iter+4)],'g');
        %     plot([cdpOnset,cdpOnset],[-3000,2500*(iter+4)],'b');
        set(gca,'ytick',[-500,500:500:500*(iter)]);
        set(gca,'yticklabel',y_label);
        ylim([-1000,500*(iter+1)])
        xlim([-2,10])
        title('Spike');
        xlabel('Time [ms]');
        clear y_label
    end
    saveas(gcf,[saveName '_LFP_align_' num2str(i) '.fig']);
    saveas(gcf,[saveName '_LFP_align_' num2str(i) '.bmp']);
    pause(2)
end

%---Spike---%
if spikeDetectionCh > 0
    plotTime = time_window(1)+1000/fs_AP:1000/fs_AP:time_window(2);
    ch = spikeDetectionCh;
    for trial = 1:num_maxTrig
        f = figure('position',[screenSize(1)+screenSize(3)*1/10,screenSize(2)+screenSize(4)*1/10,screenSize(3)*1/4,screenSize(4)*2/3]);
        set(f,'name',['stim ' num2str(trial) ', spike ' num2str(spikeAllData(ch).spikeData(trial).spikeNum)])
        subplot(4,1,1)
        hold on
        plot(plotTime,data_AP_epoch_art_mean(ch,:),'r','linewidth',1);
        plot(plotTime,squeeze(data_AP_epoch_art(ch,:,trial)),'k');
        plot([0,0],[-3000,1000],'g');
        xlim([-10,20])
        ylim([-2500,800])
        ylim([-800,800])
        title([num2str(ch) ' ch, stim ' num2str(trial) ', raw data']);
        xlabel('Time [ms]');
        ylabel('Amplitude [uV]');
        subplot(4,1,2)
        hold on
        plot(plotTime,squeeze(detrendData_AP_epoch_art(ch,:,trial)),'k');
        plot([0,0],[-1000,1000],'g');
        xlim([-10,20])
        ylim([-600,600])
        title([num2str(ch) ' ch, stim ' num2str(trial) ', detrend']);
        xlabel('Time [ms]');
        ylabel('Amplitude [uV]');
        subplot(4,1,3)
        hold on
        plot(plotTime,squeeze(filtData_AP_epoch_art(ch,:,trial)),'k');
        plot([0,0],[-1000,1000],'g');
        xlim([-10,20])
        ylim([-600,600])
        title([num2str(ch) ' ch, stim ' num2str(trial) ' filter']);
        xlabel('Time [ms]');
        ylabel('Amplitude [uV]');
        subplot(4,1,4)
        hold on
        plot(plotTime,squeeze(filtData_AP_epoch_art(ch,:,trial)),'k');
        plot([-100,100],[-spikeAllData(ch).spikeData(trial).spikeTh,-spikeAllData(ch).spikeData(trial).spikeTh],'r');
        plot([-100,100],[spikeAllData(ch).spikeData(trial).spikeTh,spikeAllData(ch).spikeData(trial).spikeTh],'r');
        for n = 1:spikeAllData(ch).spikeData(trial).spikeNum
            idx = spikeAllData(ch).spikeData(trial).locs(n);
            scatter(plotTime(idx),filtData_AP_epoch_art(ch,idx,trial),20,'g','filled');
        end
        plot([0,0],[-1000,1000],'g');
        xlim([-10,20]);
        ylim([-600,600])
        title([num2str(ch) ' ch, stim ' num2str(trial) ', spike']);
        xlabel('Time [ms]');
        ylabel('Amplitude [uV]');
        saveas(gcf,[saveName '_spike_ch' num2str(ch) '_stim' num2str(trial) '.fig']);
        saveas(gcf,[saveName '_spike_ch' num2str(ch) '_stim' num2str(trial) '.bmp']);
        pause(5)
        close(gcf);
    end
    
    %---Spike all---%
    f = figure('position',screenSize);
    set(f,'name',['ch ' num2str(ch)])
    for trial = 1:num_maxTrig
        subplot(8,7,trial)
        hold on
        plot(plotTime,squeeze(filtData_AP_epoch_trig(ch,:,trial)),'k');
        plot([-100,100],[-spikeAllData(ch).spikeData(trial).spikeTh,-spikeAllData(ch).spikeData(trial).spikeTh],'r');
        plot([-100,100],[spikeAllData(ch).spikeData(trial).spikeTh,spikeAllData(ch).spikeData(trial).spikeTh],'r');
        for n = 1:spikeAllData(ch).spikeData(trial).spikeNum
            idx = spikeAllData(ch).spikeData(trial).locs(n);
            scatter(plotTime(idx),filtData_AP_epoch_trig(ch,idx,trial),20,'g','filled');
        end
        plot([0,0],[-1500,1500],'g');
        xlim([-10,20]);
        ylim([-1500,1500])
        title([num2str(ch) ' ch, stim ' num2str(trial) ', spike ' num2str(spikeAllData(ch).spikeData(trial).spikeNum)]);
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
    subplot(4,1,1)
    hold on
    for trial = 1:num_maxTrig
        plot(plotTime,squeeze(data_AP_epoch_trig(ch,:,trial)),'k');
        for n = 1:spikeAllData(ch).spikeData(trial).spikeNum
            idx = spikeAllData(ch).spikeData(trial).locs(n);
            scatter(plotTime(idx),data_AP_epoch_trig(ch,idx,trial),20,'g','filled');
        end
    end
    plot(plotTime,data_AP_epoch_trig_mean(ch,:),'r','linewidth',1);
    plot([0,0],[-2500,2500],'g');
    plot([cdpOnset,cdpOnset],[-2500,1000],'b');
    xlim([-10,20]);
    ylim([-2500,1500])
    title([num2str(ch) ' ch, raw data']);
    xlabel('Time [ms]');
    ylabel('Amplitude [uV]');
    subplot(4,1,2)
    hold on
    for trial = 1:num_maxTrig
        plot(plotTime,squeeze(filtData_AP_epoch_trig(ch,:,trial)),'k');
        for n = 1:spikeAllData(ch).spikeData(trial).spikeNum
            idx = spikeAllData(ch).spikeData(trial).locs(n);
            scatter(plotTime(idx),filtData_AP_epoch_trig(ch,idx,trial),20,'g','filled');
        end
    end
    plot([0,0],[-1500,1500],'g');
    plot([cdpOnset,cdpOnset],[-1000,1000],'b');
    xlim([-10,20]);
    ylim([-1500,1500])
    title([num2str(ch) ' ch, filter']);
    xlabel('Time [ms]');
    ylabel('Amplitude [uV]');
    subplot(4,1,3)
    hold on
    plot([0,0],[0,num_maxTrig+5],'g');
    plot([cdpOnset,cdpOnset],[0,num_maxTrig+5],'b');
    for trial = 1:num_maxTrig
        for n = 1:spikeAllData(ch).spikeData(trial).spikeNum
            idx = spikeAllData(ch).spikeData(trial).locs;
            scatter(plotTime(idx),trial,20,'k','filled');
        end
    end
    xlim([-10,20]);
    ylim([-0,num_maxTrig+1]);
    title([num2str(ch) ' ch, raster plot']);
    xlabel('Time [ms]');
    ylabel('Stim No.');
    subplot(4,1,4)
    hold on
    plot([0,0],[0,num_maxTrig+5],'g');
    plot([cdpOnset,cdpOnset],[0,num_maxTrig+5],'b');
    g = histogram(spikeAllData(ch).allSpikeTime,'BinWidth',0.5);
    set(g,'facecolor','w');
    title([num2str(ch) ' ch, PSTH']);
    xlim([-10,20]);
    ylim([0,25]);
    ylim([0,10]);
    xlabel('Time [ms]');
    ylabel('Spike counts');
    saveas(gcf,[saveName '_spike_ch' num2str(ch) '_raster.fig']);
    saveas(gcf,[saveName '_spike_ch' num2str(ch) '_raster.bmp']);
    pause(3)
    close(gcf);
end

%---Amplitude---%
figure
hold on
plot(1:1:num_channels,amp,'k','linewidth',2)
plot([0,385],[0,0],'k--','linewidth',1)
% plot([0,385],[-50,-50],'b','linewidth',1)
% plot([0,385],[50,50],'r','linewidth',1)
title('Average LFP amplitude from 1.5 to 3 ms');
xlabel('Number of channel');
ylabel('Amplitude [uV]');
xlim([0,385]);
ylim([-200,200]);
set(gca,'fontsize',15);

saveas(gcf,[saveName '_amp.fig']);
saveas(gcf,[saveName '_amp.bmp']);


%% save

saveData.profile.fileName = saveName;
saveData.profile.fs = fs_AP;
saveData.data.time = plotTime;
saveData.data.data_AP_epoch_trig = data_AP_epoch_trig;
saveData.data.data_AP_epoch_trig_mean = data_AP_epoch_trig_mean;

save([saveName '_LFP.mat'],'saveData');
