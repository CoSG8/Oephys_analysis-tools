% Neuropixel analysis for peripheral nerve stimulation
% Programmed by Akito Kosugi
% v.1.4 07.02.2023

clc

%% Initialization

time_window = [-100,100]; % [ms]
trigTh = 100;
cdpPeakTh = 5;
cdpGain = 0.2;

trig = data_daq(trigCh,:);
data_cdp = data_daq(cdpCh,:);


%% Pre-processing

filtData_cdp = dataTemp;
passband = [15,10000]/(fs_daq/2);
dataTemp = data_cdp./(fs_AP/2);
[b,a] = butter(2,passband,'bandpass');
filtData_cdp = filtfilt(b,a,dataTemp)*(fs_daq/2);


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

%---devided by trigger---%
iter = 0;
for i = 2:length(trig)
    if trig(i) > trigTh && trig(i-1) < trigTh
        iter = iter+1;
        trigIdx(iter) = i;
    end
end
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

%---Epoching using the nearest synchronization timing---%
k = 0;
for n = firstTrig:num_maxTrig+firstTrig-1

    trigTime = recStartTime_d3+trigIdx(n)/fs_daq;      
    idxTemp = find(syncTrigTime_d3-trigTime <0);
    idx = max(idxTemp);
    if syncTrigDiff_all(idx)<0.0002
        k = k+1;
        %---AP---%
        delay = recStartDiff-syncTrigDiff_all(idx)-n*(1/fs_AP); %[s]
        temp = trigIdx(n)+floor(delay*fs_AP)+fs_AP*time_window(1)/1000+1:trigIdx(n)+floor(delay*fs_AP)+fs_AP*time_window(2)/1000;
        data_AP_epoch_trig(:,:,k) = filtData_AP(:,temp);
        %---CDP---%
        temp = trigIdx(n)+fs_daq*time_window(1)/1000+1:trigIdx(n)+fs_daq*time_window(2)/1000;
        cdp_epoch_trig(:,k) = filtData_cdp(temp);
    end
end
data_AP_epoch_trig_mean = squeeze(mean(data_AP_epoch_trig,3));
cdp_epoch_trig_mean = squeeze(mean(cdp_epoch_trig,2));

num_maxTrig = k;


%% Spike detection

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
%         [pks,locs] = findpeaks(abs(ff),'minPeakHeight',spikeTh,'minPeakDistance',intervalTh*fs_AP/1000);
        [pks,locs] = findpeaks(-ff,'minPeakHeight',spikeTh,'minPeakDistance',intervalTh*fs_AP/1000);
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


%% CDP analysis

time_window = [-100,100]; % [ms]
time_window_analysis = [0.3,1]; % [ms]

analysisIdx = (time_window_analysis(1)-time_window(1))*fs_daq/1000+1:(time_window_analysis(2)-time_window(1))*fs_daq/1000;
[pks,locs] = findpeaks(-cdp_epoch_trig_mean(analysisIdx),'minPeakHeight',cdpPeakTh);
cdpOnset = time_window_analysis(1)+(locs)*1000/fs_daq;


%% Plot

%---Spinal LFP---%
time_window = [-100,100]; % [ms]
plotTime = time_window(1)+1000/fs_AP:1000/fs_AP:time_window(2);

switch columnIdx
    case 1
        for i = 1:4
            figure('position',screenSize)
            subplot(1,3,1)
            hold on
            iter = 0;
            y_label{1} = 'CDP';
            for ch = i:4:384
                iter = iter+1;
                plot(plotTime,data_AP_epoch_trig_mean(ch,:).*gain+2500*iter,'k','linewidth',1);
                y_label{iter+1} = [num2str(ch) ' ch'];
            end
            plot(plotTime,cdp_epoch_trig_mean.*cdpGain-2500,'k','linewidth',1);
            plot([0,0],[-5000,2500*(iter+4)],'g');
        %     plot([cdpOnset,cdpOnset],[-3000,2500*(iter+4)],'b');
            set(gca,'ytick',[-2500,2500:2500:2500*iter]);    
            set(gca,'yticklabel',y_label);
            ylim([-5000,2500*(iter+1)])
            xlim([-2,10])
            title('Spinal LFP');
            xlabel('Time [ms]');
            clear y_label

            subplot(1,3,2)
            plotCh = i:4:384;
            imagesc(plotTime,1:1:iter,data_AP_epoch_trig_mean(plotCh,:));
            axis xy
            colormap('jet')
            caxis([round(-2500/gain),round(2500/gain)]);
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
            xlabel('Time [ms]');
            title('Spinal LFP, amplitude');
            clear y_label

            subplot(1,3,3)
            hold on
            iter = 0;
            y_label{1} = 'CDP';
            for ch = plotCh
                iter = iter+1;
                plot(plotTime,squeeze(filtData_AP_epoch_trig(ch,:,:)).*3+2500*iter,'k','linewidth',1);
                y_label{iter+1} = [num2str(ch) ' ch'];
            end
            plot(plotTime,cdp_epoch_trig_mean.*cdpGain-2500,'k','linewidth',1);
            plot([0,0],[-5000,2500*(iter+4)],'g');
        %     plot([cdpOnset,cdpOnset],[-3000,2500*(iter+4)],'b');
            set(gca,'ytick',[-2500,2500:2500:2500*(iter)]);    
            set(gca,'yticklabel',y_label);
            ylim([-5000,2500*(iter+1)])
            xlim([-2,10])
            title('Spike');
            xlabel('Time [ms]');
            clear y_label
            
        saveas(gcf,[saveName '_LFP_align_' num2str(i) '.fig']);
        saveas(gcf,[saveName '_LFP_align_' num2str(i) '.bmp']);
        pause(1)
        close(gcf);                    
        end
        
    case 2      
        for i = 1:2
            figure('position',screenSize)
            subplot(1,3,1)
            hold on
            iter = 0;
            y_label{1} = 'CDP';
            switch i
                case 1
                for ch = 1:4:381
                    iter = iter+1;
                    plot(plotTime,data_AP_epoch_trig_mean(ch,:).*gain+1200*iter,'k','linewidth',1);
                    y_label{iter+1} = [num2str(ch) ' ch'];  
                end
                for ch = 3:4:383
                    iter = iter+1;
                    plot(plotTime,data_AP_epoch_trig_mean(ch,:).*gain+1200*iter,'k','linewidth',1);
                    y_label{iter+1} = [num2str(ch) ' ch'];
                end                    
                case 2
                for ch = 4:4:384
                    iter = iter+1;
                    plot(plotTime,data_AP_epoch_trig_mean(ch,:).*gain+1200*iter,'k','linewidth',1);
                    y_label{iter+1} = [num2str(ch) ' ch'];
                end                                        
                for ch = 2:4:382
                    iter = iter+1;
                    plot(plotTime,data_AP_epoch_trig_mean(ch,:).*gain+1200*iter,'k','linewidth',1);
                    y_label{iter+1} = [num2str(ch) ' ch'];
                end                    
            end
            plot(plotTime,cdp_epoch_trig_mean.*cdpGain-1200,'k','linewidth',1);
            plot([0,0],[-5000,1200*(iter+4)],'g');
        %     plot([cdpOnset,cdpOnset],[-3000,2500*(iter+4)],'b');
            set(gca,'ytick',[-1200,1200:1200:1200*iter]);    
            set(gca,'yticklabel',y_label);
            ylim([-5000,1200*(iter+1)])
            xlim([-2,10])
            title('Spinal LFP');
            xlabel('Time [ms]');
            clear y_label

            subplot(1,3,2)
            hold on
            iter = 0;
            switch i
                case 1
                for ch = 1:4:381
                    iter = iter+1;
                    imagesc(plotTime,iter,data_AP_epoch_trig_mean(ch,:));
                    y_label{iter} = [num2str(ch) ' ch'];
                end
                for ch = 3:4:383
                    iter = iter+1;
                    imagesc(plotTime,iter,data_AP_epoch_trig_mean(ch,:));
                    y_label{iter} = [num2str(ch) ' ch'];
                end
                case 2
                for ch = 4:4:384
                    iter = iter+1;
                    imagesc(plotTime,iter,data_AP_epoch_trig_mean(ch,:));
                    y_label{iter} = [num2str(ch) ' ch'];
                end
                for ch = 2:4:382
                    iter = iter+1;
                    imagesc(plotTime,iter,data_AP_epoch_trig_mean(ch,:));
                    y_label{iter} = [num2str(ch) ' ch'];
                end  
            end
            axis xy
            colormap('jet')
            caxis([round(-1200/gain),round(1200/gain)]);
            xlim([-2,10])
            ylim([-1,iter+1])
            set(gca,'ytick',[1:1:iter]);    
            set(gca,'yticklabel',y_label);
            xlabel('Time [ms]');
            title('Spinal LFP, amplitude');
            clear y_label

            subplot(1,3,3)
            hold on
            iter = 0;
            y_label{1} = 'CDP';
            switch i
                case 1
                    for ch = 1:4:381            
                        iter = iter+1;
                        plot(plotTime,squeeze(filtData_AP_epoch_trig(ch,:,:)).*3+1200*iter,'k','linewidth',1);
                        y_label{iter+1} = [num2str(ch) ' ch'];
                    end
                    for ch = 3:4:383
                        iter = iter+1;
                        plot(plotTime,squeeze(filtData_AP_epoch_trig(ch,:,:)).*3+1200*iter,'k','linewidth',1);
                        y_label{iter+1} = [num2str(ch) ' ch'];
                    end
                case 2
                    for ch = 4:4:384         
                        iter = iter+1;
                        plot(plotTime,squeeze(filtData_AP_epoch_trig(ch,:,:)).*3+1200*iter,'k','linewidth',1);
                        y_label{iter+1} = [num2str(ch) ' ch'];
                    end
                    for ch = 2:4:382
                        iter = iter+1;
                        plot(plotTime,squeeze(filtData_AP_epoch_trig(ch,:,:)).*3+1200*iter,'k','linewidth',1);
                        y_label{iter+1} = [num2str(ch) ' ch'];
                    end
            end
            plot(plotTime,cdp_epoch_trig_mean.*cdpGain-1200,'k','linewidth',1);
            plot([0,0],[-5000,1200*(iter+4)],'g');
        %     plot([cdpOnset,cdpOnset],[-3000,2500*(iter+4)],'b');
            set(gca,'ytick',[-1200,1200:1200:1200*(iter)]);    
            set(gca,'yticklabel',y_label);
            ylim([-5000,1200*(iter+1)])
            xlim([-2,10])
            title('Spike');
            xlabel('Time [ms]');
            clear y_label
            
            saveas(gcf,[saveName '_LFP_align_' num2str(i) '.fig']);
            saveas(gcf,[saveName '_LFP_align_' num2str(i) '.bmp']);
            pause(1)
            close(gcf);
        end
end


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
% plot([cdpOnset,cdpOnset],[-2500,1000],'b');
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
% plot([cdpOnset,cdpOnset],[-1000,1000],'b');
xlim([-10,20]);
ylim([-1500,1500])
title([num2str(ch) ' ch, filter']);
xlabel('Time [ms]');
ylabel('Amplitude [uV]');
subplot(4,1,3)
hold on
plot([0,0],[0,num_maxTrig+5],'g');
% plot([cdpOnset,cdpOnset],[0,num_maxTrig+5],'b');
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
% plot([cdpOnset,cdpOnset],[0,num_maxTrig+5],'b');
g = histogram(spikeAllData(ch).allSpikeTime,'BinWidth',0.5);
set(g,'facecolor','w');
title([num2str(ch) ' ch, PSTH']);
xlim([-10,20]);
% ylim([0,25]);
ylim([0,10]);
xlabel('Time [ms]');
ylabel('Spike counts');
saveas(gcf,[saveName '_spike_ch' num2str(ch) '_raster.fig']);
saveas(gcf,[saveName '_spike_ch' num2str(ch) '_raster.bmp']);
pause(3)
close(gcf);
