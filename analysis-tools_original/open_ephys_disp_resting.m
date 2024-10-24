% Neuropixel analysis for noiose check
% Programmed by Akito Kosugi
% v3.0 07.30.2024

close all


%% Initialization

time_window = [0,analysisTime]; %[s]
resol_fft = 1;
analysisFreq = 50;


%% Epoching

%--- devided by time---%
plotIdx_rest = time_window(1)*fs_AP+1:time_window(2)*fs_AP;
data_AP_epoch_time = data_AP(:,plotIdx_rest);


%% Resting analysis

%---FFT---%
noverlap = 0;
nfft = fs_AP/resol_fft;
w = hanning(nfft);
[power,freq] = pwelch(detrend(data_AP_epoch_time)',w,noverlap,nfft,fs_AP);
normPower = power./sum(power,1).*100;

rms = mean(data_AP_epoch_time.^2,2).^0.5;


%% Plot

%---Resting---%
time_window = [analysisTime-2,analysisTime-1]; %[s]
plotTime = time_window(1):1/fs_AP:time_window(2)-1/fs_AP;
plotIdx = time_window(1)*fs_AP:time_window(2)*fs_AP-1;


for i = 1:columnNum
    figure('position',screenSize)
    subplot(1,3,1)
    hold on
    iter = 0;
    for ch = i:4:num_channels
        iter = iter+1;
        plot(plotTime,data_AP_epoch_time(ch,plotIdx).*gain+1000*iter,'k','linewidth',1);
        y_label{iter} = [num2str(ch) ' ch'];
    end
    set(gca,'ytick',[1000:1000:1000*iter]);
    set(gca,'yticklabel',y_label);
    ylim([-1000,1000*(iter+1)])
    xlim([time_window(1),time_window(2)])
    title('Raw  data');
    xlabel('Time [ms]')
    
    subplot(1,3,2)
    hold on
    iter = 0;
    for ch = i:4:num_channels
        iter = iter+1;
        if optionIdx == 1
            plot(freq,normPower(:,ch).*gain+1*iter,'k','linewidth',1);
        else
            plot(freq,power(:,ch).*gain+20*iter,'k','linewidth',1);
        end
        y_label{iter} = [num2str(ch) ' ch'];
    end
    plot([50,50],[-1000,10000],'r');
    plot([100,100],[-1000,10000],'r');
    xlim([0,110])
    if optionIdx == 1
        title('Normalized power spectrum density');
        ylim([-1,1*(iter+1)])
        set(gca,'ytick',[1:1:1*iter]);
    else
        title('Power spectrum density');
        ylim([-20,20*(iter+1)])
        set(gca,'ytick',[20:20:20*iter]);
    end
    set(gca,'yticklabel',y_label);
    xlabel('Frequency [Hz]');

    subplot(1,12,10)
    plotCh = i:4:num_channels;
    imagesc(1,1:1:iter,rms(plotCh,:));
    axis xy
    colormap('jet')
    caxis([10,100/gain]);
    ylim([-1,iter+1])
    clear y_label
    iter = 0;
    for ch = plotCh
        iter = iter+1;
        y_label{iter} = [num2str(ch) ' ch'];
    end
    title('Root mean square');
    set(gca,'xtick',[]);
    set(gca,'ytick',[1:1:iter]);
    set(gca,'yticklabel',y_label);

    subplot(1,12,12)
    plotCh = i:4:num_channels;
    if optionIdx == 1
        imagesc(1,1:1:iter,normPower(find(freq == analysisFreq),plotCh)');
    else
        imagesc(1,1:1:iter,power(find(freq == analysisFreq),plotCh)');
    end
    axis xy
    colormap('jet')
    set(gca,'xtick',[]);
    ylim([-1,iter+1])
    clear y_label
    iter = 0;
    plotCh = i:4:num_channels;
    for ch = plotCh
        iter = iter+1;
        y_label{iter} = [num2str(ch) ' ch'];
    end
    if optionIdx == 1
        title(['Normalized power at ' num2str(analysisFreq) ' Hz']);
        caxis([0,2/gain]);
    else
        title(['Power at ' num2str(analysisFreq) ' Hz']);
        caxis([0,10/gain]);
    end
    set(gca,'ytick',[1:1:iter]);
    set(gca,'yticklabel',y_label);
    saveas(gcf,[saveName '_resting_align_' num2str(i) '.fig']);
    saveas(gcf,[saveName '_resting_align_' num2str(i) '.bmp']);
    pause(1)
end

