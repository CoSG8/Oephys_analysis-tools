% Neuropixel analysis for noiose check
% Programmed by Akito Kosugi
% v.1.2 07.01.2023

close all


%% Initialization

time_window = [0,analysisTime]; %[s]
resol_fft = 1;
analysisFreq = 50;


%% Epoching

%--- devided by time---%
plotIdx_rest = time_window(1)*fs_AP+1:time_window(2)*fs_AP;
filtData_AP_epoch_time = filtData_AP(:,plotIdx_rest);


%% Resting analysis

%---FFT---%
noverlap = 0;
nfft = fs_AP/resol_fft;
w = hanning(nfft);
[power,freq] = pwelch(detrend(filtData_AP_epoch_time)',w,noverlap,nfft,fs_AP);
normPower = power./sum(power,1).*100;

rms = mean(filtData_AP_epoch_time.^2,2).^0.5;


%% Plot

%---Resting---%
time_window = [analysisTime-2,analysisTime-1]; %[s]
plotTime = time_window(1):1/fs_AP:time_window(2)-1/fs_AP;
plotIdx = time_window(1)*fs_AP:time_window(2)*fs_AP-1;
switch columnIdx
    case 1
        for i = 1:4
            figure('position',screenSize)
            subplot(1,3,1)
            hold on
            iter = 0;
            for ch = i:4:384
                iter = iter+1;
                plot(plotTime,filtData_AP_epoch_time(ch,plotIdx).*gain*2+1000*iter,'k','linewidth',1);
                y_label{iter} = [num2str(ch) ' ch'];
            end
            set(gca,'ytick',[1000:1000:1000*iter]);    
            set(gca,'yticklabel',y_label);
            ylim([-1000,1000*(iter+1)])
            xlim([time_window(1),time_window(2)])
            title('Raw  data');
            xlabel('Time [ms]');

            subplot(1,3,2)
            hold on
            iter = 0;
            for ch = i:4:384
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
            plotCh = i:4:384;
            imagesc(1,1:1:iter,rms(plotCh,:));
            axis xy
            colormap('jet')
            caxis([50,300/gain]);
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
            plotCh = i:4:384;
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
            plotCh = i:4:384;
            for ch = plotCh
                iter = iter+1;
                y_label{iter} = [num2str(ch) ' ch'];
            end
            if optionIdx == 1
                title(['Normalized power at ' num2str(analysisFreq) ' Hz']);
                caxis([0,2/gain]);
            else
                title(['Power at ' num2str(analysisFreq) ' Hz']);
                caxis([0,60/gain]);
            end
            set(gca,'ytick',[1:1:iter]);
            set(gca,'yticklabel',y_label);
            saveas(gcf,[saveName '_resting_align_' num2str(i) '.fig']);
            saveas(gcf,[saveName '_resting_align_' num2str(i) '.bmp']);
            pause(1)
            close(gcf)
        end
    case 2
        for i = 1:2
            figure('position',screenSize)
            subplot(1,3,1)
            hold on
            iter = 0;
            switch i
                case 1
                    for ch = 1:4:381
                        iter = iter+1;
                        plot(plotTime,filtData_AP_epoch_time(ch,plotIdx).*gain+500*iter,'k','linewidth',1);
                        y_label{iter} = [num2str(ch) ' ch'];
                    end
                    for ch = 3:4:383
                        iter = iter+1;
                        plot(plotTime,filtData_AP_epoch_time(ch,plotIdx).*gain+500*iter,'k','linewidth',1);
                        y_label{iter} = [num2str(ch) ' ch'];
                    end
                case 2
                    for ch = 4:4:384
                        iter = iter+1;
                        plot(plotTime,filtData_AP_epoch_time(ch,plotIdx).*gain+500*iter,'k','linewidth',1);
                        y_label{iter} = [num2str(ch) ' ch'];
                    end
                    for ch = 2:4:382
                        iter = iter+1;
                        plot(plotTime,filtData_AP_epoch_time(ch,plotIdx).*gain+500*iter,'k','linewidth',1);
                        y_label{iter} = [num2str(ch) ' ch'];
                    end
            end
            set(gca,'ytick',[500:500:500*iter]);    
            set(gca,'yticklabel',y_label);
            ylim([-1000,500*(iter+1)])
            xlim([time_window(1),time_window(2)])
            title('Raw  data');
            xlabel('Time [ms]');
            clear y_label
            
            subplot(1,3,2)
            hold on
            iter = 0;
            switch i
                case 1
                    for ch = 1:4:381            
                        iter = iter+1;
                        if optionIdx == 1
                            plot(freq,normPower(:,ch).*gain+1*iter,'k','linewidth',1);
                        else
                            plot(freq,power(:,ch).*gain+20*iter,'k','linewidth',1);
                        end
                        y_label{iter} = [num2str(ch) ' ch'];
                    end
                    for ch = 3:4:383            
                        iter = iter+1;
                        if optionIdx == 1
                            plot(freq,normPower(:,ch).*gain+1*iter,'k','linewidth',1);
                        else
                            plot(freq,power(:,ch).*gain+20*iter,'k','linewidth',1);
                        end
                        y_label{iter} = [num2str(ch) ' ch'];
                    end
                case 2
                    for ch = 4:4:384            
                        iter = iter+1;
                        if optionIdx == 1
                            plot(freq,normPower(:,ch).*gain+1*iter,'k','linewidth',1);
                        else
                            plot(freq,power(:,ch).*gain+20*iter,'k','linewidth',1);
                        end
                        y_label{iter} = [num2str(ch) ' ch'];
                    end
                    for ch = 2:4:382            
                        iter = iter+1;
                        if optionIdx == 1
                            plot(freq,normPower(:,ch).*gain+1*iter,'k','linewidth',1);
                        else
                            plot(freq,power(:,ch).*gain+20*iter,'k','linewidth',1);
                        end
                        y_label{iter} = [num2str(ch) ' ch'];
                    end                    
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
            clear y_label
            
            subplot(1,12,10)
            hold on
            iter = 0;
            switch i
                case 1
                    for ch = 1:4:381
                        iter = iter+1;
                        imagesc(1,iter,rms(ch,:));
                        y_label{iter} = [num2str(ch) ' ch'];
                    end
                    for ch = 3:4:383
                        iter = iter+1;
                        imagesc(1,iter,rms(ch,:));
                        y_label{iter} = [num2str(ch) ' ch'];
                    end
                case 2
                    for ch = 4:4:384
                        iter = iter+1;
                        imagesc(1,iter,rms(ch,:));
                        y_label{iter} = [num2str(ch) ' ch'];
                    end
                    for ch = 2:4:382
                        iter = iter+1;
                        imagesc(1,iter,rms(ch,:));
                        y_label{iter} = [num2str(ch) ' ch'];
                    end                                        
            end
            axis xy
            colormap('jet')
            caxis([50,300/gain]);
            ylim([-1,iter+1])
            title('Root mean square');
            set(gca,'xtick',[]);
            set(gca,'ytick',[1:1:iter]);    
            set(gca,'yticklabel',y_label);
            clear y_label
            
            subplot(1,12,12)
            hold on
            iter = 0;
            switch i
                case 1
                    for ch = 1:4:381
                        iter = iter+1;
                        if optionIdx == 1
                            imagesc(1,1:1:iter,normPower(find(freq == analysisFreq),plotCh)');
                        else
                            imagesc(1,1:1:iter,power(find(freq == analysisFreq),plotCh)');
                        end
                        y_label{iter} = [num2str(ch) ' ch'];
                    end
                    for ch = 3:4:383
                        iter = iter+1;
                        if optionIdx == 1
                            imagesc(1,1:1:iter,normPower(find(freq == analysisFreq),plotCh)');
                        else
                            imagesc(1,1:1:iter,power(find(freq == analysisFreq),plotCh)');
                        end
                        y_label{iter} = [num2str(ch) ' ch'];
                    end
                case 2
                    for ch = 4:4:384
                        iter = iter+1;                        
                        if optionIdx == 1
                            imagesc(1,1:1:iter,normPower(find(freq == analysisFreq),plotCh)');
                        else
                            imagesc(1,1:1:iter,power(find(freq == analysisFreq),plotCh)');
                        end
                        y_label{iter} = [num2str(ch) ' ch'];
                    end
                    for ch = 2:4:382
                        iter = iter+1;                        
                        if optionIdx == 1
                            imagesc(1,1:1:iter,normPower(find(freq == analysisFreq),plotCh)');
                        else
                            imagesc(1,1:1:iter,power(find(freq == analysisFreq),plotCh)');
                        end
                        y_label{iter} = [num2str(ch) ' ch'];
                    end
            end
            axis xy
            colormap('jet')
            set(gca,'xtick',[]);
            ylim([-1,iter+1])
            if optionIdx == 1
                title(['Normalized power at ' num2str(analysisFreq) ' Hz']);
                caxis([0,2/gain]);
            else
                title(['Power at ' num2str(analysisFreq) ' Hz']);
                caxis([0,60/gain]);
            end
            set(gca,'ytick',[1:1:iter]);
            set(gca,'yticklabel',y_label);
            clear y_label            
            saveas(gcf,[saveName '_resting_align_' num2str(i) '.fig']);
            saveas(gcf,[saveName '_resting_align_' num2str(i) '.bmp']);
            pause(1)
            close(gcf)
        end        
end
