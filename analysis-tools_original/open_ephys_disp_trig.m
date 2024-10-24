% for check trigger timing in Neuropixels and Open Ephys
% Programmed by Akito Kosugi
% v.3.0. 07.30.2024

clc

%% Initialization

screenSize = get(0,'ScreenSize');


%% Data loading

recStartTime = d_daq.Timestamps(1);
leverTime = d_daq.Timestamps-recStartTime;
leverPos = d_daq.Data(leverCh,:)'.*bit_volts_daq ;

syncTrigIdx = find(d_daq_2.Data == syncTrigCh);
syncTrigTime = d_daq_2.Timestamps(syncTrigIdx);
syncTrigTime = syncTrigTime-recStartTime;

leverTrigIdx = find(d_daq_2.Data == leverTrigCh);
leverTrigTime = d_daq_2.Timestamps(leverTrigIdx)-recStartTime;
cueTrigIdx = find(d_daq_2.Data == cueTrigCh);
cueTrigTime = d_daq_2.Timestamps(cueTrigIdx)-recStartTime;

messageTime = d_event.Timestamps-recStartTime;


%% Procession

k = 0;
for i = 1:length(cueTrigTime)
    temp =  find((leverTrigTime - cueTrigTime(i)) >0);
    if isempty(temp) == 0
        temp2 = leverTrigTime(min(temp))-cueTrigTime(i);
        if temp2 < taskDuration
            k = k+1;
            leverTrigIdx_after_cue(k) = min(temp);        
            leverTrigTime_after_cue(k) = leverTrigTime(min(temp));
        end 
    end
end
leverTrigIdx_after_cue = leverTrigIdx_after_cue';
leverTrigTime_after_cue = leverTrigTime_after_cue';


%% Plot

f = figure('position',[screenSize(1)+screenSize(3)*1/10 screenSize(2)+screenSize(4)*2/10 screenSize(3)*6/10 screenSize(4)*5/10]);
subplot(5,1,1)
hold on
for i = 1:length(messageTime)
    plot([messageTime(i),messageTime(i)],[0,1],'b','linewidth',1);
end
title('Word input')
xlim([min(leverTime),max(leverTime)])
set(gca,'ytick',[])
set(gca,'fontsize',14);

subplot(5,1,2)
hold on
for i = 1:length(cueTrigTime)
    plot([cueTrigTime(i),cueTrigTime(i)],[0,1],'r','linewidth',1);
end
title('Cue trigger')
xlim([min(leverTime),max(leverTime)])
set(gca,'ytick',[])
set(gca,'fontsize',14);

subplot(5,1,3)
hold on
for i = 1:length(leverTrigTime_after_cue)
    plot([leverTrigTime_after_cue(i),leverTrigTime_after_cue(i)],[0,1],'g','linewidth',1);
end
set(gca,'ytick',[])
xlim([min(leverTime),max(leverTime)])
title('Lever trigger after cue')
set(gca,'fontsize',14);

subplot(5,1,4)
hold on
for i = 1:length(leverTrigTime)
    plot([leverTrigTime(i),leverTrigTime(i)],[0,1],'g','linewidth',1);
end
set(gca,'ytick',[])
xlim([min(leverTime),max(leverTime)])
title('All lever trigger')
set(gca,'fontsize',14);

subplot(5,1,5)
plot(leverTime,leverPos,'k');
xlim([min(leverTime),max(leverTime)])
ylim([-6,6]);
title('Lever position')
xlabel('Time from recording start [s]');
set(gca,'fontsize',14);

saveas(gcf,[saveName '.fig']);
saveas(gcf,[saveName '.bmp']);


%% Save

saveData.recStartTime = recStartTime;
saveData.syncTrigTime = syncTrigTime;
saveData.leverTrigTime = leverTrigTime;
saveData.leverTrigIdx_after_cue = leverTrigIdx_after_cue;
saveData.leverTrigTime_after_cue = leverTrigTime_after_cue;
saveData.cueTrigTime = cueTrigTime;
saveData.messageTime = messageTime;
save([saveName '.mat'],'saveData');

dataLength_max = max([length(messageTime),length(cueTrigTime),length(leverTrigIdx_after_cue),length(leverTrigTime_after_cue),length(leverTrigTime)]);
saveData_csv = nan(dataLength_max,5);
saveData_csv(1:length(messageTime),1) = messageTime;
saveData_csv(1:length(cueTrigTime),2) = cueTrigTime;
saveData_csv(1:length(leverTrigIdx_after_cue),3) = leverTrigIdx_after_cue;
saveData_csv(1:length(leverTrigTime_after_cue),4) = leverTrigTime_after_cue;
saveData_csv(1:length(leverTrigTime),5) = leverTrigTime;
saveMatrix = ["Word input time" "Cue trigger time" "Lever trigger after cue index" "Lever trigger after cue time" "All lever trigger time"; saveData_csv];
writematrix(saveMatrix,[saveName '.csv']);
