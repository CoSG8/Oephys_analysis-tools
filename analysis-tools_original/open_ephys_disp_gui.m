% Programmed by Akito Kosugi
% v.1.4 07.02.2023

close all
clc

%% Initialization

screenSize = get(0,'screensize');

uiIndex = 1;
param_np_ui1.PCIdx = 1;
param_np_ui1.probeIdx = 1;
param_np_ui1.columnIdx = 1;
param_np_ui3.analysisTime = 20;
param_np_ui3.gain = 1;
param_np_ui3.optionIdx = 1;
param_np_ui4.trigCh = 4;
param_np_ui4.cdpCh = 6;
param_np_ui4.spikeDetectionCh = 101;
param_np_ui4.num_maxTrig = 20;
param_np_ui4.gain = 5;
param_np_ui5.trigCh = 7;
param_np_ui5.leverCh = 8;
param_np_ui5.spikeDetectionCh = 101;
param_np_ui5.num_maxTrig = 30;


%% Set GUI 1

%---GUI 1: Select analysis step---%
if exist('param_np_ui1.mat','file') > 0
    load('param_np_ui1.mat')
end

h = figure(100);
set(h,'position',[screenSize(1)+screenSize(3)*2/10 screenSize(2)+screenSize(4)*2/10 screenSize(3)*6/10 screenSize(4)*6/10]);
set(gcf,'name','Neuropixels brief analysis');

ui1 =uibuttongroup(h,'title','Setting',...
            'unit','normalize',...
            'position',[0.05 0.3 0.5 0.62],...
            'fontsize',20,...
            'fontweight','bold');    
        uicontrol(ui1,'style','text',...
            'unit','normalize',...
            'position',[0.05 0.7 0.6 0.15],...
            'string','PC:',...
            'fontsize',15,...
            'fontweight','bold',...
            'HorizontalAlignment','left');        
        ui1_1 = uicontrol(ui1,'style','popupmenu',...
            'unit','normalize',...
            'position',[0.4 0.80 0.15 0.05],...
            'string',{'Windows','Mac','Linux'},...
            'value',param_np_ui1.PCIdx);
        uicontrol(ui1,'style','text',...
            'unit','normalize',...
            'position',[0.05 0.55 0.6 0.15],...
            'string','Probe:',...
            'fontsize',15,...
            'fontweight','bold',...
            'HorizontalAlignment','left');            
        ui1_2 = uicontrol(ui1,'style','popupmenu',...
            'unit','normalize',...
            'position',[0.4 0.65 0.15 0.05],...
            'string',{'A','B','C'},...
            'value',1);
        uicontrol(ui1,'style','text',...
            'unit','normalize',...
            'position',[0.05 0.4 0.6 0.15],...
            'string','Number of Columns:',...
            'fontsize',15,...
            'fontweight','bold',...
            'HorizontalAlignment','left');              
        ui1_3 = uicontrol(ui1,'style','popupmenu',...
            'unit','normalize',...
            'position',[0.4 0.5 0.15 0.05],...
            'string',{'4','2'},...
            'value',1);
                
ui2 =uibuttongroup(h,'title','Purpose',...
            'unit','normalize',...
            'position',[0.35 0.3 0.3 0.62],...
            'fontsize',20,...
            'fontweight','bold');    

ui2_1 = uicontrol(ui2,'style','radiobutton',...
            'unit','normalize',...
            'position',[0.1 0.8 0.8 0.1],...
            'string','Noise check',...
            'fontsize',15,...
            'fontweight','bold',...
            'callback', 'uiIndex = 1;');
        
ui2_2=uicontrol(ui2,'style','radiobutton',...
            'unit','normalize',...
            'position',[0.1 0.65 0.8 0.1],...
            'string','LFP mapping',...
            'fontsize',15,...
            'fontweight','bold',...
            'callback', 'uiIndex = 2;');        

ui2_2=uicontrol(ui2,'style','radiobutton',...
            'unit','normalize',...
            'position',[0.1 0.5 0.8 0.1],...
            'string','Lever task',...
            'fontsize',15,...
            'fontweight','bold',...
            'callback', 'uiIndex = 3;');        

uicontrol(h,'style','push',...
            'unit','normalize',...
            'position',[0.78 0.11 0.2 0.1],...
            'string','OK',...
            'fontsize',25,...
            'fontweight','bold',...
            'callback','PCIdx=get(ui1_1,''value'');probeIdx=get(ui1_2,''value'');columnIdx=get(ui1_3,''value'');close(100);');    
                
uiwait

%% Set GUI 2

%---GUI 2: Input experiment parameters---%

switch uiIndex
    case 1
        if exist('param_np_ui3.mat','file') > 0
            load('param_np_ui3.mat')
        end
    case 2
        if exist('param_np_ui4.mat','file') > 0
            load('param_np_ui4.mat')
        end
    case 3
        if exist('param_np_ui5.mat','file') > 0
            load('param_np_ui5.mat')
        end        
end
switch uiIndex
    case 1
        h = figure(200);
        set(h,'position',[screenSize(1)+screenSize(3)*2/10 screenSize(2)+screenSize(4)*2/10 screenSize(3)*1/5 screenSize(4)*6/10]);
        set(gcf,'name','Input parameters')
        
        ui3=uibuttongroup(h,'title','Analysis parameters',...
            'unit','normalize',...
            'position',[0.05 0.3 0.9 0.62],...
            'fontsize',20,...
            'fontweight','bold');
        
        uicontrol(ui3,'style','text',...
            'unit','normalize',...
            'position',[0.05 0.81 0.6 0.15],...
            'string','Analysis time:',...
            'fontsize',15,...
            'fontweight','bold');
        ui3_1_1 = uicontrol(ui3,'style','text',...
            'unit','normalize',...
            'position',[0.9 0.905 0.2 0.05],...
            'string','[s]',...
            'fontsize',12,...            
            'HorizontalAlignment','left');
        ui3_1_2 = uicontrol(ui3,'style','edit',...
            'unit','normalize',...
            'position',[0.7 0.90 0.2 0.05],...
            'string',num2str(param_np_ui1.analysisTime),...
            'HorizontalAlignment','left');
        uicontrol(ui3,'style','text',...
            'unit','normalize',...
            'position',[0.05 0.65 0.6 0.15],...
            'string','Gain:',...
            'fontsize',15,...
            'fontweight','bold');
        ui3_2_1 = uicontrol(ui3,'style','text',...
            'unit','normalize',...
            'position',[0.67 0.745 0.2 0.05],...
            'string','x',...
            'fontsize',12,...            
            'HorizontalAlignment','left');
        ui3_2_2 = uicontrol(ui3,'style','edit',...
            'unit','normalize',...
            'position',[0.7 0.74 0.2 0.05],...
            'string',num2str(param_np_ui1.gain),...
            'HorizontalAlignment','left');
        uicontrol(ui3,'style','text',...
            'unit','normalize',...
            'position',[0.05 0.5 0.6 0.15],...
            'string','Display option:',...
            'fontsize',15,...
            'fontweight','bold');
        ui3_3_1 = uicontrol(ui3,'style','popupmenu',...
            'unit','normalize',...
            'position',[0.7 0.59 0.23 0.05],...
            'string',{'Normalize','Absolute'},...
            'value',param_np_ui1.optionIdx,...            
            'HorizontalAlignment','left');
        
        uicontrol(h,'style','push',...
            'unit','normalize',...
            'position',[0.78 0.11 0.2 0.1],...
            'string','OK',...
            'fontsize',25,...
            'fontweight','bold',...
            'callback','analysisTime=str2num(get(ui3_1_2,''string''));gain=str2num(get(ui3_2_2,''string''));optionIdx=get(ui3_3_1,''value'');close(200);');

        uiwait

    case 2
        h = figure(200);
        set(h,'position',[screenSize(1)+screenSize(3)*2/10 screenSize(2)+screenSize(4)*2/10 screenSize(3)*1/5 screenSize(4)*6/10]);
        set(gcf,'name','Input parameters')
        
        ui4 =uibuttongroup(h,'title','Analysis parameters',...
            'unit','normalize',...
            'position',[0.05 0.3 0.9 0.62],...
            'fontsize',20,...
            'fontweight','bold');        
        uicontrol(ui4,'style','text',...
            'unit','normalize',...
            'position',[0.05 0.81 0.6 0.15],...
            'string','Trigger channel:',...
            'fontsize',15,...
            'fontweight','bold');
        ui4_1_1 = uicontrol(ui4,'style','text',...
            'unit','normalize',...
            'position',[0.65 0.90 0.2 0.05],...
            'string','AI',...
            'fontsize',12,...            
            'HorizontalAlignment','left');
        ui4_1_2 = uicontrol(ui4,'style','edit',...
            'unit','normalize',...
            'position',[0.7 0.90 0.2 0.05],...
            'string',num2str(param_np_ui4.trigCh-1),...
            'HorizontalAlignment','left');        
        uicontrol(ui4,'style','text',...
            'unit','normalize',...
            'position',[0.05 0.65 0.6 0.15],...
            'string','CDP channel:',...
            'fontsize',15,...
            'fontweight','bold');
        ui4_2_1 = uicontrol(ui4,'style','text',...
            'unit','normalize',...
            'position',[0.65 0.74 0.2 0.05],...
            'string','AI',...
            'fontsize',12,...            
            'HorizontalAlignment','left');
        ui4_2_2 = uicontrol(ui4,'style','edit',...
            'unit','normalize',...
            'position',[0.7 0.74 0.2 0.05],...
            'string',num2str(param_np_ui4.cdpCh-1),...
            'HorizontalAlignment','left');
        uicontrol(ui4,'style','text',...
            'unit','normalize',...
            'position',[0.05 0.5 0.6 0.15],...
            'string','Spike detection chanel:',...
            'fontsize',15,...
            'fontweight','bold');
        ui4_3_1 = uicontrol(ui4,'style','text',...
            'unit','normalize',...
            'position',[0.9 0.59 0.2 0.05],...
            'string','ch',...
            'fontsize',12,...            
            'HorizontalAlignment','left');
        ui4_3_2 = uicontrol(ui4,'style','edit',...
            'unit','normalize',...
            'position',[0.7 0.59 0.2 0.05],...
            'string',num2str(param_np_ui4.spikeDetectionCh),...
            'HorizontalAlignment','left');        
        uicontrol(ui4,'style','text',...
            'unit','normalize',...
            'position',[0.05 0.35 0.6 0.15],...
            'string','Number of trigger:',...
            'fontsize',15,...
            'fontweight','bold');
        ui4_4_1 = uicontrol(ui4,'style','edit',...
            'unit','normalize',...
            'position',[0.7 0.44 0.2 0.05],...
            'string',num2str(param_np_ui4.num_maxTrig),...
            'HorizontalAlignment','left');
        uicontrol(ui4,'style','text',...
            'unit','normalize',...
            'position',[0.05 0.2 0.6 0.15],...
            'string','Gain:',...
            'fontsize',15,...
            'fontweight','bold');
        ui4_5_1 = uicontrol(ui4,'style','text',...
            'unit','normalize',...
            'position',[0.67 0.295 0.2 0.05],...
            'string','x',...
            'fontsize',12,...            
            'HorizontalAlignment','left');
        ui4_5_2 = uicontrol(ui4,'style','edit',...
            'unit','normalize',...
            'position',[0.7 0.29 0.2 0.05],...
            'string',num2str(param_np_ui4.gain),...
            'HorizontalAlignment','left');        
        uicontrol(h,'style','push',...
            'unit','normalize',...
            'position',[0.78 0.11 0.2 0.1],...
            'string','OK',...
            'fontsize',25,...
            'fontweight','bold',...
            'callback','trigCh=str2num(get(ui4_1_2,''string''))+1;cdpCh=str2num(get(ui4_2_2,''string''))+1;spikeDetectionCh=str2num(get(ui4_3_2,''string''));num_maxTrig=str2num(get(ui4_4_1,''string''));gain=str2num(get(ui4_5_2,''string''));close(200);');
        uiwait

    case 3
        h = figure(200);
        set(h,'position',[screenSize(1)+screenSize(3)*2/10 screenSize(2)+screenSize(4)*2/10 screenSize(3)*1/5 screenSize(4)*6/10]);
        set(gcf,'name','Input parameters')
        
        ui5 =uibuttongroup(h,'title','Analysis parameters',...
            'unit','normalize',...
            'position',[0.05 0.3 0.9 0.62],...
            'fontsize',20,...
            'fontweight','bold');
        uicontrol(ui5,'style','text',...
            'unit','normalize',...
            'position',[0.05 0.81 0.6 0.15],...
            'string','Trigger channel:',...
            'fontsize',15,...
            'fontweight','bold');
        ui5_1_1 = uicontrol(ui5,'style','text',...
            'unit','normalize',...
            'position',[0.65 0.90 0.2 0.05],...
            'string','DI',...
            'fontsize',12,...            
            'HorizontalAlignment','left');
        ui5_1_2 = uicontrol(ui5,'style','edit',...
            'unit','normalize',...
            'position',[0.7 0.90 0.2 0.05],...
            'string',num2str(param_np_ui5.trigCh-1),...
            'HorizontalAlignment','left');
        uicontrol(ui5,'style','text',...
            'unit','normalize',...
            'position',[0.05 0.65 0.6 0.15],...
            'string','Lever channel:',...
            'fontsize',15,...
            'fontweight','bold');
        ui5_2_1 = uicontrol(ui5,'style','text',...
            'unit','normalize',...
            'position',[0.65 0.74 0.2 0.05],...
            'string','AI',...
            'fontsize',12,...            
            'HorizontalAlignment','left');
        ui5_2_2 = uicontrol(ui5,'style','edit',...
            'unit','normalize',...
            'position',[0.7 0.74 0.2 0.05],...
            'string',num2str(param_np_ui5.leverCh-1),...
            'HorizontalAlignment','left');        
        uicontrol(ui5,'style','text',...
            'unit','normalize',...
            'position',[0.05 0.5 0.6 0.15],...
            'string','Spike detection channel:',...
            'fontsize',15,...
            'fontweight','bold');
        ui5_3_1 = uicontrol(ui5,'style','text',...
            'unit','normalize',...
            'position',[0.9 0.59 0.2 0.05],...
            'string','Ch',...
            'fontsize',12,...            
            'HorizontalAlignment','left');
        ui5_3_2 = uicontrol(ui5,'style','edit',...
            'unit','normalize',...
            'position',[0.7 0.59 0.2 0.05],...
            'string',num2str(param_np_ui5.spikeDetectionCh),...
            'HorizontalAlignment','left');
        uicontrol(ui5,'style','text',...
            'unit','normalize',...
            'position',[0.05 0.35 0.6 0.15],...
            'string','Number of trigger:',...
            'fontsize',15,...
            'fontweight','bold');
        ui5_4_1 = uicontrol(ui5,'style','edit',...
            'unit','normalize',...
            'position',[0.7 0.44 0.2 0.05],...
            'string',num2str(param_np_ui5.num_maxTrig),...
            'HorizontalAlignment','left');        
        uicontrol(h,'style','push',...
            'unit','normalize',...
            'position',[0.78 0.11 0.2 0.1],...
            'string','OK',...
            'fontsize',25,...
            'fontweight','bold',...
            'callback','trigCh=str2num(get(ui5_1_2,''string''))+1;leverCh=str2num(get(ui5_2_2,''string''))+1;spikeDetectionCh=str2num(get(ui5_3_2,''string''));num_maxTrig=str2num(get(ui5_4_1,''string''));close(200);');        
        uiwait

end


%% Save

param_np_ui1.PCIdx = PCIdx;
param_np_ui1.probeIdx = probeIdx;
param_np_ui1.columnIdx = columnIdx;
save('param_np_ui1.mat','param_np_ui1');

switch uiIndex
    case 1
        param_np_ui3.analysisTime = analysisTime;
        param_np_ui3.gain = gain;
        param_np_ui3.optionIdx = optionIdx;
        save('param_np_ui3.mat','param_np_ui3');
    case 2
        param_np_ui4.trigCh = trigCh;
        param_np_ui4.cdpCh = cdpCh;
        param_np_ui4.spikeDetectionCh = spikeDetectionCh;
        param_np_ui4.num_maxTrig = num_maxTrig;        
        param_np_ui4.gain = gain;
        save('param_np_ui4.mat','param_np_ui4');
    case 3
        param_np_ui5.trigCh = trigCh;
        param_np_ui5.leverCh = leverCh;
        param_np_ui5.spikeDetectionCh = spikeDetectionCh;
        param_np_ui5.num_maxTrig = num_maxTrig;        
        save('param_np_ui5.mat','param_np_ui5');        
end