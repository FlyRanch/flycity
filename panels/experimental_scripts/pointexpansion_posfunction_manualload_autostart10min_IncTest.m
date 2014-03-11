clear
clc

% cd('D:\Data\flytracker\triggerlog');
cd('E:\Dropbox\Projects\flytracker\triggerlog');

trigger_time = [];
ERROR_time = [];

t=0;
t_error=0;
t_start=clock;

min_old=t_start(5);
fname = [num2str(t_start(1)),'-',num2str(t_start(2)),'-',num2str(t_start(3)),'_',num2str(t_start(4)),num2str(t_start(5)),'h.mat'];

%% settings
% pattern_nr = 1; % POSFUNC point expansion shifted to 0 location (in cali)
% posfunc_nr = 1; % 0 to 165deg, then 300 frames steady

expansion_location = 1 % ONLY mid panel
% expansion_location = 20 % MID panel expansion
% expansion_location = 1 % TOP panel expansion
% expansion_location = 40 % BOTTOM panel expansion

%% panel speed
% % slow exp
% expansion_freq = 198 
% samplerate=1000;
% samplespertrigger=2000;

% fast exp
expansion_freq = 565 
samplerate=1000;
samplespertrigger=1000;

% % medium speed exp
% expansion_freq = 340 
% samplerate=1000;
% samplespertrigger=1000;

%% initiate daq
AI = analoginput('mcc',0);
hwch = [0];
manes = cell(size(hwch));
names{1} = 'SWtrigger';
chans = addchannel(AI,hwch,'ames');
set(AI,'SampleRate',samplerate);
set(AI,'SamplesPerTrigger',samplespertrigger);
% input('start')

% initiate panel
Panel_com('stop');

% % SUPERFAST MODE: load pattern to panels
% Panel_com('load_pattern_2panels',pattern_nr); 
% pause(10)
% 
% % CONVENTIONAL MODE
% Panel_com('set_pattern_id',pattern_nr); 
% pause(1)


% % NO posfunc
% Panel_com('send_gain_bias',[0 0 expansion_freq 0]);

% enable ext trigger
Panel_com('enable_extern_trig');

Panel_com('set_position', [expansion_location 1]);

% YES posfunc
% Panel_com('set_posfunc_id',posfunc_nr); 
% Panel_com('set_funcx_freq',0); 
Panel_com('set_funcy_freq',expansion_freq); 

% test expansion
Panel_com('start');
start(AI)
% pause(samplespertrigger/samplerate); % wait untill sampling is done
data = getdata(AI);
Panel_com('stop')
Panel_com('set_position', [expansion_location 1]);
% hold on
plot(data,'*')
grid on

%% loop
stop = 0;
while stop == 0    
    time_now = clock;
    min_now = time_now(5);
    
    start(AI)
    pause(0.05);
    
%     data1 = (getdata(AI));
%     data = max(data1);

% % panels OFF
%     data = max(getdata(AI));
%     if data>1
        
% panels ON
    data = min(getdata(AI));    
    if data<0.8
        
        Panel_com('stop');
        Panel_com('set_position', [expansion_location 1]);
        
        t=t+1;
%         clc
        trigger_time(t,:)=time_now
%         save(fname,'trigger_time','ERROR_time','pattern_nr','posfunc_nr','expansion_freq','expansion_location');
        save(fname,'trigger_time','ERROR_time','expansion_freq','expansion_location');
        
        pause(samplespertrigger/samplerate);

    elseif floor(min_now/10) > floor(min_old/10) % every 10 min
%     elseif floor(min_now/5) > floor(min_old/5) % every 5 min
        autostart=clock
        
        Panel_com('start')
        
        start(AI)
        pause(.05);
        data = min(getdata(AI));
        pause(.05);
        
        Panel_com('stop');
        Panel_com('set_position', [expansion_location 1]);
        
%         if data>-2
%         if data<2
        if data>1
            
            t_error=t_error+1;
%             clc
            
            ERROR = 'NO RESPONSE'
            ERROR_time(t_error,:)=autostart
            
%             save(fname,'trigger_time','ERROR_time','pattern_nr','posfunc_nr','expansion_freq','expansion_location');
            save(fname,'trigger_time','ERROR_time','expansion_freq','expansion_location');
        end
    end
            min_old = min_now;
end
        
delete(AI);
clear chans AI

Panel_com('stop');
Panel_com('set_position', [expansion_location 1]);

