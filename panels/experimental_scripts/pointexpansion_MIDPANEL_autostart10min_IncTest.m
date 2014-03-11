clear
clc

% cd('D:\Data\flytracker\triggerlog');
cd('C:\Users\florian\Dropbox\Projects\flytracker\triggerlog');

trigger_time = [];
ERROR_time = [];

t=0;
t_error=0;
t_start=clock;

min_old=t_start(5);
fname = [num2str(t_start(1)),'-',num2str(t_start(2)),'-',num2str(t_start(3)),'_',num2str(t_start(4)),num2str(t_start(5)),'h.mat'];

% settings
pattern_nr = 2; % point expansion shifted to 0 location (in cali)
expansion_freq = 48 % max freq without frame dropping (ON: Bus2 & 2 panels from bus1 & 3 panels from bus3)
expansion_location = 20 % MID panel expansion
% expansion_location = 1 % TOP panel expansion
% expansion_location = 40 % BOTTOM panel expansion

samplerate=100;
samplespertrigger=100;

% initiate daq
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
Panel_com('all_on');
% Panel_com('all_off');
Panel_com('set_pattern_id',pattern_nr); 
Panel_com('send_gain_bias',[0 0 expansion_freq 0]);
Panel_com('set_position', [expansion_location 1]);
% Panel_com('all_off');

stop = 0;

while stop == 0    
    time_now = clock;
    min_now = time_now(5);
    
    start(AI)
%     pause(1/samplerate);
    pause(0.05);
    data = mean(getdata(AI));
    
%     if data<-2
    if data>2
        t=t+1;
%         clc
        trigger_time(t,:)=time_now
        save(fname,'trigger_time','ERROR_time','pattern_nr','expansion_freq','expansion_location');
        
%         pause(1/samplerate);
        pause(0.05);
        
        Panel_com('stop');
Panel_com('all_on');
% Panel_com('all_off');
        Panel_com('set_pattern_id',pattern_nr); 
        Panel_com('send_gain_bias',[0 0 expansion_freq 0]);
        Panel_com('set_position', [expansion_location 1]);
        % Panel_com('all_off');

    elseif floor(min_now/10) > floor(min_old/10) % every 10 min
%     elseif floor(min_now/5) > floor(min_old/5) % every 5 min
        autostart=clock
        Panel_com('start')
        
        start(AI)
        pause(0.05);
        data = getdata(AI);
        data = mean(data);
%         data = mean(getdata(AI));
        pause(.05);
        
%         if data>-2
        if data<2
            
            t_error=t_error+1;
%             clc
            
            ERROR = 'NO RESPONSE'
            ERROR_time(t_error,:)=autostart
            
            save(fname,'trigger_time','ERROR_time','pattern_nr','expansion_freq','expansion_location');

        end
        
        Panel_com('stop');
Panel_com('all_on');
% Panel_com('all_off');
        Panel_com('set_pattern_id',pattern_nr); 
        Panel_com('send_gain_bias',[0 0 expansion_freq 0]);
        Panel_com('set_position', [expansion_location 1]);
        % Panel_com('all_off');
    end
            min_old = min_now;
end
        
delete(AI);
clear chans AI

Panel_com('stop');
Panel_com('all_on');
% Panel_com('all_off');
Panel_com('set_pattern_id',pattern_nr); 
Panel_com('send_gain_bias',[0 0 expansion_freq 0]);
Panel_com('set_position', [expansion_location 1]);
% Panel_com('all_off');
