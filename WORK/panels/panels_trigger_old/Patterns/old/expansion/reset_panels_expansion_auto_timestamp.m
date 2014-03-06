clear
clc
t=0;
t_start=clock;
fname = [num2str(t_start(1)),'-',num2str(t_start(2)),'-',num2str(t_start(3)),'_',num2str(t_start(4)),'-',num2str(t_start(5)),'.mat'];

samplerate=100;
samplespertrigger=100;

AI = analoginput('mcc',0);
hwch = [0];
manes = cell(size(hwch));
names{1} = 'SWtrigger';
chans = addchannel(AI,hwch,'ames');
set(AI,'SampleRate',samplerate);
set(AI,'SamplesPerTrigger',samplespertrigger);
% input('start')

stop = 0;

while stop == 0    
    start(AI)
%     pause(1/samplerate);
    pause(1);
    data = mean(getdata(AI));
    
    if data<-1
        t=t+1;
        trigger_time(t,:)=clock;
        save(fname,'trigger_time');
        
%         pause(1/samplerate);
        pause(1);
        Panel_com('stop')
        Panel_com('set_pattern_id',1);
        Panel_com('send_gain_bias',[0 0 50 0]);
        Panel_com('set_position', [117 1]);
        Panel_com('all_off');
    end
end
        
delete(AI);
clear chans AI
        
        
        
        




Panel_com('stop');
Panel_com('set_pattern_id',1);
Panel_com('send_gain_bias',[0 0 50 0]);
Panel_com('set_position', [117 1]);
% Panel_com('all_off');
% Panel_com('start');