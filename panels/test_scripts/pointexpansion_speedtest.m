clear
clc

cd('D:\Data\flytracker\paneltestdata');

t=0;
t_start=clock;
min_old=t_start(5);
fname = [num2str(t_start(1)),'-',num2str(t_start(2)),'-',num2str(t_start(3)),'_',num2str(t_start(4)),num2str(t_start(5)),'h.mat'];

% settings
samplerate=1000;
samplespertrigger=1000;

pattern_nr = 6; % 4 panel expansion
% pattern_nr = 7; % mirror expansion
% pattern_nr = 5; % point expansion DOUBLESPEED

expansion_location = 1; % mid panel
expansion_freq = 10
gain = 10
bias = 5
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
Panel_com('set_pattern_id',pattern_nr); 
% Panel_com('send_gain_bias',[0 0 expansion_freq 0]);
Panel_com('send_gain_bias',[0 0 gain bias]);
Panel_com('set_position', [expansion_location 1]);
% Panel_com('all_on');
% Panel_com('all_off');



tic;
Panel_com('start');
start(AI)

time = toc;
pause(1);
data = getdata(AI);

Panel_com('stop')
delete(AI);
clear chans AI

% hold on
plot(data)
% grid on