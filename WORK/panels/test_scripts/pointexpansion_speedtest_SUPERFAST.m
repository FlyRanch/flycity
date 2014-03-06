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

% pattern_nr = 6; % mirror expansion
pattern_nr = 5; % 4 panel expansion
expansion_location = 1; % only 1

% pattern_nr = 4; % point expansion DOUBLESPEED
% expansion_location = 20; % mid panel

% gain = 127
% bias = 0
if exist('gain')==1 && exist('bias')==1
    expansion_freq = gain + bias*2.5
else
    expansion_freq = 48
end

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

% % use conventional mode
% Panel_com('set_pattern_id',pattern_nr); 

% SUPERFAST MODE: load pattern to panels
Panel_com('load_pattern_2panels',pattern_nr); 



% Panel_com('send_gain_bias',[0 0 expansion_freq 0]);
if exist('gain')==1 && exist('bias')==1
    Panel_com('send_gain_bias',[0 0 gain bias]);
else
    Panel_com('send_gain_bias',[0 0 expansion_freq 0]);
end

% Panel_com('all_on');
% Panel_com('all_off');
Panel_com('set_position', [expansion_location 1]);


tic;
Panel_com('start');
start(AI)
pause(samplerate/samplespertrigger); % wait untill sampling is done
time = toc

data = getdata(AI);

Panel_com('stop')
delete(AI);
clear chans AI

% hold on
plot(data)
% grid on