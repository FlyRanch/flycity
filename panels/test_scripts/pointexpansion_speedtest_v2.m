clear
clc
close all

cd('D:\Data\flytracker\paneltestdata');

t=0;
t_start=clock;
min_old=t_start(5);
fname = [num2str(t_start(1)),'-',num2str(t_start(2)),'-',num2str(t_start(3)),'_',num2str(t_start(4)),num2str(t_start(5)),'h.mat'];

% settings
samplerate=20000;
samplespertrigger=20000;

pattern_nr = 1;
config_nr = 1;
expansion_location = 1; % only 1
posfunc_nr = [1 1];

posfunc_freq = 198 % slow speed
% posfunc_freq = 565 % high speed

% gain = 125
% bias = 22
if exist('gain')==1 && exist('bias')==1
    expansion_freq = gain + bias*2.5
else
    expansion_freq = 565
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

% % initiate panel
% Panel_com('set_config_id',config_nr);
% pause(2)
Panel_com('stop');
Panel_com('all_on');


% % normal mode: load pattern to controller
% Panel_com('set_pattern_id',pattern_nr); 

% % SUPERFAST MODE: load pattern to panels
% Panel_com('load_pattern_2panels',pattern_nr); 
% pause(10)

% Panel_com('send_gain_bias',[0 0 expansion_freq 0]);
if exist('gain')==1 && exist('bias')==1
    Panel_com('send_gain_bias',[0 0 gain bias]);
else
    Panel_com('send_gain_bias',[0 0 expansion_freq 0]);
end

% % POSITION FUNCTION SUPERFAST MODE: load pattern to panels
% Panel_com('load_pattern_2panels',pattern_nr); 
% pause(10)
% Panel_com('set_posfunc_id',posfunc_nr); 
% Panel_com('set_funcx_freq',0); 
% Panel_com('set_funcx_freq',posfunc_freq); 

% Panel_com('all_on');
% Panel_com('all_off');
Panel_com('set_position', [expansion_location 1]);

pause(1)

Panel_com('start');
tic;
start(AI)
% pause(samplespertrigger/samplerate); % wait untill sampling is done
data = getdata(AI);
time = toc

Panel_com('stop')
delete(AI);
clear chans AI

% hold on
plot(data)
grid on