close all

% settings
samplerate=20000;
samplespertrigger=20000;

% initiate daq
AI = analoginput('mcc',0);
hwch = [0];
manes = cell(size(hwch));
names{1} = 'SWtrigger';
chans = addchannel(AI,hwch,'ames');
set(AI,'SampleRate',samplerate);
set(AI,'SamplesPerTrigger',samplespertrigger);

% Panel_com('start');

start(AI)
% pause(samplespertrigger/samplerate); % wait untill sampling is done
data = getdata(AI);
time = toc

% Panel_com('stop')

delete(AI);
clear chans AI

% hold on
plot(data)
grid on