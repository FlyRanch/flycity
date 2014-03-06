
ACQUISITIONTIME = 1;
SAMPLERATE = 25000; %MAX
AI = analoginput('mcc',0);
hwch = [0 1];
names = cell(size(hwch));
names{1} = 'SW';
names{2} = 'none';
chans = addchannel(AI,hwch,'names');
set(AI,'SampleRate',SAMPLERATE);
set(AI,'SamplesPerTrigger',round(ACQUISITIONTIME*get(AI,'SampleRate')))
input('start?')
tic
start(AI)
while islogging(AI)
    pause(.1);
end
[data,time,abstime,events] = getdata(AI);
delete(AI);
clear chans AI
toc

outFileName = datestr(clock,'yyyymmdd_HHMMSS')
save(outFileName);