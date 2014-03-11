clear
% clc
SAMPLERATE = 25000 %MAX
SAMPLESPERTRIGGER = 25000
AI = analoginput('mcc',0);
hwch = [0];
names = cell(size(hwch));
names{1} = 'SWtrigger';
chans = addchannel(AI,hwch,'names');
set(AI,'SampleRate',SAMPLERATE);
set(AI,'SamplesPerTrigger',SAMPLESPERTRIGGER)
input('start?')
tic

for i=1:1
    start(AI)
%     while islogging(AI)
%         pause(1/SAMPLERATE);
%     end
    data = getdata(AI);
%     time(i,1)=toc;
end
toc

% i=1;
% data(i)=-5;
% while data(i)<-2.5
%     i=i+1;
%     start(AI)
% %     while islogging(AI)
% %         pause(1/SAMPLERATE);
% %     end
%     data(i,1) = getdata(AI);
% %     time(i,1)=toc;
% end
    



delete(AI);
clear chans AI
% toc
% 
% outFileName = datestr(clock,'yyyymmdd_HHMMSS')
% save(outFileName);