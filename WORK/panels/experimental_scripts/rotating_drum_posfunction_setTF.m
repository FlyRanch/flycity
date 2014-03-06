clear
clc

% cd('D:\Data\flytracker\triggerlog');
cd('E:\Dropbox\Projects\flytracker\triggerlog');
mfilename

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

start_frame = 1 % ONLY mid panel
% start_frame = 20 % MID panel expansion
% start_frame = 1 % TOP panel expansion
% start_frame = 40 % BOTTOM panel expansion

%% ON/OFF
% done = 0;
% while done == 0
%     panels_off_input = input('panels OFF after stimulus? Y/N [Y]: ', 's');
%     done = 1;
%     if isempty(panels_off_input)
%         panels_off_input = 'Y';
%     end
%     
%     if double(panels_off_input) == 89
%         panels_off = 1;
%     elseif double(panels_off_input) == 121
%         panels_off = 1;
%     elseif double(panels_off_input) == 78
%         panels_off = 0;
%     elseif double(panels_off_input) == 110
%         panels_off = 0;
%     else
%         done = 0;
%     end
% end

%% temp freq
done = 0;
while done == 0
    TF = input('temporal frequency?: ');
    done = 1;
    if isempty(TF)
        done = 0;
    end
end

%% panel speed

pattern_cycle = 16; % number of frames for a cycle
desired_freq = TF; % frequency in hertz


% Drum Rotation
frame_freq = desired_freq * pattern_cycle;
samplerate=1000;
samplespertrigger=1000;

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
% Panel_com('send_gain_bias',[0 0 frame_freq 0]);

% enable ext trigger
% Panel_com('enable_extern_trig');
Panel_com('set_position', [start_frame 1]);

% YES posfunc
% Panel_com('set_posfunc_id',posfunc_nr); 
% Panel_com('set_funcx_freq',0); 
Panel_com('set_funcy_freq',frame_freq); 

% test expansion
Panel_com('start');

start(AI)
% pause(samplespertrigger/samplerate); % wait untill sampling is done
data = getdata(AI);
Panel_com('stop')
Panel_com('set_position', [start_frame 1]);
% hold on
plot(data,'*')
grid on

save(fname,'TF');