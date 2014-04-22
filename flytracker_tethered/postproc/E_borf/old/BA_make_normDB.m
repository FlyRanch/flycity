% make normDB from WBmod files

clear
clc
close all

name = dir('WBdataset_all_steadyNmods*')
% name = dir('WBmod_accelbased_PitchAccel*')
name = name.name;

% list = dir('WBmod*')
% for i=1:length(list)
%     if list(i).isdir == 0
%         name = list(i).name;
%     end
% end

load(name)

% % mod locations
% Fenhanses = [0 .5 1]'
% rollaccels = [0 2.5 5]'
% % pitchaccels = [0 1.5 3]'
% pitchaccels = [0 3]'
% 
% % save norm data
% save('norm_data.mat','rollaccel_norm','pitchaccel_norm','yawaccel_norm','Fenhance_norm',...
%     'f_wb_steady','Fenhanses','rollaccels','pitchaccels')

% save norm data
save('norm_data.mat','rollaccel_norm','pitchaccel_norm','yawaccel_norm','Fenhance_norm',...
    'f_wb_steady')

