%% make MOD&normDB
clear
clc
close all

% load norm data
load('norm_data.mat')

% % mod locations
% Fenhanses = [0 .6]'
% rollaccels = [0 3]'
% % pitchaccels = [0 1.5 3]'
% pitchaccels = [0 2]'
% yawaccels = [0 2]'
% 
% FenhMods = Fenhanses/Fenhance_norm
% RollMods = rollaccels/rollaccel_norm
% % PitchMods = pitchaccels/pitchaccel_norm
% PitchMods = [-max(pitchaccels) 0 max(pitchaccels)]'/pitchaccel_norm
% YawMods = yawaccels/yawaccel_norm

FDB = 'FnMqs_Fenhance.mat'
rollDB = 'FnMqs_RollAccel.mat'
pitchDB = 'FnMqs_PitchAccel.mat'
yawDB = 'FnMqs_YawAccel.mat'

save('MOD_norm_data_QSmodel.mat')
