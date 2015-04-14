%% make MOD&normDB
clear
clc
close all

% load norm data
load('norm_data_torque.mat')

% % mod locations
% Fenhanses = [0 .6]'
% rolltorques = [0 3]'
% % pitchtorques = [0 1.5 3]'
% pitchtorques = [0 2]'
% yawtorques = [0 2]'
% 
% FenhMods = Fenhanses/Fenhance_norm
% RollMods = rolltorques/rolltorque_norm
% % PitchMods = pitchtorques/pitchtorque_norm
% PitchMods = [-max(pitchtorques) 0 max(pitchtorques)]'/pitchtorque_norm
% YawMods = yawtorques/yawtorque_norm

FDB = 'FnMqs_Fenhance.mat'
rollDB = 'FnMqs_RollTorque.mat'
pitchDB = 'FnMqs_PitchTorque.mat'
yawDB = 'FnMqs_YawTorque.mat'
RaxisDB = 'FnMqs_TorqueAxisR.mat'

% save('MOD_norm_data_QSmodel.mat')
save('norm_data_torque.mat')
