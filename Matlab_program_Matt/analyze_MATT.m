clear all, close all
global seqs_plot seqs

warning('off')
if isdir('/home/matt/Dropbox/Analyze') == 1
addpath('/home/matt/Dropbox/Analyze/');
addpath('/home/matt/Dropbox/Analyze/functions/');
addpath('/home/matt/Dropbox/Analyze/variables/');
addpath('/home/matt/Dropbox/Analyze/plotting');
addpath('/home/matt/Dropbox/flytracker_analysis/');
addpath('/home/matt/Dropbox/flytracker_analysis/Plot/');
addpath('/home/matt/Dropbox/flytracker_analysis/Plot/Smoothing_plots/');
end

if isdir('D:\Dropbox\Analyze') == 1
cd('D:\Dropbox\Analyze')
addpath('D:\Dropbox\Analyze');
addpath('D:\Dropbox\Analyze\functions\');
addpath('D:\Dropbox\Analyze\variables\');
addpath('D:\Dropbox\Analyze\plotting');
addpath('D:\Dropbox\flytracker_analysis\');
addpath('D:\Dropbox\flytracker_analysis\Plot\');
addpath('D:\Dropbox\flytracker_analysis\Plot\Smoothing_plots\');
end

% Load pathDB
load_variables
seqs = length(pathDB.x(1,:));
seqs_plot = fliplr(1:81);

% Load list of trimming parameters
load('trims.mat');
trim2smooth; % determines which need to be smoothed
trim2nan;

%%% INPUTS: seqs, pathDB.phi_R, *smooth_peaks, *pathDB.phi_L
wingbeats % set wingbeat defintions based on right wing
%%% OUTPUTS: offset, wingbeats_def

%%% INPUTS: seqs, pathDB.phi_L
wingbeats_L % determine wingbeat defintions from left wing
%%% OUTPUTS: wingbeats_def_L

%%% INPUTS: wingbeats_def, seqs, trims(:,6)
wingbeats_location % determine wingbeat location in reference to stimulus
%%% OUTPUTS: wingbeats_def(updated), wingbeats_loc

%%% INPUTS: winbeats_def_L, seqs, trims(:,6)
wingbeats_location_L % determine left wingbeat location in reference to stimulus
%%% OUTPUTS: wingbeats_def_L(updated), wingbeats_loc_L

%%% INPUTS: pathDB.*, wingbeats_def, wingbeats_def_L, wingbeats_loc
turning_define % use the difference between [theta_R,phi_R ; theta_L,phi_L] to determine turning
%%% OUTPUTS: turning_def

%%% INPUTS: pathDB.*, wingbeats_def, seqs
param_define % parameters trimed to wingbeats, same struture as wingbeats_def
%%% OUTPUTS: stroke_wb_*_MATT, dev_wb_*_MATT, pitch_wb_*_MATT

% plot_surface

% turning_tolerance = 10;
bin_generator
timeseries_generator

% hold on
% plot(pathDB.u_body(422:683,8))
% x = [422:683];
% plot(16*cos(x/length(x)+422*2*pi))
% find(pathDB.u_body(:,8) == max(pathDB.u_body(650:750,8)))

%% Plotting
% plot_gen_steady
% plot_gen_turning
plot_gen_mod

% figure_generator

% x = [1:200];

% for j = 1:8
% for i = 1:length(pitch_wb_L_MATT_bins(1,:,j))
% figure(j)
% hold on
% plot(x,(pitch_wb_L_MATT_bins(:,i,j)-pitch_wb_R_MATT_bins(:,i,j)))
% % plot(x,(dev_wb_L_MATT_bins(:,i,112)),'-g')
% % plot(x,(dev_wb_R_MATT_bins(:,i,112)),'-r')
% hold off
% 
% % figure(2)
% % hold on
% % plot(x,(stroke_wb_L_MATT_bins(:,i,112)))
% % hold off
% % 
% % figure(3)
% % hold on
% % plot(x,(stroke_wb_R_MATT_bins(:,i,112)))
% % hold off
% 
% x = [1:200]' + max(x);
% end
% end
% 
% 
% save
