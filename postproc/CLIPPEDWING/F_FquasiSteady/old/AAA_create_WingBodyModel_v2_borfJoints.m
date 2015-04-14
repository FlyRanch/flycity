clc;
clear;
close all
warning off

addpath('/home/florian/Dropbox/WORK/flytracker/flytracker');
addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/mex/');
addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/core/');
addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/results/');

% load('Wing_kinematics_data_strokeplane47.5deg.mat')
% w_length = nanmean(wing_length);

%% constants
% 8 fem & 7 male, starved for 24h, by FTM. 50%male/female
Mfly = 1.8e-6;
Mg_fly = Mfly*9.81;

w_length = 2.9860;  % [mm] mean wing-length all wingbeats
stroke_angle = 47.5;    % [deg]

rho = 1.225 / 1e9; % [kg/mm^3]
nr_sect = 20;
nr_timepoints = 200;

const.Mfly = Mfly;
const.Mg_fly = Mg_fly;
const.stroke_angle = stroke_angle;

body_model.Mg_fly = Mg_fly;

wing_model.rho = rho;
wing_model.nr_sect = nr_sect;
wing_model.nr_timepoints = nr_timepoints;

%% load model
xh = [0 0 0 0 0 0 1 0 0 0 1 0 0 0 1];
[x_mod,y_mod,z_mod,mod_fit] = load_body_model_basic(xh);

%% scale to average winglength
y_mod_R = y_mod{3};
w_length_basic           = max(y_mod_R(:,1))-min(y_mod_R(:,1));
Rwing = w_length/w_length_basic;

body_model.x_mod = x_mod{1} * Rwing;    % [mm]
body_model.y_mod = y_mod{1} * Rwing;
body_model.z_mod = z_mod{1} * Rwing;

wing_model.x_mod_L = x_mod{2} * Rwing;
wing_model.y_mod_L = y_mod{2} * Rwing;
wing_model.z_mod_L = z_mod{2} * Rwing;

wing_model.x_mod_R = x_mod{3} * Rwing;
wing_model.y_mod_R = y_mod{3} * Rwing;
wing_model.z_mod_R = z_mod{3} * Rwing;

body_model.length           = max(body_model.x_mod(:,1))-min(body_model.x_mod(:,1));
wing_model.length           = max(wing_model.y_mod_R(:,1))-min(wing_model.y_mod_R(:,1));

body_model.Joint_left       = body_model.length.*([0.2021 -0.1055 -0.1477])';
body_model.Joint_right      = body_model.length.*([0.2021 0.1055 -0.1477])';

djoint = 5.5/23 * wing_model.length;
body_model.Joint_left(2)        = -djoint;
body_model.Joint_right(2)       = djoint;

% steady cg
body_model.cg_b = [0; 0; 0];
% Strokeplane2body rotation
body_model.R_strk = [cosd(stroke_angle) 0 sind(stroke_angle); 0 1 0; -sind(stroke_angle) 0 cosd(stroke_angle)];

%% wing sections
y_sect = zeros(3,nr_sect);
y_sect(2,:) = (0.5*(w_length/nr_sect)):(w_length/nr_sect):(w_length-0.5*(w_length/nr_sect));

x_cont_w = -body_model.Joint_right(1) + wing_model.x_mod_R(:,1);
% y_cont_w = -body_model.Joint_right(2) + wing_model.y_mod_R(:,1);
y_cont_w = -body_model.length.*0.1055 + wing_model.y_mod_R(:,1);

x_sect_1_chords = interp1(y_cont_w(1:((end+1)/2)),x_cont_w(1:((end+1)/2)),y_sect(2,:));
x_sect_2_chords = interp1(y_cont_w(((end+1)/2):end),x_cont_w(((end+1)/2):end),y_sect(2,:));

chords = x_sect_1_chords-x_sect_2_chords;

wing_model.y_sect_L             = -y_sect;
wing_model.chords_L             = chords';
wing_model.y_sect_R             = y_sect;
wing_model.chords_R             = chords';

save('BodyWingModel_WingSections_TotalMean.mat','wing_model','body_model','const')