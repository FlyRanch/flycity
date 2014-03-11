clear all
close all
clc

%--------------------------------------------------------------------------
%
%   Execution script: Flycity Analysis Program
%
%   By: Johan Melis
%   Organization: Dickinson Fly Lab
%   Date: 13-08-2013
%
%--------------------------------------------------------------------------


% Submit program locations to path:

% add_paths = { 'C:/Users/Johan/Desktop/Flycity_Analysis_Program'; ...
%               'C:/Users/Johan/Desktop/Flycity_Analysis_Program/Main'; ...
%               'C:/Users/Johan/Desktop/Flycity_Analysis_Program/Filter'; ...
%               'C:/Users/Johan/Desktop/Flycity_Analysis_Program/Body_and_wing_model'; ...
%               'C:/Users/Johan/Desktop/Flycity_Analysis_Program/Select_wingbeat'; ...
%               'C:/Users/Johan/Desktop/Flycity_Analysis_Program/Body_and_wing_model/Body_model'; ...
%               'C:/Users/Johan/Desktop/Flycity_Analysis_Program/Body_and_wing_model/Body_model/mex'; ...
%               'C:/Users/Johan/Desktop/Flycity_Analysis_Program/Polynomial_Fit'; ...
%               'C:/Users/Johan/Desktop/Flycity_Analysis_Program/Dynamic_Model'; ...
%               'C:/Users/Johan/Desktop/Flycity_Analysis_Program/Dynamic_Model/Stationary_simulation'; ...
%               'C:/Users/Johan/Desktop/Flycity_Analysis_Program/Dynamic_Model/RK_simulation'; ...
%              };

add_paths = { 'C:/Users/Johan/Desktop/Flycity_Analysis_Program'; ...
              'C:/Users/Johan/Desktop/Flycity_Analysis_Program/Main'; ...
              'C:/Users/Johan/Desktop/Flycity_Analysis_Program/Filter'; ...
              'C:/Users/Johan/Desktop/Flycity_Analysis_Program/Body_and_wing_model'; ...
              'C:/Users/Johan/Desktop/Flycity_Analysis_Program/Select_wingbeat'; ...
              'C:/Users/Johan/Desktop/Flycity_Analysis_Program/Body_and_wing_model/Body_model'; ...
              'C:/Users/Johan/Desktop/Flycity_Analysis_Program/Body_and_wing_model/Body_model/mex'; ...
              'C:/Users/Johan/Desktop/Flycity_Analysis_Program/Polynomial_Fit'; ...
              'C:/Users/Johan/Desktop/Flycity_Analysis_Program/Dynamic_Model'; ...
              'C:/Users/Johan/Desktop/Flycity_Analysis_Program/Dynamic_Model/Stationary_simulation'; ...
              'C:/Users/Johan/Desktop/Flycity_Analysis_Program/Dynamic_Model/RK_simulation'; ...
              'C:/Users/Johan/Desktop/Flycity_Analysis_Program/Control_analysis'; ...
              'C:/Users/Johan/Desktop/Flycity_Analysis_Program/Robofly_experiments'; ...
              'C:/Users/Johan/Desktop/Flycity_Analysis_Program/Plot'; ...
             };
         
addpath(char(add_paths(1)))
addpath(char(add_paths(2)))
addpath(char(add_paths(3)))
addpath(char(add_paths(4)))
addpath(char(add_paths(5)))
addpath(char(add_paths(6)))
addpath(char(add_paths(7)))
addpath(char(add_paths(8)))
addpath(char(add_paths(9)))
addpath(char(add_paths(10)))
addpath(char(add_paths(11)))
addpath(char(add_paths(12)))
addpath(char(add_paths(13)))
addpath(char(add_paths(14)))


% Input:

% pathDB:

    pathDB = {};


% Create empty structure settings:

    settings = {};

    % Data location:

    settings.data_loc = {'C:/Users/Johan/Documents/Thesis_data/Sequences_14_06_2013'};
%     settings.data_loc = {'C:/Users/Johan/Documents/Thesis_data/test_seq'};
%     settings.data_loc = {'C:/Users/Johan/Documents/Thesis_data/test_seq2'};

    % Plot location:
    
    settings.plot_loc = {'C:/Users/Johan/Documents/Thesis_data/Plots'};

    % High-speed camera settings:

    settings.fps = 7500;
    settings.frame_end = 5588;
    settings.trigger_frame = 2795;    
    
    % Set time vector:

    frames = [1:settings.frame_end]';
    pathDB.t = (frames-settings.trigger_frame)/settings.fps; % [s]
    pathDB.dt = pathDB.t(2)-pathDB.t(1); % [s]

    % Set Kalman Filter and Extended Kalman filter settings:

    settings.filter_set.xyz = [pathDB.dt^2 pathDB.dt^2 pathDB.dt^2 1 1 1 1/pathDB.dt^2 1/pathDB.dt^2 1/pathDB.dt^2];
    settings.filter_set.qB1 = [ 0 0 0 0 0 0 0.003 0.003 0.003 0.003 ];
    settings.filter_set.qB2 = [ 0.0001/pathDB.dt^2 0.0001/pathDB.dt^2 0.0001/pathDB.dt^2 0.0001 0.0001 0.0001 1 1 1 1 ];
    settings.filter_set.qW1 = [ 0 0 0 0 0 0 1 1 1 1 ];
    settings.filter_set.qW2 = [ 0 0 0 1/pathDB.dt^2 1/pathDB.dt^2 1/pathDB.dt^2 1 1 1 1 ];
    
    % Set pitch angle (w.r.t. body axis) of strokeplane reference frame:
    
    settings.beta_strk = -(55/180)*pi; % [rad]
    
    % Set the number of chord sections on a wing:
    
    settings.nr_chord_sect = 20;
    
    % Give air density, gravitation constant, body mass density and wing
    % mass density:
    
    settings.rho_air            = 1.188e-9; % [kg/mm^3]
    settings.C_mfly             = 1.85e-6;
    settings.rho_cuticle        = 1200e-9; % [kg/mm^3]
    settings.h_wing             = 5.4e-4; % Non dimensional wing thickness [ t_w / L_w ]
    settings.g                  = 9810; % [mm/s^2]
    
    % Set the order of the polynomials used:
    
    settings.n_pol_theta =  12;    % Order of used polynomials for theta
    settings.n_pol_eta =    14;    % Order of used polynomials for eta
    settings.n_pol_phi =    10;    % Order of used polynomials for phi
    
    
    % Set maneuver treshold for acceleration and angular velocity:
    
    settings.w_tresh = 6*pi; % rad/s
    settings.a_tresh = 4000; % mm/s^2
    
    % Give sequence which enhance a trigger shift:
    
    settings.trigger_shift = {  '20130205_S0004'; ...
                                '20130205_S0008'; ...
                                '20130206_S0003'; ...
                                '20130206_S0006'; ...
                                '20130206_S0007'; ...
                                '20130206_S0008'; ...
                                '20130207_S0002'; ...
                                '20130208_S0005' };
    

%--------------------------------------------------------------------------

cd(char(settings.data_loc))

% Load raw data into pathDB:

if exist('pathDB1.mat')==2
    temp = load('pathDB1.mat');
    
        pathDB.raw.xyz = temp.xyz;
        pathDB.raw.qB = temp.qB;
        pathDB.raw.qL = temp.qL;
        pathDB.raw.qR = temp.qR;
        
        settings.sequence_names = temp.sequence_names;
        settings.nr_of_seq = temp.nr_of_seq;
        settings.start_stop = temp.start_stop;
        
    'loaded'
    clear temp
else
    read_raw_data(settings,pathDB);
    temp = load('pathDB1.mat');
    
        pathDB.raw.xyz = temp.xyz;
        pathDB.raw.qB = temp.qB;
        pathDB.raw.qL = temp.qL;
        pathDB.raw.qR = temp.qR;
        
        settings.sequence_names = temp.sequence_names;
        settings.nr_of_seq = temp.nr_of_seq;
        settings.start_stop = temp.start_stop;
    
    'created + loaded'
    clear temp
end


% Filter body and wing kinematics:

if exist('pathDB2.mat')==2
    temp = load('pathDB2.mat');
    
        pathDB.filt.xyz = temp.xyz;
        pathDB.filt.uvw = temp.uvw;
        pathDB.filt.a_xyz = temp.a_xyz;
        pathDB.filt.qB = temp.qB;
        pathDB.filt.wB = temp.wB;
        pathDB.filt.qL = temp.qL;
        pathDB.filt.wL = temp.wL;
        pathDB.filt.qR = temp.qR;
        pathDB.filt.wR = temp.wR;
        
    'loaded'
    clear temp
else
    filter_data(settings,pathDB);
    temp = load('pathDB2.mat');
    
        pathDB.filt.xyz = temp.xyz;
        pathDB.filt.uvw = temp.uvw;
        pathDB.filt.a_xyz = temp.a_xyz;
        pathDB.filt.qB = temp.qB;
        pathDB.filt.wB = temp.wB;
        pathDB.filt.qL = temp.qL;
        pathDB.filt.wL = temp.wL;
        pathDB.filt.qR = temp.qR;
        pathDB.filt.wR = temp.wR;
    
    'created + loaded'
    clear temp
end


% Create rotation matrices, body data & model and wing kinematics:

if exist('pathDB3.mat')==2
    temp = load('pathDB3.mat');
    
        pathDB.rot_mat = temp.rot_mat;
        pathDB.body_model = temp.body_model;
        pathDB.wing_model = temp.wing_model;
        pathDB.strkpln_kin = temp.strkpln_kin;
        pathDB.wing_kin = temp.wing_kin;
        
    'loaded'
    clear temp
else
    body_and_wing_model(settings,pathDB);
    temp = load('pathDB3.mat');
    
        pathDB.rot_mat = temp.rot_mat;
        pathDB.body_model = temp.body_model;
        pathDB.wing_model = temp.wing_model;
        pathDB.strkpln_kin = temp.strkpln_kin;
        pathDB.wing_kin = temp.wing_kin;
    
    'created + loaded'
    clear temp
end


% Find wingbeat duration, downstroke and upstroke duration:

% Create rotation matrices, body data & model and wing kinematics:

if exist('pathDB4.mat')==2
    temp = load('pathDB4.mat');
    
        pathDB.wingbeats = temp.wingbeats;
        
    'loaded'
    clear temp
else
    Select_wingbeat(settings,pathDB)
    temp = load('pathDB4.mat');
    
        pathDB.wingbeats = temp.wingbeats;
    
    'created + loaded'
    clear temp
end


% Find polynomial fits for the wingbeats of the sequences and save them per
% sequence and per random accessible wingbeat:

if exist('pathDB5.mat')==2
    temp = load('pathDB5.mat');
    
        pathDB.poly_fit = temp.poly_fit;
        pathDB.rand_wbs = temp.rand_wbs;
        pathDB.maneuver = temp.maneuver;
        
    'loaded'
    clear temp
else
    Polynomial_Fit( settings, pathDB )
    temp = load('pathDB5.mat');
    
        pathDB.poly_fit = temp.poly_fit;
        pathDB.rand_wbs = temp.rand_wbs;
        pathDB.maneuver = temp.maneuver;
    
    'created + loaded'
    clear temp
end

% Find maneuvering wing kinematics by regression of the deviation
% coefficients vs the non-dimensional forces and moments:

if exist('pathDB6.mat')==2
    temp = load('pathDB6.mat');
    
        pathDB.maneuver.c_fit_ax = temp.c_fit_ax;
        pathDB.maneuver.c_fit_ay = temp.c_fit_ay;
        pathDB.maneuver.c_fit_az = temp.c_fit_az;
        pathDB.maneuver.c_fit_wx = temp.c_fit_wx;
        pathDB.maneuver.c_fit_wy = temp.c_fit_wy;
        pathDB.maneuver.c_fit_wz = temp.c_fit_wz;
        
    'loaded'
    clear temp
else
    control_analysis( settings, pathDB )
    temp = load('pathDB6.mat');
    
        pathDB.maneuver.c_fit_ax = temp.c_fit_ax;
        pathDB.maneuver.c_fit_ay = temp.c_fit_ay;
        pathDB.maneuver.c_fit_az = temp.c_fit_az;
        pathDB.maneuver.c_fit_wx = temp.c_fit_wx;
        pathDB.maneuver.c_fit_wy = temp.c_fit_wy;
        pathDB.maneuver.c_fit_wz = temp.c_fit_wz;
    
    'created + loaded'
    clear temp
end

% Create the wingkinematics for the Robofly experiments:

if exist('pathDB7.mat')==2
    temp = load('pathDB7.mat');

        pathDB.Robofly.Fx_forward   = temp.Fx_forward;
        pathDB.Robofly.Fx_back      = temp.Fx_back;
        pathDB.Robofly.Fy           = temp.Fy;
        pathDB.Robofly.Fz_down      = temp.Fz_down;
        pathDB.Robofly.Fz_up        = temp.Fz_up;
        pathDB.Robofly.Mx           = temp.Mx;
        pathDB.Robofly.My_up        = temp.My_up;
        pathDB.Robofly.My_down      = temp.My_down;
        pathDB.Robofly.Mz           = temp.Mz;
        
    'loaded'
    clear temp
else
    Robofly_wing_kinematics( settings, pathDB )
    temp = load('pathDB7.mat');
    
        pathDB.Robofly.Fx_forward   = temp.Fx_forward;
        pathDB.Robofly.Fx_back      = temp.Fx_back;
        pathDB.Robofly.Fy           = temp.Fy;
        pathDB.Robofly.Fz_down      = temp.Fz_down;
        pathDB.Robofly.Fz_up        = temp.Fz_up;
        pathDB.Robofly.Mx           = temp.Mx;
        pathDB.Robofly.My_up        = temp.My_up;
        pathDB.Robofly.My_down      = temp.My_down;
        pathDB.Robofly.Mz           = temp.Mz;
    
    'created + loaded'
    clear temp
end

% Robofly_wing_kinematics( settings, pathDB )

Conservation_test_old( settings, pathDB )

% Conservation_test( settings, pathDB )

% RK_simulation_test( settings, pathDB )

% Test_dynamic_model( settings, pathDB )

