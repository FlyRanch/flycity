clear
clc
close all

addpath('/home/florian/Dropbox/WORK/flytracker/flightpathDB')
addpath('/home/florian/Dropbox/WORK/flytracker/flightpathDB/kalman')
addpath('/home/florian/Dropbox/WORK/flytracker/flightpathDB/kalman/KalmanAll/Kalman')
addpath('/home/florian/Dropbox/WORK/flytracker/flightpathDB/kalman/KalmanAll/KPMstats')
addpath('/home/florian/Dropbox/WORK/flytracker/flightpathDB/kalman/KalmanAll/KPMtools')

savefile = 'flightpathDB_posq_abNEW.mat';

settings.fps = 7500;
settings.frame_end = 5588;
settings.trigger_frame = 2795;

settings.kalman_pos.x = 0.0001;
settings.kalman_pos.u = 0.00001;
settings.kalman_pos.a = 0.000001;
settings.kalman_pos.R = 10000;
settings.kalman_pos.ss = 9; % state size
settings.kalman_pos.os = 3; % observation size
settings.kalman_pos.F = [1 0 0 1 0 0 0 0 0; 0 1 0 0 1 0 0 0 0; 0 0 1 0 0 1 0 0 0; 0 0 0 1 0 0 1 0 0; 0 0 0 0 1 0 0 1 0; 0 0 0 0 0 1 0 0 1; 0 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 0 1]; 
settings.kalman_pos.H = [1 0 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0 0];

settings.kalman_q.q = 0.1;
settings.kalman_q.q_dot = 0.01;
settings.kalman_q.R = 10000;
settings.kalman_q.ss = 8; % state size
settings.kalman_q.os = 4; % observation size
settings.kalman_q.F = [1 0 0 0 1 0 0 0; 0 1 0 0 0 1 0 0; 0 0 1 0 0 0 1 0; 0 0 0 1 0 0 0 1; 0 0 0 0 1 0 0 0; 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 1]; 
settings.kalman_q.H = [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0];

% initialize DB variable
frames = [1:settings.frame_end]';
pathDB.t = (frames-settings.trigger_frame)/settings.fps;

pathDB.pos_obs=[];
pathDB.qbody_obs = [];
pathDB.qwingL_obs = [];
pathDB.qwingR_obs = [];

root = cd;
dirs = dir;
for d = 3:length(dirs)

    if dirs(d).isdir == 1
        if strcmp(dirs(d).name,'slow_mid')
            expansion.speed = 0;
            expansion.VerPos = 0;
            expansion.HorPos = 0;
            expansion.stepwise = 1;
        elseif strcmp(dirs(d).name,'slow_hi')
            expansion.speed = 0;
            expansion.VerPos = 1;
            expansion.HorPos = 0;
            expansion.stepwise = 1;
        elseif strcmp(dirs(d).name,'slow_low')
            expansion.speed = 0;
            expansion.VerPos = -1;
            expansion.HorPos = 0;
            expansion.stepwise = 1;
        elseif strcmp(dirs(d).name,'fast_mid')
            expansion.speed = 1;
            expansion.VerPos = 0;
            expansion.HorPos = 0;
            expansion.stepwise = 1;
        elseif strcmp(dirs(d).name,'fast_mid_rev')
            expansion.speed = 1;
            expansion.VerPos = 0;
            expansion.HorPos = 180;
            expansion.stepwise = 1;
        elseif strcmp(dirs(d).name,'slow_1pix')
            expansion.speed = 0;
            expansion.VerPos = 0;
            expansion.HorPos = 0;
            expansion.stepwise = 0;
        elseif strcmp(dirs(d).name,'fast_1pix')
            expansion.speed = 1;
            expansion.VerPos = 0;
            expansion.HorPos = 0;
            expansion.stepwise = 0;
        else
            expansion.speed = nan;
            expansion.VerPos = nan;
            expansion.HorPos = nan;
            expansion.stepwise = nan;
        end
        
        cd(dirs(d).name)
        [pathDB,settings] = make_flightpathDB_v6_posq(pathDB,settings,expansion);
        cd(root)
    end
end

mkdir('track_figs')
cd('track_figs')
[pathDB.pos,pathDB.vel,pathDB.accel]     = batch_kalman_pos(pathDB.pos_obs,pathDB.t,settings);
[pathDB.qbody,pathDB.qbody_dot] = batch_kalman_q(pathDB.qbody_obs,pathDB.t,settings);
% figure
% qbody_test
cd(root)

mkdir('wingL_figs')
cd('wingL_figs')
[pathDB.qwingL,pathDB.qwingL_dot] = batch_kalman_q(pathDB.qwingL_obs,pathDB.t,settings);
cd(root)

mkdir('wingR_figs')
cd('wingR_figs')
[pathDB.qwingR,pathDB.qwingR_dot] = batch_kalman_q(pathDB.qwingR_obs,pathDB.t,settings);
cd(root)

[pathDB.At,pathDB.An,pathDB.alpha_dot,pathDB.At_hor,pathDB.An_hor,pathDB.alpha_dot_hor] = path_accel_incHor(pathDB.vel,pathDB.accel,pathDB.t);

[pathDB.roll,pathDB.pitch_global,pathDB.yaw_global] = qbody2angles(pathDB.qbody,pathDB.vel);
[pathDB.pitch_incroll, pathDB.yaw_incroll] = aerdyn_ref_frame_DB_incroll(pathDB.qbody, pathDB.vel);
[pathDB.pitch_NOroll, pathDB.yaw_NOroll] = aerdyn_ref_frame_DB_exclroll(pathDB.qbody, pathDB.vel);

save(savefile,'pathDB','settings')


