clear
clc
close all

addpath('/home/florian/Dropbox/WORK/toolbox')
addpath('/home/florian/Dropbox/WORK/toolbox/kalman')
addpath('/home/florian/Dropbox/WORK/toolbox/kalman/EKF')
addpath('/home/florian/Dropbox/WORK/toolbox/kalman/KalmanAll/Kalman')
addpath('/home/florian/Dropbox/WORK/toolbox/kalman/KalmanAll/KPMstats')
addpath('/home/florian/Dropbox/WORK/toolbox/kalman/KalmanAll/KPMtools')

root = cd;
dirs = dir;

% load camera data
load('seq_infoDB.mat')

% load pattern data
load('pattern_teta_temp.mat')

savefile = 'flightpathDB_pos_qbodyEKF_NOroll.mat';

% add quarternions?
incq = 0;
incq = 1;

% add wing quaternions?
    calc_qwings = 0;
%     calc_qwings = 1;

% plot data?
plotting = 0;
% plotting = 1;

% save plots?
save_plots = 0;
% save_plots = 1;

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

pathDB.pos_obs=[];
patternDB.teta=[];

if incq == 1
%     settings.kalman_q.q = 0.1;
%     settings.kalman_q.q_dot = 0.01;
%     settings.kalman_q.R = 10000;
%     settings.kalman_q.ss = 8; % state size
%     settings.kalman_q.os = 4; % observation size
%     settings.kalman_q.F = [1 0 0 0 1 0 0 0; 0 1 0 0 0 1 0 0; 0 0 1 0 0 0 1 0; 0 0 0 1 0 0 0 1; 0 0 0 0 1 0 0 0; 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 1]; 
%     settings.kalman_q.H = [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0];

    pathDB.qbody_obs = [];
    pathDB.qwingL_obs = [];
    pathDB.qwingR_obs = [];
end

% initialize DB variable
frames = [1:settings.frame_end]';
pathDB.t = (frames-settings.trigger_frame)/settings.fps;

for d = 3:length(dirs)

    if dirs(d).isdir == 1
        if strcmp(dirs(d).name,'slow_mid')
            expansion.speed = 0;
            expansion.VerPos = 0;
            expansion.HorPos = 0;
            expansion.stepwise = 1;
            expansion.maxangle = 165;
            expansion.arista = 0;
            expansion.OFF = 0;
            teta = step_slow.teta_cam;
            t_teta = step_slow.t_cam;
        elseif strcmp(dirs(d).name,'slow_hi')
            expansion.speed = 0;
            expansion.VerPos = 1;
            expansion.HorPos = 0;
            expansion.stepwise = 1;
            expansion.maxangle = 165;
            expansion.arista = 0;
            expansion.OFF = 0;
            teta = step_slow.teta_cam;
            t_teta = step_slow.t_cam;
        elseif strcmp(dirs(d).name,'slow_low')
            expansion.speed = 0;
            expansion.VerPos = -1;
            expansion.HorPos = 0;
            expansion.stepwise = 1;
            expansion.maxangle = 165;
            expansion.arista = 0;
            expansion.OFF = 0;
            teta = step_slow.teta_cam;
            t_teta = step_slow.t_cam;
        elseif strcmp(dirs(d).name,'fast_mid')
            expansion.speed = 1;
            expansion.VerPos = 0;
            expansion.HorPos = 0;
            expansion.stepwise = 1;
            expansion.maxangle = 165;
            expansion.arista = 0;
            expansion.OFF = 0;
            teta = step_fast.teta_cam;
            t_teta = step_fast.t_cam;
        elseif strcmp(dirs(d).name,'fast_mid_rev')
            expansion.speed = 1;
            expansion.VerPos = 0;
            expansion.HorPos = 180;
            expansion.stepwise = 1;
            expansion.maxangle = 165;
            expansion.arista = 0;
            expansion.OFF = 0;
            teta = step_fast.teta_cam;
            t_teta = step_fast.t_cam;
        elseif strcmp(dirs(d).name,'slow_1pix_max165')
            expansion.speed = 0;
            expansion.VerPos = 0;
            expansion.HorPos = 0;
            expansion.stepwise = 0;
            expansion.maxangle = 165;
            expansion.arista = 0;
            expansion.OFF = 0;
            teta = cont_slow_165.teta_cam;
            t_teta = cont_slow_165.t_cam;
       elseif strcmp(dirs(d).name,'fast_1pix_max165')
            expansion.speed = 2;
            expansion.VerPos = 0;
            expansion.HorPos = 0;
            expansion.stepwise = 0;
            expansion.maxangle = 165;
            expansion.arista = 0;
            expansion.OFF = 0;
            teta = cont_fast_165.teta_cam;
            t_teta = cont_fast_165.t_cam;
       elseif strcmp(dirs(d).name,'fast_1pix_max165_ARISTACLIPPED')
            expansion.speed = 2;
            expansion.VerPos = 0;
            expansion.HorPos = 0;
            expansion.stepwise = 0;
            expansion.maxangle = 165;
            expansion.arista = 1;
            expansion.OFF = 0;
            teta = cont_fast_165.teta_cam;
            t_teta = cont_fast_165.t_cam;
       elseif strcmp(dirs(d).name,'medium_1pix_max64')
            expansion.speed = 1;
            expansion.VerPos = 0;
            expansion.HorPos = 0;
            expansion.stepwise = 0;
            expansion.maxangle = 64;
            expansion.arista = 0;
            expansion.OFF = 0;
            teta = cont_medium_64.teta_cam;
            t_teta = cont_medium_64.t_cam;
       elseif strcmp(dirs(d).name,'fast_1pix_max64')
            expansion.speed = 2;
            expansion.VerPos = 0;
            expansion.HorPos = 0;
            expansion.stepwise = 0;
            expansion.maxangle = 64;
            expansion.arista = 0;
            expansion.OFF = 0;
            teta = cont_fast_64.teta_cam;
            t_teta = cont_fast_64.t_cam;
       elseif strcmp(dirs(d).name,'fast_1pix_max64OFF')
            expansion.speed = 2;
            expansion.VerPos = 0;
            expansion.HorPos = 0;
            expansion.stepwise = 0;
            expansion.maxangle = 64;
            expansion.arista = 0;
            expansion.OFF = 1;
            teta = cont_fast_64_off.teta_cam;
            t_teta = cont_fast_64_off.t_cam;
       elseif strcmp(dirs(d).name,'fast_1pix_max32')
            expansion.speed = 2;
            expansion.VerPos = 0;
            expansion.HorPos = 0;
            expansion.stepwise = 0;
            expansion.maxangle = 32;
            expansion.arista = 0;
            expansion.arista = 0;
            expansion.OFF = 0;
            teta = cont_fast_32.teta_cam;
            t_teta = cont_fast_32.t_cam;
       elseif strcmp(dirs(d).name,'fast_1pix_max32OFF')
            expansion.speed = 2;
            expansion.VerPos = 0;
            expansion.HorPos = 0;
            expansion.stepwise = 0;
            expansion.maxangle = 32;
            expansion.arista = 0;
            expansion.OFF = 1;
            teta = cont_fast_32_off.teta_cam;
            t_teta = cont_fast_32_off.t_cam;
       elseif strcmp(dirs(d).name,'fast_1pix_max16')
            expansion.speed = 2;
            expansion.VerPos = 0;
            expansion.HorPos = 0;
            expansion.stepwise = 0;
            expansion.maxangle = 16;
            expansion.arista = 0;
            expansion.OFF = 0;
            teta = cont_fast_16.teta_cam;
            t_teta = cont_fast_16.t_cam;
       elseif strcmp(dirs(d).name,'fast_1pix_max16OFF')
            expansion.speed = 2;
            expansion.VerPos = 0;
            expansion.HorPos = 0;
            expansion.stepwise = 0;
            expansion.maxangle = 16;
            expansion.arista = 0;
            expansion.OFF = 1;
            teta = cont_fast_16_off.teta_cam;
            t_teta = cont_fast_16_off.t_cam;
       else
            expansion.speed = nan;
            expansion.VerPos = nan;
            expansion.HorPos = nan;
            expansion.stepwise = nan;
            expansion.maxangle = nan;
            expansion.arista = nan;
            expansion.OFF = nan;
            teta = nan;
            t_teta = nan;
        end
        
        cd(dirs(d).name)
        [pathDB,patternDB,settings] = make_flightpathDB_v11_NOroll(pathDB,patternDB,settings,expansion,incq,seq_info,teta);
        cd(root)
    end
end

cd(root)
% mkdir('track_figs')
% cd('track_figs')
[pathDB.pos,pathDB.vel,pathDB.accel]     = batch_kalman_pos(pathDB.pos_obs,pathDB.t,settings,plotting,save_plots);
% cd(root)

% [pathDB.At,pathDB.An,pathDB.alpha_dot,pathDB.At_hor,pathDB.An_hor,pathDB.alpha_dot_hor] = path_accel_incHorPosNeg(pathDB.vel,pathDB.accel,pathDB.t);

if incq == 1
    [pathDB.qbody,pathDB.omega] = batch_kalman_qEKF(pathDB.qbody_obs,pathDB.t,settings,plotting);

    [pathDB.roll,pathDB.pitch_global,pathDB.yaw_global] = qbody2angles_manual(pathDB.qbody);
    [pathDB.roll_obs,pathDB.pitch_global_obs,pathDB.yaw_global_obs] = qbody2angles_manual(pathDB.qbody_obs);
    
%     [pathDB.pitch_incroll, pathDB.yaw_incroll] = aerdyn_ref_frame_DB_incroll(pathDB.qbody, pathDB.vel);
%     [pathDB.pitch_NOroll, pathDB.yaw_NOroll] = aerdyn_ref_frame_DB_exclroll(pathDB.qbody, pathDB.vel);
    
    if calc_qwings == 1

        mkdir('wingL_figs')
        cd('wingL_figs')
        [pathDB.qwingL,pathDB.omegaL] = batch_kalman_qEKF(pathDB.qwingL_obs,pathDB.t,settings,plotting);
        cd(root)

        mkdir('wingR_figs')
        cd('wingR_figs')
        [pathDB.qwingR,pathDB.omegaR] = batch_kalman_qEKF(pathDB.qwingR_obs,pathDB.t,settings,plotting);
        cd(root)
    end
end

save(savefile,'pathDB','patternDB','settings')


