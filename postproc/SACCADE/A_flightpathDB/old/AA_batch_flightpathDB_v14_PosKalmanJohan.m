clear
clc
close all

addpath('/home/florian/Dropbox/WORK/toolbox')
addpath('/home/florian/Dropbox/WORK/toolbox/kalman')
addpath('/home/florian/Dropbox/WORK/toolbox/kalman/EKF')
addpath('/home/florian/Dropbox/WORK/toolbox/kalman/KalmanAll/Kalman')
addpath('/home/florian/Dropbox/WORK/toolbox/kalman/KalmanAll/KPMstats')
addpath('/home/florian/Dropbox/WORK/toolbox/kalman/KalmanAll/KPMtools')

addpath('/home/florian/Dropbox/Projects/flytracker/FlightPathDB/looming_WBkin')

root = cd;
dirs = dir;

% load camera data
load('seq_infoDB.mat')

% load pattern data
load('pattern_teta_temp.mat')

% add or remove roll
roll_on = 1
% roll_on = 0

if roll_on == 1
    savefile = 'flightpathDB_pos_qbodyEKF_INCroll.mat';
else
    savefile = 'flightpathDB_pos_qbodyEKF_NOroll.mat';
end

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

settings.kalman_pos.x = 0.001;
settings.kalman_pos.u = 0.0001;
settings.kalman_pos.a = 0.00001;
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
%             teta = cont_slow_165.teta_cam;
%             t_teta = cont_slow_165.t_cam;
            teta = nan;
            t_teta = nan;
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
        
        if roll_on == 1
            [pathDB,patternDB,settings] = make_flightpathDB_INCroll(pathDB,patternDB,settings,expansion,incq,seq_info,teta);
        else
            [pathDB,patternDB,settings] = make_flightpathDB_NOroll(pathDB,patternDB,settings,expansion,incq,seq_info,teta);
        end
        
        cd(root)
    end
end

cd(root)

if save_plots == 1
    mkdir('track_figs')
    cd('track_figs')
end

[pathDB.pos,pathDB.vel,pathDB.accel]     = batch_kalman_pos_Johan(pathDB.pos_obs,pathDB.t,settings,plotting,save_plots);
cd(root)


%% filter test
% figure
% subplot(3,1,1)
% plot(pathDB.pos(:,83,1))
% hold on
% plot(pathDB.pos_obs(:,83,1)/1000,'r')
% subplot(3,1,2)
% plot(pathDB.pos(:,83,2))
% hold on
% plot(pathDB.pos_obs(:,83,2)/1000,'r')
% subplot(3,1,3)
% plot(pathDB.pos(:,83,3))
% hold on
% plot(pathDB.pos_obs(:,83,3)/1000,'r')
% 
% subplot(3,1,1)
% hold on
% plot(pathDB.pos(:,83,1),'--g')
% subplot(3,1,2)
% hold on
% plot(pathDB.pos(:,83,2),'--g')
% subplot(3,1,3)
% hold on
% plot(pathDB.pos(:,83,3),'--g')
% 
% subplot(3,1,1)
% hold on
% plot(pathDB.pos(:,83,1),'--k')
% subplot(3,1,2)
% hold on
% plot(pathDB.pos(:,83,2),'--k')
% subplot(3,1,3)
% hold on
% plot(pathDB.pos(:,83,3),'--k')
%%

% [pathDB.At,pathDB.An,pathDB.alpha_dot,pathDB.At_hor,pathDB.An_hor,pathDB.alpha_dot_hor] = path_accel_incHorPosNeg(pathDB.vel,pathDB.accel,pathDB.t);

if incq == 1
    [pathDB.qbody,pathDB.omega,settings_EKF] = batch_kalman_qEKF(pathDB.qbody_obs,pathDB.t,settings,plotting,save_plots);

    settings.kalman_body.a = settings_EKF.a;
    settings.kalman_body.b = settings_EKF.b;
    settings.kalman_body.c = settings_EKF.c;
    settings.kalman_body.d = settings_EKF.d;

    [pathDB.roll_global,pathDB.pitch_global,pathDB.yaw_global] = qbody2angles_manual(pathDB.qbody);
    [pathDB.roll_global_obs,pathDB.pitch_global_obs,pathDB.yaw_global_obs] = qbody2angles_manual(pathDB.qbody_obs);
    
    % adjust pitch & roll (global z-axis up, body z-axis down)
    pathDB.pitch = - pathDB.pitch_global;
    roll = pathDB.roll_global;
    roll_temp = pathDB.roll_global;
    roll(roll_temp<0) = roll_temp(roll_temp<0) +180;
    roll(roll_temp>0) = roll_temp(roll_temp>0) -180;
    pathDB.roll = roll;
    pathDB.yaw = pathDB.yaw_global;
    
%     [pathDB.pitch_incroll, pathDB.yaw_incroll] = aerdyn_ref_frame_DB_incroll(pathDB.qbody, pathDB.vel);
%     [pathDB.pitch_NOroll, pathDB.yaw_NOroll] = aerdyn_ref_frame_DB_exclroll(pathDB.qbody, pathDB.vel);

    if calc_qwings == 1

        if saving == 1
            mkdir('wingL_figs')
            cd('wingL_figs')
        end
        [pathDB.qwingL,pathDB.omegaL,settings_EKF] = batch_kalman_qEKF(pathDB.qwingL_obs,pathDB.t,settings,plotting,save_plots);
        settings.kalman_wingL.a = settings_EKF.a;
        settings.kalman_wingL.b = settings_EKF.b;
        settings.kalman_wingL.c = settings_EKF.c;
        settings.kalman_wingL.d = settings_EKF.d;
        cd(root)

        if saving == 1
            mkdir('wingR_figs')
            cd('wingR_figs')
        end
        [pathDB.qwingR,pathDB.omegaR,settings_EKF] = batch_kalman_qEKF(pathDB.qwingR_obs,pathDB.t,settings,plotting,save_plots);
        settings.kalman_wingR.a = settings_EKF.a;
        settings.kalman_wingR.b = settings_EKF.b;
        settings.kalman_wingR.c = settings_EKF.c;
        settings.kalman_wingR.d = settings_EKF.d;
        cd(root)
    end
end

save(savefile,'pathDB','patternDB','settings')


