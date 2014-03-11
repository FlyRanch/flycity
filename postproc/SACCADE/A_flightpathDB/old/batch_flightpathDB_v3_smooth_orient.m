clear
clc
close all

addpath('/home/florian/Dropbox/WORK/flytracker/flightpathDB')
addpath('/home/florian/Dropbox/WORK/flytracker/flightpathDB/kalman')
addpath('/home/florian/Dropbox/WORK/flytracker/flightpathDB/kalman/KalmanAll/Kalman')
addpath('/home/florian/Dropbox/WORK/flytracker/flightpathDB/kalman/KalmanAll/KPMstats')
addpath('/home/florian/Dropbox/WORK/flytracker/flightpathDB/kalman/KalmanAll/KPMtools')

savefile = 'flightpathDB.mat';

settings.fps = 7500;
settings.frame_end = 5588;
settings.trigger_frame = 2795;

settings.kalman_filter.x = 0.0001;
settings.kalman_filter.u = 0.00001;
settings.kalman_filter.a = 0.000001;
settings.kalman_filter.R = 10000;
settings.kalman_filter.ss = 9; % state size
settings.kalman_filter.os = 3; % observation size
settings.kalman_filter.F = [1 0 0 1 0 0 0 0 0; 0 1 0 0 1 0 0 0 0; 0 0 1 0 0 1 0 0 0; 0 0 0 1 0 0 1 0 0; 0 0 0 0 1 0 0 1 0; 0 0 0 0 0 1 0 0 1; 0 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 0 1]; 
settings.kalman_filter.H = [1 0 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0 0];

setting.smooth_param = .99999;

% initialize DB variable
frames = [1:settings.frame_end]';
pathDB.t = (frames-settings.trigger_frame)/settings.fps;

pathDB.x_obs=[];
pathDB.y_obs=[];
pathDB.z_obs=[];

pathDB.yaw_obs = [];
pathDB.pitch_obs = [];
pathDB.roll_obs = [];

pathDB.wingL1_obs = [];
pathDB.wingL2_obs = [];
pathDB.wingL3_obs = [];

pathDB.wingR1_obs = [];
pathDB.wingR2_obs = [];
pathDB.wingR3_obs = [];

root = cd;
dirs = dir;
for d = 3:length(dirs)

    if dirs(d).isdir == 1
        if strcmp(dirs(d).name,'slow_mid')
            expansion.speed = 0;
            expansion.VerPos = 0;
            expansion.HorPos = 0;
        elseif strcmp(dirs(d).name,'slow_hi')
            expansion.speed = 0;
            expansion.VerPos = 1;
            expansion.HorPos = 0;
        elseif strcmp(dirs(d).name,'slow_low')
            expansion.speed = 0;
            expansion.VerPos = -1;
            expansion.HorPos = 0;
        elseif strcmp(dirs(d).name,'fast_mid')
            expansion.speed = 1;
            expansion.VerPos = 0;
            expansion.HorPos = 0;
        elseif strcmp(dirs(d).name,'fast_mid_rev')
            expansion.speed = 1;
            expansion.VerPos = 0;
            expansion.HorPos = 180;
        else
            expansion.speed = nan;
            expansion.pos = nan;
        end
        
        cd(dirs(d).name)
        [settings,pathDB] = make_flightpathDB_v5_NEWorientations(settings,pathDB,expansion);
        cd(root)
    end
end

mkdir('track_figs')
cd('track_figs')
[settings,pathDB] = batch_kalman_track_subsets(settings,pathDB)
cd(root)
% save(savefile,'pathDB','settings')        
% 
% response_results_kalmandata_v2_subsets
% plot_alltracks_v3_subsets

U = sqrt(pathDB.u.^2 + pathDB.v.^2 + pathDB.w.^2 );
A = sqrt(pathDB.ax.^2 + pathDB.ay.^2 + pathDB.az.^2 );
ascend = atand(pathDB.w ./ sqrt(pathDB.u.^2 + pathDB.v.^2));


t = pathDB.t;
yaw = pathDB.yaw_obs;
pitch = pathDB.pitch_obs;
roll = pathDB.roll_obs;

wingL1 = pathDB.wingL1_obs;
wingL2 = pathDB.wingL2_obs;
wingL3 = pathDB.wingL3_obs;

wingR1 = pathDB.wingR1_obs;
wingR2 = pathDB.wingR2_obs;
wingR3 = pathDB.wingR3_obs;

for i=1:size(pitch,2)
    yaw_now = yaw(:,i);
    pitch_now = pitch(:,i);
    roll_now = roll(:,i);

    pp = csaps(t,yaw_now,.99999);
    yaw_now(isnan(yaw_now)==0) = fnval(pp,t(isnan(yaw_now)==0));
    yaw(:,i) = yaw_now;    

    pp = csaps(t,pitch_now,.99999);
    pitch_now(isnan(pitch_now)==0) = fnval(pp,t(isnan(pitch_now)==0));
    pitch(:,i) = pitch_now;    
    
    pp = csaps(t,roll_now,.99999);
    roll_now(isnan(roll_now)==0) = fnval(pp,t(isnan(roll_now)==0));
    roll(:,i) = roll_now;    
end

pitch_true = pitch - ascend;

pathDB.yaw = yaw;
pathDB.pitch = pitch;
pathDB.roll = roll;

for i=1:length(pathDB.u)
    for j=1:size(pathDB.u,2)
        U = [pathDB.u(i,j) pathDB.v(i,j) pathDB.w(i,j)]';
        A = [pathDB.ax(i,j) pathDB.ay(i,j) pathDB.az(i,j)]';

        At(i,j) = dot(U,A) / norm(U);
        An(i,j) = sqrt(norm(A)^2 - At(i,j)^2);
    end

end

pathDB.At = At;
pathDB.An = An;

save(savefile,'pathDB','settings')        