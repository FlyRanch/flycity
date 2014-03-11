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


% initialize DB variable
frames = [1:settings.frame_end]';
pathDB.t = (frames-settings.trigger_frame)/settings.fps;

pathDB.x_obs=[];
pathDB.y_obs=[];
pathDB.z_obs=[];

root = cd;
dirs = dir;
for d = 3:length(dirs)

    if dirs(d).isdir == 1
        if strcmp(dirs(d).name,'slow_mid')
            expansion.speed = 0;
            expansion.pos = 0;
        elseif strcmp(dirs(d).name,'slow_hi')
            expansion.speed = 0;
            expansion.pos = 1;
        elseif strcmp(dirs(d).name,'slow_low')
            expansion.speed = 0;
            expansion.pos = -1;
        elseif strcmp(dirs(d).name,'fast_mid')
            expansion.speed = 1;
            expansion.pos = 0;
        else
            expansion.speed = nan;
            expansion.pos = nan;
        end
        
        cd(dirs(d).name)
        [settings,pathDB] = make_flightpathDB_v3_subsets(settings,pathDB,expansion);
        cd(root)
    end
end

mkdir('track_figs')
cd('track_figs')
[settings,pathDB] = batch_kalmar_track_subsets(settings,pathDB)
cd(root)
save(savefile,'pathDB','settings')        

response_results_kalmandata_v2_subsets
plot_alltracks_v3_subsets

