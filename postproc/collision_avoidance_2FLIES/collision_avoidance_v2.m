% collision avoidance 
clear
clc
close all

plotting = 1;
saving_plots = 0;
savefile = 'flightpathDB_pos.mat';

addpath('/home/florian/Dropbox/WORK/toolbox/kalman')
addpath('/home/florian/Dropbox/WORK/toolbox/kalman/EKF')
addpath('/home/florian/Dropbox/WORK/toolbox/kalman/KalmanAll/Kalman')
addpath('/home/florian/Dropbox/WORK/toolbox/kalman/KalmanAll/KPMstats')
addpath('/home/florian/Dropbox/WORK/toolbox/kalman/KalmanAll/KPMtools')
addpath('/home/florian/Dropbox/WORK/flytracker/postproc/A_flightpathDB')

load('kine_20130122_S0006_bodywidth.mat')

b1 = (data.kine.FLYBODY1.data.length(data.kine.FLYBODY1.data.length~=0));
b2 = (data.kine.FLYBODY2.data.length(data.kine.FLYBODY2.data.length~=0));
pathDB.b = mean([b1 b2]) /1000;


load('kine_20130122_S0006_2flies.mat')
t = [0:5587]'/7500;

% frames for length calc
frames1 = find(data.kine.FLYBODY1.data.length~=0)';
frames2 = find(data.kine.FLYBODY2.data.length~=0)';
l1 = mean(data.kine.FLYBODY1.data.length(frames1));
l2 = mean(data.kine.FLYBODY2.data.length(frames2));

% frames for pos calc
x1(:,:) = data.kine.FLYBODY1.data.coords(1,1,:);
x2(:,:) = data.kine.FLYBODY2.data.coords(1,1,:);
frames1 = find(x1~=0);
frames2 = find(x2~=0);

if length(frames1) ~= length(frames2)
    something = 'wrong'
    break
end
frames = frames1;
pathDB.t = t(frames) - min(t(frames));

% observed pos matrix (head & tail of fly1 & fly2)
pos1 = data.kine.FLYBODY1.data.coords(:,:,frames);
pos2 = data.kine.FLYBODY2.data.coords(:,:,frames);
pathDB.pos_obs = [permute(pos1,[3 1 2]) permute(pos2,[3 1 2])];

% kalman filter
settings.kalman_pos.x = 1;
settings.kalman_pos.u = 0.1;
settings.kalman_pos.a = 0.01;
settings.kalman_pos.R = 1;
settings.kalman_pos.ss = 9; % state size
settings.kalman_pos.os = 3; % observation size
settings.kalman_pos.F = [1 0 0 1 0 0 0 0 0; 0 1 0 0 1 0 0 0 0; 0 0 1 0 0 1 0 0 0; 0 0 0 1 0 0 1 0 0; 0 0 0 0 1 0 0 1 0; 0 0 0 0 0 1 0 0 1; 0 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 0 1]; 
settings.kalman_pos.H = [1 0 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0 0];

[pathDB.pos,pathDB.vel,pathDB.accel]     = batch_kalman_pos(pathDB.pos_obs,pathDB.t,settings,plotting,saving_plots);
[pathDB.At,pathDB.An,pathDB.alpha_dot,pathDB.At_hor,pathDB.An_hor,pathDB.alpha_dot_hor] = path_accel_incHorPosNeg(pathDB.vel,pathDB.accel,pathDB.t);

% distance between heads
pathDB.d = sqrt( (pathDB.pos(:,1,1) - pathDB.pos(:,3,1)).^2 + (pathDB.pos(:,1,2) - pathDB.pos(:,3,2)).^2 + (pathDB.pos(:,1,3) - pathDB.pos(:,3,3)).^2 );
pathDB.teta = 2* atand(pathDB.b/2./pathDB.d);

figure
plot(pathDB.t,pathDB.teta)

save(savefile,'pathDB','settings')







