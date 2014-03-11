clc
clear
% close all

addpath('/home/florian/Dropbox/WORK/toolbox')
addpath('/home/florian/Dropbox/WORK/toolbox/plot2svg')
addpath('/home/florian/Dropbox/WORK/toolbox/CircStat')

% load WBmod data
list = dir('WBmod*.mat')
for i = 1:length(list)
    load(list(i).name)
end
steady_name = dir('WBdataset_steady_*')
% steady_name = dir('WBdataset_torquebased_steady_*')
load(steady_name.name)

all_name = dir('WBdataset_all_*')
load(all_name.name)

% loadname=dir('kinflightpathDB_pos_qbodyEKF_INCroll_9clusters*')
loadname=dir('kinflightpathDB_pos_qbodyEKF_INCroll*')
% if exist(loadname) ~= 2
% %     loadname=('flightpathDB_pos_qbodyEKF_NOroll_9clusters_2.5n-3.3n3_startframe2945.mat')
%     loadname=('flightpathDB_pos_qbodyEKF_INCroll_9clusters_2.5n-3.3n3_startframe2945_strokeplane45deg.mat')
%     loadname=('flightpathDB_pos_qbodyEKF_INCroll_9clusters_2.5n-3.3n3_startframe2945_strokeplane55deg.mat')
% end
load(loadname.name)

mkdir('response_figs')
cd('response_figs')

plot_timelines = 0
plot_timelines = 1

cluster_on = 0
% cluster_on = 1

%% settings
linewidth_timelines = 1.5;
skip = 50;

% heatmap resolution
nx = 1000;
ny = 100;

cmap_180 = jet(180);
cmap_bw = 1-colormap(gray);

% polyfit & 95% cof int settings
order = 3;
Nbin = 51;
Mbin = 3;
dn=20   % datapoints in bin
dm=20   % bin shift
csaps_filt = .00001;

%% data import
strokeplane_angle = settings.strokeplane_angle;
IDX = pathDB.IDX;
cmap_k = settings.cmap_k;
cmap_360 = settings.cmap_360;

t = pathDB.t;
V = pathDB.V;
dt = t(2)-t(1);

t_all = t;
for i=1:size(IDX,2)-1
    t_all = [t_all t];
end

stim_angle_vel = pathDB.stim_angle_vel;
stim_angle_accel = pathDB.stim_angle_accel;
stim_angle_F = pathDB.stim_angle_F;
stim_angle_yaw = pathDB.stim_angle_yaw;
stim_angle_spn = pathDB.stim_angle_spn;

accel_angle_hor_vel = pathDB.accel_angle_hor_vel;
accel_angle_hor_body = pathDB.accel_angle_hor_body;

F_angle_hor_vel = pathDB.F_angle_hor_vel;
F_angle_hor_body = pathDB.F_angle_hor_body;

slip = -pathDB.slip;

slip_global = -pathDB.slip;
pitch_global = pathDB.pitch;
roll_global = pathDB.roll;

dyaw = pathDB.dyaw_sp;
dpitch = pathDB.dpitch_sp;
droll = pathDB.droll_sp;

yaw_dot = pathDB.yaw_dot_sp;
pitch_dot = pathDB.pitch_dot_sp;
roll_dot = pathDB.roll_dot_sp;

dyaw_body = pathDB.dyaw_body;
dpitch_body = pathDB.dpitch_body;
droll_body = pathDB.droll_body;

yaw_dot_body = pathDB.yaw_dot_body;
pitch_dot_body = pathDB.pitch_dot_body;
roll_dot_body = pathDB.roll_dot_body;

slip_body = pathDB.slip_body;
pitch_body = pathDB.pitch_body;
roll_body = pathDB.roll;

slip_aero = pathDB.slip_aero;
pitch_aero = pathDB.pitch_aero;
slip_aero_hor = pathDB.slip_aero_hor;
pitch_aero_hor = pathDB.pitch_aero_hor;

A = pathDB.A;
An = pathDB.An;
At = pathDB.At;

A_hor = pathDB.A_hor;
An_hor = pathDB.An_hor;
At_hor = pathDB.At_hor;

F = pathDB.F;
F_ver = pathDB.Fz;
F_hor = pathDB.F_hor;
Fn_hor = pathDB.Fn_hor;
Ft_hor = pathDB.Ft_hor;

Fsp_pitch = pathDB.Fsp_pitch;
Fsp_roll = pathDB.Fsp_roll;

Fb_pitch = pathDB.Fb_pitch;
Fb_roll = pathDB.Fb_roll;

teta = patternDB.teta;

n_first = responseDB.n_first;
n_resp = responseDB.n_resp;
n_resp_end = responseDB.n_resp_end;
n_steady_end = responseDB.n_steady_end;

n_turn_start = responseDB.n_turn_start;
n_turn_stop = responseDB.n_turn_stop;
n_turn_max = responseDB.n_turn_max;

n_accel_start = responseDB.n_accel_start;
n_accel_stop = responseDB.n_accel_stop;
n_accel_max = responseDB.n_accel_max; 

n_decel_start = responseDB.n_decel_start;
n_decel_stop = responseDB.n_decel_stop;
n_decel_min = responseDB.n_decel_min;

%% n_pre & n_post
n_pre = min([n_turn_start n_accel_start n_decel_start]')';
n_post = max([n_turn_stop n_accel_stop n_decel_stop]')'-1;
for i = 1:length(n_post)
    if isnan(n_post(i)) == 1
        V_temp = V(:,i);
        n_post(i,1) = find(isnan(V_temp)==0, 1, 'last' );
    end
end

%% angular accels
yaw_dot_dot = nan(size(yaw_dot));
pitch_dot_dot = nan(size(pitch_dot));
roll_dot_dot = nan(size(roll_dot));

yaw_dot_dot(2:end-1,:) = yaw_dot(3:end,:) - yaw_dot(1:end-2,:);
pitch_dot_dot(2:end-1,:) = pitch_dot(3:end,:) - pitch_dot(1:end-2,:);
roll_dot_dot(2:end-1,:) = roll_dot(3:end,:) - roll_dot(1:end-2,:);

yaw_dot_dot = yaw_dot_dot/2/dt;
pitch_dot_dot = pitch_dot_dot/2/dt;
roll_dot_dot = roll_dot_dot/2/dt;

%% WBmeans
freq_steady = f_wb_steady_meanCIstd(1);
NsteadyWB = round(1/freq_steady/dt);

% yaw_dot_WBmean = nan(size(yaw_dot));
% pitch_dot_WBmean = nan(size(pitch_dot));
% roll_dot_WBmean = nan(size(roll_dot));
% F_WBmean = nan(size(F));
% for i = NsteadyWB/2+1:size(yaw_dot,1)-NsteadyWB/2
%     yaw_dot_WBmean(i,:) = nanmean(yaw_dot(i-NsteadyWB/2:i+NsteadyWB/2,:));
%     pitch_dot_WBmean(i,:) = nanmean(pitch_dot(i-NsteadyWB/2:i+NsteadyWB/2,:));
%     roll_dot_WBmean(i,:) = nanmean(roll_dot(i-NsteadyWB/2:i+NsteadyWB/2,:));
% end
%     
% yaw_dot_dot_WBmean = nan(size(yaw_dot_dot));
% pitch_dot_dot_WBmean = nan(size(pitch_dot_dot));
% roll_dot_dot_WBmean = nan(size(roll_dot_dot));
% F_WBmean = nan(size(F));
% for i = NsteadyWB/2+1:size(yaw_dot_dot,1)-NsteadyWB/2
%     yaw_dot_dot_WBmean(i,:) = nanmean(yaw_dot_dot(i-NsteadyWB/2:i+NsteadyWB/2,:));
%     pitch_dot_dot_WBmean(i,:) = nanmean(pitch_dot_dot(i-NsteadyWB/2:i+NsteadyWB/2,:));
%     roll_dot_dot_WBmean(i,:) = nanmean(roll_dot_dot(i-NsteadyWB/2:i+NsteadyWB/2,:));
%     F_WBmean(i,:) = nanmean(F(i-NsteadyWB/2:i+NsteadyWB/2,:));
% end

%% WBmod values
% yaw_dot_dot_norm = yaw_dot_dot / yawaccel_norm;
% pitch_dot_dot_norm = pitch_dot_dot / pitchaccel_norm;
% roll_dot_dot_norm = roll_dot_dot / rollaccel_norm;
% F_norm = (F-1) / Fenhance_norm;
% 
% yaw_dot_dot_WBmean_norm = yaw_dot_dot_WBmean / yawaccel_norm;
% pitch_dot_dot_WBmean_norm = pitch_dot_dot_WBmean / pitchaccel_norm;
% roll_dot_dot_WBmean_norm = roll_dot_dot_WBmean / rollaccel_norm;
% F_WBmean_norm = (F_WBmean-1) / Fenhance_norm;


%% dV & attitudes from omega INC initial conditions
roll = nan(size(droll));
yaw = nan(size(dyaw));
pitch = nan(size(dpitch));
dV =  nan(size(V));

% initial conditions
V_trig2resp = calc_mean_value(V,n_first,n_resp);
roll_global_first = calc_value(roll_global,n_first);
slip_global_first = calc_value(slip_global,n_first);
pitch_global_first = calc_value(pitch_global,n_first);

pitch_global_first_minsp = pitch_global_first + strokeplane_angle;

% clear dV
for i = 1:length(V_trig2resp)
    dV(:,i) = V(:,i) - V_trig2resp(i);
    if roll_global_first(i) < 135 && roll_global_first(i) > -135
        roll(:,i) = roll_global_first(i) + droll(:,i);
    else
        roll(:,i) = droll(:,i);
    end
%     yaw(:,i) = slip_global_first(i) + dyaw(:,i);
    yaw(:,i) = dyaw(:,i);
%     pitch(:,i) = pitch_global_first(i) + dpitch(:,i);
    pitch(:,i) = pitch_global_first_minsp(i) + dpitch(:,i);
end


%% MIRROR TURN based on An_hor_max
An_hor_max_mirror = calc_value(An_hor,n_turn_max);

stim_angle_vel_mirror = stim_angle_vel;
stim_angle_accel_mirror = stim_angle_accel;
stim_angle_yaw_mirror = stim_angle_yaw;
stim_angle_F_mirror = stim_angle_F;
stim_angle_spn_mirror = stim_angle_spn;
An_hor_mirror = An_hor;
Fn_hor_mirror = Fn_hor;

for i=1:length(An_hor_max_mirror)
    if An_hor_max_mirror(i) < 0
        stim_angle_vel_mirror(:,i) = -stim_angle_vel_mirror(:,i);
        stim_angle_accel_mirror(:,i) = -stim_angle_accel_mirror(:,i);
        stim_angle_yaw_mirror(:,i) = -stim_angle_yaw_mirror(:,i);
        stim_angle_F_mirror(:,i) = -stim_angle_F_mirror(:,i);
        stim_angle_spn_mirror(:,i) = -stim_angle_spn_mirror(:,i);
        An_hor_mirror(:,i) = -An_hor_mirror(:,i);
        Fn_hor_mirror(:,i) = -Fn_hor_mirror(:,i);
    end
end

%% MIRROR ATTITUDES based on roll_dot_dot_max
roll_dot_pre_mirror = calc_value(roll_dot,n_pre);
roll_dot_dot_pre_mirror = calc_value(roll_dot_dot,n_pre);

slip_global_mirror = slip_global;
roll_global_mirror = roll_global;
yaw_mirror = yaw;
roll_mirror = roll;
dyaw_mirror = dyaw;
droll_mirror = droll;
yaw_dot_mirror = yaw_dot;
roll_dot_mirror = roll_dot;
yaw_dot_dot_mirror = yaw_dot_dot;
roll_dot_dot_mirror = roll_dot_dot;
% yaw_dot_WBmean_mirror = yaw_dot_WBmean;
% roll_dot_WBmean_mirror = roll_dot_WBmean;
% yaw_dot_dot_WBmean_mirror = yaw_dot_dot_WBmean;
% roll_dot_dot_WBmean_mirror = roll_dot_dot_WBmean;
% yaw_dot_dot_norm_mirror = yaw_dot_dot_norm;
% roll_dot_dot_norm_mirror = roll_dot_dot_norm;
% yaw_dot_dot_WBmean_norm_mirror = yaw_dot_dot_WBmean_norm;
% roll_dot_dot_WBmean_norm_mirror = roll_dot_dot_WBmean_norm;
Fsp_roll_mirror = Fsp_roll;

for i=1:length(roll_dot_pre_mirror)
    if roll_dot_pre_mirror(i) < 0
% for i=1:length(roll_dot_dot_pre_mirror)
%     if roll_dot_dot_pre_mirror(i) < 0
        slip_global_mirror(:,i) = -slip_global_mirror(:,i);
        roll_global_mirror(:,i) = -roll_global_mirror(:,i);
        yaw_mirror(:,i) = -yaw_mirror(:,i);
        roll_mirror(:,i) = -roll_mirror(:,i);
        dyaw_mirror(:,i) = -dyaw_mirror(:,i);
        droll_mirror(:,i) = -droll_mirror(:,i);
        yaw_dot_mirror(:,i) = -yaw_dot_mirror(:,i);
        roll_dot_mirror(:,i) = -roll_dot_mirror(:,i);
        yaw_dot_dot_mirror(:,i) = -yaw_dot_dot_mirror(:,i);
        roll_dot_dot_mirror(:,i) = -roll_dot_dot_mirror(:,i);
%         yaw_dot_WBmean_mirror(:,i) = -yaw_dot_WBmean_mirror(:,i);
%         roll_dot_WBmean_mirror(:,i) = -roll_dot_WBmean_mirror(:,i);
%         yaw_dot_dot_WBmean_mirror(:,i) = -yaw_dot_dot_WBmean_mirror(:,i);
%         roll_dot_dot_WBmean_mirror(:,i) = -roll_dot_dot_WBmean_mirror(:,i);
%         yaw_dot_dot_norm_mirror(:,i) = -yaw_dot_dot_norm_mirror(:,i);
%         roll_dot_dot_norm_mirror(:,i) = -roll_dot_dot_norm_mirror(:,i);
%         yaw_dot_dot_WBmean_norm_mirror(:,i) = -yaw_dot_dot_WBmean_norm_mirror(:,i);
%         roll_dot_dot_WBmean_norm_mirror(:,i) = -roll_dot_dot_WBmean_norm_mirror(:,i);
        Fsp_roll_mirror(:,i) = -Fsp_roll_mirror(:,i);
    end
end

%% Torques (based on mean flap freq)
l_wing = kinDB.wing_length;
f_mean = mean(f_wb_L);

Mroll_accel  = nan(size(roll_dot_dot));
Mpitch_accel = nan(size(roll_dot_dot));
Myaw_accel   = nan(size(roll_dot_dot));

Mroll_damp  = nan(size(roll_dot_dot));
Mpitch_damp = nan(size(roll_dot_dot));
Myaw_damp   = nan(size(roll_dot_dot));

for i=1:length(l_wing)
    Mroll_accel(:,i)  = l_wing(i).^5 .* f_mean.^2 .* Iroll  .* roll_dot_dot_mirror(:,i);
    Mpitch_accel(:,i) = l_wing(i).^5 .* f_mean.^2 .* Ipitch .* pitch_dot_dot(:,i);
    Myaw_accel(:,i)   = l_wing(i).^5 .* f_mean.^2 .* Iyaw   .* yaw_dot_dot_mirror(:,i);

    Mroll_damp(:,i)  = l_wing(i).^5 .* f_mean.^2 .* Croll  .* roll_dot_mirror(:,i);
    Mpitch_damp(:,i) = l_wing(i).^5 .* f_mean.^2 .* Cpitch .* pitch_dot(:,i);
    Myaw_damp(:,i)   = l_wing(i).^5 .* f_mean.^2 .* Cyaw   .* yaw_dot_mirror(:,i);
end

Mroll  = Mroll_accel  + Mroll_damp;
Mpitch = Mpitch_accel + Mpitch_damp;
Myaw   = Myaw_accel   + Myaw_damp;

%% WBmod values
Myaw_Norm = Myaw / Myaw_norm;
Mpitch_Norm = Mpitch / Mpitch_norm;
Mroll_Norm = Mroll / Mroll_norm;
Fenhance_Norm = (F-1) / Fenhance_norm;

Myaw_WBmean_Norm = Myaw_mean_wb / Myaw_norm;
Mpitch_WBmean_Norm = Mpitch_mean_wb / Mpitch_norm;
Mroll_WBmean_Norm = Mroll_mean_wb / Mroll_norm;
F_WBmean_norm = (F_mean_wb-1) / Fenhance_norm;

%% WBmod kinematics
% RollAccel
loc = find(abs(DstrokeMOD_wb_RollAccel_bins_meanCIstd(:,1)) == max(abs(DstrokeMOD_wb_RollAccel_bins_meanCIstd(:,1))));
DstrokeMOD_RollAccel_max = sign(stroke_wb_steady_bins_meanCIstd(loc,1)) * DstrokeMOD_wb_RollAccel_bins_meanCIstd(loc,1);

loc = find(abs(DpitchMOD_wb_RollAccel_bins_meanCIstd(:,1)) == max(abs(DpitchMOD_wb_RollAccel_bins_meanCIstd(:,1))));
DpitchMOD_RollAccel_max = sign(pitch_wb_steady_bins_meanCIstd(loc,1)) * DpitchMOD_wb_RollAccel_bins_meanCIstd(loc,1);

loc = find(abs(DdevMOD_wb_RollAccel_bins_meanCIstd(:,1)) == max(abs(DdevMOD_wb_RollAccel_bins_meanCIstd(:,1))));
DdevMOD_RollAccel_max = sign(dev_wb_steady_bins_meanCIstd(loc,1)) * DdevMOD_wb_RollAccel_bins_meanCIstd(loc,1);

% YawAccel
loc = find(abs(DstrokeMOD_wb_YawAccel_bins_meanCIstd(:,1)) == max(abs(DstrokeMOD_wb_YawAccel_bins_meanCIstd(:,1))));
DstrokeMOD_YawAccel_max = sign(stroke_wb_steady_bins_meanCIstd(loc,1)) * DstrokeMOD_wb_YawAccel_bins_meanCIstd(loc,1);

loc = find(abs(DpitchMOD_wb_YawAccel_bins_meanCIstd(:,1)) == max(abs(DpitchMOD_wb_YawAccel_bins_meanCIstd(:,1))));
DpitchMOD_YawAccel_max = sign(pitch_wb_steady_bins_meanCIstd(loc,1)) * DpitchMOD_wb_YawAccel_bins_meanCIstd(loc,1);

loc = find(abs(DdevMOD_wb_YawAccel_bins_meanCIstd(:,1)) == max(abs(DdevMOD_wb_YawAccel_bins_meanCIstd(:,1))));
DdevMOD_YawAccel_max = sign(dev_wb_steady_bins_meanCIstd(loc,1)) * DdevMOD_wb_YawAccel_bins_meanCIstd(loc,1);

% PitchAccel
loc = find(abs(strokeMOD_wb_PitchAccel_bins_meanCIstd(:,1)) == max(abs(strokeMOD_wb_PitchAccel_bins_meanCIstd(:,1))));
strokeMOD_PitchAccel_max = sign(stroke_wb_steady_bins_meanCIstd(loc,1)) * strokeMOD_wb_PitchAccel_bins_meanCIstd(loc,1);

loc = find(abs(pitchMOD_wb_PitchAccel_bins_meanCIstd(:,1)) == max(abs(pitchMOD_wb_PitchAccel_bins_meanCIstd(:,1))));
pitchMOD_PitchAccel_max = sign(pitch_wb_steady_bins_meanCIstd(loc,1)) * pitchMOD_wb_PitchAccel_bins_meanCIstd(loc,1);

loc = find(abs(devMOD_wb_PitchAccel_bins_meanCIstd(:,1)) == max(abs(devMOD_wb_PitchAccel_bins_meanCIstd(:,1))));
devMOD_PitchAccel_max = sign(dev_wb_steady_bins_meanCIstd(loc,1)) * devMOD_wb_PitchAccel_bins_meanCIstd(loc,1);

% Fenhance
loc = find(abs(strokeMOD_wb_Fenhance_bins_meanCIstd(:,1)) == max(abs(strokeMOD_wb_Fenhance_bins_meanCIstd(:,1))));
strokeMOD_Fenhance_max = sign(stroke_wb_steady_bins_meanCIstd(loc,1)) * strokeMOD_wb_Fenhance_bins_meanCIstd(loc,1);

loc = find(abs(pitchMOD_wb_Fenhance_bins_meanCIstd(:,1)) == max(abs(pitchMOD_wb_Fenhance_bins_meanCIstd(:,1))));
pitchMOD_Fenhance_max = sign(pitch_wb_steady_bins_meanCIstd(loc,1)) * pitchMOD_wb_Fenhance_bins_meanCIstd(loc,1);

loc = find(abs(devMOD_wb_Fenhance_bins_meanCIstd(:,1)) == max(abs(devMOD_wb_Fenhance_bins_meanCIstd(:,1))));
devMOD_Fenhance_max = sign(dev_wb_steady_bins_meanCIstd(loc,1)) * devMOD_wb_Fenhance_bins_meanCIstd(loc,1);

freqMOD_Fenhance = freqMOD_wb_Fenhance_meanCIstd(1);

%% calculate variables pre post max etc
% calc_var_pre_post_etc_N

%% calculate nMAX
% calc_nMAX
calc_nMAX_torque

%% maneuver timing
% accel based timing
dt_RollAccel_max_min = (n_rolldot_dotmirrormin - n_rolldot_dotmirrormax)*dt;
dt_YawAccel_max_min = (n_yawdot_dotmirrormin - n_yawdot_dotmirrormax)*dt;
dt_PitchAccel_max_min = (n_pitchdot_dotmin - n_pitchdot_dotmax)*dt;

NwbSteady_RollAccel_max_min = dt_RollAccel_max_min*freq_steady;
NwbSteady_YawAccel_max_min = dt_YawAccel_max_min*freq_steady;
NwbSteady_PitchAccel_max_min = dt_PitchAccel_max_min*freq_steady;

% dt_RollAccel_WBmean_max_min = (n_rolldot_dot_WBmeanmirrormin - n_rolldot_dot_WBmeanmirrormax)*dt;
% dt_YawAccel_WBmean_max_min = (n_yawdot_dot_WBmeanmirrormin - n_yawdot_dot_WBmeanmirrormax)*dt;
% dt_PitchAccel_WBmean_max_min = (n_pitchdot_dot_WBmeanmin - n_pitchdot_dot_WBmeanmax)*dt;
% 
% NwbSteady_RollAccel_WBmean_max_min = dt_RollAccel_WBmean_max_min*freq_steady;
% NwbSteady_YawAccel_WBmean_max_min = dt_YawAccel_WBmean_max_min*freq_steady;
% NwbSteady_PitchAccel_WBmean_max_min = dt_PitchAccel_WBmean_max_min*freq_steady;

% velocity based timing
dt_RollVel_max_min = (n_rolldotmirrormin - n_rolldotmirrormax)*dt;
dt_YawVel_max_min = (n_yawdotmirrormin - n_yawdotmirrormax)*dt;
dt_PitchVel_max_min = (n_pitchdotmin - n_pitchdotmax)*dt;

NwbSteady_RollVel_max_min = dt_RollVel_max_min*freq_steady;
NwbSteady_YawVel_max_min = dt_YawVel_max_min*freq_steady;
NwbSteady_PitchVel_max_min = dt_PitchVel_max_min*freq_steady;

% dt_RollVel_WBmean_max_min = (n_rolldot_WBmeanmirrormin - n_rolldot_WBmeanmirrormax)*dt;
% dt_YawVel_WBmean_max_min = (n_yawdot_WBmeanmirrormin - n_yawdot_WBmeanmirrormax)*dt;
% dt_PitchVel_WBmean_max_min = (n_pitchdot_WBmeanmin - n_pitchdot_WBmeanmax)*dt;
% 
% NwbSteady_RollVel_WBmean_max_min = dt_RollVel_WBmean_max_min*freq_steady;
% NwbSteady_YawVel_WBmean_max_min = dt_YawVel_WBmean_max_min*freq_steady;
% NwbSteady_PitchVel_WBmean_max_min = dt_PitchVel_WBmean_max_min*freq_steady;

% TORQUE based timing
dt_Mroll_max_min = (n_Mroll_min - n_Mroll_max)*dt;
dt_Myaw_max_min = (n_Myaw_min - n_Myaw_max)*dt;
dt_Mpitch_max_min = (n_Mpitch_min - n_Mpitch_max)*dt;

NwbSteady_Roll_max_min = dt_Mroll_max_min*freq_steady;
NwbSteady_Yaw_max_min = dt_Myaw_max_min*freq_steady;
NwbSteady_Pitch_max_min = dt_Mpitch_max_min*freq_steady;

% dt_Roll_WBmean_max_min = (n_rolldot_WBmeanmirrormin - n_rolldot_WBmeanmirrormax)*dt;
% dt_Yaw_WBmean_max_min = (n_yawdot_WBmeanmirrormin - n_yawdot_WBmeanmirrormax)*dt;
% dt_Pitch_WBmean_max_min = (n_pitchdot_WBmeanmin - n_pitchdot_WBmeanmax)*dt;
% 
% NwbSteady_Roll_WBmean_max_min = dt_Roll_WBmean_max_min*freq_steady;
% NwbSteady_Yaw_WBmean_max_min = dt_Yaw_WBmean_max_min*freq_steady;
% NwbSteady_Pitch_WBmean_max_min = dt_Pitch_WBmean_max_min*freq_steady;

%% cal WBmod kin data
% cal_WBmod_kin_data
cal_WBmod_kin_data_torque

%% calc yaw & vel turns !!@!! FIX by using heading pre&post
% calc_turn_vectors
% calc_turn_vectors_headingyaw_pre_post
calc_heading_turn_pre_post_mirror
% calc_heading_turn_pre_post


%% TIMELINES
if plot_timelines == 1
%% timelines NO skip 

IDX_plot = IDX;
t_plot = t;

V_plot = V;
dV_plot = dV;

An_hor_plot = An_hor;
At_hor_plot = At_hor;
A_hor_plot = A_hor;

F_plot = F;
F_ver_plot = F_ver;
F_hor_plot = F_hor;
Fn_hor_plot = Fn_hor;
Ft_hor_plot = Ft_hor;

stim_angle_vel_plot = stim_angle_vel;
stim_angle_accel_plot = stim_angle_accel;
stim_angle_yaw_plot = stim_angle_yaw;
stim_angle_F_plot = stim_angle_F;
stim_angle_spn_plot = stim_angle_spn;

yaw_plot = yaw;
pitch_plot = pitch;
roll_plot = roll;

dyaw_plot = dyaw;
dpitch_plot = dpitch;
droll_plot = droll;

yaw_dot_plot = yaw_dot;
pitch_dot_plot = pitch_dot;
roll_dot_plot = roll_dot;

yaw_dot_dot_plot = yaw_dot_dot;
pitch_dot_dot_plot = pitch_dot_dot;
roll_dot_dot_plot = roll_dot_dot;

% yaw_dot_dot_norm_plot = yaw_dot_dot_norm;
% pitch_dot_dot_norm_plot = pitch_dot_dot_norm;
% roll_dot_dot_norm_plot = roll_dot_dot_norm;

slip_global_plot = slip_global;
pitch_global_plot = pitch_global;
roll_global_plot = roll_global;

Fsp_pitch_plot = Fsp_pitch;
Fsp_roll_plot = Fsp_roll;

F_plot = F;
F_norm_plot = Fenhance_Norm;

stim_angle_vel_temp = stim_angle_vel;
stim_angle_accel_temp = stim_angle_accel;
stim_angle_yaw_temp = stim_angle_yaw;
stim_angle_F_temp = stim_angle_F;
stim_angle_spn_temp = stim_angle_spn;

yaw_temp = yaw;
pitch_temp = pitch;
roll_temp = roll;

dyaw_temp = dyaw;
dpitch_temp = dpitch;
droll_temp = droll;

slip_global_temp = slip_global;
pitch_global_temp = pitch_global;
roll_global_temp = roll_global;

% remove jumps from plots (+/-180deg)
for i=1:size(stim_angle_vel_plot,2)
    for j=2:size(stim_angle_vel_plot,1)
        if abs(stim_angle_vel_temp(j,i) - stim_angle_vel_temp(j-1,i)) > 90
            stim_angle_vel_plot(j-1,i) = nan;
        end
        if abs(stim_angle_accel_temp(j,i) - stim_angle_accel_temp(j-1,i)) > 90
            stim_angle_accel_plot(j-1,i) = nan;
        end
        if abs(stim_angle_yaw_temp(j,i) - stim_angle_yaw_temp(j-1,i)) > 90
            stim_angle_yaw_plot(j-1,i) = nan;
        end
        if abs(stim_angle_F_temp(j,i) - stim_angle_F_temp(j-1,i)) > 90
            stim_angle_F_plot(j-1,i) = nan;
        end
        if abs(stim_angle_spn_temp(j,i) - stim_angle_spn_temp(j-1,i)) > 90
            stim_angle_spn_plot(j-1,i) = nan;
        end
        if abs(yaw_temp(j,i) - yaw_temp(j-1,i)) > 90
            yaw_plot(j-1,i) = nan;
        end
        if abs(pitch_temp(j,i) - pitch_temp(j-1,i)) > 90
            pitch_plot(j-1,i) = nan;
        end
        if abs(roll_temp(j,i) - roll_temp(j-1,i)) > 90
            roll_plot(j-1,i) = nan;
        end
        if abs(dyaw_temp(j,i) - dyaw_temp(j-1,i)) > 90
            dyaw_plot(j-1,i) = nan;
        end
        if abs(dpitch_temp(j,i) - dpitch_temp(j-1,i)) > 90
            dpitch_plot(j-1,i) = nan;
        end
        if abs(droll_temp(j,i) - droll_temp(j-1,i)) > 90
            droll_plot(j-1,i) = nan;
        end
        if abs(slip_global_temp(j,i) - slip_global_temp(j-1,i)) > 90
            slip_global_plot(j-1,i) = nan;
        end
        if abs(pitch_global_temp(j,i) - pitch_global_temp(j-1,i)) > 90
            pitch_global_plot(j-1,i) = nan;
        end
        if abs(roll_global_temp(j,i) - roll_global_temp(j-1,i)) > 90
            roll_global_plot(j-1,i) = nan;
        end
    end
end

t_pre = calc_t_value(t,n_pre);
t_shift = t_pre;
% t_shift = zeros(size(t_pre));
t_hist = [];
for i=1:size(V,2)    
    t_hist = [t_hist;t-t_shift(i)];
end
t_start = -.025;
t_stop = .05;

cmap_plot = cmap_360;
color_var = round(stim_angle_vel_pre)+180;

%% timelines F&SPN NO MIRROR TURN
plot_flightpath_timeline_tshift_headingstart_FnSP

saveas(gca,'flightpaths_FnSP_headingstart.fig')
saveas(gca,'flightpaths_FnSP_headingstart.png')
plot2svg('flightpaths_FnSP_headingstart.svg')

%% heatmap hist norm F&SPN NO MIRROR
plot_flightpath_timeline_histograms_norm_FnSP

saveas(gca,'flightpaths_FnSP_hist_norm_tresp.fig')
saveas(gca,'flightpaths_FnSP_hist_norm_tresp.png')
plot2svg('flightpaths_FnSP_hist_norm_tresp.svg')

%% timelines ATTITUDEnDOT NO MIRROR TURN
% plot_flightpath_timeline_tshift_headingstart_FnSPnAtt
plot_flightpath_timeline_tshift_headingstart_ATTnDOTs
% plot_flightpath_timeline_tshift_headingstart_ATTnDOTnGLOB

saveas(gca,'flightpaths_ATTnDOTs_headingstart.fig')
saveas(gca,'flightpaths_ATTnDOTs_headingstart.png')
plot2svg('flightpaths_ATTnDOTs_headingstart.svg')

%% timelines ATTITUDEnDOT NO MIRROR TURN
% plot_flightpath_timeline_tshift_headingstart_FnSPnAtt
% plot_flightpath_timeline_tshift_headingstart_ATTnDOT
plot_flightpath_timeline_tshift_headingstart_ATTnGLOB_subsets

saveas(gca,'flightpaths_ATTnGLOB_headingstart.fig')
saveas(gca,'flightpaths_ATTnGLOB_headingstart.png')
plot2svg('flightpaths_ATTnGLOB_headingstart.svg')

%% heatmap hist norm ATTITUDEnDOTs NO MIRROR
plot_flightpath_timeline_histograms_norm_ATTnDOTs

% subplot(3,2,2)
% axis([t_start t_stop -90 90])
% set(gca,'YTick',[-90;90],'fontsize',8)
% set(gca,'YTickLabel',[90;-90],'fontsize',8)
% subplot(3,2,5)
% axis([t_start t_stop -2000 2000])
% set(gca,'YTick',[-2000;0;2000],'fontsize',8)
% set(gca,'YTickLabel',[2000;0;-2000],'fontsize',8)
% subplot(3,2,6)
% axis([t_start t_stop -90 90])
% set(gca,'YTick',[-90;0;90],'fontsize',8)
% set(gca,'YTickLabel',[90;0;-90],'fontsize',8)

saveas(gca,'flightpaths_ATTnDOTs_hist_norm_tresp.fig')
saveas(gca,'flightpaths_ATTnDOTs_hist_norm_tresp.png')
plot2svg('flightpaths_ATTnDOTs_hist_norm_tresp.svg')

%% heatmap hist norm ATTITUDEnGLOB NO MIRROR
plot_flightpath_timeline_histograms_norm_ATTnGLOB
% 
% subplot(3,3,7)
% axis([t_start t_stop -2000 2000])
% set(gca,'YTick',[-2000;0;2000],'fontsize',8)
% set(gca,'YTickLabel',[2000;0;-2000],'fontsize',8)
% subplot(3,3,9)
% axis([t_start t_stop -90 90])
% set(gca,'YTick',[-90;0;90],'fontsize',8)
% set(gca,'YTickLabel',[90;0;-90],'fontsize',8)

saveas(gca,'flightpaths_ATTnGLOB_hist_norm_tresp.fig')
saveas(gca,'flightpaths_ATTnGLOB_hist_norm_tresp.png')
plot2svg('flightpaths_ATTnGLOB_hist_norm_tresp.svg')

%% timelines NO skip NO MIRROR TURN
plot_flightpath_timeline_tshift_headingstart_AnormNdir

saveas(gca,'flightpaths_ACCEL_headingstart.fig')
saveas(gca,'flightpaths_ACCEL_headingstart.png')
plot2svg('flightpaths_ACCEL_headingstart.svg')

%% heatmap hist norm NO MIRROR
plot_flightpath_timeline_histograms_norm_AnormNdir

saveas(gca,'flightpaths_ACCEL_hist_norm_tresp.fig')
saveas(gca,'flightpaths_ACCEL_hist_norm_tresp.png')
plot2svg('flightpaths_ACCEL_hist_norm_tresp.svg')

%% MIRROR TURN
for i=1:size(stim_angle_vel_plot,2)
    if An_hor_max_mirror(i) < 0
        stim_angle_vel_plot(:,i) = -stim_angle_vel_plot(:,i);
        stim_angle_accel_plot(:,i) = -stim_angle_accel_plot(:,i);
        stim_angle_yaw_plot(:,i) = -stim_angle_yaw_plot(:,i);
        stim_angle_F_plot(:,i) = -stim_angle_F_plot(:,i);
        stim_angle_spn_plot(:,i) = -stim_angle_spn_plot(:,i);
        An_hor_plot(:,i) = -An_hor_plot(:,i);
        Fn_hor_plot(:,i) = -Fn_hor_plot(:,i);
    end
    if roll_dot_pre_mirror(i) < 0
%     if roll_dot_dot_pre_mirror(i) < 0
        yaw_plot(:,i) = -yaw_plot(:,i);
        roll_plot(:,i) = -roll_plot(:,i);
        dyaw_plot(:,i) = -dyaw_plot(:,i);
        droll_plot(:,i) = -droll_plot(:,i);
        yaw_dot_plot(:,i) = -yaw_dot_plot(:,i);
        roll_dot_plot(:,i) = -roll_dot_plot(:,i);
        yaw_dot_dot_plot(:,i) = -yaw_dot_dot_plot(:,i);
        roll_dot_dot_plot(:,i) = -roll_dot_dot_plot(:,i);
%         yaw_dot_dot_norm_plot(:,i) = -yaw_dot_dot_norm_plot(:,i);
%         roll_dot_dot_norm_plot(:,i) = -roll_dot_dot_norm_plot(:,i);
        Fsp_roll_plot(:,i) = -Fsp_roll_plot(:,i);
        slip_global_plot(:,i) = -slip_global_plot(:,i);
        roll_global_plot(:,i) = -roll_global_plot(:,i);
    end
end

cmap_plot = cmap_180;
color_var = round(abs(stim_angle_vel_mirror_pre));

%% timelines F&SPN INC MIRROR TURN
plot_flightpath_timeline_tshift_headingstart_FnSP

saveas(gca,'flightpaths_FnSP_headingstart_mirror.fig')
saveas(gca,'flightpaths_FnSP_headingstart_mirror.png')
plot2svg('flightpaths_FnSP_headingstart_mirror.svg')

%% heatmap hist norm F&SPN INC MIRROR TURN
plot_flightpath_timeline_histograms_norm_FnSP

saveas(gca,'flightpaths_FnSP_hist_norm_tresp_mirror.fig')
saveas(gca,'flightpaths_FnSP_hist_norm_tresp_mirror.png')
plot2svg('flightpaths_FnSP_hist_norm_tresp_mirror.svg')

%% timelines WB mod values
plot_flightpath_timeline_tshift_headingstart_WBmodValues

saveas(gca,'flightpaths_WBmodValues_headingstart_mirror.fig')
saveas(gca,'flightpaths_WBmodValues_headingstart_mirror.png')
plot2svg('flightpaths_WBmodValues_headingstart_mirror.svg')

%% timelines ATTITUDEnDOTs INC MIRROR TURN
plot_flightpath_timeline_tshift_headingstart_ATTnDOTs

% subplot(3,2,2)
% axis([t_start t_stop -45 135])
% subplot(3,2,5)
% axis([t_start t_stop -1000 3000])
% set(gca,'YTick',[-1000;0;3000],'fontsize',8)
% subplot(3,2,6)
% axis([t_start t_stop -45 135])

saveas(gca,'flightpaths_ATTnDOTs_headingstart_mirror.fig')
saveas(gca,'flightpaths_ATTnDOTs_headingstart_mirror.png')
plot2svg('flightpaths_ATTnDOTs_headingstart_mirror.svg')

%% timelines ATTITUDEnGLOB INC MIRROR TURN
plot_flightpath_timeline_tshift_headingstart_ATTnGLOB

% subplot(3,3,2)
% axis([t_start t_stop -45 135])
% subplot(3,3,3)
% axis([t_start t_stop -45 135])
% subplot(3,3,7)
% axis([t_start t_stop -1000 3000])
% set(gca,'YTick',[-1000;0;3000],'fontsize',8)
% subplot(3,3,8)
% axis([t_start t_stop -90 90])
% subplot(3,3,9)
% axis([t_start t_stop -90 90])

saveas(gca,'flightpaths_AttSPnEULER_headingstart_mirror.fig')
saveas(gca,'flightpaths_AttSPnEULER_headingstart_mirror.png')
plot2svg('flightpaths_AttSPnEULER_headingstart_mirror.svg')

%% heatmap hist norm ATTITUDEnDOT INC MIRROR TURN
plot_flightpath_timeline_histograms_norm_ATTnDOTs

% subplot(3,2,2)
% axis([t_start t_stop -135 45])
% set(gca,'YTick',[-90;0],'fontsize',8)
% set(gca,'YTickLabel',[90;0],'fontsize',8)
% subplot(3,2,5)
% axis([t_start t_stop -3000 1000])
% set(gca,'YTick',[-3000;0;1000],'fontsize',8)
% set(gca,'YTickLabel',[1000;0;-3000],'fontsize',8)
% subplot(3,2,6)
% axis([t_start t_stop -135 45])
% set(gca,'YTick',[-90;0],'fontsize',8)
% set(gca,'YTickLabel',[90;0],'fontsize',8)

saveas(gca,'flightpaths_ATTnDOTs_hist_norm_tresp_mirror.fig')
saveas(gca,'flightpaths_ATTnDOTs_hist_norm_tresp_mirror.png')
plot2svg('flightpaths_ATTnDOTs_hist_norm_tresp_mirror.svg')

%% heatmap hist norm ATTITUDEnGLOB INC MIRROR TURN
plot_flightpath_timeline_histograms_norm_ATTnGLOB

% subplot(3,3,7)
% axis([t_start t_stop -3000 1000])
% set(gca,'YTick',[-3000;0;1000],'fontsize',8)
% set(gca,'YTickLabel',[1000;0;-3000],'fontsize',8)
% subplot(3,3,2)
% axis([t_start t_stop -135 45])
% set(gca,'YTick',[-90;0],'fontsize',8)
% set(gca,'YTickLabel',[90;0],'fontsize',8)
% subplot(3,3,3)
% axis([t_start t_stop -135 45])
% set(gca,'YTick',[-90;0],'fontsize',8)
% set(gca,'YTickLabel',[90;0],'fontsize',8)
% subplot(3,3,5)
% axis([t_start t_stop -135 45])
% set(gca,'YTick',[-90;0],'fontsize',8)
% set(gca,'YTickLabel',[90;0],'fontsize',8)
% subplot(3,3,6)
% axis([t_start t_stop -135 45])
% set(gca,'YTick',[-90;0],'fontsize',8)
% set(gca,'YTickLabel',[90;0],'fontsize',8)
% subplot(3,3,8)
% axis([t_start t_stop -90 90])
% set(gca,'YTick',[-90;0;90],'fontsize',8)
% set(gca,'YTickLabel',[90;0;-90],'fontsize',8)
% subplot(3,3,9)
% axis([t_start t_stop -90 90])
% set(gca,'YTick',[-90;0;90],'fontsize',8)
% set(gca,'YTickLabel',[90;0;-90],'fontsize',8)
% 
saveas(gca,'flightpaths_ATTnGLOB_hist_norm_tresp_mirror.fig')
saveas(gca,'flightpaths_ATTnGLOB_hist_norm_tresp_mirror.png')
plot2svg('flightpaths_ATTnGLOB_hist_norm_tresp_mirror.svg')

%% timelines NO skip INC MIRROR TURN
plot_flightpath_timeline_tshift_headingstart_AnormNdir
% 
% subplot(3,2,4)
% axis([t_start t_stop -10 20])
% set(gca,'YTick',[-10;0;20],'fontsize',8)
% subplot(3,2,6)
% axis([t_start t_stop -10 20])
% set(gca,'YTick',[-10;0;20],'fontsize',8)

saveas(gca,'flightpaths_ACCEL_headingstart_mirror.fig')
saveas(gca,'flightpaths_ACCEL_headingstart_mirror.png')
plot2svg('flightpaths_ACCEL_headingstart_mirror.svg')

%% heatmap hist norm INC MIRROR TURN
plot_flightpath_timeline_histograms_norm_AnormNdir

subplot(3,2,4)
axis([t_start t_stop -20 10])
set(gca,'YTick',[-20;0;10],'YTicklabel',[20;0;-10],'fontsize',8) 
subplot(3,2,6)
axis([t_start t_stop -20 10])
set(gca,'YTick',[-20;0;10],'YTicklabel',[20;0;-10],'fontsize',8) 

saveas(gca,'flightpaths_ACCEL_hist_norm_tresp_mirror.fig')
saveas(gca,'flightpaths_ACCEL_hist_norm_tresp_mirror.png')
plot2svg('flightpaths_ACCEL_hist_norm_tresp_mirror.svg')

%% cluster plots
if cluster_on == 1
%% timelines INC skip

IDX_plot = IDX(1:skip:end,:);
t_plot = t(1:skip:end,:);

V_plot = V(1:skip:end,:);
dV_plot = dV(1:skip:end,:);
An_hor_plot = An_hor(1:skip:end,:);
At_hor_plot = At_hor(1:skip:end,:);
A_hor_plot = A_hor(1:skip:end,:);

stim_angle_vel_plot = stim_angle_vel(1:skip:end,:);
stim_angle_accel_plot = stim_angle_accel(1:skip:end,:);
stim_angle_yaw_plot = stim_angle_yaw(1:skip:end,:);
stim_angle_F_plot = stim_angle_F(1:skip:end,:);
stim_angle_spn_plot = stim_angle_spn(1:skip:end,:);
yaw_plot = yaw(1:skip:end,:);
pitch_plot = pitch(1:skip:end,:);
roll_plot = roll(1:skip:end,:);

stim_angle_vel_temp = stim_angle_vel(1:skip:end,:);
stim_angle_accel_temp = stim_angle_accel(1:skip:end,:);
stim_angle_yaw_temp = stim_angle_yaw(1:skip:end,:);
stim_angle_F_temp = stim_angle_F(1:skip:end,:);
stim_angle_spn_temp = stim_angle_spn(1:skip:end,:);
yaw_temp = yaw(1:skip:end,:);
pitch_temp = pitch(1:skip:end,:);
roll_temp = roll(1:skip:end,:);

% remove jumps from plots (+/-180deg)
for i=1:size(stim_angle_vel_plot,2)
    for j=2:size(stim_angle_vel_plot,1)
        if abs(stim_angle_vel_temp(j,i) - stim_angle_vel_temp(j-1,i)) > 90
            stim_angle_vel_plot(j-1,i) = nan;
        end
        if abs(stim_angle_accel_temp(j,i) - stim_angle_accel_temp(j-1,i)) > 90
            stim_angle_accel_plot(j-1,i) = nan;
        end
        if abs(stim_angle_yaw_temp(j,i) - stim_angle_yaw_temp(j-1,i)) > 90
            stim_angle_yaw_plot(j-1,i) = nan;
        end
        if abs(stim_angle_F_temp(j,i) - stim_angle_F_temp(j-1,i)) > 90
            stim_angle_F_plot(j-1,i) = nan;
        end
        if abs(stim_angle_spn_temp(j,i) - stim_angle_spn_temp(j-1,i)) > 90
            stim_angle_spn_plot(j-1,i) = nan;
        end
        if abs(yaw_temp(j,i) - yaw_temp(j-1,i)) > 90
            yaw_plot(j-1,i) = nan;
        end
        if abs(pitch_temp(j,i) - pitch_temp(j-1,i)) > 90
            pitch_plot(j-1,i) = nan;
        end
        if abs(roll_temp(j,i) - roll_temp(j-1,i)) > 90
            roll_plot(j-1,i) = nan;
        end
    end
end

%% clusters NO MIRROR
plot_flightpath_timeline_tshift_cluster_AnormNdir

saveas(gca,'flightpaths_ACCEL_clusters_tresp.fig')
saveas(gca,'flightpaths_ACCEL_clusters_tresp.png')
plot2svg('flightpaths_ACCEL_clusters_tresp.svg')


%% MIRROR TURN
for i=1:size(stim_angle_vel_plot,2)
    if An_hor_max_mirror(i) < 0
        stim_angle_vel_plot(:,i) = -stim_angle_vel_plot(:,i);
        stim_angle_accel_plot(:,i) = -stim_angle_accel_plot(:,i);
        stim_angle_yaw_plot(:,i) = -stim_angle_yaw_plot(:,i);
        stim_angle_F_plot(:,i) = -stim_angle_F_plot(:,i);
        stim_angle_spn_plot(:,i) = -stim_angle_spn_plot(:,i);
        An_hor_plot(:,i) = -An_hor_plot(:,i);
        Fn_hor_plot(:,i) = -Fn_hor_plot(:,i);
    end
    if roll_dot_pre_mirror(i) < 0
%     if roll_dot_dot_pre_mirror(i) < 0
        yaw_plot(:,i) = -yaw_plot(:,i);
        roll_plot(:,i) = -roll_plot(:,i);
        yaw_dot_plot(:,i) = -yaw_dot_plot(:,i);
        roll_dot_plot(:,i) = -roll_dot_plot(:,i);
        yaw_dot_dot_plot(:,i) = -yaw_dot_dot_plot(:,i);
        roll_dot_dot_plot(:,i) = -roll_dot_dot_plot(:,i);
        yaw_dot_dot_norm_plot(:,i) = -yaw_dot_dot_norm_plot(:,i);
        roll_dot_dot_norm_plot(:,i) = -roll_dot_dot_norm_plot(:,i);
        Fsp_roll_plot(:,i) = -Fsp_roll_plot(:,i);
    end
end

%% clusters INC MIRROR
plot_flightpath_timeline_tshift_cluster_AnormNdir

subplot(3,2,4)
axis([t_start t_stop -10 20])
set(gca,'YTick',[-10;0;20],'fontsize',8)
subplot(3,2,6)
axis([t_start t_stop -10 20])
set(gca,'YTick',[-10;0;20],'fontsize',8)

saveas(gca,'flightpaths_ACCEL_clusters_tresp_mirror.fig')
saveas(gca,'flightpaths_ACCEL_clusters_tresp_mirror.png')
plot2svg('flightpaths_ACCEL_clusters_tresp_mirror.svg')

end
end

%% At_max & At_min vs heading_pre
figure
plot(stim_angle_vel_pre,An_hor_mirror_max,'ok','markerfacecolor','y','markersize',5)
hold on
plot(stim_angle_vel_pre,At_hor_max,'ok','markerfacecolor','r','markersize',5)
plot(stim_angle_vel_pre,At_hor_min,'ok','markerfacecolor','b','markersize',5)
legend('turn','accel','decel','location','sw')
grid on
axis([0 180 -20 20])
set(gca,'XTick',0:90:180) 
set(gca,'YTick',-20:20:20,'fontsize',8)
xlabel('initial heading','fontsize',10) 
ylabel('Amax','fontsize',10) 

saveas(gca,'Ahormax_VS_heading_pre.fig')
saveas(gca,'Ahormax_VS_heading_pre.png')
plot2svg('Ahormax_VS_heading_pre.svg')

%% heatmap hist all An & At
figure

binx = 0:.5:20;
biny = -20:.5:20;

x_hist = An_hor(:);
y_hist = At_hor(:);

yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);
imagesc(binx,-biny,yx_hist_log)

axis equal
axis([0 20 -20 20])
set(gca,'XTick',0:10:20) 
set(gca,'YTick',-20:10:20,'YTicklabel',20:-10:-20,'fontsize',8)
xlabel('An','fontsize',10) 
ylabel('At','fontsize',10) 
colorbar

saveas(gca,'AnAt_hist_all.fig')
saveas(gca,'AnAt_hist_all.png')
plot2svg('AnAt_hist_all.svg')

%% heatmap hist attitudes global mirror
roll_hist = roll_global_mirror(:);
pitch_hist = pitch_global(:);
yaw_hist = slip_global_mirror(:);

binx = -180:2.5:180;
biny = -180:5:180;

figure
plot_attitudes_subplots

subplot(2,2,1)
axis([-90 90 -180 90])
xlabel('slip','fontsize',10) 
subplot(2,2,2)
axis([-45 90 -180 90])
subplot(2,2,4)
axis([-45 90 -90 90])
ylabel('slip','fontsize',10) 

saveas(gca,'attitudes_global_mirror_hist_all.fig')
saveas(gca,'attitudes_global_mirror_hist_all.png')
plot2svg('attitudes_global_mirror_hist_all.svg')


%% heatmap hist attitudes mirror
roll_hist = roll_mirror(:);
pitch_hist = pitch(:);
yaw_hist = yaw_mirror(:);

binx = -180:2.5:180;
biny = -180:2.5:180;

figure
plot_attitudes_subplots

subplot(2,2,1)
axis([-45 135 -135 45])
subplot(2,2,2)
axis([-45 135 -135 45])
subplot(2,2,4)
axis([-45 135 -135 45])

saveas(gca,'attitudes_mirror_hist_all.fig')
saveas(gca,'attitudes_mirror_hist_all.png')
plot2svg('attitudes_mirror_hist_all.svg')


%% heatmap hist attitudes DOT mirror
roll_hist = roll_dot_mirror(:);
pitch_hist = pitch_dot(:);
yaw_hist = yaw_dot_mirror(:);

binx = -4000:100:4000;
biny = -4000:100:4000;

figure
plot_attitudes_subplots

subplot(2,2,1)
% axis equal
axis([-2000 4000 -4000 4000])
xlabel('yaw rate','fontsize',10) 
ylabel('roll rate','fontsize',10) 

subplot(2,2,2)
% axis equal
axis([-2000 4000 -4000 4000])
xlabel('pitch rate','fontsize',10) 
ylabel('roll rate','fontsize',10) 

subplot(2,2,4)
% axis equal
axis([-2000 4000 -4000 2000])
xlabel('pitch rate','fontsize',10) 
ylabel('yaw rate','fontsize',10) 

saveas(gca,'attitudes_dot_mirror_hist_all.fig')
saveas(gca,'attitudes_dot_mirror_hist_all.png')
plot2svg('attitudes_dot_mirror_hist_all.svg')

%% A vs attitude_global_mirror heatmaps
An_hist = An_hor_mirror(:);
At_hist = At_hor(:);
A_hist = A_hor(:);
% An_hist = An_hor(IDX~=5);
% At_hist = At_hor(IDX~=5);
% A_hist = A_hor(IDX~=5);

roll_hist = roll_global_mirror(:);
pitch_hist = pitch_global(:);
yaw_hist = slip_global_mirror(:);

binx = -180:5:180;
biny =  -20:.5:20;

figure
normy_on=0;
plot_accel_vs_attitude_subplots

subplot(3,3,1)
axis([-90 90 -20 0])
subplot(3,3,2)
axis([-90 180 -20 0])
subplot(3,3,3)
axis([-45 90 -20 0])

subplot(3,3,4)
axis([-90 90 -20 10])
subplot(3,3,5)
axis([-90 180 -20 10])
subplot(3,3,6)
axis([-45 90 -20 10])

subplot(3,3,7)
axis([-90 90 -20 10])
xlabel('slip','fontsize',10) 
subplot(3,3,8)
axis([-90 180 -20 10])
subplot(3,3,9)
axis([-45 90 -20 10])

saveas(gca,'accels_vs_attitudes_global_mirror_hist.fig')
saveas(gca,'accels_vs_attitudes_global_mirror_hist.png')
plot2svg('accels_vs_attitudes_global_mirror_hist.svg')


%% A vs attitude mirror heatmaps
An_hist = An_hor_mirror(:);
At_hist = At_hor(:);
A_hist = A_hor(:);
% An_hist = An_hor(IDX~=5);
% At_hist = At_hor(IDX~=5);
% A_hist = A_hor(IDX~=5);

roll_hist = roll_mirror(:);
pitch_hist = pitch(:);
yaw_hist = yaw_mirror(:);

binx = -180:5:180;
biny =  -20:.5:20;

figure
normy_on=0;
plot_accel_vs_attitude_subplots

subplot(3,3,1)
axis([-45 90 -20 0])
subplot(3,3,2)
axis([-90 90 -20 0])
subplot(3,3,3)
axis([-45 135 -20 0])

subplot(3,3,4)
axis([-45 90 -20 10])
subplot(3,3,5)
axis([-90 90 -20 10])
subplot(3,3,6)
axis([-45 135 -20 10])

subplot(3,3,7)
axis([-45 90 -20 10])
subplot(3,3,8)
axis([-90 90 -20 10])
subplot(3,3,9)
axis([-45 135 -20 10])

saveas(gca,'accels_vs_attitudes_mirror_hist.fig')
saveas(gca,'accels_vs_attitudes_mirror_hist.png')
plot2svg('accels_vs_attitudes_mirror_hist.svg')

%% A vs attitude DOT mirror heatmaps
An_hist = An_hor_mirror(:);
At_hist = At_hor(:);
A_hist = A_hor(:);
% An_hist = An_hor(IDX~=5);
% At_hist = At_hor(IDX~=5);
% A_hist = A_hor(IDX~=5);

roll_hist = roll_dot_mirror(:);
pitch_hist = pitch_dot(:);
yaw_hist = yaw_dot_mirror(:);

binx = -4000:100:4000;
biny =  -20:.5:20;

figure
normy_on=0;
plot_accel_vs_attitude_subplots

subplot(3,3,1)
axis([-2000 4000 -20 0])
subplot(3,3,2)
axis([-4000 4000 -20 0])
subplot(3,3,3)
axis([-2000 4000 -20 0])

subplot(3,3,4)
axis([-2000 4000 -20 10])
subplot(3,3,5)
axis([-4000 4000 -20 10])
subplot(3,3,6)
axis([-2000 4000 -20 10])

subplot(3,3,7)
axis([-2000 4000 -20 10])
xlabel('yaw rate','fontsize',10) 
subplot(3,3,8)
axis([-4000 4000 -20 10])
xlabel('roll rate','fontsize',10) 
subplot(3,3,9)
axis([-2000 4000 -20 10])
xlabel('pitch rate','fontsize',10) 

saveas(gca,'accels_vs_attitudes_mirror_hist.fig')
saveas(gca,'accels_vs_attitudes_mirror_hist.png')
plot2svg('accels_vs_attitudes_mirror_hist.svg')

%% plot turn vs accel and decel reaction time
t_turn_start = calc_value(t_all,n_turn_start);
t_accel_start = calc_value(t_all,n_accel_start);
t_decel_start = calc_value(t_all,n_decel_start);

figure
hold on

plotcolor = 'r';
plot(t_turn_start,t_accel_start,'ok','MarkerFaceColor',plotcolor,'MarkerSize',5)

plotcolor = 'b';
plot(t_turn_start,t_decel_start,'ok','MarkerFaceColor',plotcolor,'MarkerSize',5)

plot([0,1],[0,1],'--k')
legend('accel','decel')

axis equal
grid on
xlabel('turn response time','fontsize',10) 
ylabel('acceleration response time','fontsize',10) 
set(gca,'xlim',[0 .15],'ylim',[0 .15])
set(gca,'XTick',[0:.05:.2])
set(gca,'YTick',[0:.05:.2],'fontsize',8)

saveas(gca,'response_time_tAn_vs_tAt.fig')
saveas(gca,'response_time_tAn_vs_tAt.png')
plot2svg('response_time_tAn_vs_tAt.svg')


%% attitudes vs headingPre @Amax CIRC MEAN
plot_attitudes_vs_headingPRE

% saveas(gca,'attitudesAmax_vs_headingPre.fig')
% saveas(gca,'attitudesAmax_vs_headingPre.png')
% plot2svg('attitudesAmax_vs_headingPre.svg')

%% attitudes vs headingPre MIRROR@Amax CIRC MEAN
attitudesAmax_vs_headingPre_mirror

% saveas(gca,'attitudesAmax_vs_headingPre_mirror.fig')
% saveas(gca,'attitudesAmax_vs_headingPre_mirror.png')
% plot2svg('attitudesAmax_vs_headingPre_mirror.svg')

attitudesAmax_vs_headingPre_mirror_left

saveas(gca,'attitudesAmax_vs_headingPre_mirror_LnR.fig')
saveas(gca,'attitudesAmax_vs_headingPre_mirror_LnR.png')
plot2svg('attitudesAmax_vs_headingPre_mirror_LnR.svg')

%% attitudes vs stim angle pre @Amax CIRC MEAN
plot_attitudes_vs_StimAnglePre

% saveas(gca,'attitudesAmax_vs_StimAnglePre.fig')
% saveas(gca,'attitudesAmax_vs_StimAnglePre.png')
% plot2svg('attitudesAmax_vs_StimAnglePre.svg')

%% attitudes vs stim angle pre MIRROR@Amax CIRC MEAN
attitudesAmax_vs_StimAnglePre_mirror

% saveas(gca,'attitudesAmax_vs_StimAnglePre_mirror.fig')
% saveas(gca,'attitudesAmax_vs_StimAnglePre_mirror.png')
% plot2svg('attitudesAmax_vs_StimAnglePre_mirror.svg')

attitudesAmax_vs_StimAnglePre_mirror_left

saveas(gca,'attitudesAmax_vs_StimAnglePre_mirror_LnR.fig')
saveas(gca,'attitudesAmax_vs_StimAnglePre_mirror_LnR.png')
plot2svg('attitudesAmax_vs_StimAnglePre_mirror_LnR.svg')

%% attitudes Accel vs stim angle pre MIRROR@Amax CIRC MEAN
% attitudeAccelMax_vs_StimAnglePre_mirror
attitudeAccelMaxMin_vs_StimAnglePre_mirror

% saveas(gca,'attitudesAmax_vs_StimAnglePre_mirror.fig')
% saveas(gca,'attitudesAmax_vs_StimAnglePre_mirror.png')
% plot2svg('attitudesAmax_vs_StimAnglePre_mirror.svg')

attitudeAccelMaxMin_vs_StimAnglePre_mirror_left

saveas(gca,'attitudeAccels_vs_StimAnglePre_mirror_LnR.fig')
saveas(gca,'attitudeAccels_vs_StimAnglePre_mirror_LnR.png')
plot2svg('attitudeAccels_vs_StimAnglePre_mirror_LnR.svg')

%% attitude Accels & #WBs4maneuver vs vel TURN angle MIRROR
AttAccelsNmaneuverWBs_vs_VelTurnAngle_mirror
% attitudesAmaxNaccels_vs_StimAnglePre_mirror
% attitudesAmaxNaccelsWBmean_vs_StimAnglePre_mirror

saveas(gca,'AttitudeAccelsNmaneuverWBs_vs_VelTurnAngle_mirror_LnR.fig')
saveas(gca,'AttitudeAccelsNmaneuverWBs_vs_VelTurnAngle_mirror_LnR.png')
plot2svg('AttitudeAccelsNmaneuverWBs_vs_VelTurnAngle_mirror_LnR.svg')

%% attitudes@Amax & Accel & #WBs4maneuver vs stim angle pre MIRROR@Amax CIRC MEAN
attitudesAmaxNaccelsNmaneuverWBs_vs_StimAnglePre_mirror
% attitudesAmaxNaccels_vs_StimAnglePre_mirror
% attitudesAmaxNaccelsWBmean_vs_StimAnglePre_mirror

saveas(gca,'attitudesAmaxNaccelsNmaneuverWBs_vs_StimAnglePre_mirror_LnR.fig')
saveas(gca,'attitudesAmaxNaccelsNmaneuverWBs_vs_StimAnglePre_mirror_LnR.png')
plot2svg('attitudesAmaxNaccelsNmaneuverWBs_vs_StimAnglePre_mirror_LnR.svg')

%% attitudesMEAN vs headingPre MIRROR CIRC MEAN
attitudesMEAN_vs_headingPre_mirror


% saveas(gca,'attitudesMEAN_vs_headingPre_mirror.fig')
% saveas(gca,'attitudesMEAN_vs_headingPre_mirror.png')
% plot2svg('attitudesMEAN_vs_headingPre_mirror.svg')

attitudesMEAN_vs_headingPre_mirror_left

saveas(gca,'attitudesMEAN_vs_headingPre_mirror_LnR.fig')
saveas(gca,'attitudesMEAN_vs_headingPre_mirror_LnR.png')
plot2svg('attitudesMEAN_vs_headingPre_mirror_LnR.svg')

%% attitudesMEAN vs StimAnglePre MIRROR CIRC MEAN
attitudesMEAN_vs_StimAnglePre_mirror


% saveas(gca,'attitudesMEAN_vs_StimAnglePre_mirror.fig')
% saveas(gca,'attitudesMEAN_vs_StimAnglePre_mirror.png')
% plot2svg('attitudesMEAN_vs_StimAnglePre_mirror.svg')

attitudesMEAN_vs_StimAnglePre_mirror_left

saveas(gca,'attitudesMEAN_vs_StimAnglePre_mirror_LnR.fig')
saveas(gca,'attitudesMEAN_vs_StimAnglePre_mirror_LnR.png')
plot2svg('attitudesMEAN_vs_StimAnglePre_mirror_LnR.svg')


%% headingNyaw post vs pre CIRC MEAN
plot_headingNyaw_postVSpre

saveas(gca,'headingNyawpost_vs_headingNyawpre.fig')
saveas(gca,'headingNyawpost_vs_headingNyawpre.png')
plot2svg('headingNyawpost_vs_headingNyawpre.svg')


%% headingNyaw post vs pre AnMIRROR CIRC MEAN
plot_headingNyaw_postVSpre_MIRROR

% saveas(gca,'headingNyawpost_vs_headingNyawpre_Anmirror.fig')
% saveas(gca,'headingNyawpost_vs_headingNyawpre_Anmirror.png')
% plot2svg('headingNyawpost_vs_headingNyawpre_Anmirror.svg')

%% headingNyaw post vs pre AnMIRROR CIRCMEAN leftNright
plot_headingNyaw_postVSpre_MIRROR_left

saveas(gca,'headingNyawpost_vs_headingNyawpre_Anmirror_leftNright.fig')
saveas(gca,'headingNyawpost_vs_headingNyawpre_Anmirror_leftNright.png')
plot2svg('headingNyawpost_vs_headingNyawpre_Anmirror_leftNright.svg')

%% turns vs headingNyaw pre CIRC MEAN
% plot_turnVSheadingNyawPRE
% 
% saveas(gca,'turns_vs_headingNyawpre.fig')
% saveas(gca,'turns_vs_headingNyawpre.png')
% plot2svg('turns_vs_headingNyawpre.svg')

%% turns vs headingNyaw pre AnMIRROR CIRC MEAN
plot_turnVSheadingNyawPRE_MIRROR

% saveas(gca,'turns_vs_headingNyawpre_AnMIRROR.fig')
% saveas(gca,'turns_vs_headingNyawpre_AnMIRROR.png')
% plot2svg('turns_vs_headingNyawpre_AnMIRROR.svg')

%% turns vs headingNyaw pre AnMIRROR CIRCMEAN leftNright
plot_turnVSheadingNyawPRE_MIRROR_left

saveas(gca,'turns_vs_headingNyawpre_AnMIRROR_leftNright.fig')
saveas(gca,'turns_vs_headingNyawpre_AnMIRROR_leftNright.png')
plot2svg('turns_vs_headingNyawpre_AnMIRROR_leftNright.svg')

%% AdirNspNmean vs headingNyaw pre CIRC MEAN
% plot_AdirSPnMeanVSheadingNyawPRE
% 
% saveas(gca,'AdirNspN_vs_headingNyawpre.fig')
% saveas(gca,'AdirNspN_vs_headingNyawpre.png')
% plot2svg('AdirNspN_vs_headingNyawpre.svg')

%% AdirNspNmean vs headingNyaw pre AnMIRROR CIRC MEAN
plot_AdirSPnMeanVSheadingNyawPRE_mirror

% saveas(gca,'AdirNspN_vs_headingNyawpre_mirror.fig')
% saveas(gca,'AdirNspN_vs_headingNyawpre_mirror.png')
% plot2svg('AdirNspN_vs_headingNyawpre_mirror.svg')

%% AdirNspNmean vs headingNyaw pre AnMIRROR CIRCMEAN leftNright
plot_AdirSPnMeanVSheadingNyawPRE_mirror_left

saveas(gca,'AdirNspN_vs_headingNyawpre_mirror_LnR.fig')
saveas(gca,'AdirNspN_vs_headingNyawpre_mirror_LnR.png')
plot2svg('AdirNspN_vs_headingNyawpre_mirror_LnR.svg')

%% headingAdirNspNmean vs headingPre AnMIRROR CIRC MEAN
plot_headingAdirSPnMeanVSheadingPRE_mirror

% saveas(gca,'headingAdirNspNmean_vs_headingPre_mirror.fig')
% saveas(gca,'headingAdirNspNmean_vs_headingPre_mirror.png')
% plot2svg('headingAdirNspNmean_vs_headingPre_mirror.svg')

%% headingAdirNspNmean vs headingPre AnMIRROR CIRC MEAN LnR
plot_headingAdirSPnMeanVSheadingPRE_mirror_left

saveas(gca,'headingAdirNspNmean_vs_headingPre_mirror_LnR.fig')
saveas(gca,'headingAdirNspNmean_vs_headingPre_mirror_LnR.png')
plot2svg('headingAdirNspNmean_vs_headingPre_mirror_LnR.svg')

%% stim angle data
plot_response_vs_StimAnglePre

saveas(gca,'response_vs_StimAnglePre_mirror_LnR.fig')
saveas(gca,'response_vs_StimAnglePre_mirror_LnR.png')
plot2svg('response_vs_StimAnglePre_mirror_LnR.svg')

%% WBmod values min & max MIRROR
% plot_WBmodValues_vs_StimAnglePre
% plot_WBmodValues_vs_StimAnglePre_mirror
plot_WBmodValuesWBmean_vs_StimAnglePre_mirror

saveas(gca,'WBmodValues_vs_StimAnglePre_mirror_LnR.fig')
saveas(gca,'WBmodValues_vs_StimAnglePre_mirror_LnR.png')
plot2svg('WBmodValues_vs_StimAnglePre_mirror_LnR.svg')

%% WBmod kinematics @ min & max Angle Accel MIRROR
% plot_WBmodKin_AngleAccel_vs_StimAnglePre_mirror
plot_WBmeanModKin_AngleAccel_vs_StimAnglePre_mirror

saveas(gca,'WBmodKin_AngleAccel_vs_StimAnglePre_mirror_LnR.fig')
saveas(gca,'WBmodKin_AngleAccel_vs_StimAnglePre_mirror_LnR.png')
plot2svg('WBmodKin_AngleAccel_vs_StimAnglePre_mirror_LnR.svg')

%% WBmod kinematics @ Fmax MIRROR
% plot_WBmodKin_Fmax_vs_StimAnglePre_mirror
plot_WBmeanModKin_Fmax_vs_StimAnglePre_mirror

saveas(gca,'WBmodKin_AngleAccel_vs_StimAnglePre_mirror_LnR.fig')
saveas(gca,'WBmodKin_AngleAccel_vs_StimAnglePre_mirror_LnR.png')
plot2svg('WBmodKin_AngleAccel_vs_StimAnglePre_mirror_LnR.svg')

%% maneuver timing t(max) - t(min) MIRROR
% plot_maneuver_AccelTiming_vs_StimAnglePre
% 
% saveas(gca,'maneuver_Acceltiming_vs_StimAnglePre_mirror_LnR.fig')
% saveas(gca,'maneuver_Acceltiming_vs_StimAnglePre_mirror_LnR.png')
% plot2svg('maneuver_Acceltiming_vs_StimAnglePre_mirror_LnR.svg')

%% maneuver wingbeats t(max) - t(min) MIRROR
% plot_maneuverWBs_vs_StimAnglePre
% 
% saveas(gca,'maneuverWBs_vs_StimAnglePre_mirror_LnR.fig')
% saveas(gca,'maneuverWBs_vs_StimAnglePre_mirror_LnR.png')
% plot2svg('maneuverWBs_vs_StimAnglePre_mirror_LnR.svg')

%% rel heading-yaw score VS headingNyaw pre
dm = 5;
dn = 20;

figure
plot(360,1,'ok','MarkerFaceColor','r','markersize',5)
hold on
plot(360,1,'ok','MarkerFaceColor','b','markersize',5)
legend('heading','body angle','location','se')
xlabel('initial angle','fontsize',10) 

stim_angle_velVSyaw = (abs(stim_angle_yaw_post) - abs(stim_angle_vel_post)) ./ (abs(stim_angle_yaw_post) + abs(stim_angle_vel_post));
stim_angle_velVSyaw = [stim_angle_velVSyaw;stim_angle_velVSyaw];

% subplot(1,2,1)
% yaw
angle_pre = [stim_angle_vel_pre;-stim_angle_vel_pre];
plotcolor = 'b';
% plot_headingVSyaw_nowrap_csaps
plot_headingVSyaw_wrap_csaps

% xlabel('initial stim_angle_vel','fontsize',10) 
% ylabel('heading vs yaw score','fontsize',10) 
% ylabel('|body angle|-|heading| / |body angle|+|heading|','fontsize',10) 
% set(gca,'xlim',[-180 180],'ylim',[-1 1])
% set(gca,'XTick',[-180 -90 0 90 180])
% set(gca,'YTick',[-1 -.5 0 .5 1],'fontsize',8)
% grid on

% subplot(1,2,2)
% yaw
angle_pre = [stim_angle_yaw_pre;stim_angle_yaw_pre];
plotcolor = 'r';
% plot_headingVSyaw_nowrap_csaps
plot_headingVSyaw_wrap_csaps

% xlabel('initial yaw','fontsize',10) 
% ylabel('heading vs yaw score','fontsize',10) 
ylabel('|body angle|-|heading| / |body angle|+|heading|','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-1 1])
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-1 -.5 0 .5 1],'fontsize',8)
grid on


saveas(gca,'heading_yaw_score.fig')
saveas(gca,'heading_yaw_score.png')
plot2svg('heading_yaw_score.svg')

%% plot pre&post velocity vectors
V_pre = calc_value(V,n_pre);
V_post = calc_value(V,n_post);

v_pre = V_pre .* cosd(stim_angle_vel_pre);
u_pre = V_pre .* sind(stim_angle_vel_pre);

v_post = V_post .* cosd(stim_angle_vel_post);
u_post = V_post .* sind(stim_angle_vel_post);

figure
% maxHistogramValue = ceil(max(sqrt(u_post.^2 + v_post.^2)));
maxHistogramValue = (max(sqrt(u_post.^2 + v_post.^2)));

subplot(1,2,1)
polar(0, maxHistogramValue,'-k')
% plot(0,0)
% hold on
% plot(0,0,'r')
% legend('pre','post')
% Z = compass_zeroup(u_post,v_post,'r');
Z = compass(u_post,v_post,'r');
for i=1:length(Z)
    set(Z(i),'color','r','linewidth',1.5)
end
hold on
Y = compass(u_pre,v_pre,'k');
for i=1:length(Y)
    set(Y(i),'color','k','linewidth',1.5)
end 

% figure
subplot(1,2,2)
Z = rose(deg2rad(stim_angle_vel_post+90),36);
for i=1:length(Z)
    set(Z(i),'color','k','linewidth',1)
end
x = get(Z,'Xdata');
y = get(Z,'Ydata');
g=patch(x,y,'r');
alpha(.5)

hold on
Y = rose(deg2rad(stim_angle_vel_pre+90),36);
for i=1:length(Y)
    set(Y(i),'color','k','linewidth',1)
end 
x = get(Y,'Xdata');
y = get(Y,'Ydata');
g=patch(x,y,[.5 .5 .5]);
alpha(.5)

compass_zeroup(0, 0.01,'-k');

% title('heading histogram','fontsize',10) 
saveas(gca,'heading_vector_rose_pre_post.fig')
saveas(gca,'heading_rose_pre_post.png')
plot2svg('heading_vector_rose_pre_post.svg')

%% rotation axis VS stim angle & turn angle

plot_rotaxisangle_vs_stimangleNturnangle_atAmax
saveas(gca,'rotaxisangleAmax_vs_stimangleNturnangle_mirror_left.fig')
saveas(gca,'rotaxisangleAmax_vs_stimangleNturnangle_mirror_left.png')
plot2svg('rotaxisangleAmax_vs_stimangleNturnangle_mirror_left.svg')

plot_rotaxisangle_vs_stimangleNturnangle_mean
saveas(gca,'rotaxisangleMEAN_vs_stimangleNturnangle_mirror_left.fig')
saveas(gca,'rotaxisangleMEAN_vs_stimangleNturnangle_mirror_left.png')
plot2svg('rotaxisangleMEAN_vs_stimangleNturnangle_mirror_left.svg')

%% rotation axis roseplot INC stim angle color
plot_rose_rotaxisangle_stimanglecolor_atAmax

saveas(gca,'rotationaxis_angle_rose_mirror_left_stimangle_Amax.fig')
saveas(gca,'rotationaxis_angle_rose_mirror_left_stimangle_Amax.png')
plot2svg('rotationaxis_angle_rose_mirror_left_stimangle_Amax.svg')

%% plot pre&post velocity vectors INC stim angle pre color
V_pre = calc_value(V,n_pre);
V_post = calc_value(V,n_post);

Vdir_pre = stim_angle_vel_mirror_pre;
% Vdir_pre = stim_angle_vel_pre;
Vdir_post = stim_angle_vel_mirror_post;

% v_pre = [V_pre.*cosd(Vdir_pre);V_pre.*cosd(-Vdir_pre)];
% u_pre = [V_pre.*sind(Vdir_pre);V_pre.*sind(-Vdir_pre)];
v_pre = [V_pre.*cosd(Vdir_pre)];
u_pre = [V_pre.*sind(Vdir_pre)];

v_post = V_post .* cosd(Vdir_post);
u_post = V_post .* sind(Vdir_post);

v_post_l = V_post .* cosd(-Vdir_post);
u_post_l = V_post .* sind(-Vdir_post);

cmap_plot = cmap_180;
color_var = round(abs(stim_angle_yaw_mirror_pre));

plot_vector_rose_Vprepost_mirror_left_coloryawpre

% title('heading histogram','fontsize',10) 
saveas(gca,'heading_vector_rose_pre_post_mirror_left_stimangle.fig')
saveas(gca,'heading_vector_rose_pre_post_mirror_left_stimangle.png')
plot2svg('heading_vector_rose_pre_post_mirror_left_stimangle.svg')

%% plot pre&post velocity vectors
V_pre = calc_value(V,n_pre);
V_post = calc_value(V,n_post);

% Vdir_pre = stim_angle_vel_mirror_pre;
Vdir_pre = stim_angle_vel_pre;
Vdir_post = stim_angle_vel_mirror_post;

% v_pre = [V_pre.*cosd(Vdir_pre);V_pre.*cosd(-Vdir_pre)];
% u_pre = [V_pre.*sind(Vdir_pre);V_pre.*sind(-Vdir_pre)];
v_pre = [V_pre.*cosd(Vdir_pre)];
u_pre = [V_pre.*sind(Vdir_pre)];

v_post = V_post .* cosd(Vdir_post);
u_post = V_post .* sind(Vdir_post);

v_post_l = V_post .* cosd(-Vdir_post);
u_post_l = V_post .* sind(-Vdir_post);

plot_vector_rose_Vprepost_mirror_LnR

% title('heading histogram','fontsize',10) 
saveas(gca,'heading_vector_rose_pre_post_mirror_LnR.fig')
saveas(gca,'heading_vector_rose_pre_post_mirror_LnR.png')
plot2svg('heading_vector_rose_pre_post_mirror_LnR.svg')

%% plot Amax vectors
A_hor_first = calc_value(A_hor,n_first);
A_hor_post = calc_value(A_hor,n_post);
A_hor_mean = calc_mean_value(A_hor,n_pre,n_post);
% [A_hor_max,n_Ahormax] = calc_max_value(A_hor,n_pre,n_post);

Adir_mirror_first = calc_value(stim_angle_accel_mirror,n_first);
Adir_first = calc_value(stim_angle_accel,n_first);
Adir_mirror_post = calc_value(stim_angle_accel_mirror,n_post);
Adir_mirror_mean = calc_circ_mean_value(stim_angle_accel_mirror,n_pre,n_post);
Adir_mirror_Ahormax = calc_value(stim_angle_accel_mirror,n_Ahormax);

Vdir_pre = Adir_mirror_first;
Vdir_post = Adir_mirror_Ahormax;

% v_pre = [A_hor_first.*cosd(Adir_mirror_first);A_hor_first.*cosd(-Adir_mirror_first)];
% u_pre = [A_hor_first.*sind(Adir_mirror_first);A_hor_first.*sind(-Adir_mirror_first)];
v_pre = [A_hor_first.*cosd(Adir_first);A_hor_first.*cosd(-Adir_first)];
u_pre = [A_hor_first.*sind(Adir_first);A_hor_first.*sind(-Adir_first)];

v_post = A_hor_max .* cosd(Adir_mirror_Ahormax);
u_post = A_hor_max .* sind(Adir_mirror_Ahormax);

v_post_l = A_hor_max .* cosd(-Adir_mirror_Ahormax);
u_post_l = A_hor_max .* sind(-Adir_mirror_Ahormax);

plot_vector_rose_Vprepost_mirror_LnR

% title('heading histogram','fontsize',10) 
saveas(gca,'AdirMAX_vector_rose_pre_post_mirror_LnR.fig')
saveas(gca,'AdirMAX_vector_rose_pre_post_mirror_LnR.png')
plot2svg('AdirMAX_vector_rose_pre_post_mirror_LnR.svg')

%% plot Amean vectors
Vdir_pre = Adir_mirror_first;
Vdir_post = Adir_mirror_mean;

% v_pre = [A_hor_first.*cosd(Adir_mirror_first);A_hor_first.*cosd(-Adir_mirror_first)];
% u_pre = [A_hor_first.*sind(Adir_mirror_first);A_hor_first.*sind(-Adir_mirror_first)];
v_pre = [A_hor_first.*cosd(Adir_first);A_hor_first.*cosd(-Adir_first)];
u_pre = [A_hor_first.*sind(Adir_first);A_hor_first.*sind(-Adir_first)];

v_post = A_hor_mean .* cosd(Adir_mirror_mean);
u_post = A_hor_mean .* sind(Adir_mirror_mean);

v_post_l = A_hor_mean .* cosd(-Adir_mirror_mean);
u_post_l = A_hor_mean .* sind(-Adir_mirror_mean);

plot_vector_rose_Vprepost_mirror_LnR

% title('heading histogram','fontsize',10) 
saveas(gca,'AdirMEAN_vector_rose_pre_post_mirror_LnR.fig')
saveas(gca,'AdirMEAN_vector_rose_pre_post_mirror_LnR.png')
plot2svg('AdirMEAN_vector_rose_pre_post_mirror_LnR.svg')

%% plot pre&post Adir vectors MIRROR
A_hor_first = calc_value(A_hor,n_first);
A_hor_post = calc_value(A_hor,n_post);
A_hor_mean = calc_mean_value(A_hor,n_pre,n_post);
[A_hor_max,n_Ahormax] = calc_max_value(A_hor,n_pre,n_post);

Adir_first = calc_value(stim_angle_accel_mirror,n_first);
Adir_post = calc_value(stim_angle_accel_mirror,n_post);
Adir_mean = calc_circ_mean_value(stim_angle_accel_mirror,n_pre,n_post);
Adir_Ahormax = calc_value(stim_angle_accel_mirror,n_Ahormax);

% maxHistogramValue = (max(A_hor_max));
maxHistogramValue = 20;

plot_vector_rose_Adata_mirror

% saveas(gca,'Adir_vector_rose_first_max_mean_mirror.fig')
% saveas(gca,'Adir_vector_rose_first_max_mean_mirror.png')
% plot2svg('Adir_vector_rose_first_max_mean_mirror.svg')

%% plot pre&post Adir vectors MIRROR LnR
A_hor_first = calc_value(A_hor,n_first);
A_hor_post = calc_value(A_hor,n_post);
A_hor_mean = calc_mean_value(A_hor,n_pre,n_post);
[A_hor_max,n_Ahormax] = calc_max_value(A_hor,n_pre,n_post);

Adir_first = -calc_value(stim_angle_accel_mirror,n_first);
Adir_post = -calc_value(stim_angle_accel_mirror,n_post);
Adir_mean = -calc_circ_mean_value(stim_angle_accel_mirror,n_pre,n_post);
Adir_Ahormax = -calc_value(stim_angle_accel_mirror,n_Ahormax);

% maxHistogramValue = (max(A_hor_max));
maxHistogramValue = 20;

plot_vector_rose_Adata_mirror_left

saveas(gca,'Adir_vector_rose_first_max_mean_mirror_LnR.fig')
saveas(gca,'Adir_vector_rose_first_max_mean_mirror_LnR.png')
plot2svg('Adir_vector_rose_first_max_mean_mirror_LnR.svg')




%% plot FORCE data at Amax NO MIRROR
% stim_angle_F_now = calc_value(stim_angle_F,n_Amax);
% Fsp_roll_now = calc_value(Fsp_roll,n_Amax);
% Fsp_pitch_now = calc_value(Fsp_pitch,n_Amax);

% roll_now = calc_value(roll,n_Amax);
% pitch_now = calc_value(pitch,n_Amax);
% yaw_now = calc_value(yaw,n_Amax);

% F_now = calc_value(F,n_Amax);
% Fn_hor_now = calc_value(Fn_hor,n_Amax);
% Ft_hor_now = calc_value(Ft_hor,n_Amax);
% 
% figure
% maxHistogramValue = 10;
% n_bins_rose = 90;
% plot_color = 'r';
% dh = .1;
% 
% plot_rose_hist_data
% 
% saveas(gca,'Force_data_atAmax.fig')
% saveas(gca,'Force_data_atAmax.png')
% plot2svg('Force_data_atAmax.svg')

%% plot MEAN FORCE data NO MIRROR
% stim_angle_F_now = calc_circ_mean_value(stim_angle_F,n_pre,n_post);
% Fsp_roll_now = calc_circ_mean_value(Fsp_roll,n_pre,n_post);
% Fsp_pitch_now = calc_circ_mean_value(Fsp_pitch,n_pre,n_post);

% roll_now = calc_circ_mean_value(roll,n_pre,n_post);
% pitch_now = calc_circ_mean_value(pitch,n_pre,n_post);
% yaw_now = calc_circ_mean_value(yaw,n_pre,n_post);

% F_now = calc_mean_value(F,n_pre,n_post);
% Fn_hor_now = calc_mean_value(Fn_hor,n_pre,n_post);
% Ft_hor_now = calc_mean_value(Ft_hor,n_pre,n_post);
% 
% figure
% maxHistogramValue = 10;
% n_bins_rose = 90;
% plot_color = 'r';
% dh = .1;
% 
% plot_rose_hist_data
% 
% saveas(gca,'Force_data_MEAN.fig')
% saveas(gca,'Force_data_MEAN.png')
% plot2svg('Force_data_MEAN.svg')

%% plot FORCE data @Amax MIRROR LnR
stim_angle_F_now = calc_value(stim_angle_F_mirror,n_Amax);
Fsp_roll_now = calc_value(Fsp_roll_mirror,n_Amax);
Fsp_pitch_now = calc_value(Fsp_pitch,n_Amax);

roll_now = calc_value(roll_mirror,n_Amax);
pitch_now = calc_value(pitch,n_Amax);
yaw_now = calc_value(yaw_mirror,n_Amax);

F_now = calc_value(F,n_Amax);
Fn_hor_now = calc_value(Fn_hor_mirror,n_Amax);
Ft_hor_now = calc_value(Ft_hor,n_Amax);

figure
maxHistogramValue = 20;
n_bins_rose = 90;
plot_color = 'b';
dh = .1;

plot_rose_hist_data

subplot(3,3,1)
compass_zeroup(0, 0.01,'-k');
subplot(3,3,4)
compass_zeroup(0, 0.01,'-k');
% subplot(3,3,7)
% compass_zeroup(0, 0.01,'-k');
subplot(3,3,2)
compass_zeroup(0, 0.01,'-k');
% subplot(3,3,5)
% compass_zeroup(0, 0.01,'-k');
subplot(3,3,8)
compass_zeroup(0, 0.01,'-k');

% 
% saveas(gca,'Force_data_atAmax_mirror.fig')
% saveas(gca,'Force_data_atAmax_mirror.png')
% plot2svg('Force_data_atAmax_mirror.svg')

stim_angle_F_now = -calc_value(stim_angle_F_mirror,n_Amax);
Fsp_roll_now = -calc_value(Fsp_roll_mirror,n_Amax);
Fsp_pitch_now = calc_value(Fsp_pitch,n_Amax);

roll_now = -calc_value(roll_mirror,n_Amax);
pitch_now = calc_value(pitch,n_Amax);
yaw_now = -calc_value(yaw_mirror,n_Amax);

F_now = calc_value(F,n_Amax);
Fn_hor_now = -calc_value(Fn_hor_mirror,n_Amax);
Ft_hor_now = calc_value(Ft_hor,n_Amax);

maxHistogramValue = 10;
n_bins_rose = 90;
plot_color = 'r';
dh = .1;

plot_rose_hist_data

saveas(gca,'Force_data_atAmax_mirror_LnR.fig')
saveas(gca,'Force_data_atAmax_mirror_LnR.png')
plot2svg('Force_data_atAmax_mirror_LnR.svg')


%% plot FORCE data MEAN MIRROR LnR
stim_angle_F_now = calc_circ_mean_value(stim_angle_F_mirror,n_pre,n_post);
Fsp_roll_now = calc_circ_mean_value(Fsp_roll_mirror,n_pre,n_post);
Fsp_pitch_now = calc_circ_mean_value(Fsp_pitch,n_pre,n_post);

roll_now = calc_circ_mean_value(roll_mirror,n_pre,n_post);
pitch_now = calc_circ_mean_value(pitch,n_pre,n_post);
yaw_now = calc_circ_mean_value(yaw_mirror,n_pre,n_post);

F_now = calc_mean_value(F,n_pre,n_post);
Fn_hor_now = calc_mean_value(Fn_hor_mirror,n_pre,n_post);
Ft_hor_now = calc_mean_value(Ft_hor,n_pre,n_post);

figure
maxHistogramValue = 20;
n_bins_rose = 90;
plot_color = 'b';
dh = .1;

plot_rose_hist_data

subplot(3,3,1)
compass_zeroup(0, 0.01,'-k');
subplot(3,3,4)
compass_zeroup(0, 0.01,'-k');
subplot(3,3,7)
compass_zeroup(0, 0.01,'-k');
subplot(3,3,2)
compass_zeroup(0, 0.01,'-k');
subplot(3,3,5)
compass_zeroup(0, 0.01,'-k');
subplot(3,3,8)
compass_zeroup(0, 0.01,'-k');

% 
% saveas(gca,'Force_data_MEAN_mirror.fig')
% saveas(gca,'Force_data_MEAN_mirror.png')
% plot2svg('Force_data_MEAN_mirror.svg')

stim_angle_F_now = -calc_circ_mean_value(stim_angle_F_mirror,n_pre,n_post);
Fsp_roll_now = -calc_circ_mean_value(Fsp_roll_mirror,n_pre,n_post);
Fsp_pitch_now = calc_circ_mean_value(Fsp_pitch,n_pre,n_post);

roll_now = -calc_circ_mean_value(roll_mirror,n_pre,n_post);
pitch_now = calc_circ_mean_value(pitch,n_pre,n_post);
yaw_now = -calc_circ_mean_value(yaw_mirror,n_pre,n_post);

F_now = calc_mean_value(F,n_pre,n_post);
Fn_hor_now = -calc_mean_value(Fn_hor_mirror,n_pre,n_post);
Ft_hor_now = calc_mean_value(Ft_hor,n_pre,n_post);

maxHistogramValue = 10;
n_bins_rose = 90;
plot_color = 'r';
dh = .1;

plot_rose_hist_data

saveas(gca,'Force_data_MEAN_mirror_LnR.fig')
saveas(gca,'Force_data_MEAN_mirror_LnR.png')
plot2svg('Force_data_MEAN_mirror_LnR.svg')






%% rose plots heading yaw POST
% 
% figure
% % subplot(1,2,1)
% Y = rose(deg2rad(heading_post),36);
% for i=1:length(Y)
%     set(Y(i),'color','k','linewidth',1)
% end 
% x = get(Y,'Xdata');
% y = get(Y,'Ydata');
% g=patch(x,y,'r');
% % alpha(.5)
% hold on
% 
% Z = rose(deg2rad(yaw_post),36);
% for i=1:length(Z)
%     set(Z(i),'color','k','linewidth',1)
% end
% x = get(Z,'Xdata');
% y = get(Z,'Ydata');
% g=patch(x,y,'b');
% alpha(.5)
% 
% legend('','heading','','yaw')
% 
% 
% saveas(gca,'rose_heading_yaw_post.fig')
% saveas(gca,'rose_heading_yaw_post.png')
% plot2svg('rose_heading_yaw_post.svg')
% 
% 
%% rose plots heading yaw PRE&POST
% 
% 
% figure
% subplot(1,2,1)
% Z = rose(deg2rad(heading_post),36);
% for i=1:length(Z)
%     set(Z(i),'color','r','linewidth',2)
% end
% % x = get(Z,'Xdata');
% % y = get(Z,'Ydata');
% % g=patch(x,y,'r');
% 
% hold on
% Y = rose(deg2rad(heading_pre),36);
% for i=1:length(Y)
%     set(Y(i),'color','k','linewidth',2)
% end 
% % x = get(Y,'Xdata');
% % y = get(Y,'Ydata');
% % g=patch(x,y,'b');
% title('heading','fontsize',10) 
% 
% % figure
% subplot(1,2,2)
% Z = rose(deg2rad(yaw_post),36);
% for i=1:length(Z)
%     set(Z(i),'color','b','linewidth',2)
% end
% % x = get(Z,'Xdata');
% % y = get(Z,'Ydata');
% % g=patch(x,y,'r');
% 
% hold on
% Y = rose(deg2rad(yaw_pre),36);
% for i=1:length(Y)
%     set(Y(i),'color','k','linewidth',2)
% end 
% % x = get(Y,'Xdata');
% % y = get(Y,'Ydata');
% % g=patch(x,y,'b');
% title('yaw','fontsize',10) 
% 
% saveas(gca,'rose_heading_yaw_pre_post.fig')
% saveas(gca,'rose_heading_yaw_pre_post.png')
% plot2svg('rose_heading_yaw_pre_post.svg')
% 
%% plot rose teta_max ON/OFF
% 
% maxHistogramValue = 16;
% figure
% 
% plot_rose_heading_tetamax_ONOFF
% 
% saveas(gca,'heading_rose_post_turn.fig')
% saveas(gca,'heading_rose_pre_post.png')
% plot2svg('heading_rose_pre_post.svg')
% 
%% plot pre&post vertical velocity vectors
% % figure
% % compass_zeroright(sqrt(u_post.^2 + v_post.^2),w_post,'r')
% % hold on
% % compass(sqrt(u_pre.^2 + v_pre.^2),w_pre,'b')
% subplot(1,2,2)
% Z = compass(sqrt(u_post.^2 + v_post.^2),w_post);
% for i=1:length(Z)
%     set(Z(i),'color',cmap(end,:))
% end
% hold on
% Y = compass(sqrt(u_pre.^2 + v_pre.^2),w_pre);
% for i=1:length(Y)
%     set(Y(i),'color',cmap(1,:))
% end 
% title('vertical flight speed vectors')
% saveas(gca,'V_PreVsPost.fig')
% saveas(gca,'V_PreVsPost.png')
% 
%% plot histograms heading&yaw POST
% 
% figure
% % subplot(3,2,1)
% 
% dh = 360/(36)
% bins = [-180+dh/2:dh:180-dh/2];
% 
% h1 = hist(heading_post,bins);
% h2 = hist(yaw_post,bins);
% 
% bar(bins,h1,'FaceColor','r','EdgeColor','k')
% hold on
% bar(bins,h2,'FaceColor','b','EdgeColor','k')
% % bar(bins,[h1;h2]')
% 
% legend('heading','yaw')
% alpha(.5)
% set(gca,'xlim',[-180 180])
% set(gca,'XTick',[-180:90:180])
% set(gca,'YTick',[0:25:25],'fontsize',8) 
% title('escape direction','fontsize',10) 
% 
% saveas(gca,'hist_heading_yaw_post.fig')
% saveas(gca,'hist_heading_yaw_post.png')
% plot2svg('hist_heading_yaw_post.svg')
% 
% set(gca,'xlim',[-90 90])
% set(gca,'XTick',[-180:90:180])
% 
% saveas(gca,'hist_heading_yaw_post_min90to90.fig')
% saveas(gca,'hist_heading_yaw_post_min90to90.png')
% plot2svg('hist_heading_yaw_post_min90to90.svg')
% 
% figure
% subplot(2,1,1)
% bar(bins,h1,'FaceColor',[.5 .5 .5],'EdgeColor','k')
% hold on
% set(gca,'xlim',[-180 180])
% set(gca,'XTick',[-180:90:180])
% set(gca,'YTick',[0:25:25],'fontsize',8) 
% title('escape heading angle','fontsize',10) 
% 
% subplot(2,1,2)
% bar(bins,h2,'FaceColor',[.5 .5 .5],'EdgeColor','k')
% % bar(bins,[h1;h2]')
% 
% % legend('heading','yaw')
% % alpha(.5)
% set(gca,'xlim',[-180 180])
% set(gca,'XTick',[-180:90:180])
% set(gca,'YTick',[0:25:25],'fontsize',8) 
% title('escape body angle','fontsize',10) 
% 
% saveas(gca,'hist_heading_yaw_post_sep.fig')
% saveas(gca,'hist_heading_yaw_post_sep.png')
% plot2svg('hist_heading_yaw_post_sep.svg')
% 
% % -90 to 90 deg
% subplot(2,1,1)
% set(gca,'xlim',[-90 90])
% set(gca,'XTick',[-180:90:180])
% 
% subplot(2,1,2)
% set(gca,'xlim',[-90 90])
% set(gca,'XTick',[-180:90:180])
% 
% saveas(gca,'hist_heading_yaw_post_sep_min90to90.fig')
% saveas(gca,'hist_heading_yaw_post_sep_min90to90.png')
% plot2svg('hist_heading_yaw_post_sep_min90to90.svg')
% 
%% heading vs yaw post (heading pre colorcode)
% figure
% hold on
% for i = 1:length(heading_post)
%     plot(heading_post(i),yaw_post(i),'ok','MarkerFaceColor',cmap_360(round(heading_pre(i))+180,:),'linewidth',1)
% end
% axis equal
% grid on
% xlabel('escape heading','fontsize',10) 
% ylabel('escape body angle','fontsize',10) 
% set(gca,'xlim',[-180 180],'ylim',[-180 180])
% set(gca,'XTick',[-180:90:180],'YTick',[-180:90:180],'fontsize',8) 
% 
% colormap(cmap_360(1:2:end,:))
% colorbar('YTick',[1:45:181],'YTickLabel',[-180:90:180],'fontsize',8) 
% %title('escape body angle','fontsize',10) 
% 
% 
% saveas(gca,'heading_yaw_post_headingpre_colorcode.fig')
% saveas(gca,'heading_yaw_post_headingpre_colorcode.png')
% plot2svg('heading_yaw_post_headingpre_colorcode.svg')

%% plot histograms
% 
% % figure
% subplot(3,2,1)
% dh = 360/(36)
% heading_noresp = heading(1:Nnoresp,:);
% heading_noresp = heading_noresp(isnan(heading_noresp)==0);
% hist(heading_noresp,-180+dh/2:dh:180-dh/2)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor',cmap(1,:),'EdgeColor','w')
% 
% 
% title('heading pre')
% set(gca,'ylim',[0 10000])
% 
% set(gca,'xlim',[-180 180])
% % set(gca,'XTick',[-180 -90 0 90 180])
% % saveas(gca,'hist_preheading.fig')
% % saveas(gca,'hist_preheading.png')
% % plot2svg
% 
% % figure
% subplot(3,2,2)
% heading_resp = heading(Nresp:end,:);
% heading_resp = heading_resp(isnan(heading_resp)==0);
% hist(heading_resp,-180+dh/2:dh:180-dh/2)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor',cmap(1,:),'EdgeColor','w')
% 
% title('heading post')
% % set(gca,'ylim',[0 10000])
% 
% set(gca,'xlim',[-180 180])
% set(gca,'XTick',[-180 -90 0 90 180])
% % saveas(gca,'hist_postheading.fig')
% % saveas(gca,'hist_postheading.png')
% % plot2svg
% 
% % figure
% subplot(3,2,3)
% dV = 1/25
% V_noresp = V(1:Nnoresp,:);
% V_noresp = V_noresp(isnan(V_noresp)==0);
% hist(V_noresp,0+dV/2:dV:1-dV/2)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor',cmap(1,:),'EdgeColor','w')
% 
% title('velocity pre')
% % set(gca,'ylim',[0 50000])
% 
% set(gca,'xlim',[0 1])
% set(gca,'XTick',[0 .5 1])
% % saveas(gca,'hist_preV.fig')
% % saveas(gca,'hist_preV.png')
% % plot2svg
% 
% % figure
% subplot(3,2,4)
% V_resp = V(Nresp:end,:);
% V_resp = V_resp(isnan(V_resp)==0);
% hist(V_resp,0+dV/2:dV:1-dV/2)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor',cmap(1,:),'EdgeColor','w')
% 
% title('velocity post')
% % set(gca,'ylim',[0 10000])
% 
% set(gca,'xlim',[0 1])
% set(gca,'XTick',[0 .5 1])
% % saveas(gca,'hist_postV.fig')
% % saveas(gca,'hist_postV.png')
% % plot2svg
% 
% % figure
% subplot(3,2,5)
% dA = 30/25
% A_noresp = A(1:Nnoresp,:);
% A_noresp = A_noresp(isnan(A_noresp)==0);
% hist(A_noresp,0+dA/2:dA:30-dA/2)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor',cmap(1,:),'EdgeColor','w')
% 
% title('accel pre')
% % set(gca,'ylim',[0 100000])
% 
% set(gca,'xlim',[0 30])
% set(gca,'XTick',[0 10 20 30])
% % saveas(gca,'hist_preA.fig')
% % saveas(gca,'hist_preA.png')
% % plot2svg
% 
% % figure
% subplot(3,2,6)
% A_resp = A(Nresp:end,:);
% A_resp = A_resp(isnan(A_resp)==0);
% hist(A_resp,0+dA/2:dA:30-dA/2)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor',cmap(1,:),'EdgeColor','w')
% 
% title('accel post')
% % set(gca,'ylim',[0 10000])
% 
% set(gca,'xlim',[0 30])
% set(gca,'XTick',[0 10 20 30])
% % saveas(gca,'hist_postA.fig')
% % saveas(gca,'hist_postA.png')
% % plot2svg
% 
% saveas(gca,'histograms_pre_post.fig')
% saveas(gca,'histograms_pre_post.png')
% % plot2svg
% 
%% plot top view flight path vectors
% figure
% colormap(cmap)
% caxis([Nstart Nstop])
% hold on
% for j=1:size(x,2)
%     plot(x(Nstart:Nstop,j),y(Nstart:Nstop,j),'-','color',cmap(1,:))
% end
% 
% vec_col = 0;
% for i=Nstart:dN:Nstop
%     vec_col = vec_col + size(cmap,1)/Nstop;
%     for j=1:size(x,2)
% %         quiverc(x(i,:),y(i,:),u(i,:),v(i,:),'color',cmap(1,:))
%         if isnan(x(i,j)) == 0
%             quivert(x(i,j),y(i,j),u(i,j),v(i,j),i,'as',1/500,'ahr',[1 1],'nt')
%         end
%     end
%     
% end
%     
% axis equal
% colorbar
% 
% % set(gca,'xlim',[0 30])
% % set(gca,'XTick',[0 10 20 30])
% saveas(gca,'flighttracks.fig')
% saveas(gca,'flighttracks.png')
% % plot2svg

%% ON/OFF escape performance: heading, yaw, resp time, U, An, At
% 
% plot_escape_perf_ONOFF_NO165tetamax
% 
% % title('acceleration based response time')
% saveas(gca,'escape performance_ONOFF.fig')
% saveas(gca,'escape performance_ONOFF.png')
% plot2svg('escape performance_ONOFF.svg')

%% teta_max escape performance: heading, yaw, resp time, U, An, At
% 
% plot_escape_perf_tetamax
% 
% % title('acceleration based response time')
% saveas(gca,'escape performance_tetamax.fig')
% saveas(gca,'escape performance_tetamax.png')
% plot2svg('escape performance_tetamax.svg')
% 
cd ..