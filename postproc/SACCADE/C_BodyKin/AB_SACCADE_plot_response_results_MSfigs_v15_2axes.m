clc
clear
close all

%% rotation & counter-rotation time windows
t_rot1_min = -.025;
t_rot1_max = -.005;

t_rot2_min = .03;
t_rot2_max = .06;

t_start = -.05;
t_stop = .06;
t_cut = .0125;

%% load data

% load WBmod data
% list = dir('WBmod*.mat')
% for i = 1:length(list)
%     load(list(i).name)
% end
% steady_name = dir('WBdataset_steady_*')
% load(steady_name.name)
% 
% all_name = dir('WBdataset_all_9*')
% load(all_name.name)

% loadname = 'flightpathDB_pos_qbodyEKF_INCroll_strokeplane47.5deg_rotaxis27.9deg_startframe100.mat';
% loadname = 'flightpathDB_pos_qbodyEKF_INCroll_9clusters_1.4n-1.7n1.9_strokeplane47.5deg_startframe100.mat';
loadname = 'flightpathDB_pos_qbodyEKF_INCroll_strokeplane47.5deg_rotaxis22deg_startframe100.mat';
load(loadname)
% 
% % loadname=dir('kinflightpathDB_pos_qbodyEKF_INCroll_9clusters*')
% % % if exist(loadname) ~= 2
% % % %     loadname=('flightpathDB_pos_qbodyEKF_NOroll_9clusters_2.5n-3.3n3_startframe2945.mat')
% % %     loadname=('flightpathDB_pos_qbodyEKF_INCroll_9clusters_2.5n-3.3n3_startframe2945_strokeplane45deg.mat')
% % %     loadname=('flightpathDB_pos_qbodyEKF_INCroll_9clusters_2.5n-3.3n3_startframe2945_strokeplane55deg.mat')
% % % end
% % load(loadname.name)

% load('norm_data.mat')
load('norm_data_steady.mat') % from 1603 wbs

% load('flyVars_addedmassSphere.mat')
% load('flyVars_addedmassDisc.mat')
% load('flyVars_addedmassDiscSphere.mat')
const_name=dir('flyVars*')
const_name=const_name.name;
load(const_name)

%%
% mkdir('MSfigs_BodyDyn')
% cd('MSfigs_BodyDyn')

plot_timelines = 0
plot_timelines = 1

plot_paths = 0
plot_paths = 1

make_mov = 0
% make_mov = 1

cluster_on = 0
% cluster_on = 1

%% subset
% % subset_seqs=[7;65;4;13;11]; % old seq set
% % subset_seqs=[4;62;1;10;8]; % updated seq set
% % subset_seqs = [59;26;48;66;30;92;22]; % new seq set
% subset_seqs = [59;1;48;30;66;76;19;22]; % merged sets
% 
% seqs_NOmirror = [7;64;21;4;13;11];

%% CLUSTERS

% 2 MEAN clusters SACCADE
subset_cut = [90];
subset_mid = round([45 135]);

% % 3 MEAN clusters SACCADE
% subset_cut = [75 105];
% subset_mid = round([45 90 135]);

% % 5 MEAN clusters
% subset_cut = [36 72 108 144];
% subset_mid = [18 54 90 126 162];

% % 4 MEAN clusters
% subset_cut = [45 90 135];
% subset_mid = round([22.5 67.5 112.5 157.5]);

% % 5 MEAN clusters ON/OFF
% subset_cut = [20 40 60 80 100 120 140 160];
% subset_mid = round([10 50 90 130 170]);
% 
% % 4 MEAN clusters ONOFF
% dc = 180/7;
% subset_cut = [dc 2*dc 3*dc 4*dc 5*dc 6*dc];
% subset_mid = round([.5*dc 2.5*dc 4.5*dc 6.5*dc]);
% 
% % 3 MEAN clusters ONOFF
% dc = 180/5;
% subset_cut = [dc 2*dc 3*dc 4*dc];
% subset_mid = round([.5*dc 2.5*dc 4.5*dc]);

% % 4 MEAN clusters ONOFF
% dc = 180/22;
% subset_cut = [4*dc 6*dc 10*dc 12*dc 16*dc 18*dc];
% subset_mid = round([2*dc 8*dc 14*dc 20*dc]);


%% !!! mirror exception !!!
% except = 60;
except = [];

%% settings
rot_axis_angle = settings.rot_axis_angle;
g=settings.g;

linewidth_timelines = 1;
skip = 50;

% heatmap resolution
nx = 1000;
ny = 100;
% nx = 2000;
% ny = 200;

cmap_180 = jet(180);
cmap_bw = 1-colormap(gray);

grey_color = [.65 .65 .65];

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

pos_x = pathDB.pos(:,:,1);
pos_y = pathDB.pos(:,:,2);
pos_z = pathDB.pos(:,:,3);

stim_angle_vel = pathDB.stim_angle_vel;
stim_angle_accel = pathDB.stim_angle_accel;
stim_angle_F = pathDB.stim_angle_F;
stim_angle_yaw = pathDB.stim_angle_yaw;
stim_angle_spn = pathDB.stim_angle_spn;

angle_vel = rad2deg(unwrap(deg2rad(pathDB.angle_vel)));
angle_accel = rad2deg(unwrap(deg2rad(pathDB.angle_accel)));
angle_F = rad2deg(unwrap(deg2rad(pathDB.angle_F)));
angle_yaw = rad2deg(unwrap(deg2rad(pathDB.angle_yaw)));
angle_spn = rad2deg(unwrap(deg2rad(pathDB.angle_spn)));

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

drot_body_L = pathDB.drot_body_L;
drot_body_R = pathDB.drot_body_R;

drot_L = pathDB.drot_sp_L;
drot_R = pathDB.drot_sp_R;

rot_dot_body_L = pathDB.rot_dot_body_L;
rot_dot_body_R = pathDB.rot_dot_body_R;

rot_dot_L = pathDB.rot_dot_sp_L;
rot_dot_R = pathDB.rot_dot_sp_R;

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
    if isnan(n_pre(i)) == 1
        V_temp = V(:,i);
        n_pre(i,1) = find(isnan(V_temp)==0, 1, 'first' );
    end
end

%% angular accels
yaw_dot_dot = nan(size(yaw_dot));
pitch_dot_dot = nan(size(pitch_dot));
roll_dot_dot = nan(size(roll_dot));

rot_dot_dot_L = nan(size(rot_dot_L));
rot_dot_dot_R = nan(size(rot_dot_R));

yaw_dot_dot(2:end-1,:) = yaw_dot(3:end,:) - yaw_dot(1:end-2,:);
pitch_dot_dot(2:end-1,:) = pitch_dot(3:end,:) - pitch_dot(1:end-2,:);
roll_dot_dot(2:end-1,:) = roll_dot(3:end,:) - roll_dot(1:end-2,:);

rot_dot_dot_L(2:end-1,:) = rot_dot_L(3:end,:) - rot_dot_L(1:end-2,:);
rot_dot_dot_R(2:end-1,:) = rot_dot_R(3:end,:) - rot_dot_R(1:end-2,:);

yaw_dot_dot = yaw_dot_dot/2/dt;
pitch_dot_dot = pitch_dot_dot/2/dt;
roll_dot_dot = roll_dot_dot/2/dt;

rot_dot_dot_L = rot_dot_dot_L/2/dt;
rot_dot_dot_R = rot_dot_dot_R/2/dt;

%% WBmeans
freq_steady = f_wb_steady;
NsteadyWB = round(1/freq_steady/dt);

% freq_steady = f_wb_steady_meanCIstd(1);
% NsteadyWB = round(1/freq_steady/dt);
% 
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

% dAngles
dAngle_vel = nan(size(V));
dAngle_accel = nan(size(V));
dAngle_F = nan(size(V));
dAngle_yaw = nan(size(V));
dAngle_spn = nan(size(V));

% initial conditions
V_trig2resp = calc_mean_value(V,n_first,n_resp);

angle_vel_trig2resp = calc_mean_value(angle_vel,n_first,n_resp);
angle_accel_trig2resp = calc_mean_value(angle_accel,n_first,n_resp);
angle_F_trig2resp = calc_mean_value(angle_F,n_first,n_resp);
angle_yaw_trig2resp = calc_mean_value(angle_yaw,n_first,n_resp);
angle_spn_trig2resp = calc_mean_value(angle_spn,n_first,n_resp);

roll_global_first = calc_value(roll_global,n_first);
slip_global_first = calc_value(slip_global,n_first);
pitch_global_first = calc_value(pitch_global,n_first);

pitch_global_first_minsp = pitch_global_first + strokeplane_angle;

% clear dV
for i = 1:length(V_trig2resp)
    dV(:,i) = V(:,i) - V_trig2resp(i);
    
    dAngle_vel(:,i) = angle_vel(:,i) - angle_vel_trig2resp(i);
    
%     dAngle_accel(:,i) = angle_accel(:,i) - angle_accel_trig2resp(i);
%     dAngle_F(:,i) = angle_F(:,i) - angle_F_trig2resp(i);
%     dAngle_yaw(:,i) = angle_yaw(:,i) - angle_yaw_trig2resp(i);
%     dAngle_spn(:,i) = angle_spn(:,i) - angle_spn_trig2resp(i);

    dAngle_accel(:,i) = angle_accel(:,i) - angle_vel_trig2resp(i);  % !!! initial angle based on vel !!!
    dAngle_F(:,i) = angle_F(:,i) - angle_vel_trig2resp(i);  % !!! initial angle based on vel !!!
    dAngle_yaw(:,i) = angle_yaw(:,i) - angle_vel_trig2resp(i);
    dAngle_spn(:,i) = angle_spn(:,i) - angle_vel_trig2resp(i);
    
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

dV_post = calc_value(dV,n_post);

%% !!!!!!!!!!!!!! SACCADE stim_angle = dAngle !!!!!!!!!!!!!!!!!!!
stim_angle_vel_global = stim_angle_vel;
stim_angle_vel = dAngle_vel;
stim_angle_accel = dAngle_accel;
stim_angle_F = dAngle_F;
stim_angle_yaw = dAngle_yaw;
stim_angle_spn = dAngle_spn;

%% MIRROR TURN based on YAW_MAX
dyaw_post = calc_value(dyaw,n_post);

pos_y_mirror = pos_y;
stim_angle_vel_global_mirror = stim_angle_vel_global;

stim_angle_vel_mirror = stim_angle_vel;
stim_angle_accel_mirror = stim_angle_accel;
stim_angle_yaw_mirror = stim_angle_yaw;
stim_angle_F_mirror = stim_angle_F;
stim_angle_spn_mirror = stim_angle_spn;

dAngle_vel_mirror = dAngle_vel;
dAngle_accel_mirror = dAngle_accel;
dAngle_yaw_mirror = dAngle_yaw;
dAngle_F_mirror = dAngle_F;
dAngle_spn_mirror = dAngle_spn;

An_hor_mirror = An_hor;
Fn_hor_mirror = Fn_hor;

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
Fsp_roll_mirror = Fsp_roll;
% stroke_L_mirror = stroke_L;
% stroke_R_mirror = stroke_R;
% pitch_L_mirror = pitch_L;
% pitch_R_mirror = pitch_R;
% dev_L_mirror = dev_L;
% dev_R_mirror = dev_R;
% U_L_mirror = U_L;
% U_R_mirror = U_R;
% aoa_L_mirror = aoa_L;
% aoa_R_mirror = aoa_R;
% 
% drot_L_mirror = drot_L;
% drot_R_mirror = drot_R;
% rot_dot_L_mirror = rot_dot_L;
% rot_dot_R_mirror = rot_dot_R;
% rot_dot_dot_L_mirror = rot_dot_dot_L;
% rot_dot_dot_R_mirror = rot_dot_dot_R;
% 
% drot_L_temp = drot_L;
% drot_R_temp = drot_R;
% rot_dot_L_temp = rot_dot_L;
% rot_dot_R_temp = rot_dot_R;
% rot_dot_dot_L_temp = rot_dot_dot_L;
% rot_dot_dot_R_temp = rot_dot_dot_R;


for i=1:length(dyaw_post)
    if dyaw_post(i) < 0
        
        pos_y_mirror(:,i) = -pos_y(:,i);
        stim_angle_vel_global_mirror(:,i) = -stim_angle_vel_global_mirror(:,i);

        stim_angle_vel_mirror(:,i) = -stim_angle_vel_mirror(:,i);
        stim_angle_accel_mirror(:,i) = -stim_angle_accel_mirror(:,i);
        stim_angle_yaw_mirror(:,i) = -stim_angle_yaw_mirror(:,i);
        stim_angle_F_mirror(:,i) = -stim_angle_F_mirror(:,i);
        stim_angle_spn_mirror(:,i) = -stim_angle_spn_mirror(:,i);
        
        dAngle_vel_mirror(:,i) = -dAngle_vel_mirror(:,i);
        dAngle_accel_mirror(:,i) = -dAngle_accel_mirror(:,i);
        dAngle_yaw_mirror(:,i) = -dAngle_yaw_mirror(:,i);
        dAngle_F_mirror(:,i) = -dAngle_F_mirror(:,i);
        dAngle_spn_mirror(:,i) = -dAngle_spn_mirror(:,i);
        
        An_hor_mirror(:,i) = -An_hor_mirror(:,i);
        Fn_hor_mirror(:,i) = -Fn_hor_mirror(:,i);

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
        Fsp_roll_mirror(:,i) = -Fsp_roll_mirror(:,i);
%         stroke_L_mirror(:,i) = stroke_R(:,i);
%         stroke_R_mirror(:,i) = stroke_L(:,i);
%         pitch_L_mirror(:,i) = pitch_R(:,i);
%         pitch_R_mirror(:,i) = pitch_L(:,i);
%         dev_L_mirror(:,i) = dev_R(:,i);
%         dev_R_mirror(:,i) = dev_L(:,i);
%         U_L_mirror(:,i) = U_R(:,i);
%         U_R_mirror(:,i) = U_L(:,i);
%         aoa_L_mirror(:,i) = aoa_R(:,i);
%         aoa_R_mirror(:,i) = aoa_L(:,i);
% 
%         drot_L_mirror(:,i) = -drot_R_temp(:,i);
%         drot_R_mirror(:,i) = -drot_L_temp(:,i);
% 
%         rot_dot_L_mirror(:,i) = -rot_dot_R_temp(:,i);
%         rot_dot_R_mirror(:,i) = -rot_dot_L_temp(:,i);
% 
%         rot_dot_dot_L_mirror(:,i) = -rot_dot_dot_R_temp(:,i);
%         rot_dot_dot_R_mirror(:,i) = -rot_dot_dot_L_temp(:,i);
    end
end

%% recalc rotation axis angles after mirroring (off axes cant be mirrored)
drot_L_mirror = droll_mirror*sind(rot_axis_angle) - dpitch*cosd(rot_axis_angle);
drot_R_mirror = droll_mirror*cosd(rot_axis_angle) + dpitch*sind(rot_axis_angle);

rot_dot_L_mirror = roll_dot_mirror*sind(rot_axis_angle) - pitch_dot*cosd(rot_axis_angle);
rot_dot_R_mirror = roll_dot_mirror*cosd(rot_axis_angle) + pitch_dot*sind(rot_axis_angle);

rot_dot_dot_L_mirror = roll_dot_dot_mirror*sind(rot_axis_angle) - pitch_dot_dot*cosd(rot_axis_angle);
rot_dot_dot_R_mirror = roll_dot_dot_mirror*cosd(rot_axis_angle) + pitch_dot_dot*sind(rot_axis_angle);

%% calc rot1 axis data: initial rotation
drot_axis1normal_mirror = droll_mirror*sind(rot1_axis_angle) - dpitch*cosd(rot1_axis_angle);
drot_axis1_mirror = droll_mirror*cosd(rot1_axis_angle) + dpitch*sind(rot1_axis_angle);

rot_dot_axis1normal_mirror = roll_dot_mirror*sind(rot1_axis_angle) - pitch_dot*cosd(rot1_axis_angle);
rot_dot_axis1_mirror = roll_dot_mirror*cosd(rot1_axis_angle) + pitch_dot*sind(rot1_axis_angle);

rot_dot_dot_axis1normal_mirror = roll_dot_dot_mirror*sind(rot1_axis_angle) - pitch_dot_dot*cosd(rot1_axis_angle);
rot_dot_dot_axis1_mirror = roll_dot_dot_mirror*cosd(rot1_axis_angle) + pitch_dot_dot*sind(rot1_axis_angle);

%% calc rot2 axis data: counter-rotation
drot_axis2normal_mirror = droll_mirror*sind(rot2_axis_angle) - dpitch*cosd(rot2_axis_angle);
drot_axis2_mirror = droll_mirror*cosd(rot2_axis_angle) + dpitch*sind(rot2_axis_angle);

rot_dot_axis2normal_mirror = roll_dot_mirror*sind(rot2_axis_angle) - pitch_dot*cosd(rot2_axis_angle);
rot_dot_axis2_mirror = roll_dot_mirror*cosd(rot2_axis_angle) + pitch_dot*sind(rot2_axis_angle);

rot_dot_dot_axis2normal_mirror = roll_dot_dot_mirror*sind(rot2_axis_angle) - pitch_dot_dot*cosd(rot2_axis_angle);
rot_dot_dot_axis2_mirror = roll_dot_dot_mirror*cosd(rot2_axis_angle) + pitch_dot_dot*sind(rot2_axis_angle);

%% Torques (based on I & C NO nondim, CONVERT ANGLES TO RAD)
Mroll_accel  = Iroll  .* roll_dot_dot_mirror .*pi/180;
Mpitch_accel = Ipitch .* pitch_dot_dot .*pi/180;
Myaw_accel   = Iyaw  .* yaw_dot_dot_mirror .*pi/180;

Mroll_damp  = Croll  .* roll_dot_mirror .*pi/180;
Mpitch_damp = Cpitch .* pitch_dot .*pi/180;
Myaw_damp   = Cyaw  .* yaw_dot_mirror .*pi/180;

Mroll  = Mroll_accel  + Mroll_damp;
Mpitch = Mpitch_accel + Mpitch_damp;
Myaw   = Myaw_accel   + Myaw_damp;

M_L_accel = Mroll_accel*sind(rot_axis_angle) - Mpitch_accel*cosd(rot_axis_angle);
M_R_accel = Mroll_accel*cosd(rot_axis_angle) + Mpitch_accel*sind(rot_axis_angle);

M_L_damp = Mroll_damp*sind(rot_axis_angle) - Mpitch_damp*cosd(rot_axis_angle);
M_R_damp = Mroll_damp*cosd(rot_axis_angle) + Mpitch_damp*sind(rot_axis_angle);

M_L = M_L_accel + M_L_damp;
M_R = M_R_accel + M_R_damp;

% initial rot1
M_axis1normal_accel = Mroll_accel*sind(rot1_axis_angle) - Mpitch_accel*cosd(rot1_axis_angle);
M_axis1_accel = Mroll_accel*cosd(rot1_axis_angle) + Mpitch_accel*sind(rot1_axis_angle);

M_axis1normal_damp = Mroll_damp*sind(rot1_axis_angle) - Mpitch_damp*cosd(rot1_axis_angle);
M_axis1_damp = Mroll_damp*cosd(rot1_axis_angle) + Mpitch_damp*sind(rot1_axis_angle);

M_axis1 = M_axis1_accel + M_axis1_damp;
M_axis1normal = M_axis1normal_accel + M_axis1normal_damp;

% counter rot2
M_axis2normal_accel = Mroll_accel*sind(rot2_axis_angle) - Mpitch_accel*cosd(rot2_axis_angle);
M_axis2_accel = Mroll_accel*cosd(rot2_axis_angle) + Mpitch_accel*sind(rot2_axis_angle);

M_axis2normal_damp = Mroll_damp*sind(rot2_axis_angle) - Mpitch_damp*cosd(rot2_axis_angle);
M_axis2_damp = Mroll_damp*cosd(rot2_axis_angle) + Mpitch_damp*sind(rot2_axis_angle);

M_axis2 = M_axis2_accel + M_axis2_damp;
M_axis2normal = M_axis2normal_accel + M_axis2normal_damp;

% normalized torques
Mroll_accel_norm = Mroll_accel ./ Mg_fly ./ Lwing;
Mpitch_accel_norm = Mpitch_accel ./ Mg_fly ./ Lwing;
Myaw_accel_norm = Myaw_accel ./ Mg_fly ./ Lwing;
M_L_accel_norm = M_L_accel ./ Mg_fly ./ Lwing;
M_R_accel_norm = M_R_accel ./ Mg_fly ./ Lwing;
M_axis1_accel_norm = M_axis1_accel ./ Mg_fly ./ Lwing;
M_axis1normal_accel_norm = M_axis1normal_accel ./ Mg_fly ./ Lwing;
M_axis2_accel_norm = M_axis2_accel ./ Mg_fly ./ Lwing;
M_axis2normal_accel_norm = M_axis2normal_accel ./ Mg_fly ./ Lwing;

Mroll_damp_norm = Mroll_damp ./ Mg_fly ./ Lwing;
Mpitch_damp_norm = Mpitch_damp ./ Mg_fly ./ Lwing;
Myaw_damp_norm = Myaw_damp ./ Mg_fly ./ Lwing;
M_L_damp_norm = M_L_damp ./ Mg_fly ./ Lwing;
M_R_damp_norm = M_R_damp ./ Mg_fly ./ Lwing;
M_axis1_damp_norm = M_axis1_damp ./ Mg_fly ./ Lwing;
M_axis1normal_damp_norm = M_axis1normal_damp ./ Mg_fly ./ Lwing;
M_axis2_damp_norm = M_axis2_damp ./ Mg_fly ./ Lwing;
M_axis2normal_damp_norm = M_axis2normal_damp ./ Mg_fly ./ Lwing;

Mroll_norm = Mroll ./ Mg_fly ./ Lwing;
Mpitch_norm = Mpitch ./ Mg_fly ./ Lwing;
Myaw_norm = Myaw ./ Mg_fly ./ Lwing;
M_L_norm = M_L ./ Mg_fly ./ Lwing;
M_R_norm = M_R ./ Mg_fly ./ Lwing;
M_axis1_norm = M_axis1 ./ Mg_fly ./ Lwing;
M_axis1normal_norm = M_axis1normal ./ Mg_fly ./ Lwing;
M_axis2_norm = M_axis2 ./ Mg_fly ./ Lwing;
M_axis2normal_norm = M_axis2normal ./ Mg_fly ./ Lwing;

%% Torques (based on mean flap freq)
% l_wing = kinDB.wing_length;
% f_mean = mean(f_wb_L);
% 
% Mroll_accel  = nan(size(roll_dot_dot));
% Mpitch_accel = nan(size(roll_dot_dot));
% Myaw_accel   = nan(size(roll_dot_dot));
% 
% Mroll_damp  = nan(size(roll_dot_dot));
% Mpitch_damp = nan(size(roll_dot_dot));
% Myaw_damp   = nan(size(roll_dot_dot));
% 
% for i=1:length(l_wing)
%     Mroll_accel(:,i)  = l_wing(i).^5 .* f_mean.^2 .* Iroll  .* roll_dot_dot_mirror(:,i);
%     Mpitch_accel(:,i) = l_wing(i).^5 .* f_mean.^2 .* Ipitch .* pitch_dot_dot(:,i);
%     Myaw_accel(:,i)   = l_wing(i).^5 .* f_mean.^2 .* Iyaw   .* yaw_dot_dot_mirror(:,i);
% 
%     Mroll_damp(:,i)  = l_wing(i).^5 .* f_mean.^2 .* Croll  .* roll_dot_mirror(:,i);
%     Mpitch_damp(:,i) = l_wing(i).^5 .* f_mean.^2 .* Cpitch .* pitch_dot(:,i);
%     Myaw_damp(:,i)   = l_wing(i).^5 .* f_mean.^2 .* Cyaw   .* yaw_dot_mirror(:,i);
% end
% 
% Mroll  = Mroll_accel  + Mroll_damp;
% Mpitch = Mpitch_accel + Mpitch_damp;
% Myaw   = Myaw_accel   + Myaw_damp;

%% WBmod values
% % Myaw_Norm = Myaw / Myaw_norm;
% % Mpitch_Norm = Mpitch / Mpitch_norm;
% % Mroll_Norm = Mroll / Mroll_norm;
% Fenhance_Norm = (F-1) / Fenhance_norm;
% % 
% % Myaw_WBmean_Norm = Myaw_mean_wb / Myaw_norm;
% % Mpitch_WBmean_Norm = Mpitch_mean_wb / Mpitch_norm;
% % Mroll_WBmean_Norm = Mroll_mean_wb / Mroll_norm;
% F_WBmean_norm = (F_mean_wb-1) / Fenhance_norm;

%% WBmod kinematics
% % RollAccel
% loc = find(abs(DstrokeMOD_wb_RollAccel_bins_meanCIstd(:,1)) == max(abs(DstrokeMOD_wb_RollAccel_bins_meanCIstd(:,1))));
% DstrokeMOD_RollAccel_max = sign(stroke_wb_steady_bins_meanCIstd(loc,1)) * DstrokeMOD_wb_RollAccel_bins_meanCIstd(loc,1);
% 
% loc = find(abs(DpitchMOD_wb_RollAccel_bins_meanCIstd(:,1)) == max(abs(DpitchMOD_wb_RollAccel_bins_meanCIstd(:,1))));
% DpitchMOD_RollAccel_max = sign(pitch_wb_steady_bins_meanCIstd(loc,1)) * DpitchMOD_wb_RollAccel_bins_meanCIstd(loc,1);
% 
% loc = find(abs(DdevMOD_wb_RollAccel_bins_meanCIstd(:,1)) == max(abs(DdevMOD_wb_RollAccel_bins_meanCIstd(:,1))));
% DdevMOD_RollAccel_max = sign(dev_wb_steady_bins_meanCIstd(loc,1)) * DdevMOD_wb_RollAccel_bins_meanCIstd(loc,1);
% 
% % YawAccel
% loc = find(abs(DstrokeMOD_wb_YawAccel_bins_meanCIstd(:,1)) == max(abs(DstrokeMOD_wb_YawAccel_bins_meanCIstd(:,1))));
% DstrokeMOD_YawAccel_max = sign(stroke_wb_steady_bins_meanCIstd(loc,1)) * DstrokeMOD_wb_YawAccel_bins_meanCIstd(loc,1);
% 
% loc = find(abs(DpitchMOD_wb_YawAccel_bins_meanCIstd(:,1)) == max(abs(DpitchMOD_wb_YawAccel_bins_meanCIstd(:,1))));
% DpitchMOD_YawAccel_max = sign(pitch_wb_steady_bins_meanCIstd(loc,1)) * DpitchMOD_wb_YawAccel_bins_meanCIstd(loc,1);
% 
% loc = find(abs(DdevMOD_wb_YawAccel_bins_meanCIstd(:,1)) == max(abs(DdevMOD_wb_YawAccel_bins_meanCIstd(:,1))));
% DdevMOD_YawAccel_max = sign(dev_wb_steady_bins_meanCIstd(loc,1)) * DdevMOD_wb_YawAccel_bins_meanCIstd(loc,1);
% 
% % PitchAccel
% loc = find(abs(strokeMOD_wb_PitchAccel_bins_meanCIstd(:,1)) == max(abs(strokeMOD_wb_PitchAccel_bins_meanCIstd(:,1))));
% strokeMOD_PitchAccel_max = sign(stroke_wb_steady_bins_meanCIstd(loc,1)) * strokeMOD_wb_PitchAccel_bins_meanCIstd(loc,1);
% 
% loc = find(abs(pitchMOD_wb_PitchAccel_bins_meanCIstd(:,1)) == max(abs(pitchMOD_wb_PitchAccel_bins_meanCIstd(:,1))));
% pitchMOD_PitchAccel_max = sign(pitch_wb_steady_bins_meanCIstd(loc,1)) * pitchMOD_wb_PitchAccel_bins_meanCIstd(loc,1);
% 
% loc = find(abs(devMOD_wb_PitchAccel_bins_meanCIstd(:,1)) == max(abs(devMOD_wb_PitchAccel_bins_meanCIstd(:,1))));
% devMOD_PitchAccel_max = sign(dev_wb_steady_bins_meanCIstd(loc,1)) * devMOD_wb_PitchAccel_bins_meanCIstd(loc,1);
% 
% % Fenhance
% loc = find(abs(strokeMOD_wb_Fenhance_bins_meanCIstd(:,1)) == max(abs(strokeMOD_wb_Fenhance_bins_meanCIstd(:,1))));
% strokeMOD_Fenhance_max = sign(stroke_wb_steady_bins_meanCIstd(loc,1)) * strokeMOD_wb_Fenhance_bins_meanCIstd(loc,1);
% 
% loc = find(abs(pitchMOD_wb_Fenhance_bins_meanCIstd(:,1)) == max(abs(pitchMOD_wb_Fenhance_bins_meanCIstd(:,1))));
% pitchMOD_Fenhance_max = sign(pitch_wb_steady_bins_meanCIstd(loc,1)) * pitchMOD_wb_Fenhance_bins_meanCIstd(loc,1);
% 
% loc = find(abs(devMOD_wb_Fenhance_bins_meanCIstd(:,1)) == max(abs(devMOD_wb_Fenhance_bins_meanCIstd(:,1))));
% devMOD_Fenhance_max = sign(dev_wb_steady_bins_meanCIstd(loc,1)) * devMOD_wb_Fenhance_bins_meanCIstd(loc,1);
% 
% freqMOD_Fenhance = freqMOD_wb_Fenhance_meanCIstd(1);

%% calculate variables pre post max etc
% calc_var_pre_post_etc_N

%% calculate nMAX
% % calc_nMAX
% % calc_nMAX_torque
% calc_nMAX_accel
% calc_nMAX_accel_NOwb
calc_nMAX_accel_NOwb_2axes

%% maneuver timing
% % accel based timing
% dt_RollAccel_max_min = (n_rolldot_dotmirrormin - n_rolldot_dotmirrormax)*dt;
% dt_YawAccel_max_min = (n_yawdot_dotmirrormin - n_yawdot_dotmirrormax)*dt;
% dt_PitchAccel_max_min = (n_pitchdot_dotmin - n_pitchdot_dotmax)*dt;
% 
% NwbSteady_RollAccel_max_min = dt_RollAccel_max_min*freq_steady;
% NwbSteady_YawAccel_max_min = dt_YawAccel_max_min*freq_steady;
% NwbSteady_PitchAccel_max_min = dt_PitchAccel_max_min*freq_steady;
% 
% % dt_RollAccel_WBmean_max_min = (n_rolldot_dot_WBmeanmirrormin - n_rolldot_dot_WBmeanmirrormax)*dt;
% % dt_YawAccel_WBmean_max_min = (n_yawdot_dot_WBmeanmirrormin - n_yawdot_dot_WBmeanmirrormax)*dt;
% % dt_PitchAccel_WBmean_max_min = (n_pitchdot_dot_WBmeanmin - n_pitchdot_dot_WBmeanmax)*dt;
% % 
% % NwbSteady_RollAccel_WBmean_max_min = dt_RollAccel_WBmean_max_min*freq_steady;
% % NwbSteady_YawAccel_WBmean_max_min = dt_YawAccel_WBmean_max_min*freq_steady;
% % NwbSteady_PitchAccel_WBmean_max_min = dt_PitchAccel_WBmean_max_min*freq_steady;
% 
% % velocity based timing
% dt_RollVel_max_min = (n_rolldotmirrormin - n_rolldotmirrormax)*dt;
% dt_YawVel_max_min = (n_yawdotmirrormin - n_yawdotmirrormax)*dt;
% dt_PitchVel_max_min = (n_pitchdotmin - n_pitchdotmax)*dt;
% 
% NwbSteady_RollVel_max_min = dt_RollVel_max_min*freq_steady;
% NwbSteady_YawVel_max_min = dt_YawVel_max_min*freq_steady;
% NwbSteady_PitchVel_max_min = dt_PitchVel_max_min*freq_steady;
% 
% % dt_RollVel_WBmean_max_min = (n_rolldot_WBmeanmirrormin - n_rolldot_WBmeanmirrormax)*dt;
% % dt_YawVel_WBmean_max_min = (n_yawdot_WBmeanmirrormin - n_yawdot_WBmeanmirrormax)*dt;
% % dt_PitchVel_WBmean_max_min = (n_pitchdot_WBmeanmin - n_pitchdot_WBmeanmax)*dt;
% % 
% % NwbSteady_RollVel_WBmean_max_min = dt_RollVel_WBmean_max_min*freq_steady;
% % NwbSteady_YawVel_WBmean_max_min = dt_YawVel_WBmean_max_min*freq_steady;
% % NwbSteady_PitchVel_WBmean_max_min = dt_PitchVel_WBmean_max_min*freq_steady;
% 
% % TORQUE based timing
% % dt_Mroll_max_min = (n_Mroll_min - n_Mroll_max)*dt;
% % dt_Myaw_max_min = (n_Myaw_min - n_Myaw_max)*dt;
% % dt_Mpitch_max_min = (n_Mpitch_min - n_Mpitch_max)*dt;
% % 
% % NwbSteady_Roll_max_min = dt_Mroll_max_min*freq_steady;
% % NwbSteady_Yaw_max_min = dt_Myaw_max_min*freq_steady;
% % NwbSteady_Pitch_max_min = dt_Mpitch_max_min*freq_steady;
% % 
% % % dt_Roll_WBmean_max_min = (n_rolldot_WBmeanmirrormin - n_rolldot_WBmeanmirrormax)*dt;
% % % dt_Yaw_WBmean_max_min = (n_yawdot_WBmeanmirrormin - n_yawdot_WBmeanmirrormax)*dt;
% % % dt_Pitch_WBmean_max_min = (n_pitchdot_WBmeanmin - n_pitchdot_WBmeanmax)*dt;
% % % 
% % % NwbSteady_Roll_WBmean_max_min = dt_Roll_WBmean_max_min*freq_steady;
% % % NwbSteady_Yaw_WBmean_max_min = dt_Yaw_WBmean_max_min*freq_steady;
% % % NwbSteady_Pitch_WBmean_max_min = dt_Pitch_WBmean_max_min*freq_steady;

%% cal WBmod kin data
% % cal_WBmod_kin_data
% cal_WBmod_kin_data_torque

%% calc yaw & vel turns !!@!! FIX by using heading pre&post
% calc_turn_vectors
% calc_turn_vectors_headingyaw_pre_post
calc_heading_turn_pre_post_mirror
% calc_heading_turn_pre_post

%% construct flight paths & plot figures & make movies
% calc_flightpath_INCplotNmovie
calc_flightpath_INCplotNmovie_3Dview

%% TIMELINES
if plot_timelines == 1
%% timelines NO skip MIRROR

IDX_plot = IDX;
t_plot = t;

V_plot = V;
dV_plot = dV;

An_hor_plot = An_hor_mirror;
At_hor_plot = At_hor;
A_hor_plot = A_hor;

A_norm_plot = A/g;
Ahor_norm_plot = A_hor_plot/g;
Az_norm_plot = F_ver -1;

F_plot = F;
F_ver_plot = F_ver;
F_hor_plot = F_hor;
Fn_hor_plot = Fn_hor;
Ft_hor_plot = Ft_hor;

stim_angle_vel_plot = stim_angle_vel_mirror;
stim_angle_accel_plot = stim_angle_accel_mirror;
stim_angle_yaw_plot = stim_angle_yaw_mirror;
stim_angle_F_plot = stim_angle_F_mirror;
stim_angle_spn_plot = stim_angle_spn_mirror;

stim_angle_vel_plot = stim_angle_vel_mirror;
stim_angle_accel_plot = stim_angle_accel_mirror;
stim_angle_yaw_plot = stim_angle_yaw_mirror;
stim_angle_F_plot = stim_angle_F_mirror;
stim_angle_spn_plot = stim_angle_spn_mirror;

yaw_plot = yaw_mirror;
pitch_plot = pitch;
roll_plot = roll_mirror;

dyaw_plot = dyaw_mirror;
dpitch_plot = dpitch;
droll_plot = droll_mirror;

yaw_dot_plot = yaw_dot_mirror;
pitch_dot_plot = pitch_dot;
roll_dot_plot = roll_dot_mirror;

yaw_dot_dot_plot = yaw_dot_dot_mirror;
pitch_dot_dot_plot = pitch_dot_dot;
roll_dot_dot_plot = roll_dot_dot_mirror;

% yaw_dot_dot_norm_plot = yaw_dot_dot_norm;
% pitch_dot_dot_norm_plot = pitch_dot_dot_norm;
% roll_dot_dot_norm_plot = roll_dot_dot_norm;

slip_global_plot = slip_global_mirror;
pitch_global_plot = pitch_global;
roll_global_plot = roll_global_mirror;

drot_L_plot = drot_L_mirror;
drot_R_plot = drot_R_mirror;

rot_dot_L_plot = rot_dot_L_mirror;
rot_dot_R_plot = rot_dot_R_mirror;

rot_dot_dot_L_plot = rot_dot_dot_L_mirror;
rot_dot_dot_R_plot = rot_dot_dot_R_mirror;

% Fsp_pitch_plot = Fsp_pitch;
% Fsp_roll_plot = Fsp_roll_mirror;

Fsp_pitch_plot = rad2deg(unwrap(deg2rad(Fsp_pitch)));
Fsp_roll_plot = rad2deg(unwrap(deg2rad(Fsp_roll_mirror)));

F_plot = F;
% F_norm_plot = Fenhance_Norm;

stim_angle_vel_temp = stim_angle_vel_mirror;
stim_angle_accel_temp = stim_angle_accel_mirror;
stim_angle_yaw_temp = stim_angle_yaw_mirror;
stim_angle_F_temp = stim_angle_F_mirror;
stim_angle_spn_temp = stim_angle_spn_mirror;

yaw_temp = yaw_mirror;
pitch_temp = pitch;
roll_temp = roll_mirror;

dyaw_temp = dyaw_mirror;
dpitch_temp = dpitch;
droll_temp = droll_mirror;

slip_global_temp = slip_global_mirror;
pitch_global_temp = pitch_global;
roll_global_temp = roll_global_mirror;

% wrap Fhor dir for <90
stim_angle_accel_plot_wrap = stim_angle_accel_plot;
stim_angle_accel_plot_wrap(stim_angle_accel_plot_wrap<-90) = stim_angle_accel_plot_wrap(stim_angle_accel_plot_wrap<-90) +360;
stim_angle_accel_temp_wrap = stim_angle_accel_plot_wrap;

% remove jumps from plots (+/-180deg)
for i=1:size(stim_angle_vel_plot,2)
    for j=2:size(stim_angle_vel_plot,1)
        if abs(stim_angle_vel_temp(j,i) - stim_angle_vel_temp(j-1,i)) > 90
            stim_angle_vel_plot(j-1,i) = nan;
        end
        if abs(stim_angle_accel_temp(j,i) - stim_angle_accel_temp(j-1,i)) > 90
            stim_angle_accel_plot(j-1,i) = nan;
        end
        if abs(stim_angle_accel_temp_wrap(j,i) - stim_angle_accel_temp_wrap(j-1,i)) > 90
            stim_angle_accel_plot_wrap(j-1,i) = nan;
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

% time window
for i=1:length(t_shift)
        t_window(:,i) = t-t_shift(i);
end

% escape heading max at t=tstop
for j = 1:size(t_window,2)
    for i = 1:size(t_window,1)
        if isnan(stim_angle_vel_plot(i,j))==0 && t_window(i,j)>0 && t_window(i,j)<t_stop
            stim_angle_vel_escape(j,1) = stim_angle_vel_plot(i,j);
            V_escape(j,1) = V(i,j);
        end
    end
end


cmap_plot = cmap_180;
color_var = round(abs(stim_angle_vel_mirror_pre));

%% calc omega around LnR axes
% calc_rotLnR
% calc_rotVelAxis
% calc_rotAxes
% calc_rotAxes_INCtorque
calc_rotAxes_INCtorque_2axes

%% calc cluster means
% cluster_means_5c
% cluster_means_4c

% cluster_means_2c_turn
cluster_means_2c_turn_2axes
% cluster_means_3c_turn

% cluster_means_5cONOFF
% cluster_means_4cONOFF
% cluster_means_3cONOFF

%% save attitude data for statistics
% cd ..
% attitudes_atAmaxNmean_4statistics_turn_dRoll
% attitudes_atAmaxNmean_4statistics_turn_rotAxes
% attitudes_atAmaxNmean_4statistics_turn_rotAxes_plot
% calc_stats_saccade_torques
calc_stats_saccade_torques_2axes

% save stats
save('stats_saccadic_torque.mat','turn_angle_vel_mirror',...
    'M_R_norm_max_startstop','M_L_norm_max_startstop','Myaw_norm_max_startstop',...
    'M_R_norm_min_0stop','M_L_norm_min_0stop',...
    'M_R_norm_mean_startstop','M_L_norm_mean_startstop','Myaw_norm_mean_startstop',...
    'M_axis1_norm_max_startcut','M_axis1normal_norm_max_startcut',...
    'M_axis1_norm_mean_startcut','M_axis1normal_norm_mean_startcut',...
    'M_axis2_norm_min_cutstop','M_axis2normal_norm_min_cutstop',...
    'M_axis2_norm_mean_cutstop','M_axis2normal_norm_mean_cutstop',...
    'stats_M_R_norm_max_startstop','stats_M_L_norm_max_startstop','stats_Myaw_norm_max_startstop',...
    'stats_M_R_norm_min_0stop','stats_M_L_norm_min_0stop',...
    'stats_M_R_norm_mean_startstop','stats_M_L_norm_mean_startstop','stats_Myaw_norm_mean_startstop',...
    'stats_M_axis1_norm_max_startcut','stats_M_axis1normal_norm_max_startcut',...
    'stats_M_axis1_norm_mean_startcut','stats_M_axis1normal_norm_mean_startcut',...
    'stats_M_axis2_norm_min_cutstop','stats_M_axis2normal_norm_min_cutstop',...
    'stats_M_axis2_norm_mean_cutstop','stats_M_axis2normal_norm_mean_cutstop')

%% figures
mkdir('MSfigs_BodyDyn')
cd('MSfigs_BodyDyn')

% stats fig
% colormap(cmap_plot)
% colorbar
saveas(gca,'MSfigs_statistics_saccadic_torque.fig')
saveas(gca,'MSfigs_statistics_saccadic_torque.png')
plot2svg('MSfigs_statistics_saccadic_torque.svg')

% colormap(cmap_plot)
% colorbar
% saveas(gca,'MSfigs_statisticsNcolorbar.fig')
% saveas(gca,'MSfigs_statisticsNcolorbar.png')
% plot2svg('MSfigs_statisticsNcolorbar.svg')

%% fig2: saccade directionality
% % heatmap
% MSfig2_escape_dir_heatmap_notick
% 
% saveas(gca,'MSfig2_saccade_dir_heatmap.fig')
% saveas(gca,'MSfig2_saccade_dir_heatmap.png')
% plot2svg('MSfig2_saccade_dir_heatmap.svg')

% % clusters
% MSfig2_escape_dir_colorcode_MEANclusters_2c_wrap_accel
% % alpha_on_3x3
% 
% saveas(gca,'MSfig2_saccade_dir_colorcode_clusters_quartiles.fig')
% saveas(gca,'MSfig2_saccade_dir_colorcode_clusters_quartiles.png')
% plot2svg('MSfig2_saccade_dir_colorcode_clusters_quartiles.svg')

% NO clusters
MSfig2_escape_dir_colorcode_MEANclusters_1c_wrap_accel_Zup

% alpha_on_3x3

saveas(gca,'MSfig2_saccade_dir_colorcode_NOclusters_ste.fig')
saveas(gca,'MSfig2_saccade_dir_colorcode_NOclusters_ste.png')
plot2svg('MSfig2_saccade_dir_colorcode_NOclusters_ste.svg')

% % 2N1 clusters
% MSfig2_escape_dir_colorcode_MEANclusters_2cN1c_wrap_accel
% % alpha_on_3x3
% 
% saveas(gca,'MSfig2_saccade_dir_colorcode_clusters_quartiles_Fmean.fig')
% saveas(gca,'MSfig2_saccade_dir_colorcode_clusters_quartiles_Fmean.png')
% plot2svg('MSfig2_saccade_dir_colorcode_clusters_quartiles_Fmean.svg')

% % histograms
% MSfig2_hists_windowdata
% 
% saveas(gca,'MSfig2_saccade_data_histograms.fig')
% saveas(gca,'MSfig2_saccade_data_histograms.png')
% plot2svg('MSfig2_saccade_data_histograms.svg')

%% fig3: body dynamics
% % heatmap
% MSfig3_body_dyn_heatmap_notick
% saveas(gca,'MSfig3_body_dyn_heatmap.fig')
% saveas(gca,'MSfig3_body_dyn_heatmap.png')
% plot2svg('MSfig3_body_dyn_heatmap.svg')

% % clusters
% MSfig3_body_dyn_colorcode_MEANclusters_2c_dRollPitch
% % alpha_on_3x3
% 
% saveas(gca,'MSfig3_body_dyn_colorcode_clusters.fig')
% saveas(gca,'MSfig3_body_dyn_colorcode_clusters.png')
% plot2svg('MSfig3_body_dyn_colorcode_clusters.svg')

% NO clusters
MSfig3_body_dyn_colorcode_MEANclusters_1c_dRollPitch
% alpha_on_3x3

saveas(gca,'MSfig3_body_dyn_colorcode_NOclusters.fig')
saveas(gca,'MSfig3_body_dyn_colorcode_NOclusters.png')
plot2svg('MSfig3_body_dyn_colorcode_NOclusters.svg')

%% rot axes

% NO clusters
MSfig3_body_dyn_colorcode_MEANclusters_1c_rotAxesAngle
% alpha_on_3x3

saveas(gca,'MSfig3_body_dyn_colorcode_NOclusters_rotAxes.fig')
saveas(gca,'MSfig3_body_dyn_colorcode_NOclusters_rotAxes.png')
plot2svg('MSfig3_body_dyn_colorcode_NOclusters_rotAxes.svg')

% % 2 clusters
% % MSfig3_body_dyn_colorcode_MEANclusters_2c_rotAxes
% MSfig3_body_dyn_colorcode_MEANclusters_2c_rotAxesAngle
% % alpha_on_3x3
% 
% saveas(gca,'MSfig3_body_dyn_colorcode_clusters_rotAxes.fig')
% saveas(gca,'MSfig3_body_dyn_colorcode_clusters_rotAxes.png')
% plot2svg('MSfig3_body_dyn_colorcode_clusters_rotAxes.svg')
% 
% % quartiles
% MSfig3_body_dyn_colorcode_MEANclusters_2n1c_rotAxesAngle
% % alpha_on_3x3
% 
% saveas(gca,'MSfig3_body_dyn_colorcode_clusters_rotAxesMean.fig')
% saveas(gca,'MSfig3_body_dyn_colorcode_clusters_rotAxesMean.png')
% plot2svg('MSfig3_body_dyn_colorcode_clusters_rotAxesMean.svg')

%% fig4: body dynamics torques

% % clusters
% % MSfig3_bodydyn_torque_colorcode_MEANclusters_2c
% MSfig3_bodydyn_torque_colorcode_MEANclusters_2c_NONdim
% % alpha_on_3x3
% 
% saveas(gca,'MSfig4_bodydyn_torque_colorcode_clusters.fig')
% saveas(gca,'MSfig4_bodydyn_torque_colorcode_clusters.png')
% plot2svg('MSfig4_bodydyn_torque_colorcode_clusters.svg')

% NO clusters
MSfig3_bodydyn_torque_colorcode_MEANclusters_1c_NONdim
% alpha_on_3x3

saveas(gca,'MSfig4_bodydyn_torque_colorcode_NOclusters.fig')
saveas(gca,'MSfig4_bodydyn_torque_colorcode_NOclusters.png')
plot2svg('MSfig4_bodydyn_torque_colorcode_NOclusters.svg')

%% rot axes

% NO clusters
MSfig3_bodydyn_torqueaxes_colorcode_MEANclusters_1c_norm
% alpha_on_3x3

saveas(gca,'MSfig5_bodydyn_TorqueAxes_colorcode_NOclusters.fig')
saveas(gca,'MSfig5_bodydyn_TorqueAxes_colorcode_NOclusters.png')
plot2svg('MSfig5_bodydyn_TorqueAxes_colorcode_NOclusters.svg')

% % 2 clusters
% MSfig3_bodydyn_torqueaxes_colorcode_MEANclusters_2c_norm
% % alpha_on_3x3
% 
% saveas(gca,'MSfig5_bodydyn_TorqueAxes_colorcode_clusters.fig')
% saveas(gca,'MSfig5_bodydyn_TorqueAxes_colorcode_clusters.png')
% plot2svg('MSfig5_bodydyn_TorqueAxes_colorcode_clusters.svg')
% 
% % rot axes
% MSfig3_bodydyn_torqueaxes_colorcode_MEANclusters_2n1c_norm
% % alpha_on_3x3
% 
% saveas(gca,'MSfig5_bodydyn_TorqueAxes_colorcode_clusters_rotAxesMean.fig')
% saveas(gca,'MSfig5_bodydyn_TorqueAxes_colorcode_clusters_rotAxesMean.png')
% plot2svg('MSfig5_bodydyn_TorqueAxes_colorcode_clusters_rotAxesMean.svg')

%% initial rot axes rot1

% NO clusters
MSfig3_bodydyn_torqueaxes1_colorcode_MEANclusters_1c_norm
% alpha_on_3x3

saveas(gca,'MSfig5_bodydyn_TorqueAxes1_colorcode_NOclusters.fig')
saveas(gca,'MSfig5_bodydyn_TorqueAxes1_colorcode_NOclusters.png')
plot2svg('MSfig5_bodydyn_TorqueAxes1_colorcode_NOclusters.svg')

%% initial rot axes rot2

% NO clusters
MSfig3_bodydyn_torqueaxes2_colorcode_MEANclusters_1c_norm
% alpha_on_3x3

saveas(gca,'MSfig5_bodydyn_TorqueAxes2_colorcode_NOclusters.fig')
saveas(gca,'MSfig5_bodydyn_TorqueAxes2_colorcode_NOclusters.png')
plot2svg('MSfig5_bodydyn_TorqueAxes2_colorcode_NOclusters.svg')


%% min&max Torque dyn
% figure
% % MSfig3_TorqueMAXnMIN_vs_turn_TorqueStartStop
% MSfig3_TorqueMAXnMIN_vs_turn_TorqueStartStop_cyanorange
% % alpha_on_3x3
% 
% saveas(gca,'MSfig4_TorqueMINnMAX_VSturn.fig')
% saveas(gca,'MSfig4_TorqueMINnMAX_VSturn.png')
% plot2svg('MSfig4_TorqueMINnMAX_VSturn.svg')

%% min&max Torque dyn L&Raxes
% figure
% % MSfig3_TorqueMAXnMIN_vs_turn_TorqueStartStop_LnRaxes
% MSfig3_TorqueMAXnMIN_vs_turn_TorqueStartStop_LnRaxes_cyanorange
% % alpha_on_3x3
% 
% saveas(gca,'MSfig5_TorqueMINnMAX_VSturn_LnR.fig')
% saveas(gca,'MSfig5_TorqueMINnMAX_VSturn_LnR.png')
% plot2svg('MSfig5_TorqueMINnMAX_VSturn_LnR.svg')

%% min&max&mean Torque dyn L&Raxes
% figure
% % MSfig3_TorqueMAXnMINnMEAN_vs_turn_TorqueStartStop_LnRaxes
% MSfig3_TorqueMAXnMINnMEAN_vs_turn_TorqueStartStop_LnRaxes_CnO
% % alpha_on_3x3
% 
% saveas(gca,'MSfig5_TorqueMINnMAXnMEAN_VSturn_LnR.fig')
% saveas(gca,'MSfig5_TorqueMINnMAXnMEAN_VSturn_LnR.png')
% plot2svg('MSfig5_TorqueMINnMAXnMEAN_VSturn_LnR.svg')

%% min&max&mean Torque dyn axes1
% figure
% MSfig3_TorqueMAXnMINnMEAN_vs_turn_TorqueStartStop_axes1_CnO
% % alpha_on_3x3
% 
% saveas(gca,'MSfig5_TorqueMINnMAXnMEAN_VSturn_axis1.fig')
% saveas(gca,'MSfig5_TorqueMINnMAXnMEAN_VSturn_axis1.png')
% plot2svg('MSfig5_TorqueMINnMAXnMEAN_VSturn_axis1.svg')

%% min&max&mean Torque dyn axes2
% figure
% MSfig3_TorqueMAXnMINnMEAN_vs_turn_TorqueStartStop_axes2_CnO
% % alpha_on_3x3
% 
% saveas(gca,'MSfig5_TorqueMINnMAXnMEAN_VSturn_axis2.fig')
% saveas(gca,'MSfig5_TorqueMINnMAXnMEAN_VSturn_axis2.png')
% plot2svg('MSfig5_TorqueMINnMAXnMEAN_VSturn_axis2.svg')

%% min&max&mean Torque dyn axes1&2
figure
MSfig3_TorqueMAXnMINnMEAN_vs_turn_TorqueStartStop_axes1n2_CnO
% alpha_on_3x3

saveas(gca,'MSfig5_TorqueMINnMAXnMEAN_VSturn_axis1n2.fig')
saveas(gca,'MSfig5_TorqueMINnMAXnMEAN_VSturn_axis1n2.png')
plot2svg('MSfig5_TorqueMINnMAXnMEAN_VSturn_axis1n2.svg')

%% mean body dyn
% figure
% MSfig3_attitudesMEAN_vs_turn_dRoll
% MSfig3_attitudesAtAmax_vs_turn_dRoll
% MSfig3_rotaxisangle_vs_turnangle_atAmax_dRoll
% MSfig3_rotaxisangle_vs_turnangle_MEAN_dRoll
% % alpha_on_3x3
% 
% saveas(gca,'MSfig3_body_dyn_meansNatAmax.fig')
% saveas(gca,'MSfig3_body_dyn_meansNatAmax.png')
% plot2svg('MSfig3_body_dyn_meansNatAmax.svg')
% 
% figure
% MSfig3_attitudesMEAN_vs_turn_rotAxes
% MSfig3_attitudesAtAmax_vs_turn_rotAxes
% MSfig3_rotaxisangle_vs_turnangle_atAmax_dRoll
% MSfig3_rotaxisangle_vs_turnangle_MEAN_dRoll
% % alpha_on_3x3
% 
% saveas(gca,'MSfig3_body_dyn_meansNatAmax_rotAxes.fig')
% saveas(gca,'MSfig3_body_dyn_meansNatAmax_rotAxes.png')
% plot2svg('MSfig3_body_dyn_meansNatAmax_rotAxes.svg')
% 
% % histograms
% MSfig3_hists_windowdata_MEAN
% 
% saveas(gca,'MSfig3_body_dyn_histograms.fig')
% saveas(gca,'MSfig3_body_dyn_histograms.png')
% plot2svg('MSfig3_body_dyn_histograms.svg')

%% fig4: torque & kinematics
% % heatmap
% MSfig4_torques_heatmap_notick
% 
% saveas(gca,'MSfig4_torques_heatmap_notick.fig')
% saveas(gca,'MSfig4_torques_heatmap_notick.png')
% plot2svg('MSfig4_torques_heatmap_notick.svg')
% 
% % subset
% MSfig4_torques_colorcode_subset
% 
% saveas(gca,'MSfig4_torques_colorcode_subset.fig')
% saveas(gca,'MSfig4_torques_colorcode_subset.png')
% plot2svg('MSfig4_torques_colorcode_subset.svg')

end


%% INITIAL&COUNTER rotation axis roseplot rot1 Accel
% % maxHistogramValue = 60;
% plot_rose_rotAxis_rot1n2_accel
% 
% saveas(gca,'MSfig5_rotationaxis_mirror_rose_initialNcounterRotation.fig')
% saveas(gca,'MSfig5_rotationaxis_mirror_rose_initialNcounterRotation.png')
% plot2svg('MSfig5_rotationaxis_mirror_rose_initialNcounterRotation.svg')

%% INITIAL&COUNTER Torque axis roseplot rot1 Accel
% maxHistogramValue = 60;
% plot_rose_TorqueAxis_rot1n2_accel
plot_rose_TorqueAxis_rot1n2_accel_cyanorange

saveas(gca,'MSfig5_TorqueAxis_mirror_rose_initialNcounterTorque.fig')
saveas(gca,'MSfig5_TorqueAxis_mirror_rose_initialNcounterTorque.png')
plot2svg('MSfig5_TorqueAxis_mirror_rose_initialNcounterTorque.svg')

%% rotation axis roseplot INC turn angle color @Amax
% maxHistogramValue = 60;
% plot_rose_rotaxisangle_turnanglecolor_atAmax
% 
% saveas(gca,'MSfig3_rotationaxis_mirror_rose_color_turnAngle_Amax.fig')
% saveas(gca,'MSfig3_rotationaxis_mirror_rose_color_turnAngle_Amax.png')
% plot2svg('MSfig3_rotationaxis_mirror_rose_color_turnAngle_Amax.svg')

%% rotation axis roseplot INC turn angle color MEAN
% maxHistogramValue = 60;
% plot_rose_rotaxisangle_turnanglecolor_MEAN
% 
% saveas(gca,'MSfig3_rotationaxis_mirror_rose_color_turnAngle_MEAN.fig')
% saveas(gca,'MSfig3_rotationaxis_mirror_rose_color_turnAngle_MEAN.png')
% plot2svg('MSfig3_rotationaxis_mirror_rose_color_turnAngle_MEAN.svg')
% 
%% plot pre&post velocity vectors INC turn angle pre color

% pre data
% t_min = t_start - dt;
% t_max = t_start + dt;
% n_window_start = find(t_window>t_min & t_window<t_max);
% V_pre = V(n_window_start);
% Vdir_pre = stim_angle_vel_plot(n_window_start);
V_pre = calc_value(V,n_pre);
Vdir_pre = stim_angle_vel_mirror_pre;

% post data
% V_post = V_escape;
% Vdir_post = stim_angle_vel_escape;
V_post = calc_value(V,n_post);
Vdir_post = stim_angle_vel_mirror_post;

% plotting
% v_pre = [V_pre.*cosd(Vdir_pre);V_pre.*cosd(-Vdir_pre)];
% u_pre = [V_pre.*sind(Vdir_pre);V_pre.*sind(-Vdir_pre)];
v_pre = [V_pre.*cosd(Vdir_pre)];
u_pre = [V_pre.*sind(Vdir_pre)];

V_pre_mean = mean(V_pre);
Vdir_pre_mean = circ_mean_deg_nonan(Vdir_pre);
v_pre_mean = [V_pre_mean.*cosd(Vdir_pre_mean)];
u_pre_mean = [V_pre_mean.*sind(Vdir_pre_mean)];

v_post = V_post .* cosd(Vdir_post);
u_post = V_post .* sind(Vdir_post);

v_post_l = V_post .* cosd(-Vdir_post);
u_post_l = V_post .* sind(-Vdir_post);

V_post_mean = mean(V_post);
Vdir_post_mean = circ_mean_deg_nonan(Vdir_post);
v_post_mean = [V_post_mean.*cosd(Vdir_post_mean)];
u_post_mean = [V_post_mean.*sind(Vdir_post_mean)];

dVdir_post = Vdir_post-Vdir_pre;
dv_post = V_post .* cosd(dVdir_post);
du_post = V_post .* sind(dVdir_post);
dVdir_post_mean = circ_mean_deg_nonan(dVdir_post);
dv_post_mean = V_post_mean .* cosd(dVdir_post_mean);
du_post_mean = V_post_mean .* sind(dVdir_post_mean);

cmap_plot = cmap_180;
color_var = round(abs(turn_angle_vel_mirror));

x_origin = zeros(size(v_pre));
y_origin = zeros(size(v_pre));

% % Vpre blue
% figure
% maxHistogramValue = floor(max(sqrt(u_post.^2 + v_post.^2)));
% maxHistogramValue = (max(sqrt(u_post.^2 + v_post.^2)));
% % maxHistogramValue = .6;
% polar(0, maxHistogramValue,'-k')
% hold on
% 
% for i=1:length(v_pre)
%     plot([x_origin(i) u_pre(i)],[y_origin(i) v_pre(i)],'color',[0,0,1],'linewidth',2)
% end
% compass_zeroup(0, 0.01,'-k');
% 
% % title('heading histogram','fontsize',10) 
% saveas(gca,'MSfig2_heading_vector_mirror_rose_pre_black.fig')
% saveas(gca,'MSfig2_heading_vector_mirror_rose_pre_black.png')
% plot2svg('MSfig2_heading_vector_mirror_rose_pre_black.svg')
% 
% % Vpre color headingpre
% figure
% % maxHistogramValue = floor(max(sqrt(u_post.^2 + v_post.^2)))-1;
% % maxHistogramValue = .6;
% polar(0, maxHistogramValue,'-k')
% hold on
% 
% for i=1:length(v_pre)
%     plot([x_origin(i) u_pre(i)],[y_origin(i) v_pre(i)],'color',cmap_plot(color_var(i),:),'linewidth',2)
% end
% compass_zeroup(0, 0.01,'-k');
% 
% % title('heading histogram','fontsize',10) 
% saveas(gca,'MSfig2_heading_vector_mirror_rose_pre_colorpre.fig')
% saveas(gca,'MSfig2_heading_vector_mirror_rose_pre_colorpre.png')
% plot2svg('MSfig2_heading_vector_mirror_rose_pre_colorpre.svg')
% 
% % Vpost color orange
% figure
% % maxHistogramValue = floor(max(sqrt(u_post.^2 + v_post.^2)))-1;
% % maxHistogramValue = .6;
% polar(0, maxHistogramValue,'-k')
% hold on
% 
% for i=1:length(v_post)
%     plot([x_origin(i) u_post(i)],[y_origin(i) v_post(i)],'color',[1,.5,0],'linewidth',2)
% end
% plot([nanmean(x_origin) u_post_mean],[nanmean(y_origin) v_post_mean],'--k','linewidth',2)
% 
% compass_zeroup(0, 0.01,'-k');
% 
% % title('heading histogram','fontsize',10) 
% saveas(gca,'MSfig2_heading_vector_mirror_rose_post_colorpre.fig')
% saveas(gca,'MSfig2_heading_vector_mirror_rose_post_colorpre.png')
% plot2svg('MSfig2_heading_vector_mirror_rose_post_colorpre.svg')

% % dVpost: dir_post-dir_pre color orange
% figure
% % maxHistogramValue = floor(max(sqrt(u_post.^2 + v_post.^2)))-1;
% % maxHistogramValue = .6;
% polar(0, maxHistogramValue,'-k')
% hold on
% 
% for i=1:length(v_post)
%     plot([x_origin(i) du_post(i)],[y_origin(i) dv_post(i)],'color',[1,0,0],'linewidth',2)
% end
% plot([nanmean(x_origin) du_post_mean],[nanmean(y_origin) dv_post_mean],'--k','linewidth',2)
% plot([nanmean(x_origin) 0],[nanmean(y_origin) V_pre_mean],'--','color',[0,0,1],'linewidth',2)
% 
% compass_zeroup(0, 0.01,'-k');
% 
% % title('heading histogram','fontsize',10) 
% saveas(gca,'MSfig2_Dheading_vector_mirror_rose_preMeanNpost.fig')
% saveas(gca,'MSfig2_Dheading_vector_mirror_rose_preMeanNpost.png')
% plot2svg('MSfig2_Dheading_vector_mirror_rose_preMeanNpost.svg')

% dVpost: dir_post-dir_pre color orange
figure
maxHistogramValue = (max(sqrt(u_post.^2 + v_post.^2)));
% maxHistogramValue = floor(max(sqrt(u_post.^2 + v_post.^2)))-1;
% maxHistogramValue = .6;
polar(0, maxHistogramValue,'-k')
hold on

for i=1:length(v_post)
    plot([x_origin(i) du_post(i)],[y_origin(i) dv_post(i)],'color',[.5,.5,.5],'linewidth',2)
end
plot([nanmean(x_origin) du_post_mean],[nanmean(y_origin) dv_post_mean],'--','color',[1,.5,0],'linewidth',2)
plot([nanmean(x_origin) 0],[nanmean(y_origin) V_pre_mean],'--','color',[0,.5,1],'linewidth',2)

compass_zeroup(0, 0.01,'-k');

% title('heading histogram','fontsize',10) 
saveas(gca,'MSfig2_Dheading_vector_mirror_rose_preMeanNpost2.fig')
saveas(gca,'MSfig2_Dheading_vector_mirror_rose_preMeanNpost2.png')
plot2svg('MSfig2_Dheading_vector_mirror_rose_preMeanNpost2.svg')


% % Vpost color headingpre
% figure
% % maxHistogramValue = floor(max(sqrt(u_post.^2 + v_post.^2)))-1;
% % maxHistogramValue = .6;
% polar(0, maxHistogramValue,'-k')
% hold on
% 
% for i=1:length(v_post)
%     plot([x_origin(i) u_post(i)],[y_origin(i) v_post(i)],'color',cmap_plot(color_var(i),:),'linewidth',2)
% end
% compass_zeroup(0, 0.01,'-k');
% 
% % title('heading histogram','fontsize',10) 
% saveas(gca,'MSfig2_heading_vector_mirror_rose_post_colorpre.fig')
% saveas(gca,'MSfig2_heading_vector_mirror_rose_post_colorpre.png')
% plot2svg('MSfig2_heading_vector_mirror_rose_post_colorpre.svg')
% 

cd ..


