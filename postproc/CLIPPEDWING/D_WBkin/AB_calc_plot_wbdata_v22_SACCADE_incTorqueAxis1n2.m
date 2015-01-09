clc
clear
close all

% % LOOMING
%     loadname = 'kinflightpathDB_pos_qbodyEKF_INCroll_2clusters_Ahor2.75mps2_strokeplane47.5deg_startframe2945.mat';

% % SACCADE (startframe 100, kalman filter issues at start of seq)
% loadname='kinflightpathDB_pos_qbodyEKF_INCroll_9clusters_1.4n-1.7n1.9_strokeplane47.5deg_startframe100.mat'
loadname = 'kinflightpathDB_pos_qbodyEKF_INCroll_strokeplane47.5deg_rotaxis22deg_startframe100.mat';

load(loadname)

const_name=dir('flyVars*')
const_name=const_name.name;
load(const_name)

%% constants

% % Yaw data from Dickinson et al 2010 & converted to Dh scale
% % roll and pitch damping scaled from yaw acordoing to Cheng et al 2009 (table 2)
% 
% % inertial moment coeff !!! [mN mm] !!!
% Iyaw   = 1.87e-12;
% Iroll  = Iyaw;
% Ipitch = 2*Iyaw;
% 
% % damping coeff !!! [mN mm] !!!
% Cyaw   = -1.14e-10;
% Croll  = 0.55 * Iyaw;
% Cpitch = 0.19 * Iyaw;
% 
% % !!! Johan's OLD !!!
% % inertial moment coeff LINEAR NO GYRO INTERACTIONS !!! [mN mm] !!!
% Iroll  = 1.21e-13 *pi/180;
% Ipitch = 1.64e-13 *pi/180;
% Iyaw   = 7.71e-14 *pi/180;
% 
% % damping coeff LINEAR  !!! [mN mm] !!!
% Croll  = 3.47e-12 *pi/180;
% Cpitch = 1.27e-12 *pi/180;
% Cyaw   = 6.75e-12 *pi/180;
% 
% f_ref = 188;
% l_wing_ref = 2.97;   % !!! [mm] !!!
% 
% % inertial moments LINEAR NO GYRO INTERACTIONS !!! [mN mm] !!!
% Iroll  = l_wing_ref^5 * f_ref^2 * 1.21e-13;
% Ipitch = l_wing_ref^5 * f_ref^2 * 1.64e-13;
% Iyaw   = l_wing_ref^5 * f_ref^2 * 7.71e-14;
% 
% % damping coeff LINEAR  !!! [mN mm] !!!
% Croll  = l_wing_ref^5 * f_ref^2 * (-3.49e-12);
% Cpitch = l_wing_ref^5 * f_ref^2 * ( 4.16e-13);
% Cyaw   = l_wing_ref^5 * f_ref^2 * (-2.83e-12);

%% settings
n=200; % bins
fr_max = 50;

fps = settings.fps;
cmap_180 = jet(180);

%% data import
%% wingkin data
n_down_start_L = kinDB.n_down_start_L;
n_up_start_L = kinDB.n_up_start_L;
wingtip_path_L = kinDB.wingtip_path_L;
joint_pos_L = kinDB.joint_pos_L;
stroke_L = kinDB.stroke_L;
pitch_L = kinDB.pitch_L;
dev_L = kinDB.dev_L;

n_down_start_R = kinDB.n_down_start_R;
n_up_start_R = kinDB.n_up_start_R;
wingtip_path_R = kinDB.wingtip_path_R;
joint_pos_R = kinDB.joint_pos_R;
stroke_R = kinDB.stroke_R;
pitch_R = kinDB.pitch_R;
dev_R = kinDB.dev_R;

wing_length = kinDB.wing_length;

% U&aoa@75% wingspan
U_L(:,:) = kinDB.Uwing_L(:,3,:);
U_R(:,:) = kinDB.Uwing_R(:,3,:);
aoa_L(:,:) = kinDB.aoa_L(:,3,:);
aoa_R(:,:) = kinDB.aoa_R(:,3,:);

%% reverse wingkin data with reverse pattern
% stroke_L_temp = stroke_L;
% stroke_R_temp = stroke_R;
% pitch_L_temp = pitch_L;
% pitch_R_temp = pitch_R;
% dev_L_temp = dev_L;
% dev_R_temp = dev_R;
% U_L_temp = U_L;
% U_R_temp = U_R;
% aoa_L_temp = aoa_L;
% aoa_R_temp = aoa_R;
% 
% for i = 1:size(stroke_L,2)
%     if settings.expansion.HorPos(i) == 180
%         stroke_L(:,i) = stroke_R_temp(:,i);
%         stroke_R(:,i) = stroke_L_temp(:,i);
%         pitch_L(:,i) = pitch_R_temp(:,i);
%         pitch_R(:,i) = pitch_L_temp(:,i);
%         dev_L(:,i) = dev_R_temp(:,i);
%         dev_R(:,i) = dev_L_temp(:,i);
%         U_L(:,i) = U_R_temp(:,i);
%         U_R(:,i) = U_L_temp(:,i);
%         aoa_L(:,i) = aoa_R_temp(:,i);
%         aoa_R(:,i) = aoa_L_temp(:,i);
%     end
% end

%% BODY KIN DATA
% rot_axis_angle = settings.rot_axis_angle;
% strokeplane_angle = settings.strokeplane_angle;
strokeplane_angle = strokeplane_WBkin;

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

% angular accels
yaw_dot_dot = nan(size(yaw_dot));
pitch_dot_dot = nan(size(pitch_dot));
roll_dot_dot = nan(size(roll_dot));

yaw_dot_dot(2:end-1,:) = yaw_dot(3:end,:) - yaw_dot(1:end-2,:);
pitch_dot_dot(2:end-1,:) = pitch_dot(3:end,:) - pitch_dot(1:end-2,:);
roll_dot_dot(2:end-1,:) = roll_dot(3:end,:) - roll_dot(1:end-2,:);

yaw_dot_dot = yaw_dot_dot/2/dt;
pitch_dot_dot = pitch_dot_dot/2/dt;
roll_dot_dot = roll_dot_dot/2/dt;

%% n_pre & n_post
n_pre = min([n_turn_start n_accel_start n_decel_start]')';
n_post = max([n_turn_stop n_accel_stop n_decel_stop]')'-1;
for i = 1:length(n_post)
    if isnan(n_post(i)) == 1
        V_temp = V(:,i);
        n_post(i,1) = find(isnan(V_temp)==0, 1, 'last' );
    end
end

%% dV & dAngle & attitudes from omega INC initial conditions
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
    dAngle_accel(:,i) = angle_accel(:,i) - angle_accel_trig2resp(i);
    dAngle_F(:,i) = angle_F(:,i) - angle_F_trig2resp(i);
    dAngle_yaw(:,i) = angle_yaw(:,i) - angle_yaw_trig2resp(i);
    dAngle_spn(:,i) = angle_spn(:,i) - angle_spn_trig2resp(i);
    
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


%% MIRROR TURN based on YAW_MAX
dyaw_post = calc_value(dyaw,n_post);

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
stroke_L_mirror = stroke_L;
stroke_R_mirror = stroke_R;
pitch_L_mirror = pitch_L;
pitch_R_mirror = pitch_R;
dev_L_mirror = dev_L;
dev_R_mirror = dev_R;
U_L_mirror = U_L;
U_R_mirror = U_R;
aoa_L_mirror = aoa_L;
aoa_R_mirror = aoa_R;


for i=1:length(dyaw_post)
    if dyaw_post(i) < 0
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
        stroke_L_mirror(:,i) = stroke_R(:,i);
        stroke_R_mirror(:,i) = stroke_L(:,i);
        pitch_L_mirror(:,i) = pitch_R(:,i);
        pitch_R_mirror(:,i) = pitch_L(:,i);
        dev_L_mirror(:,i) = dev_R(:,i);
        dev_R_mirror(:,i) = dev_L(:,i);
        U_L_mirror(:,i) = U_R(:,i);
        U_R_mirror(:,i) = U_L(:,i);
        aoa_L_mirror(:,i) = aoa_R(:,i);
        aoa_R_mirror(:,i) = aoa_L(:,i);
    end
end

%% recalc rotation axis angles after mirroring (off axes cant be mirrored)
drot_L_mirror = droll_mirror*sind(rot_axis_angle) - dpitch*cosd(rot_axis_angle);
drot_R_mirror = droll_mirror*cosd(rot_axis_angle) + dpitch*sind(rot_axis_angle);

rot_dot_L_mirror = roll_dot_mirror*sind(rot_axis_angle) - pitch_dot*cosd(rot_axis_angle);
rot_dot_R_mirror = roll_dot_mirror*cosd(rot_axis_angle) + pitch_dot*sind(rot_axis_angle);

rot_dot_dot_L_mirror = roll_dot_dot_mirror*sind(rot_axis_angle) - pitch_dot_dot*cosd(rot_axis_angle);
rot_dot_dot_R_mirror = roll_dot_dot_mirror*cosd(rot_axis_angle) + pitch_dot_dot*sind(rot_axis_angle);

%% calc & save wb dataset
% calc_wb_data_subsets_INCrotAxes_InCnorm
calc_wb_data_subsets_INCrotAxes1n2_InCnoNorm

