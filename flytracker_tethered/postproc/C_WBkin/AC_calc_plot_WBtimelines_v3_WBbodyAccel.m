clc
clear
% close all

addpath('/home/florian/Dropbox/WORK/toolbox')
addpath('/home/florian/Dropbox/WORK/toolbox/plot2svg')
addpath('/home/florian/Dropbox/WORK/toolbox/CircStat')
addpath('/home/florian/Dropbox/WORK/toolbox/flytracker')

% loadname=dir('kinflightpathDB*')
loadname=dir('kinflightpathDB_pos_qbodyEKF_INCroll_2clusters_Ahor2.75mps2_strokeplane47.5deg_startframe2945*')
load(loadname.name)


mkdir('MSfigs_BodyDyn')
cd('MSfigs_BodyDyn')

plot_meanWB_on = 0
% plot_meanWB_on = 1
if plot_meanWB_on == 1
    % plot_WB_NOmirror = 1
    % plot_WB_ALL = 1
    % plot_WB_SUBSETS = 1
    % plot_WB_SUBSETS_FILTER = 1

    plot_WB_NOmirror = 0
    plot_WB_ALL = 0
    plot_WB_SUBSETS = 0
    plot_WB_SUBSETS_FILTER = 0
end

% plot_timelines = 0
% cluster_on = 0
plot_timelines = 1
cluster_on = 1

%% subset
% subset_seqs=[7;65;4;13;11]; % old seq set
% subset_seqs=[4;62;1;10;8]; % updated seq set
% subset_seqs = [59;26;48;66;30;92;22]; % new seq set
subset_seqs = [59;1;48;30;66;76;19;22]; % merged sets

seqs_NOmirror = [7;64;21;4;13;11];
%% settings
fps = settings.fps;

linewidth_timelines = 1;
skip = 50;

% heatmap resolution
nx = 1000;
ny = 100;

cmap_180 = jet(180);
cmap_bw = 1-colormap(gray);

grey_color = [.65 .65 .65];

% polyfit & 95% cof int settings
order = 3;
dn=20   % datapoints in bin
dm=20   % bin shift
csaps_filt = .00001;

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
stroke_L_temp = stroke_L;
stroke_R_temp = stroke_R;
pitch_L_temp = pitch_L;
pitch_R_temp = pitch_R;
dev_L_temp = dev_L;
dev_R_temp = dev_R;
U_L_temp = U_L;
U_R_temp = U_R;
aoa_L_temp = aoa_L;
aoa_R_temp = aoa_R;

for i = 1:size(stroke_L,2)
    if settings.expansion.HorPos(i) == 180
        stroke_L(:,i) = stroke_R_temp(:,i);
        stroke_R(:,i) = stroke_L_temp(:,i);
        pitch_L(:,i) = pitch_R_temp(:,i);
        pitch_R(:,i) = pitch_L_temp(:,i);
        dev_L(:,i) = dev_R_temp(:,i);
        dev_R(:,i) = dev_L_temp(:,i);
        U_L(:,i) = U_R_temp(:,i);
        U_R(:,i) = U_L_temp(:,i);
        aoa_L(:,i) = aoa_R_temp(:,i);
        aoa_R(:,i) = aoa_L_temp(:,i);
    end
end

%% BODY KIN DATA
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
    if isnan(n_pre(i)) == 1
        V_temp = V(:,i);
        n_pre(i,1) = find(isnan(V_temp)==0, 1, 'first' );
    end
end

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

%% MIRROR ATTITUDES based on roll_max
roll_dot_pre_mirror = calc_value(roll_dot,n_pre);

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

for i=1:length(roll_dot_pre_mirror)
    if roll_dot_pre_mirror(i) < 0
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

if plot_meanWB_on == 1
    plot_meanWB
end

%% TIMELINES
if plot_timelines == 1
    
% calc wb values
calc_wb_data_allWB
% calc_wb_data

% calc normalized accel
roll_dot_dot_mean_wb_norm = roll_dot_dot_mean_wb./(f_wb_L.^2);
pitch_dot_dot_mean_wb_norm = pitch_dot_dot_mean_wb./(f_wb_L.^2);
yaw_dot_dot_mean_wb_norm = yaw_dot_dot_mean_wb./(f_wb_L.^2);

roll_dot_mean_wb_norm = roll_dot_mean_wb./f_wb_L;
pitch_dot_mean_wb_norm = pitch_dot_mean_wb./f_wb_L;
yaw_dot_mean_wb_norm = yaw_dot_mean_wb./f_wb_L;

% % calculate nMAX
% calc_nMAX_accel

% calc yaw & vel turns !!@!! FIX by using heading pre&post
calc_heading_turn_pre_post

% pre & post headings
stim_angle_vel_mirror_pre = calc_value(stim_angle_vel_mirror,n_pre);
stim_angle_vel_mirror_post = calc_value(stim_angle_vel_mirror,n_post);
stim_angle_yaw_mirror_pre = calc_value(stim_angle_yaw_mirror,n_pre);
stim_angle_yaw_mirror_post = calc_value(stim_angle_yaw_mirror,n_post);

stim_angle_accel_mean = calc_circ_mean_value(stim_angle_accel,n_pre,n_post);
stim_angle_accel_mirror_mean = calc_circ_mean_value(stim_angle_accel_mirror,n_pre,n_post);
stim_angle_spn_mean = calc_circ_mean_value(stim_angle_spn,n_pre,n_post);
stim_angle_spn_mirror_mean = calc_circ_mean_value(stim_angle_spn_mirror,n_pre,n_post);

%% timelines NO skip 

t_pre = calc_t_value(t,n_pre);
t_shift = t_pre;
% t_shift = zeros(size(t_pre));

t_start = -.02;
t_stop = .04;

% frame times
t_hist = [];
for i=1:size(V,2)    
    t_hist = [t_hist;t-t_shift(i)];
end

%% FIG3: WBmean timeline BODY ROLL&PITCH&YAW &VEL&ACCEL NORM headingstart
cmap_plot = cmap_180;
color_var = round(abs(stim_angle_vel_mirror_pre));
MSfig3_body_AttNorm_colorcode_subsetNgrey
saveas(gca,'MSfig3_body_AttNorm_colorcode_subset.fig')
saveas(gca,'MSfig3_body_AttNorm_colorcode_subset.png')
plot2svg('MSfig3_body_AttNorm_colorcode_subset.svg')

%% FIG3: WBmean timeline BODY ROLL&PITCH&YAW ACCEL NORM headingstart
% cmap_plot = cmap_180;
% color_var = round(abs(stim_angle_vel_mirror_pre));
% MSfig3_body_accelNorm_colorcode_subsetNgrey
% saveas(gca,'MSfig3_body_accelNorm_colorcode_subset.fig')
% saveas(gca,'MSfig3_body_accelNorm_colorcode_subset.png')
% plot2svg('MSfig3_body_accelNorm_colorcode_subset.svg')
% 
% MSfig3_body_accelNorm_colorcode_subsetNgrey2
% saveas(gca,'MSfig3_body_accelNorm_colorcode_subset2.fig')
% saveas(gca,'MSfig3_body_accelNorm_colorcode_subset2.png')
% plot2svg('MSfig3_body_accelNorm_colorcode_subset2.svg')

%% timelines BODY ROLL&PITCH&YAW INC ACCEL headingstart
% cmap_plot = cmap_360;
% color_var = round(stim_angle_vel_pre)+180;
% plot_WingKin_timeline_tshift_bodyrollNpitchNyaw_headingstart
% saveas(gca,'WBtimelines_rollNpitchNyaw_headingstart.fig')
% saveas(gca,'WBtimelines_rollNpitchNyaw_headingstart.png')
% plot2svg('WBtimelines_rollNpitchNyaw_headingstart.svg')
% 
% %% HEATMAP BODY ROLL&PITCH&YAW INC ACCEL
% plot_WingKin_timeline_histograms_BodyAtt
% saveas(gca,'WBtimelines_rollNpitchNyaw_heatmap.fig')
% saveas(gca,'WBtimelines_rollNpitchNyaw_heatmap.png')
% plot2svg('WBtimelines_rollNpitchNyaw_heatmap.svg')
% 
% %% wingbeat times
% t_hist_wb = [];
% for i = 1:length(t_shift)
%     t_hist_wb(:,i) = t_wb_L(:,i)-t_shift(i);
% end
% 
% % frame times
% t_hist_fr = [];
% for i=1:size(V,2)    
%     t_hist_fr(:,i) = t-t_shift(i);
% end
% 
% %% timelines DEVIATION headingstart
% % cmap_plot = cmap_360;
% % color_var = round(stim_angle_vel_pre)+180;
% % plot_WingKin_timeline_tshift_DEV
% % saveas(gca,'WBtimelines_deviation_headingstart.fig')
% % saveas(gca,'WBtimelines_deviation_headingstart.png')
% % plot2svg('WBtimelines_deviation_headingstart.svg')
% 
% %% timelines DEVIATION MAXROLLdot
% cmap_plot = jet(200);
% color_var = roll_dot_mirror_max;
% color_var_mean = nanmean(color_var(:));
% color_var_std = nanstd(color_var(:));
% color_var_limit = abs(color_var_mean) + 2*color_var_std;
% color_var = 100*color_var/color_var_limit+100;
% color_var(color_var<1) = 1;
% color_var(color_var>200) = 200;
% color_var = round(color_var);
% 
% plot_WingKin_timeline_tshift_DEV
% saveas(gca,'WBtimelines_deviation_rolldotcolor.fig')
% saveas(gca,'WBtimelines_deviation_rolldotcolor.png')
% plot2svg('WBtimelines_deviation_rolldotcolor.svg')
% 
% %% HEATMAP DEVIATION
% plot_WingKin_timeline_histograms_DEV_interp
% saveas(gca,'WBtimelines_deviation_heatmap.fig')
% saveas(gca,'WBtimelines_deviation_heatmap.png')
% plot2svg('WBtimelines_deviation_heatmap.svg')
% 
% 
% %% timelines WINGSTROKE headingstart
% % cmap_plot = cmap_360;
% % color_var = round(stim_angle_vel_pre)+180;
% % plot_WingKin_timeline_tshift_STROKE
% % saveas(gca,'WBtimelines_wingstroke_headingstart.fig')
% % saveas(gca,'WBtimelines_wingstroke_headingstart.png')
% % plot2svg('WBtimelines_wingstroke_headingstart.svg')
% 
% %% timelines WINGSTROKE MAXFORCE
% cmap_plot = jet(200);
% color_var = F_max-1;
% color_var_mean = nanmean(color_var(:));
% color_var_std = nanstd(color_var(:));
% color_var_limit = abs(color_var_mean) + 1*color_var_std;
% color_var = 100*color_var/color_var_limit+100;
% color_var(color_var<1) = 1;
% color_var(color_var>200) = 200;
% color_var = round(color_var);
% 
% plot_WingKin_timeline_tshift_STROKE
% saveas(gca,'WBtimelines_wingstroke_Fcolor.fig')
% saveas(gca,'WBtimelines_wingstroke_Fcolor.png')
% plot2svg('WBtimelines_wingstroke_Fcolor.svg')
% 
% %% HEATMAP WINGSTROKE
% plot_WingKin_timeline_histograms_STROKE_interp
% saveas(gca,'WBtimelines_wingstroke_heatmap.fig')
% saveas(gca,'WBtimelines_wingstroke_heatmap.png')
% plot2svg('WBtimelines_wingstroke_heatmap.svg')
% 
% 
% %% timelines WINGPITCH headingstart
% % cmap_plot = cmap_360;
% % color_var = round(stim_angle_vel_pre)+180;
% % plot_WingKin_timeline_tshift_WINGPITCH
% % saveas(gca,'WBtimelines_wingpitch_headingstart.fig')
% % saveas(gca,'WBtimelines_wingpitch_headingstart.png')
% % plot2svg('WBtimelines_wingpitch_headingstart.svg')
% 
% %% timelines WINGPITCH MAXROLLDOT
% cmap_plot = jet(200);
% color_var = roll_dot_mirror_max;
% color_var_mean = nanmean(color_var(:));
% color_var_std = nanstd(color_var(:));
% color_var_limit = abs(color_var_mean) + 2*color_var_std;
% color_var = 100*color_var/color_var_limit+100;
% color_var(color_var<1) = 1;
% color_var(color_var>200) = 200;
% color_var = round(color_var);
% 
% plot_WingKin_timeline_tshift_WINGPITCH
% saveas(gca,'WBtimelines_wingpitch_rolldotcolor.fig')
% saveas(gca,'WBtimelines_wingpitch_rolldotcolor.png')
% plot2svg('WBtimelines_wingpitch_rolldotcolor.svg')
% 
% %% HEATMAP WINGPITCH
% plot_WingKin_timeline_histograms_WINGPITCH_interp
% saveas(gca,'WBtimelines_wingpitch_heatmap.fig')
% saveas(gca,'WBtimelines_wingpitch_heatmap.png')
% plot2svg('WBtimelines_wingpitch_heatmap.svg')
% 
% %% timelines FLAPFREQ headingstart
% % cmap_plot = cmap_360;
% % color_var = round(stim_angle_vel_pre)+180;
% % plot_WingKin_timeline_tshift_FLAPFREQ
% % saveas(gca,'WBtimelines_flapfreq_headingstart.fig')
% % saveas(gca,'WBtimelines_flapfreq_headingstart.png')
% % plot2svg('WBtimelines_flapfreq_headingstart.svg')
% 
% %% timelines FLAPFREQ MAXROLLDOT
% cmap_plot = jet(200);
% color_var = roll_dot_mirror_max;
% color_var_mean = nanmean(color_var(:));
% color_var_std = nanstd(color_var(:));
% color_var_limit = abs(color_var_mean) + 2*color_var_std;
% color_var = 100*color_var/color_var_limit+100;
% color_var(color_var<1) = 1;
% color_var(color_var>200) = 200;
% color_var = round(color_var);
% 
% plot_WingKin_timeline_tshift_FLAPFREQ
% saveas(gca,'WBtimelines_flapfreq_rolldotcolor.fig')
% saveas(gca,'WBtimelines_flapfreq_rolldotcolor.png')
% plot2svg('WBtimelines_flapfreq_rolldotcolor.svg')
% 
% %% HEATMAP FLAPFREQ
% plot_WingKin_timeline_histograms_FLAPFREQ_interp
% saveas(gca,'WBtimelines_flapfreq_heatmap.fig')
% saveas(gca,'WBtimelines_flapfreq_heatmap.png')
% plot2svg('WBtimelines_flapfreq_heatmap.svg')
cd ..
end