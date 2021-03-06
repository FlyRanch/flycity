clc
clear
% close all

addpath('/home/florian/Dropbox/WORK/toolbox')
addpath('/home/florian/Dropbox/WORK/toolbox/plot2svg')
addpath('/home/florian/Dropbox/WORK/toolbox/CircStat')

% load('flightpathDB_pos_INCq.mat')
load('flightpathDB_pos_qbodyEKF_INCroll_9clusters_2.5n-3.3n3.mat')
% load('flightpathDB_pos_qbodyEKF_NOroll_9clusters_2.5n-3.3n3.mat')
mkdir('response_figs')
cd('response_figs')

plot_timelines = 0
% plot_timelines = 1

cluster_on = 0
% cluster_on = 1

%% settings
linewidth_timelines = 1.5;
skip = 50;

% heatmap resolution
nx = 1000;
ny = 100;

cmap_180 = jet(180);

% polyfit & 95% cof int settings
order = 3;
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
    roll(:,i) = roll_global_first(i) + droll(:,i);
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
Fsp_roll_mirror = Fsp_roll;

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
        Fsp_roll_mirror(:,i) = -Fsp_roll_mirror(:,i);
    end
end

%% calculate variables pre post max etc
% calc_var_pre_post_etc_N

%% calculate nMAX
calc_nMAX

%% calc yaw & vel turns !!@!! FIX by using heading pre&post
% calc_turn_vectors
% calc_turn_vectors_headingyaw_pre_post
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

yaw_dot_plot = yaw_dot;
pitch_dot_plot = pitch_dot;
roll_dot_plot = roll_dot;

Fsp_pitch_plot = Fsp_pitch;
Fsp_roll_plot = Fsp_roll;

stim_angle_vel_temp = stim_angle_vel;
stim_angle_accel_temp = stim_angle_accel;
stim_angle_yaw_temp = stim_angle_yaw;
stim_angle_F_temp = stim_angle_F;
stim_angle_spn_temp = stim_angle_spn;

yaw_temp = yaw;
pitch_temp = pitch;
roll_temp = roll;

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
plot_flightpath_timeline_tshift_headingstart_ATTnDOT

saveas(gca,'flightpaths_ATTnDOT_headingstart.fig')
saveas(gca,'flightpaths_ATTnDOT_headingstart.png')
plot2svg('flightpaths_ATTnDOT_headingstart.svg')

%% heatmap hist norm ATTITUDEnDOT NO MIRROR
plot_flightpath_timeline_histograms_norm_ATTnDOT

subplot(3,2,2)
axis([t_start t_stop -90 90])
set(gca,'YTick',[-90;90],'fontsize',8)
set(gca,'YTickLabel',[90;-90],'fontsize',8)
subplot(3,2,5)
axis([t_start t_stop -2000 2000])
set(gca,'YTick',[-2000;0;2000],'fontsize',8)
set(gca,'YTickLabel',[2000;0;-2000],'fontsize',8)
subplot(3,2,6)
axis([t_start t_stop -90 90])
set(gca,'YTick',[-90;0;90],'fontsize',8)
set(gca,'YTickLabel',[90;0;-90],'fontsize',8)

saveas(gca,'flightpaths_ATTnDOT_hist_norm_tresp.fig')
saveas(gca,'flightpaths_ATTnDOT_hist_norm_tresp.png')
plot2svg('flightpaths_ATTnDOT_hist_norm_tresp.svg')

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
        yaw_plot(:,i) = -yaw_plot(:,i);
        roll_plot(:,i) = -roll_plot(:,i);
        yaw_dot_plot(:,i) = -yaw_dot_plot(:,i);
        roll_dot_plot(:,i) = -roll_dot_plot(:,i);
        Fsp_roll_plot(:,i) = -Fsp_roll_plot(:,i);
    end
end

cmap_plot = cmap_180;
color_var = round(abs(stim_angle_vel_mirror_pre));

%% timelines F&SPN INC MIRROR TURN
plot_flightpath_timeline_tshift_headingstart_FnSP

saveas(gca,'flightpaths_FnSP_headingstart_mirror.fig')
saveas(gca,'flightpaths_FnSP_headingstart_mirror.png')
plot2svg('flightpaths_FnSP_headingstart_mirror.svg')

%% heatmap hist norm F&SPN NO MIRROR
plot_flightpath_timeline_histograms_norm_FnSP

saveas(gca,'flightpaths_FnSP_hist_norm_tresp_mirror.fig')
saveas(gca,'flightpaths_FnSP_hist_norm_tresp_mirror.png')
plot2svg('flightpaths_FnSP_hist_norm_tresp_mirror.svg')

%% timelines ATTITUDEnDOT NO MIRROR TURN
plot_flightpath_timeline_tshift_headingstart_ATTnDOT

subplot(3,2,2)
axis([t_start t_stop -45 135])
subplot(3,2,5)
axis([t_start t_stop -1000 3000])
set(gca,'YTick',[-1000;0;3000],'fontsize',8)
subplot(3,2,6)
axis([t_start t_stop -45 135])

saveas(gca,'flightpaths_ATTnDOT_headingstart_mirror.fig')
saveas(gca,'flightpaths_ATTnDOT_headingstart_mirror.png')
plot2svg('flightpaths_ATTnDOT_headingstart_mirror.svg')

%% heatmap hist norm ATTITUDEnDOT NO MIRROR
plot_flightpath_timeline_histograms_norm_ATTnDOT

subplot(3,2,2)
axis([t_start t_stop -135 45])
set(gca,'YTick',[-90;0],'fontsize',8)
set(gca,'YTickLabel',[90;0],'fontsize',8)
subplot(3,2,5)
axis([t_start t_stop -3000 1000])
set(gca,'YTick',[-3000;0;1000],'fontsize',8)
set(gca,'YTickLabel',[1000;0;-3000],'fontsize',8)
subplot(3,2,6)
axis([t_start t_stop -135 45])
set(gca,'YTick',[-90;0],'fontsize',8)
set(gca,'YTickLabel',[90;0],'fontsize',8)

saveas(gca,'flightpaths_ATTnDOT_hist_norm_tresp_mirror.fig')
saveas(gca,'flightpaths_ATTnDOT_hist_norm_tresp_mirror.png')
plot2svg('flightpaths_ATTnDOT_hist_norm_tresp_mirror.svg')

%% timelines NO skip INC MIRROR TURN
plot_flightpath_timeline_tshift_headingstart_AnormNdir

subplot(3,2,4)
axis([t_start t_stop -10 20])
set(gca,'YTick',[-10;0;20],'fontsize',8)
subplot(3,2,6)
axis([t_start t_stop -10 20])
set(gca,'YTick',[-10;0;20],'fontsize',8)

saveas(gca,'flightpaths_ACCEL_headingstart_mirror.fig')
saveas(gca,'flightpaths_ACCEL_headingstart_mirror.png')
plot2svg('flightpaths_ACCEL_headingstart_mirror.svg')

%% heatmap hist norm INC MIRROR
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
        yaw_plot(:,i) = -yaw_plot(:,i);
        roll_plot(:,i) = -roll_plot(:,i);
        yaw_dot_plot(:,i) = -yaw_dot_plot(:,i);
        roll_dot_plot(:,i) = -roll_dot_plot(:,i);
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
figure
dm=20
dn=20
plotcolor = 'b'
n_now = n_Amax;
angle_pre = stim_angle_vel_pre;
angle_pre_min = -180;
angle_pre_max = 180;

subplot(3,1,1)
angle_post = calc_value(roll,n_now);
plot_angle_pre_post_circmean_doublemirrorsection
% xlabel('initial heading','fontsize',10) 
ylabel('roll','fontsize',10) 
set(gca,'ylim',[-90 90])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
grid on
% axis equal
title('attitudes @ Amax')

subplot(3,1,2)
angle_post = calc_value(pitch,n_now);
plot_angle_pre_post_circmean_mirrorsection
% xlabel('initial heading','fontsize',10) 
ylabel('pitch','fontsize',10) 
set(gca,'ylim',[-90 90])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
grid on
% axis equal

subplot(3,1,3)
angle_post = calc_value(yaw,n_now);
plot_angle_pre_post_circmean_doublemirrorsection
% xlabel('initial heading','fontsize',10) 
ylabel('yaw','fontsize',10) 
set(gca,'ylim',[-90 90])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
grid on
% axis equal

saveas(gca,'attitudesAmax_vs_headingPre.fig')
saveas(gca,'attitudesAmax_vs_headingPre.png')
plot2svg('attitudesAmax_vs_headingPre.svg')

%% attitudes vs headingPre MIRROR@Amax CIRC MEAN
n_now = n_Amax;

figure
subplot(3,1,1)
title('attitudes @ Amax')
hold on

% left turn
plotcolor = 'b'
angle_pre = stim_angle_vel_mirror_pre;
angle_pre_min = -225;
angle_pre_max = 180;

subplot(3,1,1)
angle_post = calc_value(roll_mirror,n_now);
plot_angle_pre_post_circmean_extendedsection
% xlabel('initial heading','fontsize',10) 
ylabel('roll','fontsize',10) 
set(gca,'ylim',[-90 90])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
set(gca,'xlim',[-225 45])
set(gca,'XTick',[-180:90:180],'fontsize',8) 
grid on
% axis equal

subplot(3,1,2)
angle_post = calc_value(dpitch,n_now);
plot_angle_pre_post_circmean_mirrorsection
% xlabel('initial heading','fontsize',10) 
ylabel('pitch','fontsize',10) 
set(gca,'ylim',[-90 90])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
set(gca,'xlim',[-225 45])
set(gca,'XTick',[-180:90:180],'fontsize',8) 
grid on
% axis equal

subplot(3,1,3)
angle_post = calc_value(yaw,n_now);
plot_angle_pre_post_circmean_doublemirrorsection
% xlabel('initial heading','fontsize',10) 
ylabel('yaw','fontsize',10) 
set(gca,'ylim',[-90 90])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
set(gca,'xlim',[-225 45])
set(gca,'XTick',[-180:90:180],'fontsize',8) 
grid on
% axis equal

saveas(gca,'attitudesAmax_vs_headingPre_mirror.fig')
saveas(gca,'attitudesAmax_vs_headingPre_mirror.png')
plot2svg('attitudesAmax_vs_headingPre_mirror.svg')

% right turn
plotcolor = 'r'
angle_pre = -stim_angle_vel_mirror_pre;
angle_pre_min = -180;
angle_pre_max = 180;

subplot(3,1,1)
angle_post = -calc_value(roll_mirror,n_now);
plot_angle_pre_post_circmean_extendedsection
% xlabel('initial heading','fontsize',10) 

subplot(3,1,2)
set(gca,'xlim',[angle_pre_min angle_pre_max])
subplot(3,1,3)
set(gca,'xlim',[angle_pre_min angle_pre_max])


saveas(gca,'attitudesAmax_vs_headingPre_mirror_LnR.fig')
saveas(gca,'attitudesAmax_vs_headingPre_mirror_LnR.png')
plot2svg('attitudesAmax_vs_headingPre_mirror_LnR.svg')

%% attitudesMEAN vs headingPre MIRROR CIRC MEAN

figure
subplot(3,1,1)
title('mean attitudes')
hold on

% left turn
plotcolor = 'b'
angle_pre = stim_angle_vel_mirror_pre;
angle_pre_min = -225;
angle_pre_max = 180;

subplot(3,1,1)
angle_post = calc_circ_mean_value(roll_mirror,n_pre,n_post);
plot_angle_pre_post_circmean_extendedsection
% xlabel('initial heading','fontsize',10) 
ylabel('roll','fontsize',10) 
set(gca,'ylim',[-90 90])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
set(gca,'xlim',[-225 45])
set(gca,'XTick',[-180:90:180],'fontsize',8) 
grid on
% axis equal

subplot(3,1,2)
angle_post = calc_circ_mean_value(pitch,n_pre,n_post);
plot_angle_pre_post_circmean_mirrorsection
% xlabel('initial heading','fontsize',10) 
ylabel('pitch','fontsize',10) 
set(gca,'ylim',[-90 90])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
set(gca,'xlim',[-225 45])
set(gca,'XTick',[-180:90:180],'fontsize',8) 
grid on
% axis equal

subplot(3,1,3)
angle_post = calc_circ_mean_value(yaw,n_pre,n_post);
plot_angle_pre_post_circmean_doublemirrorsection
% xlabel('initial heading','fontsize',10) 
ylabel('yaw','fontsize',10) 
set(gca,'ylim',[-90 90])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
set(gca,'xlim',[-225 45])
set(gca,'XTick',[-180:90:180],'fontsize',8) 
grid on
% axis equal

saveas(gca,'attitudesMEAN_vs_headingPre_mirror.fig')
saveas(gca,'attitudesMEAN_vs_headingPre_mirror.png')
plot2svg('attitudesMEAN_vs_headingPre_mirror.svg')

% right turn
plotcolor = 'r'
angle_pre = -stim_angle_vel_mirror_pre;
angle_pre_min = -180;
angle_pre_max = 180;

subplot(3,1,1)
angle_post = -calc_circ_mean_value(roll_mirror,n_pre,n_post);
plot_angle_pre_post_circmean_extendedsection
% xlabel('initial heading','fontsize',10) 

subplot(3,1,2)
set(gca,'xlim',[angle_pre_min angle_pre_max])
subplot(3,1,3)
set(gca,'xlim',[angle_pre_min angle_pre_max])

saveas(gca,'attitudesMEAN_vs_headingPre_mirror_LnR.fig')
saveas(gca,'attitudesMEAN_vs_headingPre_mirror_LnR.png')
plot2svg('attitudesMEAN_vs_headingPre_mirror_LnR.svg')



%% headingNyaw post vs pre CIRC MEAN
plot_headingNyaw_postVSpre

saveas(gca,'headingNyawpost_vs_headingNyawpre.fig')
saveas(gca,'headingNyawpost_vs_headingNyawpre.png')
plot2svg('headingNyawpost_vs_headingNyawpre.svg')


%% headingNyaw post vs pre AnMIRROR CIRC MEAN
plot_headingNyaw_postVSpre_MIRROR

saveas(gca,'headingNyawpost_vs_headingNyawpre_Anmirror.fig')
saveas(gca,'headingNyawpost_vs_headingNyawpre_Anmirror.png')
plot2svg('headingNyawpost_vs_headingNyawpre_Anmirror.svg')

%% headingNyaw post vs pre AnMIRROR CIRCMEAN leftNright
plot_headingNyaw_postVSpre_MIRROR_left

saveas(gca,'headingNyawpost_vs_headingNyawpre_Anmirror_leftNright.fig')
saveas(gca,'headingNyawpost_vs_headingNyawpre_Anmirror_leftNright.png')
plot2svg('headingNyawpost_vs_headingNyawpre_Anmirror_leftNright.svg')

%% turns vs headingNyaw pre CIRC MEAN
plot_turnVSheadingNyawPRE

saveas(gca,'turns_vs_headingNyawpre.fig')
saveas(gca,'turns_vs_headingNyawpre.png')
plot2svg('turns_vs_headingNyawpre.svg')

%% turns vs headingNyaw pre AnMIRROR CIRC MEAN
plot_turnVSheadingNyawPRE_MIRROR

saveas(gca,'turns_vs_headingNyawpre_AnMIRROR.fig')
saveas(gca,'turns_vs_headingNyawpre_AnMIRROR.png')
plot2svg('turns_vs_headingNyawpre_AnMIRROR.svg')

%% turns vs headingNyaw pre AnMIRROR CIRCMEAN leftNright
plot_turnVSheadingNyawPRE_MIRROR_left

saveas(gca,'turns_vs_headingNyawpre_AnMIRROR_leftNright.fig')
saveas(gca,'turns_vs_headingNyawpre_AnMIRROR_leftNright.png')
plot2svg('turns_vs_headingNyawpre_AnMIRROR_leftNright.svg')

%% AdirNspNmean vs headingNyaw pre CIRC MEAN
plot_AdirSPnMeanVSheadingNyawPRE

saveas(gca,'AdirNspN_vs_headingNyawpre.fig')
saveas(gca,'AdirNspN_vs_headingNyawpre.png')
plot2svg('AdirNspN_vs_headingNyawpre.svg')

%% AdirNspNmean vs headingNyaw pre AnMIRROR CIRC MEAN
plot_AdirSPnMeanVSheadingNyawPRE_mirror

saveas(gca,'AdirNspN_vs_headingNyawpre_mirror.fig')
saveas(gca,'AdirNspN_vs_headingNyawpre_mirror.png')
plot2svg('AdirNspN_vs_headingNyawpre_mirror.svg')

%% AdirNspNmean vs headingNyaw pre AnMIRROR CIRCMEAN leftNright
plot_AdirSPnMeanVSheadingNyawPRE_mirror_left

saveas(gca,'AdirNspN_vs_headingNyawpre_mirror_LnR.fig')
saveas(gca,'AdirNspN_vs_headingNyawpre_mirror_LnR.png')
plot2svg('AdirNspN_vs_headingNyawpre_mirror_LnR.svg')

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
Z = rose(deg2rad(stim_angle_vel_post),36);
for i=1:length(Z)
    set(Z(i),'color','r','linewidth',2)
end
% x = get(Z,'Xdata');
% y = get(Z,'Ydata');
% g=patch(x,y,'r');

hold on
Y = rose(deg2rad(stim_angle_vel_pre),36);
for i=1:length(Y)
    set(Y(i),'color','k','linewidth',2)
end 
% x = get(Y,'Xdata');
% y = get(Y,'Ydata');
% g=patch(x,y,'b');

% title('heading histogram','fontsize',10) 
saveas(gca,'heading_vector_rose_pre_post.fig')
saveas(gca,'heading_rose_pre_post.png')
plot2svg('heading_vector_rose_pre_post.svg')


%% plot pre&post Adir vectors
A_hor_first = calc_value(A_hor,n_first);
A_hor_post = calc_value(A_hor,n_post);
A_hor_mean = calc_mean_value(A_hor,n_pre,n_post);
[A_hor_max,n_Ahormax] = calc_max_value(A_hor,n_pre,n_post);

Adir_first = calc_value(stim_angle_accel,n_first);
Adir_post = calc_value(stim_angle_accel,n_post);
Adir_mean = calc_circ_mean_value(stim_angle_accel,n_pre,n_post);
Adir_Ahormax = calc_value(stim_angle_accel,n_Ahormax);

Ax_first = A_hor_first .* cosd(Adir_first);
Ay_first = A_hor_first .* sind(Adir_first);

Ax_post = A_hor_post .* cosd(Adir_post);
Ay_post = A_hor_post .* sind(Adir_post);

Ax_mean = A_hor_mean .* cosd(Adir_mean);
Ay_mean = A_hor_mean .* sind(Adir_mean);

Ax_max = A_hor_max .* cosd(Adir_Ahormax);
Ay_max = A_hor_max .* sind(Adir_Ahormax);

figure
maxHistogramValue = (max(A_hor_max));
maxHistogramValue = 20;

subplot(2,3,1)
% polar(0, maxHistogramValue,'-k')
% hold on
compass_zeroup(-Ay_first,Ax_first,'k');
title('Astart','fontsize',10) 

subplot(2,3,2)
polar(0, maxHistogramValue,'-k')
hold on
compass_zeroup(-Ay_max,Ax_max,'r');
title('Amax','fontsize',10) 

subplot(2,3,3)
polar(0, maxHistogramValue,'-k')
hold on
compass_zeroup(-Ay_mean,Ax_mean,'b');
title('Amean','fontsize',10) 

maxHistogramValue = 20;

% figure
subplot(2,3,4)
polar(0, maxHistogramValue,'-k')
hold on
Y = rose(deg2rad(Adir_first),36);
% Y = rose(deg2rad(Adir_body_first),36);
for i=1:length(Y)
    set(Y(i),'color','k','linewidth',1)
end
x = get(Y,'Xdata');
y = get(Y,'Ydata');
g=patch(x,y,[.5 .5 .5]);

subplot(2,3,5)
polar(0, maxHistogramValue,'-k')
hold on
X = rose(deg2rad(Adir_Ahormax),36);
for i=1:length(X)
    set(X(i),'color','k','linewidth',1)
end
x = get(X,'Xdata');
y = get(X,'Ydata');
g=patch(x,y,'r');

subplot(2,3,6)
polar(0, maxHistogramValue,'-k')
hold on
X = rose(deg2rad(Adir_mean),36);
for i=1:length(X)
    set(X(i),'color','k','linewidth',1)
end
x = get(X,'Xdata');
y = get(X,'Ydata');
g=patch(x,y,'b');

saveas(gca,'Adir_vector_rose_first_max_mean.fig')
saveas(gca,'Adir_vector_rose_first_max_mean.png')
plot2svg('Adir_vector_rose_first_max_mean.svg')




%% plot data at Amax NO MIRROR
stim_angle_F_now = calc_value(stim_angle_F,n_Amax);
Fsp_roll_now = calc_value(Fsp_roll,n_Amax);
Fsp_pitch_now = calc_value(Fsp_pitch,n_Amax);
roll_now = calc_value(roll,n_Amax);
pitch_now = calc_value(pitch,n_Amax);
yaw_now = calc_value(yaw,n_Amax);
F_now = calc_value(F,n_Amax);
Fn_hor_now = calc_value(Fn_hor,n_Amax);
Ft_hor_now = calc_value(Ft_hor,n_Amax);

figure
maxHistogramValue = 10;
n_bins_rose = 90;
plot_color = 'r';
dh = .1;

plot_rose_hist_data

saveas(gca,'Amax_variables.fig')
saveas(gca,'Amax_variables.png')
plot2svg('Amax_variables.svg')

%% plot data at Amax NO MIRROR
stim_angle_F_now = calc_circ_mean_value(stim_angle_F,n_pre,n_post);
Fsp_roll_now = calc_circ_mean_value(Fsp_roll,n_pre,n_post);
Fsp_pitch_now = calc_circ_mean_value(Fsp_pitch,n_pre,n_post);
roll_now = calc_circ_mean_value(roll,n_pre,n_post);
pitch_now = calc_circ_mean_value(pitch,n_pre,n_post);
yaw_now = calc_circ_mean_value(yaw,n_pre,n_post);
F_now = calc_mean_value(F,n_pre,n_post);
Fn_hor_now = calc_mean_value(Fn_hor,n_pre,n_post);
Ft_hor_now = calc_mean_value(Ft_hor,n_pre,n_post);

figure
maxHistogramValue = 10;
n_bins_rose = 90;
plot_color = 'r';
dh = .1;

plot_rose_hist_data

saveas(gca,'PrePostMean_variables.fig')
saveas(gca,'PrePostMean_variables.png')
plot2svg('PrePostMean_variables.svg')







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