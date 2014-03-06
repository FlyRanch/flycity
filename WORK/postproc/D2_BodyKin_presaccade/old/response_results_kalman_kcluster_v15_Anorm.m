clc
clear
close all

addpath('/home/florian/Dropbox/WORK/toolbox')
addpath('/home/florian/Dropbox/WORK/toolbox/plot2svg')
addpath('/home/florian/Dropbox/WORK/toolbox/CircStat')

% load('flightpathDB_pos_INCq.mat')
load('flightpathDB_pos_qbodyEKF_INCroll_9clusters_2.5n-3.3n3_response.mat')
% load('flightpathDB_pos_qbodyEKF_NOroll_9clusters_2.5n-3.3n3_response.mat')
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

% polyfit & 95% cof int settings
order = 3;
dn=20   % datapoints in bin
dm=5   % bin shift
csaps_filt = .00001;

%% data import
IDX = pathDB.IDX;
cmap_k = settings.cmap_k;
cmap_360 = settings.cmap_360;

t = pathDB.t;
V = pathDB.V;
An_hor = pathDB.An_hor;
An_hor_abs = abs(pathDB.An_hor);
At_hor = pathDB.At_hor;
A_hor = pathDB.A_hor;

stim_angle_vel = pathDB.stim_angle_vel;
stim_angle_accel = pathDB.stim_angle_accel;
stim_angle_yaw = pathDB.stim_angle_yaw;

accel_angle_hor_vel = pathDB.accel_angle_hor_vel;
accel_angle_hor_body = pathDB.accel_angle_hor_body;

slip = pathDB.slip;
pitch = pathDB.pitch;
roll = pathDB.roll;
% roll = pathDB.roll_global;

teta = patternDB.teta;

t_Ahormax = responseDB.t_Ahormax;
V_Ahormax = responseDB.V_Ahormax;
An_hor_Ahormax = responseDB.An_hor_Ahormax;
At_hor_Ahormax = responseDB.At_hor_Ahormax;

stim_angle_vel_Ahormax = responseDB.stim_angle_vel_Ahormax;
stim_angle_accel_Ahormax = responseDB.stim_angle_accel_Ahormax;
stim_angle_yaw_Ahormax = responseDB.stim_angle_yaw_Ahormax;

accel_angle_hor_vel_Ahormax = responseDB.accel_angle_hor_vel_Ahormax;
accel_angle_hor_body_Ahormax = responseDB.accel_angle_hor_body_Ahormax;

slip_Ahormax = responseDB.slip_Ahormax;
pitch_Ahormax = responseDB.pitch_Ahormax;
roll_Ahormax = responseDB.roll_Ahormax;

t_Amax = responseDB.t_Amax;
V_Amax = responseDB.V_Amax;
An_hor_Amax = responseDB.An_hor_Amax;
At_hor_Amax = responseDB.At_hor_Amax;

stim_angle_vel_Amax = responseDB.stim_angle_vel_Amax;
stim_angle_accel_Amax = responseDB.stim_angle_accel_Amax;
stim_angle_yaw_Amax = responseDB.stim_angle_yaw_Amax;

accel_angle_hor_vel_Amax = responseDB.accel_angle_hor_vel_Amax;
accel_angle_hor_body_Amax = responseDB.accel_angle_hor_body_Amax;

slip_Amax = responseDB.slip_Amax;
pitch_Amax = responseDB.pitch_Amax;
roll_Amax = responseDB.roll_Amax;

stim_angle_vel_pre_resp = responseDB.stim_angle_vel_pre_resp;
stim_angle_vel_pre_turn = responseDB.stim_angle_vel_pre_turn;

stim_angle_vel_post_resp = responseDB.stim_angle_vel_post_resp;
stim_angle_vel_post_turn = responseDB.stim_angle_vel_post_turn;

stim_angle_accel_pre_resp = responseDB.stim_angle_accel_pre_resp;
stim_angle_accel_pre_turn = responseDB.stim_angle_accel_pre_turn;

stim_angle_accel_post_resp = responseDB.stim_angle_accel_post_resp;
stim_angle_accel_post_turn = responseDB.stim_angle_accel_post_turn;

stim_angle_accel_resp_mean = responseDB.stim_angle_accel_resp_mean;
stim_angle_accel_resp_std = responseDB.stim_angle_accel_resp_std;
stim_angle_accel_resp_min = responseDB.stim_angle_accel_resp_min;
stim_angle_accel_resp_max = responseDB.stim_angle_accel_resp_max;

stim_angle_yaw_pre_resp = responseDB.stim_angle_yaw_pre_resp;
stim_angle_yaw_pre_turn = responseDB.stim_angle_yaw_pre_turn;

stim_angle_yaw_post_resp = responseDB.stim_angle_yaw_post_resp;
stim_angle_yaw_post_turn = responseDB.stim_angle_yaw_post_turn;

accel_angle_hor_vel_pre_resp = responseDB.accel_angle_hor_vel_pre_resp;
accel_angle_hor_vel_pre_turn = responseDB.accel_angle_hor_vel_pre_turn;
accel_angle_hor_vel_pre_accel = responseDB.accel_angle_hor_vel_pre_accel;
        
accel_angle_hor_vel_post_turn = responseDB.accel_angle_hor_vel_post_turn;
accel_angle_hor_vel_post_resp = responseDB.accel_angle_hor_vel_post_resp;
accel_angle_hor_vel_post_accel = responseDB.accel_angle_hor_vel_post_accel;
        
accel_angle_hor_body_pre_resp = responseDB.accel_angle_hor_body_pre_resp;
accel_angle_hor_body_pre_turn = responseDB.accel_angle_hor_body_pre_turn;
accel_angle_hor_body_pre_accel = responseDB.accel_angle_hor_body_pre_accel;
        
accel_angle_hor_body_post_turn = responseDB.accel_angle_hor_body_post_turn;
accel_angle_hor_body_post_resp = responseDB.accel_angle_hor_body_post_resp;
accel_angle_hor_body_post_accel = responseDB.accel_angle_hor_body_post_accel;

accel_angle_hor_vel_resp_mean = responseDB.accel_angle_hor_vel_resp_mean;
accel_angle_hor_body_resp_mean = responseDB.accel_angle_hor_body_resp_mean;

slip_pre_resp = responseDB.slip_pre_resp;
slip_pre_turn = responseDB.slip_pre_turn;

slip_post_resp = responseDB.slip_post_resp;
slip_post_turn = responseDB.slip_post_turn;

pitch_pre_resp = responseDB.pitch_pre_resp;
pitch_pre_turn = responseDB.pitch_pre_turn;

pitch_post_resp = responseDB.pitch_post_resp;
pitch_post_turn = responseDB.pitch_post_turn;

roll_pre_resp = responseDB.roll_pre_resp;
roll_pre_turn = responseDB.roll_pre_turn;

roll_post_resp = responseDB.roll_post_resp;
roll_post_turn = responseDB.roll_post_turn;

V_trig2resp = responseDB.V_trig2resp;
V_pre_resp = responseDB.V_pre_resp;
V_post_resp = responseDB.V_post_resp;

V_pre_accel = responseDB.V_pre_accel;
V_post_accel = responseDB.V_post_accel;

An_hor_max = responseDB.An_hor_max;
At_hor_max = responseDB.At_hor_max;
At_hor_min = responseDB.At_hor_min;

A_hor_pre_resp = responseDB.A_hor_pre_resp;
A_hor_pre_turn = responseDB.A_hor_pre_turn;
A_hor_pre_accel = responseDB.A_hor_pre_accel;

A_hor_post_resp = responseDB.A_hor_post_resp;
A_hor_post_turn = responseDB.A_hor_post_turn;
A_hor_post_accel = responseDB.A_hor_post_accel;

A_hor_max = responseDB.A_hor_max;
A_hor_mean = responseDB.A_hor_mean;

A_max = responseDB.A_max;
A_mean = responseDB.A_mean;

teta_max = settings.expansion.maxangle;
OFF = settings.expansion.OFF;

% t_pre = responseDB.t_resp;
% t_post = responseDB.t_turn_stop;

clear dV
for i = 1:length(V_trig2resp)
    dV(:,i) = V(:,i) - V_trig2resp(i);
end

%% calc An_max At_max At_min
calc_start_stop_data

%% calc_heading_pre_post
calc_heading_pre_post

%% calc yaw & vel turns !!@!! FIX by using heading pre&post
% calc_turn_vectors
calc_turn_vectors_headingyaw_pre_post

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

stim_angle_vel_plot = stim_angle_vel;
stim_angle_accel_plot = stim_angle_accel;
stim_angle_yaw_plot = stim_angle_yaw;
slip_plot = slip;
pitch_plot = pitch;
roll_plot = roll;

stim_angle_vel_temp = stim_angle_vel;
stim_angle_accel_temp = stim_angle_accel;
stim_angle_yaw_temp = stim_angle_yaw;
slip_temp = slip;
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
        if abs(slip_temp(j,i) - slip_temp(j-1,i)) > 90
            slip_plot(j-1,i) = nan;
        end
        if abs(pitch_temp(j,i) - pitch_temp(j-1,i)) > 90
            pitch_plot(j-1,i) = nan;
        end
        if abs(roll_temp(j,i) - roll_temp(j-1,i)) > 90
            roll_plot(j-1,i) = nan;
        end
    end
end

%% timelines NO skip NO MIRROR TURN

t_shift = responseDB.t_resp;
cmap_plot = cmap_360;
color_var = round(heading_pre)+180;

% plot_flightpath_timeline_tshift_headingstart_roll
% plot_flightpath_timeline_tshift_headingstart_accel
% plot_flightpath_timeline_tshift_headingstart_pitch
plot_flightpath_timeline_tshift_headingstart_Anorm

saveas(gca,'flightpaths_headingstart.fig')
saveas(gca,'flightpaths_headingstart.png')
plot2svg('flightpaths_headingstart.svg')

%% heatmap hist norm NO MIRROR

t_shift = responseDB.t_resp;
cmap_plot = cmap_360;
color_var = round(heading_pre)+180;

t_hist = [];
for i=1:size(V,2)    
    t_hist = [t_hist;t-t_shift(i)];
end
% plot_flightpath_timeline_histograms_norm_roll
% plot_flightpath_timeline_histograms_norm_accel
% plot_flightpath_timeline_histograms_norm_pitch
plot_flightpath_timeline_histograms_norm_Anorm

saveas(gca,'flightpath_hist_norm_tresp.fig')
saveas(gca,'flightpath_hist_norm_tresp.png')
plot2svg('flightpath_hist_norm_tresp.svg')

%% MIRROR TURN
for i=1:size(stim_angle_vel_plot,2)
    if An_hor_max(i) < 0
        stim_angle_vel_plot(:,i) = -stim_angle_vel_plot(:,i);
        stim_angle_accel_plot(:,i) = -stim_angle_accel_plot(:,i);
        stim_angle_yaw_plot(:,i) = -stim_angle_yaw_plot(:,i);
        slip_plot(:,i) = -slip_plot(:,i);
        roll_plot(:,i) = -roll_plot(:,i);
        An_hor_plot(:,i) = -An_hor_plot(:,i);
    end
end

%% timelines NO skip INC MIRROR TURN
t_shift = responseDB.t_resp;
cmap_plot = cmap_180;
color_var = abs(round(heading_pre));

% plot_flightpath_timeline_tshift_headingstart_roll
% plot_flightpath_timeline_tshift_headingstart_accel
% plot_flightpath_timeline_tshift_headingstart_pitch
plot_flightpath_timeline_tshift_headingstart_Anorm

subplot(3,3,2)
axis([-.025 .05 -45 180])
set(gca,'YTick',[-45;0;90;180],'fontsize',8) 

subplot(3,3,6)
axis([-.025 .05 -10 20])
set(gca,'YTick',[-10:10:20],'fontsize',8) 

saveas(gca,'flightpaths_headingstart_mirror.fig')
saveas(gca,'flightpaths_headingstart_mirror.png')
plot2svg('flightpaths_headingstart_mirror.svg')

%% heatmap hist norm INC MIRROR
% t_resp timeshift
t_shift = responseDB.t_resp;
cmap_plot = cmap_180;
color_var = abs(round(heading_pre));

t_hist = [];
for i=1:size(V,2)    
    t_hist = [t_hist;t-t_shift(i)];
end
% plot_flightpath_timeline_histograms_norm_roll
% plot_flightpath_timeline_histograms_norm_accel
% plot_flightpath_timeline_histograms_norm_pitch
plot_flightpath_timeline_histograms_norm_Anorm

subplot(3,3,2)
axis([-.025 .05 -180 45])
set(gca,'YTick',[-180;-90;0;45],'YTicklabel',[180;90;0;-45],'fontsize',8) 

subplot(3,3,6)
axis([-.025 .05 -20 10])
set(gca,'YTick',[-20:10:10],'YTicklabel',[20:-10:-10],'fontsize',8) 

saveas(gca,'flightpath_hist_norm_tresp_mirror.fig')
saveas(gca,'flightpath_hist_norm_tresp_mirror.png')
plot2svg('flightpath_hist_norm_tresp_mirror.svg')


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
slip_plot = slip(1:skip:end,:);
pitch_plot = pitch(1:skip:end,:);
roll_plot = roll(1:skip:end,:);

stim_angle_vel_temp = stim_angle_vel(1:skip:end,:);
stim_angle_accel_temp = stim_angle_accel(1:skip:end,:);
stim_angle_yaw_temp = stim_angle_yaw(1:skip:end,:);
slip_temp = slip(1:skip:end,:);
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
        if abs(slip_temp(j,i) - slip_temp(j-1,i)) > 90
            slip_plot(j-1,i) = nan;
        end
        if abs(pitch_temp(j,i) - pitch_temp(j-1,i)) > 90
            pitch_plot(j-1,i) = nan;
        end
        if abs(roll_temp(j,i) - roll_temp(j-1,i)) > 90
            roll_plot(j-1,i) = nan;
        end
    end
end

%% cluster plots
if cluster_on == 1
%% clusters NO MIRROR
t_shift = responseDB.t_resp;

% plot_flightpath_timeline_tshift_cluster_slip_roll
% plot_flightpath_timeline_tshift_cluster_slip_accel
% plot_flightpath_timeline_tshift_cluster_slip_pitch
plot_flightpath_timeline_tshift_cluster_slip_Anorm

saveas(gca,'flightpaths_clusters_tresp.fig')
saveas(gca,'flightpaths_clusters_tresp.png')
plot2svg('flightpaths_clusters_tresp.svg')


%% MIRROR TURN
for i=1:size(stim_angle_vel_plot,2)
    if An_hor_max(i) < 0
        stim_angle_vel_plot(:,i) = -stim_angle_vel_plot(:,i);
        stim_angle_accel_plot(:,i) = -stim_angle_accel_plot(:,i);
        stim_angle_yaw_plot(:,i) = -stim_angle_yaw_plot(:,i);
        slip_plot(:,i) = -slip_plot(:,i);
        roll_plot(:,i) = -roll_plot(:,i);
        An_hor_plot(:,i) = -An_hor_plot(:,i);
    end
end

%% clusters INC MIRROR
t_shift = responseDB.t_resp;

% plot_flightpath_timeline_tshift_cluster_slip_roll
% plot_flightpath_timeline_tshift_cluster_slip_accel
% plot_flightpath_timeline_tshift_cluster_slip_pitch
plot_flightpath_timeline_tshift_cluster_slip_Anorm

subplot(3,3,2)
axis([-.025 .05 -45 180])
set(gca,'YTick',[-45;0;90;180],'fontsize',8) 

subplot(3,3,6)
axis([-.025 .05 -10 20])
set(gca,'YTick',[-10:10:20],'fontsize',8) 

saveas(gca,'flightpaths_clusters_tresp_mirror.fig')
saveas(gca,'flightpaths_clusters_tresp_mirror.png')
plot2svg('flightpaths_clusters_tresp_mirror.svg')

end
end

%% all An vs At, heading_pre colorcoded
% figure
% hold on
% 
% color_code = round(heading_pre)+180;
% 
% for i = 1:size(An_hor,2)
%     plot(An_hor(:,i),At_hor(:,i),'-','color',cmap_360(color_code(i),:),'linewidth',linewidth_timelines)
% end
%     
% figure
% hold on
% 
% color_code = ceil(abs(heading_pre));
% 
% for i = 1:size(An_hor_abs,2)
%     plot(An_hor_abs(:,i),At_hor(:,i),'-','color',cmap_180(color_code(i),:),'linewidth',linewidth_timelines)
% end

%% At_max & At_min vs heading_pre

An_hor_max_all = abs(An_hor_max);
At_hor_max_all = At_hor_max;
At_hor_min_all = At_hor_min;

for i = 1:length(An_hor_max_all)
    if isnan(An_hor_max_all(i)) == 1 && isnan(n_pre(i)) == 0
        An_hor_max_all(i) = max(An_hor(n_pre(i):n_post(i),i));
    end
    
    if isnan(At_hor_max_all(i)) == 1 && isnan(n_pre(i)) == 0
        At_hor_max_all(i) = max(At_hor(n_pre(i):n_post(i),i));
    end
    
    if isnan(At_hor_min_all(i)) == 1 && isnan(n_pre(i)) == 0
        At_hor_min_all(i) = min(At_hor(n_pre(i):n_post(i),i));
    end
end
    
figure
plot(heading_pre,An_hor_max_all,'ok','markerfacecolor','y','markersize',5)
hold on
plot(heading_pre,At_hor_max_all,'ok','markerfacecolor','r','markersize',5)
plot(heading_pre,At_hor_min_all,'ok','markerfacecolor','b','markersize',5)
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

%% heatmap hist attitudes
figure

binx = -180:5:180;
biny = -180:5:180;

%% heatmap hist roll vs slip
x_hist = slip(:);
y_hist = roll(:);

subplot(2,2,1)
yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);
imagesc(binx,-biny,yx_hist_log)

% axis equal
axis([-180 180 -180 180])
set(gca,'XTick',-180:90:180) 
set(gca,'YTick',-180:90:180,'YTicklabel',180:-90:-180,'fontsize',8)
xlabel('slip','fontsize',10) 
ylabel('roll','fontsize',10) 
% colorbar

%% heatmap hist roll vs pitch
x_hist = pitch(:);
y_hist = roll(:);

subplot(2,2,2)
yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);
imagesc(binx,-biny,yx_hist_log)

% axis equal
axis([-45 90 -180 180])
set(gca,'XTick',-180:45:180) 
set(gca,'YTick',-180:90:180,'YTicklabel',180:-90:-180,'fontsize',8)
xlabel('pitch','fontsize',10) 
ylabel('roll','fontsize',10) 
% colorbar

%% heatmap hist pitch vs slip
x_hist = pitch(:);
y_hist = slip(:);

subplot(2,2,4)
yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);
imagesc(binx,-biny,yx_hist_log)

% axis equal
axis([-45 90 -180 180])
set(gca,'XTick',-180:45:180) 
set(gca,'YTick',-180:90:180,'YTicklabel',180:-90:-180,'fontsize',8)
xlabel('pitch','fontsize',10) 
ylabel('slip','fontsize',10) 
% colorbar

saveas(gca,'attitudes_hist_all.fig')
saveas(gca,'attitudes_hist_all.png')
plot2svg('attitudes_hist_all.svg')

%% A vs attitude heatmaps
figure

plot_accel_vs_attitude_subplots

saveas(gca,'accels_vs_attitudes_hist.fig')
saveas(gca,'accels_vs_attitudes_hist.png')
plot2svg('accels_vs_attitudes_hist.svg')

%% A vs attitude heatmaps NORM
figure

plot_accel_vs_attitude_subplots_norm

saveas(gca,'accels_vs_attitudes_hist_NORM.fig')
saveas(gca,'accels_vs_attitudes_hist_NORM.png')
plot2svg('accels_vs_attitudes_hist_NORM.svg')

%% plot turn vs accel and decel reaction time
t_Anresp = responseDB.t_turn_start;

figure
hold on

plotcolor = 'b';
t_Atresp = responseDB.t_decel_start;
plot(t_Anresp,t_Atresp,'ok','MarkerFaceColor',plotcolor,'MarkerSize',5)

plotcolor = 'r';
t_Atresp = responseDB.t_accel_start;
plot(t_Anresp,t_Atresp,'ok','MarkerFaceColor',plotcolor,'MarkerSize',5)

plot([0,1],[0,1],'--k')
legend('decel','accel')

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


%% headingNyaw post vs pre CIRC MEAN
figure

plot_headingNyaw_postVSpre

saveas(gca,'headingNyawpost_vs_headingNyawpre.fig')
saveas(gca,'headingNyawpost_vs_headingNyawpre.png')
plot2svg('headingNyawpost_vs_headingNyawpre.svg')


%% headingNyaw post vs pre AnMIRROR CIRC MEAN
figure

plot_headingNyaw_postVSpre_MIRROR

saveas(gca,'headingNyawpost_vs_headingNyawpre_Anmirror.fig')
saveas(gca,'headingNyawpost_vs_headingNyawpre_Anmirror.png')
plot2svg('headingNyawpost_vs_headingNyawpre_Anmirror.svg')

%% headingNyaw post vs pre AnMIRROR CIRCMEAN leftNright
figure

plot_headingNyaw_postVSpre_MIRRORleftNright

saveas(gca,'headingNyawpost_vs_headingNyawpre_Anmirror_leftNright.fig')
saveas(gca,'headingNyawpost_vs_headingNyawpre_Anmirror_leftNright.png')
plot2svg('headingNyawpost_vs_headingNyawpre_Anmirror_leftNright.svg')

%% turns vs headingNyaw pre CIRC MEAN
figure

plot_turnVSheadingNyawPRE

saveas(gca,'turns_vs_headingNyawpre.fig')
saveas(gca,'turns_vs_headingNyawpre.png')
plot2svg('turns_vs_headingNyawpre.svg')

%% turns vs headingNyaw pre AnMIRROR CIRC MEAN
figure

plot_turnVSheadingNyawPRE_MIRROR

saveas(gca,'turns_vs_headingNyawpre_AnMIRROR.fig')
saveas(gca,'turns_vs_headingNyawpre_AnMIRROR.png')
plot2svg('turns_vs_headingNyawpre_AnMIRROR.svg')

%% turns vs headingNyaw pre AnMIRROR CIRCMEAN leftNright
figure

plot_turnVSheadingNyawPRE_MIRRORleftNright

saveas(gca,'turns_vs_headingNyawpre_AnMIRROR_leftNright.fig')
saveas(gca,'turns_vs_headingNyawpre_AnMIRROR_leftNright.png')
plot2svg('turns_vs_headingNyawpre_AnMIRROR_leftNright.svg')

%% Adir vs headingNyaw pre CIRC MEAN
figure

plot_AdirPost_vs_headingNyawPre

saveas(gca,'Adir_vs_headingNyawpre.fig')
saveas(gca,'Adir_vs_headingNyawpre.png')
plot2svg('Adir_vs_headingNyawpre.svg')

%% Adir vs headingNyaw pre AnMIRROR CIRC MEAN
figure

plot_AdirPost_vs_headingNyawPre_MIRROR

saveas(gca,'Adir_vs_headingNyawpre_AnMIRROR.fig')
saveas(gca,'Adir_vs_headingNyawpre_AnMIRROR.png')
plot2svg('Adir_vs_headingNyawpre_AnMIRROR.svg')

%% Adir vs headingNyaw pre AnMIRROR CIRCMEAN leftNright
figure

plot_AdirPost_vs_headingNyawPre_MIRRORleftNright

saveas(gca,'Adir_vs_headingNyawpre_AnMIRROR_leftNright.fig')
saveas(gca,'Adir_vs_headingNyawpre_AnMIRROR_leftNright.png')
plot2svg('Adir_vs_headingNyawpre_AnMIRROR_leftNright.svg')

%% AdirMEAN vs headingNyaw pre CIRC MEAN
figure

plot_AdirMean_vs_headingNyawPre

saveas(gca,'AdirMean_vs_headingNyawpre.fig')
saveas(gca,'AdirMean_vs_headingNyawpre.png')
plot2svg('AdirMean_vs_headingNyawpre.svg')

%% AdirMEAN vs headingNyaw pre AnMIRROR CIRC MEAN
figure

plot_AdirMean_vs_headingNyawPre_MIRROR

saveas(gca,'AdirMean_vs_headingNyawpre_AnMIRROR.fig')
saveas(gca,'AdirMean_vs_headingNyawpre_AnMIRROR.png')
plot2svg('AdirMean_vs_headingNyawpre_AnMIRROR.svg')

%% AdirMEAN vs headingNyaw pre AnMIRROR CIRCMEAN leftNright
figure

plot_AdirMean_vs_headingNyawPre_MIRRORleftNright

saveas(gca,'AdirMean_vs_headingNyawpre_AnMIRROR_leftNright.fig')
saveas(gca,'AdirMean_vs_headingNyawpre_AnMIRROR_leftNright.png')
plot2svg('AdirMean_vs_headingNyawpre_AnMIRROR_leftNright.svg')

%% escape heading & AdirMEAN vs headingPre AnMIRROR CIRC MEAN ONE PANEL
figure

plot_headingPostNAdirMean_vs_headingPre_MIRROR_1panel

saveas(gca,'headingPostNAdirMean_vs_headingPre_AnMIRROR_1panel.fig')
saveas(gca,'headingPostNAdirMean_vs_headingPre_AnMIRROR_1panel.png')
plot2svg('headingPostNAdirMean_vs_headingPre_AnMIRROR_1panel.svg')

%% escape heading & AdirMEAN vs headingPre AnMIRROR CIRC MEAN
figure

plot_headingPostNAdirMean_vs_headingPre_MIRROR

saveas(gca,'headingPostNAdirMean_vs_headingPre_AnMIRROR.fig')
saveas(gca,'headingPostNAdirMean_vs_headingPre_AnMIRROR.png')
plot2svg('headingPostNAdirMean_vs_headingPre_AnMIRROR.svg')

%% escape heading & AdirMEAN vs headingPre AnMIRROR CIRC MEAN leftNright
figure

plot_headingPostNAdirMean_vs_headingPre_MIRRORleftNright

saveas(gca,'headingPostNAdirMean_vs_headingPre_AnMIRROR_leftNright.fig')
saveas(gca,'headingPostNAdirMean_vs_headingPre_AnMIRROR_leftNright.png')
plot2svg('headingPostNAdirMean_vs_headingPre_AnMIRROR_leftNright.svg')

%% rel heading-yaw score VS headingNyaw pre

figure
plot(360,1,'ok','MarkerFaceColor','r','markersize',5)
hold on
plot(360,1,'ok','MarkerFaceColor','b','markersize',5)
legend('heading','yaw','location','se')
xlabel('initial angle','fontsize',10) 

% subplot(1,2,1)
% yaw
angle_pre = heading_pre;
headingVSyaw = (abs(yaw_post) - abs(heading_post)) ./ (abs(yaw_post) + abs(heading_post));
plotcolor = 'b';
% plot_headingVSyaw_nowrap_csaps
plot_headingVSyaw_wrap_csaps

% xlabel('initial heading','fontsize',10) 
% ylabel('heading vs yaw score','fontsize',10) 
ylabel('|yaw|-|heading| / |yaw|+|heading|','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-1 1])
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-1 -.5 0 .5 1],'fontsize',8)
grid on

% subplot(1,2,2)
% yaw
angle_pre = yaw_pre;
headingVSyaw = (abs(yaw_post) - abs(heading_post)) ./ (abs(yaw_post) + abs(heading_post));
plotcolor = 'r';
% plot_headingVSyaw_nowrap_csaps
plot_headingVSyaw_wrap_csaps

% xlabel('initial yaw','fontsize',10) 
% ylabel('heading vs yaw score','fontsize',10) 
ylabel('|yaw|-|heading| / |yaw|+|heading|','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-1 1])
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-1 -.5 0 .5 1],'fontsize',8)
grid on


saveas(gca,'heading_yaw_score.fig')
saveas(gca,'heading_yaw_score.png')
plot2svg('heading_yaw_score.svg')

%% rel heading-yaw score VS ABS(headingNyaw pre)

figure
plot(360,1,'ok','MarkerFaceColor','r','markersize',5)
hold on
plot(360,1,'ok','MarkerFaceColor','b','markersize',5)
legend('heading','yaw','location','se')
xlabel('initial angle','fontsize',10) 

% subplot(1,2,1)
% yaw
angle_pre = abs(heading_pre);
plotcolor = 'b';
% plot_headingVSyaw_nowrap_csaps
plot_headingVSyaw_wrap_csaps

% xlabel('initial heading','fontsize',10) 
% ylabel('heading vs yaw score','fontsize',10) 
ylabel('|yaw|-|heading| / |yaw|+|heading|','fontsize',10) 
set(gca,'xlim',[0 180],'ylim',[-1 1])
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-1 -.5 0 .5 1],'fontsize',8)
grid on

% subplot(1,2,2)
% yaw
angle_pre = abs(yaw_pre);
plotcolor = 'r';
% plot_headingVSyaw_nowrap_csaps
plot_headingVSyaw_wrap_csaps

% xlabel('initial yaw','fontsize',10) 
% ylabel('heading vs yaw score','fontsize',10) 
ylabel('|yaw|-|heading| / |yaw|+|heading|','fontsize',10) 
set(gca,'xlim',[0 180],'ylim',[-1 1])
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-1 -.5 0 .5 1],'fontsize',8)
grid on

saveas(gca,'heading_yaw_score_abspre.fig')
saveas(gca,'heading_yaw_score_abspre.png')
plot2svg('heading_yaw_score_abspre.svg')

%% slip angle effect
% figure
% subplot(3,2,1)
% plot(slip_pre,heading_post,'.r')
% xlabel('initial slip')
% ylabel('escape heading')
% set(gca,'xlim',[-180 180],'ylim',[-180 180])
% set(gca,'XTick',[-180 0 180])
% set(gca,'YTick',[-180 0 180])
% grid on
% 
% subplot(3,2,3)
% plot(slip_pre,yaw_post,'.b')
% xlabel('initial slip')
% ylabel('escape yaw')
% set(gca,'xlim',[-180 180],'ylim',[-180 180])
% set(gca,'XTick',[-180 0 180])
% set(gca,'YTick',[-180 0 180])
% grid on
% 
% subplot(3,2,2)
% plot(slip_post,heading_post,'.r')
% xlabel('escape slip')
% ylabel('escape heading')
% set(gca,'xlim',[-180 180],'ylim',[-180 180])
% set(gca,'XTick',[-180 0 180])
% set(gca,'YTick',[-180 0 180])
% grid on
% 
% subplot(3,2,4)
% plot(slip_post,yaw_post,'.b')
% xlabel('escape slip')
% ylabel('escape yaw')
% set(gca,'xlim',[-180 180],'ylim',[-180 180])
% set(gca,'XTick',[-180 0 180])
% set(gca,'YTick',[-180 0 180])
% grid on
% 
% subplot(3,2,5)
% plot(heading_pre,slip_post,'.r')
% xlabel('heading pre')
% ylabel('slip post')
% set(gca,'xlim',[-180 180],'ylim',[-180 180])
% set(gca,'XTick',[-180 0 180])
% set(gca,'YTick',[-180 0 180])
% grid on
% 
% subplot(3,2,6)
% plot(yaw_pre,slip_post,'.b')
% xlabel('yaw pre')
% ylabel('slip post')
% set(gca,'xlim',[-180 180],'ylim',[-180 180])
% set(gca,'XTick',[-180 0 180])
% set(gca,'YTick',[-180 0 180])
% grid on
% 
% saveas(gca,'sideslip.fig')
% saveas(gca,'sideslip.png')
% plot2svg('sideslip.svg')
% 
% figure
% plot(slip_pre,slip_post,'.r')
% xlabel('slip pre')
% ylabel('slip post')
% set(gca,'xlim',[-180 180],'ylim',[-180 180])
% set(gca,'XTick',[-180 0 180])
% set(gca,'YTick',[-180 0 180])
% grid on
% 
% saveas(gca,'slip_pre_post.fig')
% saveas(gca,'slip_pre_post.png')
% plot2svg('slip_pre_post.svg')

%% plot pre&post velocity vectors
v_pre = V_pre .* cosd(heading_pre);
u_pre = V_pre .* sind(heading_pre);

v_post = V_post .* cosd(heading_post);
u_post = V_post .* sind(heading_post);

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
Z = rose(deg2rad(heading_post),36);
for i=1:length(Z)
    set(Z(i),'color','r','linewidth',2)
end
% x = get(Z,'Xdata');
% y = get(Z,'Ydata');
% g=patch(x,y,'r');

hold on
Y = rose(deg2rad(heading_pre),36);
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


%% plot pre&post Adir vectors in body coord
Ax_pre = A_pre .* cosd(Adir_body_pre);
Ay_pre = A_pre .* sind(Adir_body_pre);

Ax_post = A_post .* cosd(Adir_body_post);
Ay_post = A_post .* sind(Adir_body_post);

Ax_mean = A_hor_mean .* cosd(Adir_body_mean);
Ay_mean = A_hor_mean .* sind(Adir_body_mean);

Ax_max = A_hor_max .* cosd(accel_angle_hor_body_Ahormax);
Ay_max = A_hor_max .* sind(accel_angle_hor_body_Ahormax);

figure
maxHistogramValue = (max(A_hor_max));
maxHistogramValue = 20;

subplot(2,3,1)
% polar(0, maxHistogramValue,'-k')
% hold on
compass_zeroup(-Ay_pre,Ax_pre,'k');
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

% figure
subplot(2,3,4)
% polar(0, maxHistogramValue,'-k')
% hold on
Y = rose(deg2rad(Adir_body_pre),36);
% Y = rose(deg2rad(Adir_body_pre),36);
for i=1:length(Y)
    set(Y(i),'color','k','linewidth',1)
end
x = get(Y,'Xdata');
y = get(Y,'Ydata');
g=patch(x,y,[.5 .5 .5]);

subplot(2,3,5)
% polar(0, maxHistogramValue,'-k')
% hold on
X = rose(deg2rad(accel_angle_hor_body_Ahormax),36);
for i=1:length(X)
    set(X(i),'color','k','linewidth',1)
end
x = get(X,'Xdata');
y = get(X,'Ydata');
g=patch(x,y,'r');

subplot(2,3,6)
% polar(0, maxHistogramValue,'-k')
% hold on
X = rose(deg2rad(Adir_body_mean),36);
for i=1:length(X)
    set(X(i),'color','k','linewidth',1)
end
x = get(X,'Xdata');
y = get(X,'Ydata');
g=patch(x,y,'b');

saveas(gca,'Adir_body_vector_rose_pre_max_mean.fig')
saveas(gca,'Adir_body_vector_rose_pre_max_mean.png')
plot2svg('Adir_body_vector_rose_pre_max_mean.svg')


%% plot pre&post Adir vectors relative to stimulus
Ax_pre = A_pre .* cosd(Adir_pre);
Ay_pre = A_pre .* sind(Adir_pre);

Ax_post = A_post .* cosd(Adir_post);
Ay_post = A_post .* sind(Adir_post);

Ax_mean = A_hor_mean .* cosd(Adir_mean);
Ay_mean = A_hor_mean .* sind(Adir_mean);

Ax_max = A_hor_max .* cosd(stim_angle_accel_Ahormax);
Ay_max = A_hor_max .* sind(stim_angle_accel_Ahormax);

figure
maxHistogramValue = (max(A_hor_max));
maxHistogramValue = 20;

subplot(2,3,1)
% polar(0, maxHistogramValue,'-k')
% hold on
compass_zeroup(-Ay_pre,Ax_pre,'k');
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

maxHistogramValue = 25;

% figure
subplot(2,3,4)
polar(0, maxHistogramValue,'-k')
hold on
Y = rose(deg2rad(Adir_pre),36);
% Y = rose(deg2rad(Adir_body_pre),36);
for i=1:length(Y)
    set(Y(i),'color','k','linewidth',1)
end
x = get(Y,'Xdata');
y = get(Y,'Ydata');
g=patch(x,y,[.5 .5 .5]);

subplot(2,3,5)
polar(0, maxHistogramValue,'-k')
hold on
X = rose(deg2rad(stim_angle_accel_Ahormax),36);
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

saveas(gca,'Adir_stim_vector_rose_pre_max_mean.fig')
saveas(gca,'Adir_stim_vector_rose_pre_max_mean.png')
plot2svg('Adir_stim_vector_rose_pre_max_mean.svg')

%% plot Ahormax vectors in body & stim coord system
figure
maxHistogramValue = (max(A_hor_max));
maxHistogramValue = 20;

Ax_max = A_hor_max .* cosd(accel_angle_hor_body_Ahormax);
Ay_max = A_hor_max .* sind(accel_angle_hor_body_Ahormax);

subplot(2,2,1)
polar(0, maxHistogramValue,'-k')
hold on
compass_zeroup(-Ay_max,Ax_max,'k');
title('angle relative to body','fontsize',10) 

% figure
subplot(2,2,3)
% polar(0, maxHistogramValue,'-k')
% hold on
X = rose(deg2rad(accel_angle_hor_body_Ahormax),36);
for i=1:length(X)
    set(X(i),'color','k','linewidth',1)
end
x = get(X,'Xdata');
y = get(X,'Ydata');
g=patch(x,y,[.5 .5 .5]);

Ax_max = A_hor_max .* cosd(stim_angle_accel_Ahormax);
Ay_max = A_hor_max .* sind(stim_angle_accel_Ahormax);

subplot(2,2,2)
polar(0, maxHistogramValue,'-k')
hold on
compass_zeroup(-Ay_max,Ax_max,'k');
title('angle relative to stimulus','fontsize',10) 

% figure
subplot(2,2,4)
% polar(0, maxHistogramValue,'-k')
% hold on
X = rose(deg2rad(stim_angle_accel_Ahormax),36);
for i=1:length(X)
    set(X(i),'color','k','linewidth',1)
end
x = get(X,'Xdata');
y = get(X,'Ydata');
g=patch(x,y,[.5 .5 .5]);

saveas(gca,'Ahormax_vector_rose_stimNbody.fig')
saveas(gca,'Ahormax_vector_rose_stimNbody.png')
plot2svg('Ahormax_vector_rose_stimNbody.svg')


%% plot data at Ahormax MIRROR

% MIRROR TURN
for i=1:size(A_hor_max,1)
    if An_hor_Ahormax(i) < 0
        stim_angle_vel_Ahormax(i) = -stim_angle_vel_Ahormax(i);
        stim_angle_accel_Ahormax(i) = -stim_angle_accel_Ahormax(i);
        stim_angle_yaw_Ahormax(i) = -stim_angle_yaw_Ahormax(i);
        
        accel_angle_hor_vel_Ahormax(i) = -accel_angle_hor_vel_Ahormax(i);
        accel_angle_hor_body_Ahormax(i) = -accel_angle_hor_body_Ahormax(i);
        
        slip_Ahormax(i) = -slip_Ahormax(i);
        roll_Ahormax(i) = -roll_Ahormax(i);
        
        An_hor_Ahormax(i) = -An_hor_Ahormax(i);
    end
end

figure
maxHistogramValue = (max(A_hor_max));
maxHistogramValue = 10;
n_bins_rose = 90;
plot_color = 'r';
dh = 1;

plot_Ahormax_data

saveas(gca,'Ahormax_variables_mirror.fig')
saveas(gca,'Ahormax_variables_mirror.png')
plot2svg('Ahormax_variables_mirror.svg')



%% plot data at Ahormax MIRROR left n right

% MIRROR TURN
for i=1:size(A_hor_max,1)
%     if An_hor_Ahormax(i) < 0
        stim_angle_vel_Ahormax(i) = -stim_angle_vel_Ahormax(i);
        stim_angle_accel_Ahormax(i) = -stim_angle_accel_Ahormax(i);
        stim_angle_yaw_Ahormax(i) = -stim_angle_yaw_Ahormax(i);
        
        accel_angle_hor_vel_Ahormax(i) = -accel_angle_hor_vel_Ahormax(i);
        accel_angle_hor_body_Ahormax(i) = -accel_angle_hor_body_Ahormax(i);
        
        slip_Ahormax(i) = -slip_Ahormax(i);
        roll_Ahormax(i) = -roll_Ahormax(i);
        
        An_hor_Ahormax(i) = -An_hor_Ahormax(i);
%     end
end

% figure
maxHistogramValue = (max(A_hor_max));
maxHistogramValue = 10;
n_bins_rose = 90;
plot_color = 'b';
dh = 1;

plot_Ahormax_data

saveas(gca,'Ahormax_variables_mirror_leftright.fig')
saveas(gca,'Ahormax_variables_mirror_leftright.png')
plot2svg('Ahormax_variables_mirror_leftright.svg')




%% plot data at Amax MIRROR

% MIRROR TURN
for i=1:size(A_max,1)
    if An_hor_Amax(i) < 0
        stim_angle_vel_Amax(i) = -stim_angle_vel_Amax(i);
        stim_angle_accel_Amax(i) = -stim_angle_accel_Amax(i);
        stim_angle_yaw_Amax(i) = -stim_angle_yaw_Amax(i);
        
        accel_angle_hor_vel_Amax(i) = -accel_angle_hor_vel_Amax(i);
        accel_angle_hor_body_Amax(i) = -accel_angle_hor_body_Amax(i);
        
        slip_Amax(i) = -slip_Amax(i);
        roll_Amax(i) = -roll_Amax(i);
        
        An_hor_Amax(i) = -An_hor_Amax(i);
    end
end

figure
maxHistogramValue = (max(A_max));
maxHistogramValue = 10;
n_bins_rose = 90;
plot_color = 'r';
dh = 1;

plot_Amax_data

saveas(gca,'Amax_variables_mirror.fig')
saveas(gca,'Amax_variables_mirror.png')
plot2svg('Amax_variables_mirror.svg')



%% plot data at Amax MIRROR left n right

% MIRROR TURN
for i=1:size(A_max,1)
%     if An_hor_Amax(i) < 0
        stim_angle_vel_Amax(i) = -stim_angle_vel_Amax(i);
        stim_angle_accel_Amax(i) = -stim_angle_accel_Amax(i);
        stim_angle_yaw_Amax(i) = -stim_angle_yaw_Amax(i);
        
        accel_angle_hor_vel_Amax(i) = -accel_angle_hor_vel_Amax(i);
        accel_angle_hor_body_Amax(i) = -accel_angle_hor_body_Amax(i);
        
        slip_Amax(i) = -slip_Amax(i);
        roll_Amax(i) = -roll_Amax(i);
        
        An_hor_Amax(i) = -An_hor_Amax(i);
%     end
end

% figure
maxHistogramValue = (max(A_max));
maxHistogramValue = 10;
n_bins_rose = 90;
plot_color = 'b';
dh = 1;

plot_Amax_data

saveas(gca,'Amax_variables_mirror_leftright.fig')
saveas(gca,'Amax_variables_mirror_leftright.png')
plot2svg('Amax_variables_mirror_leftright.svg')




%% rose plots heading yaw POST

figure
% subplot(1,2,1)
Y = rose(deg2rad(heading_post),36);
for i=1:length(Y)
    set(Y(i),'color','k','linewidth',1)
end 
x = get(Y,'Xdata');
y = get(Y,'Ydata');
g=patch(x,y,'r');
% alpha(.5)
hold on

Z = rose(deg2rad(yaw_post),36);
for i=1:length(Z)
    set(Z(i),'color','k','linewidth',1)
end
x = get(Z,'Xdata');
y = get(Z,'Ydata');
g=patch(x,y,'b');
alpha(.5)

legend('','heading','','yaw')


saveas(gca,'rose_heading_yaw_post.fig')
saveas(gca,'rose_heading_yaw_post.png')
plot2svg('rose_heading_yaw_post.svg')


%% rose plots heading yaw PRE&POST


figure
subplot(1,2,1)
Z = rose(deg2rad(heading_post),36);
for i=1:length(Z)
    set(Z(i),'color','r','linewidth',2)
end
% x = get(Z,'Xdata');
% y = get(Z,'Ydata');
% g=patch(x,y,'r');

hold on
Y = rose(deg2rad(heading_pre),36);
for i=1:length(Y)
    set(Y(i),'color','k','linewidth',2)
end 
% x = get(Y,'Xdata');
% y = get(Y,'Ydata');
% g=patch(x,y,'b');
title('heading','fontsize',10) 

% figure
subplot(1,2,2)
Z = rose(deg2rad(yaw_post),36);
for i=1:length(Z)
    set(Z(i),'color','b','linewidth',2)
end
% x = get(Z,'Xdata');
% y = get(Z,'Ydata');
% g=patch(x,y,'r');

hold on
Y = rose(deg2rad(yaw_pre),36);
for i=1:length(Y)
    set(Y(i),'color','k','linewidth',2)
end 
% x = get(Y,'Xdata');
% y = get(Y,'Ydata');
% g=patch(x,y,'b');
title('yaw','fontsize',10) 

saveas(gca,'rose_heading_yaw_pre_post.fig')
saveas(gca,'rose_heading_yaw_pre_post.png')
plot2svg('rose_heading_yaw_pre_post.svg')

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

figure
% subplot(3,2,1)

dh = 360/(36)
bins = [-180+dh/2:dh:180-dh/2];

h1 = hist(heading_post,bins);
h2 = hist(yaw_post,bins);

bar(bins,h1,'FaceColor','r','EdgeColor','k')
hold on
bar(bins,h2,'FaceColor','b','EdgeColor','k')
% bar(bins,[h1;h2]')

legend('heading','yaw')
alpha(.5)
set(gca,'xlim',[-180 180])
set(gca,'XTick',[-180:90:180])
set(gca,'YTick',[0:25:25],'fontsize',8) 
title('escape direction','fontsize',10) 

saveas(gca,'hist_heading_yaw_post.fig')
saveas(gca,'hist_heading_yaw_post.png')
plot2svg('hist_heading_yaw_post.svg')

set(gca,'xlim',[-90 90])
set(gca,'XTick',[-180:90:180])

saveas(gca,'hist_heading_yaw_post_min90to90.fig')
saveas(gca,'hist_heading_yaw_post_min90to90.png')
plot2svg('hist_heading_yaw_post_min90to90.svg')

figure
subplot(2,1,1)
bar(bins,h1,'FaceColor',[.5 .5 .5],'EdgeColor','k')
hold on
set(gca,'xlim',[-180 180])
set(gca,'XTick',[-180:90:180])
set(gca,'YTick',[0:25:25],'fontsize',8) 
title('escape heading angle','fontsize',10) 

subplot(2,1,2)
bar(bins,h2,'FaceColor',[.5 .5 .5],'EdgeColor','k')
% bar(bins,[h1;h2]')

% legend('heading','yaw')
% alpha(.5)
set(gca,'xlim',[-180 180])
set(gca,'XTick',[-180:90:180])
set(gca,'YTick',[0:25:25],'fontsize',8) 
title('escape body angle','fontsize',10) 

saveas(gca,'hist_heading_yaw_post_sep.fig')
saveas(gca,'hist_heading_yaw_post_sep.png')
plot2svg('hist_heading_yaw_post_sep.svg')

% -90 to 90 deg
subplot(2,1,1)
set(gca,'xlim',[-90 90])
set(gca,'XTick',[-180:90:180])

subplot(2,1,2)
set(gca,'xlim',[-90 90])
set(gca,'XTick',[-180:90:180])

saveas(gca,'hist_heading_yaw_post_sep_min90to90.fig')
saveas(gca,'hist_heading_yaw_post_sep_min90to90.png')
plot2svg('hist_heading_yaw_post_sep_min90to90.svg')

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

%% heading vs slip post (heading pre colorcode)
% figure
% hold on
% for i = 1:length(heading_post)
%     plot(heading_post(i),slip_post(i),'ok','MarkerFaceColor',cmap_360(round(heading_pre(i))+180,:),'linewidth',1)
% end
% axis equal
% grid on
% xlabel('escape heading','fontsize',10) 
% ylabel('escape slip','fontsize',10) 
% set(gca,'xlim',[-180 180],'ylim',[-90 90])
% set(gca,'XTick',[-180:90:180],'YTick',[-180:90:180],'fontsize',8) 
% 
% colormap(cmap_360(1:2:end,:))
% colorbar('YTick',[1:45:181],'YTickLabel',[-180:90:180],'fontsize',8) 
% %title('escape body angle','fontsize',10) 
% 
% 
% saveas(gca,'heading_slip_post_headingpre_colorcode.fig')
% saveas(gca,'heading_slip_post_headingpre_colorcode.png')
% plot2svg('heading_slip_post_headingpre_colorcode.svg')

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