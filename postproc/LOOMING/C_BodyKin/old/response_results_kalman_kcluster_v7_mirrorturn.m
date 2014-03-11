clc
clear
close all

addpath('/home/florian/Dropbox/WORK/toolbox')
addpath('/home/florian/Dropbox/WORK/toolbox/plot2svg')

% load('flightpathDB_pos_INCq.mat')
load('flightpathDB_pos_qbodyEKF_NOroll_9clusters_2.5n-3.3n3_response.mat')
mkdir('response_figs')
cd('response_figs')

%% settings
linewidth_timelines = 1;
skip = 50;

% heatmap resolution
nx = 1000;
ny = 100;

cmap_180 = jet(180);

%% data import
IDX = pathDB.IDX;
cmap_k = settings.cmap_k;
cmap_360 = settings.cmap_360;

t = pathDB.t;
V = pathDB.V;
An_hor = pathDB.An_hor;
An_hor_abs = abs(pathDB.An_hor);
At_hor = pathDB.At_hor;

stim_angle_vel = pathDB.stim_angle_vel;
stim_angle_yaw = pathDB.stim_angle_yaw;
slip = pathDB.slip;

teta = patternDB.teta;

stim_angle_vel_pre_resp = responseDB.stim_angle_vel_pre_resp;
stim_angle_vel_pre_turn = responseDB.stim_angle_vel_pre_turn;

stim_angle_vel_post_resp = responseDB.stim_angle_vel_post_resp;
stim_angle_vel_post_turn = responseDB.stim_angle_vel_post_turn;

stim_angle_yaw_pre_resp = responseDB.stim_angle_yaw_pre_resp;
stim_angle_yaw_pre_turn = responseDB.stim_angle_yaw_pre_turn;

stim_angle_yaw_post_resp = responseDB.stim_angle_yaw_post_resp;
stim_angle_yaw_post_turn = responseDB.stim_angle_yaw_post_turn;

slip_pre_resp = responseDB.slip_pre_resp;
slip_pre_turn = responseDB.slip_pre_turn;

slip_post_resp = responseDB.slip_post_resp;
slip_post_turn = responseDB.slip_post_turn;

V_trig2resp = responseDB.V_trig2resp;
V_pre_resp = responseDB.V_pre_resp;
V_post_resp = responseDB.V_post_resp;

V_pre_accel = responseDB.V_pre_accel;
V_post_accel = responseDB.V_post_accel;

An_hor_max = responseDB.An_hor_max;
At_hor_max = responseDB.At_hor_max;
At_hor_min = responseDB.At_hor_min;

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

%% timelines NO skip 

IDX_plot = IDX;
t_plot = t;

V_plot = V;
dV_plot = dV;
An_hor_plot = An_hor;
At_hor_plot = At_hor;

stim_angle_vel_plot = stim_angle_vel;
stim_angle_yaw_plot = stim_angle_yaw;
slip_plot = slip;

stim_angle_vel_temp = stim_angle_vel;
stim_angle_yaw_temp = stim_angle_yaw;
slip_temp = slip;

% remove jumps from plots (+/-180deg)
for i=1:size(stim_angle_vel_plot,2)
    for j=2:size(stim_angle_vel_plot,1)
        if abs(stim_angle_vel_temp(j,i) - stim_angle_vel_temp(j-1,i)) > 90
            stim_angle_vel_plot(j-1,i) = nan;
        end
        if abs(stim_angle_yaw_temp(j,i) - stim_angle_yaw_temp(j-1,i)) > 90
            stim_angle_yaw_plot(j-1,i) = nan;
        end
        if abs(slip_temp(j,i) - slip_temp(j-1,i)) > 90
            slip_plot(j-1,i) = nan;
        end
    end
end

%% timelines NO skip NO MIRROR TURN

t_shift = responseDB.t_resp;
cmap_plot = cmap_360;
color_var = round(heading_pre)+180;

% plot_flightpath_timeline_tshift_cluster_slip
% plot_flightpath_timeline_tshift_cluster_slip_long
% plot_flightpath_timeline_tshift_headingstart_posangles
plot_flightpath_timeline_tshift_headingstart_subset

subplot(3,2,4)
axis([-.025 .05 -20 20])
set(gca,'YTick',-20:20:20,'fontsize',12)

saveas(gca,'flightpaths_headingUA_headingstart.fig')
saveas(gca,'flightpaths_headingUA_headingstart.png')
plot2svg('flightpaths_headingUA_headingstart.svg')

%% heatmap hist norm NO MIRROR

t_shift = responseDB.t_resp;
cmap_plot = cmap_360;
color_var = round(heading_pre)+180;

t_hist = [];
for i=1:size(V,2)    
    t_hist = [t_hist;t-t_shift(i)];
end
plot_flightpath_timeline_histograms_norm_subset

subplot(3,2,4)
axis([-.025 .05 -20 20])
set(gca,'YTick',[-20;0;20],'YTicklabel',[20;0;-20],'fontsize',12) 

saveas(gca,'flightpath_hist_norm_slip_tresp.fig')
saveas(gca,'flightpath_hist_norm_slip_tresp.png')
plot2svg('flightpath_hist_norm_slip_tresp.svg')


%% timelines NO skip INC MIRROR TURN
t_shift = responseDB.t_resp;
cmap_plot = cmap_180;
color_var = abs(round(heading_pre));

% MIRROR TURN
% if An_hor_max < 0: mirror all
for i=1:size(stim_angle_vel_plot,2)
    if An_hor_max(i) < 0
        stim_angle_vel_plot(:,i) = -stim_angle_vel_plot(:,i);
        stim_angle_yaw_plot(:,i) = -stim_angle_yaw_plot(:,i);
        slip_plot(:,i) = -slip_plot(:,i);
        An_hor_plot(:,i) = -An_hor_plot(:,i);
    end
end

% % MIRROR TURN INC REV
% % if An_hor_max < 0: mirror all
% for i=1:size(stim_angle_vel_plot,2)
%     if An_hor_max(i) < 0 && settings.expansion.HorPos(i) == 0
%         stim_angle_vel_plot(:,i) = -stim_angle_vel_plot(:,i);
%         stim_angle_yaw_plot(:,i) = -stim_angle_yaw_plot(:,i);
%         slip_plot(:,i) = -slip_plot(:,i);
%         An_hor_plot(:,i) = -An_hor_plot(:,i);
%     elseif An_hor_max(i) > 0 && settings.expansion.HorPos(i) == 1
%         stim_angle_vel_plot(:,i) = -stim_angle_vel_plot(:,i);
%         stim_angle_yaw_plot(:,i) = -stim_angle_yaw_plot(:,i);
%         slip_plot(:,i) = -slip_plot(:,i);
%         An_hor_plot(:,i) = -An_hor_plot(:,i);
%     end
% end

% plot_flightpath_timeline_tshift_cluster_slip
% plot_flightpath_timeline_tshift_cluster_slip_long
% plot_flightpath_timeline_tshift_headingstart_posangles
plot_flightpath_timeline_tshift_headingstart_subset

saveas(gca,'flightpaths_headingUA_headingstart_mirror.fig')
saveas(gca,'flightpaths_headingUA_headingstart_mirror.png')
plot2svg('flightpaths_headingUA_headingstart_mirror.svg')

%% heatmap hist norm INC MIRROR
% t_resp timeshift
t_shift = responseDB.t_resp;
cmap_plot = cmap_180;
color_var = abs(round(heading_pre));

t_hist = [];
for i=1:size(V,2)    
    t_hist = [t_hist;t-t_shift(i)];
end
plot_flightpath_timeline_histograms_norm_subset

saveas(gca,'flightpath_hist_norm_slip_tresp.fig')
saveas(gca,'flightpath_hist_norm_slip_tresp.png')
plot2svg('flightpath_hist_norm_slip_tresp.svg')


%% timelines INC skip

IDX_plot = IDX(1:skip:end,:);
t_plot = t(1:skip:end,:);

V_plot = V(1:skip:end,:);
dV_plot = dV(1:skip:end,:);
An_hor_plot = An_hor(1:skip:end,:);
At_hor_plot = At_hor(1:skip:end,:);

stim_angle_vel_plot = stim_angle_vel(1:skip:end,:);
stim_angle_yaw_plot = stim_angle_yaw(1:skip:end,:);
slip_plot = slip(1:skip:end,:);

stim_angle_vel_temp = stim_angle_vel(1:skip:end,:);
stim_angle_yaw_temp = stim_angle_yaw(1:skip:end,:);
slip_temp = slip(1:skip:end,:);

% remove jumps from plots (+/-180deg)
for i=1:size(stim_angle_vel_plot,2)
    for j=2:size(stim_angle_vel_plot,1)
        if abs(stim_angle_vel_temp(j,i) - stim_angle_vel_temp(j-1,i)) > 90
            stim_angle_vel_plot(j-1,i) = nan;
        end
        if abs(stim_angle_yaw_temp(j,i) - stim_angle_yaw_temp(j-1,i)) > 90
            stim_angle_yaw_plot(j-1,i) = nan;
        end
        if abs(slip_temp(j,i) - slip_temp(j-1,i)) > 90
            slip_plot(j-1,i) = nan;
        end
    end
end

%% clusters NO MIRROR
t_shift = responseDB.t_resp;

% plot_flightpath_timeline_tshift_cluster_slip
% plot_flightpath_timeline_tshift_cluster_slip_axisfix
% plot_flightpath_timeline_tshift_cluster_slip_long
plot_flightpath_timeline_tshift_cluster_slip_subset

subplot(3,2,4)
axis([-.025 .05 -20 20])

saveas(gca,'flightpaths_headingUA_tresp.fig')
saveas(gca,'flightpaths_headingUA_tresp.png')
plot2svg('flightpaths_headingUA_tresp.svg')

%% INC slip INC MIRROR
t_shift = responseDB.t_resp;

% MIRROR TURN
% if An_hor_max < 0: mirror all
for i=1:size(stim_angle_vel_plot,2)
    if An_hor_max(i) < 0
        stim_angle_vel_plot(:,i) = -stim_angle_vel_plot(:,i);
        stim_angle_yaw_plot(:,i) = -stim_angle_yaw_plot(:,i);
        slip_plot(:,i) = -slip_plot(:,i);
        An_hor_plot(:,i) = -An_hor_plot(:,i);
    end
end

% plot_flightpath_timeline_tshift_cluster_slip
% plot_flightpath_timeline_tshift_cluster_slip_axisfix
% plot_flightpath_timeline_tshift_cluster_slip_long
plot_flightpath_timeline_tshift_cluster_slip_subset

saveas(gca,'flightpaths_headingUA_tresp_mirror.fig')
saveas(gca,'flightpaths_headingUA_tresp_mirror.png')
plot2svg('flightpaths_headingUA_tresp_mirror.svg')


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

%% At_max & At_min vs heading_start

figure
plot(heading_start,abs(An_hor_max),'ok','markerfacecolor','y','markersize',5)
hold on
plot(heading_pre(At_hor_max>0),At_hor_max(At_hor_max>0),'ok','markerfacecolor','r','markersize',5)
plot(heading_pre(At_hor_min<0),At_hor_min(At_hor_min<0),'ok','markerfacecolor','b','markersize',5)
legend('turn','accel','decel','location','sw')
grid on
axis([0 180 -20 20])
set(gca,'XTick',0:90:180) 
set(gca,'YTick',-20:20:20,'fontsize',12)
xlabel('initial heading','fontsize',18) 
ylabel('Amax','fontsize',18) 

saveas(gca,'Amax_VS_heading_pre.fig')
saveas(gca,'Amax_VS_heading_pre.png')
plot2svg('Amax_VS_heading_pre.svg')

%% heatmap hist all An & At
binx = 0:.5:20;
biny = -20:.5:20;

An_hist = An_hor(:);
At_hist = At_hor(:);

AnAt_hist = hist3([At_hist,An_hist], {biny binx});
AnAt_hist_log = log(AnAt_hist);

% imagesc(binx,-biny,AnAt_hist)
figure
imagesc(binx,-biny,AnAt_hist_log)
axis equal
axis([0 20 -20 20])
set(gca,'XTick',0:10:20) 
set(gca,'YTick',-20:10:20,'YTicklabel',20:-10:-20,'fontsize',12)
xlabel('An','fontsize',18) 
ylabel('At','fontsize',18) 
colorbar

saveas(gca,'AnAt_hist_all.fig')
saveas(gca,'AnAt_hist_all.png')
plot2svg('AnAt_hist_all.svg')

%% heatmap hist An & At t0<t<tstop
% binx = 0:20;
% biny = -20:20;
% 
% t_stop = max([responseDB.t_turn_stop responseDB.t_accel_stop]')';
% 
% An_hist = [];
% At_hist = [];
% 
% for i = 1:length(t_stop)
%     An_hist_temp = An_hist
%     An_hist(end:end+length(An_hor(:))) = An_hor(:);
% end
% 
% AnAt_hist = hist3([At_hist,An_hist], {biny binx});
% AnAt_hist_log = log(AnAt_hist);
% 
% % imagesc(binx,-biny,AnAt_hist)
% figure
% imagesc(binx,-biny,AnAt_hist_log)
% axis equal
% axis([0 20 -20 20])
% set(gca,'XTick',0:10:20) 
% set(gca,'YTick',-20:10:20,'YTicklabel',20:-10:-20,'fontsize',12)
% xlabel('An','fontsize',18) 
% ylabel('At','fontsize',18) 
% 
% saveas(gca,'AnAt_hist_t0totstop.fig')
% saveas(gca,'AnAt_hist_t0totstop.png')
% plot2svg('AnAt_hist_t0totstop.svg')
% 

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
xlabel('turn response time','fontsize',18) 
ylabel('acceleration response time','fontsize',18) 
set(gca,'xlim',[0 .15],'ylim',[0 .15])
set(gca,'XTick',[0:.05:.2])
set(gca,'YTick',[0:.05:.2],'fontsize',12)

saveas(gca,'response_time_tAn_vs_tAt.fig')
saveas(gca,'response_time_tAn_vs_tAt.png')
plot2svg('response_time_tAn_vs_tAt.svg')


%% turn vs heading pre
figure
plotcolor = 'r'
angle_pre = responseDB.stim_angle_vel_pre_resp;
turn = turn_angle_vel
% plot_angle_pre_turn_vectors_nowrap
% plot_angle_pre_turn_vectors
% plot_angle_pre_turn_vectors_csaps
plot_angle_pre_turn_vectors_csaps_wrap

xlabel('heading pre')
ylabel('turn angle')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180])

saveas(gca,'turnangle_vs_preheading.fig')
saveas(gca,'turnangle_vs_preheading.png')

%% yaw turn vs yaw pre
figure
plotcolor = 'b'
angle_pre = responseDB.stim_angle_yaw_pre_resp;
turn = turn_angle_yaw
% plot_angle_pre_turn_vectors_nowrap
% plot_angle_pre_turn_vectors
plot_angle_pre_turn_vectors_csaps
% plot_angle_pre_turn_vectors_csaps_wrap

xlabel('heading pre')
ylabel('turn angle')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180])

saveas(gca,'yawpost_vs_yawpre.fig')
saveas(gca,'yawpost_vs_yawpre.png')

%% yaw turn vs yaw pre
figure
plotcolor = 'b'
angle_pre = responseDB.stim_angle_yaw_pre_resp;
turn = turn_angle_yaw
% plot_angle_pre_turn_vectors_nowrap
% plot_angle_pre_turn_vectors
plot_angle_pre_turn_vectors_csaps
% plot_angle_pre_turn_vectors_csaps_wrap

xlabel('heading pre')
ylabel('turn angle')
set(gca,'xlim',[-180 180],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180])

saveas(gca,'yawpost_vs_yawpre.fig')
saveas(gca,'yawpost_vs_yawpre.png')

%% heading turn n yaw turn vs heading pre
figure
plot(360,360,'ok','MarkerFaceColor','r','markersize',5)
hold on
plot(360,360,'ok','MarkerFaceColor','b','markersize',5)
legend('heading','yaw')

order = 3;
dn=20
dm=20

% yaw
angle_pre = responseDB.stim_angle_yaw_pre_resp;
turn = turn_angle_yaw;
plotcolor = 'b';
plot_angle_pre_turn_vectors_nowrap
% plot_angle_pre_turn_vectors

% heading
% angle_pre = responseDB.stim_angle_vel_pre_resp;
angle_pre = heading_pre;
turn = turn_angle_vel;
plotcolor = 'r';
% plot_angle_pre_turn_vectors_nowrap
plot_angle_pre_turn_vectors

xlabel('heading pre','fontsize',18) 
ylabel('turn angle','fontsize',18) 
set(gca,'xlim',[-180 180],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',12)

saveas(gca,'heading_vs_yaw.fig')
saveas(gca,'heading_vs_yaw.png')
plot2svg('heading_vs_yaw.svg')

%% rel heading-yaw score VS heading pre

figure

order = 3;
dn=20
dm=20

% yaw
angle_pre = heading_pre;
headingVSyaw = (abs(yaw_post) - abs(heading_post)) ./ (abs(yaw_post) + abs(heading_post));
plotcolor = 'b';
plot_headingVSyaw_nowrap
% plot_angle_pre_turn_vectors

xlabel('heading pre','fontsize',18) 
ylabel('heading vs yaw score','fontsize',18) 
ylabel('|yaw|-|heading| / |yaw|+|heading|','fontsize',18) 
set(gca,'xlim',[-180 180],'ylim',[-1 1])
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-1 -.5 0 .5 1],'fontsize',12)
grid on

saveas(gca,'heading_yaw_ratio.fig')
saveas(gca,'heading_yaw_ratio.png')
plot2svg('heading_yaw_ratio.svg')

%% rel heading-yaw score VS ABS(heading pre)

figure

order = 3;
dn=20
dm=20

% yaw
angle_pre = abs(heading_pre);
headingVSyaw = (abs(yaw_post) - abs(heading_post)) ./ (abs(yaw_post) + abs(heading_post));
plotcolor = 'b';
plot_headingVSyaw_nowrap
% plot_angle_pre_turn_vectors

xlabel('heading pre','fontsize',18) 
% ylabel('heading vs yaw score','fontsize',18) 
ylabel('|yaw|-|heading| / |yaw|+|heading|','fontsize',18) 
set(gca,'xlim',[0 180],'ylim',[-1 1])
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-1 -.5 0 .5 1],'fontsize',12)
grid on

saveas(gca,'heading_yaw_ratio_absheadingpre.fig')
saveas(gca,'heading_yaw_ratio_absheadingpre.png')
plot2svg('heading_yaw_ratio_absheadingpre.svg')
%% escape angle vs heading
% angle_pre = responseDB.stim_angle_vel_pre_resp;
% angle_post = responseDB.stim_angle_vel_post_resp;
% 
% figure
% plotcolor = 'k'
% plot_angle_pre_turn
% 
% xlabel('heading pre')
% ylabel('turn angle')
% set(gca,'xlim',[-180 180],'ylim',[-180 180])
% set(gca,'XTick',[-180 -90 0 90 180])
% set(gca,'YTick',[-180 -90 0 90 180])
% saveas(gca,'escapeangle_vs_preheading.fig')
% saveas(gca,'escapeangle_vs_preheading.png')

%% slip angle effect
% figure
% plot(responseDB.slip_pre_resp,responseDB.stim_angle_vel_post_turn,'.r')
% hold on
% plot(responseDB.slip_pre_resp,responseDB.stim_angle_yaw_post_turn,'.b')
% xlabel('slip pre')
% ylabel('heading and yaw post')
% set(gca,'xlim',[-180 180],'ylim',[-180 180])
% set(gca,'XTick',[-180 0 180])
% set(gca,'YTick',[-180 0 180])
% grid on
% 
% figure
% plot(responseDB.slip_pre_resp,responseDB.slip_post_resp,'.r')
% xlabel('slip pre')
% ylabel('slip post')
% set(gca,'xlim',[-180 180],'ylim',[-180 180])
% set(gca,'XTick',[-180 0 180])
% set(gca,'YTick',[-180 0 180])
% grid on

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

% title('heading histogram','fontsize',18) 
saveas(gca,'heading_vector_rose_pre_post.fig')
saveas(gca,'heading_rose_pre_post.png')
plot2svg('heading_vector_rose_pre_post.svg')


%% rose plots heading yaw POST

figure
subplot(1,2,1)
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
Z = rose(deg2rad(yaw_post),36);
for i=1:length(Z)
    set(Z(i),'color','r','linewidth',2)
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
title('yaw','fontsize',18) 

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
title('heading','fontsize',18) 

saveas(gca,'rose_heading_yaw_pre_post.fig')
saveas(gca,'rose_heading_yaw_pre_post.png')
plot2svg('rose_heading_yaw_pre_post.svg')

%% plot rose teta_max ON/OFF

maxHistogramValue = 5;
figure

plot_rose_heading_tetamax_ONOFF

saveas(gca,'heading_rose_post_turn.fig')
saveas(gca,'heading_rose_pre_post.png')
plot2svg('heading_rose_pre_post.svg')

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
set(gca,'YTick',[0:25:25],'fontsize',12) 
title('escape direction','fontsize',18) 

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
set(gca,'YTick',[0:25:25],'fontsize',12) 
title('escape heading angle','fontsize',18) 

subplot(2,1,2)
bar(bins,h2,'FaceColor',[.5 .5 .5],'EdgeColor','k')
% bar(bins,[h1;h2]')

% legend('heading','yaw')
% alpha(.5)
set(gca,'xlim',[-180 180])
set(gca,'XTick',[-180:90:180])
set(gca,'YTick',[0:25:25],'fontsize',12) 
title('escape body angle','fontsize',18) 

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
figure
hold on
for i = 1:length(heading_post)
    plot(heading_post(i),yaw_post(i),'ok','MarkerFaceColor',cmap_360(round(heading_pre(i))+180,:),'linewidth',1)
end
axis equal
grid on
xlabel('escape heading','fontsize',18) 
ylabel('escape body angle','fontsize',18) 
set(gca,'xlim',[-180 180],'ylim',[-180 180])
set(gca,'XTick',[-180:90:180],'YTick',[-180:90:180],'fontsize',12) 

colormap(cmap_360(1:2:end,:))
colorbar('YTick',[1:45:181],'YTickLabel',[-180:90:180],'fontsize',12) 
%title('escape body angle','fontsize',18) 


saveas(gca,'heading_yaw_post_headingpre_colorcode.fig')
saveas(gca,'heading_yaw_post_headingpre_colorcode.png')
plot2svg('heading_yaw_post_headingpre_colorcode.svg')

%% heading vs slip post (heading pre colorcode)
figure
hold on
for i = 1:length(heading_post)
    plot(heading_post(i),slip_post(i),'ok','MarkerFaceColor',cmap_360(round(heading_pre(i))+180,:),'linewidth',1)
end
axis equal
grid on
xlabel('escape heading','fontsize',18) 
ylabel('escape slip','fontsize',18) 
set(gca,'xlim',[-180 180],'ylim',[-90 90])
set(gca,'XTick',[-180:90:180],'YTick',[-180:90:180],'fontsize',12) 

colormap(cmap_360(1:2:end,:))
colorbar('YTick',[1:45:181],'YTickLabel',[-180:90:180],'fontsize',12) 
%title('escape body angle','fontsize',18) 


saveas(gca,'heading_slip_post_headingpre_colorcode.fig')
saveas(gca,'heading_slip_post_headingpre_colorcode.png')
plot2svg('heading_slip_post_headingpre_colorcode.svg')

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

plot_escape_perf_ONOFF_NO165tetamax

% title('acceleration based response time')
saveas(gca,'escape performance_ONOFF.fig')
saveas(gca,'escape performance_ONOFF.png')
plot2svg('escape performance_ONOFF.svg')

%% teta_max escape performance: heading, yaw, resp time, U, An, At

plot_escape_perf_tetamax

% title('acceleration based response time')
saveas(gca,'escape performance_tetamax.fig')
saveas(gca,'escape performance_tetamax.png')
plot2svg('escape performance_tetamax.svg')

cd ..