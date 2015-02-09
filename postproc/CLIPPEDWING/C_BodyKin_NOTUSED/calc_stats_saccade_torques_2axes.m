xp = [min(turn_angle_vel_mirror),max(turn_angle_vel_mirror)];
xp = [0,180];

% max&min
stats_M_R_norm_max_startstop = regstats(M_R_norm_max_startstop,turn_angle_vel_mirror,'linear');
stats_M_R_norm_min_0stop = regstats(M_R_norm_min_0stop,turn_angle_vel_mirror,'linear');
stats_M_L_norm_max_startstop = regstats(M_L_norm_max_startstop,turn_angle_vel_mirror,'linear');
stats_M_L_norm_min_0stop = regstats(M_L_norm_min_0stop,turn_angle_vel_mirror,'linear');

stats_M_axis1_norm_max_startcut = regstats(M_axis1_norm_max_startcut,turn_angle_vel_mirror,'linear');
stats_M_axis1normal_norm_max_startcut = regstats(M_axis1normal_norm_max_startcut,turn_angle_vel_mirror,'linear');

stats_M_axis2_norm_min_cutstop = regstats(M_axis2_norm_min_cutstop,turn_angle_vel_mirror,'linear');
stats_M_axis2normal_norm_min_cutstop = regstats(M_axis2normal_norm_min_cutstop,turn_angle_vel_mirror,'linear');

stats_Myaw_norm_max_startstop = regstats(Myaw_norm_max_startstop,turn_angle_vel_mirror,'linear');

% means
stats_M_R_norm_mean_startstop = regstats(M_R_norm_mean_startstop,turn_angle_vel_mirror,'linear');
stats_M_L_norm_mean_startstop = regstats(M_L_norm_mean_startstop,turn_angle_vel_mirror,'linear');

stats_M_axis1_norm_mean_startcut = regstats(M_axis1_norm_mean_startcut,turn_angle_vel_mirror,'linear');
stats_M_axis1normal_norm_mean_startcut = regstats(M_axis1normal_norm_mean_startcut,turn_angle_vel_mirror,'linear');

stats_M_axis2_norm_mean_cutstop = regstats(M_axis2_norm_mean_cutstop,turn_angle_vel_mirror,'linear');
stats_M_axis2normal_norm_mean_cutstop = regstats(M_axis2normal_norm_mean_cutstop,turn_angle_vel_mirror,'linear');

stats_Myaw_norm_mean_startstop = regstats(Myaw_norm_mean_startstop,turn_angle_vel_mirror,'linear');

%% plot
% figure

%% Mean
% M_axis1_mean & M_axis2_mean
subplot(3,2,1)
hold on

yp = polyval([stats_M_axis1_norm_mean_startcut.beta(2) stats_M_axis1_norm_mean_startcut.beta(1)],xp);
plot(xp,yp)
plot(turn_angle_vel_mirror,M_axis1_norm_mean_startcut,'.')

yp = polyval([stats_M_axis2_norm_mean_cutstop.beta(2) stats_M_axis2_norm_mean_cutstop.beta(1)],xp);
plot(xp,yp,'r')
plot(turn_angle_vel_mirror,M_axis2_norm_mean_cutstop,'.r')

axis([0 180 -.045 .045])
ylabel('Maxis1&2')

% M_axis1normal_mean & M_axis2normal_mean
subplot(3,2,3)
hold on

yp = polyval([stats_M_axis1normal_norm_mean_startcut.beta(2) stats_M_axis1normal_norm_mean_startcut.beta(1)],xp);
plot(xp,yp)
plot(turn_angle_vel_mirror,M_axis1normal_norm_mean_startcut,'.')

yp = polyval([stats_M_axis2normal_norm_mean_cutstop.beta(2) stats_M_axis2normal_norm_mean_cutstop.beta(1)],xp);
plot(xp,yp,'r')
plot(turn_angle_vel_mirror,M_axis2normal_norm_mean_cutstop,'.r')

axis([0 180 -.045 .045])
ylabel('Maxis1&2normal')

% Myaw mean
subplot(3,2,5)
hold on

yp = polyval([stats_Myaw_norm_mean_startstop.beta(2) stats_Myaw_norm_mean_startstop.beta(1)],xp);
plot(xp,yp)
plot(turn_angle_vel_mirror,Myaw_norm_mean_startstop,'.')

axis([0 180 0 .09])
xlabel('turn angle')
ylabel('Myaw')

%% Max & Min
% M_axis1_max & M_axis2_min
subplot(3,2,2)
hold on

yp = polyval([stats_M_axis1_norm_max_startcut.beta(2) stats_M_axis1_norm_max_startcut.beta(1)],xp);
plot(xp,yp)
plot(turn_angle_vel_mirror,M_axis1_norm_max_startcut,'.')

yp = polyval([stats_M_axis2_norm_min_cutstop.beta(2) stats_M_axis2_norm_min_cutstop.beta(1)],xp);
plot(xp,yp,'r')
plot(turn_angle_vel_mirror,M_axis2_norm_min_cutstop,'.r')

axis([0 180 -.045 .045])

% M_axis1normal_max & M_axis2normal_min
subplot(3,2,4)
hold on

yp = polyval([stats_M_axis1normal_norm_max_startcut.beta(2) stats_M_axis1normal_norm_max_startcut.beta(1)],xp);
plot(xp,yp)
plot(turn_angle_vel_mirror,M_axis1normal_norm_max_startcut,'.')

yp = polyval([stats_M_axis2normal_norm_min_cutstop.beta(2) stats_M_axis2normal_norm_min_cutstop.beta(1)],xp);
plot(xp,yp,'r')
plot(turn_angle_vel_mirror,M_axis2normal_norm_min_cutstop,'.r')

axis([0 180 -.045 .045])

% Myaw max
subplot(3,2,6)
hold on

yp = polyval([stats_Myaw_norm_max_startstop.beta(2) stats_Myaw_norm_max_startstop.beta(1)],xp);
plot(xp,yp)
plot(turn_angle_vel_mirror,Myaw_norm_max_startstop,'.')

axis([0 180 0 .09])
xlabel('turn angle')

