%% attitude @ Amax
n_now = n_Amax;

figure
subplot(3,3,1)
title('attitudes @ Amax','fontsize',10) 
hold on

% left turn
plotcolor = 'b';
angle_pre = turn_angle_Amean2Vpre_mirror;
angle_pre_min = 0;
angle_pre_max = 180;
mean_min = 0;
mean_max = 180;

subplot(3,3,1)
angle_post = calc_value(roll_mirror,n_now);
plot_angle_pre_post_circmean_extendedsection_meanlim
% xlabel('initial heading','fontsize',10) 
ylabel('roll','fontsize',10) 
set(gca,'ylim',[-90 90])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
grid on

subplot(3,3,4)
angle_post = calc_value(dpitch,n_now);
% plot_angle_pre_post_circmean_mirrorsection
plot_angle_pre_post_circmean_extendedsection_meanlim
set(gca,'ylim',[-90 90])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
% xlabel('initial heading','fontsize',10) 
ylabel('pitch','fontsize',10) 
grid on
% axis equal

subplot(3,3,7)
angle_post = calc_value(yaw_mirror,n_now);
% plot_angle_pre_post_circmean_doublemirrorsection
plot_angle_pre_post_circmean_extendedsection_meanlim
set(gca,'ylim',[-90 90])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
% xlabel('initial heading','fontsize',10) 
ylabel('yaw','fontsize',10) 
grid on
% axis equal

% right turn
plotcolor = 'r';
angle_pre = -turn_angle_Amean2Vpre_mirror;
angle_pre_min = -180;
angle_pre_max = 180;
mean_min = -180;
mean_max = 0;

subplot(3,3,1)
angle_post = -calc_value(roll_mirror,n_now);
plot_angle_pre_post_circmean_extendedsection_meanlim

subplot(3,3,4)
angle_post = calc_value(dpitch,n_now);
plot_angle_pre_post_circmean_extendedsection_meanlim
% set(gca,'xlim',[angle_pre_min angle_pre_max])

subplot(3,3,7)
angle_post = -calc_value(yaw_mirror,n_now);
plot_angle_pre_post_circmean_extendedsection_meanlim
% set(gca,'xlim',[angle_pre_min angle_pre_max])
xlabel('turn angle','fontsize',10) 

%% Attitude Accel max & min

subplot(3,3,2)
title('attitude accelerations','fontsize',10) 
hold on

% left turn
angle_pre = turn_angle_Amean2Vpre_mirror;
angle_pre_min = 0;
angle_pre_max = 180;
mean_min = 0;
mean_max = 180;

subplot(3,3,2)
plotcolor = 'b';
angle_post = roll_dot_dot_mirror_max;
plot_angle_pre_post_mean_extendedsection_meanlim
plotcolor = 'r';
angle_post = roll_dot_dot_mirror_min;
plot_angle_pre_post_mean_extendedsection_meanlim
set(gca,'ylim',[-2e5 2e5])
set(gca,'YTick',[-2e5:2e5:2e5],'fontsize',8) 
% xlabel('initial heading','fontsize',10) 
% axis equal
ylabel('roll accel','fontsize',10) 
grid on

subplot(3,3,5)
plotcolor = 'b';
angle_post = pitch_dot_dot_max;
plot_angle_pre_post_mean_extendedsection_meanlim
plotcolor = 'r';
angle_post = pitch_dot_dot_min;
plot_angle_pre_post_mean_extendedsection_meanlim
set(gca,'ylim',[-2e5 2e5])
set(gca,'YTick',[-2e5:2e5:2e5],'fontsize',8) 
% xlabel('initial heading','fontsize',10) 
% axis equal
ylabel('pitch accel','fontsize',10) 
grid on
% axis equal

subplot(3,3,8)
plotcolor = 'b';
angle_post = yaw_dot_dot_mirror_max;
plot_angle_pre_post_mean_extendedsection_meanlim
plotcolor = 'r';
angle_post = yaw_dot_dot_mirror_min;
plot_angle_pre_post_mean_extendedsection_meanlim
set(gca,'ylim',[-2e5 2e5])
set(gca,'YTick',[-2e5:2e5:2e5],'fontsize',8) 
% axis equal
ylabel('yaw accel','fontsize',10) 
grid on
% axis equal
% set(gca,'xlim',[angle_pre_min angle_pre_max])
xlabel('turn angle','fontsize',10) 

angle_pre = -turn_angle_Amean2Vpre_mirror;
angle_pre_min = -180;
angle_pre_max = 180;
mean_min = -180;
mean_max = 0;

subplot(3,3,2)
plotcolor = 'b';
angle_post = roll_dot_dot_mirror_max;
plot_angle_pre_post_mean_extendedsection_meanlim
plotcolor = 'r';
angle_post = roll_dot_dot_mirror_min;
plot_angle_pre_post_mean_extendedsection_meanlim

subplot(3,3,5)
plotcolor = 'b';
angle_post = pitch_dot_dot_max;
plot_angle_pre_post_mean_extendedsection_meanlim
plotcolor = 'r';
angle_post = pitch_dot_dot_min;
plot_angle_pre_post_mean_extendedsection_meanlim

subplot(3,3,8)
plotcolor = 'b';
angle_post = yaw_dot_dot_mirror_max;
plot_angle_pre_post_mean_extendedsection_meanlim
plotcolor = 'r';
angle_post = yaw_dot_dot_mirror_min;
plot_angle_pre_post_mean_extendedsection_meanlim

%% # of attitude accel wingbeats t(max)-t(min)
subplot(3,3,3)
title('manuever time')
hold on

% left turn
plotcolor = 'b';
angle_pre = turn_angle_Amean2Vpre_mirror;
angle_pre_min = 0;
angle_pre_max = 180;
mean_min = 0;
mean_max = 180;


subplot(3,3,3)
angle_post = abs(NwbSteady_RollAccel_max_min);
plot_angle_pre_post_circmean_extendedsection_meanlim
set(gca,'ylim',[0 10])
set(gca,'YTick',[0:5:10],'fontsize',8) 
ylabel('roll wingbeats','fontsize',10) 
grid on

subplot(3,3,6)
angle_post = abs(NwbSteady_PitchAccel_max_min);
plot_angle_pre_post_circmean_extendedsection_meanlim
set(gca,'ylim',[0 10])
set(gca,'YTick',[0:5:10],'fontsize',8) 
ylabel('pitch wingbeats','fontsize',10) 
grid on

subplot(3,3,9)
angle_post = abs(NwbSteady_YawAccel_max_min);
plot_angle_pre_post_circmean_extendedsection_meanlim
set(gca,'ylim',[0 10])
set(gca,'YTick',[0:5:10],'fontsize',8) 
ylabel('yaw wingbeats','fontsize',10) 
grid on


% right turn
plotcolor = 'r';
angle_pre = -turn_angle_Amean2Vpre_mirror;
angle_pre_min = -180;
angle_pre_max = 180;
mean_min = -180;
mean_max = 0;

subplot(3,3,3)
angle_post = abs(NwbSteady_RollAccel_max_min);
plot_angle_pre_post_circmean_extendedsection_meanlim

subplot(3,3,6)
angle_post = abs(NwbSteady_PitchAccel_max_min);
plot_angle_pre_post_circmean_extendedsection_meanlim

subplot(3,3,9)
angle_post = abs(NwbSteady_YawAccel_max_min);
plot_angle_pre_post_circmean_extendedsection_meanlim
xlabel('turn angle','fontsize',10) 






