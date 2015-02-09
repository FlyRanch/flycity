% figure
% subplot(3,3,1)
% % title('mean attitudes')
% hold on

angle_pre = turn_angle_vel_mirror;
angle_pre_min = 0;
angle_pre_max = 180;
mean_min = 0;
mean_max = 180;

%% MAX & MIN

% Max Torque axis1 & Min Torque axis2
subplot(3,3,1)

plotcolor = [0 .5 1];
angle_post = M_axis1_norm_max_startcut;
plot_angle_pre_post_mean_extendedsection_meanlim 

plotcolor = [1 .5 0];
angle_post = M_axis2_norm_min_cutstop;
plot_angle_pre_post_mean_extendedsection_meanlim

set(gca,'ylim',[-.045 .045])
set(gca,'YTick',[-.045:.015:.09],'fontsize',8) 
set(gca,'xlim',[angle_pre_min angle_pre_max])
set(gca,'XTick',[-180:45:180],'fontsize',8) 

% Max Torque axis1normal & Min Torque axis2normal
subplot(3,3,2)

plotcolor = [0 .5 1];
angle_post = M_axis1normal_norm_max_startcut;
plot_angle_pre_post_mean_extendedsection_meanlim 

plotcolor = [1 .5 0];
angle_post = M_axis2normal_norm_min_cutstop;
plot_angle_pre_post_mean_extendedsection_meanlim

set(gca,'ylim',[-.045 .045])
set(gca,'YTick',[-.045:.015:.09],'fontsize',8) 
set(gca,'xlim',[angle_pre_min angle_pre_max])
set(gca,'XTick',[-180:45:180],'fontsize',8) 

% Max Torque yaw
subplot(3,3,3)
plotcolor = [.5 .5 .5];
angle_post = Myaw_norm_max_startstop;
plot_angle_pre_post_mean_extendedsection_meanlim 

set(gca,'ylim',[0 .09])
set(gca,'YTick',[-.045:.015:.09],'fontsize',8) 
set(gca,'xlim',[angle_pre_min angle_pre_max])
set(gca,'XTick',[-180:45:180],'fontsize',8) 


%% MEAN

% MEAN Torque axis1 & MEAN Torque axis2
subplot(3,3,4)

plotcolor = [0 .5 1];
angle_post = M_axis1_norm_mean_startcut;
plot_angle_pre_post_mean_extendedsection_meanlim 

plotcolor = [1 .5 0];
angle_post = M_axis2_norm_mean_cutstop;
plot_angle_pre_post_mean_extendedsection_meanlim

set(gca,'ylim',[-.03 .03])
set(gca,'YTick',[-.045:.015:.09],'fontsize',8) 
set(gca,'xlim',[angle_pre_min angle_pre_max])
set(gca,'XTick',[-180:45:180],'fontsize',8) 

%% mean Torque axis1normal
subplot(3,3,5)
plotcolor = [0 .5 1];
angle_post = M_axis1normal_norm_mean_startcut;
plot_angle_pre_post_mean_extendedsection_meanlim 

set(gca,'ylim',[-.03 .03])
set(gca,'YTick',[-.045:.015:.09],'fontsize',8) 
set(gca,'xlim',[angle_pre_min angle_pre_max])
set(gca,'XTick',[-180:45:180],'fontsize',8) 

%% mean Torque axis2normal
subplot(3,3,8)
plotcolor = [1 .5 0];
angle_post = M_axis2normal_norm_mean_cutstop;
plot_angle_pre_post_mean_extendedsection_meanlim 

set(gca,'ylim',[-.03 .03])
set(gca,'YTick',[-.045:.015:.09],'fontsize',8) 
set(gca,'xlim',[angle_pre_min angle_pre_max])
set(gca,'XTick',[-180:45:180],'fontsize',8) 

%% mean Torque yaw
subplot(3,3,6)
plotcolor = [.5 .5 .5];
angle_post = Myaw_norm_mean_startstop;
plot_angle_pre_post_mean_extendedsection_meanlim 

set(gca,'ylim',[-.0 .06])
set(gca,'YTick',[-.045:.015:.09],'fontsize',8) 
set(gca,'xlim',[angle_pre_min angle_pre_max])
set(gca,'XTick',[-180:45:180],'fontsize',8) 

