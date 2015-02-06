figure
% subplot(2,2,1)
% title('@ max roll accel')
% hold on
% subplot(2,2,2)
% title('@ max yaw accel')
% hold on
% subplot(2,2,3)
% title('@ max pitch accel')
% hold on

%% left turn
angle_pre = stim_angle_yaw_mirror_pre;
angle_pre_min = -225;
angle_pre_max = 180;
mean_min = -180;
mean_max = 0;

% stroke
subplot(2,2,1)
plotcolor = 'b'
angle_post = stroke_Fenhance_WBmean_max;
plot_angle_pre_post_circmean_extendedsection_meanlim
ylabel('Dstroke','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-30 30])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-30 0 30],'fontsize',8) 
grid on

% pitch
subplot(2,2,2)
plotcolor = 'b'
angle_post = pitch_Fenhance_WBmean_max;
plot_angle_pre_post_circmean_extendedsection_meanlim
ylabel('Dpitch','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-30 30])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-30 0 30],'fontsize',8) 
grid on

% dev
subplot(2,2,3)
plotcolor = 'b'
angle_post = dev_Fenhance_WBmean_max;
plot_angle_pre_post_circmean_extendedsection_meanlim
ylabel('Ddev','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-30 30])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-30 0 30],'fontsize',8) 
grid on

% freq
subplot(2,2,4)
plotcolor = 'b'
angle_post = freq_Fenhance_WBmean_max;
plot_angle_pre_post_circmean_extendedsection_meanlim
ylabel('Dfreq','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[0 60])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[0 30 60],'fontsize',8) 
grid on

%% right turn
angle_pre = -stim_angle_yaw_mirror_pre;
angle_pre_min = -180;
angle_pre_max = 180;
mean_min = 0;
mean_max = 180;

% stroke
subplot(2,2,1)
plotcolor = 'b'
angle_post = stroke_Fenhance_WBmean_max;
plot_angle_pre_post_circmean_extendedsection_meanlim

% pitch
subplot(2,2,2)
plotcolor = 'b'
angle_post = pitch_Fenhance_WBmean_max;
plot_angle_pre_post_circmean_extendedsection_meanlim

% dev
subplot(2,2,3)
plotcolor = 'b'
angle_post = dev_Fenhance_WBmean_max;
plot_angle_pre_post_circmean_extendedsection_meanlim

% freq
subplot(2,2,4)
plotcolor = 'b'
angle_post = freq_Fenhance_WBmean_max;
plot_angle_pre_post_circmean_extendedsection_meanlim

