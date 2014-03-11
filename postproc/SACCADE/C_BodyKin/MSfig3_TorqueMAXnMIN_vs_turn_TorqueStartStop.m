% figure
% subplot(3,3,1)
% % title('mean attitudes')
% hold on

% left turn
% plotcolor = 'b'
plotcolor = [.5 .5 .5];
angle_pre = turn_angle_vel_mirror;
angle_pre_min = 0;
angle_pre_max = 180;
mean_min = 0;
mean_max = 180;

%% Max Torque
subplot(3,3,1)
angle_post = Mroll_norm_max_startstop;
% angle_post = Mroll_norm_max;
% angle_post = max(Mroll_norm)';
plot_angle_pre_post_mean_extendedsection_meanlim 
% xlabel('initial heading','fontsize',10) 
% axis equal
ylabel('Troll/Mgl','fontsize',10) 
% title('Max Torque','fontsize',10) 
set(gca,'ylim',[-.03 .03])
set(gca,'YTick',[-.03:.03:.03],'fontsize',8) 
set(gca,'xlim',[angle_pre_min angle_pre_max])
set(gca,'XTick',[-180:45:180],'fontsize',8) 
% grid on
% axis equal

subplot(3,3,4)
angle_post = Mpitch_norm_max_startstop;
% angle_post = Mpitch_norm_max;
% angle_post = max(Mpitch_norm)';
plot_angle_pre_post_mean_extendedsection_meanlim 
% xlabel('initial heading','fontsize',10) 
% axis equal
ylabel('Tpitch/Mgl','fontsize',10) 
set(gca,'ylim',[-.03 .03])
set(gca,'YTick',[-.03:.03:.03],'fontsize',8) 
set(gca,'xlim',[angle_pre_min angle_pre_max])
set(gca,'XTick',[-180:45:180],'fontsize',8) 
% grid on
% axis equal

subplot(3,3,7)
angle_post = Myaw_norm_max_startstop;
% angle_post = Myaw_norm_max;
% angle_post = max(Myaw_norm)';
plot_angle_pre_post_mean_extendedsection_meanlim 
xlabel('turn angle','fontsize',10) 
% axis equal
ylabel('Tyaw/Mgl','fontsize',10) 
set(gca,'ylim',[-0 .06])
set(gca,'YTick',[-.3:.3:.6],'fontsize',8) 
set(gca,'xlim',[angle_pre_min angle_pre_max])
set(gca,'XTick',[-180:45:180],'fontsize',8) 
% grid on
% axis equal

%% Min Torque
subplot(3,3,2)
angle_post = Mroll_norm_min_0stop;
% angle_post = Mroll_norm_min;
% angle_post = min(Mroll_norm)';
plot_angle_pre_post_mean_extendedsection_meanlim 
% xlabel('initial heading','fontsize',10) 
% axis equal
% ylabel('Troll/Mgl','fontsize',10) 
% title('Min Torque','fontsize',10) 
set(gca,'ylim',[-.03 .03])
set(gca,'YTick',[-.03:.03:.03],'fontsize',8) 
set(gca,'xlim',[angle_pre_min angle_pre_max])
set(gca,'XTick',[-180:45:180],'fontsize',8) 
% grid on
% axis equal

subplot(3,3,5)
angle_post = Mpitch_norm_min_0stop;
% angle_post = Mpitch_norm_min;
% angle_post = min(Mpitch_norm)';
plot_angle_pre_post_mean_extendedsection_meanlim 
% xlabel('initial heading','fontsize',10) 
% axis equal
% ylabel('Tpitch/Mgl','fontsize',10) 
set(gca,'ylim',[-.03 .03])
set(gca,'YTick',[-.03:.03:.03],'fontsize',8) 
set(gca,'xlim',[angle_pre_min angle_pre_max])
set(gca,'XTick',[-180:45:180],'fontsize',8) 
% grid on
% axis equal

subplot(3,3,8)
angle_post = Myaw_norm_min_0stop;
% angle_post = Myaw_norm_min;
% angle_post = min(Myaw_norm)';
plot_angle_pre_post_mean_extendedsection_meanlim 
xlabel('turn angle','fontsize',10) 
% axis equal
% ylabel('Tyaw/Mgl','fontsize',10) 
set(gca,'ylim',[-0 .06])
set(gca,'YTick',[-.3:.3:.6],'fontsize',8) 
set(gca,'xlim',[angle_pre_min angle_pre_max])
set(gca,'XTick',[-180:45:180],'fontsize',8) 
% grid on
% axis equal

%% Max Torque BLUE
plotcolor = [0 0 1];

subplot(3,3,3)
angle_post = Mroll_norm_max_startstop;
% angle_post = Mroll_norm_max;
% angle_post = max(Mroll_norm)';
plot_angle_pre_post_mean_extendedsection_meanlim 
% xlabel('initial heading','fontsize',10) 
% axis equal
% ylabel('Troll/Mgl','fontsize',10) 
% % title('Max Torque','fontsize',10) 
% set(gca,'ylim',[-0 .06])
% set(gca,'YTick',[-.03:.03:.03],'fontsize',8) 
% set(gca,'xlim',[angle_pre_min angle_pre_max])
% set(gca,'XTick',[-180:45:180],'fontsize',8) 
% grid on
% axis equal

subplot(3,3,6)
angle_post = Mpitch_norm_max_startstop;
% angle_post = Mpitch_norm_max;
% angle_post = max(Mpitch_norm)';
plot_angle_pre_post_mean_extendedsection_meanlim 
% xlabel('initial heading','fontsize',10) 
% axis equal
% ylabel('Tpitch/Mgl','fontsize',10) 
% set(gca,'ylim',[-.03 .03])
% set(gca,'YTick',[-.03:.03:.03],'fontsize',8) 
% set(gca,'xlim',[angle_pre_min angle_pre_max])
% set(gca,'XTick',[-180:45:180],'fontsize',8) 
% grid on
% axis equal

subplot(3,3,9)
angle_post = Myaw_norm_max_startstop;
% angle_post = Myaw_norm_max;
% angle_post = max(Myaw_norm)';
plot_angle_pre_post_mean_extendedsection_meanlim 
% xlabel('turn angle','fontsize',10) 
% axis equal
% ylabel('Tyaw/Mgl','fontsize',10) 
% set(gca,'ylim',[-0 .06])
% set(gca,'YTick',[-180:30:180],'fontsize',8) 
% set(gca,'xlim',[angle_pre_min angle_pre_max])
% set(gca,'XTick',[-180:45:180],'fontsize',8) 
% grid on
% axis equal

%% Min Torque RED
plotcolor = [1 0 0];

subplot(3,3,3)
angle_post = Mroll_norm_min_0stop;
% angle_post = Mroll_norm_min;
% angle_post = min(Mroll_norm)';
plot_angle_pre_post_mean_extendedsection_meanlim
% xlabel('initial heading','fontsize',10) 
% axis equal
% ylabel('Troll/Mgl','fontsize',10) 
% title('Max & Min Torque','fontsize',10) 
set(gca,'ylim',[-.03 .03])
set(gca,'YTick',[-.03:.03:.03],'fontsize',8) 
set(gca,'xlim',[angle_pre_min angle_pre_max])
set(gca,'XTick',[-180:45:180],'fontsize',8) 
% grid on
% axis equal

subplot(3,3,6)
angle_post = Mpitch_norm_min_0stop;
% angle_post = Mpitch_norm_min;
% angle_post = min(Mpitch_norm)';
plot_angle_pre_post_mean_extendedsection_meanlim 
% xlabel('initial heading','fontsize',10) 
% axis equal
% ylabel('Tpitch/Mgl','fontsize',10) 
set(gca,'ylim',[-.03 .03])
set(gca,'YTick',[-.03:.03:.03],'fontsize',8) 
set(gca,'xlim',[angle_pre_min angle_pre_max])
set(gca,'XTick',[-180:45:180],'fontsize',8) 
% grid on
% axis equal

subplot(3,3,9)
angle_post = Myaw_norm_min_0stop;
% angle_post = Myaw_norm_min;
% angle_post = min(Myaw_norm)';
plot_angle_pre_post_mean_extendedsection_meanlim 
xlabel('turn angle','fontsize',10) 
% axis equal
% ylabel('Tyaw/Mgl','fontsize',10) 
set(gca,'ylim',[-0 .06])
set(gca,'YTick',[-.03:.03:.06],'fontsize',8) 
set(gca,'xlim',[angle_pre_min angle_pre_max])
set(gca,'XTick',[-180:45:180],'fontsize',8) 
% grid on
% axis equal

