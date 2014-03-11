% figure
% subplot(3,3,1)
% title('mean attitudes')
% hold on

% left turn
% plotcolor = 'b'
plotcolor = [.5 .5 .5];
angle_pre = turn_angle_vel_mirror;
angle_pre_min = 0;
angle_pre_max = 180;
mean_min = 0;
mean_max = 180;

subplot(3,3,2)
angle_post = calc_value(droll_mirror,n_Amax);
plot_angle_pre_post_circmean_extendedsection_meanlim_turn
% xlabel('initial heading','fontsize',10) 
% axis equal
ylabel('droll @Amax','fontsize',10) 
set(gca,'ylim',[-15 45])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
set(gca,'xlim',[angle_pre_min angle_pre_max])
set(gca,'XTick',[-180:45:180],'fontsize',8) 
% grid on
% axis equal

subplot(3,3,5)
angle_post = calc_value(dpitch,n_Amax);
plot_angle_pre_post_circmean_extendedsection_meanlim
% plot_angle_pre_post_circmean_mirrorsection
% xlabel('initial heading','fontsize',10) 
% axis equal
ylabel('dpitch @Amax','fontsize',10) 
set(gca,'ylim',[-15 45])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
set(gca,'xlim',[angle_pre_min angle_pre_max])
set(gca,'XTick',[-180:45:180],'fontsize',8) 
% grid on
% axis equal

subplot(3,3,8)
angle_post = calc_value(yaw_mirror,n_Amax);
plot_angle_pre_post_circmean_extendedsection_meanlim
% plot_angle_pre_post_circmean_doublemirrorsection
% xlabel('initial heading','fontsize',10) 
% axis equal
xlabel('turn angle','fontsize',10) 
ylabel('yaw @Amax','fontsize',10) 
set(gca,'ylim',[0 60])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
set(gca,'xlim',[angle_pre_min angle_pre_max])
set(gca,'XTick',[-180:45:180],'fontsize',8) 
% grid on
% axis equal

