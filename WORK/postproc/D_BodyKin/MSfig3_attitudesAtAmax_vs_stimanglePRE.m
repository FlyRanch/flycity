% figure
% subplot(3,3,1)
% title('mean attitudes')
% hold on

% left turn
% plotcolor = 'b'
plotcolor = [.5 .5 .5];
angle_pre = stim_angle_yaw_mirror_pre;
angle_pre_min = -180;
angle_pre_max = 0;
mean_min = -180;
mean_max = 0;

subplot(3,3,2)
angle_post = calc_value(roll_mirror,n_Amax);
plot_angle_pre_post_circmean_extendedsection_meanlim
% xlabel('initial heading','fontsize',10) 
axis equal
ylabel('roll','fontsize',10) 
set(gca,'ylim',[-45 90])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
set(gca,'xlim',[-180 0])
set(gca,'XTick',[-180:45:180],'fontsize',8) 
% grid on
% axis equal

subplot(3,3,5)
angle_post = calc_value(dpitch,n_Amax);
plot_angle_pre_post_circmean_extendedsection_meanlim
% plot_angle_pre_post_circmean_mirrorsection
% xlabel('initial heading','fontsize',10) 
axis equal
ylabel('pitch','fontsize',10) 
set(gca,'ylim',[-45 90])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
set(gca,'xlim',[-180 0])
set(gca,'XTick',[-180:45:180],'fontsize',8) 
% grid on
% axis equal

subplot(3,3,8)
angle_post = calc_value(yaw_mirror,n_Amax);
plot_angle_pre_post_circmean_extendedsection_meanlim
% plot_angle_pre_post_circmean_doublemirrorsection
% xlabel('initial heading','fontsize',10) 
axis equal
xlabel('stimulus angle','fontsize',10) 
ylabel('yaw','fontsize',10) 
set(gca,'ylim',[-45 90])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
set(gca,'xlim',[-180 0])
set(gca,'XTick',[-180:45:180],'fontsize',8) 
% grid on
% axis equal

