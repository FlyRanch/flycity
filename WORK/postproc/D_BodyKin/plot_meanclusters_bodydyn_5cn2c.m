%% MEAN clusters
angle_pre_min = -.02;
angle_pre_max = .04;
mean_min = -.02;
mean_max = .04;

% heading
subplot(3,3,1)
angle_pre = t_c1(:);
angle_post = stim_angle_vel_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = stim_angle_vel_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c3(:);
angle_post = stim_angle_vel_c3(:);
plotcolor = cmap_plot(subset_mid(3),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c4(:);
angle_post = stim_angle_vel_c4(:);
plotcolor = cmap_plot(subset_mid(4),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c5(:);
angle_post = stim_angle_vel_c5(:);
plotcolor = cmap_plot(subset_mid(5),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

% dV
subplot(3,3,4)
angle_pre = t_c1(:);
angle_post = dV_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = dV_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c3(:);
angle_post = dV_c3(:);
plotcolor = cmap_plot(subset_mid(3),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c4(:);
angle_post = dV_c4(:);
plotcolor = cmap_plot(subset_mid(4),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c5(:);
angle_post = dV_c5(:);
plotcolor = cmap_plot(subset_mid(5),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

% Fhor
subplot(3,3,7)
angle_pre = t_c1(:);
angle_post = stim_angle_accel_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = stim_angle_accel_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c3(:);
angle_post = stim_angle_accel_c3(:);
plotcolor = cmap_plot(subset_mid(3),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c4(:);
angle_post = stim_angle_accel_c4(:);
plotcolor = cmap_plot(subset_mid(4),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c5(:);
angle_post = stim_angle_accel_c5(:);
plotcolor = cmap_plot(subset_mid(5),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

% Force
subplot(3,3,2)
angle_pre = t_c1(:);
angle_post = F_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = F_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c3(:);
angle_post = F_c3(:);
plotcolor = cmap_plot(subset_mid(3),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c4(:);
angle_post = F_c4(:);
plotcolor = cmap_plot(subset_mid(4),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c5(:);
angle_post = F_c5(:);
plotcolor = cmap_plot(subset_mid(5),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

%% mean & ci ALL
subplot(3,3,5)
angle_pre = t_plot(:);
angle_post = Fsp_roll_plot(:);
plotcolor = [.5 .5 .5];
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

subplot(3,3,8)
angle_pre = t_plot(:);
angle_post = Fsp_pitch_plot(:);
plotcolor = [.5 .5 .5];
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

