%% MEAN clusters
angle_pre_min = t_start;
angle_pre_max = t_stop;
mean_min = t_start;
mean_max = t_stop;

%% heading
subplot(3,3,1)
angle_pre = t_c1(:);
angle_post = stim_angle_vel_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = stim_angle_vel_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

% angle_pre = t_c3(:);
% angle_post = stim_angle_vel_c3(:);
% plotcolor = cmap_plot(subset_mid(3),:);
% plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c4(:);
% angle_post = stim_angle_vel_c4(:);
% plotcolor = cmap_plot(subset_mid(4),:);
% plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c5(:);
% angle_post = stim_angle_vel_c5(:);
% plotcolor = cmap_plot(subset_mid(5),:);
% plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

%% dV
subplot(3,3,4)
angle_pre = t_c1(:);
angle_post = dV_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = dV_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

% angle_pre = t_c3(:);
% angle_post = dV_c3(:);
% plotcolor = cmap_plot(subset_mid(3),:);
% plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c4(:);
% angle_post = dV_c4(:);
% plotcolor = cmap_plot(subset_mid(4),:);
% plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c5(:);
% angle_post = dV_c5(:);
% plotcolor = cmap_plot(subset_mid(5),:);
% plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

%% Fhor
subplot(3,3,7)
angle_pre = t_c1(:);
angle_post = stim_angle_accel_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
% plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp
plot_timelines_dt_circmean_NOextendedsection_NOdp_min90to270

angle_pre = t_c2(:);
angle_post = stim_angle_accel_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
% plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp
plot_timelines_dt_circmean_NOextendedsection_NOdp_min90to270

% angle_pre = t_c3(:);
% angle_post = stim_angle_accel_c3(:);
% plotcolor = cmap_plot(subset_mid(3),:);
% plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c4(:);
% angle_post = stim_angle_accel_c4(:);
% plotcolor = cmap_plot(subset_mid(4),:);
% plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c5(:);
% angle_post = stim_angle_accel_c5(:);
% plotcolor = cmap_plot(subset_mid(5),:);
% plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

%% Force
subplot(3,3,2)
angle_pre = t_c1(:);
angle_post = F_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = F_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

% % angle_pre = t_c3(:);
% % angle_post = F_c3(:);
% % plotcolor = cmap_plot(subset_mid(3),:);
% % plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp
% % 
% % angle_pre = t_c4(:);
% % angle_post = F_c4(:);
% % plotcolor = cmap_plot(subset_mid(4),:);
% % plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp
% % 
% % angle_pre = t_c5(:);
% % angle_post = F_c5(:);
% % plotcolor = cmap_plot(subset_mid(5),:);
% % plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

%% Force roll
subplot(3,3,5)
angle_pre = t_c1(:);
angle_post = Fsp_roll_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = Fsp_roll_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

% % angle_pre = t_c3(:);
% % angle_post = Fsp_roll_c3(:);
% % plotcolor = cmap_plot(subset_mid(3),:);
% % plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp
% % 
% % angle_pre = t_c4(:);
% % angle_post = Fsp_roll_c4(:);
% % plotcolor = cmap_plot(subset_mid(4),:);
% % plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp
% % 
% % angle_pre = t_c5(:);
% % angle_post = Fsp_roll_c5(:);
% % plotcolor = cmap_plot(subset_mid(5),:);
% % plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

%% Force pitch
subplot(3,3,8)
angle_pre = t_c1(:);
angle_post = Fsp_pitch_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = Fsp_pitch_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

% % angle_pre = t_c3(:);
% % angle_post = Fsp_pitch_c3(:);
% % plotcolor = cmap_plot(subset_mid(3),:);
% % plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp
% % 
% % angle_pre = t_c4(:);
% % angle_post = Fsp_pitch_c4(:);
% % plotcolor = cmap_plot(subset_mid(4),:);
% % plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp
% % 
% % angle_pre = t_c5(:);
% % angle_post = Fsp_pitch_c5(:);
% % plotcolor = cmap_plot(subset_mid(5),:);
% % plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

%% a/g
subplot(3,3,3)
angle_pre = t_c1(:);
angle_post = A_norm_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = A_norm_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

%% a/g hor
subplot(3,3,6)
angle_pre = t_c1(:);
angle_post = Ahor_norm_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = Ahor_norm_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

%% a/g ver
subplot(3,3,9)
angle_pre = t_c1(:);
angle_post = Az_norm_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = Az_norm_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

%% mean & ci ALL
% subplot(3,3,5)
% angle_pre = t_plot(:);
% angle_post = Fsp_roll_plot(:);
% plotcolor = [.5 .5 .5];
% plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp
% 
% subplot(3,3,8)
% angle_pre = t_plot(:);
% angle_post = Fsp_pitch_plot(:);
% plotcolor = [.5 .5 .5];
% plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

%% mean & quartiles ALL
% subplot(3,3,2)
% angle_pre = t_plot(:);
% angle_post = F_plot(:);
% plotcolor = [.5 .5 .5];
% plot_timelines_dt_circmean_NOext_meanlim_NOdp_quart
% 
% subplot(3,3,5)
% angle_pre = t_plot(:);
% angle_post = Fsp_roll_plot(:);
% plotcolor = [.5 .5 .5];
% plot_timelines_dt_circmean_NOext_meanlim_NOdp_quart
% 
% subplot(3,3,8)
% angle_pre = t_plot(:);
% angle_post = Fsp_pitch_plot(:);
% plotcolor = [.5 .5 .5];
% plot_timelines_dt_circmean_NOext_meanlim_NOdp_quart

