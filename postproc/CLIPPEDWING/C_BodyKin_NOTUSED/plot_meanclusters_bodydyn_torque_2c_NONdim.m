%% MEAN clusters
angle_pre_min = t_start;
angle_pre_max = t_stop;
mean_min = t_start;
mean_max = t_stop;

% Mroll_accel_norm
subplot(3,3,1)
angle_pre = t_c1(:);
angle_post = Mroll_accel_norm_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = Mroll_accel_norm_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% angle_pre = t_c3(:);
% angle_post = Mroll_accel_norm_c3(:);
% plotcolor = cmap_plot(subset_mid(3),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c4(:);
% angle_post = Mroll_accel_norm_c4(:);
% plotcolor = cmap_plot(subset_mid(4),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c5(:);
% angle_post = Mroll_accel_norm_c5(:);
% plotcolor = cmap_plot(subset_mid(5),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% Mpitch_accel_norm
subplot(3,3,4)
angle_pre = t_c1(:);
angle_post = Mpitch_accel_norm_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = Mpitch_accel_norm_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% angle_pre = t_c3(:);
% angle_post = Mpitch_accel_norm_c3(:);
% plotcolor = cmap_plot(subset_mid(3),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c4(:);
% angle_post = Mpitch_accel_norm_c4(:);
% plotcolor = cmap_plot(subset_mid(4),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c5(:);
% angle_post = Mpitch_accel_norm_c5(:);
% plotcolor = cmap_plot(subset_mid(5),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% Myaw_accel_norm
subplot(3,3,7)
angle_pre = t_c1(:);
angle_post = Myaw_accel_norm_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = Myaw_accel_norm_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% angle_pre = t_c3(:);
% angle_post = Myaw_accel_norm_c3(:);
% plotcolor = cmap_plot(subset_mid(3),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c4(:);
% angle_post = Myaw_accel_norm_c4(:);
% plotcolor = cmap_plot(subset_mid(4),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c5(:);
% angle_post = Myaw_accel_norm_c5(:);
% plotcolor = cmap_plot(subset_mid(5),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% Mroll_damp_norm
subplot(3,3,2)
angle_pre = t_c1(:);
angle_post = Mroll_damp_norm_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = Mroll_damp_norm_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% angle_pre = t_c3(:);
% angle_post = Mroll_damp_norm_c3(:);
% plotcolor = cmap_plot(subset_mid(3),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c4(:);
% angle_post = Mroll_damp_norm_c4(:);
% plotcolor = cmap_plot(subset_mid(4),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c5(:);
% angle_post = Mroll_damp_norm_c5(:);
% plotcolor = cmap_plot(subset_mid(5),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% Mpitch_damp_norm
subplot(3,3,5)
angle_pre = t_c1(:);
angle_post = Mpitch_damp_norm_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = Mpitch_damp_norm_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% angle_pre = t_c3(:);
% angle_post = Mpitch_damp_norm_c3(:);
% plotcolor = cmap_plot(subset_mid(3),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c4(:);
% angle_post = Mpitch_damp_norm_c4(:);
% plotcolor = cmap_plot(subset_mid(4),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c5(:);
% angle_post = Mpitch_damp_norm_c5(:);
% plotcolor = cmap_plot(subset_mid(5),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% Myaw_damp_norm
subplot(3,3,8)
angle_pre = t_c1(:);
angle_post = Myaw_damp_norm_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = Myaw_damp_norm_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% angle_pre = t_c3(:);
% angle_post = Myaw_damp_norm_c3(:);
% plotcolor = cmap_plot(subset_mid(3),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c4(:);
% angle_post = Myaw_damp_norm_c4(:);
% plotcolor = cmap_plot(subset_mid(4),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c5(:);
% angle_post = Myaw_damp_norm_c5(:);
% plotcolor = cmap_plot(subset_mid(5),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% Mroll_norm
subplot(3,3,3)
angle_pre = t_c1(:);
angle_post = Mroll_norm_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = Mroll_norm_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% angle_pre = t_c3(:);
% angle_post = Mroll_norm_c3(:);
% plotcolor = cmap_plot(subset_mid(3),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c4(:);
% angle_post = Mroll_norm_c4(:);
% plotcolor = cmap_plot(subset_mid(4),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c5(:);
% angle_post = Mroll_norm_c5(:);
% plotcolor = cmap_plot(subset_mid(5),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% Mpitch_norm
subplot(3,3,6)
angle_pre = t_c1(:);
angle_post = Mpitch_norm_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = Mpitch_norm_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% angle_pre = t_c3(:);
% angle_post = Mpitch_norm_c3(:);
% plotcolor = cmap_plot(subset_mid(3),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c4(:);
% angle_post = Mpitch_norm_c4(:);
% plotcolor = cmap_plot(subset_mid(4),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c5(:);
% angle_post = Mpitch_norm_c5(:);
% plotcolor = cmap_plot(subset_mid(5),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% Myaw_norm
subplot(3,3,9)
angle_pre = t_c1(:);
angle_post = Myaw_norm_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = Myaw_norm_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% angle_pre = t_c3(:);
% angle_post = Myaw_norm_c3(:);
% plotcolor = cmap_plot(subset_mid(3),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c4(:);
% angle_post = Myaw_norm_c4(:);
% plotcolor = cmap_plot(subset_mid(4),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c5(:);
% angle_post = Myaw_norm_c5(:);
% plotcolor = cmap_plot(subset_mid(5),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
