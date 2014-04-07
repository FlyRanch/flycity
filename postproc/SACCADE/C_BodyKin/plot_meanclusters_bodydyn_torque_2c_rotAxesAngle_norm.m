%% MEAN clusters
angle_pre_min = t_start;
angle_pre_max = t_stop;
mean_min = t_start;
mean_max = t_stop;

% M_R_accel
subplot(3,3,1)
angle_pre = t_c1(:);
angle_post = M_R_accel_norm_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = M_R_accel_norm_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

% M_L_accel
subplot(3,3,4)
angle_pre = t_c1(:);
angle_post = M_L_accel_norm_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = M_L_accel_norm_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

% M_accel axis angle
subplot(3,3,7)
angle_pre = t_c1(:);
angle_post = Maccel_axis_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_circmean_NOextendedsection_NOdp_min90to270

angle_pre = t_c2(:);
angle_post = Maccel_axis_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_circmean_NOextendedsection_NOdp_min90to270


% M_R_damp
subplot(3,3,2)
angle_pre = t_c1(:);
angle_post = M_R_damp_norm_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = M_R_damp_norm_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

% M_L_damp
subplot(3,3,5)
angle_pre = t_c1(:);
angle_post = M_L_damp_norm_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = M_L_damp_norm_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

% M_damp axis angle
subplot(3,3,8)
angle_pre = t_c1(:);
angle_post = Mdamp_axis_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_circmean_NOextendedsection_NOdp_min90to270

angle_pre = t_c2(:);
angle_post = Mdamp_axis_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_circmean_NOextendedsection_NOdp_min90to270

% M_R
subplot(3,3,3)
angle_pre = t_c1(:);
angle_post = M_R_norm_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = M_R_norm_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

% M_L
subplot(3,3,6)
angle_pre = t_c1(:);
angle_post = M_L_norm_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = M_L_norm_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

% M axis angle
subplot(3,3,9)
angle_pre = t_c1(:);
angle_post = M_axis_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_circmean_NOextendedsection_NOdp_min90to270

angle_pre = t_c2(:);
angle_post = M_axis_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_circmean_NOextendedsection_NOdp_min90to270
