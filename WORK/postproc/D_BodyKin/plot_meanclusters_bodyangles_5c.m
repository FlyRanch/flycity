%% MEAN clusters
angle_pre_min = -.02;
angle_pre_max = .04;
mean_min = -.02;
mean_max = .04;

% roll
subplot(3,3,1)
angle_pre = t_c1(:);
angle_post = roll_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = roll_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c3(:);
angle_post = roll_c3(:);
plotcolor = cmap_plot(subset_mid(3),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c4(:);
angle_post = roll_c4(:);
plotcolor = cmap_plot(subset_mid(4),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c5(:);
angle_post = roll_c5(:);
plotcolor = cmap_plot(subset_mid(5),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

% pitch
subplot(3,3,4)
angle_pre = t_c1(:);
angle_post = pitch_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = pitch_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c3(:);
angle_post = pitch_c3(:);
plotcolor = cmap_plot(subset_mid(3),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c4(:);
angle_post = pitch_c4(:);
plotcolor = cmap_plot(subset_mid(4),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c5(:);
angle_post = pitch_c5(:);
plotcolor = cmap_plot(subset_mid(5),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

% yaw
subplot(3,3,7)
angle_pre = t_c1(:);
angle_post = yaw_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = yaw_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c3(:);
angle_post = yaw_c3(:);
plotcolor = cmap_plot(subset_mid(3),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c4(:);
angle_post = yaw_c4(:);
plotcolor = cmap_plot(subset_mid(4),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

angle_pre = t_c5(:);
angle_post = yaw_c5(:);
plotcolor = cmap_plot(subset_mid(5),:);
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

% roll dot
subplot(3,3,2)
angle_pre = t_c1(:);
angle_post = roll_dot_c1(:)/freq_steady;
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = roll_dot_c2(:)/freq_steady;
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c3(:);
angle_post = roll_dot_c3(:)/freq_steady;
plotcolor = cmap_plot(subset_mid(3),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c4(:);
angle_post = roll_dot_c4(:)/freq_steady;
plotcolor = cmap_plot(subset_mid(4),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c5(:);
angle_post = roll_dot_c5(:)/freq_steady;
plotcolor = cmap_plot(subset_mid(5),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% pitch_dot
subplot(3,3,5)
angle_pre = t_c1(:);
angle_post = pitch_dot_c1(:)/freq_steady;
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = pitch_dot_c2(:)/freq_steady;
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c3(:);
angle_post = pitch_dot_c3(:)/freq_steady;
plotcolor = cmap_plot(subset_mid(3),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c4(:);
angle_post = pitch_dot_c4(:)/freq_steady;
plotcolor = cmap_plot(subset_mid(4),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c5(:);
angle_post = pitch_dot_c5(:)/freq_steady;
plotcolor = cmap_plot(subset_mid(5),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% yaw_dot
subplot(3,3,8)
angle_pre = t_c1(:);
angle_post = yaw_dot_c1(:)/freq_steady;
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = yaw_dot_c2(:)/freq_steady;
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c3(:);
angle_post = yaw_dot_c3(:)/freq_steady;
plotcolor = cmap_plot(subset_mid(3),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c4(:);
angle_post = yaw_dot_c4(:)/freq_steady;
plotcolor = cmap_plot(subset_mid(4),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c5(:);
angle_post = yaw_dot_c5(:)/freq_steady;
plotcolor = cmap_plot(subset_mid(5),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% roll dot_dot
subplot(3,3,3)
angle_pre = t_c1(:);
angle_post = roll_dot_dot_c1(:)/freq_steady^2;
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = roll_dot_dot_c2(:)/freq_steady^2;
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c3(:);
angle_post = roll_dot_dot_c3(:)/freq_steady^2;
plotcolor = cmap_plot(subset_mid(3),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c4(:);
angle_post = roll_dot_dot_c4(:)/freq_steady^2;
plotcolor = cmap_plot(subset_mid(4),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c5(:);
angle_post = roll_dot_dot_c5(:)/freq_steady^2;
plotcolor = cmap_plot(subset_mid(5),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% pitch_dot_dot
subplot(3,3,6)
angle_pre = t_c1(:);
angle_post = pitch_dot_dot_c1(:)/freq_steady^2;
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = pitch_dot_dot_c2(:)/freq_steady^2;
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c3(:);
angle_post = pitch_dot_dot_c3(:)/freq_steady^2;
plotcolor = cmap_plot(subset_mid(3),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c4(:);
angle_post = pitch_dot_dot_c4(:)/freq_steady^2;
plotcolor = cmap_plot(subset_mid(4),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c5(:);
angle_post = pitch_dot_dot_c5(:)/freq_steady^2;
plotcolor = cmap_plot(subset_mid(5),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% yaw_dot_dot
subplot(3,3,9)
angle_pre = t_c1(:);
angle_post = yaw_dot_dot_c1(:)/freq_steady^2;
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = yaw_dot_dot_c2(:)/freq_steady^2;
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c3(:);
angle_post = yaw_dot_dot_c3(:)/freq_steady^2;
plotcolor = cmap_plot(subset_mid(3),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c4(:);
angle_post = yaw_dot_dot_c4(:)/freq_steady^2;
plotcolor = cmap_plot(subset_mid(4),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c5(:);
angle_post = yaw_dot_dot_c5(:)/freq_steady^2;
plotcolor = cmap_plot(subset_mid(5),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
