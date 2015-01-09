%% MEAN clusters
angle_pre_min = t_start;
angle_pre_max = t_stop;
mean_min = t_start;
mean_max = t_stop;

% Mroll_accel
subplot(3,3,1)
angle_pre = t_c1(:);
angle_post = Mroll_accel_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = Mroll_accel_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% angle_pre = t_c3(:);
% angle_post = Mroll_accel_c3(:);
% plotcolor = cmap_plot(subset_mid(3),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c4(:);
% angle_post = Mroll_accel_c4(:);
% plotcolor = cmap_plot(subset_mid(4),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c5(:);
% angle_post = Mroll_accel_c5(:);
% plotcolor = cmap_plot(subset_mid(5),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% Mpitch_accel
subplot(3,3,4)
angle_pre = t_c1(:);
angle_post = Mpitch_accel_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = Mpitch_accel_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% angle_pre = t_c3(:);
% angle_post = Mpitch_accel_c3(:);
% plotcolor = cmap_plot(subset_mid(3),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c4(:);
% angle_post = Mpitch_accel_c4(:);
% plotcolor = cmap_plot(subset_mid(4),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c5(:);
% angle_post = Mpitch_accel_c5(:);
% plotcolor = cmap_plot(subset_mid(5),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% Myaw_accel
subplot(3,3,7)
angle_pre = t_c1(:);
angle_post = Myaw_accel_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = Myaw_accel_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% angle_pre = t_c3(:);
% angle_post = Myaw_accel_c3(:);
% plotcolor = cmap_plot(subset_mid(3),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c4(:);
% angle_post = Myaw_accel_c4(:);
% plotcolor = cmap_plot(subset_mid(4),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c5(:);
% angle_post = Myaw_accel_c5(:);
% plotcolor = cmap_plot(subset_mid(5),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% Mroll_damp
subplot(3,3,2)
angle_pre = t_c1(:);
angle_post = Mroll_damp_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = Mroll_damp_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% angle_pre = t_c3(:);
% angle_post = Mroll_damp_c3(:);
% plotcolor = cmap_plot(subset_mid(3),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c4(:);
% angle_post = Mroll_damp_c4(:);
% plotcolor = cmap_plot(subset_mid(4),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c5(:);
% angle_post = Mroll_damp_c5(:);
% plotcolor = cmap_plot(subset_mid(5),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% Mpitch_damp
subplot(3,3,5)
angle_pre = t_c1(:);
angle_post = Mpitch_damp_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = Mpitch_damp_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% angle_pre = t_c3(:);
% angle_post = Mpitch_damp_c3(:);
% plotcolor = cmap_plot(subset_mid(3),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c4(:);
% angle_post = Mpitch_damp_c4(:);
% plotcolor = cmap_plot(subset_mid(4),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c5(:);
% angle_post = Mpitch_damp_c5(:);
% plotcolor = cmap_plot(subset_mid(5),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% Myaw_damp
subplot(3,3,8)
angle_pre = t_c1(:);
angle_post = Myaw_damp_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = Myaw_damp_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% angle_pre = t_c3(:);
% angle_post = Myaw_damp_c3(:);
% plotcolor = cmap_plot(subset_mid(3),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c4(:);
% angle_post = Myaw_damp_c4(:);
% plotcolor = cmap_plot(subset_mid(4),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c5(:);
% angle_post = Myaw_damp_c5(:);
% plotcolor = cmap_plot(subset_mid(5),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% Mroll
subplot(3,3,3)
angle_pre = t_c1(:);
angle_post = Mroll_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = Mroll_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% angle_pre = t_c3(:);
% angle_post = Mroll_c3(:);
% plotcolor = cmap_plot(subset_mid(3),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c4(:);
% angle_post = Mroll_c4(:);
% plotcolor = cmap_plot(subset_mid(4),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c5(:);
% angle_post = Mroll_c5(:);
% plotcolor = cmap_plot(subset_mid(5),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% Mpitch
subplot(3,3,6)
angle_pre = t_c1(:);
angle_post = Mpitch_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = Mpitch_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% angle_pre = t_c3(:);
% angle_post = Mpitch_c3(:);
% plotcolor = cmap_plot(subset_mid(3),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c4(:);
% angle_post = Mpitch_c4(:);
% plotcolor = cmap_plot(subset_mid(4),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c5(:);
% angle_post = Mpitch_c5(:);
% plotcolor = cmap_plot(subset_mid(5),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% Myaw
subplot(3,3,9)
angle_pre = t_c1(:);
angle_post = Myaw_c1(:);
plotcolor = cmap_plot(subset_mid(1),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

angle_pre = t_c2(:);
angle_post = Myaw_c2(:);
plotcolor = cmap_plot(subset_mid(2),:);
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% angle_pre = t_c3(:);
% angle_post = Myaw_c3(:);
% plotcolor = cmap_plot(subset_mid(3),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c4(:);
% angle_post = Myaw_c4(:);
% plotcolor = cmap_plot(subset_mid(4),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c5(:);
% angle_post = Myaw_c5(:);
% plotcolor = cmap_plot(subset_mid(5),:);
% plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp
