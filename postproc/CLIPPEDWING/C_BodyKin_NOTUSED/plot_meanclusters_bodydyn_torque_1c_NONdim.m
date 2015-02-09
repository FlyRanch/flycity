%% MEAN clusters
angle_pre_min = t_start;
angle_pre_max = t_stop;
mean_min = t_start;
mean_max = t_stop;

% Mroll_accel_norm
subplot(3,3,1)
angle_pre = t_plot(:);
angle_post = Mroll_accel_norm(:);
plotcolor = [.5 .5 .5];
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% Mpitch_accel_norm
subplot(3,3,4)
angle_pre = t_plot(:);
angle_post = Mpitch_accel_norm(:);
plotcolor = [.5 .5 .5];
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp


% Myaw_accel_norm
subplot(3,3,7)
angle_pre = t_plot(:);
angle_post = Myaw_accel_norm(:);
plotcolor = [.5 .5 .5];
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp


% Mroll_damp_norm
subplot(3,3,2)
angle_pre = t_plot(:);
angle_post = Mroll_damp_norm(:);
plotcolor = [.5 .5 .5];
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp


% Mpitch_damp_norm
subplot(3,3,5)
angle_pre = t_plot(:);
angle_post = Mpitch_damp_norm(:);
plotcolor = [.5 .5 .5];
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp


% Myaw_damp_norm
subplot(3,3,8)
angle_pre = t_plot(:);
angle_post = Myaw_damp_norm(:);
plotcolor = [.5 .5 .5];
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp


% Mroll_norm
subplot(3,3,3)
angle_pre = t_plot(:);
angle_post = Mroll_norm(:);
plotcolor = [.5 .5 .5];
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp


% Mpitch_norm
subplot(3,3,6)
angle_pre = t_plot(:);
angle_post = Mpitch_norm(:);
plotcolor = [.5 .5 .5];
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp


% Myaw_norm
subplot(3,3,9)
angle_pre = t_plot(:);
angle_post = Myaw_norm(:);
plotcolor = [.5 .5 .5];
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp


