%% MEAN clusters
angle_pre_min = t_start;
angle_pre_max = t_stop;
mean_min = t_start;
mean_max = t_stop;

% droll
subplot(3,3,1)
angle_pre = t_plot(:);
angle_post = droll_plot(:);
plotcolor = [.5 .5 .5];
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

% pitch
subplot(3,3,4)
angle_pre = t_plot(:);
angle_post = dpitch_plot(:);
plotcolor = [.5 .5 .5];
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

% yaw
subplot(3,3,7)
angle_pre = t_plot(:);
angle_post = yaw_plot(:);
plotcolor = [.5 .5 .5];
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

% roll dot
subplot(3,3,2)
angle_pre = t_plot(:);
angle_post = roll_dot_plot(:)/freq_steady;
plotcolor = [.5 .5 .5];
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% pitch_dot
subplot(3,3,5)
angle_pre = t_plot(:);
angle_post = pitch_dot_plot(:)/freq_steady;
plotcolor = [.5 .5 .5];
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% yaw_dot
subplot(3,3,8)
angle_pre = t_plot(:);
angle_post = yaw_dot_plot(:)/freq_steady;
plotcolor = [.5 .5 .5];
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% roll dot_dot
subplot(3,3,3)
angle_pre = t_plot(:);
angle_post = roll_dot_dot_plot(:)/freq_steady^2;
plotcolor = [.5 .5 .5];
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% pitch_dot_dot
subplot(3,3,6)
angle_pre = t_plot(:);
angle_post = pitch_dot_dot_plot(:)/freq_steady^2;
plotcolor = [.5 .5 .5];
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% yaw_dot_dot
subplot(3,3,9)
angle_pre = t_plot(:);
angle_post = yaw_dot_dot_plot(:)/freq_steady^2;
plotcolor = [.5 .5 .5];
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp


