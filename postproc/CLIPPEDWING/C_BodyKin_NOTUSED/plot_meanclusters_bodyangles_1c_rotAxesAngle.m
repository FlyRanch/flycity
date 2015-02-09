%% MEAN clusters
angle_pre_min = t_start;
angle_pre_max = t_stop;
mean_min = t_start;
mean_max = t_stop;

% drot_R
subplot(3,3,1)
angle_pre = t_plot(:);
angle_post = drot_R_plot(:);
plotcolor = [.5 .5 .5];
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

% drot_L
subplot(3,3,4)
angle_pre = t_plot(:);
angle_post = drot_L_plot(:);
plotcolor = [.5 .5 .5];
plot_timelines_dt_circmean_NOextendedsection_meanlim_NOdp

% rotation axis angle
subplot(3,3,7)
angle_pre = t_plot(:);
angle_post = rot_axis_pos_plot(:);
plotcolor = [.5 .5 .5];
plot_timelines_dt_circmean_NOextendedsection_NOdp_min90to270

% rot_R dot
subplot(3,3,2)
angle_pre = t_plot(:);
angle_post = rot_dot_R_plot(:)/freq_steady;
plotcolor = [.5 .5 .5];
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% rot_L dot
subplot(3,3,5)
angle_pre = t_plot(:);
angle_post = rot_dot_L_plot(:)/freq_steady;
plotcolor = [.5 .5 .5];
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% rotation axis angle
subplot(3,3,8)
angle_pre = t_plot(:);
angle_post = rot_axis_vel_plot(:);
plotcolor = [.5 .5 .5];
plot_timelines_dt_circmean_NOextendedsection_NOdp_min90to270

% rot_R dot_dot
subplot(3,3,3)
angle_pre = t_plot(:);
angle_post = rot_dot_dot_R_plot(:)/freq_steady^2;
plotcolor = [.5 .5 .5];
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% rot_L dot_dot
subplot(3,3,6)
angle_pre = t_plot(:);
angle_post = rot_dot_dot_L_plot(:)/freq_steady^2;
plotcolor = [.5 .5 .5];
plot_timelines_dt_mean_NOextendedsection_meanlim_NOdp

% rotation axis angle
subplot(3,3,9)
angle_pre = t_plot(:);
angle_post = rot_axis_accel_plot(:);
plotcolor = [.5 .5 .5];
plot_timelines_dt_circmean_NOextendedsection_NOdp_min90to270


