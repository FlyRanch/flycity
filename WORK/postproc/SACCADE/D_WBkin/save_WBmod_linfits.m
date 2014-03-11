
WBmod_linfit.fit_rollaccel_ddevmax_ds = fit_rollaccel_ddevmax_ds;

WBmod_linfit.fit_rollaccel_ddevmin_ds = fit_rollaccel_ddevmin_ds;

WBmod_linfit.fit_rollaccel_ddevmax_us = fit_rollaccel_ddevmax_us;

WBmod_linfit.fit_rollaccel_ddevmin_us = fit_rollaccel_ddevmin_us;

WBmod_linfit.fit_rollaccel_dAdev_ds = fit_rollaccel_dAdev_ds;

WBmod_linfit.fit_rollaccel_dAdev_us = fit_rollaccel_dAdev_us;

WBmod_linfit.fit_rollaccel_dstrokemax = fit_rollaccel_dstrokemax;

WBmod_linfit.fit_rollaccel_dstrokemin = fit_rollaccel_dstrokemin;

WBmod_linfit.fit_rollaccel_dAstroke = fit_rollaccel_dAstroke;

WBmod_linfit.fit_rollaccel_dpitchmid_ds = fit_rollaccel_dpitchmid_ds;

WBmod_linfit.fit_rollaccel_dpitchmid_us = fit_rollaccel_dpitchmid_us;

WBmod_linfit.fit_rollaccel_dApitchmid = fit_rollaccel_dApitchmid;

WBmod_linfit.fit_rollaccel_ddt_ds = fit_rollaccel_ddt_ds;

WBmod_linfit.fit_rollaccel_ddt_us = fit_rollaccel_ddt_us;

WBmod_linfit.fit_rollaccel_dRds = fit_rollaccel_dRds;

WBmod_linfit.fit_rollaccel_range = [min(roll_dot_dot_mean_wb(:));max(roll_dot_dot_mean_wb(:))];


WBmod_linfit.fit_yawaccel_ddevmax_ds = fit_yawaccel_ddevmax_ds;

WBmod_linfit.fit_yawaccel_ddevmax_us = fit_yawaccel_ddevmax_us;

WBmod_linfit.fit_yawaccel_ddevmin_ds = fit_yawaccel_ddevmin_ds;

WBmod_linfit.fit_yawaccel_ddevmin_us = fit_yawaccel_ddevmin_us; 

WBmod_linfit.fit_yawaccel_dAdev_ds = fit_yawaccel_dAdev_ds;

WBmod_linfit.fit_yawaccel_dAdev_us = fit_yawaccel_dAdev_us;

WBmod_linfit.fit_yawaccel_dstrokemax = fit_yawaccel_dstrokemax;

WBmod_linfit.fit_yawaccel_dstrokemin = fit_yawaccel_dstrokemin;

WBmod_linfit.fit_yawaccel_dAstroke = fit_yawaccel_dAstroke;

WBmod_linfit.fit_yawaccel_dpitchmid_ds = fit_yawaccel_dpitchmid_ds;

WBmod_linfit.fit_yawaccel_dpitchmid_us = fit_yawaccel_dpitchmid_us;

WBmod_linfit.fit_yawaccel_dApitchmid = fit_yawaccel_dApitchmid;

WBmod_linfit.fit_yawaccel_ddt_ds = fit_yawaccel_ddt_ds;

WBmod_linfit.fit_yawaccel_ddt_us = fit_yawaccel_ddt_us;

WBmod_linfit.fit_yawaccel_dRds = fit_yawaccel_dRds;

WBmod_linfit.fit_yawaccel_range = [min(yaw_dot_dot_mean_wb(:));max(yaw_dot_dot_mean_wb(:))];


WBmod_linfit.fit_pitchaccel_devmax_ds = fit_pitchaccel_devmax_ds;

WBmod_linfit.fit_pitchaccel_devmax_us = fit_pitchaccel_devmax_us;

WBmod_linfit.fit_pitchaccel_devmin_ds = fit_pitchaccel_devmin_ds;

WBmod_linfit.fit_pitchaccel_devmin_us = fit_pitchaccel_devmin_us;

WBmod_linfit.fit_pitchaccel_Adev_ds = fit_pitchaccel_Adev_ds;

WBmod_linfit.fit_pitchaccel_Adev_us = fit_pitchaccel_Adev_us;

WBmod_linfit.fit_pitchaccel_strokemax = fit_pitchaccel_strokemax;

WBmod_linfit.fit_pitchaccel_strokemin = fit_pitchaccel_strokemin;

WBmod_linfit.fit_pitchaccel_Astroke = fit_pitchaccel_Astroke;

WBmod_linfit.fit_pitchaccel_pitchmid_ds = fit_pitchaccel_pitchmid_ds;

WBmod_linfit.fit_pitchaccel_pitchmid_us = fit_pitchaccel_pitchmid_us;

WBmod_linfit.fit_pitchaccel_Apitchmid = fit_pitchaccel_Apitchmid;

WBmod_linfit.fit_pitchaccel_dt_ds = fit_pitchaccel_dt_ds;

WBmod_linfit.fit_pitchaccel_dt_us = fit_pitchaccel_dt_us;

WBmod_linfit.fit_pitchaccel_freq = fit_pitchaccel_freq;

WBmod_linfit.fit_pitchaccel_range = [min(pitch_dot_dot_mean_wb(:));max(pitch_dot_dot_mean_wb(:))];


WBmod_linfit.fit_force_devmax_ds = fit_force_devmax_ds;

WBmod_linfit.fit_force_devmax_us = fit_force_devmax_us;

WBmod_linfit.fit_force_devmin_ds = fit_force_devmin_ds;

WBmod_linfit.fit_force_devmin_us = fit_force_devmin_us;

WBmod_linfit.fit_force_Adev_ds = fit_force_Adev_ds;

WBmod_linfit.fit_force_Adev_us = fit_force_Adev_us;

WBmod_linfit.fit_force_strokemax = fit_force_strokemax;

WBmod_linfit.fit_force_strokemin = fit_force_strokemin;

WBmod_linfit.fit_force_Astroke = fit_force_Astroke;

WBmod_linfit.fit_force_pitchmid_ds = fit_force_pitchmid_ds;

WBmod_linfit.fit_force_pitchmid_us = fit_force_pitchmid_us;

WBmod_linfit.fit_force_Apitchmid = fit_force_Apitchmid;

WBmod_linfit.fit_force_dt_ds = fit_force_dt_ds;

WBmod_linfit.fit_force_dt_us = fit_force_dt_us;

WBmod_linfit.fit_force_freq = fit_force_freq;

WBmod_linfit.fit_force_range = [min(F_mean_wb(:));max(F_mean_wb(:))];


save('WBmod_linear_fits.mat','WBmod_linfit')
































