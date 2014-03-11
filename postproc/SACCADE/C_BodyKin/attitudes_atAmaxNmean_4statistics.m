

% move wrap point to 90 deg
stim_angle_yaw_mirror_pre(stim_angle_yaw_mirror_pre>90) = stim_angle_yaw_mirror_pre(stim_angle_yaw_mirror_pre>90)-360;
xp = [min(stim_angle_yaw_mirror_pre),max(stim_angle_yaw_mirror_pre)];

roll_mirror_Amax = calc_value(roll_mirror,n_Amax);
yaw_mirror_Amax = calc_value(yaw_mirror,n_Amax);
Dpitch_Amax = calc_value(dpitch,n_Amax);
pitch_Amax = calc_value(pitch,n_Amax);


roll_mirror_mean = calc_circ_mean_value(roll_mirror,n_pre,n_post);
yaw_mirror_mean = calc_circ_mean_value(yaw_mirror,n_pre,n_post);
Dpitch_mean = calc_circ_mean_value(dpitch,n_pre,n_post);
pitch_mean = calc_circ_mean_value(pitch,n_pre,n_post);

% stats all
stats_roll_atAmax = regstats(roll_mirror_Amax,stim_angle_yaw_mirror_pre,'linear');
stats_Dpitch_atAmax = regstats(Dpitch_Amax,stim_angle_yaw_mirror_pre,'linear');
stats_yaw_atAmax = regstats(yaw_mirror_Amax,stim_angle_yaw_mirror_pre,'linear');

% % stats means
% stats_roll_atAmax = regstats(roll_mirror_mean,stim_angle_yaw_mirror_pre,'linear');
% stats_Dpitch_atAmax = regstats(Dpitch_mean,stim_angle_yaw_mirror_pre,'linear');
% stats_yaw_atAmax = regstats(yaw_mirror_mean,stim_angle_yaw_mirror_pre,'linear');

stats_roll_atAmax.beta
stats_roll_atAmax.tstat.pval
[h,p]=ttest(roll_mirror_Amax)
yp = polyval([stats_roll_atAmax.beta(2) stats_roll_atAmax.beta(1)],xp);

figure
subplot(1,3,1)
plot(xp,yp)
hold on
plot(stim_angle_yaw_mirror_pre,roll_mirror_Amax,'.')

stats_Dpitch_atAmax.beta
stats_Dpitch_atAmax.tstat.pval
[h,p]=ttest(Dpitch_Amax)
yp = polyval([stats_Dpitch_atAmax.beta(2) stats_Dpitch_atAmax.beta(1)],xp);

subplot(1,3,2)
plot(xp,yp)
hold on
plot(stim_angle_yaw_mirror_pre,Dpitch_Amax,'.')

stats_yaw_atAmax.beta
stats_yaw_atAmax.tstat.pval
[h,p]=ttest(yaw_mirror_Amax)
yp = polyval([stats_yaw_atAmax.beta(2) stats_yaw_atAmax.beta(1)],xp);

subplot(1,3,3)
plot(xp,yp)
hold on
plot(stim_angle_yaw_mirror_pre,yaw_mirror_Amax,'.')


save('attitude_atAmaxNmean_4stats.mat','stim_angle_yaw_mirror_pre','roll_mirror_Amax','yaw_mirror_Amax','Dpitch_Amax','pitch_Amax','roll_mirror_mean','yaw_mirror_mean','Dpitch_mean','pitch_mean','stats_roll_atAmax','stats_Dpitch_atAmax','stats_yaw_atAmax');