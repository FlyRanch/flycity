

% move wrap point to 90 deg
% turn_angle_vel_mirror(turn_angle_vel_mirror>90) = turn_angle_vel_mirror(turn_angle_vel_mirror>90)-360;
xp = [min(turn_angle_vel_mirror),max(turn_angle_vel_mirror)];

droll_mirror_Amax = calc_value(droll_mirror,n_Amax);
yaw_mirror_Amax = calc_value(yaw_mirror,n_Amax);
Dpitch_Amax = calc_value(dpitch,n_Amax);
pitch_Amax = calc_value(pitch,n_Amax);
rot_angle_Amax = atand(Dpitch_Amax./droll_mirror_Amax);

droll_mirror_mean = calc_circ_mean_value(droll_mirror,n_pre,n_post);
yaw_mirror_mean = calc_circ_mean_value(yaw_mirror,n_pre,n_post);
Dpitch_mean = calc_circ_mean_value(dpitch,n_pre,n_post);
pitch_mean = calc_circ_mean_value(pitch,n_pre,n_post);
rot_angle_mean = atand(Dpitch_mean./droll_mirror_mean);

% stats Amax
stats_droll_atAmax = regstats(droll_mirror_Amax,turn_angle_vel_mirror,'linear');
stats_Dpitch_atAmax = regstats(Dpitch_Amax,turn_angle_vel_mirror,'linear');
stats_yaw_atAmax = regstats(yaw_mirror_Amax,turn_angle_vel_mirror,'linear');
stats_rot_angle_atAmax = regstats(rot_angle_Amax,turn_angle_vel_mirror,'linear');

stats_droll_atAmax.beta
stats_droll_atAmax.tstat.pval
[h,p]=ttest(droll_mirror_Amax)

stats_Dpitch_atAmax.beta
stats_Dpitch_atAmax.tstat.pval
[h,p]=ttest(Dpitch_Amax)

stats_yaw_atAmax.beta
stats_yaw_atAmax.tstat.pval
[h,p]=ttest(yaw_mirror_Amax)

stats_rot_angle_atAmax.beta
stats_rot_angle_atAmax.tstat.pval
[h,p]=ttest(rot_angle_Amax)

% stats means
stats_droll_mean = regstats(droll_mirror_mean,turn_angle_vel_mirror,'linear');
stats_Dpitch_mean = regstats(Dpitch_mean,turn_angle_vel_mirror,'linear');
stats_yaw_mean = regstats(yaw_mirror_mean,turn_angle_vel_mirror,'linear');
stats_rot_angle_mean = regstats(rot_angle_mean,turn_angle_vel_mirror,'linear');

stats_droll_mean.beta
stats_droll_mean.tstat.pval
[h,p]=ttest(droll_mirror_mean)

stats_Dpitch_mean.beta
stats_Dpitch_mean.tstat.pval
[h,p]=ttest(Dpitch_mean)

stats_yaw_mean.beta
stats_yaw_mean.tstat.pval
[h,p]=ttest(yaw_mirror_mean)

stats_rot_angle_mean.beta
stats_rot_angle_mean.tstat.pval
[h,p]=ttest(rot_angle_mean)

% droll
figure
subplot(2,2,1)
hold on

yp = polyval([stats_droll_atAmax.beta(2) stats_droll_atAmax.beta(1)],xp);
plot(xp,yp)
plot(turn_angle_vel_mirror,droll_mirror_Amax,'.')

yp = polyval([stats_droll_mean.beta(2) stats_droll_mean.beta(1)],xp);
plot(xp,yp,'r')
plot(turn_angle_vel_mirror,droll_mirror_mean,'.r')

axis([0 180 0 45])
xlabel('turn angle')
ylabel('droll')

% pitch
subplot(2,2,2)
hold on

yp = polyval([stats_Dpitch_atAmax.beta(2) stats_Dpitch_atAmax.beta(1)],xp);
plot(xp,yp)
plot(turn_angle_vel_mirror,Dpitch_Amax,'.')

yp = polyval([stats_Dpitch_mean.beta(2) stats_Dpitch_mean.beta(1)],xp);
plot(xp,yp,'r')
plot(turn_angle_vel_mirror,Dpitch_mean,'.r')

axis([0 180 0 45])
xlabel('turn angle')
ylabel('dpitch')

% yaw
subplot(2,2,3)
hold on

yp = polyval([stats_yaw_atAmax.beta(2) stats_yaw_atAmax.beta(1)],xp);
plot(xp,yp)
plot(turn_angle_vel_mirror,yaw_mirror_Amax,'.')

yp = polyval([stats_yaw_mean.beta(2) stats_yaw_mean.beta(1)],xp);
plot(xp,yp,'r')
plot(turn_angle_vel_mirror,yaw_mirror_mean,'.r')

axis([0 180 0 90])
xlabel('turn angle')
ylabel('yaw')

% rot_angle
subplot(2,2,4)
hold on

yp = polyval([stats_rot_angle_atAmax.beta(2) stats_rot_angle_atAmax.beta(1)],xp);
plot(xp,yp)
plot(turn_angle_vel_mirror,rot_angle_Amax,'.')

yp = polyval([stats_rot_angle_mean.beta(2) stats_rot_angle_mean.beta(1)],xp);
plot(xp,yp,'r')
plot(turn_angle_vel_mirror,rot_angle_mean,'.r')

axis([0 180 0 90])
xlabel('turn angle')
ylabel('rotation axis angle')

% save data
save('attitude_atAmaxNmean_4stats.mat','turn_angle_vel_mirror',...
    'rot_angle_Amax','droll_mirror_Amax','yaw_mirror_Amax','Dpitch_Amax','pitch_Amax',...
    'rot_angle_mean','droll_mirror_mean','yaw_mirror_mean','Dpitch_mean','pitch_mean',...
    'stats_rot_angle_atAmax','stats_droll_atAmax','stats_Dpitch_atAmax','stats_yaw_atAmax',...
    'stats_rot_angle_mean','stats_droll_mean','stats_Dpitch_mean','stats_yaw_mean');



