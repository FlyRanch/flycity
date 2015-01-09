xp = [min(turn_angle_vel_mirror),max(turn_angle_vel_mirror)];
xp = [0,180];

% max&min
stats_M_R_norm_max_startstop = regstats(M_R_norm_max_startstop,turn_angle_vel_mirror,'linear');
stats_M_R_norm_min_0stop = regstats(M_R_norm_min_0stop,turn_angle_vel_mirror,'linear');
stats_M_L_norm_max_startstop = regstats(M_L_norm_max_startstop,turn_angle_vel_mirror,'linear');
stats_M_L_norm_min_0stop = regstats(M_L_norm_min_0stop,turn_angle_vel_mirror,'linear');
stats_Myaw_norm_max_startstop = regstats(Myaw_norm_max_startstop,turn_angle_vel_mirror,'linear');

% means
stats_M_R_norm_mean_startstop = regstats(M_R_norm_mean_startstop,turn_angle_vel_mirror,'linear');
stats_M_L_norm_mean_startstop = regstats(M_L_norm_mean_startstop,turn_angle_vel_mirror,'linear');
stats_Myaw_norm_mean_startstop = regstats(Myaw_norm_mean_startstop,turn_angle_vel_mirror,'linear');

%% show results
% stats_M_R_norm_max_startstop.beta
% stats_M_R_norm_max_startstop.tstat.pval
% % [h,p]=ttest(M_R_norm_max_startstop)
% 
% stats_M_L_norm_max_startstop.beta
% stats_M_L_norm_max_startstop.tstat.pval
% % [h,p]=ttest(M_L_norm_max_startstop)
% 
% stats_Myaw_norm_max_startstop.beta
% stats_Myaw_norm_max_startstop.tstat.pval
% % [h,p]=ttest(Myaw_norm_max_startstop)
% 
% stats_M_R_norm_min_0stop.beta
% stats_M_R_norm_min_0stop.tstat.pval
% % [h,p]=ttest(M_R_norm_min_0stop)
% 
% stats_M_L_norm_min_0stop.beta
% stats_M_L_norm_min_0stop.tstat.pval
% % [h,p]=ttest(M_L_norm_min_0stop)
% 
% % stats_Myaw_norm_min_0stop.beta
% % stats_Myaw_norm_min_0stop.tstat.pval
% % % [h,p]=ttest(Myaw_norm_min_0stop)
% 
% stats_M_R_norm_mean_startstop.beta
% stats_M_R_norm_mean_startstop.tstat.pval
% % [h,p]=ttest(M_R_norm_mean_startstop)
% 
% stats_M_L_norm_mean_startstop.beta
% stats_M_L_norm_mean_startstop.tstat.pval
% % [h,p]=ttest(M_L_norm_mean_startstop)
% 
% stats_Myaw_norm_mean_startstop.beta
% stats_Myaw_norm_mean_startstop.tstat.pval
% % [h,p]=ttest(Myaw_norm_mean_startstop)

%% plot
% figure

% M_R min&max
subplot(3,2,2)
hold on

yp = polyval([stats_M_R_norm_max_startstop.beta(2) stats_M_R_norm_max_startstop.beta(1)],xp);
plot(xp,yp)
plot(turn_angle_vel_mirror,M_R_norm_max_startstop,'.')

yp = polyval([stats_M_R_norm_min_0stop.beta(2) stats_M_R_norm_min_0stop.beta(1)],xp);
plot(xp,yp,'r')
plot(turn_angle_vel_mirror,M_R_norm_min_0stop,'.r')

axis([0 180 -.045 .045])
% xlabel('turn angle')
% ylabel('M_R')

% M_L min&max
subplot(3,2,4)
hold on

yp = polyval([stats_M_L_norm_max_startstop.beta(2) stats_M_L_norm_max_startstop.beta(1)],xp);
plot(xp,yp)
plot(turn_angle_vel_mirror,M_L_norm_max_startstop,'.')

yp = polyval([stats_M_L_norm_min_0stop.beta(2) stats_M_L_norm_min_0stop.beta(1)],xp);
plot(xp,yp,'r')
plot(turn_angle_vel_mirror,M_L_norm_min_0stop,'.r')

axis([0 180 -.045 .045])
% xlabel('turn angle')
% ylabel('M_L')

% Myaw max
subplot(3,2,6)
hold on

yp = polyval([stats_Myaw_norm_max_startstop.beta(2) stats_Myaw_norm_max_startstop.beta(1)],xp);
plot(xp,yp)
plot(turn_angle_vel_mirror,Myaw_norm_max_startstop,'.')

axis([0 180 0 .09])
xlabel('turn angle')
% ylabel('Myaw')

% M_R mean
subplot(3,2,1)
hold on

yp = polyval([stats_M_R_norm_mean_startstop.beta(2) stats_M_R_norm_mean_startstop.beta(1)],xp);
plot(xp,yp,'k')
plot(turn_angle_vel_mirror,M_R_norm_mean_startstop,'.k')

axis([0 180 -.045 .045])
% xlabel('turn angle')
ylabel('M_R')

% M_L min&mean
subplot(3,2,3)
hold on

yp = polyval([stats_M_L_norm_mean_startstop.beta(2) stats_M_L_norm_mean_startstop.beta(1)],xp);
plot(xp,yp,'k')
plot(turn_angle_vel_mirror,M_L_norm_mean_startstop,'.k')

axis([0 180 -.045 .045])
% xlabel('turn angle')
ylabel('M_L')

% Myaw mean
subplot(3,2,5)
hold on

yp = polyval([stats_Myaw_norm_mean_startstop.beta(2) stats_Myaw_norm_mean_startstop.beta(1)],xp);
plot(xp,yp,'k')
plot(turn_angle_vel_mirror,Myaw_norm_mean_startstop,'.k')

axis([0 180 0 .09])
xlabel('turn angle')
ylabel('Myaw')

%% save data
% save('stats_saccadic_torque.mat','turn_angle_vel_mirror',...
%     'M_R_norm_max_startstop','M_L_norm_max_startstop','Myaw_norm_max_startstop',...
%     'M_R_norm_min_0stop','M_L_norm_min_0stop',...
%     'M_R_norm_mean_startstop','M_L_norm_mean_startstop','Myaw_norm_mean_startstop',...
%     'stats_M_R_norm_max_startstop','stats_M_L_norm_max_startstop','stats_Myaw_norm_max_startstop',...
%     'stats_M_R_norm_min_0stop','stats_M_L_norm_min_0stop',...
%     'stats_M_R_norm_mean_startstop','stats_M_L_norm_mean_startstop','stats_Myaw_norm_mean_startstop');



