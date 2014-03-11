

%% torque
% Mroll, Mpitch
figure
subplot(3,3,1)
hold on
plot(t_wb_seq_pre_mean,Mroll_mean_wb_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,Mroll_mean_wb_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,Mroll_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,Mroll_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,Mroll_mean_wb_seq_mean,':k')

% xlabel('time [sec]')
ylabel('Mroll')
title('Total Torque')
ylim([-.05 .05])
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,2)
hold on
plot(t_wb_seq_pre_mean,Mroll_mean_wb_accel_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,Mroll_mean_wb_accel_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,Mroll_mean_wb_accel_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,Mroll_mean_wb_accel_seq_post_mean,'c')
plot(t_wb_seq_mean,Mroll_mean_wb_accel_seq_mean,':k')

% xlabel('time [sec]')
% ylabel('Mroll')
title('Inertial Torque')
ylim([-.05 .05])
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,3)
hold on
plot(t_wb_seq_pre_mean,Mroll_mean_wb_damp_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,Mroll_mean_wb_damp_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,Mroll_mean_wb_damp_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,Mroll_mean_wb_damp_seq_post_mean,'c')
plot(t_wb_seq_mean,Mroll_mean_wb_damp_seq_mean,':k')

% xlabel('time [sec]')
% ylabel('Mroll')
title('Damping Torque')
ylim([-.05 .05])
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,4)
hold on
plot(t_wb_seq_pre_mean,Mpitch_mean_wb_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,Mpitch_mean_wb_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,Mpitch_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,Mpitch_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,Mpitch_mean_wb_seq_mean,':k')

% xlabel('time [sec]')
ylabel('Mpitch')
ylim([-.05 .05])
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,5)
hold on
plot(t_wb_seq_pre_mean,Mpitch_mean_wb_accel_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,Mpitch_mean_wb_accel_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,Mpitch_mean_wb_accel_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,Mpitch_mean_wb_accel_seq_post_mean,'c')
plot(t_wb_seq_mean,Mpitch_mean_wb_accel_seq_mean,':k')

% xlabel('time [sec]')
% ylabel('Mpitch')
ylim([-.05 .05])
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,6)
hold on
plot(t_wb_seq_pre_mean,Mpitch_mean_wb_damp_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,Mpitch_mean_wb_damp_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,Mpitch_mean_wb_damp_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,Mpitch_mean_wb_damp_seq_post_mean,'c')
plot(t_wb_seq_mean,Mpitch_mean_wb_damp_seq_mean,':k')

% xlabel('time [sec]')
% ylabel('Mpitch')
ylim([-.05 .05])
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,7)
hold on
plot(t_wb_seq_pre_mean,Myaw_mean_wb_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,Myaw_mean_wb_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,Myaw_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,Myaw_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,Myaw_mean_wb_seq_mean,':k')

xlabel('time [sec]')
ylabel('Myaw')
ylim([-.05 .05])
set(gca,'YTick',[-.1:.025:1])

subplot(3,3,8)
hold on
plot(t_wb_seq_pre_mean,Myaw_mean_wb_accel_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,Myaw_mean_wb_accel_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,Myaw_mean_wb_accel_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,Myaw_mean_wb_accel_seq_post_mean,'c')
plot(t_wb_seq_mean,Myaw_mean_wb_accel_seq_mean,':k')

xlabel('time [sec]')
% ylabel('Myaw')
ylim([-.05 .05])
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,9)
hold on
plot(t_wb_seq_pre_mean,Myaw_mean_wb_damp_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,Myaw_mean_wb_damp_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,Myaw_mean_wb_damp_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,Myaw_mean_wb_damp_seq_post_mean,'c')
plot(t_wb_seq_mean,Myaw_mean_wb_damp_seq_mean,':k')

xlabel('time [sec]')
% ylabel('Myaw')
ylim([-.05 .05])
set(gca,'YTick',[-.1:.025:1])

%% M_R, M_L
figure
subplot(3,3,1)
hold on
plot(t_wb_seq_pre_mean,M_R_mean_wb_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,M_R_mean_wb_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,M_R_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,M_R_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,M_R_mean_wb_seq_mean,':k')

% xlabel('time [sec]')
ylabel('M_R')
title('Total Torque')
ylim([-.05 .05])
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,2)
hold on
plot(t_wb_seq_pre_mean,M_R_mean_wb_accel_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,M_R_mean_wb_accel_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,M_R_mean_wb_accel_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,M_R_mean_wb_accel_seq_post_mean,'c')
plot(t_wb_seq_mean,M_R_mean_wb_accel_seq_mean,':k')

% xlabel('time [sec]')
% ylabel('M_R')
title('Inertial Torque')
ylim([-.05 .05])
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,3)
hold on
plot(t_wb_seq_pre_mean,M_R_mean_wb_damp_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,M_R_mean_wb_damp_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,M_R_mean_wb_damp_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,M_R_mean_wb_damp_seq_post_mean,'c')
plot(t_wb_seq_mean,M_R_mean_wb_damp_seq_mean,':k')

% xlabel('time [sec]')
% ylabel('M_R')
title('Damping Torque')
ylim([-.05 .05])
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,4)
hold on
plot(t_wb_seq_pre_mean,M_L_mean_wb_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,M_L_mean_wb_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,M_L_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,M_L_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,M_L_mean_wb_seq_mean,':k')

% xlabel('time [sec]')
ylabel('M_L')
ylim([-.05 .05])
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,5)
hold on
plot(t_wb_seq_pre_mean,M_L_mean_wb_accel_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,M_L_mean_wb_accel_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,M_L_mean_wb_accel_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,M_L_mean_wb_accel_seq_post_mean,'c')
plot(t_wb_seq_mean,M_L_mean_wb_accel_seq_mean,':k')

% xlabel('time [sec]')
% ylabel('M_L')
ylim([-.05 .05])
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,6)
hold on
plot(t_wb_seq_pre_mean,M_L_mean_wb_damp_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,M_L_mean_wb_damp_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,M_L_mean_wb_damp_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,M_L_mean_wb_damp_seq_post_mean,'c')
plot(t_wb_seq_mean,M_L_mean_wb_damp_seq_mean,':k')

% xlabel('time [sec]')
% ylabel('M_L')
ylim([-.05 .05])
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,7)
hold on
plot(t_wb_seq_pre_mean,Myaw_mean_wb_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,Myaw_mean_wb_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,Myaw_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,Myaw_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,Myaw_mean_wb_seq_mean,':k')

xlabel('time [sec]')
ylabel('Myaw')
ylim([-.05 .05])
set(gca,'YTick',[-.1:.025:1])

subplot(3,3,8)
hold on
plot(t_wb_seq_pre_mean,Myaw_mean_wb_accel_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,Myaw_mean_wb_accel_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,Myaw_mean_wb_accel_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,Myaw_mean_wb_accel_seq_post_mean,'c')
plot(t_wb_seq_mean,Myaw_mean_wb_accel_seq_mean,':k')

xlabel('time [sec]')
% ylabel('Myaw')
ylim([-.05 .05])
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,9)
hold on
plot(t_wb_seq_pre_mean,Myaw_mean_wb_damp_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,Myaw_mean_wb_damp_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,Myaw_mean_wb_damp_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,Myaw_mean_wb_damp_seq_post_mean,'c')
plot(t_wb_seq_mean,Myaw_mean_wb_damp_seq_mean,':k')

xlabel('time [sec]')
% ylabel('Myaw')
ylim([-.05 .05])
set(gca,'YTick',[-.1:.025:1])

