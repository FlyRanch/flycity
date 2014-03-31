%% force
figure
subplot(3,1,1)
hold on
plot(t_wb_seq_pre_mean_all,F_mean_wb_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,F_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,F_mean_wb_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,F_mean_wb_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,F_mean_wb_seq_mean_all,':k')
ciplot(F_mean_wb_seq_mean_all-1.96*F_mean_wb_seq_ste_all,F_mean_wb_seq_mean_all+1.96*F_mean_wb_seq_ste_all,t_wb_seq_mean_all)

xlabel('time [sec]')
ylabel('F/mg')
xlim([-.05 .06]) 
ylim([.5 1.5])
set(gca,'XTick',-0.05:0.05:.05) 
set(gca,'YTick',[0:.5:2]) 

saveas(gca,['ForceVsTime_allNmean.fig'])
saveas(gca,['ForceVsTime_allNmean.png'])
plot2svg(['ForceVsTime_allNmean.svg'])

%% torque
% Mroll, Mpitch
figure
subplot(3,3,1)
hold on
plot(t_wb_seq_pre_mean_all,Mroll_mean_wb_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,Mroll_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,Mroll_mean_wb_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,Mroll_mean_wb_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,Mroll_mean_wb_seq_mean_all,':k')
ciplot(Mroll_mean_wb_seq_mean_all-1.96*Mroll_mean_wb_seq_ste_all,Mroll_mean_wb_seq_mean_all+1.96*Mroll_mean_wb_seq_ste_all,t_wb_seq_mean_all)

% xlabel('time [sec]')
ylabel('Mroll')
title('Total Torque')
xlim([-.05 .06]) 
ylim([-.05 .05])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,2)
hold on
plot(t_wb_seq_pre_mean_all,Mroll_mean_wb_accel_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,Mroll_mean_wb_accel_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,Mroll_mean_wb_accel_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,Mroll_mean_wb_accel_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,Mroll_mean_wb_accel_seq_mean_all,':k')
ciplot(Mroll_mean_wb_accel_seq_mean_all-1.96*Mroll_mean_wb_accel_seq_ste_all,Mroll_mean_wb_accel_seq_mean_all+1.96*Mroll_mean_wb_accel_seq_ste_all,t_wb_seq_mean_all)

% xlabel('time [sec]')
% ylabel('Mroll')
title('Inertial Torque')
xlim([-.05 .06]) 
ylim([-.05 .05])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,3)
hold on
plot(t_wb_seq_pre_mean_all,Mroll_mean_wb_damp_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,Mroll_mean_wb_damp_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,Mroll_mean_wb_damp_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,Mroll_mean_wb_damp_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,Mroll_mean_wb_damp_seq_mean_all,':k')
ciplot(Mroll_mean_wb_damp_seq_mean_all-1.96*Mroll_mean_wb_damp_seq_ste_all,Mroll_mean_wb_damp_seq_mean_all+1.96*Mroll_mean_wb_damp_seq_ste_all,t_wb_seq_mean_all)

% xlabel('time [sec]')
% ylabel('Mroll')
title('Damping Torque')
xlim([-.05 .06]) 
ylim([-.05 .05])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,4)
hold on
plot(t_wb_seq_pre_mean_all,Mpitch_mean_wb_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,Mpitch_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,Mpitch_mean_wb_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,Mpitch_mean_wb_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,Mpitch_mean_wb_seq_mean_all,':k')
ciplot(Mpitch_mean_wb_seq_mean_all-1.96*Mpitch_mean_wb_seq_ste_all,Mpitch_mean_wb_seq_mean_all+1.96*Mpitch_mean_wb_seq_ste_all,t_wb_seq_mean_all)

% xlabel('time [sec]')
ylabel('Mpitch')
xlim([-.05 .06]) 
ylim([-.05 .05])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,5)
hold on
plot(t_wb_seq_pre_mean_all,Mpitch_mean_wb_accel_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,Mpitch_mean_wb_accel_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,Mpitch_mean_wb_accel_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,Mpitch_mean_wb_accel_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,Mpitch_mean_wb_accel_seq_mean_all,':k')
ciplot(Mpitch_mean_wb_accel_seq_mean_all-1.96*Mpitch_mean_wb_accel_seq_ste_all,Mpitch_mean_wb_accel_seq_mean_all+1.96*Mpitch_mean_wb_accel_seq_ste_all,t_wb_seq_mean_all)

% xlabel('time [sec]')
% ylabel('Mpitch')
xlim([-.05 .06]) 
ylim([-.05 .05])
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,6)
hold on
plot(t_wb_seq_pre_mean_all,Mpitch_mean_wb_damp_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,Mpitch_mean_wb_damp_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,Mpitch_mean_wb_damp_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,Mpitch_mean_wb_damp_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,Mpitch_mean_wb_damp_seq_mean_all,':k')
ciplot(Mpitch_mean_wb_damp_seq_mean_all-1.96*Mpitch_mean_wb_damp_seq_ste_all,Mpitch_mean_wb_damp_seq_mean_all+1.96*Mpitch_mean_wb_damp_seq_ste_all,t_wb_seq_mean_all)

% xlabel('time [sec]')
% ylabel('Mpitch')
xlim([-.05 .06]) 
ylim([-.05 .05])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,7)
hold on
plot(t_wb_seq_pre_mean_all,Myaw_mean_wb_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,Myaw_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,Myaw_mean_wb_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,Myaw_mean_wb_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,Myaw_mean_wb_seq_mean_all,':k')
ciplot(Myaw_mean_wb_seq_mean_all-1.96*Myaw_mean_wb_seq_ste_all,Myaw_mean_wb_seq_mean_all+1.96*Myaw_mean_wb_seq_ste_all,t_wb_seq_mean_all)

xlabel('time [sec]')
ylabel('Myaw')
xlim([-.05 .06]) 
ylim([-.0 .1])
set(gca,'XTick',-0.05:0.05:.05)
set(gca,'YTick',[-.1:.025:1])

subplot(3,3,8)
hold on
plot(t_wb_seq_pre_mean_all,Myaw_mean_wb_accel_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,Myaw_mean_wb_accel_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,Myaw_mean_wb_accel_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,Myaw_mean_wb_accel_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,Myaw_mean_wb_accel_seq_mean_all,':k')
ciplot(Myaw_mean_wb_accel_seq_mean_all-1.96*Myaw_mean_wb_accel_seq_ste_all,Myaw_mean_wb_accel_seq_mean_all+1.96*Myaw_mean_wb_accel_seq_ste_all,t_wb_seq_mean_all)

xlabel('time [sec]')
% ylabel('Myaw')
xlim([-.05 .06]) 
ylim([-.05 .05])
set(gca,'XTick',-0.05:0.05:.05)
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,9)
hold on
plot(t_wb_seq_pre_mean_all,Myaw_mean_wb_damp_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,Myaw_mean_wb_damp_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,Myaw_mean_wb_damp_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,Myaw_mean_wb_damp_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,Myaw_mean_wb_damp_seq_mean_all,':k')
ciplot(Myaw_mean_wb_damp_seq_mean_all-1.96*Myaw_mean_wb_damp_seq_ste_all,Myaw_mean_wb_damp_seq_mean_all+1.96*Myaw_mean_wb_damp_seq_ste_all,t_wb_seq_mean_all)

xlabel('time [sec]')
% ylabel('Myaw')
xlim([-.05 .06]) 
ylim([-.0 .1])
set(gca,'XTick',-0.05:0.05:.05)
set(gca,'YTick',[-.1:.025:1])

saveas(gca,['TorqueVsTime_RollPitchYaw_allNmean.fig'])
saveas(gca,['TorqueVsTime_RollPitchYaw_allNmean.png'])
plot2svg(['TorqueVsTime_RollPitchYaw_allNmean.svg'])

%% M_R, M_L
figure
subplot(3,3,1)
hold on
plot(t_wb_seq_pre_mean_all,M_R_mean_wb_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,M_R_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,M_R_mean_wb_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,M_R_mean_wb_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,M_R_mean_wb_seq_mean_all,':k')
ciplot(M_R_mean_wb_seq_mean_all-1.96*M_R_mean_wb_seq_ste_all,M_R_mean_wb_seq_mean_all+1.96*M_R_mean_wb_seq_ste_all,t_wb_seq_mean_all)

% xlabel('time [sec]')
ylabel('M_R')
title('Total Torque')
xlim([-.05 .06]) 
ylim([-.05 .05])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,2)
hold on
plot(t_wb_seq_pre_mean_all,M_R_mean_wb_accel_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,M_R_mean_wb_accel_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,M_R_mean_wb_accel_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,M_R_mean_wb_accel_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,M_R_mean_wb_accel_seq_mean_all,':k')
ciplot(M_R_mean_wb_accel_seq_mean_all-1.96*M_R_mean_wb_accel_seq_ste_all,M_R_mean_wb_accel_seq_mean_all+1.96*M_R_mean_wb_accel_seq_ste_all,t_wb_seq_mean_all)

% xlabel('time [sec]')
% ylabel('M_','color',[0 .5 1])
title('Inertial Torque')
xlim([-.05 .06]) 
ylim([-.05 .05])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,3)
hold on
plot(t_wb_seq_pre_mean_all,M_R_mean_wb_damp_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,M_R_mean_wb_damp_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,M_R_mean_wb_damp_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,M_R_mean_wb_damp_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,M_R_mean_wb_damp_seq_mean_all,':k')
ciplot(M_R_mean_wb_damp_seq_mean_all-1.96*M_R_mean_wb_damp_seq_ste_all,M_R_mean_wb_damp_seq_mean_all+1.96*M_R_mean_wb_damp_seq_ste_all,t_wb_seq_mean_all)

% xlabel('time [sec]')
% ylabel('M_','color',[0 .5 1])
title('Damping Torque')
xlim([-.05 .06]) 
ylim([-.05 .05])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,4)
hold on
plot(t_wb_seq_pre_mean_all,M_L_mean_wb_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,M_L_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,M_L_mean_wb_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,M_L_mean_wb_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,M_L_mean_wb_seq_mean_all,':k')
ciplot(M_L_mean_wb_seq_mean_all-1.96*M_L_mean_wb_seq_ste_all,M_L_mean_wb_seq_mean_all+1.96*M_L_mean_wb_seq_ste_all,t_wb_seq_mean_all)

% xlabel('time [sec]')
ylabel('M_L')
xlim([-.05 .06]) 
ylim([-.05 .05])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,5)
hold on
plot(t_wb_seq_pre_mean_all,M_L_mean_wb_accel_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,M_L_mean_wb_accel_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,M_L_mean_wb_accel_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,M_L_mean_wb_accel_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,M_L_mean_wb_accel_seq_mean_all,':k')
ciplot(M_L_mean_wb_accel_seq_mean_all-1.96*M_L_mean_wb_accel_seq_ste_all,M_L_mean_wb_accel_seq_mean_all+1.96*M_L_mean_wb_accel_seq_ste_all,t_wb_seq_mean_all)

% xlabel('time [sec]')
% ylabel('M_L')
xlim([-.05 .06]) 
ylim([-.05 .05])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,6)
hold on
plot(t_wb_seq_pre_mean_all,M_L_mean_wb_damp_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,M_L_mean_wb_damp_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,M_L_mean_wb_damp_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,M_L_mean_wb_damp_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,M_L_mean_wb_damp_seq_mean_all,':k')
ciplot(M_L_mean_wb_damp_seq_mean_all-1.96*M_L_mean_wb_damp_seq_ste_all,M_L_mean_wb_damp_seq_mean_all+1.96*M_L_mean_wb_damp_seq_ste_all,t_wb_seq_mean_all)

% xlabel('time [sec]')
% ylabel('M_L')
xlim([-.05 .06]) 
ylim([-.05 .05])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,7)
hold on
plot(t_wb_seq_pre_mean_all,Myaw_mean_wb_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,Myaw_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,Myaw_mean_wb_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,Myaw_mean_wb_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,Myaw_mean_wb_seq_mean_all,':k')
ciplot(Myaw_mean_wb_seq_mean_all-1.96*Myaw_mean_wb_seq_ste_all,Myaw_mean_wb_seq_mean_all+1.96*Myaw_mean_wb_seq_ste_all,t_wb_seq_mean_all)

xlabel('time [sec]')
ylabel('Myaw')
xlim([-.05 .06]) 
ylim([-.0 .1])
set(gca,'XTick',-0.05:0.05:.05)
set(gca,'YTick',[-.1:.025:1])

subplot(3,3,8)
hold on
plot(t_wb_seq_pre_mean_all,Myaw_mean_wb_accel_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,Myaw_mean_wb_accel_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,Myaw_mean_wb_accel_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,Myaw_mean_wb_accel_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,Myaw_mean_wb_accel_seq_mean_all,':k')
ciplot(Myaw_mean_wb_accel_seq_mean_all-1.96*Myaw_mean_wb_accel_seq_ste_all,Myaw_mean_wb_accel_seq_mean_all+1.96*Myaw_mean_wb_accel_seq_ste_all,t_wb_seq_mean_all)

xlabel('time [sec]')
% ylabel('Myaw')
xlim([-.05 .06]) 
ylim([-.05 .05])
set(gca,'XTick',-0.05:0.05:.05)
set(gca,'YTick',[-.1:.025:1]) 

subplot(3,3,9)
hold on
plot(t_wb_seq_pre_mean_all,Myaw_mean_wb_damp_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,Myaw_mean_wb_damp_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,Myaw_mean_wb_damp_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,Myaw_mean_wb_damp_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,Myaw_mean_wb_damp_seq_mean_all,':k')
ciplot(Myaw_mean_wb_damp_seq_mean_all-1.96*Myaw_mean_wb_damp_seq_ste_all,Myaw_mean_wb_damp_seq_mean_all+1.96*Myaw_mean_wb_damp_seq_ste_all,t_wb_seq_mean_all)

xlabel('time [sec]')
% ylabel('Myaw')
xlim([-.05 .06]) 
ylim([-.0 .1])
set(gca,'XTick',-0.05:0.05:.05)
set(gca,'YTick',[-.1:.025:1])

saveas(gca,['TorqueVsTime_rotLnRnYaw_allNmean.fig'])
saveas(gca,['TorqueVsTime_rotLnRnYaw_allNmean.png'])
plot2svg(['TorqueVsTime_rotLnRnYaw_allNmean.svg'])
