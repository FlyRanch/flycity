% body kin
% roll, pitch
figure
subplot(3,3,1)
hold on
plot(t_wb_seq_pre_mean_all,roll_mean_wb_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,roll_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,roll_mean_wb_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,roll_mean_wb_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,roll_mean_wb_seq_mean_all,':k')
ciplot(roll_mean_wb_seq_mean_all-1.96*roll_mean_wb_seq_ste_all,roll_mean_wb_seq_mean_all+1.96*roll_mean_wb_seq_ste_all,t_wb_seq_mean_all)

% xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('roll angle [deg]')
ylim([-30 60])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-90:30:180])

subplot(3,3,2)
hold on
plot(t_wb_seq_pre_mean_all,roll_dot_mean_wb_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,roll_dot_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,roll_dot_mean_wb_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,roll_dot_mean_wb_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,roll_dot_mean_wb_seq_mean_all,':k')
ciplot(roll_dot_mean_wb_seq_mean_all-1.96*roll_dot_mean_wb_seq_ste_all,roll_dot_mean_wb_seq_mean_all+1.96*roll_dot_mean_wb_seq_ste_all,t_wb_seq_mean_all)

% xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('roll velocity [deg/sec]')
ylim([-1500 1500])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-1500:1500:1500])

subplot(3,3,3)
hold on
plot(t_wb_seq_pre_mean_all,roll_dot_dot_mean_wb_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,roll_dot_dot_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,roll_dot_dot_mean_wb_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,roll_dot_dot_mean_wb_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,roll_dot_dot_mean_wb_seq_mean_all,':k')
ciplot(roll_dot_dot_mean_wb_seq_mean_all-1.96*roll_dot_dot_mean_wb_seq_ste_all,roll_dot_dot_mean_wb_seq_mean_all+1.96*roll_dot_dot_mean_wb_seq_ste_all,t_wb_seq_mean_all)

% xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('roll acceleration/freq^2 [-]')
ylim([-3 2])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-3:1:2])

subplot(3,3,4)
hold on
plot(t_wb_seq_pre_mean_all,pitch_mean_wb_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,pitch_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,pitch_mean_wb_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,pitch_mean_wb_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,pitch_mean_wb_seq_mean_all,':k')
ciplot(pitch_mean_wb_seq_mean_all-1.96*pitch_mean_wb_seq_ste_all,pitch_mean_wb_seq_mean_all+1.96*pitch_mean_wb_seq_ste_all,t_wb_seq_mean_all)

% xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('pitch angle [deg]')
ylim([-30 60])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-90:30:180])

subplot(3,3,5)
hold on
plot(t_wb_seq_pre_mean_all,pitch_dot_mean_wb_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,pitch_dot_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,pitch_dot_mean_wb_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,pitch_dot_mean_wb_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,pitch_dot_mean_wb_seq_mean_all,':k')
ciplot(pitch_dot_mean_wb_seq_mean_all-1.96*pitch_dot_mean_wb_seq_ste_all,pitch_dot_mean_wb_seq_mean_all+1.96*pitch_dot_mean_wb_seq_ste_all,t_wb_seq_mean_all)

% xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('pitch velocity [deg/sec]')
ylim([-1500 1500])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-1500:1500:1500])

subplot(3,3,6)
hold on
plot(t_wb_seq_pre_mean_all,pitch_dot_dot_mean_wb_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,pitch_dot_dot_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,pitch_dot_dot_mean_wb_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,pitch_dot_dot_mean_wb_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,pitch_dot_dot_mean_wb_seq_mean_all,':k')
ciplot(pitch_dot_dot_mean_wb_seq_mean_all-1.96*pitch_dot_dot_mean_wb_seq_ste_all,pitch_dot_dot_mean_wb_seq_mean_all+1.96*pitch_dot_dot_mean_wb_seq_ste_all,t_wb_seq_mean_all)

% xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('pitch acceleration/freq^2 [-]')
ylim([-3 2])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-3:1:2])

subplot(3,3,7)
hold on
plot(t_wb_seq_pre_mean_all,yaw_mean_wb_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,yaw_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,yaw_mean_wb_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,yaw_mean_wb_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,yaw_mean_wb_seq_mean_all,':k')
ciplot(yaw_mean_wb_seq_mean_all-1.96*yaw_mean_wb_seq_ste_all,yaw_mean_wb_seq_mean_all+1.96*yaw_mean_wb_seq_ste_all,t_wb_seq_mean_all)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('yaw angle [deg]')
ylim([0 90])
set(gca,'XTick',-0.05:0.05:.05) 
set(gca,'YTick',[-90:30:180])

subplot(3,3,8)
hold on
plot(t_wb_seq_pre_mean_all,yaw_dot_mean_wb_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,yaw_dot_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,yaw_dot_mean_wb_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,yaw_dot_mean_wb_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,yaw_dot_mean_wb_seq_mean_all,':k')
ciplot(yaw_dot_mean_wb_seq_mean_all-1.96*yaw_dot_mean_wb_seq_ste_all,yaw_dot_mean_wb_seq_mean_all+1.96*yaw_dot_mean_wb_seq_ste_all,t_wb_seq_mean_all)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('yaw velocity [deg/sec]')
ylim([-1000 2000])
set(gca,'XTick',-0.05:0.05:.05) 
set(gca,'YTick',[-1000:1000:2000])

subplot(3,3,9)
hold on
plot(t_wb_seq_pre_mean_all,yaw_dot_dot_mean_wb_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,yaw_dot_dot_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,yaw_dot_dot_mean_wb_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,yaw_dot_dot_mean_wb_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,yaw_dot_dot_mean_wb_seq_mean_all,':k')
ciplot(yaw_dot_dot_mean_wb_seq_mean_all-1.96*yaw_dot_dot_mean_wb_seq_ste_all,yaw_dot_dot_mean_wb_seq_mean_all+1.96*yaw_dot_dot_mean_wb_seq_ste_all,t_wb_seq_mean_all)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('yaw acceleration/freq^2 [-]')
ylim([-2 3])
set(gca,'XTick',-0.05:0.05:.05) 
set(gca,'YTick',[-3:1:3])

saveas(gca,['BodyDynvsTime_RollPitchYaw_allNmean.fig'])
saveas(gca,['BodyDynvsTime_RollPitchYaw_allNmean.png'])
plot2svg(['BodyDynvsTime_RollPitchYaw_allNmean.svg'])

%% rotL&R
figure
subplot(3,3,1)
hold on
plot(t_wb_seq_pre_mean_all,drot_R_mean_wb_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,drot_R_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,drot_R_mean_wb_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,drot_R_mean_wb_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,drot_R_mean_wb_seq_mean_all,':k')
ciplot(drot_R_mean_wb_seq_mean_all-1.96*drot_R_mean_wb_seq_ste_all,drot_R_mean_wb_seq_mean_all+1.96*drot_R_mean_wb_seq_ste_all,t_wb_seq_mean_all)

% xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('right axis rotation angle [deg]')
ylim([-30 60])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-90:30:180])

subplot(3,3,2)
hold on
plot(t_wb_seq_pre_mean_all,rot_dot_R_mean_wb_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,rot_dot_R_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,rot_dot_R_mean_wb_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,rot_dot_R_mean_wb_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,rot_dot_R_mean_wb_seq_mean_all,':k')
ciplot(rot_dot_R_mean_wb_seq_mean_all-1.96*rot_dot_R_mean_wb_seq_ste_all,rot_dot_R_mean_wb_seq_mean_all+1.96*rot_dot_R_mean_wb_seq_ste_all,t_wb_seq_mean_all)

% xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('right axis rotation velocity [deg/sec]')
ylim([-1000 2000])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-2000:1000:2000])

subplot(3,3,3)
hold on
plot(t_wb_seq_pre_mean_all,rot_dot_dot_R_mean_wb_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,rot_dot_dot_R_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,rot_dot_dot_R_mean_wb_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,rot_dot_dot_R_mean_wb_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,rot_dot_dot_R_mean_wb_seq_mean_all,':k')
ciplot(rot_dot_dot_R_mean_wb_seq_mean_all-1.96*rot_dot_dot_R_mean_wb_seq_ste_all,rot_dot_dot_R_mean_wb_seq_mean_all+1.96*rot_dot_dot_R_mean_wb_seq_ste_all,t_wb_seq_mean_all)

% xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('right axis rotation acceleration/freq^2 [-]')
ylim([-3 2])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-3:1:3])

subplot(3,3,4)
hold on
plot(t_wb_seq_pre_mean_all,drot_L_mean_wb_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,drot_L_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,drot_L_mean_wb_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,drot_L_mean_wb_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,drot_L_mean_wb_seq_mean_all,':k')
ciplot(drot_L_mean_wb_seq_mean_all-1.96*drot_L_mean_wb_seq_ste_all,drot_L_mean_wb_seq_mean_all+1.96*drot_L_mean_wb_seq_ste_all,t_wb_seq_mean_all)

% xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('left axis rotation angle [deg]')
ylim([-60 30])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-90:30:180])

subplot(3,3,5)
hold on
plot(t_wb_seq_pre_mean_all,rot_dot_L_mean_wb_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,rot_dot_L_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,rot_dot_L_mean_wb_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,rot_dot_L_mean_wb_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,rot_dot_L_mean_wb_seq_mean_all,':k')
ciplot(rot_dot_L_mean_wb_seq_mean_all-1.96*rot_dot_L_mean_wb_seq_ste_all,rot_dot_L_mean_wb_seq_mean_all+1.96*rot_dot_L_mean_wb_seq_ste_all,t_wb_seq_mean_all)

% xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('left axis rotation velocity [deg/sec]')
ylim([-2000 1000])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-2000:1000:2000])

subplot(3,3,6)
hold on
plot(t_wb_seq_pre_mean_all,rot_dot_dot_L_mean_wb_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,rot_dot_dot_L_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,rot_dot_dot_L_mean_wb_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,rot_dot_dot_L_mean_wb_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,rot_dot_dot_L_mean_wb_seq_mean_all,':k')
ciplot(rot_dot_dot_L_mean_wb_seq_mean_all-1.96*rot_dot_dot_L_mean_wb_seq_ste_all,rot_dot_dot_L_mean_wb_seq_mean_all+1.96*rot_dot_dot_L_mean_wb_seq_ste_all,t_wb_seq_mean_all)

% xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('left axis rotation acceleration/freq^2 [-]')
ylim([-3 2])
set(gca,'XTick',-0.05:0.05:.05,'XTickLabel',[]) 
set(gca,'YTick',[-3:1:3])

subplot(3,3,7)
hold on
plot(t_wb_seq_pre_mean_all,yaw_mean_wb_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,yaw_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,yaw_mean_wb_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,yaw_mean_wb_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,yaw_mean_wb_seq_mean_all,':k')
ciplot(yaw_mean_wb_seq_mean_all-1.96*yaw_mean_wb_seq_ste_all,yaw_mean_wb_seq_mean_all+1.96*yaw_mean_wb_seq_ste_all,t_wb_seq_mean_all)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('yaw angle [deg]')
ylim([0 90])
set(gca,'XTick',-0.05:0.05:.05) 
set(gca,'YTick',[-90:30:180])

subplot(3,3,8)
hold on
plot(t_wb_seq_pre_mean_all,yaw_dot_mean_wb_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,yaw_dot_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,yaw_dot_mean_wb_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,yaw_dot_mean_wb_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,yaw_dot_mean_wb_seq_mean_all,':k')
ciplot(yaw_dot_mean_wb_seq_mean_all-1.96*yaw_dot_mean_wb_seq_ste_all,yaw_dot_mean_wb_seq_mean_all+1.96*yaw_dot_mean_wb_seq_ste_all,t_wb_seq_mean_all)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('yaw velocity [deg/sec]')
ylim([-1000 2000])
set(gca,'XTick',-0.05:0.05:.05) 
set(gca,'YTick',[-1000:1000:2000])

subplot(3,3,9)
hold on
plot(t_wb_seq_pre_mean_all,yaw_dot_dot_mean_wb_seq_pre_all,'.','color',[0 .5 1])
plot(t_wb_seq_post_mean_all,yaw_dot_dot_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean_all,yaw_dot_dot_mean_wb_seq_pre_mean_all,'color',[0 .5 1])
plot(t_wb_seq_post_mean_all,yaw_dot_dot_mean_wb_seq_post_mean_all,'color',[1 .5 0])

plot(t_wb_seq_mean_all,yaw_dot_dot_mean_wb_seq_mean_all,':k')
ciplot(yaw_dot_dot_mean_wb_seq_mean_all-1.96*yaw_dot_dot_mean_wb_seq_ste_all,yaw_dot_dot_mean_wb_seq_mean_all+1.96*yaw_dot_dot_mean_wb_seq_ste_all,t_wb_seq_mean_all)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('yaw acceleration/freq^2 [-]')
ylim([-2 3])
set(gca,'XTick',-0.05:0.05:.05) 
set(gca,'YTick',[-3:1:3])

saveas(gca,['BodyDynvsTime_rotLnRnYaw_allNmean.fig'])
saveas(gca,['BodyDynvsTime_rotLnRnYaw_allNmean.png'])
plot2svg(['BodyDynvsTime_rotLnRnYaw_allNmean.svg'])
