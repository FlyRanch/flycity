% wb kin

%% freq
figure
hold on
plot(t_wb_seq_pre_mean,f_wb_seq_pre,'.','color',[1 .5 0])
plot(t_wb_seq_post_mean,f_wb_seq_post,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean,f_wb_seq_pre_mean,'color',[1 .5 0])
plot(t_wb_seq_post_mean,f_wb_seq_post_mean,'color',[1 .5 0])

plot(t_wb_seq_mean,f_wb_seq_mean,':k')
ciplot(f_wb_seq_mean-1.96*f_wb_seq_ste,f_wb_seq_mean+1.96*f_wb_seq_ste,t_wb_seq_mean)

xlabel('time [sec]')
ylabel('freq [1/sec]')
ylim([150 230])

saveas(gca,['WBvsTimeSeqs_freq_95pCI.fig'])
saveas(gca,['WBvsTimeSeqs_freq_95pCI.png'])
plot2svg(['WBvsTimeSeqs_freq_95pCI.svg'])

%% angles ALLnMean
figure
subplot(3,1,1)
hold on
plot(t_wb_seq_bins_mean,stroke_wb_L_seq_bins_all,'b','linewidth',.5)
plot(t_wb_seq_bins_mean,stroke_wb_R_seq_bins_all,'r','linewidth',.5)

plot(t_wb_seq_bins_mean,stroke_wb_L_seq_bins_mean,'c','linewidth',2)
plot(t_wb_seq_bins_mean,stroke_wb_R_seq_bins_mean,'color',[1 .5 0],'linewidth',2)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('stroke angle [deg]')
ylim([-90 90])
set(gca,'YTick',[-90:45:90])

subplot(3,1,2)
hold on
plot(t_wb_seq_bins_mean,dev_wb_L_seq_bins_all,'b','linewidth',.5)
plot(t_wb_seq_bins_mean,dev_wb_R_seq_bins_all,'r','linewidth',.5)

plot(t_wb_seq_bins_mean,dev_wb_L_seq_bins_mean,'c','linewidth',2)
plot(t_wb_seq_bins_mean,dev_wb_R_seq_bins_mean,'color',[1 .5 0],'linewidth',2)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('deviation angle [deg]')
ylim([-30 30])
set(gca,'YTick',[-90:30:90])

subplot(3,1,3)
hold on
plot(t_wb_seq_bins_mean,pitch_wb_L_seq_bins_all,'b','linewidth',.5)
plot(t_wb_seq_bins_mean,pitch_wb_R_seq_bins_all,'r','linewidth',.5)

plot(t_wb_seq_bins_mean,pitch_wb_L_seq_bins_mean,'c','linewidth',2)
plot(t_wb_seq_bins_mean,pitch_wb_R_seq_bins_mean,'color',[1 .5 0],'linewidth',2)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('rotation angle [deg]')
ylim([0 180])
set(gca,'YTick',[-90:45:180])

saveas(gca,['WBvsTimeSeqs_WBangles_allNmean.fig'])
saveas(gca,['WBvsTimeSeqs_WBangles_allNmean.png'])
plot2svg(['WBvsTimeSeqs_WBangles_allNmean.svg'])


%% angles Mean & 95% CI
figure
subplot(3,1,1)
hold on
% plot(t_wb_seq_bins_mean,stroke_wb_L_seq_bins_all,'b','linewidth',.5)
% plot(t_wb_seq_bins_mean,stroke_wb_R_seq_bins_all,'r','linewidth',.5)

plot(t_wb_seq_bins_mean,stroke_wb_L_seq_bins_mean,'--','color',[0 .5 1],'linewidth',2)
plot(t_wb_seq_bins_mean,stroke_wb_R_seq_bins_mean,'--','color',[1 .5 0],'linewidth',2)

ciplot(stroke_wb_L_seq_bins_mean-1.96*stroke_wb_L_seq_bins_ste,...
    stroke_wb_L_seq_bins_mean+1.96*stroke_wb_L_seq_bins_ste,t_wb_seq_bins_mean,[0 .5 1])
ciplot(stroke_wb_R_seq_bins_mean-1.96*stroke_wb_R_seq_bins_ste,...
    stroke_wb_R_seq_bins_mean+1.96*stroke_wb_R_seq_bins_ste,t_wb_seq_bins_mean,[1 .5 0])

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('stroke angle [deg]')
ylim([-90 90])
set(gca,'YTick',[-90:45:90])

subplot(3,1,2)
hold on
% plot(t_wb_seq_bins_mean,dev_wb_L_seq_bins_all,'b','linewidth',.5)
% plot(t_wb_seq_bins_mean,dev_wb_R_seq_bins_all,'r','linewidth',.5)

plot(t_wb_seq_bins_mean,dev_wb_L_seq_bins_mean,'--','color',[0 .5 1],'linewidth',2)
plot(t_wb_seq_bins_mean,dev_wb_R_seq_bins_mean,'--','color',[1 .5 0],'linewidth',2)

ciplot(dev_wb_L_seq_bins_mean-1.96*dev_wb_L_seq_bins_ste,...
    dev_wb_L_seq_bins_mean+1.96*dev_wb_L_seq_bins_ste,t_wb_seq_bins_mean,[0 .5 1])
ciplot(dev_wb_R_seq_bins_mean-1.96*dev_wb_R_seq_bins_ste,...
    dev_wb_R_seq_bins_mean+1.96*dev_wb_R_seq_bins_ste,t_wb_seq_bins_mean,[1 .5 0])

xlabel('time [sec]')
xlim([-.05 .06])
ylabel('deviation angle [deg]')
ylim([-30 30])
set(gca,'YTick',[-90:30:90])

subplot(3,1,3)
hold on
% plot(t_wb_seq_bins_mean,pitch_wb_L_seq_bins_all,'b','linewidth',.5)
% plot(t_wb_seq_bins_mean,pitch_wb_R_seq_bins_all,'r','linewidth',.5)

plot(t_wb_seq_bins_mean,pitch_wb_L_seq_bins_mean,'--','color',[0 .5 1],'linewidth',2)
plot(t_wb_seq_bins_mean,pitch_wb_R_seq_bins_mean,'--','color',[1 .5 0],'linewidth',2)

ciplot(pitch_wb_L_seq_bins_mean-1.96*pitch_wb_L_seq_bins_ste,...
    pitch_wb_L_seq_bins_mean+1.96*pitch_wb_L_seq_bins_ste,t_wb_seq_bins_mean,[0 .5 1])
ciplot(pitch_wb_R_seq_bins_mean-1.96*pitch_wb_R_seq_bins_ste,...
    pitch_wb_R_seq_bins_mean+1.96*pitch_wb_R_seq_bins_ste,t_wb_seq_bins_mean,[1 .5 0])

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('rotation angle [deg]')
ylim([0 180])
set(gca,'YTick',[-90:45:180])

saveas(gca,['WBvsTimeSeqs_WBangles_meanN95pCI.fig'])
saveas(gca,['WBvsTimeSeqs_WBangles_meanN95pCI.png'])
plot2svg(['WBvsTimeSeqs_WBangles_meanN95pCI.svg'])


%% angles Only Mean
figure
subplot(3,1,1)
hold on
% plot(t_wb_seq_bins_mean,stroke_wb_L_seq_bins_all,'b','linewidth',.5)
% plot(t_wb_seq_bins_mean,stroke_wb_R_seq_bins_all,'r','linewidth',.5)

plot(t_wb_seq_bins_mean,stroke_wb_L_seq_bins_mean,'--','color',[0 .5 1],'linewidth',2)
plot(t_wb_seq_bins_mean,stroke_wb_R_seq_bins_mean,'--','color',[1 .5 0],'linewidth',2)

% ciplot(stroke_wb_L_seq_bins_mean-1.96*stroke_wb_L_seq_bins_ste,...
%     stroke_wb_L_seq_bins_mean+1.96*stroke_wb_L_seq_bins_ste,t_wb_seq_bins_mean,[0 .5 1])
% ciplot(stroke_wb_R_seq_bins_mean-1.96*stroke_wb_R_seq_bins_ste,...
%     stroke_wb_R_seq_bins_mean+1.96*stroke_wb_R_seq_bins_ste,t_wb_seq_bins_mean,[1 .5 0])

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('stroke angle [deg]')
ylim([-90 90])
set(gca,'YTick',[-90:45:90])

subplot(3,1,2)
hold on
% plot(t_wb_seq_bins_mean,dev_wb_L_seq_bins_all,'b','linewidth',.5)
% plot(t_wb_seq_bins_mean,dev_wb_R_seq_bins_all,'r','linewidth',.5)

plot(t_wb_seq_bins_mean,dev_wb_L_seq_bins_mean,'--','color',[0 .5 1],'linewidth',2)
plot(t_wb_seq_bins_mean,dev_wb_R_seq_bins_mean,'--','color',[1 .5 0],'linewidth',2)

% ciplot(dev_wb_L_seq_bins_mean-1.96*dev_wb_L_seq_bins_ste,...
%     dev_wb_L_seq_bins_mean+1.96*dev_wb_L_seq_bins_ste,t_wb_seq_bins_mean,[0 .5 1])
% ciplot(dev_wb_R_seq_bins_mean-1.96*dev_wb_R_seq_bins_ste,...
%     dev_wb_R_seq_bins_mean+1.96*dev_wb_R_seq_bins_ste,t_wb_seq_bins_mean,[1 .5 0])

xlabel('time [sec]')
xlim([-.05 .06])
ylabel('deviation angle [deg]')
ylim([-30 30])
set(gca,'YTick',[-90:30:90])

subplot(3,1,3)
hold on
% plot(t_wb_seq_bins_mean,pitch_wb_L_seq_bins_all,'b','linewidth',.5)
% plot(t_wb_seq_bins_mean,pitch_wb_R_seq_bins_all,'r','linewidth',.5)

plot(t_wb_seq_bins_mean,pitch_wb_L_seq_bins_mean,'--','color',[0 .5 1],'linewidth',2)
plot(t_wb_seq_bins_mean,pitch_wb_R_seq_bins_mean,'--','color',[1 .5 0],'linewidth',2)

% ciplot(pitch_wb_L_seq_bins_mean-1.96*pitch_wb_L_seq_bins_ste,...
%     pitch_wb_L_seq_bins_mean+1.96*pitch_wb_L_seq_bins_ste,t_wb_seq_bins_mean,[0 .5 1])
% ciplot(pitch_wb_R_seq_bins_mean-1.96*pitch_wb_R_seq_bins_ste,...
%     pitch_wb_R_seq_bins_mean+1.96*pitch_wb_R_seq_bins_ste,t_wb_seq_bins_mean,[1 .5 0])

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('rotation angle [deg]')
ylim([0 180])
set(gca,'YTick',[-90:45:180])

saveas(gca,['WBvsTimeSeqs_WBangles_mean.fig'])
saveas(gca,['WBvsTimeSeqs_WBangles_mean.png'])
plot2svg(['WBvsTimeSeqs_WBangles_mean.svg'])

%% body kin
% roll, pitch
figure
subplot(3,3,1)
hold on
plot(t_wb_seq_pre_mean,roll_mean_wb_seq_pre_all,'.c')
plot(t_wb_seq_post_mean,roll_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean,roll_mean_wb_seq_pre_mean,'c')
plot(t_wb_seq_post_mean,roll_mean_wb_seq_post_mean,'color',[1 .5 0])

plot(t_wb_seq_mean,roll_mean_wb_seq_mean,':k')
ciplot(roll_mean_wb_seq_mean-1.96*roll_mean_wb_seq_ste,roll_mean_wb_seq_mean+1.96*roll_mean_wb_seq_ste,t_wb_seq_mean)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('roll angle [deg]')
ylim([-30 60])
set(gca,'YTick',[-90:30:180])

subplot(3,3,2)
hold on
plot(t_wb_seq_pre_mean,roll_dot_mean_wb_seq_pre_all,'.c')
plot(t_wb_seq_post_mean,roll_dot_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean,roll_dot_mean_wb_seq_pre_mean,'c')
plot(t_wb_seq_post_mean,roll_dot_mean_wb_seq_post_mean,'color',[1 .5 0])

plot(t_wb_seq_mean,roll_dot_mean_wb_seq_mean,':k')
ciplot(roll_dot_mean_wb_seq_mean-1.96*roll_dot_mean_wb_seq_ste,roll_dot_mean_wb_seq_mean+1.96*roll_dot_mean_wb_seq_ste,t_wb_seq_mean)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('roll velocity [deg/sec]')
ylim([-1500 1500])
set(gca,'YTick',[-1500:1500:1500])

subplot(3,3,3)
hold on
plot(t_wb_seq_pre_mean,roll_dot_dot_mean_wb_seq_pre_all,'.c')
plot(t_wb_seq_post_mean,roll_dot_dot_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean,roll_dot_dot_mean_wb_seq_pre_mean,'c')
plot(t_wb_seq_post_mean,roll_dot_dot_mean_wb_seq_post_mean,'color',[1 .5 0])

plot(t_wb_seq_mean,roll_dot_dot_mean_wb_seq_mean,':k')
ciplot(roll_dot_dot_mean_wb_seq_mean-1.96*roll_dot_dot_mean_wb_seq_ste,roll_dot_dot_mean_wb_seq_mean+1.96*roll_dot_dot_mean_wb_seq_ste,t_wb_seq_mean)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('roll acceleration/freq^2 [-]')
ylim([-3 2])
set(gca,'YTick',[-3:1:2])

subplot(3,3,4)
hold on
plot(t_wb_seq_pre_mean,pitch_mean_wb_seq_pre_all,'.c')
plot(t_wb_seq_post_mean,pitch_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean,pitch_mean_wb_seq_pre_mean,'c')
plot(t_wb_seq_post_mean,pitch_mean_wb_seq_post_mean,'color',[1 .5 0])

plot(t_wb_seq_mean,pitch_mean_wb_seq_mean,':k')
ciplot(pitch_mean_wb_seq_mean-1.96*pitch_mean_wb_seq_ste,pitch_mean_wb_seq_mean+1.96*pitch_mean_wb_seq_ste,t_wb_seq_mean)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('pitch angle [deg]')
ylim([-30 60])
set(gca,'YTick',[-90:30:180])

subplot(3,3,5)
hold on
plot(t_wb_seq_pre_mean,pitch_dot_mean_wb_seq_pre_all,'.c')
plot(t_wb_seq_post_mean,pitch_dot_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean,pitch_dot_mean_wb_seq_pre_mean,'c')
plot(t_wb_seq_post_mean,pitch_dot_mean_wb_seq_post_mean,'color',[1 .5 0])

plot(t_wb_seq_mean,pitch_dot_mean_wb_seq_mean,':k')
ciplot(pitch_dot_mean_wb_seq_mean-1.96*pitch_dot_mean_wb_seq_ste,pitch_dot_mean_wb_seq_mean+1.96*pitch_dot_mean_wb_seq_ste,t_wb_seq_mean)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('pitch velocity [deg/sec]')
ylim([-1500 1500])
set(gca,'YTick',[-1500:1500:1500])

subplot(3,3,6)
hold on
plot(t_wb_seq_pre_mean,pitch_dot_dot_mean_wb_seq_pre_all,'.c')
plot(t_wb_seq_post_mean,pitch_dot_dot_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean,pitch_dot_dot_mean_wb_seq_pre_mean,'c')
plot(t_wb_seq_post_mean,pitch_dot_dot_mean_wb_seq_post_mean,'color',[1 .5 0])

plot(t_wb_seq_mean,pitch_dot_dot_mean_wb_seq_mean,':k')
ciplot(pitch_dot_dot_mean_wb_seq_mean-1.96*pitch_dot_dot_mean_wb_seq_ste,pitch_dot_dot_mean_wb_seq_mean+1.96*pitch_dot_dot_mean_wb_seq_ste,t_wb_seq_mean)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('pitch acceleration/freq^2 [-]')
ylim([-3 2])
set(gca,'YTick',[-3:1:2])

subplot(3,3,7)
hold on
plot(t_wb_seq_pre_mean,yaw_mean_wb_seq_pre_all,'.c')
plot(t_wb_seq_post_mean,yaw_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean,yaw_mean_wb_seq_pre_mean,'c')
plot(t_wb_seq_post_mean,yaw_mean_wb_seq_post_mean,'color',[1 .5 0])

plot(t_wb_seq_mean,yaw_mean_wb_seq_mean,':k')
ciplot(yaw_mean_wb_seq_mean-1.96*yaw_mean_wb_seq_ste,yaw_mean_wb_seq_mean+1.96*yaw_mean_wb_seq_ste,t_wb_seq_mean)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('yaw angle [deg]')
ylim([0 90])
set(gca,'YTick',[-90:30:180])

subplot(3,3,8)
hold on
plot(t_wb_seq_pre_mean,yaw_dot_mean_wb_seq_pre_all,'.c')
plot(t_wb_seq_post_mean,yaw_dot_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean,yaw_dot_mean_wb_seq_pre_mean,'c')
plot(t_wb_seq_post_mean,yaw_dot_mean_wb_seq_post_mean,'color',[1 .5 0])

plot(t_wb_seq_mean,yaw_dot_mean_wb_seq_mean,':k')
ciplot(yaw_dot_mean_wb_seq_mean-1.96*yaw_dot_mean_wb_seq_ste,yaw_dot_mean_wb_seq_mean+1.96*yaw_dot_mean_wb_seq_ste,t_wb_seq_mean)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('yaw velocity [deg/sec]')
ylim([-1000 2000])
set(gca,'YTick',[-1000:1000:2000])

subplot(3,3,9)
hold on
plot(t_wb_seq_pre_mean,yaw_dot_dot_mean_wb_seq_pre_all,'.c')
plot(t_wb_seq_post_mean,yaw_dot_dot_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean,yaw_dot_dot_mean_wb_seq_pre_mean,'c')
plot(t_wb_seq_post_mean,yaw_dot_dot_mean_wb_seq_post_mean,'color',[1 .5 0])

plot(t_wb_seq_mean,yaw_dot_dot_mean_wb_seq_mean,':k')
ciplot(yaw_dot_dot_mean_wb_seq_mean-1.96*yaw_dot_dot_mean_wb_seq_ste,yaw_dot_dot_mean_wb_seq_mean+1.96*yaw_dot_dot_mean_wb_seq_ste,t_wb_seq_mean)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('yaw acceleration/freq^2 [-]')
ylim([-2 3])
set(gca,'YTick',[-3:1:3])

saveas(gca,['BodyDynvsTime_RollPitchYaw_allNmean.fig'])
saveas(gca,['BodyDynvsTime_RollPitchYaw_allNmean.png'])
plot2svg(['BodyDynvsTime_RollPitchYaw_allNmean.svg'])

%% rotL&R
figure
subplot(3,3,1)
hold on
plot(t_wb_seq_pre_mean,drot_R_mean_wb_seq_pre_all,'.c')
plot(t_wb_seq_post_mean,drot_R_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean,drot_R_mean_wb_seq_pre_mean,'c')
plot(t_wb_seq_post_mean,drot_R_mean_wb_seq_post_mean,'color',[1 .5 0])

plot(t_wb_seq_mean,drot_R_mean_wb_seq_mean,':k')
ciplot(drot_R_mean_wb_seq_mean-1.96*drot_R_mean_wb_seq_ste,drot_R_mean_wb_seq_mean+1.96*drot_R_mean_wb_seq_ste,t_wb_seq_mean)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('right axis rotation angle [deg]')
ylim([-30 60])
set(gca,'YTick',[-90:30:180])

subplot(3,3,2)
hold on
plot(t_wb_seq_pre_mean,rot_dot_R_mean_wb_seq_pre_all,'.c')
plot(t_wb_seq_post_mean,rot_dot_R_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean,rot_dot_R_mean_wb_seq_pre_mean,'c')
plot(t_wb_seq_post_mean,rot_dot_R_mean_wb_seq_post_mean,'color',[1 .5 0])

plot(t_wb_seq_mean,rot_dot_R_mean_wb_seq_mean,':k')
ciplot(rot_dot_R_mean_wb_seq_mean-1.96*rot_dot_R_mean_wb_seq_ste,rot_dot_R_mean_wb_seq_mean+1.96*rot_dot_R_mean_wb_seq_ste,t_wb_seq_mean)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('right axis rotation velocity [deg/sec]')
ylim([-1000 2000])
set(gca,'YTick',[-2000:1000:2000])

subplot(3,3,3)
hold on
plot(t_wb_seq_pre_mean,rot_dot_dot_R_mean_wb_seq_pre_all,'.c')
plot(t_wb_seq_post_mean,rot_dot_dot_R_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean,rot_dot_dot_R_mean_wb_seq_pre_mean,'c')
plot(t_wb_seq_post_mean,rot_dot_dot_R_mean_wb_seq_post_mean,'color',[1 .5 0])

plot(t_wb_seq_mean,rot_dot_dot_R_mean_wb_seq_mean,':k')
ciplot(rot_dot_dot_R_mean_wb_seq_mean-1.96*rot_dot_dot_R_mean_wb_seq_ste,rot_dot_dot_R_mean_wb_seq_mean+1.96*rot_dot_dot_R_mean_wb_seq_ste,t_wb_seq_mean)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('right axis rotation acceleration/freq^2 [-]')
ylim([-3 2])
set(gca,'YTick',[-3:1:3])

subplot(3,3,4)
hold on
plot(t_wb_seq_pre_mean,drot_L_mean_wb_seq_pre_all,'.c')
plot(t_wb_seq_post_mean,drot_L_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean,drot_L_mean_wb_seq_pre_mean,'c')
plot(t_wb_seq_post_mean,drot_L_mean_wb_seq_post_mean,'color',[1 .5 0])

plot(t_wb_seq_mean,drot_L_mean_wb_seq_mean,':k')
ciplot(drot_L_mean_wb_seq_mean-1.96*drot_L_mean_wb_seq_ste,drot_L_mean_wb_seq_mean+1.96*drot_L_mean_wb_seq_ste,t_wb_seq_mean)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('left axis rotation angle [deg]')
ylim([-60 30])
set(gca,'YTick',[-90:30:180])

subplot(3,3,5)
hold on
plot(t_wb_seq_pre_mean,rot_dot_L_mean_wb_seq_pre_all,'.c')
plot(t_wb_seq_post_mean,rot_dot_L_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean,rot_dot_L_mean_wb_seq_pre_mean,'c')
plot(t_wb_seq_post_mean,rot_dot_L_mean_wb_seq_post_mean,'color',[1 .5 0])

plot(t_wb_seq_mean,rot_dot_L_mean_wb_seq_mean,':k')
ciplot(rot_dot_L_mean_wb_seq_mean-1.96*rot_dot_L_mean_wb_seq_ste,rot_dot_L_mean_wb_seq_mean+1.96*rot_dot_L_mean_wb_seq_ste,t_wb_seq_mean)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('left axis rotation velocity [deg/sec]')
ylim([-2000 1000])
set(gca,'YTick',[-2000:1000:2000])

subplot(3,3,6)
hold on
plot(t_wb_seq_pre_mean,rot_dot_dot_L_mean_wb_seq_pre_all,'.c')
plot(t_wb_seq_post_mean,rot_dot_dot_L_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean,rot_dot_dot_L_mean_wb_seq_pre_mean,'c')
plot(t_wb_seq_post_mean,rot_dot_dot_L_mean_wb_seq_post_mean,'color',[1 .5 0])

plot(t_wb_seq_mean,rot_dot_dot_L_mean_wb_seq_mean,':k')
ciplot(rot_dot_dot_L_mean_wb_seq_mean-1.96*rot_dot_dot_L_mean_wb_seq_ste,rot_dot_dot_L_mean_wb_seq_mean+1.96*rot_dot_dot_L_mean_wb_seq_ste,t_wb_seq_mean)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('left axis rotation acceleration/freq^2 [-]')
ylim([-3 2])
set(gca,'YTick',[-3:1:3])

subplot(3,3,7)
hold on
plot(t_wb_seq_pre_mean,yaw_mean_wb_seq_pre_all,'.c')
plot(t_wb_seq_post_mean,yaw_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean,yaw_mean_wb_seq_pre_mean,'c')
plot(t_wb_seq_post_mean,yaw_mean_wb_seq_post_mean,'color',[1 .5 0])

plot(t_wb_seq_mean,yaw_mean_wb_seq_mean,':k')
ciplot(yaw_mean_wb_seq_mean-1.96*yaw_mean_wb_seq_ste,yaw_mean_wb_seq_mean+1.96*yaw_mean_wb_seq_ste,t_wb_seq_mean)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('yaw angle [deg]')
ylim([0 90])
set(gca,'YTick',[-90:30:180])

subplot(3,3,8)
hold on
plot(t_wb_seq_pre_mean,yaw_dot_mean_wb_seq_pre_all,'.c')
plot(t_wb_seq_post_mean,yaw_dot_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean,yaw_dot_mean_wb_seq_pre_mean,'c')
plot(t_wb_seq_post_mean,yaw_dot_mean_wb_seq_post_mean,'color',[1 .5 0])

plot(t_wb_seq_mean,yaw_dot_mean_wb_seq_mean,':k')
ciplot(yaw_dot_mean_wb_seq_mean-1.96*yaw_dot_mean_wb_seq_ste,yaw_dot_mean_wb_seq_mean+1.96*yaw_dot_mean_wb_seq_ste,t_wb_seq_mean)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('yaw velocity [deg/sec]')
ylim([-1000 2000])
set(gca,'YTick',[-1000:1000:2000])

subplot(3,3,9)
hold on
plot(t_wb_seq_pre_mean,yaw_dot_dot_mean_wb_seq_pre_all,'.c')
plot(t_wb_seq_post_mean,yaw_dot_dot_mean_wb_seq_post_all,'.','color',[1 .5 0])

plot(t_wb_seq_pre_mean,yaw_dot_dot_mean_wb_seq_pre_mean,'c')
plot(t_wb_seq_post_mean,yaw_dot_dot_mean_wb_seq_post_mean,'color',[1 .5 0])

plot(t_wb_seq_mean,yaw_dot_dot_mean_wb_seq_mean,':k')
ciplot(yaw_dot_dot_mean_wb_seq_mean-1.96*yaw_dot_dot_mean_wb_seq_ste,yaw_dot_dot_mean_wb_seq_mean+1.96*yaw_dot_dot_mean_wb_seq_ste,t_wb_seq_mean)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('yaw acceleration/freq^2 [-]')
ylim([-2 3])
set(gca,'YTick',[-3:1:3])

saveas(gca,['BodyDynvsTime_rotLnRnYaw_allNmean.fig'])
saveas(gca,['BodyDynvsTime_rotLnRnYaw_allNmean.png'])
plot2svg(['BodyDynvsTime_rotLnRnYaw_allNmean.svg'])
