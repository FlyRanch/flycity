% wb kin

% freq
figure
hold on
plot(t_wb_seq_pre_mean,f_wb_seq_pre,'.r')
plot(t_wb_seq_post_mean,f_wb_seq_post,'.c')

plot(t_wb_seq_pre_mean,f_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,f_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,f_wb_seq_mean,':k')

xlabel('time [sec]')
ylabel('freq [1/sec]')
ylim([150 230])

% angles
figure
subplot(3,1,1)
hold on
plot(t_wb_seq_bins_mean,stroke_wb_L_seq_bins_all,'b','linewidth',.5)
plot(t_wb_seq_bins_mean,stroke_wb_R_seq_bins_all,'r','linewidth',.5)

plot(t_wb_seq_bins_mean,stroke_wb_L_seq_bins_mean,'c','linewidth',2)
plot(t_wb_seq_bins_mean,stroke_wb_R_seq_bins_mean,'y','linewidth',2)

xlabel('time [sec]')
ylabel('stroke angle [deg]')
ylim([-90 90])
set(gca,'YTick',[-90:45:90])

subplot(3,1,2)
hold on
plot(t_wb_seq_bins_mean,dev_wb_L_seq_bins_all,'b','linewidth',.5)
plot(t_wb_seq_bins_mean,dev_wb_R_seq_bins_all,'r','linewidth',.5)

plot(t_wb_seq_bins_mean,dev_wb_L_seq_bins_mean,'c','linewidth',2)
plot(t_wb_seq_bins_mean,dev_wb_R_seq_bins_mean,'y','linewidth',2)

xlabel('time [sec]')
ylabel('deviation angle [deg]')
ylim([-30 30])
set(gca,'YTick',[-90:30:90])

subplot(3,1,3)
hold on
plot(t_wb_seq_bins_mean,pitch_wb_L_seq_bins_all,'b','linewidth',.5)
plot(t_wb_seq_bins_mean,pitch_wb_R_seq_bins_all,'r','linewidth',.5)

plot(t_wb_seq_bins_mean,pitch_wb_L_seq_bins_mean,'c','linewidth',2)
plot(t_wb_seq_bins_mean,pitch_wb_R_seq_bins_mean,'y','linewidth',2)

xlabel('time [sec]')
ylabel('rotation angle [deg]')
ylim([0 180])
set(gca,'YTick',[-90:45:180])

%% body kin
% roll, pitch
figure
subplot(3,3,1)
hold on
plot(t_wb_seq_pre_mean,roll_mean_wb_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,roll_mean_wb_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,roll_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,roll_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,roll_mean_wb_seq_mean,':k')

xlabel('time [sec]')
ylabel('roll angle [deg]')
ylim([-30 60])
set(gca,'YTick',[-90:30:180])

subplot(3,3,2)
hold on
plot(t_wb_seq_pre_mean,roll_dot_mean_wb_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,roll_dot_mean_wb_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,roll_dot_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,roll_dot_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,roll_dot_mean_wb_seq_mean,':k')

xlabel('time [sec]')
ylabel('roll velocity [deg/sec]')
ylim([-1500 1500])
set(gca,'YTick',[-1500:1500:1500])

subplot(3,3,3)
hold on
plot(t_wb_seq_pre_mean,roll_dot_dot_mean_wb_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,roll_dot_dot_mean_wb_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,roll_dot_dot_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,roll_dot_dot_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,roll_dot_dot_mean_wb_seq_mean,':k')

xlabel('time [sec]')
ylabel('roll acceleration/freq^2 [-]')
ylim([-3 2])
set(gca,'YTick',[-3:1:2])

subplot(3,3,4)
hold on
plot(t_wb_seq_pre_mean,pitch_mean_wb_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,pitch_mean_wb_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,pitch_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,pitch_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,pitch_mean_wb_seq_mean,':k')

xlabel('time [sec]')
ylabel('pitch angle [deg]')
ylim([-30 60])
set(gca,'YTick',[-90:30:180])

subplot(3,3,5)
hold on
plot(t_wb_seq_pre_mean,pitch_dot_mean_wb_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,pitch_dot_mean_wb_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,pitch_dot_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,pitch_dot_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,pitch_dot_mean_wb_seq_mean,':k')

xlabel('time [sec]')
ylabel('pitch velocity [deg/sec]')
ylim([-1500 1500])
set(gca,'YTick',[-1500:1500:1500])

subplot(3,3,6)
hold on
plot(t_wb_seq_pre_mean,pitch_dot_dot_mean_wb_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,pitch_dot_dot_mean_wb_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,pitch_dot_dot_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,pitch_dot_dot_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,pitch_dot_dot_mean_wb_seq_mean,':k')

xlabel('time [sec]')
ylabel('pitch acceleration/freq^2 [-]')
ylim([-3 2])
set(gca,'YTick',[-3:1:2])

subplot(3,3,7)
hold on
plot(t_wb_seq_pre_mean,yaw_mean_wb_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,yaw_mean_wb_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,yaw_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,yaw_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,yaw_mean_wb_seq_mean,':k')

xlabel('time [sec]')
ylabel('yaw angle [deg]')
ylim([0 90])
set(gca,'YTick',[-90:30:180])

subplot(3,3,8)
hold on
plot(t_wb_seq_pre_mean,yaw_dot_mean_wb_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,yaw_dot_mean_wb_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,yaw_dot_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,yaw_dot_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,yaw_dot_mean_wb_seq_mean,':k')

xlabel('time [sec]')
ylabel('yaw velocity [deg/sec]')
ylim([-1000 2000])
set(gca,'YTick',[-1000:1000:2000])

subplot(3,3,9)
hold on
plot(t_wb_seq_pre_mean,yaw_dot_dot_mean_wb_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,yaw_dot_dot_mean_wb_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,yaw_dot_dot_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,yaw_dot_dot_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,yaw_dot_dot_mean_wb_seq_mean,':k')

xlabel('time [sec]')
ylabel('yaw acceleration/freq^2 [-]')
ylim([-2 3])
set(gca,'YTick',[-3:1:3])

% rotL&R
figure
subplot(3,3,1)
hold on
plot(t_wb_seq_pre_mean,drot_R_mean_wb_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,drot_R_mean_wb_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,drot_R_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,drot_R_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,drot_R_mean_wb_seq_mean,':k')

xlabel('time [sec]')
ylabel('right axis rotation angle [deg]')
ylim([-30 60])
set(gca,'YTick',[-90:30:180])

subplot(3,3,2)
hold on
plot(t_wb_seq_pre_mean,rot_dot_R_mean_wb_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,rot_dot_R_mean_wb_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,rot_dot_R_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,rot_dot_R_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,rot_dot_R_mean_wb_seq_mean,':k')

xlabel('time [sec]')
ylabel('right axis rotation velocity [deg/sec]')
ylim([-1000 2000])
set(gca,'YTick',[-2000:1000:2000])

subplot(3,3,3)
hold on
plot(t_wb_seq_pre_mean,rot_dot_dot_R_mean_wb_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,rot_dot_dot_R_mean_wb_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,rot_dot_dot_R_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,rot_dot_dot_R_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,rot_dot_dot_R_mean_wb_seq_mean,':k')

xlabel('time [sec]')
ylabel('right axis rotation acceleration/freq^2 [-]')
ylim([-3 2])
set(gca,'YTick',[-3:1:3])

subplot(3,3,4)
hold on
plot(t_wb_seq_pre_mean,drot_L_mean_wb_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,drot_L_mean_wb_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,drot_L_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,drot_L_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,drot_L_mean_wb_seq_mean,':k')

xlabel('time [sec]')
ylabel('left axis rotation angle [deg]')
ylim([-60 30])
set(gca,'YTick',[-90:30:180])

subplot(3,3,5)
hold on
plot(t_wb_seq_pre_mean,rot_dot_L_mean_wb_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,rot_dot_L_mean_wb_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,rot_dot_L_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,rot_dot_L_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,rot_dot_L_mean_wb_seq_mean,':k')

xlabel('time [sec]')
ylabel('left axis rotation velocity [deg/sec]')
ylim([-2000 1000])
set(gca,'YTick',[-2000:1000:2000])

subplot(3,3,6)
hold on
plot(t_wb_seq_pre_mean,rot_dot_dot_L_mean_wb_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,rot_dot_dot_L_mean_wb_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,rot_dot_dot_L_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,rot_dot_dot_L_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,rot_dot_dot_L_mean_wb_seq_mean,':k')

xlabel('time [sec]')
ylabel('left axis rotation acceleration/freq^2 [-]')
ylim([-3 2])
set(gca,'YTick',[-3:1:3])

subplot(3,3,7)
hold on
plot(t_wb_seq_pre_mean,yaw_mean_wb_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,yaw_mean_wb_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,yaw_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,yaw_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,yaw_mean_wb_seq_mean,':k')

xlabel('time [sec]')
ylabel('yaw angle [deg]')
ylim([0 90])
set(gca,'YTick',[-90:30:180])

subplot(3,3,8)
hold on
plot(t_wb_seq_pre_mean,yaw_dot_mean_wb_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,yaw_dot_mean_wb_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,yaw_dot_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,yaw_dot_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,yaw_dot_mean_wb_seq_mean,':k')

xlabel('time [sec]')
ylabel('yaw velocity [deg/sec]')
ylim([-1000 2000])
set(gca,'YTick',[-1000:1000:2000])

subplot(3,3,9)
hold on
plot(t_wb_seq_pre_mean,yaw_dot_dot_mean_wb_seq_pre_all,'.r')
plot(t_wb_seq_post_mean,yaw_dot_dot_mean_wb_seq_post_all,'.c')

plot(t_wb_seq_pre_mean,yaw_dot_dot_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,yaw_dot_dot_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,yaw_dot_dot_mean_wb_seq_mean,':k')

xlabel('time [sec]')
ylabel('yaw acceleration/freq^2 [-]')
ylim([-2 3])
set(gca,'YTick',[-3:1:3])

