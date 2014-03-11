% normalized time
% wb kin
figure
hold on
plot(t_norm_wb_seq_pre_mean,f_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,f_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,f_wb_seq_mean,':k')

figure
hold on
plot(t_norm_wb_seq_bins_pre_mean,stroke_wb_L_seq_bins_pre_mean,'b')
plot(t_norm_wb_seq_bins_post_mean,stroke_wb_L_seq_bins_post_mean,'c')

plot(t_norm_wb_seq_bins_pre_mean,stroke_wb_R_seq_bins_pre_mean,'r')
plot(t_norm_wb_seq_bins_post_mean,stroke_wb_R_seq_bins_post_mean,'m')

figure
hold on
plot(t_norm_wb_seq_bins_mean,stroke_wb_L_seq_bins_mean,'b')
plot(t_norm_wb_seq_bins_mean,stroke_wb_R_seq_bins_mean,'r')

figure
hold on
plot(t_norm_wb_seq_bins_pre_mean,dev_wb_L_seq_bins_pre_mean,'b')
plot(t_norm_wb_seq_bins_post_mean,dev_wb_L_seq_bins_post_mean,'c')

plot(t_norm_wb_seq_bins_pre_mean,dev_wb_R_seq_bins_pre_mean,'r')
plot(t_norm_wb_seq_bins_post_mean,dev_wb_R_seq_bins_post_mean,'m')

figure
hold on
plot(t_norm_wb_seq_bins_mean,dev_wb_L_seq_bins_mean,'b')
plot(t_norm_wb_seq_bins_mean,dev_wb_R_seq_bins_mean,'r')

figure
hold on
plot(t_norm_wb_seq_bins_pre_mean,pitch_wb_L_seq_bins_pre_mean,'b')
plot(t_norm_wb_seq_bins_post_mean,pitch_wb_L_seq_bins_post_mean,'c')

plot(t_norm_wb_seq_bins_pre_mean,pitch_wb_R_seq_bins_pre_mean,'r')
plot(t_norm_wb_seq_bins_post_mean,pitch_wb_R_seq_bins_post_mean,'m')

figure
hold on
plot(t_norm_wb_seq_bins_mean,pitch_wb_L_seq_bins_mean,'b')
plot(t_norm_wb_seq_bins_mean,pitch_wb_R_seq_bins_mean,'r')

% body kin
figure
subplot(3,3,1)
hold on
plot(t_norm_wb_seq_pre_mean,roll_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,roll_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,roll_mean_wb_seq_mean,':k')

subplot(3,3,2)
hold on
plot(t_norm_wb_seq_pre_mean,roll_dot_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,roll_dot_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,roll_dot_mean_wb_seq_mean,':k')

subplot(3,3,3)
hold on
plot(t_norm_wb_seq_pre_mean,roll_dot_dot_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,roll_dot_dot_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,roll_dot_dot_mean_wb_seq_mean,':k')


subplot(3,3,4)
hold on
plot(t_norm_wb_seq_pre_mean,pitch_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,pitch_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,pitch_mean_wb_seq_mean,':k')

subplot(3,3,5)
hold on
plot(t_norm_wb_seq_pre_mean,pitch_dot_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,pitch_dot_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,pitch_dot_mean_wb_seq_mean,':k')

subplot(3,3,6)
hold on
plot(t_norm_wb_seq_pre_mean,pitch_dot_dot_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,pitch_dot_dot_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,pitch_dot_dot_mean_wb_seq_mean,':k')


subplot(3,3,7)
hold on
plot(t_norm_wb_seq_pre_mean,yaw_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,yaw_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,yaw_mean_wb_seq_mean,':k')

subplot(3,3,8)
hold on
plot(t_norm_wb_seq_pre_mean,yaw_dot_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,yaw_dot_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,yaw_dot_mean_wb_seq_mean,':k')

subplot(3,3,9)
hold on
plot(t_norm_wb_seq_pre_mean,yaw_dot_dot_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,yaw_dot_dot_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,yaw_dot_dot_mean_wb_seq_mean,':k')

figure
subplot(3,3,1)
hold on
plot(t_norm_wb_seq_pre_mean,drot_R_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,drot_R_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,drot_R_mean_wb_seq_mean,':k')

subplot(3,3,2)
hold on
plot(t_norm_wb_seq_pre_mean,rot_dot_R_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,rot_dot_R_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,rot_dot_R_mean_wb_seq_mean,':k')

subplot(3,3,3)
hold on
plot(t_norm_wb_seq_pre_mean,rot_dot_dot_R_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,rot_dot_dot_R_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,rot_dot_dot_R_mean_wb_seq_mean,':k')

subplot(3,3,4)
hold on
plot(t_norm_wb_seq_pre_mean,drot_L_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,drot_L_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,drot_L_mean_wb_seq_mean,':k')

subplot(3,3,5)
hold on
plot(t_norm_wb_seq_pre_mean,rot_dot_L_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,rot_dot_L_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,rot_dot_L_mean_wb_seq_mean,':k')

subplot(3,3,6)
hold on
plot(t_norm_wb_seq_pre_mean,rot_dot_dot_L_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,rot_dot_dot_L_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,rot_dot_dot_L_mean_wb_seq_mean,':k')

subplot(3,3,7)
hold on
plot(t_norm_wb_seq_pre_mean,yaw_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,yaw_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,yaw_mean_wb_seq_mean,':k')

subplot(3,3,8)
hold on
plot(t_norm_wb_seq_pre_mean,yaw_dot_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,yaw_dot_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,yaw_dot_mean_wb_seq_mean,':k')

subplot(3,3,9)
hold on
plot(t_norm_wb_seq_pre_mean,yaw_dot_dot_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,yaw_dot_dot_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,yaw_dot_dot_mean_wb_seq_mean,':k')


%% MEAN DATA: real time
% wb kin
figure
hold on
plot(t_wb_seq_pre_mean,f_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,f_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,f_wb_seq_mean,':k')

figure
hold on
plot(t_wb_seq_bins_pre_mean,stroke_wb_L_seq_bins_pre_mean,'b')
plot(t_wb_seq_bins_post_mean,stroke_wb_L_seq_bins_post_mean,'c')

plot(t_wb_seq_bins_pre_mean,stroke_wb_R_seq_bins_pre_mean,'r')
plot(t_wb_seq_bins_post_mean,stroke_wb_R_seq_bins_post_mean,'m')

figure
hold on
plot(t_wb_seq_bins_mean,stroke_wb_L_seq_bins_mean,'b')
plot(t_wb_seq_bins_mean,stroke_wb_R_seq_bins_mean,'r')

figure
hold on
plot(t_wb_seq_bins_pre_mean,dev_wb_L_seq_bins_pre_mean,'b')
plot(t_wb_seq_bins_post_mean,dev_wb_L_seq_bins_post_mean,'c')

plot(t_wb_seq_bins_pre_mean,dev_wb_R_seq_bins_pre_mean,'r')
plot(t_wb_seq_bins_post_mean,dev_wb_R_seq_bins_post_mean,'m')

figure
hold on
plot(t_wb_seq_bins_mean,dev_wb_L_seq_bins_mean,'b')
plot(t_wb_seq_bins_mean,dev_wb_R_seq_bins_mean,'r')

figure
hold on
plot(t_wb_seq_bins_pre_mean,pitch_wb_L_seq_bins_pre_mean,'b')
plot(t_wb_seq_bins_post_mean,pitch_wb_L_seq_bins_post_mean,'c')

plot(t_wb_seq_bins_pre_mean,pitch_wb_R_seq_bins_pre_mean,'r')
plot(t_wb_seq_bins_post_mean,pitch_wb_R_seq_bins_post_mean,'m')

figure
hold on
plot(t_wb_seq_bins_mean,pitch_wb_L_seq_bins_mean,'b')
plot(t_wb_seq_bins_mean,pitch_wb_R_seq_bins_mean,'r')

% body kin
figure
subplot(3,3,1)
hold on
plot(t_wb_seq_pre_mean,roll_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,roll_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,roll_mean_wb_seq_mean,':k')

subplot(3,3,2)
hold on
plot(t_wb_seq_pre_mean,roll_dot_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,roll_dot_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,roll_dot_mean_wb_seq_mean,':k')

subplot(3,3,3)
hold on
plot(t_wb_seq_pre_mean,roll_dot_dot_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,roll_dot_dot_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,roll_dot_dot_mean_wb_seq_mean,':k')


subplot(3,3,4)
hold on
plot(t_wb_seq_pre_mean,pitch_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,pitch_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,pitch_mean_wb_seq_mean,':k')

subplot(3,3,5)
hold on
plot(t_wb_seq_pre_mean,pitch_dot_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,pitch_dot_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,pitch_dot_mean_wb_seq_mean,':k')

subplot(3,3,6)
hold on
plot(t_wb_seq_pre_mean,pitch_dot_dot_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,pitch_dot_dot_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,pitch_dot_dot_mean_wb_seq_mean,':k')


subplot(3,3,7)
hold on
plot(t_wb_seq_pre_mean,yaw_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,yaw_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,yaw_mean_wb_seq_mean,':k')

subplot(3,3,8)
hold on
plot(t_wb_seq_pre_mean,yaw_dot_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,yaw_dot_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,yaw_dot_mean_wb_seq_mean,':k')

subplot(3,3,9)
hold on
plot(t_wb_seq_pre_mean,yaw_dot_dot_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,yaw_dot_dot_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,yaw_dot_dot_mean_wb_seq_mean,':k')

figure
subplot(3,3,1)
hold on
plot(t_wb_seq_pre_mean,drot_R_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,drot_R_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,drot_R_mean_wb_seq_mean,':k')

subplot(3,3,2)
hold on
plot(t_wb_seq_pre_mean,rot_dot_R_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,rot_dot_R_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,rot_dot_R_mean_wb_seq_mean,':k')

subplot(3,3,3)
hold on
plot(t_wb_seq_pre_mean,rot_dot_dot_R_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,rot_dot_dot_R_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,rot_dot_dot_R_mean_wb_seq_mean,':k')

subplot(3,3,4)
hold on
plot(t_wb_seq_pre_mean,drot_L_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,drot_L_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,drot_L_mean_wb_seq_mean,':k')

subplot(3,3,5)
hold on
plot(t_wb_seq_pre_mean,rot_dot_L_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,rot_dot_L_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,rot_dot_L_mean_wb_seq_mean,':k')

subplot(3,3,6)
hold on
plot(t_wb_seq_pre_mean,rot_dot_dot_L_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,rot_dot_dot_L_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,rot_dot_dot_L_mean_wb_seq_mean,':k')

subplot(3,3,7)
hold on
plot(t_wb_seq_pre_mean,yaw_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,yaw_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,yaw_mean_wb_seq_mean,':k')

subplot(3,3,8)
hold on
plot(t_wb_seq_pre_mean,yaw_dot_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,yaw_dot_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,yaw_dot_mean_wb_seq_mean,':k')

subplot(3,3,9)
hold on
plot(t_wb_seq_pre_mean,yaw_dot_dot_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,yaw_dot_dot_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,yaw_dot_dot_mean_wb_seq_mean,':k')


