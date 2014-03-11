% define timeline variables

%% ALL

%% wb kin vars
f_wb_seq_pre = nan(max(seq_nr),wb_max);
f_wb_seq_post = nan(max(seq_nr),wb_max);

stroke_wb_L_seq_bins_pre = nan(n_bins,max(seq_nr),wb_max);
dev_wb_L_seq_bins_pre = nan(n_bins,max(seq_nr),wb_max);
pitch_wb_L_seq_bins_pre = nan(n_bins,max(seq_nr),wb_max);

stroke_wb_R_seq_bins_pre = nan(n_bins,max(seq_nr),wb_max);
dev_wb_R_seq_bins_pre = nan(n_bins,max(seq_nr),wb_max);
pitch_wb_R_seq_bins_pre = nan(n_bins,max(seq_nr),wb_max);

stroke_wb_L_seq_bins_post = nan(n_bins,max(seq_nr),wb_max);
dev_wb_L_seq_bins_post = nan(n_bins,max(seq_nr),wb_max);
pitch_wb_L_seq_bins_post = nan(n_bins,max(seq_nr),wb_max);

stroke_wb_R_seq_bins_post = nan(n_bins,max(seq_nr),wb_max);
dev_wb_R_seq_bins_post = nan(n_bins,max(seq_nr),wb_max);
pitch_wb_R_seq_bins_post = nan(n_bins,max(seq_nr),wb_max);

%% body kin vars
% pre
roll_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
roll_dot_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
roll_dot_dot_mean_wb_seq_pre = nan(max(seq_nr),wb_max);

pitch_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
pitch_dot_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
pitch_dot_dot_mean_wb_seq_pre = nan(max(seq_nr),wb_max);

yaw_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
yaw_dot_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
yaw_dot_dot_mean_wb_seq_pre = nan(max(seq_nr),wb_max);

% drot_L_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
% rot_L_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
% rot_dot_L_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
% rot_dot_dot_L_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
% 
% drot_R_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
% rot_R_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
% rot_dot_R_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
% rot_dot_dot_R_mean_wb_seq_pre = nan(max(seq_nr),wb_max);

% post
roll_mean_wb_seq_post = nan(max(seq_nr),wb_max);
roll_dot_mean_wb_seq_post = nan(max(seq_nr),wb_max);
roll_dot_dot_mean_wb_seq_post = nan(max(seq_nr),wb_max);

pitch_mean_wb_seq_post = nan(max(seq_nr),wb_max);
pitch_dot_mean_wb_seq_post = nan(max(seq_nr),wb_max);
pitch_dot_dot_mean_wb_seq_post = nan(max(seq_nr),wb_max);

yaw_mean_wb_seq_post = nan(max(seq_nr),wb_max);
yaw_dot_mean_wb_seq_post = nan(max(seq_nr),wb_max);
yaw_dot_dot_mean_wb_seq_post = nan(max(seq_nr),wb_max);

% drot_L_mean_wb_seq_post = nan(max(seq_nr),wb_max);
% rot_L_mean_wb_seq_post = nan(max(seq_nr),wb_max);
% rot_dot_L_mean_wb_seq_post = nan(max(seq_nr),wb_max);
% rot_dot_dot_L_mean_wb_seq_post = nan(max(seq_nr),wb_max);
% 
% drot_R_mean_wb_seq_post = nan(max(seq_nr),wb_max);
% rot_R_mean_wb_seq_post = nan(max(seq_nr),wb_max);
% rot_dot_R_mean_wb_seq_post = nan(max(seq_nr),wb_max);
% rot_dot_dot_R_mean_wb_seq_post = nan(max(seq_nr),wb_max);

%% Torque vars
% % pre
% Mroll_mean_wb_accel_seq_pre = nan(max(seq_nr),wb_max);
% Mroll_mean_wb_damp_seq_pre = nan(max(seq_nr),wb_max);
% Mroll_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
% 
% Mpitch_mean_wb_accel_seq_pre = nan(max(seq_nr),wb_max);
% Mpitch_mean_wb_damp_seq_pre = nan(max(seq_nr),wb_max);
% Mpitch_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
% 
% Myaw_mean_wb_accel_seq_pre = nan(max(seq_nr),wb_max);
% Myaw_mean_wb_damp_seq_pre = nan(max(seq_nr),wb_max);
% Myaw_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
% 
% M_L_mean_wb_accel_seq_pre = nan(max(seq_nr),wb_max);
% M_L_mean_wb_damp_seq_pre = nan(max(seq_nr),wb_max);
% M_L_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
% 
% M_R_mean_wb_accel_seq_pre = nan(max(seq_nr),wb_max);
% M_R_mean_wb_damp_seq_pre = nan(max(seq_nr),wb_max);
% M_R_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
% 
% % post
% Mroll_mean_wb_accel_seq_post = nan(max(seq_nr),wb_max);
% Mroll_mean_wb_damp_seq_post = nan(max(seq_nr),wb_max);
% Mroll_mean_wb_seq_post = nan(max(seq_nr),wb_max);
% 
% Mpitch_mean_wb_accel_seq_post = nan(max(seq_nr),wb_max);
% Mpitch_mean_wb_damp_seq_post = nan(max(seq_nr),wb_max);
% Mpitch_mean_wb_seq_post = nan(max(seq_nr),wb_max);
% 
% Myaw_mean_wb_accel_seq_post = nan(max(seq_nr),wb_max);
% Myaw_mean_wb_damp_seq_post = nan(max(seq_nr),wb_max);
% Myaw_mean_wb_seq_post = nan(max(seq_nr),wb_max);
% 
% M_L_mean_wb_accel_seq_post = nan(max(seq_nr),wb_max);
% M_L_mean_wb_damp_seq_post = nan(max(seq_nr),wb_max);
% M_L_mean_wb_seq_post = nan(max(seq_nr),wb_max);
% 
% M_R_mean_wb_accel_seq_post = nan(max(seq_nr),wb_max);
% M_R_mean_wb_damp_seq_post = nan(max(seq_nr),wb_max);
% M_R_mean_wb_seq_post = nan(max(seq_nr),wb_max);


%% CLUSTER 1

%% wb kin vars
f_wb_seq_pre_c1 = nan(max(seq_nr),wb_max);
f_wb_seq_post_c1 = nan(max(seq_nr),wb_max);

stroke_wb_L_seq_bins_pre_c1 = nan(n_bins,max(seq_nr),wb_max);
dev_wb_L_seq_bins_pre_c1 = nan(n_bins,max(seq_nr),wb_max);
pitch_wb_L_seq_bins_pre_c1 = nan(n_bins,max(seq_nr),wb_max);

stroke_wb_R_seq_bins_pre_c1 = nan(n_bins,max(seq_nr),wb_max);
dev_wb_R_seq_bins_pre_c1 = nan(n_bins,max(seq_nr),wb_max);
pitch_wb_R_seq_bins_pre_c1 = nan(n_bins,max(seq_nr),wb_max);

stroke_wb_L_seq_bins_post_c1 = nan(n_bins,max(seq_nr),wb_max);
dev_wb_L_seq_bins_post_c1 = nan(n_bins,max(seq_nr),wb_max);
pitch_wb_L_seq_bins_post_c1 = nan(n_bins,max(seq_nr),wb_max);

stroke_wb_R_seq_bins_post_c1 = nan(n_bins,max(seq_nr),wb_max);
dev_wb_R_seq_bins_post_c1 = nan(n_bins,max(seq_nr),wb_max);
pitch_wb_R_seq_bins_post_c1 = nan(n_bins,max(seq_nr),wb_max);

%% body kin vars
% pre
roll_mean_wb_seq_pre_c1 = nan(max(seq_nr),wb_max);
roll_dot_mean_wb_seq_pre_c1 = nan(max(seq_nr),wb_max);
roll_dot_dot_mean_wb_seq_pre_c1 = nan(max(seq_nr),wb_max);

pitch_mean_wb_seq_pre_c1 = nan(max(seq_nr),wb_max);
pitch_dot_mean_wb_seq_pre_c1 = nan(max(seq_nr),wb_max);
pitch_dot_dot_mean_wb_seq_pre_c1 = nan(max(seq_nr),wb_max);

yaw_mean_wb_seq_pre_c1 = nan(max(seq_nr),wb_max);
yaw_dot_mean_wb_seq_pre_c1 = nan(max(seq_nr),wb_max);
yaw_dot_dot_mean_wb_seq_pre_c1 = nan(max(seq_nr),wb_max);

% drot_L_mean_wb_seq_pre_c1 = nan(max(seq_nr),wb_max);
% rot_L_mean_wb_seq_pre_c1 = nan(max(seq_nr),wb_max);
% rot_dot_L_mean_wb_seq_pre_c1 = nan(max(seq_nr),wb_max);
% rot_dot_dot_L_mean_wb_seq_pre_c1 = nan(max(seq_nr),wb_max);
% 
% drot_R_mean_wb_seq_pre_c1 = nan(max(seq_nr),wb_max);
% rot_R_mean_wb_seq_pre_c1 = nan(max(seq_nr),wb_max);
% rot_dot_R_mean_wb_seq_pre_c1 = nan(max(seq_nr),wb_max);
% rot_dot_dot_R_mean_wb_seq_pre_c1 = nan(max(seq_nr),wb_max);

% post
roll_mean_wb_seq_post_c1 = nan(max(seq_nr),wb_max);
roll_dot_mean_wb_seq_post_c1 = nan(max(seq_nr),wb_max);
roll_dot_dot_mean_wb_seq_post_c1 = nan(max(seq_nr),wb_max);

pitch_mean_wb_seq_post_c1 = nan(max(seq_nr),wb_max);
pitch_dot_mean_wb_seq_post_c1 = nan(max(seq_nr),wb_max);
pitch_dot_dot_mean_wb_seq_post_c1 = nan(max(seq_nr),wb_max);

yaw_mean_wb_seq_post_c1 = nan(max(seq_nr),wb_max);
yaw_dot_mean_wb_seq_post_c1 = nan(max(seq_nr),wb_max);
yaw_dot_dot_mean_wb_seq_post_c1 = nan(max(seq_nr),wb_max);

% drot_L_mean_wb_seq_post_c1 = nan(max(seq_nr),wb_max);
% rot_L_mean_wb_seq_post_c1 = nan(max(seq_nr),wb_max);
% rot_dot_L_mean_wb_seq_post_c1 = nan(max(seq_nr),wb_max);
% rot_dot_dot_L_mean_wb_seq_post_c1 = nan(max(seq_nr),wb_max);
% 
% drot_R_mean_wb_seq_post_c1 = nan(max(seq_nr),wb_max);
% rot_R_mean_wb_seq_post_c1 = nan(max(seq_nr),wb_max);
% rot_dot_R_mean_wb_seq_post_c1 = nan(max(seq_nr),wb_max);
% rot_dot_dot_R_mean_wb_seq_post_c1 = nan(max(seq_nr),wb_max);

%% Torque vars
% % pre
% Mroll_mean_wb_accel_seq_pre_c1 = nan(max(seq_nr),wb_max);
% Mroll_mean_wb_damp_seq_pre_c1 = nan(max(seq_nr),wb_max);
% Mroll_mean_wb_seq_pre_c1 = nan(max(seq_nr),wb_max);
% 
% Mpitch_mean_wb_accel_seq_pre_c1 = nan(max(seq_nr),wb_max);
% Mpitch_mean_wb_damp_seq_pre_c1 = nan(max(seq_nr),wb_max);
% Mpitch_mean_wb_seq_pre_c1 = nan(max(seq_nr),wb_max);
% 
% Myaw_mean_wb_accel_seq_pre_c1 = nan(max(seq_nr),wb_max);
% Myaw_mean_wb_damp_seq_pre_c1 = nan(max(seq_nr),wb_max);
% Myaw_mean_wb_seq_pre_c1 = nan(max(seq_nr),wb_max);
% 
% M_L_mean_wb_accel_seq_pre_c1 = nan(max(seq_nr),wb_max);
% M_L_mean_wb_damp_seq_pre_c1 = nan(max(seq_nr),wb_max);
% M_L_mean_wb_seq_pre_c1 = nan(max(seq_nr),wb_max);
% 
% M_R_mean_wb_accel_seq_pre_c1 = nan(max(seq_nr),wb_max);
% M_R_mean_wb_damp_seq_pre_c1 = nan(max(seq_nr),wb_max);
% M_R_mean_wb_seq_pre_c1 = nan(max(seq_nr),wb_max);
% 
% % post
% Mroll_mean_wb_accel_seq_post_c1 = nan(max(seq_nr),wb_max);
% Mroll_mean_wb_damp_seq_post_c1 = nan(max(seq_nr),wb_max);
% Mroll_mean_wb_seq_post_c1 = nan(max(seq_nr),wb_max);
% 
% Mpitch_mean_wb_accel_seq_post_c1 = nan(max(seq_nr),wb_max);
% Mpitch_mean_wb_damp_seq_post_c1 = nan(max(seq_nr),wb_max);
% Mpitch_mean_wb_seq_post_c1 = nan(max(seq_nr),wb_max);
% 
% Myaw_mean_wb_accel_seq_post_c1 = nan(max(seq_nr),wb_max);
% Myaw_mean_wb_damp_seq_post_c1 = nan(max(seq_nr),wb_max);
% Myaw_mean_wb_seq_post_c1 = nan(max(seq_nr),wb_max);
% 
% M_L_mean_wb_accel_seq_post_c1 = nan(max(seq_nr),wb_max);
% M_L_mean_wb_damp_seq_post_c1 = nan(max(seq_nr),wb_max);
% M_L_mean_wb_seq_post_c1 = nan(max(seq_nr),wb_max);
% 
% M_R_mean_wb_accel_seq_post_c1 = nan(max(seq_nr),wb_max);
% M_R_mean_wb_damp_seq_post_c1 = nan(max(seq_nr),wb_max);
% M_R_mean_wb_seq_post_c1 = nan(max(seq_nr),wb_max);

%% CLUSTER 2

%% wb kin vars
f_wb_seq_pre_c2 = nan(max(seq_nr),wb_max);
f_wb_seq_post_c2 = nan(max(seq_nr),wb_max);

stroke_wb_L_seq_bins_pre_c2 = nan(n_bins,max(seq_nr),wb_max);
dev_wb_L_seq_bins_pre_c2 = nan(n_bins,max(seq_nr),wb_max);
pitch_wb_L_seq_bins_pre_c2 = nan(n_bins,max(seq_nr),wb_max);

stroke_wb_R_seq_bins_pre_c2 = nan(n_bins,max(seq_nr),wb_max);
dev_wb_R_seq_bins_pre_c2 = nan(n_bins,max(seq_nr),wb_max);
pitch_wb_R_seq_bins_pre_c2 = nan(n_bins,max(seq_nr),wb_max);

stroke_wb_L_seq_bins_post_c2 = nan(n_bins,max(seq_nr),wb_max);
dev_wb_L_seq_bins_post_c2 = nan(n_bins,max(seq_nr),wb_max);
pitch_wb_L_seq_bins_post_c2 = nan(n_bins,max(seq_nr),wb_max);

stroke_wb_R_seq_bins_post_c2 = nan(n_bins,max(seq_nr),wb_max);
dev_wb_R_seq_bins_post_c2 = nan(n_bins,max(seq_nr),wb_max);
pitch_wb_R_seq_bins_post_c2 = nan(n_bins,max(seq_nr),wb_max);

%% body kin vars
% pre
roll_mean_wb_seq_pre_c2 = nan(max(seq_nr),wb_max);
roll_dot_mean_wb_seq_pre_c2 = nan(max(seq_nr),wb_max);
roll_dot_dot_mean_wb_seq_pre_c2 = nan(max(seq_nr),wb_max);

pitch_mean_wb_seq_pre_c2 = nan(max(seq_nr),wb_max);
pitch_dot_mean_wb_seq_pre_c2 = nan(max(seq_nr),wb_max);
pitch_dot_dot_mean_wb_seq_pre_c2 = nan(max(seq_nr),wb_max);

yaw_mean_wb_seq_pre_c2 = nan(max(seq_nr),wb_max);
yaw_dot_mean_wb_seq_pre_c2 = nan(max(seq_nr),wb_max);
yaw_dot_dot_mean_wb_seq_pre_c2 = nan(max(seq_nr),wb_max);

% drot_L_mean_wb_seq_pre_c2 = nan(max(seq_nr),wb_max);
% rot_L_mean_wb_seq_pre_c2 = nan(max(seq_nr),wb_max);
% rot_dot_L_mean_wb_seq_pre_c2 = nan(max(seq_nr),wb_max);
% rot_dot_dot_L_mean_wb_seq_pre_c2 = nan(max(seq_nr),wb_max);
% 
% drot_R_mean_wb_seq_pre_c2 = nan(max(seq_nr),wb_max);
% rot_R_mean_wb_seq_pre_c2 = nan(max(seq_nr),wb_max);
% rot_dot_R_mean_wb_seq_pre_c2 = nan(max(seq_nr),wb_max);
% rot_dot_dot_R_mean_wb_seq_pre_c2 = nan(max(seq_nr),wb_max);

% post
roll_mean_wb_seq_post_c2 = nan(max(seq_nr),wb_max);
roll_dot_mean_wb_seq_post_c2 = nan(max(seq_nr),wb_max);
roll_dot_dot_mean_wb_seq_post_c2 = nan(max(seq_nr),wb_max);

pitch_mean_wb_seq_post_c2 = nan(max(seq_nr),wb_max);
pitch_dot_mean_wb_seq_post_c2 = nan(max(seq_nr),wb_max);
pitch_dot_dot_mean_wb_seq_post_c2 = nan(max(seq_nr),wb_max);

yaw_mean_wb_seq_post_c2 = nan(max(seq_nr),wb_max);
yaw_dot_mean_wb_seq_post_c2 = nan(max(seq_nr),wb_max);
yaw_dot_dot_mean_wb_seq_post_c2 = nan(max(seq_nr),wb_max);

% drot_L_mean_wb_seq_post_c2 = nan(max(seq_nr),wb_max);
% rot_L_mean_wb_seq_post_c2 = nan(max(seq_nr),wb_max);
% rot_dot_L_mean_wb_seq_post_c2 = nan(max(seq_nr),wb_max);
% rot_dot_dot_L_mean_wb_seq_post_c2 = nan(max(seq_nr),wb_max);
% 
% drot_R_mean_wb_seq_post_c2 = nan(max(seq_nr),wb_max);
% rot_R_mean_wb_seq_post_c2 = nan(max(seq_nr),wb_max);
% rot_dot_R_mean_wb_seq_post_c2 = nan(max(seq_nr),wb_max);
% rot_dot_dot_R_mean_wb_seq_post_c2 = nan(max(seq_nr),wb_max);

%% Torque vars
% % pre
% Mroll_mean_wb_accel_seq_pre_c2 = nan(max(seq_nr),wb_max);
% Mroll_mean_wb_damp_seq_pre_c2 = nan(max(seq_nr),wb_max);
% Mroll_mean_wb_seq_pre_c2 = nan(max(seq_nr),wb_max);
% 
% Mpitch_mean_wb_accel_seq_pre_c2 = nan(max(seq_nr),wb_max);
% Mpitch_mean_wb_damp_seq_pre_c2 = nan(max(seq_nr),wb_max);
% Mpitch_mean_wb_seq_pre_c2 = nan(max(seq_nr),wb_max);
% 
% Myaw_mean_wb_accel_seq_pre_c2 = nan(max(seq_nr),wb_max);
% Myaw_mean_wb_damp_seq_pre_c2 = nan(max(seq_nr),wb_max);
% Myaw_mean_wb_seq_pre_c2 = nan(max(seq_nr),wb_max);
% 
% M_L_mean_wb_accel_seq_pre_c2 = nan(max(seq_nr),wb_max);
% M_L_mean_wb_damp_seq_pre_c2 = nan(max(seq_nr),wb_max);
% M_L_mean_wb_seq_pre_c2 = nan(max(seq_nr),wb_max);
% 
% M_R_mean_wb_accel_seq_pre_c2 = nan(max(seq_nr),wb_max);
% M_R_mean_wb_damp_seq_pre_c2 = nan(max(seq_nr),wb_max);
% M_R_mean_wb_seq_pre_c2 = nan(max(seq_nr),wb_max);
% 
% % post
% Mroll_mean_wb_accel_seq_post_c2 = nan(max(seq_nr),wb_max);
% Mroll_mean_wb_damp_seq_post_c2 = nan(max(seq_nr),wb_max);
% Mroll_mean_wb_seq_post_c2 = nan(max(seq_nr),wb_max);
% 
% Mpitch_mean_wb_accel_seq_post_c2 = nan(max(seq_nr),wb_max);
% Mpitch_mean_wb_damp_seq_post_c2 = nan(max(seq_nr),wb_max);
% Mpitch_mean_wb_seq_post_c2 = nan(max(seq_nr),wb_max);
% 
% Myaw_mean_wb_accel_seq_post_c2 = nan(max(seq_nr),wb_max);
% Myaw_mean_wb_damp_seq_post_c2 = nan(max(seq_nr),wb_max);
% Myaw_mean_wb_seq_post_c2 = nan(max(seq_nr),wb_max);
% 
% M_L_mean_wb_accel_seq_post_c2 = nan(max(seq_nr),wb_max);
% M_L_mean_wb_damp_seq_post_c2 = nan(max(seq_nr),wb_max);
% M_L_mean_wb_seq_post_c2 = nan(max(seq_nr),wb_max);
% 
% M_R_mean_wb_accel_seq_post_c2 = nan(max(seq_nr),wb_max);
% M_R_mean_wb_damp_seq_post_c2 = nan(max(seq_nr),wb_max);
% M_R_mean_wb_seq_post_c2 = nan(max(seq_nr),wb_max);

%% CLUSTER 3

%% wb kin vars
f_wb_seq_pre_c3 = nan(max(seq_nr),wb_max);
f_wb_seq_post_c3 = nan(max(seq_nr),wb_max);

stroke_wb_L_seq_bins_pre_c3 = nan(n_bins,max(seq_nr),wb_max);
dev_wb_L_seq_bins_pre_c3 = nan(n_bins,max(seq_nr),wb_max);
pitch_wb_L_seq_bins_pre_c3 = nan(n_bins,max(seq_nr),wb_max);

stroke_wb_R_seq_bins_pre_c3 = nan(n_bins,max(seq_nr),wb_max);
dev_wb_R_seq_bins_pre_c3 = nan(n_bins,max(seq_nr),wb_max);
pitch_wb_R_seq_bins_pre_c3 = nan(n_bins,max(seq_nr),wb_max);

stroke_wb_L_seq_bins_post_c3 = nan(n_bins,max(seq_nr),wb_max);
dev_wb_L_seq_bins_post_c3 = nan(n_bins,max(seq_nr),wb_max);
pitch_wb_L_seq_bins_post_c3 = nan(n_bins,max(seq_nr),wb_max);

stroke_wb_R_seq_bins_post_c3 = nan(n_bins,max(seq_nr),wb_max);
dev_wb_R_seq_bins_post_c3 = nan(n_bins,max(seq_nr),wb_max);
pitch_wb_R_seq_bins_post_c3 = nan(n_bins,max(seq_nr),wb_max);

%% body kin vars
% pre
roll_mean_wb_seq_pre_c3 = nan(max(seq_nr),wb_max);
roll_dot_mean_wb_seq_pre_c3 = nan(max(seq_nr),wb_max);
roll_dot_dot_mean_wb_seq_pre_c3 = nan(max(seq_nr),wb_max);

pitch_mean_wb_seq_pre_c3 = nan(max(seq_nr),wb_max);
pitch_dot_mean_wb_seq_pre_c3 = nan(max(seq_nr),wb_max);
pitch_dot_dot_mean_wb_seq_pre_c3 = nan(max(seq_nr),wb_max);

yaw_mean_wb_seq_pre_c3 = nan(max(seq_nr),wb_max);
yaw_dot_mean_wb_seq_pre_c3 = nan(max(seq_nr),wb_max);
yaw_dot_dot_mean_wb_seq_pre_c3 = nan(max(seq_nr),wb_max);

% drot_L_mean_wb_seq_pre_c3 = nan(max(seq_nr),wb_max);
% rot_L_mean_wb_seq_pre_c3 = nan(max(seq_nr),wb_max);
% rot_dot_L_mean_wb_seq_pre_c3 = nan(max(seq_nr),wb_max);
% rot_dot_dot_L_mean_wb_seq_pre_c3 = nan(max(seq_nr),wb_max);
% 
% drot_R_mean_wb_seq_pre_c3 = nan(max(seq_nr),wb_max);
% rot_R_mean_wb_seq_pre_c3 = nan(max(seq_nr),wb_max);
% rot_dot_R_mean_wb_seq_pre_c3 = nan(max(seq_nr),wb_max);
% rot_dot_dot_R_mean_wb_seq_pre_c3 = nan(max(seq_nr),wb_max);

% post
roll_mean_wb_seq_post_c3 = nan(max(seq_nr),wb_max);
roll_dot_mean_wb_seq_post_c3 = nan(max(seq_nr),wb_max);
roll_dot_dot_mean_wb_seq_post_c3 = nan(max(seq_nr),wb_max);

pitch_mean_wb_seq_post_c3 = nan(max(seq_nr),wb_max);
pitch_dot_mean_wb_seq_post_c3 = nan(max(seq_nr),wb_max);
pitch_dot_dot_mean_wb_seq_post_c3 = nan(max(seq_nr),wb_max);

yaw_mean_wb_seq_post_c3 = nan(max(seq_nr),wb_max);
yaw_dot_mean_wb_seq_post_c3 = nan(max(seq_nr),wb_max);
yaw_dot_dot_mean_wb_seq_post_c3 = nan(max(seq_nr),wb_max);

% drot_L_mean_wb_seq_post_c3 = nan(max(seq_nr),wb_max);
% rot_L_mean_wb_seq_post_c3 = nan(max(seq_nr),wb_max);
% rot_dot_L_mean_wb_seq_post_c3 = nan(max(seq_nr),wb_max);
% rot_dot_dot_L_mean_wb_seq_post_c3 = nan(max(seq_nr),wb_max);
% 
% drot_R_mean_wb_seq_post_c3 = nan(max(seq_nr),wb_max);
% rot_R_mean_wb_seq_post_c3 = nan(max(seq_nr),wb_max);
% rot_dot_R_mean_wb_seq_post_c3 = nan(max(seq_nr),wb_max);
% rot_dot_dot_R_mean_wb_seq_post_c3 = nan(max(seq_nr),wb_max);

%% Torque vars
% % pre
% Mroll_mean_wb_accel_seq_pre_c3 = nan(max(seq_nr),wb_max);
% Mroll_mean_wb_damp_seq_pre_c3 = nan(max(seq_nr),wb_max);
% Mroll_mean_wb_seq_pre_c3 = nan(max(seq_nr),wb_max);
% 
% Mpitch_mean_wb_accel_seq_pre_c3 = nan(max(seq_nr),wb_max);
% Mpitch_mean_wb_damp_seq_pre_c3 = nan(max(seq_nr),wb_max);
% Mpitch_mean_wb_seq_pre_c3 = nan(max(seq_nr),wb_max);
% 
% Myaw_mean_wb_accel_seq_pre_c3 = nan(max(seq_nr),wb_max);
% Myaw_mean_wb_damp_seq_pre_c3 = nan(max(seq_nr),wb_max);
% Myaw_mean_wb_seq_pre_c3 = nan(max(seq_nr),wb_max);
% 
% M_L_mean_wb_accel_seq_pre_c3 = nan(max(seq_nr),wb_max);
% M_L_mean_wb_damp_seq_pre_c3 = nan(max(seq_nr),wb_max);
% M_L_mean_wb_seq_pre_c3 = nan(max(seq_nr),wb_max);
% 
% M_R_mean_wb_accel_seq_pre_c3 = nan(max(seq_nr),wb_max);
% M_R_mean_wb_damp_seq_pre_c3 = nan(max(seq_nr),wb_max);
% M_R_mean_wb_seq_pre_c3 = nan(max(seq_nr),wb_max);
% 
% % post
% Mroll_mean_wb_accel_seq_post_c3 = nan(max(seq_nr),wb_max);
% Mroll_mean_wb_damp_seq_post_c3 = nan(max(seq_nr),wb_max);
% Mroll_mean_wb_seq_post_c3 = nan(max(seq_nr),wb_max);
% 
% Mpitch_mean_wb_accel_seq_post_c3 = nan(max(seq_nr),wb_max);
% Mpitch_mean_wb_damp_seq_post_c3 = nan(max(seq_nr),wb_max);
% Mpitch_mean_wb_seq_post_c3 = nan(max(seq_nr),wb_max);
% 
% Myaw_mean_wb_accel_seq_post_c3 = nan(max(seq_nr),wb_max);
% Myaw_mean_wb_damp_seq_post_c3 = nan(max(seq_nr),wb_max);
% Myaw_mean_wb_seq_post_c3 = nan(max(seq_nr),wb_max);
% 
% M_L_mean_wb_accel_seq_post_c3 = nan(max(seq_nr),wb_max);
% M_L_mean_wb_damp_seq_post_c3 = nan(max(seq_nr),wb_max);
% M_L_mean_wb_seq_post_c3 = nan(max(seq_nr),wb_max);
% 
% M_R_mean_wb_accel_seq_post_c3 = nan(max(seq_nr),wb_max);
% M_R_mean_wb_damp_seq_post_c3 = nan(max(seq_nr),wb_max);
% M_R_mean_wb_seq_post_c3 = nan(max(seq_nr),wb_max);

%% CLUSTER 4

%% wb kin vars
f_wb_seq_pre_c4 = nan(max(seq_nr),wb_max);
f_wb_seq_post_c4 = nan(max(seq_nr),wb_max);

stroke_wb_L_seq_bins_pre_c4 = nan(n_bins,max(seq_nr),wb_max);
dev_wb_L_seq_bins_pre_c4 = nan(n_bins,max(seq_nr),wb_max);
pitch_wb_L_seq_bins_pre_c4 = nan(n_bins,max(seq_nr),wb_max);

stroke_wb_R_seq_bins_pre_c4 = nan(n_bins,max(seq_nr),wb_max);
dev_wb_R_seq_bins_pre_c4 = nan(n_bins,max(seq_nr),wb_max);
pitch_wb_R_seq_bins_pre_c4 = nan(n_bins,max(seq_nr),wb_max);

stroke_wb_L_seq_bins_post_c4 = nan(n_bins,max(seq_nr),wb_max);
dev_wb_L_seq_bins_post_c4 = nan(n_bins,max(seq_nr),wb_max);
pitch_wb_L_seq_bins_post_c4 = nan(n_bins,max(seq_nr),wb_max);

stroke_wb_R_seq_bins_post_c4 = nan(n_bins,max(seq_nr),wb_max);
dev_wb_R_seq_bins_post_c4 = nan(n_bins,max(seq_nr),wb_max);
pitch_wb_R_seq_bins_post_c4 = nan(n_bins,max(seq_nr),wb_max);

%% body kin vars
% pre
roll_mean_wb_seq_pre_c4 = nan(max(seq_nr),wb_max);
roll_dot_mean_wb_seq_pre_c4 = nan(max(seq_nr),wb_max);
roll_dot_dot_mean_wb_seq_pre_c4 = nan(max(seq_nr),wb_max);

pitch_mean_wb_seq_pre_c4 = nan(max(seq_nr),wb_max);
pitch_dot_mean_wb_seq_pre_c4 = nan(max(seq_nr),wb_max);
pitch_dot_dot_mean_wb_seq_pre_c4 = nan(max(seq_nr),wb_max);

yaw_mean_wb_seq_pre_c4 = nan(max(seq_nr),wb_max);
yaw_dot_mean_wb_seq_pre_c4 = nan(max(seq_nr),wb_max);
yaw_dot_dot_mean_wb_seq_pre_c4 = nan(max(seq_nr),wb_max);

% drot_L_mean_wb_seq_pre_c4 = nan(max(seq_nr),wb_max);
% rot_L_mean_wb_seq_pre_c4 = nan(max(seq_nr),wb_max);
% rot_dot_L_mean_wb_seq_pre_c4 = nan(max(seq_nr),wb_max);
% rot_dot_dot_L_mean_wb_seq_pre_c4 = nan(max(seq_nr),wb_max);
% 
% drot_R_mean_wb_seq_pre_c4 = nan(max(seq_nr),wb_max);
% rot_R_mean_wb_seq_pre_c4 = nan(max(seq_nr),wb_max);
% rot_dot_R_mean_wb_seq_pre_c4 = nan(max(seq_nr),wb_max);
% rot_dot_dot_R_mean_wb_seq_pre_c4 = nan(max(seq_nr),wb_max);

% post
roll_mean_wb_seq_post_c4 = nan(max(seq_nr),wb_max);
roll_dot_mean_wb_seq_post_c4 = nan(max(seq_nr),wb_max);
roll_dot_dot_mean_wb_seq_post_c4 = nan(max(seq_nr),wb_max);

pitch_mean_wb_seq_post_c4 = nan(max(seq_nr),wb_max);
pitch_dot_mean_wb_seq_post_c4 = nan(max(seq_nr),wb_max);
pitch_dot_dot_mean_wb_seq_post_c4 = nan(max(seq_nr),wb_max);

yaw_mean_wb_seq_post_c4 = nan(max(seq_nr),wb_max);
yaw_dot_mean_wb_seq_post_c4 = nan(max(seq_nr),wb_max);
yaw_dot_dot_mean_wb_seq_post_c4 = nan(max(seq_nr),wb_max);

% drot_L_mean_wb_seq_post_c4 = nan(max(seq_nr),wb_max);
% rot_L_mean_wb_seq_post_c4 = nan(max(seq_nr),wb_max);
% rot_dot_L_mean_wb_seq_post_c4 = nan(max(seq_nr),wb_max);
% rot_dot_dot_L_mean_wb_seq_post_c4 = nan(max(seq_nr),wb_max);
% 
% drot_R_mean_wb_seq_post_c4 = nan(max(seq_nr),wb_max);
% rot_R_mean_wb_seq_post_c4 = nan(max(seq_nr),wb_max);
% rot_dot_R_mean_wb_seq_post_c4 = nan(max(seq_nr),wb_max);
% rot_dot_dot_R_mean_wb_seq_post_c4 = nan(max(seq_nr),wb_max);

%% Torque vars
% % pre
% Mroll_mean_wb_accel_seq_pre_c4 = nan(max(seq_nr),wb_max);
% Mroll_mean_wb_damp_seq_pre_c4 = nan(max(seq_nr),wb_max);
% Mroll_mean_wb_seq_pre_c4 = nan(max(seq_nr),wb_max);
% 
% Mpitch_mean_wb_accel_seq_pre_c4 = nan(max(seq_nr),wb_max);
% Mpitch_mean_wb_damp_seq_pre_c4 = nan(max(seq_nr),wb_max);
% Mpitch_mean_wb_seq_pre_c4 = nan(max(seq_nr),wb_max);
% 
% Myaw_mean_wb_accel_seq_pre_c4 = nan(max(seq_nr),wb_max);
% Myaw_mean_wb_damp_seq_pre_c4 = nan(max(seq_nr),wb_max);
% Myaw_mean_wb_seq_pre_c4 = nan(max(seq_nr),wb_max);
% 
% M_L_mean_wb_accel_seq_pre_c4 = nan(max(seq_nr),wb_max);
% M_L_mean_wb_damp_seq_pre_c4 = nan(max(seq_nr),wb_max);
% M_L_mean_wb_seq_pre_c4 = nan(max(seq_nr),wb_max);
% 
% M_R_mean_wb_accel_seq_pre_c4 = nan(max(seq_nr),wb_max);
% M_R_mean_wb_damp_seq_pre_c4 = nan(max(seq_nr),wb_max);
% M_R_mean_wb_seq_pre_c4 = nan(max(seq_nr),wb_max);
% 
% % post
% Mroll_mean_wb_accel_seq_post_c4 = nan(max(seq_nr),wb_max);
% Mroll_mean_wb_damp_seq_post_c4 = nan(max(seq_nr),wb_max);
% Mroll_mean_wb_seq_post_c4 = nan(max(seq_nr),wb_max);
% 
% Mpitch_mean_wb_accel_seq_post_c4 = nan(max(seq_nr),wb_max);
% Mpitch_mean_wb_damp_seq_post_c4 = nan(max(seq_nr),wb_max);
% Mpitch_mean_wb_seq_post_c4 = nan(max(seq_nr),wb_max);
% 
% Myaw_mean_wb_accel_seq_post_c4 = nan(max(seq_nr),wb_max);
% Myaw_mean_wb_damp_seq_post_c4 = nan(max(seq_nr),wb_max);
% Myaw_mean_wb_seq_post_c4 = nan(max(seq_nr),wb_max);
% 
% M_L_mean_wb_accel_seq_post_c4 = nan(max(seq_nr),wb_max);
% M_L_mean_wb_damp_seq_post_c4 = nan(max(seq_nr),wb_max);
% M_L_mean_wb_seq_post_c4 = nan(max(seq_nr),wb_max);
% 
% M_R_mean_wb_accel_seq_post_c4 = nan(max(seq_nr),wb_max);
% M_R_mean_wb_damp_seq_post_c4 = nan(max(seq_nr),wb_max);
% M_R_mean_wb_seq_post_c4 = nan(max(seq_nr),wb_max);

