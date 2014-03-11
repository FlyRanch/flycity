
for i = 1:length(t_shift)
    t_plot(1:length(t),i) = t-t_shift(i);
end

% subset_cut = [36 72 108 144];
% subset_mid = [18 54 90 126 162];

n1 = find(abs(turn_angle_vel_mirror)<subset_cut(1));
n2 = find(abs(turn_angle_vel_mirror)>subset_cut(1));
% n3 = find(abs(turn_angle_vel_mirror)<subset_cut(3)&abs(turn_angle_vel_mirror)>subset_cut(2));
% n4 = find(abs(turn_angle_vel_mirror)>subset_cut(3));
% n4 = find(abs(turn_angle_vel_mirror)<subset_cut(4)&abs(turn_angle_vel_mirror)>subset_cut(3));
% n5 = find(abs(turn_angle_vel_mirror)>subset_cut(4));

% figure
% plot(abs(turn_angle_vel_mirror(n1)),'.')
% hold on
% plot(abs(turn_angle_vel_mirror(n2)),'r.')
% plot(abs(turn_angle_vel_mirror(n3)),'.g')
% plot(abs(turn_angle_vel_mirror(n4)),'.c')
% plot(abs(turn_angle_vel_mirror(n5)),'.m')

t_c1 = t_plot(:,n1);
t_c2 = t_plot(:,n2);
% t_c3 = t_plot(:,n3);
% t_c4 = t_plot(:,n4);
% t_c5 = t_plot(:,n5);

stim_angle_vel_c1 = stim_angle_vel_plot(:,n1);
stim_angle_vel_c2 = stim_angle_vel_plot(:,n2);
% stim_angle_vel_c3 = stim_angle_vel_plot(:,n3);
% stim_angle_vel_c4 = stim_angle_vel_plot(:,n4);
% stim_angle_vel_c5 = stim_angle_vel_plot(:,n5);

dV_c1 = dV_plot(:,n1);
dV_c2 = dV_plot(:,n2);
% dV_c3 = dV_plot(:,n3);
% dV_c4 = dV_plot(:,n4);
% dV_c5 = dV_plot(:,n5);

stim_angle_accel_c1 = stim_angle_accel_plot(:,n1);
stim_angle_accel_c2 = stim_angle_accel_plot(:,n2);
% stim_angle_accel_c3 = stim_angle_accel_plot(:,n3);
% stim_angle_accel_c4 = stim_angle_accel_plot(:,n4);
% stim_angle_accel_c5 = stim_angle_accel_plot(:,n5);

F_c1 = F_plot(:,n1);
F_c2 = F_plot(:,n2);
% F_c3 = F_plot(:,n3);
% F_c4 = F_plot(:,n4);
% F_c5 = F_plot(:,n5);

A_norm_c1 = A_norm_plot(:,n1);
A_norm_c2 = A_norm_plot(:,n2);

Ahor_norm_c1 = Ahor_norm_plot(:,n1);
Ahor_norm_c2 = Ahor_norm_plot(:,n2);

Az_norm_c1 = Az_norm_plot(:,n1);
Az_norm_c2 = Az_norm_plot(:,n2);

Fsp_roll_c1 = Fsp_roll_plot(:,n1);
Fsp_roll_c2 = Fsp_roll_plot(:,n2);
% Fsp_roll_c3 = Fsp_roll_plot(:,n3);
% Fsp_roll_c4 = Fsp_roll_plot(:,n4);
% Fsp_roll_c5 = Fsp_roll_plot(:,n5);

Fsp_pitch_c1 = Fsp_pitch_plot(:,n1);
Fsp_pitch_c2 = Fsp_pitch_plot(:,n2);
% Fsp_pitch_c3 = Fsp_pitch_plot(:,n3);
% Fsp_pitch_c4 = Fsp_pitch_plot(:,n4);
% Fsp_pitch_c5 = Fsp_pitch_plot(:,n5);

droll_c1 = droll_plot(:,n1);
droll_c2 = droll_plot(:,n2);
% droll_c3 = droll_plot(:,n3);
% droll_c4 = droll_plot(:,n4);
% droll_c5 = droll_plot(:,n5);

roll_c1 = roll_plot(:,n1);
roll_c2 = roll_plot(:,n2);
% roll_c3 = roll_plot(:,n3);
% roll_c4 = roll_plot(:,n4);
% roll_c5 = roll_plot(:,n5);

roll_dot_c1 = roll_dot_plot(:,n1);
roll_dot_c2 = roll_dot_plot(:,n2);
% roll_dot_c3 = roll_dot_plot(:,n3);
% roll_dot_c4 = roll_dot_plot(:,n4);
% roll_dot_c5 = roll_dot_plot(:,n5);

roll_dot_dot_c1 = roll_dot_dot_plot(:,n1);
roll_dot_dot_c2 = roll_dot_dot_plot(:,n2);
% roll_dot_dot_c3 = roll_dot_dot_plot(:,n3);
% roll_dot_dot_c4 = roll_dot_dot_plot(:,n4);
% roll_dot_dot_c5 = roll_dot_dot_plot(:,n5);

dpitch_c1 = dpitch_plot(:,n1);
dpitch_c2 = dpitch_plot(:,n2);
% dpitch_c3 = dpitch_plot(:,n3);
% dpitch_c4 = dpitch_plot(:,n4);
% dpitch_c5 = dpitch_plot(:,n5);

pitch_c1 = pitch_plot(:,n1);
pitch_c2 = pitch_plot(:,n2);
% pitch_c3 = pitch_plot(:,n3);
% pitch_c4 = pitch_plot(:,n4);
% pitch_c5 = pitch_plot(:,n5);

pitch_dot_c1 = pitch_dot_plot(:,n1);
pitch_dot_c2 = pitch_dot_plot(:,n2);
% pitch_dot_c3 = pitch_dot_plot(:,n3);
% pitch_dot_c4 = pitch_dot_plot(:,n4);
% pitch_dot_c5 = pitch_dot_plot(:,n5);

pitch_dot_dot_c1 = pitch_dot_dot_plot(:,n1);
pitch_dot_dot_c2 = pitch_dot_dot_plot(:,n2);
% pitch_dot_dot_c3 = pitch_dot_dot_plot(:,n3);
% pitch_dot_dot_c4 = pitch_dot_dot_plot(:,n4);
% pitch_dot_dot_c5 = pitch_dot_dot_plot(:,n5);

yaw_c1 = yaw_plot(:,n1);
yaw_c2 = yaw_plot(:,n2);
% yaw_c3 = yaw_plot(:,n3);
% yaw_c4 = yaw_plot(:,n4);
% yaw_c5 = yaw_plot(:,n5);

yaw_dot_c1 = yaw_dot_plot(:,n1);
yaw_dot_c2 = yaw_dot_plot(:,n2);
% yaw_dot_c3 = yaw_dot_plot(:,n3);
% yaw_dot_c4 = yaw_dot_plot(:,n4);
% yaw_dot_c5 = yaw_dot_plot(:,n5);

yaw_dot_dot_c1 = yaw_dot_dot_plot(:,n1);
yaw_dot_dot_c2 = yaw_dot_dot_plot(:,n2);
% yaw_dot_dot_c3 = yaw_dot_dot_plot(:,n3);
% yaw_dot_dot_c4 = yaw_dot_dot_plot(:,n4);
% yaw_dot_dot_c5 = yaw_dot_dot_plot(:,n5);

% rot axes
rot_axis_pos_c1 = rot_axis_pos_plot(:,n1);
rot_axis_pos_c2 = rot_axis_pos_plot(:,n2);

rot_axis_vel_c1 = rot_axis_vel_plot(:,n1);
rot_axis_vel_c2 = rot_axis_vel_plot(:,n2);

rot_axis_accel_c1 = rot_axis_accel_plot(:,n1);
rot_axis_accel_c2 = rot_axis_accel_plot(:,n2);

drot_L_c1 = drot_L_plot(:,n1);
drot_L_c2 = drot_L_plot(:,n2);
drot_R_c1 = drot_R_plot(:,n1);
drot_R_c2 = drot_R_plot(:,n2);

rot_dot_L_c1 = rot_dot_L_plot(:,n1);
rot_dot_L_c2 = rot_dot_L_plot(:,n2);
rot_dot_R_c1 = rot_dot_R_plot(:,n1);
rot_dot_R_c2 = rot_dot_R_plot(:,n2);

rot_dot_dot_L_c1 = rot_dot_dot_L_plot(:,n1);
rot_dot_dot_L_c2 = rot_dot_dot_L_plot(:,n2);
rot_dot_dot_R_c1 = rot_dot_dot_R_plot(:,n1);
rot_dot_dot_R_c2 = rot_dot_dot_R_plot(:,n2);

% torque
Mroll_accel_c1 = Mroll_accel(:,n1);
Mroll_accel_c2 = Mroll_accel(:,n2);

Mpitch_accel_c1 = Mpitch_accel(:,n1);
Mpitch_accel_c2 = Mpitch_accel(:,n2);

Myaw_accel_c1 = Myaw_accel(:,n1);
Myaw_accel_c2 = Myaw_accel(:,n2);

M_L_accel_c1 = M_L_accel(:,n1);
M_L_accel_c2 = M_L_accel(:,n2);

M_R_accel_c1 = M_R_accel(:,n1);
M_R_accel_c2 = M_R_accel(:,n2);

Mroll_damp_c1 = Mroll_damp(:,n1);
Mroll_damp_c2 = Mroll_damp(:,n2);

Mpitch_damp_c1 = Mpitch_damp(:,n1);
Mpitch_damp_c2 = Mpitch_damp(:,n2);

Myaw_damp_c1 = Myaw_damp(:,n1);
Myaw_damp_c2 = Myaw_damp(:,n2);

M_L_damp_c1 = M_L_damp(:,n1);
M_L_damp_c2 = M_L_damp(:,n2);

M_R_damp_c1 = M_R_damp(:,n1);
M_R_damp_c2 = M_R_damp(:,n2);

Mroll_c1 = Mroll(:,n1);
Mroll_c2 = Mroll(:,n2);

Mpitch_c1 = Mpitch(:,n1);
Mpitch_c2 = Mpitch(:,n2);

Myaw_c1 = Myaw(:,n1);
Myaw_c2 = Myaw(:,n2);

M_L_c1 = M_L(:,n1);
M_L_c2 = M_L(:,n2);

M_R_c1 = M_R(:,n1);
M_R_c2 = M_R(:,n2);

% NONdim torque
Mroll_accel_norm_c1 = Mroll_accel_norm(:,n1);
Mroll_accel_norm_c2 = Mroll_accel_norm(:,n2);

Mpitch_accel_norm_c1 = Mpitch_accel_norm(:,n1);
Mpitch_accel_norm_c2 = Mpitch_accel_norm(:,n2);

Myaw_accel_norm_c1 = Myaw_accel_norm(:,n1);
Myaw_accel_norm_c2 = Myaw_accel_norm(:,n2);

M_L_accel_norm_c1 = M_L_accel_norm(:,n1);
M_L_accel_norm_c2 = M_L_accel_norm(:,n2);

M_R_accel_norm_c1 = M_R_accel_norm(:,n1);
M_R_accel_norm_c2 = M_R_accel_norm(:,n2);

Mroll_damp_norm_c1 = Mroll_damp_norm(:,n1);
Mroll_damp_norm_c2 = Mroll_damp_norm(:,n2);

Mpitch_damp_norm_c1 = Mpitch_damp_norm(:,n1);
Mpitch_damp_norm_c2 = Mpitch_damp_norm(:,n2);

Myaw_damp_norm_c1 = Myaw_damp_norm(:,n1);
Myaw_damp_norm_c2 = Myaw_damp_norm(:,n2);

M_L_damp_norm_c1 = M_L_damp_norm(:,n1);
M_L_damp_norm_c2 = M_L_damp_norm(:,n2);

M_R_damp_norm_c1 = M_R_damp_norm(:,n1);
M_R_damp_norm_c2 = M_R_damp_norm(:,n2);

Mroll_norm_c1 = Mroll_norm(:,n1);
Mroll_norm_c2 = Mroll_norm(:,n2);

Mpitch_norm_c1 = Mpitch_norm(:,n1);
Mpitch_norm_c2 = Mpitch_norm(:,n2);

Myaw_norm_c1 = Myaw_norm(:,n1);
Myaw_norm_c2 = Myaw_norm(:,n2);

M_L_norm_c1 = M_L_norm(:,n1);
M_L_norm_c2 = M_L_norm(:,n2);

M_R_norm_c1 = M_R_norm(:,n1);
M_R_norm_c2 = M_R_norm(:,n2);

% torque axes
Maccel_axis_c1 = Maccel_axis_plot(:,n1);
Maccel_axis_c2 = Maccel_axis_plot(:,n2);

Mdamp_axis_c1 = Mdamp_axis_plot(:,n1);
Mdamp_axis_c2 = Mdamp_axis_plot(:,n2);

M_axis_c1 = M_axis_plot(:,n1);
M_axis_c2 = M_axis_plot(:,n2);

% % make means
% angle_pre_min = -.02;
% angle_pre_max = .04;
% mean_min = -.02;
% mean_max = .04;
% 
% figure
% hold on
% plot(t_c1(:),stim_angle_vel_c1(:),'.')
% plot(t_c2(:),stim_angle_vel_c2(:),'.r')
% plot(t_c3(:),stim_angle_vel_c3(:),'.g')
% plot(t_c4(:),stim_angle_vel_c4(:),'.c')
% plot(t_c5(:),stim_angle_vel_c5(:),'.m')
% axis([angle_pre_min angle_pre_max -180 180])
% 
% angle_pre = t_c1(:);
% angle_post = stim_angle_vel_c1(:);
% plotcolor = cmap_plot(subset_mid(1),:);
% plot_angle_pre_post_circmean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c2(:);
% angle_post = stim_angle_vel_c2(:);
% plotcolor = cmap_plot(subset_mid(2),:);
% plot_angle_pre_post_circmean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c3(:);
% angle_post = stim_angle_vel_c3(:);
% plotcolor = cmap_plot(subset_mid(3),:);
% plot_angle_pre_post_circmean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c4(:);
% angle_post = stim_angle_vel_c4(:);
% plotcolor = cmap_plot(subset_mid(4),:);
% plot_angle_pre_post_circmean_NOextendedsection_meanlim_NOdp
% 
% angle_pre = t_c5(:);
% angle_post = stim_angle_vel_c5(:);
% plotcolor = cmap_plot(subset_mid(5),:);
% plot_angle_pre_post_circmean_NOextendedsection_meanlim_NOdp
% 
% colorbar


