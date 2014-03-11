% heading
stim_angle_vel_pre = calc_value(stim_angle_vel,n_pre);
stim_angle_vel_post = calc_value(stim_angle_vel,n_post);
[stim_angle_vel_mean stim_angle_vel_ul stim_angle_vel_ll] = calc_circ_mean_value(stim_angle_vel,n_pre,n_post);

% body dir
stim_angle_yaw_pre = calc_value(stim_angle_yaw,n_pre);
stim_angle_yaw_post = calc_value(stim_angle_yaw,n_post);
[stim_angle_yaw_mean stim_angle_yaw_ul stim_angle_yaw_ll] = calc_circ_mean_value(stim_angle_yaw,n_pre,n_post);

% Adir
stim_angle_accel_pre = calc_value(stim_angle_accel,n_pre);
stim_angle_accel_post = calc_value(stim_angle_accel,n_post);
[stim_angle_accel_mean stim_angle_accel_ul stim_angle_accel_ll] = calc_circ_mean_value(stim_angle_accel,n_pre,n_post);

for i = 1:length(stim_angle_vel_pre)
    if isnan(stim_angle_vel_pre(i)) == 0 && isnan(stim_angle_vel_post(i)) == 0

        % 2d yaw vector pre
        x_yaw_pre = cosd(stim_angle_yaw_pre(i));
        y_yaw_pre = sind(stim_angle_yaw_pre(i));

        x_yaw_post = cosd(stim_angle_yaw_post(i));
        y_yaw_post = sind(stim_angle_yaw_post(i));

        % 2d velocity vector
        x_vel_pre = cosd(stim_angle_vel_pre(i));
        y_vel_pre = sind(stim_angle_vel_pre(i));

        x_vel_post = cosd(stim_angle_vel_post(i));
        y_vel_post = sind(stim_angle_vel_post(i));
        
        % 2d A vector
        x_A_pre = cosd(stim_angle_accel_pre(i));
        y_A_pre = sind(stim_angle_accel_pre(i));

        x_A_post = cosd(stim_angle_accel_post(i));
        y_A_post = sind(stim_angle_accel_post(i));

        x_A_mean = cosd(stim_angle_accel_mean(i));
        y_A_mean = sind(stim_angle_accel_mean(i));

        turn_angle_vel(i,1) =   atan2(x_vel_pre*y_vel_post-y_vel_pre*x_vel_post,x_vel_pre*x_vel_post+y_vel_pre*y_vel_post) *180/pi();
        turn_angle_yaw(i,1) =   atan2(x_yaw_pre*y_yaw_post-y_yaw_pre*x_yaw_post,x_yaw_pre*x_yaw_post+y_yaw_pre*y_yaw_post) *180/pi();
        turn_angle_Amean2Vpre(i,1) = atan2(x_vel_pre*y_A_mean-y_vel_pre*x_A_mean,x_vel_pre*x_A_mean+y_vel_pre*y_A_mean) *180/pi();
    else
        turn_angle_vel(i,1) =   nan;
        turn_angle_yaw(i,1) =   nan;
        turn_angle_Amean2Vpre(i,1) = nan;
    end
end

% mirror
turn_angle_vel_mirror = turn_angle_vel;
turn_angle_yaw_mirror = turn_angle_yaw;
turn_angle_Amean2Vpre_mirror = turn_angle_Amean2Vpre;
for i=1:length(An_hor_max_mirror)
    if An_hor_max_mirror(i) < 0
        turn_angle_vel_mirror(i) = -turn_angle_vel_mirror(i);
        turn_angle_yaw_mirror(i) = -turn_angle_yaw_mirror(i);
        turn_angle_Amean2Vpre_mirror(i) = -turn_angle_Amean2Vpre_mirror(i);
    end
end
