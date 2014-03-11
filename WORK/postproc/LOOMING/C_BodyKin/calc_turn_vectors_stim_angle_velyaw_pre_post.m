for i = 1:length(stim_angle_vel_pre)
    if isnan(stim_angle_vel_pre(i)) == 0 && isnan(stim_angle_vel_post(i)) == 0

        % 2d yaw vector pre
        x_yaw_pre = cosd(yaw_pre(i));
        y_yaw_pre = sind(yaw_pre(i));

        x_yaw_post = cosd(yaw_post(i));
        y_yaw_post = sind(yaw_post(i));

        % 2d velocity vector
        x_vel_pre = cosd(stim_angle_vel_pre(i));
        y_vel_pre = sind(stim_angle_vel_pre(i));

        x_vel_post = cosd(stim_angle_vel_post(i));
        y_vel_post = sind(stim_angle_vel_post(i));

        turn_angle_vel(i,1) =   atan2(x_vel_pre*y_vel_post-y_vel_pre*x_vel_post,x_vel_pre*x_vel_post+y_vel_pre*y_vel_post) *180/pi();
        turn_angle_yaw(i,1) =   atan2(x_yaw_pre*y_yaw_post-y_yaw_pre*x_yaw_post,x_yaw_pre*x_yaw_post+y_yaw_pre*y_yaw_post) *180/pi();
%         turn_angle_vel(i,1) =   atan2(x_vel_post*y_vel_pre-y_vel_post*x_vel_pre,x_vel_post*x_vel_pre+y_vel_post*y_vel_pre) *180/pi();
%         turn_angle_yaw(i,1) =   atan2(x_yaw_post*y_yaw_pre-y_yaw_post*x_yaw_pre,x_yaw_post*x_yaw_pre+y_yaw_post*y_yaw_pre) *180/pi();
    else
        turn_angle_vel(i,1) =   nan;
        turn_angle_yaw(i,1) =   nan;
    end
end
