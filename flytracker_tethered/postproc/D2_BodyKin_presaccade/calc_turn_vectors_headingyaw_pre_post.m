for i = 1:length(heading_pre)
    if isnan(heading_pre(i)) == 0 && isnan(heading_post(i)) == 0

        % 2d yaw vector pre
        x_yaw_pre = cosd(stim_angle_yaw_pre(i));
        y_yaw_pre = sind(stim_angle_yaw_pre(i));

        x_yaw_post = cosd(stim_angle_yaw_post(i));
        y_yaw_post = sind(stim_angle_yaw_post(i));

        % 2d velocity vector
        x_vel_pre = cosd(heading_pre(i));
        y_vel_pre = sind(heading_pre(i));

        x_vel_post = cosd(heading_post(i));
        y_vel_post = sind(heading_post(i));
        
        % 2d A vector
        x_A_pre = cosd(Adir_pre(i));
        y_A_pre = sind(Adir_pre(i));

        x_A_post = cosd(Adir_post(i));
        y_A_post = sind(Adir_post(i));

        x_Adir_mean = cosd(Adir_mean(i));
        y_Adir_mean = sind(Adir_mean(i));
        

        turn_angle_vel(i,1) =   atan2(x_vel_pre*y_vel_post-y_vel_pre*x_vel_post,x_vel_pre*x_vel_post+y_vel_pre*y_vel_post) *180/pi();
        turn_angle_yaw(i,1) =   atan2(x_yaw_pre*y_yaw_post-y_yaw_pre*x_yaw_post,x_yaw_pre*x_yaw_post+y_yaw_pre*y_yaw_post) *180/pi();
        turn_angle_Amean2Vpre(i,1) = atan2(x_vel_pre*y_Adir_mean-y_vel_pre*x_Adir_mean,x_vel_pre*x_Adir_mean+y_vel_pre*y_Adir_mean) *180/pi();
    else
        turn_angle_vel(i,1) =   nan;
        turn_angle_yaw(i,1) =   nan;
        turn_angle_Amean2Vpre(i,1) = nan;
    end
end
