for i = 1:length(t_pre)
    if isnan(t_pre(i)) == 0 && isnan(t_post(i)) == 0
        n_pre = find(t==t_pre(i));
        n_post = find(t==t_post(i));

        % 2d yaw vector pre
        x_yaw_pre = cosd(pathDB.yaw_global(n_pre,i));
        y_yaw_pre = sind(pathDB.yaw_global(n_pre,i));

        x_yaw_post = cosd(pathDB.yaw_global(n_post,i));
        y_yaw_post = sind(pathDB.yaw_global(n_post,i));

        % 2d velocity vector
        x_vel_pre = pathDB.vel(n_pre,i,1);
        y_vel_pre = pathDB.vel(n_pre,i,2);

        x_vel_post = pathDB.vel(n_post,i,1);
        y_vel_post = pathDB.vel(n_post,i,2);
% 
%         turn_angle_vel(i,1) =   atan2(x_vel_pre*y_vel_post-y_vel_pre*x_vel_post,x_vel_pre*x_vel_post+y_vel_pre*y_vel_post) *180/pi();
%         turn_angle_yaw(i,1) =   atan2(x_yaw_pre*y_yaw_post-y_yaw_pre*x_yaw_post,x_yaw_pre*x_yaw_post+y_yaw_pre*y_yaw_post) *180/pi();
        turn_angle_vel(i,1) =   atan2(x_vel_post*y_vel_pre-y_vel_post*x_vel_pre,x_vel_post*x_vel_pre+y_vel_post*y_vel_pre) *180/pi();
        turn_angle_yaw(i,1) =   atan2(x_yaw_post*y_yaw_pre-y_yaw_post*x_yaw_pre,x_yaw_post*x_yaw_pre+y_yaw_post*y_yaw_pre) *180/pi();
    else
        turn_angle_vel(i,1) =   nan;
        turn_angle_yaw(i,1) =   nan;
    end
end
