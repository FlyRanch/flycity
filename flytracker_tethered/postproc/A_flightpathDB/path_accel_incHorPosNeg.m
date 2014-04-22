function [At,An,alpha_dot,At_hor,An_hor,alpha_dot_hor,accel_angle_hor_vel] = path_accel_incHorPosNeg(vel,accel,t)
dt = (t(2)-t(1));
for j=1:size(vel,2)
    counter = size(vel,2) - j
    
    for i=1:size(vel,1)
        
        vel_now(:,1) = vel(i,j,:);
        accel_now(:,1) = accel(i,j,:);
        At(i,j) = dot(vel_now,accel_now) / norm(vel_now);
        An(i,j) = sqrt(norm(accel_now)^2 - At(i,j)^2);
        
        vel_hor_now(:,1) = vel(i,j,1:2);
        accel_hor_now(:,1) = accel(i,j,1:2);
        
        % 2d horizontal accel vector
        x_a = accel_hor_now(1);
        y_a = accel_hor_now(2);
        % 2d horizontal velocity vector
        x_vel = vel_hor_now(1);
        y_vel = vel_hor_now(2);
        % -180 < angle < 180
        accel_angle_now = atan2(x_a*y_vel-y_a*x_vel,x_a*x_vel+y_a*y_vel) *180/pi();
        
        At_hor(i,j) = dot(vel_hor_now,accel_hor_now) / norm(vel_hor_now);
        An_hor(i,j) = sign(accel_angle_now)*sqrt(norm(accel_hor_now)^2 - At_hor(i,j)^2);
        accel_angle_hor_vel(i,j) = accel_angle_now;

    end

    for i=2:size(vel,1)
        
        vel_pre(:,1) = vel(i-1,j,:);
        vel_now(:,1) = vel(i,j,:);
        alpha = atan2(norm(cross(vel_now,vel_pre)),dot(vel_now,vel_pre)) *180/pi();
        alpha_dot(i,j) = alpha.*180./pi()./dt;

        vel_hor_pre(:,1) = vel(i-1,j,1:2);
        vel_hor_now(:,1) = vel(i,j,1:2);
        % 2d horizontal velocity vector @ t-dt
        x1 = vel_hor_pre(1);
        y1 = vel_hor_pre(2);
        % 2d horizontal velocity vector @ t
        x2 = vel_hor_now(1);
        y2 = vel_hor_now(2);
        % -180 < angle < 180
        alpha_hor = atan2(x1*y2-y1*x2,x1*x2+y1*y2) *180/pi();
        alpha_dot_hor(i,j) = alpha_hor./dt;
%         alpha_hor(i,j) = acos(dot(vel_hor_now,vel_hor_pre) / norm(vel_hor_pre) / norm(vel_hor_now));
    end

end

% remove zeros at first row
alpha_dot_hor(1,:) = nan;

