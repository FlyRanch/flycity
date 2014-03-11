function [At,An,alpha_dot,At_hor,An_hor,alpha_dot_hor] = path_accel_incHor(vel,accel,t)

alpha_dot = nan(size(vel,1),size(vel,2));

for j=1:size(vel,2)
    counter = size(vel,2) - j
    
    for i=1:size(vel,1)
        
        vel_now(:,1) = vel(i,j,:);
        accel_now(:,1) = accel(i,j,:);
        At(i,j) = dot(vel_now,accel_now) / norm(vel_now);
        An(i,j) = sqrt(norm(accel_now)^2 - At(i,j)^2);
        
        vel_hor_now(:,1) = vel(i,j,1:2);
        accel_hor_now(:,1) = accel(i,j,1:2);
        At_hor(i,j) = dot(vel_hor_now,accel_hor_now) / norm(vel_hor_now);
        An_hor(i,j) = sqrt(norm(accel_hor_now)^2 - At_hor(i,j)^2);
    end

    for i=2:size(vel,1)
        
        vel_pre(:,1) = vel(i-1,j,:);
        vel_now(:,1) = vel(i,j,:);
        alpha(i,j) = atan2(norm(cross(vel_now,vel_pre)),dot(vel_now,vel_pre));

        vel_hor_pre(:,1) = vel(i-1,j,1:2);
        vel_hor_now(:,1) = vel(i,j,1:2);
        alpha_hor(i,j) = acos(dot(vel_hor_now,vel_hor_pre) / norm(vel_hor_pre) / norm(vel_hor_now));
    end

end

alpha_dot = rad2deg(alpha) ./ (t(2)-t(1));
alpha_dot_hor = rad2deg(alpha_hor) ./ (t(2)-t(1));
