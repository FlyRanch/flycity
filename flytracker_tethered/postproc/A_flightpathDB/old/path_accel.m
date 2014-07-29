function [At,An,alpha_dot] = path_accel(vel,accel,t)


alpha_dot = nan(size(vel,1),size(vel,2));

for j=1:size(vel,2)
    counter = size(vel,2) - j
    
    for i=1:size(vel,1)
        
        vel_now(:,1) = vel(i,j,:);
        accel_now(:,1) = accel(i,j,:);
        
        At(i,j) = dot(vel_now,accel_now) / norm(vel_now);
        An(i,j) = sqrt(norm(accel_now)^2 - At(i,j)^2);
    end

    for i=2:size(vel,1)
        
        vel_pre(:,1) = vel(i-1,j,:);
        vel_now(:,1) = vel(i,j,:);

        alpha_dot(i,j) = acosd(dot(vel_now,vel_pre) / norm(vel_pre) / norm(vel_now));
    end

end

alpha_dot = alpha_dot ./ (t(2)-t(1));

At = At;
An = An;
alpha_dot = alpha_dot;
