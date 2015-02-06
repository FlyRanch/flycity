function [roll,pitch_global,yaw_global,pitch_true,yaw_true,up,dir] = qbody2angles(qbody,vel)

for i=1:size(qbody,1)
    for j=1:size(qbody,2)
        
        q_now(:,:) = qbody(i,j,:);
        U_now(:,:) = vel(i,j,:);

    [yaw_now,pitch_now,roll_now] = quat2angle([q_now(4);q_now(1:3)]');
    yaw_global(i,j) = rad2deg(yaw_now);
    pitch_global(i,j) = rad2deg(pitch_now);
    roll(i,j) = rad2deg(roll_now);
    
    up(i,j) = atand(U_now(3) / sqrt(U_now(1)^2 +U_now(2)^2));
    dir(i,j) = atand(U_now(1) / U_now(2)) - 180;
end

% adjust pitch & roll
pitch_global = - pitch_global;
roll_temp = roll;
roll(roll_temp<0) = roll_temp(roll_temp<0) +180;
roll(roll_temp>0) = roll_temp(roll_temp>0) -180;

pitch_true = pitch_global - up;
yaw_true = yaw_global - dir;
