function [roll,pitch,yaw] = qbody2angles(qbody)

for i=1:size(qbody,1)
    for j=1:size(qbody,2)
        
        q_now(:,:) = qbody(i,j,:);

        [yaw_now pitch_now roll_now] = quat2angle([q_now(4);q_now(1:3)]');
        yaw(i,j) = rad2deg(yaw_now);
        pitch(i,j) = rad2deg(pitch_now);
        roll(i,j) = rad2deg(roll_now);
    end
end

% adjust pitch & roll
pitch = - pitch;
roll_temp = roll;
roll(roll_temp<0) = roll_temp(roll_temp<0) +180;
roll(roll_temp>0) = roll_temp(roll_temp>0) -180;

