function [roll,pitch,yaw] = qbody2angles_manual_yawrollpitch(qbody)

roll = nan(size(qbody,1),size(qbody,2));
pitch = nan(size(qbody,1),size(qbody,2));
yaw = nan(size(qbody,1),size(qbody,2));

for i=1:size(qbody,1)
    for j=1:size(qbody,2)
        
        q_now(:,:) = qbody(i,j,:);
    
%         roll(i,j) = rad2deg(atan2(-2*(q_now(3)*q_now(1)-q_now(2)*q_now(4)),1-2*(q_now(1)^2+q_now(2)^2)));
%         pitch(i,j) = rad2deg(asin( (1-2*(q_now(1)^2 + q_now(2)^2)) / (-2*(q_now(1)*q_now(2)-q_now(3)*q_now(4))*cosd(roll(i,j))) ));
%         pitch(i,j) = rad2deg(asin(-2*(q_now(1)*q_now(2)-q_now(3)*q_now(4)) * ( (1-2*(q_now(1)^2 + q_now(2)^2))/cosd(roll(i,j)) ) ));
%         pitch(i,j) = acosd( (1-2*(q_now(1)^2 + q_now(2)^2)) / cosd(yaw(i,j)));
%         yaw(i,j) = rad2deg(asin(2*(q_now(4)*q_now(1)+q_now(2)*q_now(3))));

        pitch(i,j) = atand(-2*(q_now(1)*q_now(3)-q_now(2)*q_now(4))/(1-2*(q_now(1)^2+q_now(2)^2)));
        yaw(i,j) = asind(-2*(q_now(1)*q_now(2)-q_now(3)*q_now(4)) * (1-2*(q_now(1)^2 + q_now(2)^2)) / cosd(pitch(i,j)) );
        roll(i,j) = asind(2*(q_now(1)*q_now(4) + q_now(2)*q_now(3)));
    end
end


