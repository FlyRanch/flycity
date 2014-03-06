function [roll,pitch,yaw] = qbody2angles_THAD(qbody)

roll = nan(1,size(qbody,2));
pitch = nan(1,size(qbody,2));
yaw = nan(1,size(qbody,2));

for i=1:size(qbody,2)
        q_now = qbody(:,i);
    
%         roll(i,j) = rad2deg(atan2(-2*(q_now(3)*q_now(1)-q_now(2)*q_now(4)),1-2*(q_now(1)^2+q_now(2)^2)));
%         pitch(i,j) = rad2deg(asin( (1-2*(q_now(1)^2 + q_now(2)^2)) / (-2*(q_now(1)*q_now(2)-q_now(3)*q_now(4))*cosd(roll(i,j))) ));
%         pitch(i,j) = rad2deg(asin(-2*(q_now(1)*q_now(2)-q_now(3)*q_now(4)) * ( (1-2*(q_now(1)^2 + q_now(2)^2))/cosd(roll(i,j)) ) ));
%         pitch(i,j) = acosd( (1-2*(q_now(1)^2 + q_now(2)^2)) / cosd(yaw(i,j)));
%         yaw(i,j) = rad2deg(asin(2*(q_now(4)*q_now(1)+q_now(2)*q_now(3))));

        pitch(i) = atand(-2*(q_now(1)*q_now(3)-q_now(2)*q_now(4))/(1-2*(q_now(1)^2+q_now(2)^2)));
        yaw(i) = asind(-2*(q_now(1)*q_now(2)-q_now(3)*q_now(4)) * (1-2*(q_now(1)^2 + q_now(2)^2)) / cosd(pitch(i)) );
        roll(i) = asind(2*(q_now(1)*q_now(4) + q_now(2)*q_now(3)));
    end
end


