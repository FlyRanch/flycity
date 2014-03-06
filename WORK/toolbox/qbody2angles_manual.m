function [roll,pitch,yaw] = qbody2angles_manual(qbody)

for i=1:size(qbody,1)
    for j=1:size(qbody,2)
        
        q_now(:,:) = qbody(i,j,:);
    
        roll(i,j) = rad2deg(atan2(2*(q_now(4)*q_now(1)+q_now(2)*q_now(3)),1-2*(q_now(1)^2+q_now(2)^2)));
        pitch(i,j) = rad2deg(asin(2*(q_now(4)*q_now(2)-q_now(3)*q_now(1))));
        yaw(i,j) = rad2deg(atan2(2*(q_now(4)*q_now(3)+q_now(1)*q_now(2)),1-2*(q_now(2)^2+q_now(3)^2)));
    end
end


