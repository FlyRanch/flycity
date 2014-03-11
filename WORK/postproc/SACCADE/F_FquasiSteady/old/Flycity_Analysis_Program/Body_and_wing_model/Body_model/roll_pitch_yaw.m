function [roll,pitch,yaw] = roll_pitch_yaw(qB)
    % Compute the pitch, roll and yaw angle in radians with the body
    % quaternion:
    
    roll    = atan2(2*qB(2)*qB(3)+2*qB(4)*qB(1),qB(3)^2-qB(2)^2-qB(1)^2+qB(4)^2);
    pitch   = -asin(2*qB(1)*qB(3)-qB(4)*qB(2));
    yaw     = atan2(2*qB(1)*qB(2)+2*qB(4)*qB(3),qB(1)^2+qB(4)^2-qB(3)^2-qB(2)^2);

end

