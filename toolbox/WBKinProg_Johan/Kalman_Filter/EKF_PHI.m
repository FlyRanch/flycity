function PHI = EKF_PHI(x,dt)

%    % Update quaternion:
%     
%         % Quaternion rate matrix:
%     
%         W_q = [x(7) -x(6) x(5) -x(4); ...
%                x(6) x(7) -x(4) -x(5); ...
%                -x(5) x(4) x(7) -x(6)];
%        
%         % Quaternion derivative:
%         
%         q2 = 0.5*W_q'*[x(1); x(2); x(3)];
%         
%         if q2(4) == 0
%             
%             q2(4) = 1;
%             
%         end
%         
%         q2 = q2./norm(q2);
%         
%         q1 = [0; 0; 0; 1];
%         
%         % Now Slerp is used in order to find qdot:
%         
%         if q1'*q2 < 0
%             
%             q2 = -q2;
%             
%             theta = real(acos(q1'*q2));
%             
%             q_dot = (sin((1-dt)*theta)/sin(theta))*q1 + (sin(dt*theta)/sin(theta))*q2;
%             
%             q_dot = q_dot/norm(q_dot);
%             
%         elseif q1'*q2 == 1
%             
%             q_dot = q1;
%             
%         else
%             
%             theta = real(acos(q1'*q2));
%             
%             q_dot = (sin((1-dt)*theta)/sin(theta))*q1 + (sin(dt*theta)/sin(theta))*q2;
%             
%             q_dot = q_dot/norm(q_dot);
%             
%         end
%         
% %         % Update quaternion by multiplying(quat) by the quaternion
% %         % derivative:
% %         
% %         Q = [x(7) -x(6) x(5) x(4); ...
% %              x(6) x(7) -x(4) x(5); ...
% %              -x(5) x(4) x(7) x(6); ...
% %              -x(4) -x(5) -x(6) x(7)];
% %          
% %         q_upd = Q*q_dot;
% %         
% %         q_upd = q_upd./norm(q_upd);      
% %         
% %         
% %      % Update angular rates:
% %      
% %         % Convert omega to predicted reference frame
% %         
% %         omega = 2*W_q*q_dot;



PHI_1  = [1 0 0 0 0 0 0 0 0 0];
PHI_2  = [0 1 0 0 0 0 0 0 0 0];
PHI_3  = [0 0 1 0 0 0 0 0 0 0];
PHI_4  = [dt 0 0 1 0 0 0 0 0 0];
PHI_5  = [0 dt 0 0 1 0 0 0 0 0];
PHI_6  = [0 0 dt 0 0 1 0 0 0 0];
PHI_7  = [0 0 0 x(10)*dt/2 -x(9)*dt/2 x(8)*dt/2 1 x(6)*dt/2 -x(5)*dt/2 x(4)*dt/2];
PHI_8  = [0 0 0 x(9)*dt/2 x(10)*dt/2 -x(7)*dt/2 -x(6)*dt/2 1 x(4)*dt/2 x(5)*dt/2];
PHI_9  = [0 0 0 -x(8)*dt/2 x(7)*dt/2 x(10)*dt/2 x(5)*dt/2 -x(4)*dt/2 1 x(6)*dt/2];
PHI_10 = [0 0 0 -x(7)*dt/2 -x(8)*dt/2 -x(9)*dt/2 -x(4)*dt/2 -x(5)*dt/2 -x(6)*dt/2 1];

PHI = [PHI_1; PHI_2; PHI_3; PHI_4; PHI_5; PHI_6; PHI_7; PHI_8; PHI_9; PHI_10];
        


% PHI_1  = [1 0 0 0 0 0 0 0 0 0];
% PHI_2  = [0 1 0 0 0 0 0 0 0 0];
% PHI_3  = [0 0 1 0 0 0 0 0 0 0];
% PHI_4  = [(1/(2*pi))*dt 0 0 (1/(2*pi))*1 0 0 0 0 0 0];
% PHI_5  = [0 (1/(2*pi))*dt 0 0 (1/(2*pi))*1 0 0 0 0 0];
% PHI_6  = [0 0 (1/(2*pi))*dt 0 0 (1/(2*pi))*1 0 0 0 0];
% PHI_7  = [0 0 0 x(10)*dt/2 -x(9)*dt/2 x(8)*dt/2 1 x(6)*dt/2 -x(5)*dt/2 x(4)*dt/2];
% PHI_8  = [0 0 0 x(9)*dt/2 x(10)*dt/2 -x(7)*dt/2 -x(6)*dt/2 1 x(4)*dt/2 x(5)*dt/2];
% PHI_9  = [0 0 0 -x(8)*dt/2 x(7)*dt/2 x(10)*dt/2 x(5)*dt/2 -x(4)*dt/2 1 x(6)*dt/2];
% PHI_10 = [0 0 0 -x(7)*dt/2 -x(8)*dt/2 -x(9)*dt/2 -x(4)*dt/2 -x(5)*dt/2 -x(6)*dt/2 1];
% 
% PHI = [PHI_1; PHI_2; PHI_3; PHI_4; PHI_5; PHI_6; PHI_7; PHI_8; PHI_9; PHI_10];


% PHI_1  = [0 0 0 0 0 0 0 0 0 0];
% PHI_2  = [0 0 0 0 0 0 0 0 0 0];
% PHI_3  = [0 0 0 0 0 0 0 0 0 0];
% PHI_4  = [0 0 0 1 0 0 0 0 0 0];
% PHI_5  = [0 0 0 0 1 0 0 0 0 0];
% PHI_6  = [0 0 0 0 0 1 0 0 0 0];
% PHI_7  = [0 0 0 x(10)*dt/2 -x(9)*dt/2 x(8)*dt/2 1 x(6)*dt/2 -x(5)*dt/2 x(4)*dt/2];
% PHI_8  = [0 0 0 x(9)*dt/2 x(10)*dt/2 -x(7)*dt/2 -x(6)*dt/2 1 x(4)*dt/2 x(5)*dt/2];
% PHI_9  = [0 0 0 -x(8)*dt/2 x(7)*dt/2 x(10)*dt/2 x(5)*dt/2 -x(4)*dt/2 1 x(6)*dt/2];
% PHI_10 = [0 0 0 -x(7)*dt/2 -x(8)*dt/2 -x(9)*dt/2 -x(4)*dt/2 -x(5)*dt/2 -x(6)*dt/2 1];
% 
% PHI = [PHI_1; PHI_2; PHI_3; PHI_4; PHI_5; PHI_6; PHI_7; PHI_8; PHI_9; PHI_10];






% PHI_1 = [-dt*(q_dot(4)^2+q_dot(1)^2-q_dot(2)^2-q_dot(3)^2) -dt*(2*q_dot(1)*q_dot(2)+2*q_dot(4)*q_dot(3)) -dt*(2*q_dot(1)*q_dot(3)-2*q_dot(4)*q_dot(2)) -2*dt*q_dot(4) 2*dt*q_dot(3) -2*dt*q_dot(2) 2*dt*q_dot(1)];
% PHI_2 = [-dt*(2*q_dot(1)*q_dot(2)-2*q_dot(4)*q_dot(3)) -dt*(q_dot(4)^2-q_dot(1)^2+q_dot(2)^2-q_dot(3)^2) -dt*(2*q_dot(2)*q_dot(3)+2*q_dot(4)*q_dot(1)) -2*dt*q_dot(3) -2*dt*q_dot(4) 2*dt*q_dot(1) 2*dt*q_dot(2)];
% PHI_3 = [-dt*(2*q_dot(1)*q_dot(3)-2*q_dot(4)*q_dot(3)) -dt*(2*q_dot(2)*q_dot(3)-2*q_dot(4)*q_dot(3)) -dt*(q_dot(4)^2-q_dot(1)^2-q_dot(2)^2+q_dot(3)^2) 2*dt*q_dot(2) -2*dt*q_dot(1) -2*dt*q_dot(4) 2*dt*q_dot(3)];
% PHI_4 = [x(7)*dt/2 -x(6)*dt/2 x(5)*dt/2 1 x(3)*dt/2 -x(2)*dt/2 x(1)*dt/2];
% PHI_5 = [x(6)*dt/2 x(7)*dt/2 -x(4)*dt/2 -x(3)*dt/2 1 x(1)*dt/2 x(2)*dt/2];
% PHI_6 = [-x(5)*dt/2 x(4)*dt/2 x(7)*dt/2 x(2)*dt/2 -x(1)*dt/2 1 x(3)*dt/2];
% PHI_7 = [-x(4)*dt/2 -x(5)*dt/2 -x(6)*dt/2 -x(1)*dt/2 -x(2)*dt/2 -x(3)*dt/2 1];
% 
%  
% PHI = [PHI_1; PHI_2; PHI_3; PHI_4; PHI_5; PHI_6; PHI_7];



        
    