function x = EKF_F(x,dt)
    
    % Update quaternion:
    
        % Quaternion rate matrix:
    
        W_q = [x(10) -x(9) x(8) -x(7); ...
               x(9) x(10) -x(7) -x(8); ...
               -x(8) x(7) x(10) -x(9)];
       
%         W_q = [x(10) x(9) -x(8) -x(7); ...
%                -x(9) x(10) x(7) -x(8); ...
%                x(8) -x(7) x(10) -x(9)];

        % Quaternion derivative:
        
        q2 = 0.5*W_q'*[x(4); x(5); x(6)];
        
        if q2(4) == 0
            
            q2(4) = 1;
            
        end
        
        q2 = q2./norm(q2);
        
        q1 = [0; 0; 0; 1];
        
        % Now Slerp is used in order to find qdot:
        
        c = dt;
        
        if q1'*q2 < 0
            
            q2 = -q2;
            
            theta = real(acos(q1'*q2));
            
            q_dot1 = (sin((1-c)*theta)/sin(theta))*q1 + (sin(c*theta)/sin(theta))*q2;
            
            q_dot1 = q_dot1/norm(q_dot1);
            
            q_dot2 = (sin((1-0.5*c)*theta)/sin(theta))*q1 + (sin(0.5*c*theta)/sin(theta))*q2;
            
            q_dot2 = q_dot2/norm(q_dot2);
            
        elseif q1'*q2 == 1
            
            q_dot1 = q1;
            
            q_dot2 = q1;
            
        else
            
            theta = real(acos(q1'*q2));
            
            q_dot1 = (sin((1-c)*theta)/sin(theta))*q1 + (sin(c*theta)/sin(theta))*q2;
            
            q_dot1 = q_dot1/norm(q_dot1);
            
            q_dot2 = (sin((1-0.5*c)*theta)/sin(theta))*q1 + (sin(0.5*c*theta)/sin(theta))*q2;
            
            q_dot2 = q_dot2/norm(q_dot2);
            
        end
        
        % Update quaternion by multiplying(quat) by the quaternion
        % derivative:
        
%         Q = [x(10) -x(9) x(8) x(7); ...
%              x(9) x(10) -x(7) x(8); ...
%              -x(8) x(7) x(10) x(9); ...
%              -x(7) -x(8) -x(9) x(10)];

        Q = [x(10) x(9) -x(8) x(7); ...
             -x(9) x(10) x(7) x(8); ...
             x(8) -x(7) x(10) x(9); ...
             -x(7) -x(8) -x(9) x(10)];
         
        q_upd = Q*q_dot1;

        
        q_upd = q_upd./norm(q_upd);
        
        
        % Calculate omega_dot
               
        q_upd1 = Q*q_dot1;
        
        q_upd1 = q_upd1./norm(q_upd1);
        
        q_upd2 = Q*q_dot1;
        
        q_upd2 = q_upd2./norm(q_upd2);
        
        W_q_upd1 = [q_upd1(4) -q_upd1(3) q_upd1(2) -q_upd1(1); ...
                   q_upd1(3) q_upd1(4) -q_upd1(1) -q_upd1(2); ...
                   -q_upd1(2) q_upd1(1) q_upd1(4) -q_upd1(3)];
               
%         W_q_upd1 = [q_upd1(4) q_upd1(3) -q_upd1(2) -q_upd1(1); ...
%                     -q_upd1(3) q_upd1(4) q_upd1(1) -q_upd1(2); ...
%                     q_upd1(2) -q_upd1(1) q_upd1(4) -q_upd1(3)];
               
        W_q_upd2 = [q_upd2(4) -q_upd2(3) q_upd2(2) -q_upd2(1); ...
                   q_upd2(3) q_upd2(4) -q_upd2(1) -q_upd2(2); ...
                   -q_upd2(2) q_upd2(1) q_upd2(4) -q_upd2(3)];

%         W_q_upd2 = [q_upd2(4) q_upd2(3) -q_upd2(2) -q_upd2(1); ...
%                     -q_upd2(3) q_upd2(4) q_upd2(1) -q_upd2(2); ...
%                     q_upd2(2) -q_upd2(1) q_upd2(4) -q_upd2(3)];
        
        omega1 = 2*W_q_upd1*q_dot1;
        
        omega2 = 2*W_q_upd2*q_dot2;
        
%        omega_dot = 2.*(omega1-omega2)./(dt);

%        omega_dot = 
        
        R_q_dot = [ q_dot1(4)^2+q_dot1(1)^2-q_dot1(2)^2-q_dot1(3)^2    2*q_dot1(1)*q_dot1(2)+2*q_dot1(4)*q_dot1(3)      2*q_dot1(1)*q_dot1(3)-2*q_dot1(4)*q_dot1(2); ...
                    2*q_dot1(1)*q_dot1(2)-2*q_dot1(4)*q_dot1(3)        q_dot1(4)^2-q_dot1(1)^2+q_dot1(2)^2-q_dot1(3)^2  2*q_dot1(2)*q_dot1(3)+2*q_dot1(4)*q_dot1(1); ...
                    2*q_dot1(1)*q_dot1(3)+2*q_dot1(4)*q_dot1(2)        2*q_dot1(2)*q_dot1(3)-2*q_dot1(4)*q_dot1(1)      q_dot1(4)^2-q_dot1(1)^2-q_dot1(2)^2+q_dot1(3)^2];
        
        %omega_dot = (omega1-2*omega2+[x(4); x(5); x(6)])./dt;
        
        %omega_dot = (2*W_q_upd1*q_dot1 - [x(4); x(5); x(6)])./dt;
        
        % Calculate omega
        
        
       %omega = 2*W_q_upd1*q_dot1;
       
       %omega = 2*W_q*q_dot1;
        
        %omega = (omega1+omega2)./2;

        %omega = q2 ; %+ [omega_dot(1)*dt; omega_dot(2)*dt; omega_dot(3)*dt];
        
        %omega = [x(4); x(5); x(6)] + [omega_dot(1)*dt; omega_dot(2)*dt; omega_dot(3)*dt];
        
        %omega = [x(4); x(5); x(6)] + [omega_dot(1)*dt; omega_dot(2)*dt; omega_dot(3)*dt];
        
        %omega = W_q*q_dot1 + W_q*q_dot2;

        %omega = [x(4); x(5); x(6)];
        

        

        
%         if abs(phi) > 0
%         
%         omega(1) = phi*(q_dot1(1)/sin(0.5*phi))+x(1)*dt;
%         omega(2) = phi*(q_dot1(2)/sin(0.5*phi))+x(2)*dt;
%         omega(3) = phi*(q_dot1(3)/sin(0.5*phi))+x(3)*dt;
%         
%         else
%             
%         omega(1) = x(1)*dt;
%         omega(2) = x(2)*dt;
%         omega(3) = x(3)*dt;
%         
%         end

%         phi = 2*acos(q_dot1(4));
% 
%         if abs(phi) > 0
%         
%         omega(1) = phi*(q_dot1(1)/sin(0.5*phi));
%         omega(2) = phi*(q_dot1(2)/sin(0.5*phi));
%         omega(3) = phi*(q_dot1(3)/sin(0.5*phi));
%         
%         else
%             
%         omega(1) = 0;
%         omega(2) = 0;
%         omega(3) = 0;
%         
%         end

%        q3 = q_upd-x(7:10);
%        
%        q3 = q3./norm(q3);
% 
%         phi = 2*acos(q3(4));
% 
%         if abs(phi) > 0
%         
%         omega(1) = phi*(q3(1)/sin(0.5*phi));
%         omega(2) = phi*(q3(2)/sin(0.5*phi));
%         omega(3) = phi*(q3(3)/sin(0.5*phi));
%         
%         else
%             
%         omega(1) = 0;
%         omega(2) = 0;
%         omega(3) = 0;
%         
%         end

        omega = [x(4); x(5); x(6)] + [x(1); x(2); x(3)].*dt;
        
        omega_dot = ([omega(1); omega(2); omega(3)] - [x(4); x(5); x(6)])./dt;
        
        

       
     % Return new value of x
     
     x(1) = omega_dot(1);
     x(2) = omega_dot(2);
     x(3) = omega_dot(3);
     x(4) = omega(1);
     x(5) = omega(2);
     x(6) = omega(3);
     x(7) = q_upd(1);
     x(8) = q_upd(2);
     x(9) = q_upd(3);
     x(10) = q_upd(4);

        
     
          
    
    
end
    
    