function x = EKF_F2(x,dt)
    
    % Update quaternion:
    
        % Quaternion rate matrix:
    
        W_q = [x(10) -x(9) x(8) -x(7); ...
               x(9) x(10) -x(7) -x(8); ...
               -x(8) x(7) x(10) -x(9)];
       
        % Quaternion derivative:
        
        q2 = 0.5*W_q'*[x(4)*dt; x(5)*dt; x(6)*dt];
        
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
        
        Q = [x(10) x(9) -x(8) x(7); ...
             -x(9) x(10) x(7) x(8); ...
             x(8) -x(7) x(10) x(9); ...
             -x(7) -x(8) -x(9) x(10)];
         
        q_upd = Q*q2;
        
        q_upd = q_upd./norm(q_upd);
        
        
        % Calculate omega
        
        rn = real(atan2(2*(q_upd(4).*q_upd(1)+q_upd(2).*q_upd(3)),1-2*(q_upd(1).^2+q_upd(2).^2)));
        pn = real(asin(2*(q_upd(4).*q_upd(2)-q_upd(3).*q_upd(1))));
        yn = real(atan2(2*(q_upd(4).*q_upd(3)+q_upd(1).*q_upd(2)),1-2*(q_upd(2).^2+q_upd(3).^2)));
        
        ro = real(atan2(2*(x(10).*x(7)+x(8).*x(9)),1-2*(x(7).^2+x(8).^2)));
        po = real(asin(2*(x(10).*x(8)-x(9).*x(7))));
        yo = real(atan2(2*(x(10).*x(9)+x(7).*x(8)),1-2*(x(8).^2+x(9).^2)));
        
        om_r = (rn-ro)/dt;
        om_p = (pn-po)/dt;
        om_y = (yn-yo)/dt;
        
        Tr = [1 0 0; 0 cos(rn) sin(rn); 0 -sin(rn) cos(rn)];
        Tp = [cos(pn) 0 -sin(pn); 0 1 0; sin(pn) 0 cos(pn)];
        Ty = [cos(yn) -sin(yn) 0; sin(yn) cos(yn) 0; 0 0 1];

        om_g = (Tr*Tp*[0; 0; om_y]+Tr*[0; om_p; 0] + [om_r; 0; 0]);

       
        omega= om_g;
        
        omega_dot = [x(1); x(2); x(3)];
        
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