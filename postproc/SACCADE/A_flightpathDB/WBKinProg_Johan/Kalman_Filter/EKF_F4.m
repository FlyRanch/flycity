function x = EKF_F4(x,dt)



    % Step 1: Transfer omega_i, omega_dot_i and omega_dot_i_1 from body
    % fixed reference frame to the world frame:
    
        R_i   = [ x(10)^2+x(7)^2-x(8)^2-x(9)^2    2*x(7)*x(8)+2*x(10)*x(9)      2*x(7)*x(9)-2*x(10)*x(8); ...
                  2*x(7)*x(8)-2*x(10)*x(9)        x(10)^2-x(7)^2+x(8)^2-x(9)^2  2*x(8)*x(9)+2*x(10)*x(7); ...
                  2*x(7)*x(9)+2*x(10)*x(8)        2*x(8)*x(9)-2*x(10)*x(7)      x(10)^2-x(7)^2-x(8)^2+x(9)^2];    
              
%         R_i_1 = [ x_prev(10)^2+x_prev(7)^2-x_prev(8)^2-x_prev(9)^2    2*x_prev(7)*x_prev(8)+2*x_prev(10)*x_prev(9)      2*x_prev(7)*x_prev(9)-2*x_prev(10)*x_prev(8); ...
%                   2*x_prev(7)*x_prev(8)-2*x_prev(10)*x_prev(9)        x_prev(10)^2-x_prev(7)^2+x_prev(8)^2-x_prev(9)^2  2*x_prev(8)*x_prev(9)+2*x_prev(10)*x_prev(7); ...
%                   2*x_prev(7)*x_prev(9)+2*x_prev(10)*x_prev(8)        2*x_prev(8)*x_prev(9)-2*x_prev(10)*x_prev(7)      x_prev(10)^2-x_prev(7)^2-x_prev(8)^2+x_prev(9)^2];
              
        omega_i   = R_i'*[x(4); x(5); x(6)];
        
%         omega_i_1 = R_i'*[x_prev(4); x_prev(5); x_prev(6)];
        
        omega_dot_i = R_i'*[x(1); x(2); x(3)];
        
%         omega_dot_i_1 = R_i_1'*[x_prev(4); x_prev(5); x_prev(6)];
        
     % Step 2: Calculate omega_i+1 with a two-step Adams-Bashforth method:
     
        %omega_upd = omega_i + 1.5.*dt.*omega_dot_i+0.5.*dt.*omega_dot_i_1;
        
       
     % Step 3: Calculate the average omega between t and t+dt needed to
     % calculate the quaternion rate:
     
        %omega_avg = (omega_i+omega_upd)./2;
        
        omega_avg = omega_i +0.5.*omega_dot_i.*dt;
        
     % Step 4: Calculate the quaternion rate and update q:
     
        W_q = [x(10) -x(9) x(8) -x(7); ...
               x(9) x(10) -x(7) -x(8); ...
               -x(8) x(7) x(10) -x(9)];
           
        %q_dot = 0.5*W_q'*[x(4); x(5); x(6)];

        q_dot = 0.5*W_q'*[x(4)+0.5*x(1)*dt; x(5)+0.5*x(2)*dt; x(6)+0.5*x(3)*dt];
    
        % Use Spherical Linear Interpolation (SLERP) to calculate the
        % rotation (captured in q_dot_dt) following from angular rate omega
        % in order to update the orientation.
        
        q2 = q_dot;
        
        if q2(4) == 0
            
            q2 = [0; 0; 0; 1];
            
        end
        
        q2 = q2./norm(q2);
        
        q1 = [0; 0; 0; 1];
                
        c = dt;
        
        if q1'*q2 < 0
            
            q2 = -q2;
            
            theta = real(acos(q1'*q2));
            
            q_dot1 = (sin((1-c)*theta)/sin(theta))*q1 + (sin(c*theta)/sin(theta))*q2;
            
            q_dot_dt = q_dot1/norm(q_dot1);
            
        elseif q1'*q2 == 1
            
            q_dot_dt = q1;
            
        else
            
            theta = real(acos(q1'*q2));
            
            q_dot1 = (sin((1-c)*theta)/sin(theta))*q1 + (sin(c*theta)/sin(theta))*q2;
            
            q_dot_dt = q_dot1/norm(q_dot1);
            
        end
        
        clear q1 q2 q_dot1 theta     
        
        Q = [ x(10) x(9) -x(8)  x(7); ...  %quaternion multiplication conjugate matrix
             -x(9)  x(10) x(7)  x(8); ...
              x(8) -x(7)  x(10) x(9); ...
             -x(7)  -x(8)  -x(9)  x(10)];
          
        q_upd = Q*q_dot_dt;
        
     % Step 5: Update omega_dot in the world frame:
     
                phi = 2*acos(q_upd(4));
        
                if abs(phi) > 0
                
                omega_t(1) = phi*(q_upd(1)/sin(0.5*phi));
                omega_t(2) = phi*(q_upd(2)/sin(0.5*phi));
                omega_t(3) = phi*(q_upd(3)/sin(0.5*phi));
                
                else
                    
                omega_t(1) = 0;
                omega_t(2) = 0;
                omega_t(3) = 0;
                
                end
                
        %omega_upd = omega_i; % + omega_dot_i.*dt;
                
        %omega_upd = 0.5.*[omega_t(1); omega_t(2); omega_t(3)] + 0.5.*(omega_i + 1.5.*dt.*omega_dot_i+0.5.*dt.*omega_dot_i_1);
        
        %omega_upd = (omega_i + 1.5.*dt.*omega_dot_i+0.5.*dt.*omega_dot_i_1);
                
        %omega_upd = omega_i+[omega_t(1); omega_t(2); omega_t(3)].*dt;
        
        %omega_upd = [omega_t(1); omega_t(2); omega_t(3)];
    
        %omega_dot_upd = 0.5.*omega_dot_i+0.5.*(omega_upd-2*omega_i+omega_i_1)./dt;
        
        %omega_dot_upd = (omega_upd-omega_dot_i)./dt;
        
        %omega_dot_upd = omega_dot_i;
        
        %omega_dot_upd = (omega_upd-omega_i)./(dt);
        
     % Step 6: Transfer omega_dot_upd and omega_upd to the body reference
     % frame:
     
%         R_upd = [ q_upd(4)^2+q_upd(1)^2-q_upd(2)^2-q_upd(3)^2    2*q_upd(1)*q_upd(2)+2*q_upd(4)*q_upd(3)      2*q_upd(1)*q_upd(3)-2*q_upd(4)*q_upd(2); ...
%                     2*q_upd(1)*q_upd(2)-2*q_upd(4)*q_upd(3)        q_upd(4)^2-q_upd(1)^2+q_upd(2)^2-q_upd(3)^2  2*q_upd(2)*q_upd(3)+2*q_upd(4)*q_upd(1); ...
%                     2*q_upd(1)*q_upd(3)+2*q_upd(4)*q_upd(2)        2*q_upd(2)*q_upd(3)-2*q_upd(4)*q_upd(1)      q_upd(4)^2-q_upd(1)^2-q_upd(2)^2+q_upd(3)^2];
%                 
%         omega = R_upd*omega_upd;
%         
%         omega_dot = R_upd*omega_dot_upd;
        
        omega = [x(4); x(5); x(6)] + dt.*[x(1); x(2); x(3)];
        
        omega_dot = [x(1); x(2); x(3)];


     % Step 7: Return the values:
     
        x(1)  = omega_dot(1);
        x(2)  = omega_dot(2);
        x(3)  = omega_dot(3);
        x(4)  = omega(1);
        x(5)  = omega(2);
        x(6)  = omega(3);
        x(7)  = q_upd(1);
        x(8)  = q_upd(2);
        x(9)  = q_upd(3);
        x(10) = q_upd(4);
        
        

end

