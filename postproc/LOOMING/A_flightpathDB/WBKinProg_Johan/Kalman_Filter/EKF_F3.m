function x = EKF_F3(x,x_prev,dt)

    % Calculate i+1 values for q, omega and omega_dot
    
    % Step 1: calculate q_dot in the body frame of reference
    
        W_q = [x(10) x(9) -x(8) -x(7); ...
               -x(9) x(10) x(7) -x(8); ...
               x(8) -x(7) x(10) -x(9)];
           
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
        
    % Step 2: Update the quaternion in the world frame of reference
        
        % Now update the quaternion q, q_dot_dt is a quaternion which
        % rotates q_i to q_i+1 in the body frame of reference. As the
        % quaternions represent orientations in the world frame of
        % reference q_dot_dt needs to be transferred to the world frame of
        % reference first.
        
        q_inv = [x(7); x(8); x(9); -x(10)]./norm([x(7); x(8); x(9); x(10)]);
        
        Q_inv = [ q_inv(4)  q_inv(3) -q_inv(2) q_inv(1); ... %quaternion multiplication matrix
                 -q_inv(3)  q_inv(4)  q_inv(1) q_inv(2); ...
                  q_inv(2) -q_inv(1)  q_inv(4) q_inv(3); ...
                 -q_inv(1) -q_inv(2) -q_inv(3) q_inv(4)];
             
        Q_bar = [ x(10) -x(9)   x(8)  x(7); ...  %quaternion multiplication conjugate matrix
                  x(9)   x(10) -x(7)  x(8); ...
                 -x(8)   x(7)   x(10) x(9); ...
                 -x(7)  -x(8)  -x(9)  x(10)];
          
        q_upd = Q_bar*(Q_inv*q_dot_dt);
        
   % Step 3: derive the angular velocity update in the body frame of
   % reference
   
        % The new omega can be simply calculated by adding the angular rate
        % to the old omega. However the body frame of reference changes in
        % orientation. In the world frame of reference this doesn't have
        % much effect on omega, however in the body frame of reference it
        % is important to transfer the updated vector omega to the updated
        % reference frame:
        
        omega_new = [x(4)+x(1)*dt; x(5)+x(2)*dt; x(6)+x(3)*dt];
        
        % Rotation matrix:
        
        R_q_dot = [ q_dot_dt(4)^2+q_dot_dt(1)^2-q_dot_dt(2)^2-q_dot_dt(3)^2    2*q_dot_dt(1)*q_dot_dt(2)+2*q_dot_dt(4)*q_dot_dt(3)      2*q_dot_dt(1)*q_dot_dt(3)-2*q_dot_dt(4)*q_dot_dt(2); ...
                    2*q_dot_dt(1)*q_dot_dt(2)-2*q_dot_dt(4)*q_dot_dt(3)        q_dot_dt(4)^2-q_dot_dt(1)^2+q_dot_dt(2)^2-q_dot_dt(3)^2  2*q_dot_dt(2)*q_dot_dt(3)+2*q_dot_dt(4)*q_dot_dt(1); ...
                    2*q_dot_dt(1)*q_dot_dt(3)+2*q_dot_dt(4)*q_dot_dt(2)        2*q_dot_dt(2)*q_dot_dt(3)-2*q_dot_dt(4)*q_dot_dt(1)      q_dot_dt(4)^2-q_dot_dt(1)^2-q_dot_dt(2)^2+q_dot_dt(3)^2];
                
        % Omega update
        
        omega_upd = R_q_dot*omega_new;
        
   % Step 4: Finally a new update for omega_dot is needed. Also omega_dot
   % changes in orientation and 2 values for omega are needed in order to
   % obtain a new value for omega_dot. Therefore the value of omega from
   % i-1 is necessary. This value is supplied from the EKF filter. It is
   % important to keep in mind that the reference frames have changed
   % between omega_i and omega_i-1. Therefore first omega_i, omega_i-1 and
   % omega_upd are transferred to the world coordinate system, after which
   % omega_dot is calculated. Than omega_dot is transferred to the
   % predicted attitude of the body frame of reference.
   
   
        % Rotation matrices
   
        R_upd = [ q_upd(4)^2+q_upd(1)^2-q_upd(2)^2-q_upd(3)^2    2*q_upd(1)*q_upd(2)+2*q_upd(4)*q_upd(3)      2*q_upd(1)*q_upd(3)-2*q_upd(4)*q_upd(2); ...
                    2*q_upd(1)*q_upd(2)-2*q_upd(4)*q_upd(3)        q_upd(4)^2-q_upd(1)^2+q_upd(2)^2-q_upd(3)^2  2*q_upd(2)*q_upd(3)+2*q_upd(4)*q_upd(1); ...
                    2*q_upd(1)*q_upd(3)+2*q_upd(4)*q_upd(2)        2*q_upd(2)*q_upd(3)-2*q_upd(4)*q_upd(1)      q_upd(4)^2-q_upd(1)^2-q_upd(2)^2+q_upd(3)^2];
                
        R_i   = [ x(10)^2+x(7)^2-x(8)^2-x(9)^2    2*x(7)*x(8)+2*x(10)*x(9)      2*x(7)*x(9)-2*x(10)*x(8); ...
                  2*x(7)*x(8)-2*x(10)*x(9)        x(10)^2-x(7)^2+x(8)^2-x(9)^2  2*x(8)*x(9)+2*x(10)*x(7); ...
                  2*x(7)*x(9)+2*x(10)*x(8)        2*x(8)*x(9)-2*x(10)*x(7)      x(10)^2-x(7)^2-x(8)^2+x(9)^2];
                
        R_i_1 = [ x_prev(10)^2+x_prev(7)^2-x_prev(8)^2-x_prev(9)^2    2*x_prev(7)*x_prev(8)+2*x_prev(10)*x_prev(9)      2*x_prev(7)*x_prev(9)-2*x_prev(10)*x_prev(8); ...
                  2*x_prev(7)*x_prev(8)-2*x_prev(10)*x_prev(9)        x_prev(10)^2-x_prev(7)^2+x_prev(8)^2-x_prev(9)^2  2*x_prev(8)*x_prev(9)+2*x_prev(10)*x_prev(7); ...
                  2*x_prev(7)*x_prev(9)+2*x_prev(10)*x_prev(8)        2*x_prev(8)*x_prev(9)-2*x_prev(10)*x_prev(7)      x_prev(10)^2-x_prev(7)^2-x_prev(8)^2+x_prev(9)^2];
              
        % angular velocities in world frame
   
        omega_i_1 = R_upd' * [ x_prev(4); x_prev(5); x_prev(6) ];
        omega_i   = R_i'   * [ x(4); x(5); x(6) ];
        omega_upd = R_i_1' * [ omega_upd(1); omega_upd(2); omega_upd(3)];
        
        
        % omega_dot update in world frame
        
        omega_dot = (omega_upd-2.*omega_i+omega_i_1)./(2*dt);
        
        
   
    

end

