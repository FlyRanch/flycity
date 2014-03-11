function x = EKF_F5(x, x_m, dt )

    % Update the quaternion:
    
        % Calculate the DCM matrix of current quaternion:
        
            R_i   = [   x(10)^2+x(7)^2-x(8)^2-x(9)^2    2*x(7)*x(8)+2*x(10)*x(9)      2*x(7)*x(9)-2*x(10)*x(8); ...
                        2*x(7)*x(8)-2*x(10)*x(9)        x(10)^2-x(7)^2+x(8)^2-x(9)^2  2*x(8)*x(9)+2*x(10)*x(7); ...
                        2*x(7)*x(9)+2*x(10)*x(8)        2*x(8)*x(9)-2*x(10)*x(7)      x(10)^2-x(7)^2-x(8)^2+x(9)^2];
        
        % Calculate the update rotation vector in world reference frame:
        
            d_eps = R_i'*[x(4)*dt; x(5)*dt; x(6)*dt];
            
        % Calculate the update quaternion:
            
            if norm(d_eps) == 0
                
                d_q = [0; 0; 0; 1];
                
            else
                
                d_q = [(sin(norm(d_eps)/2)/(norm(d_eps)/2)).*(d_eps./2); cos(norm(d_eps)/2)];
                
            end
            
        % Calculate the quaternion multiplication matrix:
        
            Q_d_q = [ d_q(4)  d_q(3) -d_q(2)  d_q(1); ...
                     -d_q(3)  d_q(4)  d_q(1)  d_q(2); ...
                      d_q(2) -d_q(1)  d_q(4)  d_q(3); ...
                     -d_q(1) -d_q(2) -d_q(3)  d_q(4)];
                 
%             Q_d_q-eye(4)
            
        % Update and normalize the quaternion:
        
            q_temp = Q_d_q*[x(7); x(8); x(9); x(10)];
            
            q_upd = q_temp./norm(q_temp);
            
    % Update omega:
    
        % Calculate the DCM matrix of the updated quaternion:
        
            R_ip1 = [   q_upd(4)^2+q_upd(1)^2-q_upd(2)^2-q_upd(3)^2    2*q_upd(1)*q_upd(2)+2*q_upd(4)*q_upd(3)      2*q_upd(1)*q_upd(3)-2*q_upd(4)*q_upd(2); ...
                        2*q_upd(1)*q_upd(2)-2*q_upd(4)*q_upd(3)        q_upd(4)^2-q_upd(1)^2+q_upd(2)^2-q_upd(3)^2  2*q_upd(2)*q_upd(3)+2*q_upd(4)*q_upd(1); ...
                        2*q_upd(1)*q_upd(3)+2*q_upd(4)*q_upd(2)        2*q_upd(2)*q_upd(3)-2*q_upd(4)*q_upd(1)      q_upd(4)^2-q_upd(1)^2-q_upd(2)^2+q_upd(3)^2];
                    
        % Calculate the DCM matrix of the previous quaternion:
        
            R_im1 = [   x_m(10)^2+x_m(7)^2-x_m(8)^2-x_m(9)^2    2*x_m(7)*x_m(8)+2*x_m(10)*x_m(9)      2*x_m(7)*x_m(9)-2*x_m(10)*x_m(8); ...
                        2*x_m(7)*x_m(8)-2*x_m(10)*x_m(9)        x_m(10)^2-x_m(7)^2+x_m(8)^2-x_m(9)^2  2*x_m(8)*x_m(9)+2*x_m(10)*x_m(7); ...
                        2*x_m(7)*x_m(9)+2*x_m(10)*x_m(8)        2*x_m(8)*x_m(9)-2*x_m(10)*x_m(7)      x_m(10)^2-x_m(7)^2-x_m(8)^2+x_m(9)^2];
                    
        % Calculate the derivitve of the rotation matrix between t = t and
        % t = t + dt:
        
            DR_Dt = (R_i-R_im1)./(dt);
            
        % Calculate the cross matrix of omega:
         
            W_cross = R_i'*DR_Dt;
        
        % Update omega with omega_dot and the relative rotation to the new reference frame:
        
            omega = R_i*[(W_cross(3,2)-W_cross(2,3))/2; (-W_cross(3,1)+W_cross(1,3))/2; (W_cross(2,1)-W_cross(1,2))/2];
        
    % Update omega_dot:
    
        % Update omega_dot by rotating the vector to the new reference
        % frame:
        
            omega_dot = R_ip1*R_i'*[x(1); x(2); x(3)];
        
    % Return the state vector x:
    
%     
%         x(7)-q_upd(1)
%         x(8)-q_upd(2)
%         x(9)-q_upd(3)
%         x(10)-q_upd(4)
    
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
        
        
%         pause
        

end

