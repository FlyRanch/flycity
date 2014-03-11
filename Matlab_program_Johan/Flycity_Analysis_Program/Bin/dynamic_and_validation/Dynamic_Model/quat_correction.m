function [ qb_new ] = quat_correction( Fg_b, qb )

    % Correct the body quaternion for the direction of the gravity vector
    % in the inertial reference frame:
    
    Rb = quat2mat(qb);
    
    Fg = Rb'*Fg_b;
    
    Fg_norm = Fg/norm(Fg);
    
    e_xyz = cross(Fg_norm,[0;0; -1])/norm(cross(Fg_norm,[0;0; -1]));
    
    theta = asin(norm(cross(Fg_norm,[0;0; -1]))/(norm(Fg_norm)*1));
    
    q_err_t = [e_xyz*sin(theta/2); cos(theta/2)];
    
    q_err = q_err_t/norm(q_err_t);
    
    Q = [ q_err(4) -q_err(3) q_err(2) q_err(1); ...
          q_err(3) q_err(4) -q_err(1) q_err(2); ...
          -q_err(2) q_err(1) q_err(4) q_err(3); ...
          -q_err(1) -q_err(2) -q_err(3) q_err(4)];
      
    qb_t = Q*qb;
    
    qb_new = qb_t/norm(qb_t);
    
end

