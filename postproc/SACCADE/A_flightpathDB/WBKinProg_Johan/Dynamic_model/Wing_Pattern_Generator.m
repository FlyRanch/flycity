function [q_L,w_L,w_dot_L,q_R,w_R,w_dot_R,up,down] = Wing_Pattern_Generator(Pattern,t, dt)


    % Return the wing kinematics for a certain value of time t:
    
    strokeplane_WBkin = settings.strokeplane_WBkin;

    
    frame_nr = round(t/(0.5*dt))+1;
    
    
    phi_L = Pattern.phi(frame_nr,1);
    
    theta_L = Pattern.theta(frame_nr,1);
    
    eta_L = pi-Pattern.eta(frame_nr,1);

    
    phi_R = Pattern.phi(frame_nr,2);
    
    theta_R = Pattern.theta(frame_nr,2);
    
    eta_R = pi-Pattern.eta(frame_nr,2);
    
    
    phi_L_dot = Pattern.phi_dot(frame_nr,1);
    
    theta_L_dot = Pattern.theta_dot(frame_nr,1);
    
    eta_L_dot = -Pattern.eta_dot(frame_nr,1);

    
    phi_R_dot = Pattern.phi_dot(frame_nr,2);
    
    theta_R_dot = Pattern.theta_dot(frame_nr,2);
    
    eta_R_dot = -Pattern.eta_dot(frame_nr,2);
        
    
    phi_L_ddot = Pattern.phi_ddot(frame_nr,1);
    
    theta_L_ddot = Pattern.theta_ddot(frame_nr,1);
    
    eta_L_ddot = -Pattern.eta_ddot(frame_nr,1);

    
    phi_R_ddot = Pattern.phi_ddot(frame_nr,2);
    
    theta_R_ddot = Pattern.theta_ddot(frame_nr,2);
    
    eta_R_ddot = -Pattern.eta_ddot(frame_nr,2);
    
    
    up = Pattern.up(frame_nr);
    
    down = Pattern.down(frame_nr);
    
    
    
%     beta = -(55/180)*pi;    
    beta = (strokeplane_WBkin/180)*pi;
    
    R_phi_L =  [cos(-phi_L) sin(-phi_L) 0; ...
                -sin(-phi_L) cos(-phi_L) 0; ...
                0 0 1];
    
    R_theta_L = [1 0 0; ...
                 0 cos(theta_L) sin(theta_L); ...
                 0 -sin(theta_L) cos(theta_L)];
    
    R_eta_L = [cos(eta_L) 0 -sin(eta_L); ...
               0 1 0; ...
               sin(eta_L) 0 cos(eta_L)];

    R_phi_R =  [cos(phi_R) sin(phi_R) 0; ...
                -sin(phi_R) cos(phi_R) 0; ...
                0 0 1];
           
    R_theta_R = [1 0 0; ...
                 0 cos(-theta_R) sin(-theta_R); ...
                 0 -sin(-theta_R) cos(-theta_R)];
    
    R_eta_R = [cos(eta_R) 0 -sin(eta_R); ...
               0 1 0; ...
               sin(eta_R) 0 cos(eta_R)];    
           
    R_beta = [cos(beta) 0 -sin(beta); ...
               0 1 0; ...
               sin(beta) 0 cos(beta)]; 
           
    R_rot_L = R_eta_L*R_theta_L*R_phi_L*R_beta;
           
    R_rot_R = R_eta_R*R_theta_R*R_phi_R*R_beta;  
           
    q_L = quat2matNEW(R_rot_L');
    
    q_R = quat2matNEW(R_rot_R');
    


    w_L = ([0; eta_L_dot; 0] + R_eta_L*[theta_L_dot; 0; 0] + R_eta_L*R_theta_L*[0; 0; -phi_L_dot]);
    
    w_R = ([0; eta_R_dot; 0] + R_eta_R*[-theta_R_dot; 0; 0] + R_eta_R*R_theta_R*[0; 0; phi_R_dot]);
    
    
   
    
    % Compute w_L and w_R:
    
    w_dot_L = ([0; eta_L_ddot; 0] + R_eta_L*[theta_L_ddot; 0; 0] + R_eta_L*R_theta_L*[0; 0; -phi_L_ddot]);
    
    w_dot_R = ([0; eta_R_ddot; 0] + R_eta_R*[-theta_R_ddot; 0; 0] + R_eta_R*R_theta_R*[0; 0; phi_R_ddot]);

end

