function [q_L,w_L,w_dot_L,q_R,w_R,w_dot_R,Rot_L,Rot_R] = Wingkin_2_quaternion( wing_kin )


    % Function which transforms wing kinematic angles into quaternions,
    % angular velocities and angular accelerations

    theta_L = wing_kin.theta_L;
    eta_L = wing_kin.eta_L;
    phi_L = wing_kin.phi_L;
    
    theta_R = wing_kin.theta_R;
    eta_R = wing_kin.eta_R;
    phi_R = wing_kin.phi_R;
    
    theta_dot_L = wing_kin.theta_dot_L;
    eta_dot_L = wing_kin.eta_dot_L;
    phi_dot_L = wing_kin.phi_dot_L;
    
    theta_dot_R = wing_kin.theta_dot_R;
    eta_dot_R = wing_kin.eta_dot_R;
    phi_dot_R = wing_kin.phi_dot_R;
    
    theta_ddot_L = wing_kin.theta_ddot_L;
    eta_ddot_L = wing_kin.eta_ddot_L;
    phi_ddot_L = wing_kin.phi_ddot_L;
    
    theta_ddot_R = wing_kin.theta_ddot_R;
    eta_ddot_R = wing_kin.eta_ddot_R;
    phi_ddot_R = wing_kin.phi_ddot_R;
    
    n = length(theta_L);
    
    q_L = zeros(4,n);
    q_R = zeros(4,n);
    
    w_L = zeros(3,n);
    w_R = zeros(3,n);
    
    w_dot_L = zeros(3,n);
    w_dot_R = zeros(3,n);
    
    Rot_L = zeros(3,3,n);
    
    Rot_R = zeros(3,3,n);
    
    beta = (55/180)*pi;    
    
    for i = 1:n
    
        R_phi_L =  [cos(-phi_L(i)) sin(-phi_L(i)) 0; ...
                    -sin(-phi_L(i)) cos(-phi_L(i)) 0; ...
                    0 0 1];

        R_theta_L = [1 0 0; ...
                     0 cos(theta_L(i)) sin(theta_L(i)); ...
                     0 -sin(theta_L(i)) cos(theta_L(i))];

        R_eta_L = [cos(pi-eta_L(i)) 0 -sin(pi-eta_L(i)); ...
                   0 1 0; ...
                   sin(pi-eta_L(i)) 0 cos(pi-eta_L(i))];

        R_phi_R =  [cos(phi_R(i)) sin(phi_R(i)) 0; ...
                    -sin(phi_R(i)) cos(phi_R(i)) 0; ...
                    0 0 1];

        R_theta_R = [1 0 0; ...
                     0 cos(-theta_R(i)) sin(-theta_R(i)); ...
                     0 -sin(-theta_R(i)) cos(-theta_R(i))];

        R_eta_R = [cos(pi-eta_R(i)) 0 -sin(pi-eta_R(i)); ...
                   0 1 0; ...
                   sin(pi-eta_R(i)) 0 cos(pi-eta_R(i))];    

        R_beta = [cos(-beta) 0 -sin(-beta); ...
                   0 1 0; ...
                   sin(-beta) 0 cos(-beta)]; 

        R_rot_L = R_eta_L*R_theta_L*R_phi_L*R_beta;

        R_rot_R = R_eta_R*R_theta_R*R_phi_R*R_beta;  

        q_L(:,i) = quat2matNEW(R_rot_L');

        q_R(:,i) = quat2matNEW(R_rot_R');
        
        w_L(:,i) = ([0; -eta_dot_L(i); 0] + R_eta_L*[theta_dot_L(i); 0; 0] + R_eta_L*R_theta_L*[0; 0; -phi_dot_L(i)]);
    
        w_R(:,i) = ([0; -eta_dot_R(i); 0] + R_eta_R*[-theta_dot_R(i); 0; 0] + R_eta_R*R_theta_R*[0; 0; phi_dot_R(i)]);

        w_dot_L(:,i) = ([0; -eta_ddot_L(i); 0] + R_eta_L*[theta_ddot_L(i); 0; 0] + R_eta_L*R_theta_L*[0; 0; -phi_ddot_L(i)]);
    
        w_dot_R(:,i) = ([0; -eta_ddot_R(i); 0] + R_eta_R*[-theta_ddot_R(i); 0; 0] + R_eta_R*R_theta_R*[0; 0; phi_ddot_R(i)]);
        
        Rot_L(:,:,i) = R_rot_L;
        
        Rot_R(:,:,i) = R_rot_R;
    
    end
    
end

