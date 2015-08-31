function [ kine ] = angular_velocities_polynomial2( theta_L_t, eta_L_t, phi_L_t, theta_R_t, eta_R_t, phi_R_t, f, R_strk )


    % Compute the wing kinematics:
    
    nr_points   = length(theta_L_t);
    
    dt = (7/f)/(nr_points-1);  
    
    t = 0:dt:nr_points*dt;
    
    theta_L         = theta_L_t;
    eta_L           = -eta_L_t;
    phi_L           = -phi_L_t;
    
    theta_R         = -theta_R_t;
    eta_R           = -eta_R_t;
    phi_R           = phi_R_t;

    
    theta_dot_L     = gradient(theta_L)./dt;
    eta_dot_L       = gradient(eta_L)./dt;
    phi_dot_L       = gradient(phi_L)./dt;
    
    theta_dot_R     = gradient(theta_R)./dt;
    eta_dot_R       = gradient(eta_R)./dt;
    phi_dot_R       = gradient(phi_R)./dt;
    
    theta_ddot_L    = gradient(theta_dot_L)./dt;
    eta_ddot_L      = gradient(eta_dot_L)./dt;
    phi_ddot_L      = gradient(phi_dot_L)./dt;
    
    theta_ddot_R    = gradient(theta_dot_R)./dt;
    eta_ddot_R      = gradient(eta_dot_R)./dt;
    phi_ddot_R      = gradient(phi_dot_R)./dt;
    
    R_phi_L         = nan(3,3,nr_points);
    R_theta_L       = nan(3,3,nr_points);
    R_eta_L         = nan(3,3,nr_points);
    
    R_phi_R         = nan(3,3,nr_points);
    R_theta_R       = nan(3,3,nr_points);
    R_eta_R         = nan(3,3,nr_points);
    
    R_dot_theta_L   = nan(3,3,nr_points);
    R_dot_eta_L     = nan(3,3,nr_points);
    
    R_dot_theta_R   = nan(3,3,nr_points);
    R_dot_eta_R     = nan(3,3,nr_points);
    
    RL = nan(3,3,nr_points);
    RR = nan(3,3,nr_points);
    
    qL = nan(4,nr_points);
    qR = nan(4,nr_points);
    
    wL = nan(3,nr_points);
    wR = nan(3,nr_points);
    
    w_dot_L = nan(3,nr_points);
    w_dot_R = nan(3,nr_points);
    
    wL_b = nan(3,nr_points);
    wR_b = nan(3,nr_points);
    
    w_dot_L_b = nan(3,nr_points);
    w_dot_R_b = nan(3,nr_points);
    
    wL_strk = nan(3,nr_points);
    wR_strk = nan(3,nr_points);
    
    w_dot_L_strk = nan(3,nr_points);
    w_dot_R_strk = nan(3,nr_points);
       
    for i = 1:nr_points
        
        % Matrix to switch orientation wing:
        
        R_180              = [ -1  0  0; ...
                                0  1  0; ...
                                0  0 -1];
                           
        % Wing kinematics rotation matrices:
        
        R_phi_L(:,:,i)     = [ cos(phi_L(i))  sin(phi_L(i)) 0; ...
                              -sin(phi_L(i))  cos(phi_L(i)) 0; ...
                               0              0             1];

        R_theta_L(:,:,i)   = [ 1  0               0              ; ...
                               0  cos(theta_L(i)) sin(theta_L(i)); ...
                               0 -sin(theta_L(i)) cos(theta_L(i))];

        R_eta_L(:,:,i)     = [  cos(eta_L(i)) 0 -sin(eta_L(i)); ...
                                0             1  0            ; ...
                                sin(eta_L(i)) 0  cos(eta_L(i))];
    
        R_phi_R(:,:,i)     = [ cos(phi_R(i))  sin(phi_R(i)) 0; ...
                              -sin(phi_R(i))  cos(phi_R(i)) 0; ...
                               0              0             1];

        R_theta_R(:,:,i)   = [ 1  0               0              ; ...
                               0  cos(theta_R(i)) sin(theta_R(i)); ...
                               0 -sin(theta_R(i)) cos(theta_R(i))];

        R_eta_R(:,:,i)     = [  cos(eta_R(i)) 0 -sin(eta_R(i)); ...
                                0             1  0            ; ...
                                sin(eta_R(i)) 0  cos(eta_R(i))];
                           
        R_dot_eta_L(:,:,i) = [ -eta_dot_L(i)*sin(eta_L(i)) 0 -eta_dot_L(i)*cos(eta_L(i)); ...
                                0                          0  0                         ; ...
                                eta_dot_L(i)*cos(eta_L(i)) 0 -eta_dot_L(i)*sin(eta_L(i))];
                           
        R_dot_theta_L(:,:,i) = [ 0  0                               0                             ; ...
                                 0 -theta_dot_L(i)*sin(theta_L(i))  theta_dot_L(i)*cos(theta_L(i)); ...
                                 0 -theta_dot_L(i)*cos(theta_L(i)) -theta_dot_L(i)*sin(theta_L(i))];
                           
        R_dot_eta_R(:,:,i) = [ -eta_dot_R(i)*sin(eta_R(i)) 0 -eta_dot_R(i)*cos(eta_R(i)); ...
                                0                          0  0                         ; ...
                                eta_dot_R(i)*cos(eta_R(i)) 0 -eta_dot_R(i)*sin(eta_R(i))];
                           
        R_dot_theta_R(:,:,i) = [ 0  0                               0                             ; ...
                                 0 -theta_dot_R(i)*sin(theta_R(i))  theta_dot_R(i)*cos(theta_R(i)); ...
                                 0 -theta_dot_R(i)*cos(theta_R(i)) -theta_dot_R(i)*sin(theta_R(i))];

        RL(:,:,i)    = R_180*R_eta_L(:,:,i)*R_theta_L(:,:,i)*R_phi_L(:,:,i)*R_strk;
        RR(:,:,i)    = R_180*R_eta_R(:,:,i)*R_theta_R(:,:,i)*R_phi_R(:,:,i)*R_strk;
        
        RL(:,:,i)    = RL(:,:,i)./norm(RL(:,:,i));
        RR(:,:,i)    = RR(:,:,i)./norm(RR(:,:,i));
        
        qL_t         = quat2mat(RL(:,:,i));
        qR_t         = quat2mat(RR(:,:,i));
        
        qL(:,i)      = qL_t/norm(qL_t);
        qR(:,i)      = qR_t/norm(qR_t);

        wL(:,i)      = R_180*R_eta_L(:,:,i)*R_theta_L(:,:,i)*[0; 0; phi_dot_L(i)] + R_180*R_eta_L(:,:,i)*[theta_dot_L(i); 0; 0] + R_180*[0; eta_dot_L(i); 0];
        wR(:,i)      = R_180*R_eta_R(:,:,i)*R_theta_R(:,:,i)*[0; 0; phi_dot_R(i)] + R_180*R_eta_R(:,:,i)*[theta_dot_R(i); 0; 0] + R_180*[0; eta_dot_R(i); 0];
        
        w_dot_L(:,i) = R_180*R_dot_eta_L(:,:,i)*R_theta_L(:,:,i)*[0; 0; phi_dot_L(i)] + ...
                       R_180*R_eta_L(:,:,i)*R_dot_theta_L(:,:,i)*[0; 0; phi_dot_L(i)] + ...
                       R_180*R_eta_L(:,:,i)*R_theta_L(:,:,i)*[0; 0; phi_ddot_L(i)] + ...
                       R_180*R_dot_eta_L(:,:,i)*[theta_dot_L(i); 0; 0] + ...
                       R_180*R_eta_L(:,:,i)*[theta_ddot_L(i); 0; 0] + ...
                       R_180*[0; eta_ddot_L(i); 0];
        
        w_dot_R(:,i) = R_180*R_dot_eta_R(:,:,i)*R_theta_R(:,:,i)*[0; 0; phi_dot_R(i)] + ...
                       R_180*R_eta_R(:,:,i)*R_dot_theta_R(:,:,i)*[0; 0; phi_dot_R(i)] + ...
                       R_180*R_eta_R(:,:,i)*R_theta_R(:,:,i)*[0; 0; phi_ddot_R(i)] + ...
                       R_180*R_dot_eta_R(:,:,i)*[theta_dot_R(i); 0; 0] + ...
                       R_180*R_eta_R(:,:,i)*[theta_ddot_R(i); 0; 0] + ...
                       R_180*[0; eta_ddot_R(i); 0];
                   
        wL_b(:,i) = RL(:,:,i)'*wL(:,i);
        wR_b(:,i) = RR(:,:,i)'*wR(:,i);
        
        w_dot_L_b(:,i) = RL(:,:,i)'*w_dot_L(:,i);
        w_dot_R_b(:,i) = RR(:,:,i)'*w_dot_R(:,i);
        
        wL_strk(:,i) = R_strk*wL_b(:,i);
        wR_strk(:,i) = R_strk*wR_b(:,i);
        
        w_dot_L_strk(:,i) = R_strk*w_dot_L_b(:,i);
        w_dot_R_strk(:,i) = R_strk*w_dot_R_b(:,i);
        
    end


    kine.dt             = dt;
    kine.t              = t;
    kine.theta_L        = theta_L;
    kine.eta_L          = eta_L;
    kine.phi_L          = phi_L;
    kine.theta_R        = theta_R;
    kine.eta_R          = eta_R;
    kine.phi_R          = phi_R;
    kine.theta_dot_L    = theta_dot_L;
    kine.eta_dot_L      = eta_dot_L;
    kine.phi_dot_L      = phi_dot_L;
    kine.theta_dot_R    = theta_dot_R;
    kine.eta_dot_R      = eta_dot_R;
    kine.phi_dot_R      = phi_dot_R;
    kine.theta_ddot_L   = theta_ddot_L;
    kine.eta_ddot_L     = eta_ddot_L;
    kine.phi_ddot_L     = phi_ddot_L;
    kine.theta_ddot_R   = theta_ddot_R;
    kine.eta_ddot_R     = eta_ddot_R;
    kine.phi_ddot_R     = phi_ddot_R;
    kine.RL             = RL;
    kine.RR             = RR;
    kine.qL             = qL;
    kine.qR             = qR;
    kine.wL             = wL;
    kine.wR             = wR;
    kine.w_dot_L        = w_dot_L;
    kine.w_dot_R        = w_dot_R;
    kine.wL_b           = wL_b;
    kine.wR_b           = wR_b;
    kine.w_dot_L_b      = w_dot_L_b;
    kine.w_dot_R_b      = w_dot_R_b;
    kine.wL_strk        = wL_strk;
    kine.wR_strk        = wR_strk;
    kine.w_dot_L_strk   = w_dot_L_strk;
    kine.w_dot_R_strk   = w_dot_R_strk;

    end