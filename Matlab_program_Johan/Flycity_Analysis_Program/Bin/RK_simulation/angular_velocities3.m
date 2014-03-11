function [ kine ] = angular_velocities3( a_fit )


    % Compute the wing kinematics:
    
    theta_L0    = a_fit.theta_L0;
    eta_L0      = a_fit.eta_L0;
    phi_L0      = a_fit.phi_L0;
    theta_R0    = a_fit.theta_R0;
    eta_R0      = a_fit.eta_R0;
    phi_R0      = a_fit.phi_R0;
    A_theta_L   = a_fit.A_theta_L;
    A_eta_L     = a_fit.A_eta_L;
    A_phi_L     = a_fit.A_phi_L;
    A_theta_R   = a_fit.A_theta_R;
    A_eta_R     = a_fit.A_eta_R;
    A_phi_R     = a_fit.A_phi_R;
    f           = a_fit.f;
    nr_points   = a_fit.nr_points;
    R_strk      = a_fit.R_strk;
    
    dt = (1/f)/(nr_points-1);
    
    t = 0:dt:(dt*(nr_points-1));
    
    theta_L         = theta_L0+A_theta_L*cos(4*pi*f*t);
    eta_L           = eta_L0+A_eta_L*sin(2*pi*f*t);
    phi_L           = phi_L0+A_phi_L*cos(2*pi*f*t);
    
    theta_R         = theta_R0+A_theta_R*cos(4*pi*f*t);
    eta_R           = eta_R0+A_eta_R*sin(2*pi*f*t);
    phi_R           = phi_R0+A_phi_R*cos(2*pi*f*t);
    
    theta_dot_L     = -A_theta_L*4*pi*f*sin(4*pi*f*t);
    eta_dot_L       = A_eta_L*2*pi*f*cos(2*pi*f*t);
    phi_dot_L       = A_phi_L*2*pi*f*sin(2*pi*f*t);
    
    theta_dot_R     = A_theta_R*4*pi*f*sin(4*pi*f*t);
    eta_dot_R       = A_eta_R*2*pi*f*cos(2*pi*f*t);
    phi_dot_R       = -A_phi_R*2*pi*f*sin(2*pi*f*t);
    
    theta_ddot_L    = -A_theta_L*16*pi^2*f^2*cos(4*pi*f*t);
    eta_ddot_L      = -A_eta_L*4*pi^2*f^2*sin(2*pi*f*t);
    phi_ddot_L      = A_phi_L*4*pi^2*f^2*cos(2*pi*f*t);
    
    theta_ddot_R    = A_theta_R*16*pi^2*f^2*cos(4*pi*f*t);
    eta_ddot_R      = -A_eta_R*4*pi^2*f^2*sin(2*pi*f*t);
    phi_ddot_R      = -A_phi_R*4*pi^2*f^2*cos(2*pi*f*t);
    
    R_phi_L     = nan(3,3,nr_points);
    R_theta_L   = nan(3,3,nr_points);
    R_eta_L     = nan(3,3,nr_points);
    
    R_phi_R     = nan(3,3,nr_points);
    R_theta_R   = nan(3,3,nr_points);
    R_eta_R     = nan(3,3,nr_points);
    
    RL = nan(3,3,nr_points);
    RR = nan(3,3,nr_points);
    
    dR_dt_L = nan(3,3,points);
    dR_dt_R = nan(3,3,points);
    
     for i = 1:nr_points
        
        % Matrix to switch orientation wing:
        
        R_180              = [ -1  0  0; ...
                                0  1  0; ...
                                0  0 -1];
                           
        % Wing kinematics rotation matrices:
        
        R_theta_L   = [ 1  0               0              ; ...
                        0  cos(theta_L(i)) sin(theta_L(i)); ...
                        0 -sin(theta_L(i)) cos(theta_L(i))];

        R_eta_L     = [  cos(eta_L(i)) 0 sin(eta_L(i)); ...
                         0             1 0            ; ...
                        -sin(eta_L(i)) 0 cos(eta_L(i))];
    
        R_phi_R     = [ cos(phi_R(i)) sin(phi_R(i)) 0; ...
                       -sin(phi_R(i)) cos(phi_R(i)) 0; ...
                        0             0             1];

        R_theta_R   = [ 1 0                0              ; ...
                        0 cos(theta_R(i)) -sin(theta_R(i)); ...
                        0 sin(theta_R(i))  cos(theta_R(i))];

        R_eta_R     = [  cos(eta_R(i)) 0 sin(eta_R(i)); ...
                         0             1 0            ; ...
                        -sin(eta_R(i)) 0 cos(eta_R(i))];
                    
        R_phi_L     = [ cos(phi_L(i)) -sin(phi_L(i)) 0; ...
                        sin(phi_L(i))  cos(phi_L(i)) 0; ...
                        0              0             1];
                    
        RL(:,:,i) = R_180*R_eta_L*R_theta_L*R_phi_L*R_strk;
        RR(:,:,i) = R_180*R_eta_R*R_theta_R*R_phi_R*R_strk;
        
        dR_dt_L(:,:,i) = 
        
     end
     
     % Convert the rotation matrices to quaternions:
     
     qL = zeros(4,nr_points);
     qR = zeros(4,nr_points);
     
     for i = 1:nr_points
         
         qL_t = quat2mat(RL(:,:,i));
         
         qL_t = qL_t/norm(qL_t);
         
         qR_t = quat2mat(RR(:,:,i));
         
         qR_t = qR_t/norm(qR_t);
         
         qL(:,i) = qL_t;
         qR(:,i) = qR_t;
         
         % quaternion rate:
         
         
         
     end
     
     figure()
     hold on
     subplot(4,1,1); plot(t,qL(1,:),t,qR(1,:))
     subplot(4,1,2); plot(t,qL(2,:),t,qR(2,:))
     subplot(4,1,3); plot(t,qL(3,:),t,qR(3,:))
     subplot(4,1,4); plot(t,qL(4,:),t,qR(4,:))
     hold off
    
    pause
     
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
        
        R_phi_L(:,:,i)     = [ cos(phi_L(i)) -sin(phi_L(i)) 0; ...
                               sin(phi_L(i))  cos(phi_L(i)) 0; ...
                               0              0             1];

        R_theta_L(:,:,i)   = [ 1  0               0              ; ...
                               0  cos(theta_L(i)) sin(theta_L(i)); ...
                               0 -sin(theta_L(i)) cos(theta_L(i))];

        R_eta_L(:,:,i)     = [  cos(eta_L(i)) 0 sin(eta_L(i)); ...
                                0             1 0            ; ...
                               -sin(eta_L(i)) 0 cos(eta_L(i))];
    
        R_phi_R(:,:,i)     = [ cos(phi_R(i)) sin(phi_R(i)) 0; ...
                              -sin(phi_R(i)) cos(phi_R(i)) 0; ...
                               0             0             1];

        R_theta_R(:,:,i)   = [ 1 0                0              ; ...
                               0 cos(theta_R(i)) -sin(theta_R(i)); ...
                               0 sin(theta_R(i))  cos(theta_R(i))];

        R_eta_R(:,:,i)     = [  cos(eta_R(i)) 0 sin(eta_R(i)); ...
                                0             1 0            ; ...
                               -sin(eta_R(i)) 0 cos(eta_R(i))];
                    
        RL(:,:,i) = R_180*R_eta_L(:,:,i)*R_theta_L(:,:,i)*R_phi_L(:,:,i)*R_strk;
        RR(:,:,i) = R_180*R_eta_R(:,:,i)*R_theta_R(:,:,i)*R_phi_R(:,:,i)*R_strk;
        
%         RL(:,:,i) = RL(:,:,i)./norm(RL(:,:,i));
%         RR(:,:,i) = RR(:,:,i)./norm(RR(:,:,i));

        wL(:,i) = R_180*R_eta_L(:,:,i)*R_theta_L(:,:,i)*[0; 0; phi_dot_L(i)] + R_180*R_eta_L(:,:,i)*[theta_dot_L(i); 0; 0] + R_180*[0; eta_dot_L(i); 0];
        wR(:,i) = R_180*R_eta_R(:,:,i)*R_theta_R(:,:,i)*[0; 0; phi_dot_R(i)] + R_180*R_eta_R(:,:,i)*[theta_dot_R(i); 0; 0] + R_180*[0; eta_dot_R(i); 0];
           
        w_dot_L(:,i) = R_180*R_eta_L(:,:,i)*R_theta_L(:,:,i)*[0; 0; phi_ddot_L(i)] + R_180*R_eta_L(:,:,i)*[theta_ddot_L(i); 0; 0] + R_180*[0; eta_ddot_L(i); 0];
        w_dot_R(:,i) = R_180*R_eta_R(:,:,i)*R_theta_R(:,:,i)*[0; 0; phi_ddot_R(i)] + R_180*R_eta_R(:,:,i)*[theta_ddot_R(i); 0; 0] + R_180*[0; eta_ddot_R(i); 0];

        
        wL_b(:,i) = RL(:,:,i)'*wL(:,i);
        wR_b(:,i) = RR(:,:,i)'*wR(:,i);

        w_dot_L_b(:,i) = RL(:,:,i)'*w_dot_L(:,i);
        w_dot_R_b(:,i) = RR(:,:,i)'*w_dot_R(:,i);
        
        wL_strk(:,i) = R_strk*wL_b(:,i);
        wR_strk(:,i) = R_strk*wR_b(:,i);
        
        w_dot_L_strk(:,i) = R_strk*w_dot_L_b(:,i);
        w_dot_R_strk(:,i) = R_strk*w_dot_R_b(:,i);
        
    end
    
    wL_b_0 = [mean(wL_b(1,1:(end-1))); mean(wL_b(2,1:(end-1))); mean(wL_b(3,1:(end-1)))];
    wR_b_0 = [mean(wR_b(1,1:(end-1))); mean(wR_b(2,1:(end-1))); mean(wR_b(3,1:(end-1)))];

    
%     wL_b(1,:) = wL_b(1,:)-mean(wL_b(1,1:(end-1)));
%     wL_b(2,:) = wL_b(2,:)-mean(wL_b(2,1:(end-1)));
%     wL_b(3,:) = wL_b(3,:)-mean(wL_b(3,1:(end-1)));
%     wR_b(1,:) = wR_b(1,:)-mean(wR_b(1,1:(end-1)));
%     wR_b(2,:) = wR_b(2,:)-mean(wR_b(2,1:(end-1)));
%     wR_b(3,:) = wR_b(3,:)-mean(wR_b(3,1:(end-1)));
%     w_dot_L_b(1,:) = w_dot_L_b(1,:)-mean(w_dot_L_b(1,1:(end-1)));
%     w_dot_L_b(2,:) = w_dot_L_b(2,:)-mean(w_dot_L_b(2,1:(end-1)));
%     w_dot_L_b(3,:) = w_dot_L_b(3,:)-mean(w_dot_L_b(3,1:(end-1)));
%     w_dot_R_b(1,:) = w_dot_R_b(1,:)-mean(w_dot_R_b(1,1:(end-1)));
%     w_dot_R_b(2,:) = w_dot_R_b(2,:)-mean(w_dot_R_b(2,1:(end-1)));
%     w_dot_R_b(3,:) = w_dot_R_b(3,:)-mean(w_dot_R_b(3,1:(end-1)));
    
    
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
    kine.wL_b_0         = wL_b_0;
    kine.wR_b_0         = wR_b_0;


end