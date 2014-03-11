function [ rot_mat ] = wingkin_rotation_mat( wing_kin, nr_points )

    % Compute the rotation matrices, angular velocities and angular
    % accelerations:
    
    R_strk      = wing_kin.R_strk;
    
    theta_L     = wing_kin.theta_L;
    eta_L       = wing_kin.eta_L;
    phi_L       = wing_kin.phi_L;
    
    theta_R     = wing_kin.theta_R;
    eta_R       = wing_kin.eta_R;
    phi_R       = wing_kin.phi_R;
    
    theta_dot_L     = wing_kin.theta_dot_L;
    eta_dot_L       = wing_kin.eta_dot_L;
    phi_dot_L       = wing_kin.phi_dot_L;
    
    theta_dot_R     = wing_kin.theta_dot_R;
    eta_dot_R       = wing_kin.eta_dot_R;
    phi_dot_R       = wing_kin.phi_dot_R;
    
    theta_ddot_L     = wing_kin.theta_ddot_L;
    eta_ddot_L       = wing_kin.eta_ddot_L;
    phi_ddot_L       = wing_kin.phi_ddot_L;
    
    theta_ddot_R     = wing_kin.theta_ddot_R;
    eta_ddot_R       = wing_kin.eta_ddot_R;
    phi_ddot_R       = wing_kin.phi_ddot_R;
    
    R_phi_L     = nan(3,3,nr_points);
    R_theta_L   = nan(3,3,nr_points);
    R_eta_L     = nan(3,3,nr_points);
    
    R_phi_R     = nan(3,3,nr_points);
    R_theta_R   = nan(3,3,nr_points);
    R_eta_R     = nan(3,3,nr_points);
    
    RL = nan(3,3,nr_points);
    RR = nan(3,3,nr_points);
    
    wL = nan(3,nr_points);
    wR = nan(3,nr_points);
    
    w_dot_L = nan(3,nr_points);
    w_dot_R = nan(3,nr_points);
    
       
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
        
        wL(:,i) = R_180*R_eta_L(:,:,i)*R_theta_L(:,:,i)*[0; 0; phi_dot_L(i)] + R_180*R_eta_L(:,:,i)*[theta_dot_L(i); 0; 0] + [0; eta_dot_L(i); 0];
        wR(:,i) = R_180*R_eta_R(:,:,i)*R_theta_R(:,:,i)*[0; 0; phi_dot_R(i)] + R_180*R_eta_R(:,:,i)*[theta_dot_R(i); 0; 0] + [0; eta_dot_R(i); 0];
        
        w_dot_L(:,i) = R_180*R_eta_L(:,:,i)*R_theta_L(:,:,i)*[0; 0; phi_ddot_L(i)] + R_180*R_eta_L(:,:,i)*[theta_ddot_L(i); 0; 0] + [0; eta_ddot_L(i); 0];
        w_dot_R(:,i) = R_180*R_eta_R(:,:,i)*R_theta_R(:,:,i)*[0; 0; phi_ddot_R(i)] + R_180*R_eta_R(:,:,i)*[theta_ddot_R(i); 0; 0] + [0; eta_ddot_R(i); 0];

%         RL(:,:,i) = R_eta_L(:,:,i)*R_theta_L(:,:,i)*R_phi_L(:,:,i)*R_strk;
%         
%         RR(:,:,i) = R_eta_R(:,:,i)*R_theta_R(:,:,i)*R_phi_R(:,:,i)*R_strk;
%         
%         wL(:,i) = R_eta_L(:,:,i)*R_theta_L(:,:,i)*[0; 0; phi_dot_L(i)] + R_eta_L(:,:,i)*[theta_dot_L(i); 0; 0] + [0; eta_dot_L(i); 0];
%         
%         wR(:,i) = R_eta_R(:,:,i)*R_theta_R(:,:,i)*[0; 0; phi_dot_R(i)] + R_eta_R(:,:,i)*[theta_dot_R(i); 0; 0] + [0; eta_dot_R(i); 0];
%         
%         w_dot_L(:,i) = R_eta_L(:,:,i)*R_theta_L(:,:,i)*[0; 0; phi_ddot_L(i)] + R_eta_L(:,:,i)*[theta_ddot_L(i); 0; 0] + [0; eta_ddot_L(i); 0];
%         
%         w_dot_R(:,i) = R_eta_R(:,:,i)*R_theta_R(:,:,i)*[0; 0; phi_ddot_R(i)] + R_eta_R(:,:,i)*[theta_ddot_R(i); 0; 0] + [0; eta_ddot_R(i); 0];

    end    
    

    
%     wL(1,:) = wL(1,:) - mean(wL(1,1:(nr_points-1)),2);
%     wR(1,:) = wR(1,:) - mean(wR(1,1:(nr_points-1)),2);
%     wL(2,:) = wL(2,:) - mean(wL(2,1:(nr_points-1)),2);
%     wR(2,:) = wR(2,:) - mean(wR(2,1:(nr_points-1)),2);
%     wL(3,:) = wL(3,:) - mean(wL(3,1:(nr_points-1)),2);
%     wR(3,:) = wR(3,:) - mean(wR(3,1:(nr_points-1)),2);
%     w_dot_L(1,:) = w_dot_L(1,:) - mean(w_dot_L(1,1:(nr_points-1)),2);
%     w_dot_R(1,:) = w_dot_R(1,:) - mean(w_dot_R(1,1:(nr_points-1)),2);
%     w_dot_L(2,:) = w_dot_L(2,:) - mean(w_dot_L(2,1:(nr_points-1)),2);
%     w_dot_R(2,:) = w_dot_R(2,:) - mean(w_dot_R(2,1:(nr_points-1)),2);
%     w_dot_L(3,:) = w_dot_L(3,:) - mean(w_dot_L(3,1:(nr_points-1)),2);
%     w_dot_R(3,:) = w_dot_R(3,:) - mean(w_dot_R(3,1:(nr_points-1)),2);
    
%     figure()
%     plot(1:nr_points,wL(1,:),'r',1:nr_points,wR(1,:),'b',1:nr_points,wL(1,:)+wR(1,:),'k')
%     
%     figure()
%     plot(1:nr_points,wL(2,:),'r',1:nr_points,wR(2,:),'b',1:nr_points,mean(wL(2,:)),'k')
%     
%     figure()
%     plot(1:nr_points,wL(3,:),'r',1:nr_points,wR(3,:),'b',1:nr_points,wL(3,:)+wR(3,:),'k')
%     
%         
%     wL_strk = nan(3,nr_points);
%     wR_strk = nan(3,nr_points);
%     
%     for i = 1:nr_points
%         
%         wL_strk(:,i) = R_strk*RL(:,:,i)'*wL(:,i);
%         wR_strk(:,i) = R_strk*RR(:,:,i)'*wR(:,i);
%         
%     end
    
%     figure()
%     plot(1:nr_points,wL_strk(1,:),'r',1:nr_points,wR_strk(1,:),'b',1:nr_points,wL_strk(1,:)+wR_strk(1,:),'k')
%     
%     figure()
%     plot(1:nr_points,wL_strk(2,:),'r',1:nr_points,wR_strk(2,:),'b',1:nr_points,mean(wL_strk(2,:)),'k')
%     
%     figure()
%     plot(1:nr_points,wL_strk(3,:),'r',1:nr_points,wR_strk(3,:),'b',1:nr_points,wL_strk(3,:)+wR_strk(3,:),'k')
%     
%     wL_strk(1,:) = wL_strk(1,:)-mean(wL_strk(1,1:(nr_points)));
%     wL_strk(2,:) = wL_strk(2,:)-mean(wL_strk(2,1:(nr_points)));
%     wL_strk(3,:) = wL_strk(3,:)-mean(wL_strk(3,1:(nr_points)));
%     
%     wR_strk(1,:) = wR_strk(1,:)-mean(wR_strk(1,1:(nr_points)));
%     wR_strk(2,:) = wR_strk(2,:)-mean(wR_strk(2,1:(nr_points)));
%     wR_strk(3,:) = wR_strk(3,:)-mean(wR_strk(3,1:(nr_points)));
%     
%     for i = 1:nr_points
%         
%         wL(:,i) = RL(:,:,i)*R_strk'*wL_strk(:,i);
%         wR(:,i) = RR(:,:,i)*R_strk'*wR_strk(:,i);
%         
%     end
    
%     figure()
%     plot(1:nr_points,wL_strk(1,:),'r',1:nr_points,wR_strk(1,:),'b',1:nr_points,wL_strk(1,:)+wR_strk(1,:),'k')
%     
%     figure()
%     plot(1:nr_points,wL_strk(2,:),'r',1:nr_points,wR_strk(2,:),'b',1:nr_points,mean(wL_strk(2,:)),'k')
%     
%     figure()
%     plot(1:nr_points,wL_strk(3,:),'r',1:nr_points,wR_strk(3,:),'b',1:nr_points,wL_strk(3,:)+wR_strk(3,:),'k')
%     
%     figure()
%     plot(1:nr_points,wL(1,:),'r',1:nr_points,wR(1,:),'b',1:nr_points,wL(1,:)+wR(1,:),'k')
%     
%     figure()
%     plot(1:nr_points,wL(2,:),'r',1:nr_points,wR(2,:),'b',1:nr_points,mean(wL(2,:)),'k')
%     
%     figure()
%     plot(1:nr_points,wL(3,:),'r',1:nr_points,wR(3,:),'b',1:nr_points,wL(3,:)+wR(3,:),'k')
%     
%     pause
    
    rot_mat.R_theta_L       = R_theta_L;
    rot_mat.R_eta_L         = R_eta_L;
    rot_mat.R_phi_L         = R_phi_L;
    rot_mat.R_theta_R       = R_theta_R;
    rot_mat.R_eta_R         = R_eta_R;
    rot_mat.R_phi_R         = R_phi_R;
    rot_mat.RL              = RL;
    rot_mat.RR              = RR;
    rot_mat.wL              = wL;
    rot_mat.wR              = wR;
    rot_mat.w_dot_L         = w_dot_L;
    rot_mat.w_dot_R         = w_dot_R;

end

