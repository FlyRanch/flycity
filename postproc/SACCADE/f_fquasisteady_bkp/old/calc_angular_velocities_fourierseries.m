function [ kine ] = calc_angular_velocities_fourierseries(eta_L,phi_L,theta_L,eta_R,phi_R,theta_R,f,R_strk)


    % Compute the wing kinematics:
    nr_points = length(eta_L);
    dt = (1/f)/(nr_points-1);
    t = 0:dt:(dt*(nr_points-1));  

    % diff data
    theta_dot_L     = diff(theta_L)./dt;
    eta_dot_L     = diff(eta_L)./dt;
    phi_dot_L     = diff(phi_L)./dt;
    
    theta_dot_R     = diff(theta_R)./dt;
    eta_dot_R     = diff(eta_R)./dt;
    phi_dot_R     = diff(phi_R)./dt;

    theta_ddot_L     = diff(theta_dot_L)./dt;
    eta_ddot_L     = diff(eta_dot_L)./dt;
    phi_ddot_L     = diff(phi_dot_L)./dt;
    
    theta_ddot_R     = diff(theta_dot_R)./dt;
    eta_ddot_R     = diff(eta_dot_R)./dt;
    phi_ddot_R     = diff(phi_dot_R)./dt;

    % calc R&w's
    R_phi_L     = nan(3,3,nr_points);
    R_theta_L   = nan(3,3,nr_points);
    R_eta_L     = nan(3,3,nr_points);
    
    R_phi_R     = nan(3,3,nr_points);
    R_theta_R   = nan(3,3,nr_points);
    R_eta_R     = nan(3,3,nr_points);
    
    R_dot_theta_L   = nan(3,3,nr_points);
    R_dot_eta_L     = nan(3,3,nr_points);
    
    R_dot_theta_R   = nan(3,3,nr_points);
    R_dot_eta_R     = nan(3,3,nr_points);
    
    RL = nan(3,3,nr_points);
    RR = nan(3,3,nr_points);
    
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
    
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta_L,'r',t,theta_R,'b',t,ones(nr_points,1)*mean(theta_L(1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(theta_R(1:(nr_points-1))),'b')
%     title('wing kinematics')
%     subplot(3,1,2); plot(t,eta_L,'r',t,eta_R,'b',t,ones(nr_points,1)*mean(eta_L(1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(eta_R(1:(nr_points-1))),'b')
%     subplot(3,1,3); plot(t,phi_L,'r',t,phi_R,'b',t,ones(nr_points,1)*mean(phi_L(1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(phi_R(1:(nr_points-1))),'b')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta_dot_L,'r',t,theta_dot_R,'b',t,ones(nr_points,1)*mean(theta_dot_L(1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(theta_dot_R(1:(nr_points-1))),'b')
%     title('euler angle derivatives')
%     subplot(3,1,2); plot(t,eta_dot_L,'r',t,eta_dot_R,'b',t,ones(nr_points,1)*mean(eta_dot_L(1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(eta_dot_R(1:(nr_points-1))),'b')
%     subplot(3,1,3); plot(t,phi_dot_L,'r',t,phi_dot_R,'b',t,ones(nr_points,1)*mean(phi_dot_L(1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(phi_dot_R(1:(nr_points-1))),'b')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta_ddot_L,'r',t,theta_ddot_R,'b',t,ones(nr_points,1)*mean(theta_ddot_L(1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(theta_ddot_R(1:(nr_points-1))),'b')
%     title('euler angle second derivative')
%     subplot(3,1,2); plot(t,eta_ddot_L,'r',t,eta_ddot_R,'b',t,ones(nr_points,1)*mean(eta_ddot_L(1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(eta_ddot_R(1:(nr_points-1))),'b')
%     subplot(3,1,3); plot(t,phi_ddot_L,'r',t,phi_ddot_R,'b',t,ones(nr_points,1)*mean(phi_ddot_L(1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(phi_ddot_R(1:(nr_points-1))),'b')
%     hold off
% 
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,wL(1,:),'r',t,wR(1,:),'b',t,ones(nr_points,1)*mean(wL(1,1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(wR(1,1:(nr_points-1))),'b')
%     title('wL & wR')
%     subplot(3,1,2); plot(t,wL(2,:),'r',t,wR(2,:),'b',t,ones(nr_points,1)*mean(wL(2,1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(wR(2,1:(nr_points-1))),'b')
%     subplot(3,1,3); plot(t,wL(3,:),'r',t,wR(3,:),'b',t,ones(nr_points,1)*mean(wL(3,1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(wR(3,1:(nr_points-1))),'b')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,w_dot_L(1,:),'r',t,w_dot_R(1,:),'b',t,ones(nr_points,1)*mean(w_dot_L(1,1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(w_dot_R(1,1:(nr_points-1))),'b')
%     title('w dot L & w dot R')
%     subplot(3,1,2); plot(t,w_dot_L(2,:),'r',t,w_dot_R(2,:),'b',t,ones(nr_points,1)*mean(w_dot_L(2,1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(w_dot_R(2,1:(nr_points-1))),'b')
%     subplot(3,1,3); plot(t,w_dot_L(3,:),'r',t,w_dot_R(3,:),'b',t,ones(nr_points,1)*mean(w_dot_L(3,1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(w_dot_R(3,1:(nr_points-1))),'b')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,wL_b(1,:),'r',t,wR_b(1,:),'b',t,ones(nr_points,1)*mean(wL_b(1,1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(wR_b(1,1:(nr_points-1))),'b')
%     title('wL b & wR b')
%     subplot(3,1,2); plot(t,wL_b(2,:),'r',t,wR_b(2,:),'b',t,ones(nr_points,1)*mean(wL_b(2,1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(wR_b(2,1:(nr_points-1))),'b')
%     subplot(3,1,3); plot(t,wL_b(3,:),'r',t,wR_b(3,:),'b',t,ones(nr_points,1)*mean(wL_b(3,1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(wR_b(3,1:(nr_points-1))),'b')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,w_dot_L_b(1,:),'r',t,w_dot_R_b(1,:),'b',t,ones(nr_points,1)*mean(w_dot_L_b(1,1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(w_dot_R_b(1,1:(nr_points-1))),'b')
%     title('w dot L b & w dot R b')
%     subplot(3,1,2); plot(t,w_dot_L_b(2,:),'r',t,w_dot_R_b(2,:),'b',t,ones(nr_points,1)*mean(w_dot_L_b(2,1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(w_dot_R_b(2,1:(nr_points-1))),'b')
%     subplot(3,1,3); plot(t,w_dot_L_b(3,:),'r',t,w_dot_R_b(3,:),'b',t,ones(nr_points,1)*mean(w_dot_L_b(3,1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(w_dot_R_b(3,1:(nr_points-1))),'b')
%     hold off    
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,wL_strk(1,:),'r',t,wR_strk(1,:),'b',t,ones(nr_points,1)*mean(wL_strk(1,1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(wR_strk(1,1:(nr_points-1))),'b')
%     title('wL strk & wR strk')
%     subplot(3,1,2); plot(t,wL_strk(2,:),'r',t,wR_strk(2,:),'b',t,ones(nr_points,1)*mean(wL_strk(2,1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(wR_strk(2,1:(nr_points-1))),'b')
%     subplot(3,1,3); plot(t,wL_strk(3,:),'r',t,wR_strk(3,:),'b',t,ones(nr_points,1)*mean(wL_strk(3,1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(wR_strk(3,1:(nr_points-1))),'b')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,w_dot_L_strk(1,:),'r',t,w_dot_R_strk(1,:),'b',t,ones(nr_points,1)*mean(w_dot_L_strk(1,1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(w_dot_R_strk(1,1:(nr_points-1))),'b')
%     title('w dot L strk & w dot R strk')
%     subplot(3,1,2); plot(t,w_dot_L_strk(2,:),'r',t,w_dot_R_strk(2,:),'b',t,ones(nr_points,1)*mean(w_dot_L_strk(2,1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(w_dot_R_strk(2,1:(nr_points-1))),'b')
%     subplot(3,1,3); plot(t,w_dot_L_strk(3,:),'r',t,w_dot_R_strk(3,:),'b',t,ones(nr_points,1)*mean(w_dot_L_strk(3,1:(nr_points-1))),'r',t,ones(nr_points,1)*mean(w_dot_R_strk(3,1:(nr_points-1))),'b')
%     hold off    
%     
%     Wt_L1 = [ 0.05; -1; 0];
%     Wt_L2 = [-0.05; -1 ; 0];
%     Wt_R1 = [ 0.05; 1; 0];
%     Wt_R2 = [-0.05; 1 ; 0];
%     
%     xyz_wing_L1 = zeros(3,nr_points);
%     xyz_wing_L2 = zeros(3,nr_points);
%     xyz_wing_R1 = zeros(3,nr_points);
%     xyz_wing_R2 = zeros(3,nr_points);
%     
%     for k = 1:nr_points
%         
%         xyz_wing_L1(:,k) = RL(:,:,k)'*Wt_L1;
%         xyz_wing_L2(:,k) = RL(:,:,k)'*Wt_L2;
%         xyz_wing_R1(:,k) = RR(:,:,k)'*Wt_R1;
%         xyz_wing_R2(:,k) = RR(:,:,k)'*Wt_R2;
%         
%     end
%     
%     figure()
%     hold on
%     for k = 1:nr_points
%         plot3([xyz_wing_L1(1,k) xyz_wing_L2(1,k)], [xyz_wing_L1(2,k) xyz_wing_L2(2,k)], [xyz_wing_L1(3,k) xyz_wing_L2(3,k)],'k')
%         plot3([xyz_wing_R1(1,k) xyz_wing_R2(1,k)], [xyz_wing_R1(2,k) xyz_wing_R2(2,k)], [xyz_wing_R1(3,k) xyz_wing_R2(3,k)],'k')
%     end
%     hold off
%     axis equal
%     
%     pause

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


end

