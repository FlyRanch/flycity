function [ theta_L, eta_L, phi_L, theta_R, eta_R, phi_R ] = theta_eta_phi(RL,RR,Rstr)

    % Convert the left and right wing rotation matrices into wing
    % kinematics:

    LWT = [0; -1; 0];
    LWT2 = [0.01; -sqrt(1-0.01^2); 0];
    RWT = [0; 1; 0];
    RWT2 = [0.01; sqrt(1-0.01^2); 0];
    
    wingtip_left    = Rstr*RL'*LWT;
    wingtip_left2   = Rstr*RL'*LWT2;
    
    wingtip_right   = Rstr*RR'*RWT;
    wingtip_right2  = Rstr*RR'*RWT2;
    
    phi_L       = real(atan2(-wingtip_left(1),-wingtip_left(2)));
    theta_L     = real(asin(-wingtip_left(3)/norm(wingtip_left)));
    
    phi_L2      = real(atan2(-wingtip_left2(1),-wingtip_left2(2)));   
    theta_L2    = real(asin(-wingtip_left2(3)/norm(wingtip_left2)));
    
    d_phi_L     = phi_L2-phi_L;
    d_theta_L   = theta_L2-theta_L;
    
    eta_L = real(atan2(d_theta_L,d_phi_L));
    
    if eta_L < 0 && abs(eta_L) >= (pi/2)
        
        eta_L = abs(eta_L);
        
    end
    
    phi_R       = real(atan2(-wingtip_right(1),wingtip_right(2)));
    theta_R     = real(asin(-wingtip_right(3)/norm(wingtip_right)));
    
    phi_R2      = real(atan2(-wingtip_right2(1),wingtip_right2(2)));   
    theta_R2    = real(asin(-wingtip_right2(3)/norm(wingtip_right2)));
    
    d_phi_R     = phi_R2-phi_R;
    d_theta_R   = theta_R2-theta_R;

    
    eta_R = real(atan2(d_theta_R,d_phi_R));
    
    if eta_R < 0 && abs(eta_R) >= (pi/2)
        
        eta_R = abs(eta_R);
        
%     elseif eta_R < 0 && abs(eta_R) < (pi/2)
%         
%         eta_R = 
        
    end
    
end

