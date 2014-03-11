function [ phi_L, theta_L, eta_L, phi_R, theta_R, eta_R ] = wing_stroke_plane( qL1, qL2, qL3, qL4, qR1, qR2, qR3, qR4 )

    % Calculate the wingstroke amplitude angle, deviation angle and wing
    % pitch angle for the left and right wing.
    
        %Left wingtip parameters

    Lwingtip = zeros(length(qL1),3);
    
    Lwt = [0; -1; 0];
    
    Lwingtip2 = zeros(length(qL1),3);
    
    Lwt2 = [0.01; -sqrt(1-0.01^2); 0];
    
    N = length(qL1);
    
    % Calculate wing_tip path
    
    for j = 1:N
        
        DCM_L_r = quat2matNEW([qL1(j) qL2(j) qL3(j) qL4(j)]);

        Lwingtip(j,:) = DCM_L_r*Lwt;
        
        Lwingtip2(j,:) = DCM_L_r*Lwt2;
        
    end
    
    
    
    
    % Right wingtip parameters
    
    Rwingtip = zeros(length(qR1),3);
    
    Rwt = [0; 1; 0];
    
    Rwingtip2 = zeros(length(qR1),3);
    
    Rwt2 = [0.01; sqrt(1-0.01^2); 0];
    
    
    % Calculate wing_tip path
    
    for j = 1:N
        
        DCM_R_r = quat2matNEW([qR1(j) qR2(j) qR3(j) qR4(j)]);

        Rwingtip(j,:) = DCM_R_r*Rwt;
        
        Rwingtip2(j,:) = DCM_R_r*Rwt2;
        
    end
    
    
        % Stroke plane angle is chosen to be constant at 55 degrees.
    
    beta = -(55/180)*pi;
    
    Rot_mat = [cos(beta) 0 -sin(beta); ...
                   0 1 0; ...
                   sin(beta) 0 cos(beta)];
    
    WTL_stroke = zeros(N,3);
    
    WTR_stroke = zeros(N,3);
    
    WTL_stroke2 = zeros(N,3);
    
    WTR_stroke2 = zeros(N,3);
    
    for j = 1:N
        
        % Rotate the wingtip coordinates from the body reference frame to
        % the strokesplane coordinates.
               
        WTL_stroke(j,:) = Rot_mat*Lwingtip(j,:)';
        
        WTR_stroke(j,:) = Rot_mat*Rwingtip(j,:)';
        
        WTL_stroke2(j,:) = Rot_mat*Lwingtip2(j,:)';
        
        WTR_stroke2(j,:) = Rot_mat*Rwingtip2(j,:)';
        
    end
    
    phi_L = -real(atan2(WTL_stroke(:,1),-WTL_stroke(:,2)));
    theta_L = -real(asin(WTL_stroke(:,3)./sqrt(WTL_stroke(:,1).^2+WTL_stroke(:,2).^2+WTL_stroke(:,3).^2)));
    
    phi_L2 = -real(atan2(WTL_stroke2(:,1),-WTL_stroke2(:,2)));   
    theta_L2 = -real(asin(WTL_stroke2(:,3)./sqrt(WTL_stroke2(:,1).^2+WTL_stroke2(:,2).^2+WTL_stroke2(:,3).^2)));
    
    d_phi_L = phi_L2-phi_L;
    d_theta_L = theta_L2-theta_L;

    
    eta_L = real(atan2(d_theta_L(:),d_phi_L(:)));
    
    phi_R = -real(atan2(WTR_stroke(:,1),WTR_stroke(:,2)));
    theta_R = -real(asin(WTR_stroke(:,3)./sqrt(WTR_stroke(:,1).^2+WTR_stroke(:,2).^2+WTR_stroke(:,3).^2)));
    
    phi_R2 = -real(atan2(WTR_stroke2(:,1),WTR_stroke2(:,2)));
    theta_R2 = -real(asin(WTR_stroke2(:,3)./sqrt(WTR_stroke2(:,1).^2+WTR_stroke2(:,2).^2+WTR_stroke2(:,3).^2)));
    
    d_phi_R = phi_R2-phi_R;
    d_theta_R = theta_R2-theta_R;
    
    eta_R = real(atan2(d_theta_R(:),d_phi_R(:)));
    
end