function Stroke_plane_pte_prototype(qL1, qL2, qL3, qL4, qR1, qR2, qR3, qR4, t, fig_nr1, fig_nr2)
    
    % Program which calculates the position of the stroke plane and plots
    % the angles phi, theta and eta w.r.t. this stroke plane.
    
%     r_L = real(atan2(2*(qL4.*qL1+qL2.*qL3),1-2*(qL1.^2+qL2.^2)));
%     p_L = real(asin(2*(qL4.*qL2-qL3.*qL1)));
%     y_L = real(atan2(2*(qL4.*qL3+qL1.*qL2),1-2*(qL2.^2+qL3.^2)));

%     [Gamma_downy loc_downy] = findpeaks(-y_L, 'minpeakdistance' , 25)
%     [Gamma_upy loc_upy]  = findpeaks(, 'minpeakdistance' , 25)
    
    Lwingtip = zeros(length(qL1),3);
    
    Lwt = [0; -1; 0];
    
    N = length(qL1);
        
    % Calculate wing_tip path
    
    for j = 1:N
        
        DCM_L_r = quat2matNEW([qL1(j) qL2(j) qL3(j) qL4(j)]);

        Lwingtip(j,:) = DCM_L_r*Lwt;
        
    end
    
    phi_L = -real(atan2(Lwingtip(:,3),-Lwingtip(:,2)));
    
    theta_L = real(atan2(Lwingtip(:,1),sqrt(Lwingtip(:,2).^2+Lwingtip(:,3).^2)));
    
    eta_L = -real(atan2(2*(qL4.*qL1+qL2.*qL3),1-2*(qL1.^2+qL2.^2)));
    
%     figure()
%     plot(radtodeg(phi_L))
%     
%     figure()
%     plot(radtodeg(theta_L))
%     
%     figure()
%     plot(radtodeg(eta_L))
%     
%     figure()
%     plot(Lwingtip(:,2))
    
    
    
    
    
    Rwingtip = zeros(length(qR1),3);
    
    Rwt = [0; -1; 0];
    
    
    % Calculate wing_tip path
    
    for j = 1:N
        
        DCM_R_r = quat2matNEW([qR1(j) qR2(j) qR3(j) qR4(j)]);

        Rwingtip(j,:) = DCM_R_r*Rwt;
        
    end
    
    phi_R = real(atan2(Rwingtip(:,3),-Rwingtip(:,2)));
    
    theta_R = -real(atan2(Rwingtip(:,1),sqrt(Rwingtip(:,2).^2+Rwingtip(:,3).^2)));
    
    eta_R = real(atan2(2*(qR4.*qR1+qR2.*qR3),1-2*(qR1.^2+qR2.^2)));
    
%     figure()
%     plot(radtodeg(phi_R))
%     
%     figure()
%     plot(radtodeg(theta_R))
%     
%     figure()
%     plot(radtodeg(eta_R))
%     
%     figure()
%     plot(Rwingtip(:,2))
    
    
    [GammaL locL] = findpeaks(Lwingtip(:,2), 'minpeakdistance' , 8);
   
    
    BetaL = zeros(length(locL)-1,1);
    
    for k = 1:(length(locL)-1)
        
        if GammaL(k) > GammaL(k+1)
            
            BetaL(k) = real(atan2(abs(Lwingtip(locL(k),3)-Lwingtip(locL(k+1),3)),abs(Lwingtip(locL(k),1)-Lwingtip(locL(k+1),1))));
            
        elseif GammaL(k+1) > GammaL(k)
            
            BetaL(k) = real(atan2(abs(Lwingtip(locL(k+1),3)-Lwingtip(locL(k),3)),abs(Lwingtip(locL(k+1),1)-Lwingtip(locL(k),1))));
            
        end
        
    end
    
    
    [GammaR locR] = findpeaks(Rwingtip(:,2), 'minpeakdistance' , 8);
   
    
    BetaR = zeros(length(locR)-1,1);
    
    for k = 1:(length(locR)-1)
        
        if GammaR(k) > GammaR(k+1)
            
            BetaR(k) = real(atan2(abs(Rwingtip(locR(k),3)-Rwingtip(locR(k+1),3)),abs(Rwingtip(locR(k),1)-Rwingtip(locR(k+1),1))));
            
        elseif GammaR(k+1) > GammaR(k)
            
            BetaR(k) = real(atan2(abs(Rwingtip(locR(k+1),3)-Rwingtip(locR(k),3)),abs(Rwingtip(locR(k+1),1)-Rwingtip(locR(k),1))));
            
        end
        
    end
    
    
    % Now calculate the average wing-stroke plane per half wingbeat:
    
    Beta_avg_L = zeros(length(locL)-1,1);
    
    
    for k = 2:length(locL)
        
        temp_x = zeros(locL(k)-locL(k-1)-1,1);
        temp_z = zeros(locL(k)-locL(k-1)-1,1);
        
        j = 1;
        
        for i = (locL(k-1)):(locL(k))
            
            temp_x(j) = Lwingtip(i,1);
            
            temp_z(j) = Lwingtip(i,3);
            
            j = j+1;
        end
        
        [r,m,b] = regression(temp_z,temp_x,'one');
        
        Beta_avg_L(k-1) =  (pi/2)-real(atan2(m,1));
        
        clear r m b temp_x temp_z
    end
    
    
    
    % Now calculate the average wing-stroke plane per half wingbeat:
    
    Beta_avg_R = zeros(length(locR)-1,1);
    
    
    for k = 2:length(locR)
        
        temp_x = zeros(locR(k)-locR(k-1)-1,1);
        temp_z = zeros(locR(k)-locR(k-1)-1,1);
        
        j = 1;
        
        for i = (locR(k-1)):(locR(k))
            
            temp_x(j) = Rwingtip(i,1);
            
            temp_z(j) = Rwingtip(i,3);
            
            j = j+1;
        end
        
        [r,m,b] = regression(temp_z,temp_x,'one');
        
        Beta_avg_R(k-1) =  (pi/2)-real(atan2(m,1));
        
        clear r m b temp_x temp_z
    end
    
    
    
    figure()
    plot(t(locL(1:(length(locL)-1))),radtodeg(BetaL),t(locR(1:(length(locR)-1))),radtodeg(BetaR))
    title('Left and right wing-stroke plane')
    legend('left', 'right')
    
    figure()
    plot(t(locL(1:(length(locL)-1))),radtodeg(Beta_avg_L),t(locL(1:(length(locR)-1))),radtodeg(Beta_avg_R))
  


end

