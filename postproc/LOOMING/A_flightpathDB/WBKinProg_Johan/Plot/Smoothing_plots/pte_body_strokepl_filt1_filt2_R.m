function pte_body_strokepl_filt1_filt2_R(settings, seq_nr, qR1_f1,qR2_f1,qR3_f1,qR4_f1,qR1_f2,qR2_f2,qR3_f2,qR4_f2,t,fig_nr, save_on_off)

    % Program that plots phi, theta and eta. Phi being the amplitude angle
    % within the stroke plane, theta the angle of deviation of the stroke
    % plane and eta the roll angle of the wing w.r.t to the stroke plane.
    % The stroke plane is defined as the y-z plane of the body reference
    % frame of the fly.
    
    N = length(qR1_f1);
    
    Rwingtip_f1 = zeros(N,3);
    
    Rwt = [0; -1; 0];
        
    % Calculate wing_tip path
    
    for j = 1:N
        
        DCM_R = quat2matNEW([qR1_f1(j) qR2_f1(j) qR3_f1(j) qR4_f1(j)]);

        Rwingtip_f1(j,:) = DCM_R*Rwt;
        
    end
    
    phi_f1 = real(atan2(Rwingtip_f1(:,3),-Rwingtip_f1(:,2)));
    theta_f1 = -real(atan2(Rwingtip_f1(:,1),sqrt(Rwingtip_f1(:,2).^2+Rwingtip_f1(:,3).^2)));
    eta_f1 = real(atan2(2*(qR4_f1.*qR1_f1+qR2_f1.*qR3_f1),1-2*(qR1_f1.^2+qR2_f1.^2)));
    
    
    Rwingtip_f2 = zeros(N,3);
        
    % Calculate wing_tip path
    
    for j = 1:N
        
        DCM_R = quat2matNEW([qR1_f2(j) qR2_f2(j) qR3_f2(j) qR4_f2(j)]);

        Rwingtip_f2(j,:) = DCM_R*Rwt;
        
    end
    
    phi_f2 = real(atan2(Rwingtip_f2(:,3),-Rwingtip_f2(:,2)));
    theta_f2 = -real(atan2(Rwingtip_f2(:,1),sqrt(Rwingtip_f2(:,2).^2+Rwingtip_f2(:,3).^2)));
    eta_f2 = real(atan2(2*(qR4_f2.*qR1_f2+qR2_f2.*qR3_f2),1-2*(qR1_f2.^2+qR2_f2.^2)));
    
    
    figure(fig_nr)
    subplot(3,1,1); plot(t,radtodeg(phi_f1),t,radtodeg(phi_f2))
    title('Unfiltered and filtered wing kinematics with body mounted stroke plane for the right wing')
    xlabel('t [s]')
    ylabel('phi [deg]')
    subplot(3,1,2); plot(t,radtodeg(theta_f1),t,radtodeg(theta_f2))
    xlabel('t [s]')
    ylabel('theta [deg]')
    subplot(3,1,3); plot(t,radtodeg(eta_f1),t,radtodeg(eta_f2))
    xlabel('t [s]')
    ylabel('eta [deg]')
    
    % Save plots
    
    if save_on_off == 1
    
    saveas(fig_nr, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/pte_right'], 'fig')
    
    saveas(fig_nr, [char(settings.plot_folders(2)) '/pte_right/pte_right_' int2str(seq_nr)], 'fig')
    
    end

end



