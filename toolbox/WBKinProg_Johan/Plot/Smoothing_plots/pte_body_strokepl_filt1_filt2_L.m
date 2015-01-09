function pte_body_strokepl_filt1_filt2_L(settings, seq_nr, qL1_f1,qL2_f1,qL3_f1,qL4_f1,qL1_f2,qL2_f2,qL3_f2,qL4_f2,t,fig_nr, save_on_off)

    % Program that plots phi, theta and eta. Phi being the amplitude angle
    % within the stroke plane, theta the angle of deviation of the stroke
    % plane and eta the roll angle of the wing w.r.t to the stroke plane.
    % The stroke plane is defined as the y-z plane of the body reference
    % frame of the fly.
    
    N = length(qL1_f1);
    
    Lwingtip_f1 = zeros(N,3);
    
    Lwt = [0; -1; 0];
        
    % Calculate wing_tip path
    
    for j = 1:N
        
        DCM_L = quat2matNEW([qL1_f1(j) qL2_f1(j) qL3_f1(j) qL4_f1(j)]);

        Lwingtip_f1(j,:) = DCM_L*Lwt;
        
    end
    
    phi_f1 = -real(atan2(Lwingtip_f1(:,3),-Lwingtip_f1(:,2)));
    theta_f1 = real(atan2(Lwingtip_f1(:,1),sqrt(Lwingtip_f1(:,2).^2+Lwingtip_f1(:,3).^2)));
    eta_f1 = -real(atan2(2*(qL4_f1.*qL1_f1+qL2_f1.*qL3_f1),1-2*(qL1_f1.^2+qL2_f1.^2)));
    
    clear DCM_L Lwingtip_f1
    
    Lwingtip_f2 = zeros(N,3);
        
    % Calculate wing_tip path
    
    for j = 1:N
        
        DCM_L = quat2matNEW([qL1_f2(j) qL2_f2(j) qL3_f2(j) qL4_f2(j)]);

        Lwingtip_f2(j,:) = DCM_L*Lwt;
        
    end
    
    phi_f2 = -real(atan2(Lwingtip_f2(:,3),-Lwingtip_f2(:,2)));
    theta_f2 = real(atan2(Lwingtip_f2(:,1),sqrt(Lwingtip_f2(:,2).^2+Lwingtip_f2(:,3).^2)));
    eta_f2 = -real(atan2(2*(qL4_f2.*qL1_f2+qL2_f2.*qL3_f2),1-2*(qL1_f2.^2+qL2_f2.^2)));
    
    clear DCM_L Lwingtip_f2
    
    figure(fig_nr)
    subplot(3,1,1); plot(t,radtodeg(phi_f1),t,radtodeg(phi_f2))
    title('Unfiltered and filtered wing kinematics with body mounted stroke plane for the left wing')
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
    
    saveas(fig_nr, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/pte_left'], 'fig')
    
    saveas(fig_nr, [char(settings.plot_folders(2)) '/pte_left/pte_left' int2str(seq_nr)], 'fig')
    
    end
    

end

