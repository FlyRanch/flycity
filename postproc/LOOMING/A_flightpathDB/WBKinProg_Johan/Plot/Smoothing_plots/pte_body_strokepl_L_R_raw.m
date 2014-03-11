function pte_body_strokepl_L_R_raw(settings,seq_nr,qL1_f2,qL2_f2,qL3_f2,qL4_f2,qR1_f2,qR2_f2,qR3_f2,qR4_f2,t,fig_nr,save_on_off)

    % Compares phi, theta and eta of the left wing with the phi, theta and
    % eta of the right wing for the filt2 quaternions
    
    N = length(qL1_f2);
    
    Lwt = [0; -1; 0];  
    
    Lwingtip_f2 = zeros(N,3);
        
    % Calculate wing_tip path
    
    for j = 1:N
        
        DCM_L = quat2matNEW([qL1_f2(j) qL2_f2(j) qL3_f2(j) qL4_f2(j)]);

        Lwingtip_f2(j,:) = DCM_L*Lwt;
        
    end
    
    phi_f2_L = -real(atan2(Lwingtip_f2(:,3),-Lwingtip_f2(:,2)));
    theta_f2_L = real(atan2(Lwingtip_f2(:,1),sqrt(Lwingtip_f2(:,2).^2+Lwingtip_f2(:,3).^2)));
    eta_f2_L = -real(atan2(2*(qL4_f2.*qL1_f2+qL2_f2.*qL3_f2),1-2*(qL1_f2.^2+qL2_f2.^2)));
    
    
    
    
    Rwt = [0; -1; 0];
    
    Rwingtip_f2 = zeros(N,3);
        
    % Calculate wing_tip path
    
    for j = 1:N
        
        DCM_R = quat2matNEW([qR1_f2(j) qR2_f2(j) qR3_f2(j) qR4_f2(j)]);

        Rwingtip_f2(j,:) = DCM_R*Rwt;
        
    end
    
    phi_f2_R = real(atan2(Rwingtip_f2(:,3),-Rwingtip_f2(:,2)));
    theta_f2_R = -real(atan2(Rwingtip_f2(:,1),sqrt(Rwingtip_f2(:,2).^2+Rwingtip_f2(:,3).^2)));
    eta_f2_R = real(atan2(2*(qR4_f2.*qR1_f2+qR2_f2.*qR3_f2),1-2*(qR1_f2.^2+qR2_f2.^2)));
    
    figure(fig_nr)
    subplot(3,1,1); plot(t,radtodeg(phi_f2_L),t,radtodeg(phi_f2_R))
    title('Comparison wing kinematics left and right, body mounted strokeplane, unfiltered')
    xlabel('t [s]')
    ylabel('phi [deg]')
    legend('left','right')
    subplot(3,1,2); plot(t,radtodeg(theta_f2_L),t,radtodeg(theta_f2_R))
    xlabel('t [s]')
    ylabel('theta [deg]')
    legend('left','right')
    subplot(3,1,3); plot(t,radtodeg(eta_f2_L),t,radtodeg(eta_f2_R))
    xlabel('t [s]')
    ylabel('eta [deg]')
    legend('left','right')
    
    % Save plots
    
    if save_on_off == 1
    
    saveas(fig_nr, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/pte_bref_lr_raw'], 'fig')
    
    saveas(fig_nr, [char(settings.plot_folders(2)) '/pte_bref_lr/pte_bref_lr_raw_' int2str(seq_nr)], 'fig')
    
    end
    

end