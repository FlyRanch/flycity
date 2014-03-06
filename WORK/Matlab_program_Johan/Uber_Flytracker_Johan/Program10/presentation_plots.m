function presentation_plots(settings, pathDB)


    a_avg_A = pathDB.a_avg.a_avg_40;
    
    a_dev_A = pathDB.a_dev.a_dev_40;
    
    
    down_up_avg_A = pathDB.down_up_avg.down_up_avg_40;
    
    a_theta_L_A = [a_avg_A.theta_L1; a_avg_A.theta_L2];
    a_theta_R_A = [a_avg_A.theta_R1; a_avg_A.theta_R2];
    a_eta_L_A = [a_avg_A.eta_L1; a_avg_A.eta_L2];
    a_eta_R_A = [a_avg_A.eta_R1; a_avg_A.eta_R2];
    a_phi_L_A = [a_avg_A.phi_L1; a_avg_A.phi_L2];
    a_phi_R_A = [a_avg_A.phi_R1; a_avg_A.phi_R2];
    a_theta_LR_A = [a_avg_A.theta_LR1; a_avg_A.theta_LR2];
    a_eta_LR_A = [a_avg_A.eta_LR1; a_avg_A.eta_LR2];
    a_phi_LR_A = [a_avg_A.phi_LR1; a_avg_A.phi_LR2];
    
    a_dev_theta_L_A = [a_dev_A.theta_L1(:,end-3); a_dev_A.theta_L2(:,end-3)];
    a_dev_theta_R_A = [a_dev_A.theta_R1(:,end-3); a_dev_A.theta_R2(:,end-3)];
    a_dev_eta_L_A = [a_dev_A.eta_L1(:,end-3); a_dev_A.eta_L2(:,end-3)];
    a_dev_eta_R_A = [a_dev_A.eta_R1(:,end-3); a_dev_A.eta_R2(:,end-3)];
    a_dev_phi_L_A = [a_dev_A.phi_L1(:,end-3); a_dev_A.phi_L2(:,end-3)];
    a_dev_phi_R_A = [a_dev_A.phi_R1(:,end-3); a_dev_A.phi_R2(:,end-3)];
    
    [ t, X_theta ] = Wingbeat_Legendre_matrix( 12, down_up_avg_A, 100, 0, 1, 0 );
    [ ~, X_eta ] = Wingbeat_Legendre_matrix( 14, down_up_avg_A, 100, 0, 1, 0 );
    [ ~, X_phi ] = Wingbeat_Legendre_matrix( 10, down_up_avg_A, 100, 0, 1, 0 );
    
    theta_L_A_avg = X_theta*a_theta_L_A;
    theta_R_A_avg = X_theta*a_theta_R_A;
    theta_LR_A_avg = X_theta*a_theta_LR_A;
    
    eta_L_A_avg = X_eta*a_eta_L_A;
    eta_R_A_avg = X_eta*a_eta_R_A;
    eta_LR_A_avg = X_eta*a_eta_LR_A;
    
    phi_L_A_avg = X_phi*a_phi_L_A;
    phi_R_A_avg = X_phi*a_phi_R_A;
    phi_LR_A_avg = X_phi*a_phi_LR_A;
    
    theta_L_A_dev = X_theta*a_dev_theta_L_A;
    theta_R_A_dev = X_theta*a_dev_theta_R_A;
    
    eta_L_A_dev = X_eta*a_dev_eta_L_A;
    eta_R_A_dev = X_eta*a_dev_eta_R_A;
    
    phi_L_A_dev = X_phi*a_dev_phi_L_A;
    phi_R_A_dev = X_phi*a_dev_phi_R_A;
    
    
    figure()
    hold on
    subplot(3,1,1); plot(t,radtodeg(theta_L_A_avg),'r',t,radtodeg(theta_R_A_avg),'b',t,radtodeg(theta_LR_A_avg),'m')
    title('\theta')
    xlabel('wingbeat')
    ylabel('angle [deg]')
    legend('avg left','avg right','avg symmetric')
    subplot(3,1,2); plot(t,radtodeg(eta_L_A_avg),'r',t,radtodeg(eta_R_A_avg),'b',t,radtodeg(eta_LR_A_avg),'m')
    title('\eta')
    xlabel('wingbeat')
    ylabel('angle [deg]')
    legend('avg left','avg right','avg symmetric')
    subplot(3,1,3); plot(t,radtodeg(phi_L_A_avg),'r',t,radtodeg(phi_R_A_avg),'b',t,radtodeg(phi_LR_A_avg),'m')
    title('\phi')
    xlabel('wingbeat')
    ylabel('angle [deg]')
    legend('avg left','avg right','avg symmetric')
    hold off
    
    figure()
    hold on
    subplot(2,1,1); plot(t,radtodeg(theta_LR_A_avg+theta_L_A_dev),'k',t,radtodeg(theta_LR_A_avg),'r',t,radtodeg(theta_L_A_dev),'b')
    title('\theta left wing')
    xlabel('wingbeat')
    ylabel('angle [deg]')
    legend('actual','avg','dev')
    subplot(2,1,2); plot(t,radtodeg(theta_LR_A_avg+theta_R_A_dev),'k',t,radtodeg(theta_LR_A_avg),'r',t,radtodeg(theta_R_A_dev),'b')
    title('\theta right wing')
    xlabel('wingbeat')
    ylabel('angle [deg]')
    legend('actual','avg','dev')
    hold off

end
