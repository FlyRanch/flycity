function Wingkinematics_Legendre(settings, pathDB)


    % Use a Legendre polynomial fit to describe the wingkinematics in a
    % minimal number of parameters:
    
    
    
    for i =1:size(pathDB.x,2)
    
    seq_nr = i
    
    nr_wb = find(isnan(pathDB.wingbeat_time(:,1,seq_nr))==0, 1, 'last' );


    
    
%     n_pol_theta = 13; % Order of used polynomials
%     n_pol_eta = 16; % Order of used polynomials
%     n_pol_phi = 8; % Order of used polynomials
%     
    n_pol_theta = 10; % Order of used polynomials
    n_pol_eta = 10; % Order of used polynomials
    n_pol_phi = 10; % Order of used polynomials

    n_pol = {};
    
    n_pol.theta = n_pol_theta;
    
    n_pol.eta = n_pol_eta;
    
    n_pol.phi = n_pol_phi;

    [a_fit_L,a_fit_R,a_avg_L,a_avg_R,a_avg_LR,f_avg,down_up,trigger_wb] = Standard_wingbeat( settings, pathDB, i , n_pol_theta, n_pol_eta, n_pol_phi);
    
    down_up
    
    f_avg
    
    trigger_wb
    
    t = -1:0.01:1;
    
    phi_3D_L = zeros(nr_wb,201);
    theta_3D_L = zeros(nr_wb,201);
    eta_3D_L = zeros(nr_wb,201);
    
    for j = 1:nr_wb
        
        PN_theta = Legendre_polynomial(n_pol_theta,1,t);
        PN_eta = Legendre_polynomial(n_pol_eta,1,t);
        PN_phi = Legendre_polynomial(n_pol_phi,1,t);
        eta_3D_L(j,:) = PN_eta(:,:,1)'*a_fit_L.eta(:,j);
        theta_3D_L(j,:) = PN_theta(:,:,1)'*a_fit_L.theta(:,j);
        phi_3D_L(j,:) = PN_phi(:,:,1)'*a_fit_L.phi(:,j);
        
    end
    
    % Plot left wingbeats over time:

    figure()
    surf(theta_3D_L)
    title('theta left')
    
    figure()
    surf(eta_3D_L)
    title('eta left')
    
    figure()
    surf(phi_3D_L)
    title('phi left')
    
    phi_3D_R = zeros(nr_wb,201);
    theta_3D_R = zeros(nr_wb,201);
    eta_3D_R = zeros(nr_wb,201);
    
    for j = 1:nr_wb
        
        PN_theta = Legendre_polynomial(n_pol_theta,1,t);
        PN_eta = Legendre_polynomial(n_pol_eta,1,t);
        PN_phi = Legendre_polynomial(n_pol_phi,1,t);
        eta_3D_R(j,:) = PN_eta(:,:,1)'*a_fit_R.eta(:,j);
        theta_3D_R(j,:) = PN_theta(:,:,1)'*a_fit_R.theta(:,j);
        phi_3D_R(j,:) = PN_phi(:,:,1)'*a_fit_R.phi(:,j);
        
    end
    
    % Plot right wingbeats over time:
    
    figure()
    surf(theta_3D_R)
    title('theta right')
    
    figure()
    surf(eta_3D_R)
    title('eta right')
    
    figure()
    surf(phi_3D_R)    
    title('phi right')    
    
    
    % Plot average wingbeats left
    
    
    % Plot average wingbeats right
    
    if trigger_wb > 3
    
    figure()
    plot(radtodeg(PN_theta(:,:,1)'*a_avg_L.theta),'r')
    hold on
    plot(radtodeg(PN_theta(:,:,1)'*a_avg_R.theta),'g')
    plot(radtodeg(PN_theta(:,:,1)'*a_avg_LR.theta),'b')
    hold off

    figure()
    plot(radtodeg(PN_eta(:,:,1)'*a_avg_L.eta),'r')
    hold on
    plot(radtodeg(PN_eta(:,:,1)'*a_avg_R.eta),'g')
    plot(radtodeg(PN_eta(:,:,1)'*a_avg_LR.eta),'b')
    hold off
    
    figure()
    plot(radtodeg(PN_phi(:,:,1)'*a_avg_L.phi),'r')
    hold on
    plot(radtodeg(PN_phi(:,:,1)'*a_avg_R.phi),'g')
    plot(radtodeg(PN_phi(:,:,1)'*a_avg_LR.phi),'b')
    hold off
    
    
    figure()
    hold on
    for k = 3:trigger_wb
        plot(radtodeg(PN_theta(:,:,1)'*a_fit_L.theta(:,k)),'Color',[0.5 0.5 0.5])
    end
    plot(radtodeg(PN_theta(:,:,1)'*a_avg_L.theta),'r')
    hold off
    
    figure()
    hold on
    for k = 3:trigger_wb
        plot(radtodeg(PN_theta(:,:,1)'*a_fit_R.theta(:,k)),'Color',[0.5 0.5 0.5])
    end
    plot(radtodeg(PN_theta(:,:,1)'*a_avg_R.theta),'r')
    hold off
    
    figure()
    hold on
    for k = 3:trigger_wb
        plot(radtodeg(PN_eta(:,:,1)'*a_fit_L.eta(:,k)),'Color',[0.5 0.5 0.5])
    end
    plot(radtodeg(PN_eta(:,:,1)'*a_avg_L.eta),'r')
    hold off
    
    figure()
    hold on
    for k = 3:trigger_wb
        plot(radtodeg(PN_eta(:,:,1)'*a_fit_R.eta(:,k)),'Color',[0.5 0.5 0.5])
    end
    plot(radtodeg(PN_eta(:,:,1)'*a_avg_R.eta),'r')
    hold off
    
    figure()
    hold on
    for k = 3:trigger_wb
        plot(radtodeg(PN_phi(:,:,1)'*a_fit_L.phi(:,k)),'Color',[0.5 0.5 0.5])
    end
    plot(radtodeg(PN_phi(:,:,1)'*a_avg_L.phi),'r')
    hold off
    
    figure()
    hold on
    for k = 3:trigger_wb
        plot(radtodeg(PN_phi(:,:,1)'*a_fit_R.phi(:,k)),'Color',[0.5 0.5 0.5])
    end
    plot(radtodeg(PN_phi(:,:,1)'*a_avg_R.phi),'r')
    hold off
    
    end
    
    figure()
    plot(pathDB.omegax_body_mean(1:nr_wb,1,seq_nr))
    
    figure()
    plot(pathDB.omegay_body_mean(1:nr_wb,1,seq_nr))
    
    figure()
    plot(pathDB.omegaz_body_mean(1:nr_wb,1,seq_nr))


    start = find(isnan(pathDB.x(:,seq_nr))==0, 1 );
    stop = find(isnan(pathDB.x(:,seq_nr))==0, 1, 'last' );

    omega_x_body = pathDB.b_omega1(start:stop,seq_nr);
    omega_y_body = pathDB.b_omega2(start:stop,seq_nr);
    omega_z_body = pathDB.b_omega3(start:stop,seq_nr);
    
    % Create omega_x, omega_y and omega_z for the strokeplane reference
    % frame:
    
    omega_x_strkpln = cosd(-55).*omega_x_body-sind(-55).*omega_z_body;
    omega_y_strkpln = omega_y_body;
    omega_z_strkpln = sind(-55).*omega_x_body+cosd(-55).*omega_z_body;
    
    
    figure()
    subplot(3,1,1); plot(omega_x_strkpln)
    subplot(3,1,2); plot(omega_y_strkpln)
    subplot(3,1,3); plot(omega_z_strkpln)
    
    Maneuver_wingbeat(a_fit_L,a_fit_R,a_avg_L,a_avg_R,a_avg_LR,n_pol,trigger_wb)
    
    
    pause
    
    close all
    
    end
    

    
    
    
    
    

end

