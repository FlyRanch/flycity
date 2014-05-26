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
    n_pol_theta = 13; % Order of used polynomials
    n_pol_eta = 14; % Order of used polynomials
    n_pol_phi = 11; % Order of used polynomials
%     n_pol_theta = 9; % Order of used polynomials
%     n_pol_eta = 9; % Order of used polynomials
%     n_pol_phi = 7; % Order of used polynomials

    n_pol = {};
    
    n_pol.theta = n_pol_theta;
    
    n_pol.eta = n_pol_eta;
    
    n_pol.phi = n_pol_phi;


%     [a_fit, a_avg,f_avg,down_up,trigger_wb] = Standard_wingbeat( settings, pathDB, i , n_pol_theta, n_pol_eta, n_pol_phi);

    [a_fit,a_avg,f_avg,down_up,trigger_wb,ratio_1_2,ratio_1_2_avg] = Standard_wingbeat( settings, pathDB, i, n_pol_theta, n_pol_eta, n_pol_phi);



    a_fit
    
    a_avg
    
    down_up
    
    f_avg
    
    trigger_wb
    
    pause
    
    t = -1:0.01:1;

    
    PN_theta = Legendre_polynomial(n_pol_theta,1,t);
    PN_eta = Legendre_polynomial(n_pol_eta,1,t);
    PN_phi = Legendre_polynomial(n_pol_phi,1,t);
    
    theta_fit_L = zeros(nr_wb,401);
    eta_fit_L = zeros(nr_wb,401);
    phi_fit_L = zeros(nr_wb,401);
    
    theta_fit_R = zeros(nr_wb,401);
    eta_fit_R = zeros(nr_wb,401);
    phi_fit_R = zeros(nr_wb,401);
    
    for j = 1:(2*nr_wb)
        
        if mod(j,2) == 1
            eta_fit_L(ceil(j/2),1:200) = PN_eta(:,1:200,1)'*a_fit.eta_L1(:,ceil(j/2));
            theta_fit_L(ceil(j/2),1:200) = PN_theta(:,1:200,1)'*a_fit.theta_L1(:,ceil(j/2));
            phi_fit_L(ceil(j/2),1:200) = PN_phi(:,1:200,1)'*a_fit.phi_L1(:,ceil(j/2));
            eta_fit_R(ceil(j/2),1:200) = PN_eta(:,1:200,1)'*a_fit.eta_R1(:,ceil(j/2));
            theta_fit_R(ceil(j/2),1:200) = PN_theta(:,1:200,1)'*a_fit.theta_R1(:,ceil(j/2));
            phi_fit_R(ceil(j/2),1:200) = PN_phi(:,1:200,1)'*a_fit.phi_R1(:,ceil(j/2));
        elseif mod(j,2) == 0
            eta_fit_L(j/2,201:401) = PN_eta(:,:,1)'*a_fit.eta_L2(:,j/2);
            theta_fit_L(j/2,201:401) = PN_theta(:,:,1)'*a_fit.theta_L2(:,j/2);
            phi_fit_L(j/2,201:401) = PN_phi(:,:,1)'*a_fit.phi_L2(:,j/2);
            eta_fit_R(j/2,201:401) = PN_eta(:,:,1)'*a_fit.eta_R2(:,j/2);
            theta_fit_R(j/2,201:401) = PN_theta(:,:,1)'*a_fit.theta_R2(:,j/2);
            phi_fit_R(j/2,201:401) = PN_phi(:,:,1)'*a_fit.phi_R2(:,j/2);
        end
        
    end
    
    delta_x = 0.001;
    
    x = 0:delta_x:1;
    
    x_dist = [0; 0.3; 0.4; 0.8; 1];
    
    loc_max_min = zeros(length(x_dist),1);
    
    for j = 1:length(x_dist)
        
        
        loc_max_min(j) = round(1+x_dist(j)/delta_x);
        
    end
    
    
%     Mn = Melisoids( loc_max_min, x);
%     
%     figure()
%     plot(Mn*[1; 0.5; 2; -3; 1])
%     
%     pause
    
    % Plot left wingbeats over time:

%     figure()
%     surf(eta_fit_L)
%     title('theta left')
%     
%     figure()
%     surf(theta_fit_L)
%     title('eta left')
%     
%     figure()
%     surf(phi_fit_L)
%     title('phi left')
%  
%     figure()
%     surf(eta_fit_R)
%     title('theta right')
%     
%     figure()
%     surf(theta_fit_R)
%     title('eta right')
%     
%     figure()
%     surf(phi_fit_R)
%     title('phi right')
    

 
%     ratio_1_2 = zeros(nr_wb,1);
%     
%     wb_loc_12 = a_fit.wb_loc_12;
%     
%     for j = 1:nr_wb
%         
%         ratio_1_2(j) = (wb_loc_12(j*2-1,2)-wb_loc_12(j*2-1,1)+1)/(wb_loc_12(j*2,2)-wb_loc_12(j*2,1)+1);
%         
%     end
%     
%     mean(ratio_1_2(3:trigger_wb))
    
    % Plot average wingbeats left
    
    
    % Plot average wingbeats right
    
    if trigger_wb > 3
    

        
    theta_avg_L = zeros(1,401);
    eta_avg_L = zeros(1,401);
    phi_avg_L = zeros(1,401);
    
    theta_avg_R = zeros(1,401);
    eta_avg_R = zeros(1,401);
    phi_avg_R = zeros(1,401);
    
    theta_avg_LR = zeros(1,401);
    eta_avg_LR = zeros(1,401);
    phi_avg_LR = zeros(1,401);
    

    eta_avg_L(1:200) = PN_eta(:,1:200,1)'*a_avg.eta_L1;
    theta_avg_L(1:200) = PN_theta(:,1:200,1)'*a_avg.theta_L1;
    phi_avg_L(1:200) = PN_phi(:,1:200,1)'*a_avg.phi_L1;
    eta_avg_R(1:200) = PN_eta(:,1:200,1)'*a_avg.eta_R1;
    theta_avg_R(1:200) = PN_theta(:,1:200,1)'*a_avg.theta_R1;
    phi_avg_R(1:200) = PN_phi(:,1:200,1)'*a_avg.phi_R1;
    theta_avg_LR(1:200) = PN_theta(:,1:200,1)'*a_avg.theta_LR1;
    eta_avg_LR(1:200) = PN_eta(:,1:200,1)'*a_avg.eta_LR1;
    phi_avg_LR(1:200) = PN_phi(:,1:200,1)'*a_avg.phi_LR1;

    eta_avg_L(201:401) = PN_eta(:,:,1)'*a_avg.eta_L2;
    theta_avg_L(201:401) = PN_theta(:,:,1)'*a_avg.theta_L2;
    phi_avg_L(201:401) = PN_phi(:,:,1)'*a_avg.phi_L2;
    eta_avg_R(201:401) = PN_eta(:,:,1)'*a_avg.eta_R2;
    theta_avg_R(201:401) = PN_theta(:,:,1)'*a_avg.theta_R2;
    phi_avg_R(201:401) = PN_phi(:,:,1)'*a_avg.phi_R2;
    theta_avg_LR(201:401) = PN_theta(:,:,1)'*a_avg.theta_LR2;
    eta_avg_LR(201:401) = PN_eta(:,:,1)'*a_avg.eta_LR2;
    phi_avg_LR(201:401) = PN_phi(:,:,1)'*a_avg.phi_LR2;

        
    
    figure()
    plot(radtodeg(theta_avg_L),'r')
    hold on
    plot(radtodeg(theta_avg_R),'g')
    plot(radtodeg(theta_avg_LR),'b')
    hold off        
        
    figure()
    plot(radtodeg(eta_avg_L),'r')
    hold on
    plot(radtodeg(eta_avg_R),'g')
    plot(radtodeg(eta_avg_LR),'b')
    hold off
    
    figure()
    plot(radtodeg(phi_avg_L),'r')
    hold on
    plot(radtodeg(phi_avg_R),'g')
    plot(radtodeg(phi_avg_LR),'b')
    hold off

    
    figure()
    hold on
    t_tot = 1/f_avg;
    for k = 3:trigger_wb
        t1_fit = (ratio_1_2(k)/2)*t_tot;
        t2_fit = ((2-ratio_1_2(k))/2)*t_tot;
        t_fit = [(0:(t1_fit/200):(t1_fit-(t1_fit/200))) (t1_fit:(t2_fit/200):t_tot)];
        plot(t_fit,radtodeg(theta_fit_L(k,:)),'Color',[0.5 0.5 0.5])
    end
    t1_avg = (mean(ratio_1_2(3:trigger_wb))/2)*t_tot;
    t2_avg = ((2-mean(ratio_1_2(3:trigger_wb)))/2)*t_tot;
    t_avg = [0:(t1_avg/200):(t1_avg-t1_avg/200) t1_avg:(t2_avg/200):t_tot];
    plot(t_avg,radtodeg(theta_avg_L),'r')
    hold off
    
    figure()
    hold on
    t_tot = 1/f_avg;
    for k = 3:trigger_wb
        t1_fit = (ratio_1_2(k)/2)*t_tot;
        t2_fit = ((2-ratio_1_2(k))/2)*t_tot;
        t_fit = [(0:(t1_fit/200):(t1_fit-(t1_fit/200))) (t1_fit:(t2_fit/200):t_tot)];
        plot(t_fit,radtodeg(eta_fit_L(k,:)),'Color',[0.5 0.5 0.5])
    end
    t1_avg = (mean(ratio_1_2(3:trigger_wb))/2)*t_tot;
    t2_avg = ((2-mean(ratio_1_2(3:trigger_wb)))/2)*t_tot;
    t_avg = [0:(t1_avg/200):(t1_avg-t1_avg/200) t1_avg:(t2_avg/200):t_tot];
    plot(t_avg,radtodeg(eta_avg_L),'r')
    hold off

    figure()
    hold on
    t_tot = 1/f_avg;
    for k = 3:trigger_wb
        t1_fit = (ratio_1_2(k)/2)*t_tot;
        t2_fit = ((2-ratio_1_2(k))/2)*t_tot;
        t_fit = [(0:(t1_fit/200):(t1_fit-(t1_fit/200))) (t1_fit:(t2_fit/200):t_tot)];
        plot(t_fit,radtodeg(phi_fit_L(k,:)),'Color',[0.5 0.5 0.5])
    end
    t1_avg = (mean(ratio_1_2(3:trigger_wb))/2)*t_tot;
    t2_avg = ((2-mean(ratio_1_2(3:trigger_wb)))/2)*t_tot;
    t_avg = [0:(t1_avg/200):(t1_avg-t1_avg/200) t1_avg:(t2_avg/200):t_tot];
    plot(t_avg,radtodeg(phi_avg_L),'r')
    hold off

    figure()
    hold on
    t_tot = 1/f_avg;
    for k = 3:trigger_wb
        t1_fit = (ratio_1_2(k)/2)*t_tot;
        t2_fit = ((2-ratio_1_2(k))/2)*t_tot;
        t_fit = [(0:(t1_fit/200):(t1_fit-(t1_fit/200))) (t1_fit:(t2_fit/200):t_tot)];
        plot(t_fit,radtodeg(theta_fit_R(k,:)),'Color',[0.5 0.5 0.5])
    end
    t1_avg = (mean(ratio_1_2(3:trigger_wb))/2)*t_tot;
    t2_avg = ((2-mean(ratio_1_2(3:trigger_wb)))/2)*t_tot;
    t_avg = [0:(t1_avg/200):(t1_avg-t1_avg/200) t1_avg:(t2_avg/200):t_tot];
    plot(t_avg,radtodeg(theta_avg_R),'r')
    hold off
    
    figure()
    hold on
    t_tot = 1/f_avg;
    for k = 3:trigger_wb
        t1_fit = (ratio_1_2(k)/2)*t_tot;
        t2_fit = ((2-ratio_1_2(k))/2)*t_tot;
        t_fit = [(0:(t1_fit/200):(t1_fit-(t1_fit/200))) (t1_fit:(t2_fit/200):t_tot)];
        plot(t_fit,radtodeg(eta_fit_R(k,:)),'Color',[0.5 0.5 0.5])
    end
    t1_avg = (mean(ratio_1_2(3:trigger_wb))/2)*t_tot;
    t2_avg = ((2-mean(ratio_1_2(3:trigger_wb)))/2)*t_tot;
    t_avg = [0:(t1_avg/200):(t1_avg-t1_avg/200) t1_avg:(t2_avg/200):t_tot];
    plot(t_avg,radtodeg(eta_avg_R),'r')
    hold off

    figure()
    hold on
    t_tot = 1/f_avg;
    for k = 3:trigger_wb
        t1_fit = (ratio_1_2(k)/2)*t_tot;
        t2_fit = ((2-ratio_1_2(k))/2)*t_tot;
        t_fit = [(0:(t1_fit/200):(t1_fit-(t1_fit/200))) (t1_fit:(t2_fit/200):t_tot)];
        plot(t_fit,radtodeg(phi_fit_R(k,:)),'Color',[0.5 0.5 0.5])
    end
    t1_avg = (mean(ratio_1_2(3:trigger_wb))/2)*t_tot;
    t2_avg = ((2-mean(ratio_1_2(3:trigger_wb)))/2)*t_tot;
    t_avg = [0:(t1_avg/200):(t1_avg-t1_avg/200) t1_avg:(t2_avg/200):t_tot];
    plot(t_avg,radtodeg(phi_avg_R),'r')
    hold off
    
    
    figure()
    hold on
    for k = 3:trigger_wb
        plot(-1:(1/200):1,radtodeg(theta_fit_R(k,:)),'Color',[0.5 0.5 0.5])
    end
    plot(-1:(1/200):1,radtodeg(theta_avg_R),'r')
    hold off
    
    figure()
    hold on
    for k = 3:trigger_wb
        plot(-1:(1/200):1,radtodeg(eta_fit_L(k,:)),'Color',[0.5 0.5 0.5])
    end
    plot(-1:(1/200):1,radtodeg(eta_avg_L),'r')
    hold off
    
    figure()
    hold on
    for k = 3:trigger_wb
        plot(-1:(1/200):1,radtodeg(eta_fit_R(k,:)),'Color',[0.5 0.5 0.5])
    end
    plot(-1:(1/200):1,radtodeg(eta_avg_R),'r')
    hold off
    
    figure()
    hold on
    for k = 3:trigger_wb
        plot(-1:(1/200):1,radtodeg(phi_fit_L(k,:)),'Color',[0.5 0.5 0.5])
    end
    plot(-1:(1/200):1,radtodeg(phi_avg_L),'r')
    hold off
    
    figure()
    hold on
    for k = 3:trigger_wb
        plot(-1:(1/200):1,radtodeg(phi_fit_R(k,:)),'Color',[0.5 0.5 0.5])
    end
    plot(-1:(1/200):1,radtodeg(phi_avg_R),'r')
    hold off
    
    end
    
%     figure()
%     plot(pathDB.omegax_body_mean(1:nr_wb,1,seq_nr))
%     
%     figure()
%     plot(pathDB.omegay_body_mean(1:nr_wb,1,seq_nr))
%     
%     figure()
%     plot(pathDB.omegaz_body_mean(1:nr_wb,1,seq_nr))
% 
% 
%     start = find(isnan(pathDB.x(:,seq_nr))==0, 1 );
%     stop = find(isnan(pathDB.x(:,seq_nr))==0, 1, 'last' );
% 
%     omega_x_body = pathDB.b_omega1(start:stop,seq_nr);
%     omega_y_body = pathDB.b_omega2(start:stop,seq_nr);
%     omega_z_body = pathDB.b_omega3(start:stop,seq_nr);
%     
%     % Create omega_x, omega_y and omega_z for the strokeplane reference
%     % frame:
%     
%     omega_x_strkpln = cosd(-55).*omega_x_body-sind(-55).*omega_z_body;
%     omega_y_strkpln = omega_y_body;
%     omega_z_strkpln = sind(-55).*omega_x_body+cosd(-55).*omega_z_body;
%     
%     
%     figure()
%     subplot(3,1,1); plot(omega_x_strkpln)
%     subplot(3,1,2); plot(omega_y_strkpln)
%     subplot(3,1,3); plot(omega_z_strkpln)
%     
%     Maneuver_wingbeat(a_fit_L,a_fit_R,a_avg_L,a_avg_R,a_avg_LR,n_pol,trigger_wb)
    
    
    pause
    
    close all
    
    end
    

    
    
    
    
    

end
