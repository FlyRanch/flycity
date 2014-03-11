function [ a_fit, a_avg, a_dev ] = Standard_Wingbeat( settings, pathDB, seq_nr )

    % Obtain local polynomial fit and sequence averaged fit:
    
    n_pol_theta     = settings.n_pol_theta;
    n_pol_eta       = settings.n_pol_eta;
    n_pol_phi       = settings.n_pol_phi;
    
    start   = settings.start_stop(seq_nr,1);
    stop    = settings.start_stop(seq_nr,2);
    
    t = pathDB.t;
    dt = pathDB.dt;

    theta_L     = pathDB.wing_kin.theta_L(seq_nr,:);
    eta_L       = pathDB.wing_kin.eta_L(seq_nr,:);
    phi_L       = pathDB.wing_kin.phi_L(seq_nr,:);
    theta_R     = pathDB.wing_kin.theta_R(seq_nr,:);
    eta_R       = pathDB.wing_kin.eta_R(seq_nr,:);
    phi_R       = pathDB.wing_kin.phi_R(seq_nr,:);
    
    nr_wb = pathDB.wingbeats.nr_of_wb(seq_nr);
    
    wingbeat_loc    = pathDB.wingbeats.wingbeat_loc(1:nr_wb,:,seq_nr);
    downstroke_loc  = pathDB.wingbeats.downstroke_loc(1:nr_wb,:,seq_nr);
    upstroke_loc    = pathDB.wingbeats.upstroke_loc(1:nr_wb,:,seq_nr);
    
    a_fit.theta_L   = nan(2*(n_pol_theta+1),200);
    a_fit.eta_L     = nan(2*(n_pol_eta+1),200);
    a_fit.phi_L     = nan(2*(n_pol_phi+1),200);
    a_fit.theta_R   = nan(2*(n_pol_theta+1),200);
    a_fit.eta_R     = nan(2*(n_pol_eta+1),200);
    a_fit.phi_R     = nan(2*(n_pol_phi+1),200);
    
       
    % Determine downstroke upstroke ratio:
    
    a_fit.down_up = (downstroke_loc(:,2)-downstroke_loc(:,1)+1)./(wingbeat_loc(:,2)-wingbeat_loc(:,1)+1);
    
    % Determine the wingbeat frequency:
    
    a_fit.f = 1./((wingbeat_loc(:,2)-wingbeat_loc(:,1)+1)*dt);
    
    wb_loc = nan(nr_wb*2,2);
    down_up_ratio = nan(nr_wb*2,1);
    f_piece = nan(nr_wb*2,1);
    
    for k = 1:nr_wb
        
        wb_loc(2*k-1,:)         = [ downstroke_loc(k,1) downstroke_loc(k,2)+1 ];
        wb_loc(2*k,:)           = [ upstroke_loc(k,1) upstroke_loc(k,2)+1 ];
        down_up_ratio(2*k-1)    = a_fit.down_up(k);
        down_up_ratio(2*k)      = 1-a_fit.down_up(k);
        
    end
    
    a_theta_L = Piecewise_polynomial_fit(theta_L',n_pol_theta,wb_loc,down_up_ratio);
    a_eta_L = Piecewise_polynomial_fit(eta_L',n_pol_eta,wb_loc,down_up_ratio);
    a_phi_L = Piecewise_polynomial_fit(phi_L',n_pol_phi,wb_loc,down_up_ratio);
    a_theta_R = Piecewise_polynomial_fit(theta_R',n_pol_theta,wb_loc,down_up_ratio);
    a_eta_R = Piecewise_polynomial_fit(eta_R',n_pol_eta,wb_loc,down_up_ratio);
    a_phi_R = Piecewise_polynomial_fit(phi_R',n_pol_phi,wb_loc,down_up_ratio);
    
    a_fit.theta_L(:,1:nr_wb)   = [a_theta_L(:,1:2:end); a_theta_L(:,2:2:end)];
    a_fit.eta_L(:,1:nr_wb)     = [a_eta_L(:,1:2:end); a_eta_L(:,2:2:end)];
    a_fit.phi_L(:,1:nr_wb)     = [a_phi_L(:,1:2:end); a_phi_L(:,2:2:end)];
    a_fit.theta_R(:,1:nr_wb)   = [a_theta_R(:,1:2:end); a_theta_R(:,2:2:end)];
    a_fit.eta_R(:,1:nr_wb)     = [a_eta_R(:,1:2:end); a_eta_R(:,2:2:end)];
    a_fit.phi_R(:,1:nr_wb)     = [a_phi_R(:,1:2:end); a_phi_R(:,2:2:end)];
    
    
    frames              = settings.frame_end;
    
    theta_L_fit         = nan(1,frames);
    eta_L_fit           = nan(1,frames);
    phi_L_fit           = nan(1,frames);
    
    theta_R_fit         = nan(1,frames);
    eta_R_fit           = nan(1,frames);
    phi_R_fit           = nan(1,frames);
    
    theta_dot_L_fit     = nan(1,frames);
    eta_dot_L_fit       = nan(1,frames);
    phi_dot_L_fit       = nan(1,frames);
    
    theta_dot_R_fit     = nan(1,frames);
    eta_dot_R_fit       = nan(1,frames);
    phi_dot_R_fit       = nan(1,frames);
    
    theta_ddot_L_fit    = nan(1,frames);
    eta_ddot_L_fit      = nan(1,frames);
    phi_ddot_L_fit      = nan(1,frames);
    
    theta_ddot_R_fit    = nan(1,frames);
    eta_ddot_R_fit      = nan(1,frames);
    phi_ddot_R_fit      = nan(1,frames);
    
    for k = 1:nr_wb
        
        [ ~, X_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, a_fit.down_up(k), wingbeat_loc(k,2)-wingbeat_loc(k,1)+2,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),0);
        [ ~,    X_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, a_fit.down_up(k), wingbeat_loc(k,2)-wingbeat_loc(k,1)+2,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),0);
        [ ~, X_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, a_fit.down_up(k), wingbeat_loc(k,2)-wingbeat_loc(k,1)+2,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),0);
        
        [ ~, X_theta_dot ] = Wingbeat_Legendre_matrix( n_pol_theta, a_fit.down_up(k), wingbeat_loc(k,2)-wingbeat_loc(k,1)+2,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),1);
        [ ~,    X_eta_dot ] = Wingbeat_Legendre_matrix( n_pol_eta, a_fit.down_up(k), wingbeat_loc(k,2)-wingbeat_loc(k,1)+2,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),1);
        [ ~, X_phi_dot ] = Wingbeat_Legendre_matrix( n_pol_phi, a_fit.down_up(k), wingbeat_loc(k,2)-wingbeat_loc(k,1)+2,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),1);
        
        [ ~, X_theta_ddot ] = Wingbeat_Legendre_matrix( n_pol_theta, a_fit.down_up(k), wingbeat_loc(k,2)-wingbeat_loc(k,1)+2,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),2);
        [ ~,    X_eta_ddot ] = Wingbeat_Legendre_matrix( n_pol_eta, a_fit.down_up(k), wingbeat_loc(k,2)-wingbeat_loc(k,1)+2,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),2);
        [ ~, X_phi_ddot ] = Wingbeat_Legendre_matrix( n_pol_phi, a_fit.down_up(k), wingbeat_loc(k,2)-wingbeat_loc(k,1)+2,t(wingbeat_loc(k,1)),t(wingbeat_loc(k,2)+1),2);
        
        theta_L_fit(wingbeat_loc(k,1):(wingbeat_loc(k,2)+1))         = X_theta*a_fit.theta_L(:,k);
        eta_L_fit(wingbeat_loc(k,1):(wingbeat_loc(k,2)+1))           = X_eta*a_fit.eta_L(:,k);
        phi_L_fit(wingbeat_loc(k,1):(wingbeat_loc(k,2)+1))           = X_phi*a_fit.phi_L(:,k);

        theta_R_fit(wingbeat_loc(k,1):(wingbeat_loc(k,2)+1))         = X_theta*a_fit.theta_R(:,k);
        eta_R_fit(wingbeat_loc(k,1):(wingbeat_loc(k,2)+1))           = X_eta*a_fit.eta_R(:,k);
        phi_R_fit(wingbeat_loc(k,1):(wingbeat_loc(k,2)+1))           = X_phi*a_fit.phi_R(:,k);

        theta_dot_L_fit(wingbeat_loc(k,1):(wingbeat_loc(k,2)+1))     = (X_theta_dot*a_fit.theta_L(:,k)).*(2*a_fit.f(k));
        eta_dot_L_fit(wingbeat_loc(k,1):(wingbeat_loc(k,2)+1))       = (X_eta_dot*a_fit.eta_L(:,k)).*(2*a_fit.f(k));
        phi_dot_L_fit(wingbeat_loc(k,1):(wingbeat_loc(k,2)+1))       = (X_phi_dot*a_fit.phi_L(:,k)).*(2*a_fit.f(k));

        theta_dot_R_fit(wingbeat_loc(k,1):(wingbeat_loc(k,2)+1))     = (X_theta_dot*a_fit.theta_R(:,k)).*(2*a_fit.f(k));
        eta_dot_R_fit(wingbeat_loc(k,1):(wingbeat_loc(k,2)+1))       = (X_eta_dot*a_fit.eta_R(:,k)).*(2*a_fit.f(k));
        phi_dot_R_fit(wingbeat_loc(k,1):(wingbeat_loc(k,2)+1))       = (X_phi_dot*a_fit.phi_R(:,k)).*(2*a_fit.f(k));

        theta_ddot_L_fit(wingbeat_loc(k,1):(wingbeat_loc(k,2)+1))    = (X_theta_ddot*a_fit.theta_L(:,k)).*(2*a_fit.f(k))^2;
        eta_ddot_L_fit(wingbeat_loc(k,1):(wingbeat_loc(k,2)+1))      = (X_eta_ddot*a_fit.eta_L(:,k)).*(2*a_fit.f(k))^2;
        phi_ddot_L_fit(wingbeat_loc(k,1):(wingbeat_loc(k,2)+1))      = (X_phi_ddot*a_fit.phi_L(:,k)).*(2*a_fit.f(k))^2;

        theta_ddot_R_fit(wingbeat_loc(k,1):(wingbeat_loc(k,2)+1))    = (X_theta_ddot*a_fit.theta_R(:,k)).*(2*a_fit.f(k))^2;
        eta_ddot_R_fit(wingbeat_loc(k,1):(wingbeat_loc(k,2)+1))      = (X_eta_ddot*a_fit.eta_R(:,k)).*(2*a_fit.f(k))^2;
        phi_ddot_R_fit(wingbeat_loc(k,1):(wingbeat_loc(k,2)+1))      = (X_phi_ddot*a_fit.phi_R(:,k)).*(2*a_fit.f(k))^2;
    
    end
    
    theta_dot_L_dt      = nan(1,frames);
    eta_dot_L_dt        = nan(1,frames);
    phi_dot_L_dt        = nan(1,frames);
    
    theta_dot_R_dt      = nan(1,frames);
    eta_dot_R_dt        = nan(1,frames);
    phi_dot_R_dt        = nan(1,frames);
    
    theta_ddot_L_dt     = nan(1,frames);
    eta_ddot_L_dt       = nan(1,frames);
    phi_ddot_L_dt       = nan(1,frames);
    
    theta_ddot_R_dt     = nan(1,frames);
    eta_ddot_R_dt       = nan(1,frames);
    phi_ddot_R_dt       = nan(1,frames);
    
    for i = wingbeat_loc(1,1):wingbeat_loc(end,2)
        
        if i == wingbeat_loc(1,1)
        
            theta_dot_L_dt(i)    = (theta_L_fit(i+1)-theta_L_fit(i))/dt;
            eta_dot_L_dt(i)      = (eta_L_fit(i+1)-eta_L_fit(i))/dt;
            phi_dot_L_dt(i)      = (phi_L_fit(i+1)-phi_L_fit(i))/dt;

            theta_dot_R_dt(i)    = (theta_R_fit(i+1)-theta_R_fit(i))/dt;
            eta_dot_R_dt(i)      = (eta_R_fit(i+1)-eta_R_fit(i))/dt;
            phi_dot_R_dt(i)      = (phi_R_fit(i+1)-phi_R_fit(i))/dt;

            theta_ddot_L_dt(i)  = (theta_L_fit(i+2)-2*theta_L_fit(i+1)+theta_L_fit(i))/(dt^2);
            eta_ddot_L_dt(i)    = (eta_L_fit(i+2)-2*eta_L_fit(i+1)+eta_L_fit(i))/(dt^2);
            phi_ddot_L_dt(i)    = (phi_L_fit(i+2)-2*phi_L_fit(i+1)+phi_L_fit(i))/(dt^2);

            theta_ddot_R_dt(i)  = (theta_R_fit(i+2)-2*theta_R_fit(i+1)+theta_R_fit(i))/(dt^2);
            eta_ddot_R_dt(i)    = (eta_R_fit(i+2)-2*eta_R_fit(i+1)+eta_R_fit(i))/(dt^2);
            phi_ddot_R_dt(i)    = (phi_R_fit(i+2)-2*phi_R_fit(i+1)+phi_R_fit(i))/(dt^2);
        
        elseif i == wingbeat_loc(end,2)

            theta_dot_L_dt(i)    = (theta_L_fit(i)-theta_L_fit(i-1))/dt;
            eta_dot_L_dt(i)      = (eta_L_fit(i)-eta_L_fit(i-1))/dt;
            phi_dot_L_dt(i)      = (phi_L_fit(i)-phi_L_fit(i-1))/dt;

            theta_dot_R_dt(i)    = (theta_R_fit(i)-theta_R_fit(i-1))/dt;
            eta_dot_R_dt(i)      = (eta_R_fit(i)-eta_R_fit(i-1))/dt;
            phi_dot_R_dt(i)      = (phi_R_fit(i)-phi_R_fit(i-1))/dt;

            theta_ddot_L_dt(i)  = (theta_L_fit(i)-2*theta_L_fit(i-1)+theta_L_fit(i-2))/(dt^2);
            eta_ddot_L_dt(i)    = (eta_L_fit(i)-2*eta_L_fit(i-1)+eta_L_fit(i-2))/(dt^2);
            phi_ddot_L_dt(i)    = (phi_L_fit(i)-2*phi_L_fit(i-1)+phi_L_fit(i-2))/(dt^2);

            theta_ddot_R_dt(i)  = (theta_R_fit(i)-2*theta_R_fit(i-1)+theta_R_fit(i-2))/(dt^2);
            eta_ddot_R_dt(i)    = (eta_R_fit(i)-2*eta_R_fit(i-1)+eta_R_fit(i-2))/(dt^2);
            phi_ddot_R_dt(i)    = (phi_R_fit(i)-2*phi_R_fit(i-1)+phi_R_fit(i-2))/(dt^2);
            
        else
            
            theta_dot_L_dt(i)    = (theta_L_fit(i+1)-theta_L_fit(i-1))/(2*dt);
            eta_dot_L_dt(i)      = (eta_L_fit(i+1)-eta_L_fit(i-1))/(2*dt);
            phi_dot_L_dt(i)      = (phi_L_fit(i+1)-phi_L_fit(i-1))/(2*dt);

            theta_dot_R_dt(i)    = (theta_R_fit(i+1)-theta_R_fit(i-1))/(2*dt);
            eta_dot_R_dt(i)      = (eta_R_fit(i+1)-eta_R_fit(i-1))/(2*dt);
            phi_dot_R_dt(i)      = (phi_R_fit(i+1)-phi_R_fit(i-1))/(2*dt);

            theta_ddot_L_dt(i)  = (theta_L_fit(i+1)-2*theta_L_fit(i)+theta_L_fit(i-1))/(dt^2);
            eta_ddot_L_dt(i)    = (eta_L_fit(i+1)-2*eta_L_fit(i)+eta_L_fit(i-1))/(dt^2);
            phi_ddot_L_dt(i)    = (phi_L_fit(i+1)-2*phi_L_fit(i)+phi_L_fit(i-1))/(dt^2);

            theta_ddot_R_dt(i)  = (theta_R_fit(i+1)-2*theta_R_fit(i)+theta_R_fit(i-1))/(dt^2);
            eta_ddot_R_dt(i)    = (eta_R_fit(i+1)-2*eta_R_fit(i)+eta_R_fit(i-1))/(dt^2);
            phi_ddot_R_dt(i)    = (phi_R_fit(i+1)-2*phi_R_fit(i)+phi_R_fit(i-1))/(dt^2);
            
        end
        
    end

 
    
%     figure()
%     hold on
%     plot(t,theta_L,'b')
%     plot(t,theta_L_fit,'r')
%     hold off
%     title('Piecewise polynomial fit \theta_L')
%     xlabel('Time [s]')
%     ylabel('\theta_L [rad]')
%     legend('raw','fit')
%     
%     figure()
%     hold on
%     plot(t,eta_L,'b')
%     plot(t,eta_L_fit,'r')
%     hold off
%     title('Piecewise polynomial fit \eta_L')
%     xlabel('Time [s]')
%     ylabel('\eta_L [rad]')
%     legend('raw','fit')
%     
%     figure()
%     hold on
%     plot(t,phi_L,'b')
%     plot(t,phi_L_fit,'r')
%     hold off
%     title('Piecewise polynomial fit \phi_L')
%     xlabel('Time [s]')
%     ylabel('\phi_L [rad]')
%     legend('raw','fit')
%     
%     figure()
%     hold on
%     plot(t,theta_R,'b')
%     plot(t,theta_R_fit,'g')
%     hold off
%     title('Piecewise polynomial fit \theta_R')
%     xlabel('Time [s]')
%     ylabel('\theta_R [rad]')
%     legend('raw','fit')
%     
%     figure()
%     hold on
%     plot(t,eta_R,'b')
%     plot(t,eta_R_fit,'g')
%     hold off
%     title('Piecewise polynomial fit \eta_R')
%     xlabel('Time [s]')
%     ylabel('\eta_R [rad]')
%     legend('raw','fit')
%     
%     figure()
%     hold on
%     plot(t,phi_R,'b')
%     plot(t,phi_R_fit,'g')
%     hold off
%     title('Piecewise polynomial fit \phi_R')
%     xlabel('Time [s]')
%     ylabel('\phi_R [rad]')
%     legend('raw','fit')
%     
%     figure()
%     hold on
%     plot(t,theta_dot_L_dt,'b')
%     plot(t,theta_dot_L_fit,'r')
%     hold off
%     title('Piecewise polynomial fit d\theta_L / dt')
%     xlabel('Time [s]')
%     ylabel('d\theta_L / dt [rad/s]')
%     legend('raw','fit')
%     
%     figure()
%     hold on
%     plot(t,eta_dot_L_dt,'b')
%     plot(t,eta_dot_L_fit,'r')
%     hold off
%     title('Piecewise polynomial fit d\eta_L / dt')
%     xlabel('Time [s]')
%     ylabel('d\eta_L / dt [rad/s]')
%     legend('raw','fit')
%     
%     figure()
%     hold on
%     plot(t,phi_dot_L_dt,'b')
%     plot(t,phi_dot_L_fit,'r')
%     hold off
%     title('Piecewise polynomial fit d\phi_L / dt')
%     xlabel('Time [s]')
%     ylabel('d\phi_L / dt [rad/s]')
%     legend('raw','fit')
%     
%     figure()
%     hold on
%     plot(t,theta_dot_R_dt,'b')
%     plot(t,theta_dot_R_fit,'g')
%     hold off
%     title('Piecewise polynomial fit d\theta_R / dt')
%     xlabel('Time [s]')
%     ylabel('d\theta_R / dt [rad/s]')
%     legend('raw','fit')
%     
%     figure()
%     hold on
%     plot(t,eta_dot_R_dt,'b')
%     plot(t,eta_dot_R_fit,'g')
%     hold off
%     title('Piecewise polynomial fit d\eta_R / dt')
%     xlabel('Time [s]')
%     ylabel('d\eta_R / dt [rad/s]')
%     legend('raw','fit')
%     
%     figure()
%     hold on
%     plot(t,phi_dot_R_dt,'b')
%     plot(t,phi_dot_R_fit,'g')
%     hold off
%     title('Piecewise polynomial fit d\phi_R / dt')
%     xlabel('Time [s]')
%     ylabel('d\phi_R / dt [rad/s]')
%     legend('raw','fit')
% 
%     figure()
%     hold on
%     plot(t,theta_ddot_L_dt,'b')
%     plot(t,theta_ddot_L_fit,'r')
%     hold off
%     title('Piecewise polynomial fit d^2\theta_L / dt^2')
%     xlabel('Time [s]')
%     ylabel('d^2\theta_L / dt^2 [rad/s^2]')
%     legend('raw','fit')
%     
%     figure()
%     hold on
%     plot(t,eta_ddot_L_dt,'b')
%     plot(t,eta_ddot_L_fit,'r')
%     hold off
%     title('Piecewise polynomial fit d^2\eta_L / dt^2')
%     xlabel('Time [s]')
%     ylabel('d^2\eta_L / dt^2 [rad/s^2]')
%     legend('raw','fit')
%     
%     figure()
%     hold on
%     plot(t,phi_ddot_L_dt,'b')
%     plot(t,phi_ddot_L_fit,'r')
%     hold off
%     title('Piecewise polynomial fit d^2\phi_L / dt^2')
%     xlabel('Time [s]')
%     ylabel('d^2\phi_L / dt^2 [rad/s^2]')
%     legend('raw','fit')
%     
%     figure()
%     hold on
%     plot(t,theta_ddot_R_dt,'b')
%     plot(t,theta_ddot_R_fit,'g')
%     hold off
%     title('Piecewise polynomial fit d^2\theta_R / dt^2')
%     xlabel('Time [s]')
%     ylabel('d^2\theta_R / dt^2 [rad/s^2]')
%     legend('raw','fit')
%     
%     figure()
%     hold on
%     plot(t,eta_ddot_R_dt,'b')
%     plot(t,eta_ddot_R_fit,'g')
%     hold off
%     title('Piecewise polynomial fit d^2\eta_R / dt^2')
%     xlabel('Time [s]')
%     ylabel('d^2\eta_R / dt^2 [rad/s^2]')
%     legend('raw','fit')
%     
%     figure()
%     hold on
%     plot(t,phi_ddot_R_dt,'b')
%     plot(t,phi_ddot_R_fit,'g')
%     hold off
%     title('Piecewise polynomial fit d^2\phi_R / dt^2')
%     xlabel('Time [s]')
%     ylabel('d^2\phi_R / dt^2 [rad/s^2]')
%     legend('raw','fit')
%     
%     pause


    % Find the average wingbeat:
    
    % find the wingbeat which includes the trigger:
    
    trigger_wb = 0;
    
    seq_trigger_shift = settings.trigger_shift;
    
    trig_shift = sum(strcmp(char(settings.sequence_names(seq_nr)),seq_trigger_shift)==1);
                          
    if trig_shift ~= 0
        trigger_frame = 2683;
    else
        trigger_frame = 2795;
    end
    
    for j = 1:nr_wb
        
        if wingbeat_loc(j,1)<=trigger_frame && wingbeat_loc(j,2)>trigger_frame
            
            trigger_wb = j;
            
        elseif wingbeat_loc(j,2) == trigger_frame
        
            trigger_wb = j;
            
        end
    end
    
    % Find the average wingbeat:
    
    a_avg.theta_L    = nan((n_pol_theta+1)*2,1);
    a_avg.eta_L      = nan((n_pol_eta+1)*2,1);
    a_avg.phi_L      = nan((n_pol_phi+1)*2,1);
    
    a_avg.theta_R    = nan((n_pol_theta+1)*2,1);
    a_avg.eta_R      = nan((n_pol_eta+1)*2,1);
    a_avg.phi_R      = nan((n_pol_phi+1)*2,1);
    
    a_avg.theta_LR   = nan((n_pol_theta+1)*2,1);
    a_avg.eta_LR     = nan((n_pol_eta+1)*2,1);
    a_avg.phi_LR     = nan((n_pol_phi+1)*2,1);
    
    a_avg.down_up    = nan;
    a_avg.f          = nan;
    a_avg.trigger_wb = trigger_wb;
    a_avg.m_fly      = nan;
    a_avg.wing_l     = nan;
    
    a_dev.theta_L = nan(2*(n_pol_theta+1),200);
    a_dev.eta_L = nan(2*(n_pol_eta+1),200);
    a_dev.phi_L = nan(2*(n_pol_phi+1),200);
    
    a_dev.theta_R = nan(2*(n_pol_theta+1),200);
    a_dev.eta_R = nan(2*(n_pol_eta+1),200);
    a_dev.phi_R = nan(2*(n_pol_phi+1),200);
    
    if trigger_wb >= 5  % minimum of 3 wingbeats in order to determine the average
        
        a_avg.down_up   = mean(a_fit.down_up(3:trigger_wb));  % exclude the first 2 wingbeats from analysis
    
        a_avg.f         = mean(a_fit.f(3:trigger_wb));  % exclude the first 2 wingbeats from analysis
    
        a_avg.m_fly     = pathDB.body_model.mass_fly(seq_nr);
        
        a_avg.wing_l    = pathDB.wing_model.length(seq_nr);
        
        % Determine weights for the weighted average:
        
        [ W_theta_L1, W_theta_L2, W_theta_R1, W_theta_R2, W_theta_LR1, W_theta_LR2 ] = weight(  a_fit.theta_L(1:(n_pol_theta+1),1:nr_wb), ...
                                                                                                a_fit.theta_L((n_pol_theta+2):(2*(n_pol_theta+1)),1:nr_wb), ...
                                                                                                a_fit.theta_R(1:(n_pol_theta+1),1:nr_wb), ... 
                                                                                                a_fit.theta_R((n_pol_theta+2):(2*(n_pol_theta+1)),1:nr_wb), ... 
                                                                                                n_pol_theta, trigger_wb);
        
        [ W_eta_L1, W_eta_L2, W_eta_R1, W_eta_R2, W_eta_LR1, W_eta_LR2 ] = weight(  a_fit.eta_L(1:(n_pol_eta+1),1:nr_wb), ... 
                                                                                    a_fit.eta_L((n_pol_eta+2):(2*(n_pol_eta+1)),1:nr_wb), ... 
                                                                                    a_fit.eta_R(1:(n_pol_eta+1),1:nr_wb), ... 
                                                                                    a_fit.eta_R((n_pol_eta+2):(2*(n_pol_eta+1)),1:nr_wb), ...
                                                                                    n_pol_eta, trigger_wb);
                                                                                
        [ W_phi_L1, W_phi_L2, W_phi_R1, W_phi_R2, W_phi_LR1, W_phi_LR2 ] = weight(  a_fit.phi_L(1:(n_pol_phi+1),1:nr_wb), ... 
                                                                                    a_fit.phi_L((n_pol_phi+2):(2*(n_pol_phi+1)),1:nr_wb), ...
                                                                                    a_fit.phi_R(1:(n_pol_phi+1),1:nr_wb), ...
                                                                                    a_fit.phi_R((n_pol_phi+2):(2*(n_pol_phi+1)),1:nr_wb), ...
                                                                                    n_pol_phi, trigger_wb);
        
        % Use function average_fit to obtain the average for the left and
        % right wing seperately:

        [a_avg_theta_L1, a_avg_theta_L2] = average_fit( a_fit.theta_L((1:(n_pol_theta+1)),3:trigger_wb), ...
                                                        a_fit.theta_L((n_pol_theta+2):(2*(n_pol_theta+1)),3:trigger_wb), ...
                                                        n_pol_theta,W_theta_L1,W_theta_L2,a_avg.down_up);
                                                    
        [a_avg_theta_R1, a_avg_theta_R2] = average_fit( a_fit.theta_R((1:(n_pol_theta+1)),3:trigger_wb), ...
                                                        a_fit.theta_R((n_pol_theta+2):(2*(n_pol_theta+1)),3:trigger_wb), ...
                                                        n_pol_theta,W_theta_R1,W_theta_R2,a_avg.down_up);
                                                    
        [a_avg_eta_L1, a_avg_eta_L2] = average_fit( a_fit.eta_L(1:(n_pol_eta+1),3:trigger_wb), ...
                                                    a_fit.eta_L((n_pol_eta+2):(2*(n_pol_eta+1)),3:trigger_wb), ...
                                                    n_pol_eta,W_eta_L1,W_eta_L2,a_avg.down_up);
                                                
        [a_avg_eta_R1, a_avg_eta_R2] = average_fit( a_fit.eta_R(1:(n_pol_eta+1),3:trigger_wb), ...
                                                    a_fit.eta_R((n_pol_eta+2):(2*(n_pol_eta+1)),3:trigger_wb), ...
                                                    n_pol_eta,W_eta_R1,W_eta_R2,a_avg.down_up);
                                                
        [a_avg_phi_L1, a_avg_phi_L2] = average_fit( a_fit.phi_L(1:(n_pol_phi+1),3:trigger_wb), ...
                                                    a_fit.phi_L((n_pol_phi+2):(2*(n_pol_phi+1)),3:trigger_wb), ...
                                                    n_pol_phi,W_phi_L1,W_phi_L2,a_avg.down_up);
                                                
        [a_avg_phi_R1, a_avg_phi_R2] = average_fit( a_fit.phi_R(1:(n_pol_phi+1),3:trigger_wb), ...
                                                    a_fit.phi_R((n_pol_phi+2):(2*(n_pol_phi+1)),3:trigger_wb), ...
                                                    n_pol_phi,W_phi_R1,W_phi_R2,a_avg.down_up);     

        % Use function average_fit to obtain the average for the left and
        % right wing combined:

        [a_avg_theta_LR1, a_avg_theta_LR2] = average_fit(   [a_fit.theta_L((1:(n_pol_theta+1)),3:trigger_wb) ...
                                                             a_fit.theta_R((1:(n_pol_theta+1)),3:trigger_wb)], ...
                                                            [a_fit.theta_L((n_pol_theta+2):(2*(n_pol_theta+1)),3:trigger_wb) ... 
                                                             a_fit.theta_R((n_pol_theta+2):(2*(n_pol_theta+1)),3:trigger_wb)], ...
                                                             n_pol_theta,W_theta_LR1,W_theta_LR2,a_avg.down_up);
                                                        
        [a_avg_eta_LR1, a_avg_eta_LR2] = average_fit(   [a_fit.eta_L((1:(n_pol_eta+1)),3:trigger_wb) ...
                                                         a_fit.eta_R((1:(n_pol_eta+1)),3:trigger_wb)], ...
                                                        [a_fit.eta_L((n_pol_eta+2):(2*(n_pol_eta+1)),3:trigger_wb) ... 
                                                         a_fit.eta_R((n_pol_eta+2):(2*(n_pol_eta+1)),3:trigger_wb)], ...
                                                         n_pol_eta,W_eta_LR1,W_eta_LR2,a_avg.down_up);
                                                         
        [a_avg_phi_LR1, a_avg_phi_LR2] = average_fit(   [a_fit.phi_L((1:(n_pol_phi+1)),3:trigger_wb) ...
                                                         a_fit.phi_R((1:(n_pol_phi+1)),3:trigger_wb)], ...
                                                        [a_fit.phi_L((n_pol_phi+2):(2*(n_pol_phi+1)),3:trigger_wb) ... 
                                                         a_fit.phi_R((n_pol_phi+2):(2*(n_pol_phi+1)),3:trigger_wb)], ...
                                                         n_pol_phi,W_phi_LR1,W_phi_LR2,a_avg.down_up);

        
        % Store average coefficients in a_avg
        
        a_avg.theta_L   = [a_avg_theta_L1; a_avg_theta_L2];
        a_avg.eta_L     = [a_avg_eta_L1; a_avg_eta_L2];
        a_avg.phi_L     = [a_avg_phi_L1; a_avg_phi_L2];

        a_avg.theta_R   = [a_avg_theta_R1; a_avg_theta_R2];
        a_avg.eta_R     = [a_avg_eta_R1; a_avg_eta_R2];
        a_avg.phi_R     = [a_avg_phi_R1; a_avg_phi_R2];     
        
        a_avg.theta_LR = [a_avg_theta_LR1; a_avg_theta_LR2];
        a_avg.eta_LR   = [a_avg_eta_LR1; a_avg_eta_LR2];
        a_avg.phi_LR   = [a_avg_phi_LR1; a_avg_phi_LR2];
        
        % Determine the deviation from the average wingbeat for the left and
        % right wings seperately:
        
        for k = 1:nr_wb

            a_dev.theta_L(:,k)   = a_fit.theta_L(:,k)-a_avg.theta_L;
            a_dev.eta_L(:,k)     = a_fit.eta_L(:,k)-a_avg.eta_L;
            a_dev.phi_L(:,k)     = a_fit.phi_L(:,k)-a_avg.phi_L;

            a_dev.theta_R(:,k)   = a_fit.theta_R(:,k)-a_avg.theta_R;
            a_dev.eta_R(:,k)     = a_fit.eta_R(:,k)-a_avg.eta_R;
            a_dev.phi_R(:,k)     = a_fit.phi_R(:,k)-a_avg.phi_R;
            
        end

%             [ t_avg , X_theta_avg ] = Wingbeat_Legendre_matrix( n_pol_theta, a_avg.down_up,1000,0,1/a_avg.f,0);
%             [ ~,    X_eta_avg ] = Wingbeat_Legendre_matrix( n_pol_eta, a_avg.down_up,1000,0,1/a_avg.f,0);
%             [ ~, X_phi_avg ] = Wingbeat_Legendre_matrix( n_pol_phi, a_avg.down_up,1000,0,1/a_avg.f,0);
% 
%             [ ~, X_theta_dot_avg ] = Wingbeat_Legendre_matrix( n_pol_theta, a_avg.down_up, 1000,0,1/a_avg.f,1);
%             [ ~,    X_eta_dot_avg ] = Wingbeat_Legendre_matrix( n_pol_eta, a_avg.down_up, 1000,0,1/a_avg.f,1);
%             [ ~, X_phi_dot_avg ] = Wingbeat_Legendre_matrix( n_pol_phi, a_avg.down_up, 1000,0,1/a_avg.f,1);
% 
%             [ ~, X_theta_ddot_avg ] = Wingbeat_Legendre_matrix( n_pol_theta, a_avg.down_up, 1000,0,1/a_avg.f,2);
%             [ ~,    X_eta_ddot_avg ] = Wingbeat_Legendre_matrix( n_pol_eta, a_avg.down_up, 1000,0,1/a_avg.f,2);
%             [ ~, X_phi_ddot_avg ] = Wingbeat_Legendre_matrix( n_pol_phi, a_avg.down_up, 1000,0,1/a_avg.f,2);
% 
%             theta_L_avg         = X_theta_avg*a_avg.theta_L;
%             eta_L_avg           = X_eta_avg*a_avg.eta_L;
%             phi_L_avg           = X_phi_avg*a_avg.phi_L;
% 
%             theta_R_avg         = X_theta_avg*a_avg.theta_R;
%             eta_R_avg           = X_eta_avg*a_avg.eta_R;
%             phi_R_avg           = X_phi_avg*a_avg.phi_R;
%             
%             theta_LR_avg         = X_theta_avg*a_avg.theta_LR;
%             eta_LR_avg           = X_eta_avg*a_avg.eta_LR;
%             phi_LR_avg           = X_phi_avg*a_avg.phi_LR;
% 
%             theta_dot_L_avg     = (X_theta_dot_avg*a_avg.theta_L).*(2*a_avg.f);
%             eta_dot_L_avg       = (X_eta_dot_avg*a_avg.eta_L).*(2*a_avg.f);
%             phi_dot_L_avg       = (X_phi_dot_avg*a_avg.phi_L).*(2*a_avg.f);
% 
%             theta_dot_R_avg     = (X_theta_dot_avg*a_avg.theta_R).*(2*a_avg.f);
%             eta_dot_R_avg       = (X_eta_dot_avg*a_avg.eta_R).*(2*a_avg.f);
%             phi_dot_R_avg       = (X_phi_dot_avg*a_avg.phi_R).*(2*a_avg.f);
%             
%             theta_dot_LR_avg     = (X_theta_dot_avg*a_avg.theta_LR).*(2*a_avg.f);
%             eta_dot_LR_avg       = (X_eta_dot_avg*a_avg.eta_LR).*(2*a_avg.f);
%             phi_dot_LR_avg       = (X_phi_dot_avg*a_avg.phi_LR).*(2*a_avg.f);
% 
%             theta_ddot_L_avg    = (X_theta_ddot_avg*a_avg.theta_L).*(2*a_avg.f)^2;
%             eta_ddot_L_avg      = (X_eta_ddot_avg*a_avg.eta_L).*(2*a_avg.f)^2;
%             phi_ddot_L_avg      = (X_phi_ddot_avg*a_avg.phi_L).*(2*a_avg.f)^2;
% 
%             theta_ddot_R_avg    = (X_theta_ddot_avg*a_avg.theta_R).*(2*a_avg.f)^2;
%             eta_ddot_R_avg      = (X_eta_ddot_avg*a_avg.eta_R).*(2*a_avg.f)^2;
%             phi_ddot_R_avg      = (X_phi_ddot_avg*a_avg.phi_R).*(2*a_avg.f)^2;
%             
%             theta_ddot_LR_avg    = (X_theta_ddot_avg*a_avg.theta_LR).*(2*a_avg.f)^2;
%             eta_ddot_LR_avg      = (X_eta_ddot_avg*a_avg.eta_LR).*(2*a_avg.f)^2;
%             phi_ddot_LR_avg      = (X_phi_ddot_avg*a_avg.phi_LR).*(2*a_avg.f)^2;
%             
%             
%             d_theta_L_dt    = nan(1,1000);
%             d_eta_L_dt      = nan(1,1000);
%             d_phi_L_dt      = nan(1,1000);
% 
%             d_theta_R_dt    = nan(1,1000);
%             d_eta_R_dt      = nan(1,1000);
%             d_phi_R_dt      = nan(1,1000);
%             
%             d_theta_LR_dt    = nan(1,1000);
%             d_eta_LR_dt      = nan(1,1000);
%             d_phi_LR_dt      = nan(1,1000);
% 
%             d2_theta_L_dt2  = nan(1,1000);
%             d2_eta_L_dt2    = nan(1,1000);
%             d2_phi_L_dt2    = nan(1,1000);
% 
%             d2_theta_R_dt2  = nan(1,1000);
%             d2_eta_R_dt2    = nan(1,1000);
%             d2_phi_R_dt2    = nan(1,1000);
%             
%             d2_theta_LR_dt2  = nan(1,1000);
%             d2_eta_LR_dt2    = nan(1,1000);
%             d2_phi_LR_dt2    = nan(1,1000);
%             
%             dt2 = (1/a_avg.f)/999;
% 
%             for i = 1:1000
% 
%                 if i == 1
% 
%                     d_theta_L_dt(i)    = (theta_L_avg(i+1)-theta_L_avg(i))/dt2;
%                     d_eta_L_dt(i)      = (eta_L_avg(i+1)-eta_L_avg(i))/dt2;
%                     d_phi_L_dt(i)      = (phi_L_avg(i+1)-phi_L_avg(i))/dt2;
% 
%                     d_theta_R_dt(i)    = (theta_R_avg(i+1)-theta_R_avg(i))/dt2;
%                     d_eta_R_dt(i)      = (eta_R_avg(i+1)-eta_R_avg(i))/dt2;
%                     d_phi_R_dt(i)      = (phi_R_avg(i+1)-phi_R_avg(i))/dt2;
%                     
%                     d_theta_LR_dt(i)    = (theta_LR_avg(i+1)-theta_LR_avg(i))/dt2;
%                     d_eta_LR_dt(i)      = (eta_LR_avg(i+1)-eta_LR_avg(i))/dt2;
%                     d_phi_LR_dt(i)      = (phi_LR_avg(i+1)-phi_LR_avg(i))/dt2;
% 
%                     d2_theta_L_dt2(i)  = (theta_L_avg(i+2)-2*theta_L_avg(i+1)+theta_L_avg(i))/(dt2^2);
%                     d2_eta_L_dt2(i)    = (eta_L_avg(i+2)-2*eta_L_avg(i+1)+eta_L_avg(i))/(dt2^2);
%                     d2_phi_L_dt2(i)    = (phi_L_avg(i+2)-2*phi_L_avg(i+1)+phi_L_avg(i))/(dt2^2);
% 
%                     d2_theta_R_dt2(i)  = (theta_R_avg(i+2)-2*theta_R_avg(i+1)+theta_R_avg(i))/(dt2^2);
%                     d2_eta_R_dt2(i)    = (eta_R_avg(i+2)-2*eta_R_avg(i+1)+eta_R_avg(i))/(dt2^2);
%                     d2_phi_R_dt2(i)    = (phi_R_avg(i+2)-2*phi_R_avg(i+1)+phi_R_avg(i))/(dt2^2);
%                     
%                     d2_theta_LR_dt2(i)  = (theta_LR_avg(i+2)-2*theta_LR_avg(i+1)+theta_LR_avg(i))/(dt2^2);
%                     d2_eta_LR_dt2(i)    = (eta_LR_avg(i+2)-2*eta_LR_avg(i+1)+eta_LR_avg(i))/(dt2^2);
%                     d2_phi_LR_dt2(i)    = (phi_LR_avg(i+2)-2*phi_LR_avg(i+1)+phi_LR_avg(i))/(dt2^2);
% 
%                 elseif i == 1000
% 
%                     d_theta_L_dt(i)    = (theta_L_avg(i)-theta_L_avg(i-1))/dt2;
%                     d_eta_L_dt(i)      = (eta_L_avg(i)-eta_L_avg(i-1))/dt2;
%                     d_phi_L_dt(i)      = (phi_L_avg(i)-phi_L_avg(i-1))/dt2;
% 
%                     d_theta_R_dt(i)    = (theta_R_avg(i)-theta_R_avg(i-1))/dt2;
%                     d_eta_R_dt(i)      = (eta_R_avg(i)-eta_R_avg(i-1))/dt2;
%                     d_phi_R_dt(i)      = (phi_R_avg(i)-phi_R_avg(i-1))/dt2;
%                     
%                     d_theta_LR_dt(i)    = (theta_LR_avg(i)-theta_LR_avg(i-1))/dt2;
%                     d_eta_LR_dt(i)      = (eta_LR_avg(i)-eta_LR_avg(i-1))/dt2;
%                     d_phi_LR_dt(i)      = (phi_LR_avg(i)-phi_LR_avg(i-1))/dt2;
% 
%                     d2_theta_L_dt2(i)  = (theta_L_avg(i)-2*theta_L_avg(i-1)+theta_L_avg(i-2))/(dt2^2);
%                     d2_eta_L_dt2(i)    = (eta_L_avg(i)-2*eta_L_avg(i-1)+eta_L_avg(i-2))/(dt2^2);
%                     d2_phi_L_dt2(i)    = (phi_L_avg(i)-2*phi_L_avg(i-1)+phi_L_avg(i-2))/(dt2^2);
% 
%                     d2_theta_R_dt2(i)  = (theta_R_avg(i)-2*theta_R_avg(i-1)+theta_R_avg(i-2))/(dt2^2);
%                     d2_eta_R_dt2(i)    = (eta_R_avg(i)-2*eta_R_avg(i-1)+eta_R_avg(i-2))/(dt2^2);
%                     d2_phi_R_dt2(i)    = (phi_R_avg(i)-2*phi_R_avg(i-1)+phi_R_avg(i-2))/(dt2^2);
%                     
%                     d2_theta_LR_dt2(i)  = (theta_LR_avg(i)-2*theta_LR_avg(i-1)+theta_LR_avg(i-2))/(dt2^2);
%                     d2_eta_LR_dt2(i)    = (eta_LR_avg(i)-2*eta_LR_avg(i-1)+eta_LR_avg(i-2))/(dt2^2);
%                     d2_phi_LR_dt2(i)    = (phi_LR_avg(i)-2*phi_LR_avg(i-1)+phi_LR_avg(i-2))/(dt2^2);
% 
%                 else
% 
%                     d_theta_L_dt(i)    = (theta_L_avg(i+1)-theta_L_avg(i-1))/(2*dt2);
%                     d_eta_L_dt(i)      = (eta_L_avg(i+1)-eta_L_avg(i-1))/(2*dt2);
%                     d_phi_L_dt(i)      = (phi_L_avg(i+1)-phi_L_avg(i-1))/(2*dt2);
% 
%                     d_theta_R_dt(i)    = (theta_R_avg(i+1)-theta_R_avg(i-1))/(2*dt2);
%                     d_eta_R_dt(i)      = (eta_R_avg(i+1)-eta_R_avg(i-1))/(2*dt2);
%                     d_phi_R_dt(i)      = (phi_R_avg(i+1)-phi_R_avg(i-1))/(2*dt2);
%                     
%                     d_theta_LR_dt(i)    = (theta_LR_avg(i+1)-theta_LR_avg(i-1))/(2*dt2);
%                     d_eta_LR_dt(i)      = (eta_LR_avg(i+1)-eta_LR_avg(i-1))/(2*dt2);
%                     d_phi_LR_dt(i)      = (phi_LR_avg(i+1)-phi_LR_avg(i-1))/(2*dt2);
% 
%                     d2_theta_L_dt2(i)  = (theta_L_avg(i+1)-2*theta_L_avg(i)+theta_L_avg(i-1))/(dt2^2);
%                     d2_eta_L_dt2(i)    = (eta_L_avg(i+1)-2*eta_L_avg(i)+eta_L_avg(i-1))/(dt2^2);
%                     d2_phi_L_dt2(i)    = (phi_L_avg(i+1)-2*phi_L_avg(i)+phi_L_avg(i-1))/(dt2^2);
% 
%                     d2_theta_R_dt2(i)  = (theta_R_avg(i+1)-2*theta_R_avg(i)+theta_R_avg(i-1))/(dt2^2);
%                     d2_eta_R_dt2(i)    = (eta_R_avg(i+1)-2*eta_R_avg(i)+eta_R_avg(i-1))/(dt2^2);
%                     d2_phi_R_dt2(i)    = (phi_R_avg(i+1)-2*phi_R_avg(i)+phi_R_avg(i-1))/(dt2^2);
%                     
%                     d2_theta_LR_dt2(i)  = (theta_LR_avg(i+1)-2*theta_LR_avg(i)+theta_LR_avg(i-1))/(dt2^2);
%                     d2_eta_LR_dt2(i)    = (eta_LR_avg(i+1)-2*eta_LR_avg(i)+eta_LR_avg(i-1))/(dt2^2);
%                     d2_phi_LR_dt2(i)    = (phi_LR_avg(i+1)-2*phi_LR_avg(i)+phi_LR_avg(i-1))/(dt2^2);
% 
%                 end
% 
%             end    
%             
%         figure()
%         hold on
%         plot(t_avg,theta_L_avg,'r')
%         plot(t_avg,theta_R_avg,'g')
%         hold off
% 
%         figure()
%         hold on
%         plot(t_avg,eta_L_avg,'r')
%         plot(t_avg,eta_R_avg,'g')
%         hold off
% 
%         figure()
%         hold on
%         plot(t_avg,phi_L_avg,'r')
%         plot(t_avg,phi_R_avg,'g')
%         hold off
% 
%         figure()
%         hold on
%         plot(t_avg,theta_dot_L_avg,'b')
%         plot(t_avg,d_theta_L_dt,'r')
%         hold off
% 
%         figure()
%         hold on
%         plot(t_avg,eta_dot_L_avg,'b')
%         plot(t_avg,d_eta_L_dt,'r')
%         hold off
% 
%         figure()
%         hold on
%         plot(t_avg,phi_dot_L_avg,'b')
%         plot(t_avg,d_phi_L_dt,'r')
%         hold off
% 
%         figure()
%         hold on
%         plot(t_avg,theta_dot_R_avg,'b')
%         plot(t_avg,d_theta_R_dt,'g')
%         hold off
% 
%         figure()
%         hold on
%         plot(t_avg,eta_dot_R_avg,'b')
%         plot(t_avg,d_eta_R_dt,'g')
%         hold off
% 
%         figure()
%         hold on
%         plot(t_avg,phi_dot_R_avg,'b')
%         plot(t_avg,d_phi_R_dt,'g')
%         hold off
% 
%         figure()
%         hold on
%         plot(t_avg,theta_ddot_L_avg,'b')
%         plot(t_avg,d2_theta_L_dt2,'r')
%         hold off
% 
%         figure()
%         hold on
%         plot(t_avg,eta_ddot_L_avg,'b')
%         plot(t_avg,d2_eta_L_dt2,'r')
%         hold off
% 
%         figure()
%         hold on
%         plot(t_avg,phi_ddot_L_avg,'b')
%         plot(t_avg,d2_phi_L_dt2,'r')
%         hold off
% 
%         figure()
%         hold on
%         plot(t_avg,theta_ddot_R_avg,'b')
%         plot(t_avg,d2_theta_R_dt2,'g')
%         hold off
% 
%         figure()
%         hold on
%         plot(t_avg,eta_ddot_R_avg,'b')
%         plot(t_avg,d2_eta_R_dt2,'g')
%         hold off
% 
%         figure()
%         hold on
%         plot(t_avg,phi_ddot_R_avg,'b')
%         plot(t_avg,d2_phi_R_dt2,'g')
%         hold off
%         
%         figure()
%         plot(t_avg,theta_LR_avg,'r')
% 
%         figure()
%         plot(t_avg,eta_LR_avg,'r')
% 
%         figure()
%         plot(t_avg,phi_LR_avg,'r')
% 
%         figure()
%         hold on
%         plot(t_avg,theta_dot_LR_avg,'b')
%         plot(t_avg,d_theta_LR_dt,'r')
%         hold off
% 
%         figure()
%         hold on
%         plot(t_avg,eta_dot_LR_avg,'b')
%         plot(t_avg,d_eta_LR_dt,'r')
%         hold off
% 
%         figure()
%         hold on
%         plot(t_avg,phi_dot_LR_avg,'b')
%         plot(t_avg,d_phi_LR_dt,'r')
%         hold off
% 
%         figure()
%         hold on
%         plot(t_avg,theta_ddot_LR_avg,'b')
%         plot(t_avg,d2_theta_LR_dt2,'r')
%         hold off
% 
%         figure()
%         hold on
%         plot(t_avg,eta_ddot_LR_avg,'b')
%         plot(t_avg,d2_eta_LR_dt2,'r')
%         hold off
% 
%         figure()
%         hold on
%         plot(t_avg,phi_ddot_LR_avg,'b')
%         plot(t_avg,d2_phi_LR_dt2,'r')
%         hold off
% 
%         
%         pause
                
    end    


end

