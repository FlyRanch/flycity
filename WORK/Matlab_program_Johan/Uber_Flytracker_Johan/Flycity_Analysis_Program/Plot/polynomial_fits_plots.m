function polynomial_fits_plots( settings, pathDB )

    % Create figures for the chapter on polynomial fits:
    
    seq_nr = 5;
    
    % Legendre polynomials upon order 5:
    
    x = -1:(2/99):1;
    
    Px_0 = ones(100,1);
    Px_1 = x;
    Px_2 = 0.5*(3*x.^2-1);
    Px_3 = 0.5*(5*x.^3-3*x);
    Px_4 = 0.125*(35*x.^4-30*x.^2+3);
    Px_5 = 0.125*(63*x.^5-70*x.^3+15*x);
    
    figure(1)
    hFig = figure(1);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 1200 1000]);
    hold on
    plot(x,Px_0,'k')
    plot(x,Px_1,'b')
    plot(x,Px_2,'r')
    plot(x,Px_3,'g')
    plot(x,Px_4,'m')
    plot(x,Px_5,'c')
    hold off
    title('Legendre polynomials')
    xlabel('x')
    ylabel('y')
    ylim([-1.1 1.1])
    legend('P_0','P_1','P_2','P_3','P_4','P_5')
    saveas(hFig,[char(settings.plot_loc) '/Polynomial_fits/Legendre_polynomials'],'fig')
    
    
    % Raw wing kinematics and polynomial fit:
    
    nr_wb = pathDB.wingbeats.nr_of_wb(seq_nr);
    
    wb_loc = pathDB.wingbeats.wingbeat_loc(1:nr_wb,:,seq_nr);
    
    start = settings.start_stop(seq_nr,1);
    stop = settings.start_stop(seq_nr,2);
    
    t = pathDB.t;
    
    theta_L_raw     = pathDB.wing_kin.theta_L(seq_nr,:);
    eta_L_raw       = pathDB.wing_kin.eta_L(seq_nr,:);
    phi_L_raw       = pathDB.wing_kin.phi_L(seq_nr,:);
    
    n_pol_theta     = settings.n_pol_theta;
    n_pol_eta       = settings.n_pol_eta;
    n_pol_phi       = settings.n_pol_phi;
    
    a_theta_L       = pathDB.poly_fit.a_fit.theta_L(:,1:nr_wb,seq_nr);
    a_eta_L         = pathDB.poly_fit.a_fit.eta_L(:,1:nr_wb,seq_nr);
    a_phi_L         = pathDB.poly_fit.a_fit.phi_L(:,1:nr_wb,seq_nr);
    
    down_up         = pathDB.poly_fit.a_fit.down_up(1:nr_wb,seq_nr);
    f               = pathDB.poly_fit.a_fit.f(1:nr_wb,seq_nr);
    
    
    theta_L_fit     = nan(1,5588);
    eta_L_fit       = nan(1,5588);
    phi_L_fit       = nan(1,5588);
    
    for i = 1:nr_wb
        
        data_points = wb_loc(i,2)-wb_loc(i,1)+2;
        
        [~, X_theta]     = Wingbeat_Legendre_matrix( n_pol_theta, down_up(i), data_points, 0, 1, 0 );
        [~, X_eta]       = Wingbeat_Legendre_matrix( n_pol_eta, down_up(i), data_points, 0, 1, 0 );
        [~, X_phi]       = Wingbeat_Legendre_matrix( n_pol_phi, down_up(i), data_points, 0, 1, 0 );
        
        theta_L_fit(wb_loc(i,1):(wb_loc(i,2)+1))    = X_theta*a_theta_L(:,i);
        eta_L_fit(wb_loc(i,1):(wb_loc(i,2)+1))      = X_eta*a_eta_L(:,i);
        phi_L_fit(wb_loc(i,1):(wb_loc(i,2)+1))      = X_phi*a_phi_L(:,i);
        
    end
    
    figure(2)
    hFig = figure(2);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 1400 800]);
    hold on
    subplot(3,3,[1 2]); hold on
    plot(t(start:stop),theta_L_raw(start:stop),'b')
    plot(t(start:stop),theta_L_fit(start:stop),'r')
    ylabel('\theta_L [rad]')
    xlim([t(start) t(stop)])
    hold off
    subplot(3,3,3); hold on
    plot(t(wb_loc(20,1):wb_loc(22,2)),theta_L_raw(wb_loc(20,1):wb_loc(22,2)),'b')
    plot(t(wb_loc(20,1):wb_loc(22,2)),theta_L_fit(wb_loc(20,1):wb_loc(22,2)),'r')
    xlim([t(wb_loc(20,1)) t(wb_loc(22,2))])
    hold off
    subplot(3,3,[4 5]); hold on
    plot(t(start:stop),eta_L_raw(start:stop),'b')
    plot(t(start:stop),eta_L_fit(start:stop),'r')
    ylabel('\eta_L [rad]')
    xlim([t(start) t(stop)])
    hold off
    hold off
    subplot(3,3,6); hold on
    plot(t(wb_loc(20,1):wb_loc(22,2)),eta_L_raw(wb_loc(20,1):wb_loc(22,2)),'b')
    plot(t(wb_loc(20,1):wb_loc(22,2)),eta_L_fit(wb_loc(20,1):wb_loc(22,2)),'r')
    xlim([t(wb_loc(20,1)) t(wb_loc(22,2))])
    hold off
    subplot(3,3,[7 8]); hold on
    plot(t(start:stop),phi_L_raw(start:stop),'b')
    plot(t(start:stop),phi_L_fit(start:stop),'r')
    xlabel('t [s]')
    ylabel('\phi_L [rad]')
    xlim([t(start) t(stop)])
    hold off
    subplot(3,3,9); hold on
    plot(t(wb_loc(20,1):wb_loc(22,2)),phi_L_raw(wb_loc(20,1):wb_loc(22,2)),'b')
    plot(t(wb_loc(20,1):wb_loc(22,2)),phi_L_fit(wb_loc(20,1):wb_loc(22,2)),'r')
    xlabel('t [s]')
    xlim([t(wb_loc(20,1)) t(wb_loc(22,2))])
    hold off
    hold off
    [~,h1] = suplabel('Polynomial fit left wing, seq 5', 't');
    set(h1,'FontSize',10)
    saveas(hFig,[char(settings.plot_loc) '/Polynomial_fits/Raw_poly_fit'],'fig')
    
    
    
    % Average left, average right and average symmetric wingbeat:
    
    a_avg_theta_L   = pathDB.poly_fit.a_avg.theta_L(:,seq_nr);
    a_avg_eta_L     = pathDB.poly_fit.a_avg.eta_L(:,seq_nr);
    a_avg_phi_L     = pathDB.poly_fit.a_avg.phi_L(:,seq_nr);
    
    a_avg_theta_R   = pathDB.poly_fit.a_avg.theta_R(:,seq_nr);
    a_avg_eta_R     = pathDB.poly_fit.a_avg.eta_R(:,seq_nr);
    a_avg_phi_R     = pathDB.poly_fit.a_avg.phi_R(:,seq_nr);
    
    a_avg_theta_LR  = pathDB.poly_fit.a_avg.theta_LR(:,seq_nr);
    a_avg_eta_LR    = pathDB.poly_fit.a_avg.eta_LR(:,seq_nr);
    a_avg_phi_LR    = pathDB.poly_fit.a_avg.phi_LR(:,seq_nr);
    
    down_up_avg     = pathDB.poly_fit.a_avg.down_up(seq_nr);
    f_avg           = pathDB.poly_fit.a_avg.f(seq_nr);    
    
    [t_avg, X_theta_avg]    = Wingbeat_Legendre_matrix( n_pol_theta, down_up_avg, 101, 0, 1/f_avg, 0 );
    [~, X_eta_avg]          = Wingbeat_Legendre_matrix( n_pol_eta, down_up_avg, 101, 0, 1/f_avg, 0 );
    [~, X_phi_avg]          = Wingbeat_Legendre_matrix( n_pol_phi, down_up_avg, 101, 0, 1/f_avg, 0 );
    
    theta_L_avg     = X_theta_avg*a_avg_theta_L;
    eta_L_avg       = X_eta_avg*a_avg_eta_L;
    phi_L_avg       = X_phi_avg*a_avg_phi_L;
    
    theta_R_avg     = X_theta_avg*a_avg_theta_R;
    eta_R_avg       = X_eta_avg*a_avg_eta_R;
    phi_R_avg       = X_phi_avg*a_avg_phi_R;
    
    theta_LR_avg    = X_theta_avg*a_avg_theta_LR;
    eta_LR_avg      = X_eta_avg*a_avg_eta_LR;
    phi_LR_avg      = X_phi_avg*a_avg_phi_LR;
    
    figure(3)
    hFig = figure(3);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 1400 800]);
    hold on
    subplot(3,1,1); plot(t_avg,theta_L_avg,'r',t_avg,theta_R_avg,'b',t_avg,theta_LR_avg,'m')
    ylabel('\theta [rad]')
    subplot(3,1,2); plot(t_avg,eta_L_avg,'r',t_avg,eta_R_avg,'b',t_avg,eta_LR_avg,'m')
    ylabel('\eta [rad]')
    subplot(3,1,3); plot(t_avg,phi_L_avg,'r',t_avg,phi_R_avg,'b',t_avg,phi_LR_avg,'m')
    xlabel('t [s]')
    ylabel('\phi [rad]')
    legend('left','right','symmetric')
    hold off
    [~,h1] = suplabel('Average wingbeats, seq 5', 't');
    set(h1,'FontSize',10)
    saveas(hFig,[char(settings.plot_loc) '/Polynomial_fits/avg_poly_fit'],'fig')
    
    
end

