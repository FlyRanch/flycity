function [ c_fit ] = symmetric_motion_regression( maneuver, down_up_glob)

    % Relate deviation coefficients to the non-dimensional forces and
    % non-dimensional moments:
    
    n_pol_theta         = (size(maneuver.a_dev_theta_L,1)-2)/2;
    n_pol_eta           = (size(maneuver.a_dev_eta_L,1)-2)/2;
    n_pol_phi           = (size(maneuver.a_dev_phi_L,1)-2)/2;
    
    nr_wb               = size(maneuver.a_dev_theta_L,2);
    
    F_star              = maneuver.F_star;
    M_star              = maneuver.M_star;
    
    a_dev_theta_L       = maneuver.a_dev_theta_L;
    a_dev_eta_L         = maneuver.a_dev_eta_L;
    a_dev_phi_L         = maneuver.a_dev_phi_L;
    
    a_dev_theta_R       = maneuver.a_dev_theta_R;
    a_dev_eta_R         = maneuver.a_dev_eta_R;
    a_dev_phi_R         = maneuver.a_dev_phi_R;
    
    a_dev_theta         = nan(n_pol_theta*2+2,nr_wb);
    a_dev_eta           = nan(n_pol_eta*2+2,nr_wb);
    a_dev_phi           = nan(n_pol_phi*2+2,nr_wb);
    
    for i = 1:nr_wb

        a_dev_theta(:,i)    = mean([a_dev_theta_L(:,i) a_dev_theta_R(:,i)],2);
        a_dev_eta(:,i)      = mean([a_dev_eta_L(:,i) a_dev_eta_R(:,i)],2);
        a_dev_phi(:,i)      = mean([a_dev_phi_L(:,i) a_dev_phi_R(:,i)],2);
        
    end
    
    
    % Split positive and negative values of F_sort and M_sort and a_dev
    % accordingly:
    
    b_theta_n    = zeros(n_pol_theta*2+2,3,6);
    b_eta_n      = zeros(n_pol_eta*2+2,3,6);
    b_phi_n      = zeros(n_pol_phi*2+2,3,6);
    b_theta_p    = zeros(n_pol_theta*2+2,3,6);
    b_eta_p      = zeros(n_pol_eta*2+2,3,6);
    b_phi_p      = zeros(n_pol_phi*2+2,3,6);

    
    for k = 1:6
        
        if k <= 3
            
            [ b, FM_n_t, FM_p_t, N_n_t, N_p_t, N_t ] = dev_fit_symmetric( F_star(k,:), a_dev_theta, a_dev_eta, a_dev_phi, down_up_glob );

            b_theta_n(:,:,k)   = b.theta_n;
            b_eta_n(:,:,k)     = b.eta_n;
            b_phi_n(:,:,k)     = b.phi_n;

            b_theta_p(:,:,k)   = b.theta_p;
            b_eta_p(:,:,k)     = b.eta_p;
            b_phi_p(:,:,k)     = b.phi_p;
            
            FM.(['FM_' char(num2str(k))])   = [FM_n_t FM_p_t];
            FM_n.(['FM_' char(num2str(k))]) = FM_n_t;
            FM_p.(['FM_' char(num2str(k))]) = FM_p_t;
            N_n.(['FM_' char(num2str(k))])  = N_n_t;
            N_p.(['FM_' char(num2str(k))])  = N_p_t;
            N.(['FM_' char(num2str(k))])    = N_t;

        else
            
            [ b, FM_n_t, FM_p_t, N_n_t, N_p_t, N_t ] = dev_fit_symmetric( M_star(k-3,:), a_dev_theta, a_dev_eta, a_dev_phi, down_up_glob );

            b_theta_n(:,:,k)   = b.theta_n;
            b_eta_n(:,:,k)     = b.eta_n;
            b_phi_n(:,:,k)     = b.phi_n;

            b_theta_p(:,:,k)   = b.theta_p;
            b_eta_p(:,:,k)     = b.eta_p;
            b_phi_p(:,:,k)     = b.phi_p;


            FM.(['FM_' char(num2str(k))])   = [FM_n_t FM_p_t];
            FM_n.(['FM_' char(num2str(k))]) = FM_n_t;
            FM_p.(['FM_' char(num2str(k))]) = FM_p_t;
            N_n.(['FM_' char(num2str(k))])  = N_n_t;
            N_p.(['FM_' char(num2str(k))])  = N_p_t;
            N.(['FM_' char(num2str(k))])    = N_t;

       
        end   

    end

    c_fit.b_theta_n  = b_theta_n;
    c_fit.b_eta_n    = b_eta_n;
    c_fit.b_phi_n    = b_phi_n;
    c_fit.b_theta_p  = b_theta_p;
    c_fit.b_eta_p    = b_eta_p;
    c_fit.b_phi_p    = b_phi_p;
    
    c_fit.FM         = FM;
    c_fit.FM_n       = FM_n;
    c_fit.FM_p       = FM_p;
    c_fit.N_n        = N_n;
    c_fit.N_p        = N_p;
    c_fit.N          = N;
    

end

