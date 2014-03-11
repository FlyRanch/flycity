function [ c_fit ] = asymmetric_motion_regression( maneuver, down_up_glob )

    % Relate deviation coefficients to the non-dimensional forces and
    % non-dimensional moments:
    
    n_pol_theta     = (size(maneuver.a_dev_theta_L,1)-2)/2;
    n_pol_eta       = (size(maneuver.a_dev_eta_L,1)-2)/2;
    n_pol_phi       = (size(maneuver.a_dev_phi_L,1)-2)/2;
    
    nr_wb           = size(maneuver.a_dev_theta_L,2);
    
    F_star          = maneuver.F_star;
    M_star          = maneuver.M_star;
    
    a_dev_theta_L      = maneuver.a_dev_theta_L;
    a_dev_eta_L        = maneuver.a_dev_eta_L;
    a_dev_phi_L        = maneuver.a_dev_phi_L;
    
    a_dev_theta_R      = maneuver.a_dev_theta_R;
    a_dev_eta_R        = maneuver.a_dev_eta_R;
    a_dev_phi_R        = maneuver.a_dev_phi_R;
     

    
    % Split positive and negative values of F_sort and M_sort and a_dev
    % accordingly:
    
    b_theta_L     = nan(n_pol_theta*2+2,3,6);
    b_eta_L       = nan(n_pol_eta*2+2,3,6);
    b_phi_L       = nan(n_pol_phi*2+2,3,6);
    
    b_theta_R     = nan(n_pol_theta*2+2,3,6);
    b_eta_R       = nan(n_pol_eta*2+2,3,6);
    b_phi_R       = nan(n_pol_phi*2+2,3,6);  
    
    b_theta_L_0     = nan(n_pol_theta*2+2,3,6);
    b_eta_L_0       = nan(n_pol_eta*2+2,3,6);
    b_phi_L_0       = nan(n_pol_phi*2+2,3,6);
    
    b_theta_R_0     = nan(n_pol_theta*2+2,3,6);
    b_eta_R_0       = nan(n_pol_eta*2+2,3,6);
    b_phi_R_0       = nan(n_pol_phi*2+2,3,6);  
    
    for k = 1:6
        
        if k <= 3
            
            [b,FM_LR_t,N_t] = dev_fit_asymmetric(F_star(k,:),a_dev_theta_L,a_dev_eta_L,a_dev_phi_L,a_dev_theta_R,a_dev_eta_R,a_dev_phi_R,down_up_glob);

            b_theta_L(:,:,k)   = b.theta_L;
            b_eta_L(:,:,k)     = b.eta_L;
            b_phi_L(:,:,k)     = b.phi_L;

            b_theta_R(:,:,k)   = b.theta_R;
            b_eta_R(:,:,k)     = b.eta_R;
            b_phi_R(:,:,k)     = b.phi_R;
            
            b_theta_L_0(:,:,k)   = b.theta_L_0;
            b_eta_L_0(:,:,k)     = b.eta_L_0;
            b_phi_L_0(:,:,k)     = b.phi_L_0;

            b_theta_R_0(:,:,k)   = b.theta_R_0;
            b_eta_R_0(:,:,k)     = b.eta_R_0;
            b_phi_R_0(:,:,k)     = b.phi_R_0;

            FM.(['FM_' char(num2str(k))])   = FM_LR_t;
            N.(['FM_' char(num2str(k))])    = N_t;
            
        else

            [b,FM_LR_t,N_t] = dev_fit_asymmetric(M_star(k-3,:),a_dev_theta_L,a_dev_eta_L,a_dev_phi_L,a_dev_theta_R,a_dev_eta_R,a_dev_phi_R,down_up_glob);

            b_theta_L(:,:,k)   = b.theta_L;
            b_eta_L(:,:,k)     = b.eta_L;
            b_phi_L(:,:,k)     = b.phi_L;

            b_theta_R(:,:,k)   = b.theta_R;
            b_eta_R(:,:,k)     = b.eta_R;
            b_phi_R(:,:,k)     = b.phi_R;
            
            b_theta_L_0(:,:,k)   = b.theta_L_0;
            b_eta_L_0(:,:,k)     = b.eta_L_0;
            b_phi_L_0(:,:,k)     = b.phi_L_0;

            b_theta_R_0(:,:,k)   = b.theta_R_0;
            b_eta_R_0(:,:,k)     = b.eta_R_0;
            b_phi_R_0(:,:,k)     = b.phi_R_0;
            
            FM.(['FM_' char(num2str(k))])   = FM_LR_t;
            N.(['FM_' char(num2str(k))])    = N_t;

            
        end
        
    end
    
    c_fit.b_theta_L_0   = b_theta_L_0;
    c_fit.b_eta_L_0     = b_eta_L_0;
    c_fit.b_phi_L_0     = b_phi_L_0;

    c_fit.b_theta_R_0   = b_theta_R_0;
    c_fit.b_eta_R_0     = b_eta_R_0;
    c_fit.b_phi_R_0     = b_phi_R_0;
    
    c_fit.b_theta_L   = b_theta_L;
    c_fit.b_eta_L     = b_eta_L;
    c_fit.b_phi_L     = b_phi_L;

    c_fit.b_theta_R   = b_theta_R;
    c_fit.b_eta_R     = b_eta_R;
    c_fit.b_phi_R     = b_phi_R;
    
    c_fit.N           = N;
    c_fit.FM          = FM;


end

