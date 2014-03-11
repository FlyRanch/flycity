function [ man ] = FM_aero_maneuver(maneuver,R_strk,wb_names,nr_wb,n_pol_theta,n_pol_eta,n_pol_phi)


    man.a_dev_theta_L   = nan(n_pol_theta*2+2,nr_wb);
    man.a_dev_eta_L     = nan(n_pol_eta*2+2,nr_wb);
    man.a_dev_phi_L     = nan(n_pol_phi*2+2,nr_wb);
    man.a_dev_theta_R   = nan(n_pol_theta*2+2,nr_wb);
    man.a_dev_eta_R     = nan(n_pol_eta*2+2,nr_wb);
    man.a_dev_phi_R     = nan(n_pol_phi*2+2,nr_wb);
    man.FM_aero_man     = nan(6,nr_wb);
    man.FM_aero_avg     = nan(6,nr_wb);
    man.F_star          = nan(3,nr_wb);
    man.M_star          = nan(3,nr_wb);
    man.f               = nan(1,nr_wb);
    man.down_up         = nan(1,nr_wb);
    man.wing_L          = nan(1,nr_wb);
    man.mass_fly        = nan(1,nr_wb);
    
    nr_points           = 201;
    
    nr_wb
    
    for i = 1:nr_wb
        
        i
        
        % Obtain a_avg and a_dev
        
        a_avg_theta     = maneuver.(char(wb_names(i))).a_avg.theta_LR;
        a_avg_eta       = maneuver.(char(wb_names(i))).a_avg.eta_LR;
        a_avg_phi       = maneuver.(char(wb_names(i))).a_avg.phi_LR;
        
        a_dev_theta_L   = maneuver.(char(wb_names(i))).a_dev.theta_L;
        a_dev_eta_L     = maneuver.(char(wb_names(i))).a_dev.eta_L;
        a_dev_phi_L     = maneuver.(char(wb_names(i))).a_dev.phi_L;
        a_dev_theta_R   = maneuver.(char(wb_names(i))).a_dev.theta_R;
        a_dev_eta_R     = maneuver.(char(wb_names(i))).a_dev.eta_R;
        a_dev_phi_R     = maneuver.(char(wb_names(i))).a_dev.phi_R;
        
        f               = maneuver.(char(wb_names(i))).a_avg.f;
        down_up         = maneuver.(char(wb_names(i))).a_avg.down_up;
        
        % Obtain body model:
        
        body_model      = maneuver.(char(wb_names(i))).body_model;
                        
        % Obtain wing model:
        
        wing_model          = maneuver.(char(wb_names(i))).wing_model;
        wing_model.y_sect_L = wing_model.y_sect_L';
        wing_model.y_sect_R = wing_model.y_sect_R';
        wing_model.chords_L = wing_model.chords_L';
        wing_model.chords_R = wing_model.chords_R';
        
        m_fly           = body_model.mass_fly;
        wing_L          = wing_model.length;
        
        % Set up the wing kinematics for the maneuver:
        
        a_theta_L       = a_avg_theta+a_dev_theta_L;
        a_eta_L         = a_avg_eta+a_dev_eta_L;
        a_phi_L         = a_avg_phi+a_dev_phi_L;
        a_theta_R       = a_avg_theta+a_dev_theta_R;
        a_eta_R         = a_avg_eta+a_dev_eta_R;
        a_phi_R         = a_avg_phi+a_dev_phi_R;

        a_fit.a_theta_L = a_theta_L;
        a_fit.a_eta_L   = a_eta_L;
        a_fit.a_phi_L   = a_phi_L;
        a_fit.a_theta_R = a_theta_R;
        a_fit.a_eta_R   = a_eta_R;
        a_fit.a_phi_R   = a_phi_R;
        a_fit.f         = f;
        a_fit.down_up   = down_up;
        a_fit.nr_points = nr_points;
        a_fit.R_strk    = R_strk;
        
        [ kine_t ] = angular_velocities( a_fit );
        
        kine.u_strk          = zeros(3,nr_points);
        kine.w_strk          = zeros(3,nr_points);
        kine.wL              = kine_t.wL(:,1:(nr_points));
        kine.wR              = kine_t.wR(:,1:(nr_points));
        kine.RL              = kine_t.RL(:,:,1:(nr_points));
        kine.RR              = kine_t.RR(:,:,1:(nr_points));
        kine.R_strk          = R_strk;
        
        clear kine_t
        
        wb.wb_loc            = [ 1 nr_points ];
        wb.down_loc          = [ 1 round(down_up*nr_points) ];
        wb.up_loc            = [ round(down_up*nr_points)+1 nr_points ];
        wb.dt                = (1/f)/(nr_points-1);
        
        % Compute aerodynamic forces for the maneuver:
        
        [ FM_strkpln_man, ~, ~, ~, ~, ~, ~, ~, ~ ] = Aerodynamic_forces( kine, body_model, wing_model, wb, 1);
        
        clear kine
        
        % Set up the wing kinematics for the average wingbeat:
        
        a_fit.a_theta_L = a_avg_theta;
        a_fit.a_eta_L   = a_avg_eta;
        a_fit.a_phi_L   = a_avg_phi;
        a_fit.a_theta_R = a_avg_theta;
        a_fit.a_eta_R   = a_avg_eta;
        a_fit.a_phi_R   = a_avg_phi;
        
        [ kine_t ] = angular_velocities( a_fit );
        
        kine.u_strk          = zeros(3,nr_points);
        kine.w_strk          = zeros(3,nr_points);
        kine.wL              = kine_t.wL(:,1:(nr_points));
        kine.wR              = kine_t.wR(:,1:(nr_points));
        kine.RL              = kine_t.RL(:,:,1:(nr_points));
        kine.RR              = kine_t.RR(:,:,1:(nr_points));
        kine.R_strk          = R_strk;
        
        clear kine_t a_fit
        
        % Compute aerodynamic forces for the average wingbeat:
        
        [ FM_strkpln_avg, ~, ~, ~, ~, ~, ~, ~, ~ ] = Aerodynamic_forces( kine, body_model, wing_model, wb, 1);
        
        clear kine
        
        % Non-dimensionalize
        
        F_star = (mean(FM_strkpln_man(1:3,:),2)-mean(FM_strkpln_avg(1:3,:),2))/(m_fly*wing_L*f^2);
        
        M_star = (mean(FM_strkpln_man(4:6,:),2)-mean(FM_strkpln_avg(4:6,:),2))/(m_fly*wing_L^2*f^2);
        
        % Save the data
        
        man.a_dev_theta_L(:,i)  = a_dev_theta_L;
        man.a_dev_eta_L(:,i)    = a_dev_eta_L;
        man.a_dev_phi_L(:,i)    = a_dev_phi_L;
        man.a_dev_theta_R(:,i)  = a_dev_theta_R;
        man.a_dev_eta_R(:,i)    = a_dev_eta_R;
        man.a_dev_phi_R(:,i)    = a_dev_phi_R;
        man.FM_aero_man(:,i)    = mean(FM_strkpln_man,2);
        man.FM_aero_avg(:,i)    = mean(FM_strkpln_avg,2);
        man.F_star(:,i)         = F_star;
        man.M_star(:,i)         = M_star;
        man.f(i)                = f;
        man.down_up(i)          = down_up;
        man.wing_L(i)           = wing_L;
        man.mass_fly(i)         = m_fly;
        
        clear FM_strkpln_man FM_strkpln_avg F_star M_star body_model wing_model wb

    end
    
end

