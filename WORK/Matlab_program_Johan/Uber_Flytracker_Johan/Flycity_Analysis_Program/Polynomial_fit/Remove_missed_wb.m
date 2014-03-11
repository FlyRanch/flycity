function [ maneuver_star, weird_wb ] = Remove_missed_wb( settings, maneuver )

    wb_names = fieldnames( maneuver );
    
    nr_wb_man = length(wb_names);
    
    n_pol_theta = settings.n_pol_theta;
    n_pol_eta   = settings.n_pol_eta;
    n_pol_phi   = settings.n_pol_phi;
    
    a_dev_theta_L       = nan(n_pol_theta*2+2,nr_wb_man);
    a_dev_eta_L         = nan(n_pol_eta*2+2,nr_wb_man);
    a_dev_phi_L         = nan(n_pol_phi*2+2,nr_wb_man);
    
    a_dev_theta_R       = nan(n_pol_theta*2+2,nr_wb_man);
    a_dev_eta_R         = nan(n_pol_eta*2+2,nr_wb_man);
    a_dev_phi_R         = nan(n_pol_phi*2+2,nr_wb_man);
    
    for i = 1:nr_wb_man

        a_dev_theta_L(:,i)  = maneuver.(char(wb_names(i))).a_dev.theta_L;
        a_dev_eta_L(:,i)    = maneuver.(char(wb_names(i))).a_dev.eta_L;
        a_dev_phi_L(:,i)    = maneuver.(char(wb_names(i))).a_dev.phi_L;
        
        a_dev_theta_R(:,i)  = maneuver.(char(wb_names(i))).a_dev.theta_R;
        a_dev_eta_R(:,i)    = maneuver.(char(wb_names(i))).a_dev.eta_R;
        a_dev_phi_R(:,i)    = maneuver.(char(wb_names(i))).a_dev.phi_R;
        
    end
    
    % Remove missed wingbeats
    
    [ weird_wb ] = Remove_missed_wingbeats( a_dev_theta_L, a_dev_eta_L, a_dev_phi_L, a_dev_theta_R, a_dev_eta_R, a_dev_phi_R );
    
    maneuver_new = rmfield(maneuver,wb_names(weird_wb));
    
    % Rename the remaining wingbeats
    
    maneuver_new_names = fieldnames(maneuver_new);
    nr_wb_new = length(maneuver_new_names);
    
    maneuver_star = {};
    
    for i = 1:nr_wb_new
        
        maneuver_star.(char(['wb_' num2str(i)])) = maneuver_new.(char(maneuver_new_names(i)));
            
    end
    

end

