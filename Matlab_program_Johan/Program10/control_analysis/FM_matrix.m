function [ FM_matrix_theta, FM_matrix_eta, FM_matrix_phi ] = FM_matrix(maneuver, maneuver_names, nr_maneuver, n_pol_theta, n_pol_eta, n_pol_phi, force_on, moment_on )


    % Compute the deviation coefficients versus the symmetric forces or
    % moments:

    
    
    FM_matrix_theta = [];
    FM_matrix_eta = [];
    FM_matrix_phi = [];
    
    
    
    for i = 1:nr_maneuver
        
        temp_man = maneuver.(char(maneuver_names(i)));
        
        f_avg = temp_man.f_avg;
        down_up_avg = temp_man.down_up_avg;
        
        f = temp_man.f;
        down_up = temp_man.down_up;
        
        wing_l = temp_man.wing_l;
        mass_fly = temp_man.mass_fly;
        
        a_avg = temp_man.a_avg;
        
        a_dev = temp_man.a_dev;
        
        a_theta_L = [a_dev.theta_L1; a_dev.theta_L2];
        a_eta_L = [a_dev.eta_L1; a_dev.eta_L2];
        a_phi_L = [a_dev.phi_L1; a_dev.phi_L2];
        
        a_theta_R = [a_dev.theta_R1; a_dev.theta_R2];
        a_eta_R = [a_dev.eta_R1; a_dev.eta_R2];
        a_phi_R = [a_dev.phi_R1; a_dev.phi_R2];
        
        a_theta_LR = [a_avg.theta_LR1; a_avg.theta_LR2];
        a_eta_LR = [a_avg.eta_LR1; a_avg.eta_LR2];
        a_phi_LR = [a_avg.phi_LR1; a_avg.phi_LR2];
       
        
        if force_on == 1
        
            F_I_avg     = temp_man.FM_req.F_M_avg.F_I_mean;
            F_a_L_avg   = temp_man.FM_req.F_M_avg.F_a_L_mean;
            F_a_R_avg   = temp_man.FM_req.F_M_avg.F_a_R_mean;
            F_g_avg     = temp_man.FM_req.F_M_avg.F_g_mean;
            
            F_req_avg = F_I_avg+F_a_L_avg+F_a_R_avg+F_g_avg;
            
            wb_names = fieldnames(temp_man.FM_req);
            
            nr_wb = length(f);
            
            FM_temp_theta = zeros(nr_wb,12+(n_pol_theta+1)*6);
            FM_temp_eta = zeros(nr_wb,12+(n_pol_eta+1)*6);
            FM_temp_phi = zeros(nr_wb,12+(n_pol_phi+1)*6);
            
            for j = 1:nr_wb
               
                F_I     = temp_man.FM_req.(char(wb_names(j))).F_I_mean;
                F_a_L   = temp_man.FM_req.(char(wb_names(j))).F_a_L_mean;
                F_a_R   = temp_man.FM_req.(char(wb_names(j))).F_a_R_mean;
                F_g     = temp_man.FM_req.(char(wb_names(j))).F_g_mean;
                
                F_req = F_I+F_a_L+F_a_R+F_g;
                                
                FM_temp_theta(j,:)   = [ F_req' wing_l f(j) down_up(j) F_req_avg' f_avg down_up_avg mass_fly a_theta_L(:,j)' a_theta_R(:,j)' a_theta_LR' ];
                FM_temp_eta(j,:)     = [ F_req' wing_l f(j) down_up(j) F_req_avg' f_avg down_up_avg mass_fly a_eta_L(:,j)' a_eta_R(:,j)' a_eta_LR' ];
                FM_temp_phi(j,:)     = [ F_req' wing_l f(j) down_up(j) F_req_avg' f_avg down_up_avg mass_fly a_phi_L(:,j)' a_phi_R(:,j)' a_phi_LR' ];
                
            end
            
        FM_matrix_theta = [ FM_matrix_theta ; FM_temp_theta ];
        FM_matrix_eta = [ FM_matrix_eta ; FM_temp_eta ];
        FM_matrix_phi = [ FM_matrix_phi ; FM_temp_phi ];    
            
        end
        
        
        if moment_on == 1
            
                    
            M_I_avg     = temp_man.FM_req.F_M_avg.M_I_mean;
            M_a_L_avg   = temp_man.FM_req.F_M_avg.M_a_L_mean;
            M_a_R_avg   = temp_man.FM_req.F_M_avg.M_a_R_mean;
            
            M_req_avg = M_I_avg + M_a_L_avg + M_a_R_avg;
            
            wb_names = fieldnames(temp_man.FM_req);
            
            nr_wb = length(f);
            
            FM_temp_theta = zeros(nr_wb,12+(n_pol_theta+1)*6);
            FM_temp_eta = zeros(nr_wb,12+(n_pol_eta+1)*6);
            FM_temp_phi = zeros(nr_wb,12+(n_pol_phi+1)*6);
            
            for j = 1:nr_wb
               
                
                M_I     = temp_man.FM_req.(char(wb_names(j))).M_I_mean;
                M_a_L   = temp_man.FM_req.(char(wb_names(j))).M_a_L_mean;
                M_a_R   = temp_man.FM_req.(char(wb_names(j))).M_a_R_mean;
                
                M_req = M_I+M_a_L+M_a_R;
                                
                FM_temp_theta(j,:)   = [ M_req' wing_l f(j) down_up(j) M_req_avg' f_avg down_up_avg mass_fly a_theta_L(:,j)' a_theta_R(:,j)' a_theta_LR' ];
                FM_temp_eta(j,:)     = [ M_req' wing_l f(j) down_up(j) M_req_avg' f_avg down_up_avg mass_fly a_eta_L(:,j)' a_eta_R(:,j)' a_eta_LR' ];
                FM_temp_phi(j,:)     = [ M_req' wing_l f(j) down_up(j) M_req_avg' f_avg down_up_avg mass_fly a_phi_L(:,j)' a_phi_R(:,j)' a_phi_LR' ];
                
            end
            
        FM_matrix_theta = [ FM_matrix_theta ; FM_temp_theta ];
        FM_matrix_eta = [ FM_matrix_eta ; FM_temp_eta ];
        FM_matrix_phi = [ FM_matrix_phi ; FM_temp_phi ];    
            
        end
        
        
    end
    


end

