function Polynomial_Fit( settings, pathDB )


    % A polynomial fit is made for both the downstroke and the upstroke in
    % a wingbeat. The polynomials which are used are Legendre polynomials.
    % Once a fit is made, an average wingbeat is determined. This average
    % wingbeat is used to determine the deviation from this wingbeat of
    % both left and right wings.
    
    savefile = 'pathDB5.mat';
    
    N = settings.nr_of_seq;
    
    M = settings.frame_end;
    
    n_pol_theta     = settings.n_pol_theta;
    n_pol_eta       = settings.n_pol_eta;
    n_pol_phi       = settings.n_pol_phi;
    
    w_tresh = settings.w_tresh;
    a_tresh = settings.a_tresh;
    
    
    
    % Obtain the direct, average and deviation polynomial fits:
    
    poly_fit.a_fit.theta_L         = nan(2*(n_pol_theta+1),200,N);
    poly_fit.a_fit.eta_L           = nan(2*(n_pol_eta+1),200,N);
    poly_fit.a_fit.phi_L           = nan(2*(n_pol_phi+1),200,N);
    poly_fit.a_fit.theta_R         = nan(2*(n_pol_theta+1),200,N);
    poly_fit.a_fit.eta_R           = nan(2*(n_pol_eta+1),200,N);
    poly_fit.a_fit.phi_R           = nan(2*(n_pol_phi+1),200,N);
    poly_fit.a_fit.down_up         = nan(200,N);
    poly_fit.a_fit.f               = nan(200,N);
    
    poly_fit.a_avg.theta_L         = nan(2*(n_pol_theta+1),N);
    poly_fit.a_avg.eta_L           = nan(2*(n_pol_eta+1),N);
    poly_fit.a_avg.phi_L           = nan(2*(n_pol_phi+1),N);
    poly_fit.a_avg.theta_R         = nan(2*(n_pol_theta+1),N);
    poly_fit.a_avg.eta_R           = nan(2*(n_pol_eta+1),N);
    poly_fit.a_avg.phi_R           = nan(2*(n_pol_phi+1),N);
    poly_fit.a_avg.theta_LR        = nan(2*(n_pol_theta+1),N);
    poly_fit.a_avg.eta_LR          = nan(2*(n_pol_eta+1),N);
    poly_fit.a_avg.phi_LR          = nan(2*(n_pol_phi+1),N);
    poly_fit.a_avg.down_up         = nan(N,1);
    poly_fit.a_avg.f               = nan(N,1);
    poly_fit.a_avg.trigger_wb      = nan(N,1);
    poly_fit.a_avg.m_fly           = nan(N,1);
    poly_fit.a_avg.wing_l          = nan(N,1);
        
    poly_fit.a_dev.theta_L         = nan(2*(n_pol_theta+1),200,N);
    poly_fit.a_dev.eta_L           = nan(2*(n_pol_eta+1),200,N);
    poly_fit.a_dev.phi_L           = nan(2*(n_pol_phi+1),200,N);
    poly_fit.a_dev.theta_R         = nan(2*(n_pol_theta+1),200,N);
    poly_fit.a_dev.eta_R           = nan(2*(n_pol_eta+1),200,N);
    poly_fit.a_dev.phi_R           = nan(2*(n_pol_phi+1),200,N);
    
    count = 1;
    
    for i = 1:N
        
        nr_wb = pathDB.wingbeats.nr_of_wb(i);
    
        [ a_fit_t, a_avg_t, a_dev_t ] = Standard_Wingbeat(settings,pathDB,i);
        
        % Save per sequence:
        
        poly_fit.a_fit.theta_L(:,:,i)       = a_fit_t.theta_L;
        poly_fit.a_fit.eta_L(:,:,i)         = a_fit_t.eta_L;
        poly_fit.a_fit.phi_L(:,:,i)         = a_fit_t.phi_L;
        poly_fit.a_fit.theta_R(:,:,i)       = a_fit_t.theta_R;
        poly_fit.a_fit.eta_R(:,:,i)         = a_fit_t.eta_R;
        poly_fit.a_fit.phi_R(:,:,i)         = a_fit_t.phi_R;
        poly_fit.a_fit.down_up(1:nr_wb,i)   = a_fit_t.down_up;
        poly_fit.a_fit.f(1:nr_wb,i)         = a_fit_t.f;
        
        poly_fit.a_avg.theta_L(:,i)         = a_avg_t.theta_L;
        poly_fit.a_avg.eta_L(:,i)           = a_avg_t.eta_L;
        poly_fit.a_avg.phi_L(:,i)           = a_avg_t.phi_L;
        poly_fit.a_avg.theta_R(:,i)         = a_avg_t.theta_R;
        poly_fit.a_avg.eta_R(:,i)           = a_avg_t.eta_R;
        poly_fit.a_avg.phi_R(:,i)           = a_avg_t.phi_R;
        poly_fit.a_avg.theta_LR(:,i)        = a_avg_t.theta_LR;
        poly_fit.a_avg.eta_LR(:,i)          = a_avg_t.eta_LR;
        poly_fit.a_avg.phi_LR(:,i)          = a_avg_t.phi_LR;
        poly_fit.a_avg.down_up(i)           = a_avg_t.down_up;
        poly_fit.a_avg.f(i)                 = a_avg_t.f;
        poly_fit.a_avg.trigger_wb(i)        = a_avg_t.trigger_wb;
        poly_fit.a_avg.m_fly(i)             = a_avg_t.m_fly;
        poly_fit.a_avg.wing_l(i)            = a_avg_t.wing_l;
        
        poly_fit.a_dev.theta_L(:,:,i)       = a_dev_t.theta_L;
        poly_fit.a_dev.eta_L(:,:,i)         = a_dev_t.eta_L;
        poly_fit.a_dev.phi_L(:,:,i)         = a_dev_t.phi_L;
        poly_fit.a_dev.theta_R(:,:,i)       = a_dev_t.theta_R;
        poly_fit.a_dev.eta_R(:,:,i)         = a_dev_t.eta_R;
        poly_fit.a_dev.phi_R(:,:,i)         = a_dev_t.phi_R;
        
        % Save per wingbeat:
        
        if isnan(a_avg_t.f) == 0
        
        for j = 1:nr_wb
        
            wb_name = ['wb_' num2str(count)];
                        
            rand_wbs.(char(wb_name)).seq_name                       = char(settings.sequence_names(i,:));
            rand_wbs.(char(wb_name)).wb_nr                          = j;
            
            rand_wbs.(char(wb_name)).body_model.mass_fly            = pathDB.body_model.mass_fly(i);
            rand_wbs.(char(wb_name)).body_model.mass_body           = pathDB.body_model.mass_body(i);
            rand_wbs.(char(wb_name)).body_model.Inertia             = pathDB.body_model.Inertia(:,:,i);
            rand_wbs.(char(wb_name)).body_model.Joint_left          = pathDB.body_model.Joint_left(i,:)';
            rand_wbs.(char(wb_name)).body_model.Joint_right         = pathDB.body_model.Joint_right(i,:)';
            rand_wbs.(char(wb_name)).body_model.x_mod               = pathDB.body_model.x_mod(:,:,i);
            rand_wbs.(char(wb_name)).body_model.y_mod               = pathDB.body_model.y_mod(:,:,i);
            rand_wbs.(char(wb_name)).body_model.z_mod               = pathDB.body_model.z_mod(:,:,i);
            rand_wbs.(char(wb_name)).body_model.cg_b                = pathDB.body_model.cg(i,:)';
            
            rand_wbs.(char(wb_name)).wing_model.length              = pathDB.wing_model.length(i);
            rand_wbs.(char(wb_name)).wing_model.mass                = pathDB.wing_model.mass(i);
            rand_wbs.(char(wb_name)).wing_model.virtual_mass        = pathDB.wing_model.virtual_mass(i);
            rand_wbs.(char(wb_name)).wing_model.Inertia             = pathDB.wing_model.Inertia(:,:,i);
            rand_wbs.(char(wb_name)).wing_model.virtual_Inertia     = pathDB.wing_model.virtual_Inertia(:,:,i);
            rand_wbs.(char(wb_name)).wing_model.area                = pathDB.wing_model.area(i);
            rand_wbs.(char(wb_name)).wing_model.wing_cg_L           = pathDB.wing_model.wing_cg_L(i,:);
            rand_wbs.(char(wb_name)).wing_model.y_sect_L            = pathDB.wing_model.y_sect_L(:,:,i);
            rand_wbs.(char(wb_name)).wing_model.chords_L            = pathDB.wing_model.chords_L(i,:);
            rand_wbs.(char(wb_name)).wing_model.x_mod_L             = pathDB.wing_model.x_mod_L(:,:,i);
            rand_wbs.(char(wb_name)).wing_model.y_mod_L             = pathDB.wing_model.y_mod_L(:,:,i);
            rand_wbs.(char(wb_name)).wing_model.z_mod_L             = pathDB.wing_model.z_mod_L(:,:,i);
            rand_wbs.(char(wb_name)).wing_model.wing_cg_R           = pathDB.wing_model.wing_cg_R(i,:);
            rand_wbs.(char(wb_name)).wing_model.y_sect_R            = pathDB.wing_model.y_sect_R(:,:,i);
            rand_wbs.(char(wb_name)).wing_model.chords_R            = pathDB.wing_model.chords_R(i,:);
            rand_wbs.(char(wb_name)).wing_model.x_mod_R             = pathDB.wing_model.x_mod_R(:,:,i);
            rand_wbs.(char(wb_name)).wing_model.y_mod_R             = pathDB.wing_model.y_mod_R(:,:,i);
            rand_wbs.(char(wb_name)).wing_model.z_mod_R             = pathDB.wing_model.z_mod_R(:,:,i);
            rand_wbs.(char(wb_name)).wing_model.rho                 = settings.rho_air;
            
            wb_loc = pathDB.wingbeats.wingbeat_loc(j,:,i);
            
            rand_wbs.(char(wb_name)).wingbeats.nr_of_wb             = nr_wb;
            rand_wbs.(char(wb_name)).wingbeats.wingbeat_loc         = pathDB.wingbeats.wingbeat_loc(j,:,i);
            rand_wbs.(char(wb_name)).wingbeats.downstroke_loc       = pathDB.wingbeats.downstroke_loc(j,:,i);
            rand_wbs.(char(wb_name)).wingbeats.upstroke_loc         = pathDB.wingbeats.upstroke_loc(j,:,i);
            
            rand_wbs.(char(wb_name)).strkpln_kin.alfa               = pathDB.strkpln_kin.alfa(i,wb_loc(1):wb_loc(2));
            rand_wbs.(char(wb_name)).strkpln_kin.beta               = pathDB.strkpln_kin.beta(i,wb_loc(1):wb_loc(2));
            rand_wbs.(char(wb_name)).strkpln_kin.roll               = pathDB.strkpln_kin.roll(i,wb_loc(1):wb_loc(2));
            rand_wbs.(char(wb_name)).strkpln_kin.pitch              = pathDB.strkpln_kin.pitch(i,wb_loc(1):wb_loc(2));
            rand_wbs.(char(wb_name)).strkpln_kin.yaw                = pathDB.strkpln_kin.yaw(i,wb_loc(1):wb_loc(2));
            rand_wbs.(char(wb_name)).strkpln_kin.uvw                = pathDB.strkpln_kin.uvw(wb_loc(1):wb_loc(2),:,i);
            rand_wbs.(char(wb_name)).strkpln_kin.w                  = pathDB.strkpln_kin.w(wb_loc(1):wb_loc(2),:,i);
            rand_wbs.(char(wb_name)).strkpln_kin.a_xyz              = pathDB.strkpln_kin.a_xyz(wb_loc(1):wb_loc(2),:,i);
            rand_wbs.(char(wb_name)).strkpln_kin.Fg                 = pathDB.strkpln_kin.Fg(wb_loc(1):wb_loc(2),:,i);
            rand_wbs.(char(wb_name)).strkpln_kin.alfa_mean          = mean(pathDB.strkpln_kin.alfa(i,wb_loc(1):wb_loc(2)));
            rand_wbs.(char(wb_name)).strkpln_kin.beta_mean          = mean(pathDB.strkpln_kin.beta(i,wb_loc(1):wb_loc(2)));
            rand_wbs.(char(wb_name)).strkpln_kin.roll_mean          = mean(pathDB.strkpln_kin.roll(i,wb_loc(1):wb_loc(2)));
            rand_wbs.(char(wb_name)).strkpln_kin.pitch_mean         = mean(pathDB.strkpln_kin.pitch(i,wb_loc(1):wb_loc(2)));
            rand_wbs.(char(wb_name)).strkpln_kin.yaw_mean           = mean(pathDB.strkpln_kin.yaw(i,wb_loc(1):wb_loc(2)));
            rand_wbs.(char(wb_name)).strkpln_kin.uvw_mean           = mean(pathDB.strkpln_kin.uvw(wb_loc(1):wb_loc(2),:,i));
            rand_wbs.(char(wb_name)).strkpln_kin.w_mean             = mean(pathDB.strkpln_kin.w(wb_loc(1):wb_loc(2),:,i));
            rand_wbs.(char(wb_name)).strkpln_kin.a_xyz_mean         = mean(pathDB.strkpln_kin.a_xyz(wb_loc(1):wb_loc(2),:,i));
            rand_wbs.(char(wb_name)).strkpln_kin.Fg_mean            = mean(pathDB.strkpln_kin.Fg(wb_loc(1):wb_loc(2),:,i));
            
            rand_wbs.(char(wb_name)).wing_kin.theta_L               = pathDB.wing_kin.theta_L(i,wb_loc(1):wb_loc(2));
            rand_wbs.(char(wb_name)).wing_kin.eta_L                 = pathDB.wing_kin.eta_L(i,wb_loc(1):wb_loc(2));
            rand_wbs.(char(wb_name)).wing_kin.phi_L                 = pathDB.wing_kin.phi_L(i,wb_loc(1):wb_loc(2));
            rand_wbs.(char(wb_name)).wing_kin.theta_R               = pathDB.wing_kin.theta_R(i,wb_loc(1):wb_loc(2));
            rand_wbs.(char(wb_name)).wing_kin.eta_R                 = pathDB.wing_kin.eta_R(i,wb_loc(1):wb_loc(2));
            rand_wbs.(char(wb_name)).wing_kin.phi_R                 = pathDB.wing_kin.phi_R(i,wb_loc(1):wb_loc(2));
            
            
            rand_wbs.(char(wb_name)).a_fit.theta_L                  = a_fit_t.theta_L(:,j);
            rand_wbs.(char(wb_name)).a_fit.eta_L                    = a_fit_t.eta_L(:,j);
            rand_wbs.(char(wb_name)).a_fit.phi_L                    = a_fit_t.phi_L(:,j);
            rand_wbs.(char(wb_name)).a_fit.theta_R                  = a_fit_t.theta_R(:,j);
            rand_wbs.(char(wb_name)).a_fit.eta_R                    = a_fit_t.eta_R(:,j);
            rand_wbs.(char(wb_name)).a_fit.phi_R                    = a_fit_t.phi_R(:,j);
            rand_wbs.(char(wb_name)).a_fit.down_up                  = a_fit_t.down_up(j);
            rand_wbs.(char(wb_name)).a_fit.f                        = a_fit_t.f(j);
            
            rand_wbs.(char(wb_name)).a_avg.theta_L                  = a_avg_t.theta_L;
            rand_wbs.(char(wb_name)).a_avg.eta_L                    = a_avg_t.eta_L;
            rand_wbs.(char(wb_name)).a_avg.phi_L                    = a_avg_t.phi_L;
            rand_wbs.(char(wb_name)).a_avg.theta_R                  = a_avg_t.theta_R;
            rand_wbs.(char(wb_name)).a_avg.eta_R                    = a_avg_t.eta_R;
            rand_wbs.(char(wb_name)).a_avg.phi_R                    = a_avg_t.phi_R;
            rand_wbs.(char(wb_name)).a_avg.theta_LR                 = a_avg_t.theta_LR;
            rand_wbs.(char(wb_name)).a_avg.eta_LR                   = a_avg_t.eta_LR;
            rand_wbs.(char(wb_name)).a_avg.phi_LR                   = a_avg_t.phi_LR;
            rand_wbs.(char(wb_name)).a_avg.down_up                  = a_avg_t.down_up;
            rand_wbs.(char(wb_name)).a_avg.f                        = a_avg_t.f;
            rand_wbs.(char(wb_name)).a_avg.trigger_wb               = a_avg_t.trigger_wb;
            
            rand_wbs.(char(wb_name)).a_dev.theta_L                  = a_dev_t.theta_L(:,j);
            rand_wbs.(char(wb_name)).a_dev.eta_L                    = a_dev_t.eta_L(:,j);
            rand_wbs.(char(wb_name)).a_dev.phi_L                    = a_dev_t.phi_L(:,j);
            rand_wbs.(char(wb_name)).a_dev.theta_R                  = a_dev_t.theta_R(:,j);
            rand_wbs.(char(wb_name)).a_dev.eta_R                    = a_dev_t.eta_R(:,j);
            rand_wbs.(char(wb_name)).a_dev.phi_R                    = a_dev_t.phi_R(:,j);
        
        
            count = count+1;
            
        end
        
        end
    
    end
    
    % Determine the global wingbeat:
    
    [ a_glob ] = Global_wingbeat( poly_fit.a_avg );
    
    poly_fit.a_glob = a_glob;
    
    
    % Sort maneuvering wingbeats based on 
    
    count_ax = 1;
    count_ay = 1;
    count_az = 1;
    count_wx = 1;
    count_wy = 1;
    count_wz = 1;
    
    wb_names = fieldnames(rand_wbs);
    
    for k = 1:length(wb_names)
            
                
        a_xyz_mean = rand_wbs.(char(wb_names(k))).strkpln_kin.a_xyz_mean;
        w_xyz_mean = rand_wbs.(char(wb_names(k))).strkpln_kin.w_mean;
            
        if abs(a_xyz_mean(1)) >= a_tresh
                
            wb_name_ax = ['wb_' num2str(count_ax)];
                
            maneuver.ax.(char(wb_name_ax)) = rand_wbs.(char(wb_names(k)));
                
            count_ax = count_ax+1;
                
        end
        
        if abs(a_xyz_mean(2)) >= a_tresh
                
            wb_name_ay = ['wb_' num2str(count_ay)];
                
            maneuver.ay.(char(wb_name_ay)) = rand_wbs.(char(wb_names(k)));
                
            count_ay = count_ay+1;
                
        end
        
        if abs(a_xyz_mean(3)) >= a_tresh
                
            wb_name_az = ['wb_' num2str(count_az)];
                
            maneuver.az.(char(wb_name_az)) = rand_wbs.(char(wb_names(k)));
                
            count_az = count_az+1;
                
        end
        
        if abs(w_xyz_mean(1)) >= w_tresh
                
            wb_name_wx = ['wb_' num2str(count_wx)];
                
            maneuver.wx.(char(wb_name_wx)) = rand_wbs.(char(wb_names(k)));
                
            count_wx = count_wx+1;
                
        end
        
        if abs(w_xyz_mean(2)) >= w_tresh
                
            wb_name_wy = ['wb_' num2str(count_wy)];
                
            maneuver.wy.(char(wb_name_wy)) = rand_wbs.(char(wb_names(k)));
                
            count_wy = count_wy+1;
                
        end
        
        if abs(w_xyz_mean(3)) >= w_tresh
                
            wb_name_wz = ['wb_' num2str(count_wz)];
                
            maneuver.wz.(char(wb_name_wz)) = rand_wbs.(char(wb_names(k)));
                
            count_wz = count_wz+1;
                
        end
            
    end
    
    % Save the structures:
    
    save(savefile,'poly_fit','rand_wbs','maneuver')

end

