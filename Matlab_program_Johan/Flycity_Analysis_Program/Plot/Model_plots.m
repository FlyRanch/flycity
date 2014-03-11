function Model_plots( settings, pathDB )


    % Test the Inertia model for conservation of linear and angular
    % momentum:
    
    seq_nr      = 1;
    nr_points   = 40;
    
    % Setup the body_model and wing_model:
    
    body_model.mass_fly         = pathDB.body_model.mass_fly(seq_nr);
    body_model.mass_body        = pathDB.body_model.mass_body(seq_nr);
    body_model.Joint_left       = pathDB.body_model.Joint_left(seq_nr,:)';
    body_model.Joint_right      = pathDB.body_model.Joint_right(seq_nr,:)';
    body_model.cg_b             = pathDB.body_model.cg(seq_nr,:)';
    body_model.Inertia          = pathDB.body_model.Inertia(:,:,seq_nr);
    body_model.g                = 1e-3*settings.g;
    body_model.x_mod            = pathDB.body_model.x_mod(:,:,seq_nr);
    body_model.y_mod            = pathDB.body_model.y_mod(:,:,seq_nr);
    body_model.z_mod            = pathDB.body_model.z_mod(:,:,seq_nr);
    
    wing_model.virtual_mass     = pathDB.wing_model.virtual_mass(seq_nr);
    wing_model.wing_cg_L        = pathDB.wing_model.wing_cg_L(seq_nr,:)';
    wing_model.wing_cg_R        = pathDB.wing_model.wing_cg_R(seq_nr,:)';
    wing_model.virtual_Inertia  = pathDB.wing_model.virtual_Inertia(:,:,seq_nr);
    wing_model.y_sect_L         = pathDB.wing_model.y_sect_L(:,:,seq_nr)';
    wing_model.chords_L         = pathDB.wing_model.chords_L(seq_nr,:)';
    wing_model.y_sect_R         = pathDB.wing_model.y_sect_R(:,:,seq_nr)';
    wing_model.chords_R         = pathDB.wing_model.chords_R(seq_nr,:)';
    wing_model.rho              = settings.rho_air;
    wing_model.length           = pathDB.wing_model.length(seq_nr);
    wing_model.x_mod_L          = pathDB.wing_model.x_mod_L(:,:,seq_nr);
    wing_model.y_mod_L          = pathDB.wing_model.y_mod_L(:,:,seq_nr);
    wing_model.z_mod_L          = pathDB.wing_model.z_mod_L(:,:,seq_nr);
    wing_model.x_mod_R          = pathDB.wing_model.x_mod_R(:,:,seq_nr);
    wing_model.y_mod_R          = pathDB.wing_model.y_mod_R(:,:,seq_nr);
    wing_model.z_mod_R          = pathDB.wing_model.z_mod_R(:,:,seq_nr);
    
    f                           = pathDB.poly_fit.a_glob.f;
    dt                          = (1/f)/(nr_points-1);
    R_strk                      = pathDB.rot_mat.Rstr;
    
    % Create the wing kinematics:
    
    n_pol_theta     = settings.n_pol_theta;
    n_pol_eta       = settings.n_pol_eta;
    n_pol_phi       = settings.n_pol_phi;

    a_avg_theta    = pathDB.poly_fit.a_glob.theta;
    a_avg_eta      = pathDB.poly_fit.a_glob.eta;
    a_avg_phi      = pathDB.poly_fit.a_glob.phi;
    down_up_avg    = pathDB.poly_fit.a_glob.down_up;

    a_theta_L_avg  = a_avg_theta;
    a_eta_L_avg    = a_avg_eta;
    a_phi_L_avg    = a_avg_phi;
    a_theta_R_avg  = a_avg_theta;
    a_eta_R_avg    = a_avg_eta;
    a_phi_R_avg    = a_avg_phi;
       
    % Create the average wing kinematics:

    a_fit_avg.a_theta_L     = a_theta_L_avg;
    a_fit_avg.a_eta_L       = a_eta_L_avg;
    a_fit_avg.a_phi_L       = a_phi_L_avg;
    a_fit_avg.a_theta_R     = a_theta_R_avg;
    a_fit_avg.a_eta_R       = a_eta_R_avg;
    a_fit_avg.a_phi_R       = a_phi_R_avg;
    a_fit_avg.f             = f;
    a_fit_avg.down_up       = down_up_avg;
    a_fit_avg.nr_points     = nr_points;
    a_fit_avg.R_strk        = R_strk;
        
    [ kine_avg ] = angular_velocities_polynomial( a_fit_avg );
    
    
    RL_avg = kine_avg.RL;
    RR_avg = kine_avg.RR;
    
    minx = -3.5;
    maxx = 3.5;
    miny = -3.5;
    maxy = 3.5;
    minz = -3.5;
    maxz = 3.5;
    
    figure(1)
    g1 = figure(1);
    hold on
    Fly_model_plot( [0; 0; 0], [1 0 0; 0 -1 0; 0 0 -1], eye(3), eye(3), body_model, wing_model )
    axis([minx maxx miny maxy minz maxz])
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    hold off
    set(gcf, 'InvertHardCopy', 'off');
    set(1,'Color','k')
%     saveas(g1,[char(settings.plot_loc) '/Manuevering_wing_kin2/' man_save_name ],'fig')




end

