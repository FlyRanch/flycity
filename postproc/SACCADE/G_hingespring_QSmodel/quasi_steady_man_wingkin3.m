function [ FM_strkpln, FM_inertia, Vel_wingtip ] = quasi_steady_man_wingkin3( settings, pathDB, theta_L, eta_L, phi_L, theta_R, eta_R, phi_R, f, down_up, seq_nr, rot_on )

    % Load a_glob and maneuvering wing kinematic functions:
        
    f               = 7*f;
    
    R_strk          = pathDB.rot_mat.Rstr;
    
    nr_of_points    = length(theta_L);
    
            
%     % Seq_nr 98 is cooresponding the best to the global fly conditions:
%     
%     seq_nr = 92;
    
    % Step 2: Setup the body_model and wing_model:
    
    body_model.mass_fly         = pathDB.body_model.mass_fly(seq_nr);
    body_model.mass_body        = pathDB.body_model.mass_body(seq_nr);
    body_model.Joint_left       = pathDB.body_model.Joint_left(seq_nr,:)';
    body_model.Joint_right      = pathDB.body_model.Joint_right(seq_nr,:)';
    body_model.cg_b             = pathDB.body_model.cg(seq_nr,:)';
    body_model.Inertia          = pathDB.body_model.Inertia(:,:,seq_nr);
    body_model.x_mod            = pathDB.body_model.x_mod(:,:,seq_nr);
    body_model.y_mod            = pathDB.body_model.y_mod(:,:,seq_nr);
    body_model.z_mod            = pathDB.body_model.z_mod(:,:,seq_nr);
    body_model.g                = 1e-3*settings.g;
    
    wing_model.virtual_mass     = pathDB.wing_model.virtual_mass(seq_nr);
    wing_model.mass             = pathDB.wing_model.mass(seq_nr);
    wing_model.wing_cg_L        = pathDB.wing_model.wing_cg_L(seq_nr,:)';
    wing_model.wing_cg_R        = pathDB.wing_model.wing_cg_R(seq_nr,:)';
    wing_model.virtual_Inertia  = pathDB.wing_model.virtual_Inertia(:,:,seq_nr);
    wing_model.y_sect_L         = pathDB.wing_model.y_sect_L(:,:,seq_nr)';
    wing_model.chords_L         = pathDB.wing_model.chords_L(seq_nr,:)';
    wing_model.x_LE_L           = pathDB.wing_model.x_LE_L(seq_nr,:)';
    wing_model.y_sect_R         = pathDB.wing_model.y_sect_R(:,:,seq_nr)';
    wing_model.chords_R         = pathDB.wing_model.chords_R(seq_nr,:)';
    wing_model.x_LE_R           = pathDB.wing_model.x_LE_R(seq_nr,:)';
    wing_model.x_mod_L          = pathDB.wing_model.x_mod_L(:,:,seq_nr);
    wing_model.y_mod_L          = pathDB.wing_model.y_mod_L(:,:,seq_nr);
    wing_model.z_mod_L          = pathDB.wing_model.z_mod_L(:,:,seq_nr);
    wing_model.x_mod_R          = pathDB.wing_model.x_mod_R(:,:,seq_nr);
    wing_model.y_mod_R          = pathDB.wing_model.y_mod_R(:,:,seq_nr);
    wing_model.z_mod_R          = pathDB.wing_model.z_mod_R(:,:,seq_nr);
    wing_model.rho              = settings.rho_air;
    wing_model.nr_sect          = size(wing_model.y_sect_L,2);
    
    
    % Compute the global quasi-steady forces and moments:
   
    [ wingkin ] = angular_velocities_polynomial2( theta_L, eta_L, phi_L, theta_R, eta_R, phi_R, f, R_strk );

    kine.R_strk          = R_strk;
    kine.u_strk          = zeros(3,nr_of_points);
    kine.w_strk          = zeros(3,nr_of_points);
    kine.wL              = wingkin.wL;
    kine.wR              = wingkin.wR;
    kine.w_dot_L         = wingkin.w_dot_L;
    kine.w_dot_R         = wingkin.w_dot_R;
    kine.RL              = wingkin.RL;
    kine.RR              = wingkin.RR;

    kine.vb              = zeros(3,nr_of_points);
    kine.wb              = zeros(3,nr_of_points);
    kine.wL_b            = wingkin.wL_b;
    kine.wR_b            = wingkin.wR_b;
    kine.w_dot_L_b       = wingkin.w_dot_L_b;
    kine.w_dot_R_b       = wingkin.w_dot_R_b;
    kine.theta_dot_L     = wingkin.theta_dot_L;
    kine.eta_dot_L       = wingkin.eta_dot_L;
    kine.phi_dot_L       = wingkin.phi_dot_L;
    kine.theta_ddot_L    = wingkin.theta_ddot_L;
    kine.eta_ddot_L      = wingkin.eta_ddot_L;
    kine.phi_ddot_L      = wingkin.phi_ddot_L;
    
    FM_strkpln = {};
    
    [ FM_strkpln.FM, ~, ~ ,U_left, U_right, alfa_L, alfa_R, ~, ~ ] = Aerodynamic_forces( kine, body_model, wing_model, rot_on);
     
    [ ~, FM_inertia ] = Inertia_instantaneous( kine, body_model, wing_model );
    
    Vel_wingtip.U_left    = U_left;
    Vel_wingtip.U_right   = U_right;
    Vel_wingtip.alfa_L    = alfa_L;
    Vel_wingtip.alfa_R    = alfa_R;

    FM_strkpln.kine = kine;
    
end