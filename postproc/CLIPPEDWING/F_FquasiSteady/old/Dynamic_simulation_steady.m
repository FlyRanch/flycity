function [kine_steady, sim_data] = Dynamic_simulation_steady( settings, pathDB, body_IC, wb_name, seq_nr, case_nr )

    tic

    % Simulate an average and a maneuvering wingbeat and return the forces,
    % moments, accelerations and velocities.
    
    % Step 1: Load the maneuvering and average wingbeat.
    
    nr_points       = 1000;
    
    R_strk          = pathDB.rot_mat.Rstr;
    
    a_steady_theta_L_t   = pathDB.rand_wbs.(char(wb_name)).a_avg.theta_LR+pathDB.rand_wbs.(char(wb_name)).a_dev.theta_L;
    a_steady_eta_L_t     = pathDB.rand_wbs.(char(wb_name)).a_avg.eta_LR+pathDB.rand_wbs.(char(wb_name)).a_dev.eta_L;
    a_steady_phi_L_t     = pathDB.rand_wbs.(char(wb_name)).a_avg.phi_LR+pathDB.rand_wbs.(char(wb_name)).a_dev.phi_L;
    a_steady_theta_R_t   = pathDB.rand_wbs.(char(wb_name)).a_avg.theta_LR+pathDB.rand_wbs.(char(wb_name)).a_dev.theta_R;
    a_steady_eta_R_t     = pathDB.rand_wbs.(char(wb_name)).a_avg.eta_LR+pathDB.rand_wbs.(char(wb_name)).a_dev.eta_R;
    a_steady_phi_R_t     = pathDB.rand_wbs.(char(wb_name)).a_avg.phi_LR+pathDB.rand_wbs.(char(wb_name)).a_dev.phi_R;
    
%     a_steady_theta_L   = pathDB.rand_wbs.(char(wb_name)).a_avg.theta_L+pathDB.rand_wbs.(char(wb_name)).a_dev.theta_L;
%     a_steady_eta_L     = pathDB.rand_wbs.(char(wb_name)).a_avg.eta_L+pathDB.rand_wbs.(char(wb_name)).a_dev.eta_L;
%     a_steady_phi_L     = pathDB.rand_wbs.(char(wb_name)).a_avg.phi_L+pathDB.rand_wbs.(char(wb_name)).a_dev.phi_L;
%     a_steady_theta_R   = pathDB.rand_wbs.(char(wb_name)).a_avg.theta_R+pathDB.rand_wbs.(char(wb_name)).a_dev.theta_R;
%     a_steady_eta_R     = pathDB.rand_wbs.(char(wb_name)).a_avg.eta_R+pathDB.rand_wbs.(char(wb_name)).a_dev.eta_R;
%     a_steady_phi_R     = pathDB.rand_wbs.(char(wb_name)).a_avg.phi_R+pathDB.rand_wbs.(char(wb_name)).a_dev.phi_R;

    down_up            = pathDB.rand_wbs.(char(wb_name)).a_fit.down_up;
    f                  = pathDB.rand_wbs.(char(wb_name)).a_fit.f;
    
    % Create periodic boundary conditions:
    
    a_steady_theta_L   = periodic_bc_func( a_steady_theta_L_t, down_up );
    a_steady_eta_L     = periodic_bc_func( a_steady_eta_L_t, down_up );
    a_steady_phi_L     = periodic_bc_func( a_steady_phi_L_t, down_up );
    a_steady_theta_R   = periodic_bc_func( a_steady_theta_R_t, down_up );
    a_steady_eta_R     = periodic_bc_func( a_steady_eta_R_t, down_up );
    a_steady_phi_R     = periodic_bc_func( a_steady_phi_R_t, down_up );
    
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
    
        
    % Step 3: Create the wing kinematics:
        
    a_fit.a_theta_L     = a_steady_theta_L;
    a_fit.a_eta_L       = a_steady_eta_L;
    a_fit.a_phi_L       = a_steady_phi_L;
    a_fit.a_theta_R     = a_steady_theta_R;
    a_fit.a_eta_R       = a_steady_eta_R;
    a_fit.a_phi_R       = a_steady_phi_R;
    a_fit.f             = f;
    a_fit.down_up       = down_up;
    a_fit.nr_points     = nr_points+1;
    a_fit.R_strk        = R_strk;
        
    [ kine_steady ] = angular_velocities_polynomial( a_fit );
    
    t_func                  = kine_steady.t;
    kine_steady.R_strk      = R_strk;
    kine_steady.down_up_t   = t_func(end)*down_up;

    clear a_fit
    

    % Step 4: Simulation
        
    % Tolerance settings ode45:
    
    options = odeset('RelTol',1e-5,'AbsTol',[1e-15 1e-15 1e-15 1e-15 1e-15 1e-15]);
    
    % case_nr == 0: Inertia only
    % case_nr == 1: Inertia, translational aerodynamics, gravity
    % case_nr == 2: Inertia, translational and rotational aerodynamics,
    % gravity
    
    % Simulation
    
    
    sim_data = dyn_sim(t_func,kine_steady,body_model,wing_model,body_IC,options,case_nr);

    toc
    
end

