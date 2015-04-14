function [kine_avg, kine_man, FM_avg, FM_man] = Dynamic_simulation_plot( settings, pathDB, maneuver_type, wb_name, glob_on )


    % Simulate an average and a maneuvering wingbeat and return the forces,
    % moments, accelerations and velocities.
    
    % Step 1: Load the maneuvering and average wingbeat.
    
    seq_nr          = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).seq_nr;
    
    nr_points       = 1000;
    
    R_strk          = pathDB.rot_mat.Rstr;
    
    if glob_on == 0
    
        a_avg_theta_L   = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.theta_LR;
        a_avg_eta_L     = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.eta_LR;
        a_avg_phi_L     = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.phi_LR;
        a_avg_theta_R   = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.theta_LR;
        a_avg_eta_R     = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.eta_LR;
        a_avg_phi_R     = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.phi_LR;

        down_up         = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.down_up;
        f               = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.f;
    
    elseif glob_on == 1
        
        a_avg_theta_L   = pathDB.poly_fit.a_glob.theta;
        a_avg_eta_L     = pathDB.poly_fit.a_glob.eta;
        a_avg_phi_L     = pathDB.poly_fit.a_glob.phi;
        a_avg_theta_R   = pathDB.poly_fit.a_glob.theta;
        a_avg_eta_R     = pathDB.poly_fit.a_glob.eta;
        a_avg_phi_R     = pathDB.poly_fit.a_glob.phi;

        down_up         = pathDB.poly_fit.a_glob.down_up;
        f               = pathDB.poly_fit.a_glob.f;
        
    end
    
    a_man_theta_L   = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.theta_LR+pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_dev.theta_L;
    a_man_eta_L     = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.eta_LR+pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_dev.eta_L;
    a_man_phi_L     = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.phi_LR+pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_dev.phi_L;
    a_man_theta_R   = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.theta_LR+pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_dev.theta_R;
    a_man_eta_R     = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.eta_LR+pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_dev.eta_R;
    a_man_phi_R     = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.phi_LR+pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_dev.phi_R;
    

    
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
        
    a_fit.a_theta_L     = a_avg_theta_L;
    a_fit.a_eta_L       = a_avg_eta_L;
    a_fit.a_phi_L       = a_avg_phi_L;
    a_fit.a_theta_R     = a_avg_theta_R;
    a_fit.a_eta_R       = a_avg_eta_R;
    a_fit.a_phi_R       = a_avg_phi_R;
    a_fit.f             = f;
    a_fit.down_up       = down_up;
    a_fit.nr_points     = nr_points+1;
    a_fit.R_strk        = R_strk;
        
    [ kine_avg ] = angular_velocities_polynomial( a_fit );
    
    t_func              = kine_avg.t;
    kine_avg.R_strk     = R_strk;
    kine_avg.down_up_t  = t_func(end)*down_up;

    clear a_fit
    
    a_fit.a_theta_L     = a_man_theta_L;
    a_fit.a_eta_L       = a_man_eta_L;
    a_fit.a_phi_L       = a_man_phi_L;
    a_fit.a_theta_R     = a_man_theta_R;
    a_fit.a_eta_R       = a_man_eta_R;
    a_fit.a_phi_R       = a_man_phi_R;
    a_fit.f             = f;
    a_fit.down_up       = down_up;
    a_fit.nr_points     = nr_points+1;
    a_fit.R_strk        = R_strk;
        
    [ kine_man ] = angular_velocities_polynomial( a_fit );
    
    kine_man.R_strk     = R_strk;
    kine_man.down_up_t  = t_func(end)*down_up;
    
    clear a_fit
    
    % Average wingbeat
    
    % Body IC:
    
    body_kin.vb = [0;0;0];
    body_kin.wb = [0;0;0];
    
    wing_kin.RL     = kine_avg.RL(:,:,1);
    wing_kin.RR     = kine_avg.RR(:,:,1);
    wing_kin.wL     = kine_avg.wL(:,1);
    wing_kin.wR     = kine_avg.wR(:,1);
    wing_kin.wL_b   = kine_avg.wL_b(:,1);
    wing_kin.wR_b   = kine_avg.wR_b(:,1);
    
    [ M_matrix ] = Mass_matrix( body_model, wing_model, body_kin, wing_kin );
    
    body_IC.xyz_0   = [0; 0; 0];
    body_IC.qb_0    = quat2mat((R_strk'*[1 0 0; 0 -1 0; 0 0 -1]))';
    body_IC.vb_0    = M_matrix.vb_0;
    body_IC.wb_0    = M_matrix.wb_0;
    body_IC.F_0     = [0; 0; 0];
    body_IC.M_0     = [0; 0; 0];
    
    clear sim_data
    
    % Tolerance settings ode45:
    
    options = odeset('RelTol',1e-5,'AbsTol',[1e-15 1e-15 1e-15 1e-15 1e-15 1e-15]);
    
    % case_nr == 0: Inertia only
    % case_nr == 1: Inertia, translational aerodynamics, gravity
    % case_nr == 2: Inertia, translational and rotational aerodynamics,
    % gravity
    
    case_nr = 2;
    
    tic
    
    sim_data = dyn_sim(t_func,kine_avg,body_model,wing_model,body_IC,options,case_nr);
    
    toc
    
    n = 0;
    
    [n] = plot_dyn_sim( sim_data, case_nr, n+1 );
    
    clear body_IC case_nr options
    
    
    % Maneuvering wingbeat
    
    % Obtain start body kinematics:
    
%     xyz_strk
%     qb_strk
%     u_man_strk
    
    body_IC.xyz_0   = [0; 0; 0];
    qb_start = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).body_kin.qb(1,:)';
    body_IC.qb_0    = qb_start;
    Rb_start        = quat2mat(qb_start);
%     body_IC.qb_0    = quat2mat((R_strk'*[1 0 0; 0 -1 0; 0 0 -1]))';
    vb_start        = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).body_kin.uvw(1,:)';
    wb_start        = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).body_kin.wb(1,:)';
    body_IC.vb_0    = Rb_start*vb_start;
    body_IC.wb_0    = wb_start;
%     body_IC.vb_0    = M_matrix.vb_0;
%     body_IC.wb_0    = M_matrix.wb_0;
    body_IC.F_0     = [0; 0; 0];
    body_IC.M_0     = [0; 0; 0];
    
    clear sim_data
    
    % Tolerance settings ode45:
    
    options = odeset('RelTol',1e-5,'AbsTol',[1e-15 1e-15 1e-15 1e-15 1e-15 1e-15]);
    
    % case_nr == 0: Inertia only
    % case_nr == 1: Inertia, translational aerodynamics, gravity
    % case_nr == 2: Inertia, translational and rotational aerodynamics,
    % gravity
    
    case_nr = 2;
    
    tic
    
    sim_data = dyn_sim(t_func,kine_man,body_model,wing_model,body_IC,options,case_nr);
    
    toc
    
    [n] = plot_dyn_sim( sim_data, case_nr, n+1 );
    
    clear body_IC case_nr options
    
    sim_data
    
    % ook wingkinematics exporteren
    
end

