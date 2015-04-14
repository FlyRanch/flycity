function [kine_avg, sim_data] = Dynamic_simulation_avg( settings, pathDB, seq_nr, glob_on, case_nr )

    tic

    % Simulate an average and a maneuvering wingbeat and return the forces,
    % moments, accelerations and velocities.
    
    % Step 1: Load the maneuvering and average wingbeat.
    
    nr_points       = 1000;
    
    R_strk          = pathDB.rot_mat.Rstr;
    
    if glob_on == 0
    
        a_avg_theta_L   = pathDB.poly_fit.a_avg.theta_LR(:,seq_nr);
        a_avg_eta_L     = pathDB.poly_fit.a_avg.eta_LR(:,seq_nr);
        a_avg_phi_L     = pathDB.poly_fit.a_avg.phi_LR(:,seq_nr);
        a_avg_theta_R   = pathDB.poly_fit.a_avg.theta_LR(:,seq_nr);
        a_avg_eta_R     = pathDB.poly_fit.a_avg.eta_LR(:,seq_nr);
        a_avg_phi_R     = pathDB.poly_fit.a_avg.phi_LR(:,seq_nr);

        down_up         = pathDB.poly_fit.a_avg.down_up(seq_nr);
        f               = pathDB.poly_fit.a_avg.f(seq_nr);
    
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

%     if glob_on == 0
%     
%         a_avg_theta_L   = pathDB.poly_fit.a_avg.theta_L(:,seq_nr);
%         a_avg_eta_L     = pathDB.poly_fit.a_avg.eta_L(:,seq_nr);
%         a_avg_phi_L     = pathDB.poly_fit.a_avg.phi_L(:,seq_nr);
%         a_avg_theta_R   = pathDB.poly_fit.a_avg.theta_R(:,seq_nr);
%         a_avg_eta_R     = pathDB.poly_fit.a_avg.eta_R(:,seq_nr);
%         a_avg_phi_R     = pathDB.poly_fit.a_avg.phi_R(:,seq_nr);
% 
%         down_up         = pathDB.poly_fit.a_avg.down_up(seq_nr);
%         f               = pathDB.poly_fit.a_avg.f(seq_nr);
%     
%     elseif glob_on == 1
%         
%         a_avg_theta_L   = pathDB.poly_fit.a_glob.theta;
%         a_avg_eta_L     = pathDB.poly_fit.a_glob.eta;
%         a_avg_phi_L     = pathDB.poly_fit.a_glob.phi;
%         a_avg_theta_R   = pathDB.poly_fit.a_glob.theta;
%         a_avg_eta_R     = pathDB.poly_fit.a_glob.eta;
%         a_avg_phi_R     = pathDB.poly_fit.a_glob.phi;
% 
%         down_up         = pathDB.poly_fit.a_glob.down_up;
%         f               = pathDB.poly_fit.a_glob.f;
%         
%     end
    
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
    

    % Step 4: vb_0 and wb_0
    
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
        
    % Simulation
        
    sim_data1 = dyn_sim(t_func,kine_avg,body_model,wing_model,body_IC,options,case_nr);

    F_tot1  = sim_data1.FI_acc_b_mean;
    M_tot1  = sim_data1.MI_acc_b_mean;
    
    body_IC.xyz_0   = [0; 0; 0];
    body_IC.qb_0    = quat2mat((R_strk'*[1 0 0; 0 -1 0; 0 0 -1]))';
    body_IC.vb_0    = M_matrix.vb_0;
    body_IC.wb_0    = M_matrix.wb_0;
    body_IC.F_0     = F_tot1;
    body_IC.M_0     = M_tot1;

    sim_data2 = dyn_sim(t_func,kine_avg,body_model,wing_model,body_IC,options,case_nr);
    
    F_tot2  = sim_data2.FI_acc_b_mean;
    M_tot2  = sim_data2.MI_acc_b_mean;
    
    F_0 = F_tot1+F_tot2;
    M_0 = M_tot1+M_tot2;
    
    sim_data = sim_data2;
    sim_data.F_0 = F_0;
    sim_data.M_0 = M_0;

    toc
    
end

