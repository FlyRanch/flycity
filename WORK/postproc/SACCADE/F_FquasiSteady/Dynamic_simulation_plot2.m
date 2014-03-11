function Dynamic_simulation_plot2( settings, pathDB, maneuver_type, wb_name, sym_on )

    % Step 1: obtain wing kinematics:

    seq_nr          = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).seq_nr;
    
    wb_nr           = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).wb_nr;
    
    nr_points       = 1000;
    
    R_strk          = pathDB.rot_mat.Rstr;
    
    if sym_on == 1
        
        a_avg_theta_L   = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.theta_LR;
        a_avg_eta_L     = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.eta_LR;
        a_avg_phi_L     = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.phi_LR;
        a_avg_theta_R   = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.theta_LR;
        a_avg_eta_R     = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.eta_LR;
        a_avg_phi_R     = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.phi_LR;
    
        a_man_theta_L_t = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.theta_LR+pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_dev.theta_L;
        a_man_eta_L_t   = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.eta_LR+pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_dev.eta_L;
        a_man_phi_L_t   = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.phi_LR+pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_dev.phi_L;
        a_man_theta_R_t = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.theta_LR+pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_dev.theta_R;
        a_man_eta_R_t   = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.eta_LR+pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_dev.eta_R;
        a_man_phi_R_t   = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.phi_LR+pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_dev.phi_R;
        
        down_up         = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_fit.down_up;
        f               = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_fit.f;
        
        % Create periodic boundary conditions
    
        a_man_theta_L   = periodic_bc_func( a_man_theta_L_t, down_up );
        a_man_eta_L     = periodic_bc_func( a_man_eta_L_t, down_up );
        a_man_phi_L     = periodic_bc_func( a_man_phi_L_t, down_up );
        a_man_theta_R   = periodic_bc_func( a_man_theta_R_t, down_up );
        a_man_eta_R     = periodic_bc_func( a_man_eta_R_t, down_up );
        a_man_phi_R     = periodic_bc_func( a_man_phi_R_t, down_up );
    
    else
        
        a_avg_theta_L   = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.theta_LR;
        a_avg_eta_L     = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.eta_LR;
        a_avg_phi_L     = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.phi_LR;
        a_avg_theta_R   = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.theta_LR;
        a_avg_eta_R     = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.eta_LR;
        a_avg_phi_R     = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_avg.phi_LR;        
        
        a_man_theta_L_t = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_fit.theta_L;
        a_man_eta_L_t   = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_fit.eta_L;
        a_man_phi_L_t   = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_fit.phi_L;
        a_man_theta_R_t = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_fit.theta_R;
        a_man_eta_R_t   = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_fit.eta_R;
        a_man_phi_R_t   = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_fit.phi_R;
        
        down_up         = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_fit.down_up;
        f               = pathDB.maneuver.(char(maneuver_type)).(char(wb_name)).a_fit.f;
        
        % Create periodic boundary conditions
    
        a_man_theta_L   = periodic_bc_func( a_man_theta_L_t, down_up );
        a_man_eta_L     = periodic_bc_func( a_man_eta_L_t, down_up );
        a_man_phi_L     = periodic_bc_func( a_man_phi_L_t, down_up );
        a_man_theta_R   = periodic_bc_func( a_man_theta_R_t, down_up );
        a_man_eta_R     = periodic_bc_func( a_man_eta_R_t, down_up );
        a_man_phi_R     = periodic_bc_func( a_man_phi_R_t, down_up );
        
    end
        
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
    
    wing_model.length           = pathDB.wing_model.length(seq_nr);
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
    body_IC.fg_fact = 1;
    body_IC.cg_fact = 1;

    
    clear sim_data
    
    % Tolerance settings ode45:
    
    options = odeset('RelTol',1e-5,'AbsTol',[1e-15 1e-15 1e-15 1e-15 1e-15 1e-15]);
    
    % case_nr == 0: Inertia only
    % case_nr == 1: Inertia, translational aerodynamics, gravity
    % case_nr == 2: Inertia, translational and rotational aerodynamics,
    % gravity
    
    case_nr = 2;
    
%     tic
%     
%     sim_data = dyn_sim(t_func,kine_avg,body_model,wing_model,body_IC,options,case_nr);
%     
%     toc
    
    n = 0;
    
    tic
    
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
    
%     [n] = plot_dyn_sim( sim_data, case_nr, n+1 );
    
    clear body_IC case_nr options
    
    
    % Maneuvering wingbeat
    
    % Obtain raw body kinematics:
    
    dt = pathDB.dt;
    
    wb_range = pathDB.wingbeats.wingbeat_loc(wb_nr,1,seq_nr):pathDB.wingbeats.wingbeat_loc(wb_nr,2,seq_nr);
    wb_length = length(wb_range);
    
    xyz_raw                 = pathDB.filt.xyz(wb_range,:,seq_nr)';
    qb_raw                  = pathDB.filt.qB(wb_range,:,seq_nr)';
    Rb_raw                  = pathDB.rot_mat.RB(:,:,wb_range,seq_nr);
    v_inertial_raw          = pathDB.filt.uvw(wb_range,:,seq_nr)';
    wb_raw                  = pathDB.filt.wB(wb_range,:,seq_nr)';
    a_inertial_raw          = pathDB.filt.a_xyz(wb_range,:,seq_nr)';
    
    w_dot_b_raw             = gradient(wb_raw)./dt;
    
    vb_raw                  = zeros(3,wb_length);
    w_inertial_raw          = zeros(3,wb_length);
    ab_raw                  = zeros(3,wb_length);
    w_dot_inertial_raw      = zeros(3,wb_length);
    
    for i = 1:wb_length
        
        vb_raw(:,i)                 = Rb_raw(:,:,i)*v_inertial_raw(:,i);
        w_inertial_raw(:,i)         = Rb_raw(:,:,i)'*wb_raw(:,i);
        ab_raw(:,i)                 = Rb_raw(:,:,i)*a_inertial_raw(:,i);
        w_dot_inertial_raw(:,i)     = Rb_raw(:,:,i)'*w_dot_b_raw(:,i);
        
    end
    
    v_strk_raw          = R_strk*mean(vb_raw,2)
    w_strk_raw          = R_strk*mean(wb_raw,2)
    a_strk_raw          = R_strk*mean(ab_raw,2)
    w_dot_strk_raw      = R_strk*mean(w_dot_b_raw,2)
    
   
    % Obtain start body kinematics:
    
    body_IC.xyz_0   = [0; 0; 0];
    qb_start        = qb_raw(:,1);
    body_IC.qb_0    = qb_start;
    Rb_start        = Rb_raw(:,:,1);
    vb_start        = vb_raw(:,1);
    wb_start        = wb_raw(:,1);
    body_IC.vb_0    = Rb_start*vb_start;
    body_IC.wb_0    = wb_start;
    body_IC.F_0     = sim_data.F_0;
    body_IC.M_0     = sim_data.M_0;
%     body_IC.F_0     = [0; 0; 0];
%     body_IC.M_0     = [0; 0; 0];
%     body_IC.fg_fact = fg_fact;

    clear sim_data
    
    %Plot lay-out  maneuvering wing kinematics:

%     n = n+1;
%     
%     figure()
%     Fly_plot_3D_wingtip( body_IC.xyz_0, Rb_start, kine_man.RL, kine_man.RR, down_up, body_model, wing_model, 1 )

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
    
    n = 0;
    
    [n] = plot_dyn_sim( sim_data, case_nr, n+1 );
    
    clear case_nr options
    
    % Create average accelerations, forces and moments in inertial
    % reference frame:
    
    T_sim               = sim_data.T;
    qb_sim              = sim_data.qb;
    vb_sim              = sim_data.vb;
    wb_sim              = sim_data.wb;
    ab_sim              = sim_data.ab;
    w_dot_b_sim         = sim_data.w_dot_b;
    
    FI_vel_b_sim        = sim_data.FI_vel_b;
    MI_vel_b_sim        = sim_data.MI_vel_b;
    FI_acc_b_sim        = sim_data.FI_acc_b;
    MI_acc_b_sim        = sim_data.MI_acc_b;
    FA_b_sim            = sim_data.FA_b;
    MA_b_sim            = sim_data.MA_b;
    Fg_b_sim            = sim_data.Fg_b;
    
    v_strk_sim          = R_strk*sim_data.vb_mean
    w_strk_sim          = R_strk*sim_data.wb_mean
    a_strk_sim          = R_strk*sim_data.ab_mean
    w_dot_strk_sim      = R_strk*sim_data.w_dot_b_mean
    
    FA_strk_sim         = sim_data.FA_strk_mean
    MA_strk_sim         = sim_data.MA_strk_mean
    FI_vel_strk_sim     = sim_data.FI_vel_strk_mean
    MI_vel_strk_sim     = sim_data.MI_vel_strk_mean
    FI_acc_strk_sim     = sim_data.FI_acc_strk_mean
    MI_acc_strk_sim     = sim_data.MI_acc_strk_mean
    Fg_strk_sim         = sim_data.Fg_strk_mean

     
    sim_length = length(T_sim);
    
    % Convert to inertial reference frame:
    
    Rb_sim              = zeros(3,3,sim_length);
    
    v_inert_sim         = zeros(3,sim_length);
    w_inert_sim         = zeros(3,sim_length);
    a_inert_sim         = zeros(3,sim_length);
    w_dot_inert_sim     = zeros(3,sim_length);
    
    FI_vel_inert_sim    = zeros(3,sim_length);
    MI_vel_inert_sim    = zeros(3,sim_length);
    FI_acc_inert_sim    = zeros(3,sim_length);
    MI_acc_inert_sim    = zeros(3,sim_length);
    FA_inert_sim        = zeros(3,sim_length);
    MA_inert_sim        = zeros(3,sim_length);
    Fg_inert_sim        = zeros(3,sim_length);
    
    
    for i = 1:sim_length
        
        Rb_sim(:,:,i)       = quat2mat(qb_sim(:,i));
        
        v_inert_sim(:,i)    = Rb_sim(:,:,i)'*vb_sim(:,i);
        w_inert_sim(:,i)    = Rb_sim(:,:,i)'*wb_sim(:,i);
        a_inert_sim(:,i)    = Rb_sim(:,:,i)'*ab_sim(:,i);
        w_dot_inert_sim(:,i)= Rb_sim(:,:,i)'*w_dot_b_sim(:,i);
        
        FI_vel_inert_sim(:,i)    = Rb_sim(:,:,i)'*FI_vel_b_sim(:,i);
        MI_vel_inert_sim(:,i)    = Rb_sim(:,:,i)'*MI_vel_b_sim(:,i);
        FI_acc_inert_sim(:,i)    = Rb_sim(:,:,i)'*FI_acc_b_sim(:,i);
        MI_acc_inert_sim(:,i)    = Rb_sim(:,:,i)'*MI_acc_b_sim(:,i);
        FA_inert_sim(:,i)        = Rb_sim(:,:,i)'*(FA_b_sim(:,i)-body_IC.F_0);
        MA_inert_sim(:,i)        = Rb_sim(:,:,i)'*(MA_b_sim(:,i)-body_IC.M_0);
        Fg_inert_sim(:,i)        = Rb_sim(:,:,i)'*Fg_b_sim(:,i);
        
    end
    
%     figure()
%     plot(T_sim,a_inert_sim)
%     
%     figure()
%     plot(T_sim,w_dot_inert_sim)
%     
%     figure()
%     plot(T_sim,FI_vel_inert_sim)
%     
%     figure()
%     plot(T_sim,MI_vel_inert_sim)
%     
%     figure()
%     plot(T_sim,FI_acc_inert_sim)
%     
%     figure()
%     plot(T_sim,MI_acc_inert_sim)
    
    
    % Obtain averages in the inertial reference frame:
    
    T_end = T_sim(sim_length);
    
    t_func = 0:(T_end/99):T_end;
    
    
    v_inert_sim_mean        = FM_mean_func(T_sim,t_func,v_inert_sim);
    w_inert_sim_mean        = FM_mean_func(T_sim,t_func,w_inert_sim);
    a_inert_sim_mean        = FM_mean_func(T_sim,t_func,a_inert_sim);
    w_dot_inert_sim_mean    = FM_mean_func(T_sim,t_func,w_dot_inert_sim);
    
    FI_vel_inert_sim_mean   = FM_mean_func(T_sim,t_func,FI_vel_inert_sim);
    MI_vel_inert_sim_mean   = FM_mean_func(T_sim,t_func,MI_vel_inert_sim);
    FI_acc_inert_sim_mean   = FM_mean_func(T_sim,t_func,FI_acc_inert_sim);
    MI_acc_inert_sim_mean   = FM_mean_func(T_sim,t_func,MI_acc_inert_sim);
    FA_inert_sim_mean       = FM_mean_func(T_sim,t_func,FA_inert_sim);
    MA_inert_sim_mean       = FM_mean_func(T_sim,t_func,MA_inert_sim);
    Fg_inert_sim_mean       = FM_mean_func(T_sim,t_func,Fg_inert_sim); 
    
%     figure()
%     hold on
%     subplot(4,1,1); hold on
%     plot(0:dt:(dt*(wb_length-1)),qb_raw(1,:),'b')
%     plot(T_sim,qb_sim(1,:),'r')
%     hold off
%     subplot(4,1,2); hold on
%     plot(0:dt:(dt*(wb_length-1)),qb_raw(2,:),'b')
%     plot(T_sim,qb_sim(2,:),'r')
%     hold off
%     subplot(4,1,3); hold on
%     plot(0:dt:(dt*(wb_length-1)),qb_raw(3,:),'b')
%     plot(T_sim,qb_sim(3,:),'r')
%     hold off
%     subplot(4,1,4); hold on
%     plot(0:dt:(dt*(wb_length-1)),qb_raw(4,:),'b')
%     plot(T_sim,qb_sim(4,:),'r')
%     hold off
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); hold on
%     plot(0:dt:(dt*(wb_length-1)),w_inertial_raw(1,:),'b')
%     plot(T_sim,w_inert_sim(1,:),'r')
%     plot(0:dt:(dt*(wb_length-1)),ones(wb_length,1)*w_inert_sim_mean(1),'k')
%     hold off
%     subplot(3,1,2); hold on
%     plot(0:dt:(dt*(wb_length-1)),w_inertial_raw(2,:),'b')
%     plot(T_sim,w_inert_sim(2,:),'r')
%     plot(0:dt:(dt*(wb_length-1)),ones(wb_length,1)*w_inert_sim_mean(2),'k')
%     hold off
%     subplot(3,1,3); hold on
%     plot(0:dt:(dt*(wb_length-1)),w_inertial_raw(3,:),'b')
%     plot(T_sim,w_inert_sim(3,:),'r')
%     plot(0:dt:(dt*(wb_length-1)),ones(wb_length,1)*w_inert_sim_mean(3),'k')
%     hold off
%     hold off
    
    % Raw data:
    
    a_inertial_raw_mean     = mean(a_inertial_raw,2);
    w_dot_inertial_raw_mean = mean(w_dot_inertial_raw,2);
    
    % Export data:
    
    F_sim.FA                = FA_inert_sim_mean;
    F_sim.FI_vel            = FI_vel_inert_sim_mean;
    F_sim.FI_acc            = FI_acc_inert_sim_mean;
    F_sim.Fg                = Fg_inert_sim_mean;
    
    F_sim.MA                = MA_inert_sim_mean;
    F_sim.MI_vel            = MI_vel_inert_sim_mean;
    F_sim.MI_acc            = MI_acc_inert_sim_mean;
    
    kine_sim.a_xyz          = a_inert_sim_mean;
    kine_raw.a_xyz          = a_inertial_raw_mean;
    kine_sim.w_dot          = w_dot_inert_sim_mean;
    kine_raw.w_dot          = w_dot_inertial_raw_mean;
    
    figure(3)
    Fly_FM_plot( Rb_start, R_strk ,body_model, wing_model, F_sim, kine_sim, kine_raw )
    set(3,'Color','k')
    set(gcf, 'InvertHardCopy', 'off');
    
end

