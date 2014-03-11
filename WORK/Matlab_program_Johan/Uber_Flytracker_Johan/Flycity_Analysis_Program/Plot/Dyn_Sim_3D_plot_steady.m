function Dyn_Sim_3D_plot_steady( settings, pathDB, wb_name )


    R_strk = pathDB.rot_mat.Rstr;
    
    kine = pathDB.DynSim_steady.(wb_name).kine;
    
    sim_data = pathDB.DynSim_steady.(wb_name).sim_data;
    
    seq_nr      = pathDB.rand_wbs.(wb_name).seq_nr;
    
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
    wing_model.x_LE_L           = pathDB.wing_model.x_LE_L(seq_nr,:)';
    wing_model.y_sect_R         = pathDB.wing_model.y_sect_R(:,:,seq_nr)';
    wing_model.chords_R         = pathDB.wing_model.chords_R(seq_nr,:)';
    wing_model.x_LE_R           = pathDB.wing_model.x_LE_R(seq_nr,:)';
    wing_model.rho              = settings.rho_air;
    wing_model.length           = pathDB.wing_model.length(seq_nr);
    wing_model.x_mod_L          = pathDB.wing_model.x_mod_L(:,:,seq_nr);
    wing_model.y_mod_L          = pathDB.wing_model.y_mod_L(:,:,seq_nr);
    wing_model.z_mod_L          = pathDB.wing_model.z_mod_L(:,:,seq_nr);
    wing_model.x_mod_R          = pathDB.wing_model.x_mod_R(:,:,seq_nr);
    wing_model.y_mod_R          = pathDB.wing_model.y_mod_R(:,:,seq_nr);
    wing_model.z_mod_R          = pathDB.wing_model.z_mod_R(:,:,seq_nr);
    
    
    % Plot the start condition of the simulation:
    
    xyz_start       = sim_data.xyz(:,1);
    qb_start        = sim_data.qb(:,1);
    Rb_start        = quat2mat(qb_start)';
    vb_start        = sim_data.vb(:,1);
    wb_start        = sim_data.wb(:,1);
    ab_start        = sim_data.ab(:,1);
    w_dot_b_start   = sim_data.w_dot_b(:,1);
    FA_b_start      = sim_data.FA_b(:,1);
    MA_b_start      = sim_data.MA_b(:,1);
    FI_vel_b_start  = sim_data.FI_vel_b(:,1);
    MI_vel_b_start  = sim_data.MI_vel_b(:,1);
    FI_acc_b_start  = sim_data.FI_acc_b(:,1);
    MI_acc_b_start  = sim_data.MI_acc_b(:,1);
    Fg_b_start      = sim_data.Fg_b(:,2);
    
    norm(ab_start)
    
    ab_raw_start    = pathDB.rand_wbs.(wb_name).body_kin.a_xyz(1,:)';
    
    norm(ab_raw_start)
    
    figure()
    Fly_force_plot( Rb_start, R_strk, xyz_start, vb_start, ab_start, ab_raw_start, FA_b_start, FI_vel_b_start, FI_acc_b_start, Fg_b_start, body_model, wing_model )
    
    % Plot the average condition of the simulation:
    
    xyz_avg       = FM_mean_func(sim_data.T,kine.t,sim_data.xyz);
    qb_avg        = quat_mean( kine.t, sim_data.T, sim_data.qb );
    Rb_avg        = quat2mat(qb_avg)';
    vb_avg        = sim_data.vb_mean;
    wb_avg        = sim_data.wb_mean;
    ab_avg        = sim_data.ab_mean;
    w_dot_b_avg   = sim_data.w_dot_b_mean;
    FA_b_avg      = sim_data.FA_b_mean;
    MA_b_avg      = sim_data.MA_b_mean;
    FI_vel_b_avg  = sim_data.FI_vel_b_mean;
    MI_vel_b_avg  = sim_data.MI_vel_b_mean;
    FI_acc_b_avg  = sim_data.FI_acc_b_mean;
    MI_acc_b_avg  = sim_data.MI_acc_b_mean;
    Fg_b_avg      = sim_data.Fg_b_mean;
    
    norm(ab_avg)

    ab_raw_avg    = pathDB.rand_wbs.(wb_name).body_kin.a_xyz_mean';
    
    norm(ab_raw_avg)
    
    figure()
    Fly_force_plot( Rb_avg, R_strk, xyz_avg, vb_avg, ab_avg, ab_raw_avg, FA_b_avg, FI_vel_b_avg, FI_acc_b_avg, Fg_b_avg, body_model, wing_model )
    
end


