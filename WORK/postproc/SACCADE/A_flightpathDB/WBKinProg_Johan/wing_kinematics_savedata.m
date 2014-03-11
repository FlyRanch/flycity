function wing_kinematics_savedata(settings,pathDB)


    savefile = ['Wing_kinematics_data_strokeplane',num2str(-settings.strokeplane_WBkin),'deg.mat'];
    
    
    % Save wing kinematics data in a .mat file:
    
    seq_names = settings.sequence_names;
    
    phi_L = pathDB.phi_L;
    
    theta_L = pathDB.theta_L;
    
    eta_L = pathDB.eta_L;
    
    phi_R = pathDB.phi_R;
    
    theta_R = pathDB.theta_R;
    
    eta_R = pathDB.eta_R;
    
    wingbeat_time = pathDB.wingbeat_time;
    
    down_time_L = pathDB.down_time_L;
    
    up_time_L = pathDB.up_time_L;
    
    down_time_R = pathDB.down_time_R;
    
    up_time_R = pathDB.up_time_R;
    
    wingtip_path_L = pathDB.Lwingtip;
    
    wingtip_path_R = pathDB.Rwingtip;
    
    wing_length = pathDB.wing_l;
    
    joint_pos_L = pathDB.joint_pos_L;
    
    joint_pos_R = pathDB.joint_pos_R;
    
    u_wing_L = pathDB.u_wing_L;
    
    v_wing_L = pathDB.v_wing_L;
    
    w_wing_L = pathDB.w_wing_L;
    
    u_wing_R = pathDB.u_wing_R;
    
    v_wing_R = pathDB.v_wing_R;
    
    w_wing_R = pathDB.w_wing_R;
    
    alfa_L = pathDB.alfa_L;
    
    beta_L = pathDB.beta_L;
    
    alfa_R = pathDB.alfa_R;
    
    beta_R = pathDB.beta_R;
      
    
    
    
    save(savefile, 'seq_names', 'phi_L', 'theta_L', 'eta_L', 'phi_R', 'theta_R', 'eta_R', 'wingbeat_time', 'down_time_L', 'up_time_L', 'down_time_R', 'up_time_R' , 'wingtip_path_L', 'wingtip_path_R', 'wing_length', 'joint_pos_L', 'joint_pos_R', ...
        'u_wing_L', 'v_wing_L', 'w_wing_L', 'u_wing_R', 'v_wing_R', 'w_wing_R', 'alfa_L', 'beta_L', 'alfa_R', 'beta_R')


end

