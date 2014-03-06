% initialize DB variable

% Computation variables

frames = [1:settings.frame_end]';
pathDB.t = (frames-settings.trigger_frame)/settings.fps;

pathDB.x = [];
pathDB.y = [];
pathDB.z = [];

pathDB.x_filt = [];
pathDB.y_filt = [];
pathDB.z_filt = [];

pathDB.u_filt = [];
pathDB.v_filt = [];
pathDB.w_filt = [];

pathDB.ax_filt = [];
pathDB.ay_filt = [];
pathDB.az_filt = [];

pathDB.u_body = [];
pathDB.v_body = [];
pathDB.w_body = [];

pathDB.ax_body = [];
pathDB.ay_body = [];
pathDB.az_body = [];

pathDB.qb1 = [];
pathDB.qb2 = [];
pathDB.qb3 = [];
pathDB.qb4 = [];

pathDB.qL1 = [];
pathDB.qL2 = [];
pathDB.qL3 = [];
pathDB.qL4 = [];

pathDB.qR1 = [];
pathDB.qR2 = [];
pathDB.qR3 = [];
pathDB.qR4 = [];

pathDB.qb1_filt = [];
pathDB.qb2_filt = [];
pathDB.qb3_filt = [];
pathDB.qb4_filt = [];

pathDB.b_omega1 = [];
pathDB.b_omega2 = [];
pathDB.b_omega3 = [];

pathDB.qL1_filt1 = [];
pathDB.qL2_filt1 = [];
pathDB.qL3_filt1 = [];
pathDB.qL4_filt1 = [];

pathDB.qR1_filt1 = [];
pathDB.qR2_filt1 = [];
pathDB.qR3_filt1 = [];
pathDB.qR4_filt1 = [];

pathDB.omega1_L = [];
pathDB.omega2_L = [];
pathDB.omega3_L = [];

pathDB.qL1_filt2 = [];
pathDB.qL2_filt2 = [];
pathDB.qL3_filt2 = [];
pathDB.qL4_filt2 = [];

pathDB.omega1_R = [];
pathDB.omega2_R = [];
pathDB.omega3_R = [];

pathDB.qR1_filt2 = [];
pathDB.qR2_filt2 = [];
pathDB.qR3_filt2 = [];
pathDB.qR4_filt2 = [];

pathDB.b_roll = [];
pathDB.b_pitch = [];
pathDB.b_yaw = [];
pathDB.b_alfa = [];
pathDB.b_beta = [];

pathDB.body_l = [];
pathDB.wing_l = [];
pathDB.joint_pos_L = [];
pathDB.joint_pos_R = [];

pathDB.u_wing_L = [];
pathDB.v_wing_L = [];
pathDB.w_wing_L = [];

pathDB.u_wing_R = [];
pathDB.v_wing_R = [];
pathDB.w_wing_R = [];

pathDB.u_wing_rot_L = [];
pathDB.v_wing_rot_L = [];
pathDB.w_wing_rot_L = [];
    
pathDB.u_wing_rot_R = [];
pathDB.v_wing_rot_R = [];
pathDB.w_wing_rot_R = [];
    
pathDB.u_wing_trans_L = [];
pathDB.v_wing_trans_L = [];
pathDB.w_wing_trans_L = [];
    
pathDB.u_wing_trans_R = [];
pathDB.v_wing_trans_R = [];
pathDB.w_wing_trans_R = [];
    
pathDB.u_wing_body_rot_L = [];
pathDB.v_wing_body_rot_L = [];
pathDB.w_wing_body_rot_L = [];
    
pathDB.u_wing_body_rot_R = [];
pathDB.v_wing_body_rot_R = [];
pathDB.w_wing_body_rot_R = [];
    
pathDB.alfa_L = [];
pathDB.beta_L = [];
    
pathDB.alfa_R = [];
pathDB.beta_R = [];

pathDB.Lwingtip = [];
pathDB.Rwingtip = [];

pathDB.L_wingbeat_loc = [];
pathDB.R_wingbeat_loc = [];

pathDB.phi_L = [];
pathDB.theta_L = [];
pathDB.eta_L = [];

pathDB.phi_R = [];
pathDB.theta_R = [];
pathDB.eta_R = [];

% Histogram data:

    % Body:
    
    pathDB.phi_body_mean = [];
    pathDB.theta_body_mean = [];
    pathDB.xsi_body_mean = [];
    
    pathDB.phi_body_sd = [];
    pathDB.theta_body_sd = [];
    pathDB.xsi_body_sd = [];
    
    pathDB.omegax_body_mean = [];
    pathDB.omegay_body_mean = [];
    pathDB.omegaz_body_mean = [];
    pathDB.Omega_body_mean = [];
    
    pathDB.omegax_body_sd = [];
    pathDB.omegay_body_sd = [];
    pathDB.omegaz_body_sd = [];
    pathDB.Omega_body_sd = [];
    
    pathDB.alfa_body_mean = [];
    pathDB.beta_body_mean = [];
    
    pathDB.alfa_body_sd = [];
    pathDB.beta_body_sd = [];
    
    pathDB.u_body_mean = [];
    pathDB.v_body_mean = [];
    pathDB.w_body_mean = [];
    pathDB.U_body_mean = [];
    
    pathDB.u_body_sd = []; 
    pathDB.v_body_sd = [];
    pathDB.w_body_sd = [];
    pathDB.U_body_sd = [];
    
    pathDB.ax_body_mean = [];
    pathDB.ay_body_mean = [];
    pathDB.az_body_mean = [];
    pathDB.a_body_mean = [];
    
    pathDB.ax_body_sd = [];
    pathDB.ay_body_sd = [];
    pathDB.az_body_sd = [];
    pathDB.a_body_sd = [];
    
    % Left wing:
    
    pathDB.phi_L_up_mean = []; 
    pathDB.theta_L_up_mean = [];
    pathDB.eta_L_up_mean = [];
    
    pathDB.phi_L_up_sd = [];
    pathDB.theta_L_up_sd = [];
    pathDB.eta_L_up_sd = [];
    
    pathDB.phi_L_down_mean = [];
    pathDB.theta_L_down_mean = [];
    pathDB.eta_L_down_mean = [];
    
    pathDB.phi_L_down_sd = [];
    pathDB.theta_L_down_sd = [];
    pathDB.eta_L_down_sd = [];
    
    pathDB.alfa_L_up_mean = []; 
    pathDB.beta_L_up_mean = [];
    
    pathDB.alfa_L_up_sd = [];
    pathDB.beta_L_up_sd = [];
    
    pathDB.alfa_L_down_mean = [];
    pathDB.beta_L_down_mean = [];
    
    pathDB.alfa_L_down_sd = [];
    pathDB.beta_L_down_sd = [];
    
    pathDB.u_L_up_mean = [];
    pathDB.v_L_up_mean = [];
    pathDB.w_L_up_mean = [];
    pathDB.U_L_up_mean = [];
    
    pathDB.u_L_up_sd = [];
    pathDB.v_L_up_sd = [];
    pathDB.w_L_up_sd = [];
    pathDB.U_L_up_sd = [];
    
    pathDB.u_L_down_mean = [];
    pathDB.v_L_down_mean = [];
    pathDB.w_L_down_mean = [];
    pathDB.U_L_down_mean = [];
    
    pathDB.u_L_down_sd = []; 
    pathDB.v_L_down_sd = [];
    pathDB.w_L_down_sd = [];
    pathDB.U_L_down_sd = [];
    
    pathDB.omegax_L_up_mean = [];
    pathDB.omegay_L_up_mean = [];
    pathDB.omegaz_L_up_mean = [];
    pathDB.Omega_L_up_mean = [];
    
    pathDB.omegax_L_up_sd = [];
    pathDB.omegay_L_up_sd = [];
    pathDB.omegaz_L_up_sd = [];
    pathDB.Omega_L_up_sd = [];
    
    pathDB.omegax_L_down_mean = [];
    pathDB.omegay_L_down_mean = [];
    pathDB.omegaz_L_down_mean = [];
    pathDB.Omega_L_down_mean = [];
    
    pathDB.omegax_L_down_sd = [];
    pathDB.omegay_L_down_sd = [];
    pathDB.omegaz_L_down_sd = [];
    pathDB.Omega_L_down_sd = [];
    
    
    % Right wing:
    
    pathDB.phi_R_up_mean = [];
    pathDB.theta_R_up_mean = [];
    pathDB.eta_R_up_mean = [];
    
    pathDB.phi_R_up_sd = [];
    pathDB.theta_R_up_sd = [];
    pathDB.eta_R_up_sd = [];
    
    pathDB.phi_R_down_mean = [];
    pathDB.theta_R_down_mean = [];
    pathDB.eta_R_down_mean = [];
    
    pathDB.phi_R_down_sd = [];
    pathDB.theta_R_down_sd = [];
    pathDB.eta_R_down_sd = [];
    
    pathDB.alfa_R_up_mean = [];
    pathDB.beta_R_up_mean = [];
    
    pathDB.alfa_R_up_sd = [];
    pathDB.beta_R_up_sd = [];
    
    pathDB.alfa_R_down_mean = [];
    pathDB.beta_R_down_mean = [];
    
    pathDB.alfa_R_down_sd = [];
    pathDB.beta_R_down_sd = [];
    
    pathDB.u_R_up_mean = [];
    pathDB.v_R_up_mean = [];
    pathDB.w_R_up_mean = [];
    pathDB.U_R_up_mean = [];
    
    pathDB.u_R_up_sd = [];
    pathDB.v_R_up_sd = [];
    pathDB.w_R_up_sd = [];
    pathDB.U_R_up_sd = [];
    
    pathDB.u_R_down_mean = [];
    pathDB.v_R_down_mean = [];
    pathDB.w_R_down_mean = [];
    pathDB.U_R_down_mean = [];
    
    pathDB.u_R_down_sd = [];
    pathDB.v_R_down_sd = [];
    pathDB.w_R_down_sd = [];
    pathDB.U_R_down_sd = [];
    
    pathDB.omegax_R_up_mean = [];
    pathDB.omegay_R_up_mean = [];
    pathDB.omegaz_R_up_mean = [];
    pathDB.Omega_R_up_mean = [];
    
    pathDB.omegax_R_up_sd = [];
    pathDB.omegay_R_up_sd = [];
    pathDB.omegaz_R_up_sd = [];
    pathDB.Omega_R_up_sd = [];
    
    pathDB.omegax_R_down_mean = [];
    pathDB.omegay_R_down_mean = [];
    pathDB.omegaz_R_down_mean = [];
    pathDB.Omega_R_down_mean = [];
    
    pathDB.omegax_R_down_sd = [];
    pathDB.omegay_R_down_sd = [];
    pathDB.omegaz_R_down_sd = [];
    pathDB.Omega_R_down_sd = [];

    
    % PathDB7
   
    pathDB.F_joint_L_down = [];
    pathDB.M_joint_L_down = [];
    
    pathDB.F_joint_L_up = [];
    pathDB.M_joint_L_up = [];
    
    pathDB.F_joint_R_down = [];
    pathDB.M_joint_R_down = [];
    
    pathDB.F_joint_R_up = [];
    pathDB.M_joint_R_up = [];
    
    pathDB.F_cg = [];
    pathDB.M_cg = [];
    
    pathDB.R_L_down = [];
    pathDB.Fn_L_down = [];
    pathDB.Ft_L_down = [];
    pathDB.Mn_L_down = [];
    pathDB.Mt_L_down = [];
    
    pathDB.R_L_up = [];
    pathDB.Fn_L_up = [];
    pathDB.Ft_L_up = [];
    pathDB.Mn_L_up = [];
    pathDB.Mt_L_up = [];
    
    pathDB.R_R_down = [];
    pathDB.Fn_R_down = [];
    pathDB.Ft_R_down = [];
    pathDB.Mn_R_down = [];
    pathDB.Mt_R_down = [];
    
    pathDB.R_R_up = [];
    pathDB.Fn_R_up = [];
    pathDB.Ft_R_up = [];
    pathDB.Mn_R_up = [];
    pathDB.Mt_R_up = [];
      
    pathDB.m_est_point = [];
    
    pathDB.m_est_dynamical_model = [];
    
    pathDB.r_body_L_down = [];
    pathDB.r_body_L_up = [];
    
    pathDB.r_body_R_down = [];
    pathDB.r_body_R_up = [];
    
    pathDB.q_avg_down_body = [];
    pathDB.q_avg_up_body = [];
    
    pathDB.q_avg_down_L = [];
    pathDB.q_avg_up_L = [];
    
    pathDB.q_avg_down_R = [];
    pathDB.q_avg_up_R = [];
    
    pathDB.down_time_L = [];
    pathDB.up_time_L = [];
        
    pathDB.down_time_R= [];
    pathDB.up_time_R= [];
    
    pathDB.wingbeat_time= [];
    
    pathDB.F_L = [];
    pathDB.F_R = [];
    
    pathDB.F_L_arm = [];
    pathDB.F_R_arm = [];  