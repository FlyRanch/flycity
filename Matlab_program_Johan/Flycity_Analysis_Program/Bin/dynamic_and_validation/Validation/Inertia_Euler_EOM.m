function [ FM_b, FM_strkpln ] = Inertia_Euler_EOM( kine, body_model, wing_model )

    % Give the instantaneous inertial forces and moments:
    
    vb              = kine.vb;
    wb              = kine.wb;
    ab              = kine.ab;
    w_dot_b         = kine.
    
    theta_L         = kine.theta_L;
    eta_L           = kine.eta_L;
    phi_L           = kine.phi_L;
    theta_R         = kine.theta_R;
    eta_R           = kine.eta_R;
    phi_R           = kine.phi_R;
    theta_L         = kine.theta_dot_L;
    eta_L           = kine.eta_dot_L;
    phi_L           = kine.phi_dot_L;
    theta_R         = kine.theta_dot_R;
    eta_R           = kine.eta_dot_R;
    phi_R           = kine.phi_dot_R;
    wL_b            = kine.wL_b;
    wR_b            = kine.wR_b;
    w_dot_L_b       = kine.w_dot_L_b;
    w_dot_R_b       = kine.w_dot_R_b;
    R_strk          = kine.R_strk;
    
    JL              = body_model.
    JR              = body_model.
    
    wing_cg_L       = wing_model.wing_cg_L;
    wing_cg_R       = wing_model.wing_cg_R;


    % Compute the intertial forces and moments:
    
    F_vel_b = 1e-3*(-w_b_dL_dv_b-M_dot_12*wb-M_dot_13*wL_b- ...
                    M13*w_dot_L_b-M_dot_14*wR_b-M14*w_dot_R_b);     %[N]

    M_vel_b = 1e-3*(-v_b_dL_dv_b-w_b_dL_dw_b-M_dot_21*vb-M_dot_22*wb- ...
                    M_dot_23*wL_b-M23*w_dot_L_b-M_dot_24*wR_b-M24*w_dot_R_b); %[N*mm]

    FM_b.M_mat_b = [M11 M12; M21 M22];
    FM_b.F_I_vel = F_vel_b;
    FM_b.M_I_vel = M_vel_b;
    FM_b.LM      = M_matrix.Lin_momentum;
    FM_b.AM      = M_matrix.Ang_momentum;
    FM_b.KE      = 0.5*[vb; wb; wL_b; wR_b]'*M*[vb; wb; wL_b; wR_b];
    FM_b.KE_lin  = 0.5*vb'*[M11 M12 M13 M14]*[vb; wb; wL_b; wR_b];
    FM_b.KE_ang  = 0.5*[wb; wL_b; wR_b]'*[ M12 M22 M23 M24; M31 M32 M33 M34; M41 M42 M43 M44]*[vb; wb; wL_b; wR_b];
    FM_b.vb_0    = M_matrix.vb_0;
    FM_b.wb_0    = M_matrix.wb_0;
    
    FM_strkpln.F_I_vel = R_strk*F_vel_b;
    FM_strkpln.M_I_vel = R_strk*M_vel_b;
    
        
end

