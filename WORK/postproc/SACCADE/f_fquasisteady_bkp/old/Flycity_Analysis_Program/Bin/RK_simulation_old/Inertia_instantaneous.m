function [ FM_b, FM_strkpln ] = Inertia_instantaneous( kine, body_model, wing_model )


    % Compute the instantaneous inertial forces and moments on the center
    % of gravity in the strokeplane reference frame:
    
    RL              = kine.RL;
    RR              = kine.RR;
%     wL              = kine.wL;
%     wR              = kine.wR;
%     w_dot_L         = kine.w_dot_L;
%     w_dot_R         = kine.w_dot_R;
    wL              = kine.wL_b;
    wR              = kine.wR_b;
    w_dot_L         = kine.w_dot_L_b;
    w_dot_R         = kine.w_dot_R_b;
%     wL_b            = kine.wL_b;
%     wR_b            = kine.wR_b;
%     w_dot_L_b       = kine.w_dot_L_b;
%     w_dot_R_b       = kine.w_dot_R_b;
    vb              = kine.vb;
    wb              = kine.wb;
    ab              = kine.ab;
    w_dot_b         = kine.w_dot_b;
    R_strk          = kine.R_strk;
    
    body_kin.vb     = vb;
    body_kin.wb     = wb;
    
    wing_kin.RL     = RL;
    wing_kin.RR     = RR;
    wing_kin.wL     = wL;
    wing_kin.wR     = wR;
    
    
    % Compute the mass matrix:
    
    [ M_matrix ] = Mass_matrix_new( body_model, wing_model, body_kin, wing_kin );

    
    M11 = M_matrix.M11;
    M12 = M_matrix.M12;
    M13 = M_matrix.M13;
    M14 = M_matrix.M14;
    M21 = M_matrix.M21;
    M22 = M_matrix.M22;
    M23 = M_matrix.M23;
    M24 = M_matrix.M24;
    M31 = M_matrix.M31;
    M32 = M_matrix.M32;
    M33 = M_matrix.M33;
    M34 = M_matrix.M34;
    M41 = M_matrix.M41;
    M42 = M_matrix.M42;
    M43 = M_matrix.M43;
    M44 = M_matrix.M44;
    
    M   = M_matrix.M;
    
    M_dot_12 = M_matrix.M12_dot;
    M_dot_13 = M_matrix.M13_dot;
    M_dot_14 = M_matrix.M14_dot;
    M_dot_21 = M_matrix.M21_dot;
    M_dot_22 = M_matrix.M22_dot;
    M_dot_23 = M_matrix.M23_dot;
    M_dot_24 = M_matrix.M24_dot;
    
    w_b_dL_dv_b = M_matrix.Wb_dL_dvb;
    v_b_dL_dv_b = M_matrix.Vb_dL_dvb;
    w_b_dL_dw_b = M_matrix.Wb_dL_dwb;
    
%     M13
%     M_dot_13
%     M14
%     M_dot_14
%     
%     M23
%     M_dot_23
%     M24
%     M_dot_24
%     
%     pause
    

    % Compute the intertial forces and moments:
    
    F_acc_b = 1e-3*[M11 M12]*[ab; w_dot_b]; 
    
%     F_vel_b = 1e-3*(-w_b_dL_dv_b-M_dot_12*wb-M_dot_13*wL- ...
%                     M13*w_dot_L-M_dot_14*wR-M14*w_dot_R);

%     F_vel_b = 1e-3*(-M_dot_13*wL- ...
%                     M13*w_dot_L-M_dot_14*wR-M14*w_dot_R);

    F_vel_b = 1e-3*(-M_dot_13*wL- ...
                    M_dot_14*wR);
                
    F_vel_b_L = 1e-3*(-M_dot_13*wL-M13*w_dot_L);
    F_vel_b_R = 1e-3*(-M_dot_14*wR-M14*w_dot_R);

    M_acc_b = 1e-3*[M21 M22]*[ab; w_dot_b]; 

%     M_vel_b = 1e-3*(-v_b_dL_dv_b-w_b_dL_dw_b-M_dot_21*vb-M_dot_22*wb- ...
%                     M_dot_23*wL-M23*w_dot_L-M_dot_24*wR-M24*w_dot_R);

    M_vel_b_L = 1e-3*(-M_dot_23*wL-M23*w_dot_L);
    M_vel_b_R = 1e-3*(-M_dot_24*wR-M24*w_dot_R);
          
%     M_vel_b = 1e-3*(- ...
%                     M_dot_23*wL-M23*w_dot_L-M_dot_24*wR-M24*w_dot_R);

    M_vel_b = 1e-3*(- ...
                    M_dot_23*wL-M_dot_24*wR);



    FM_b.F_I_acc = F_acc_b;
    FM_b.M_I_acc = M_acc_b;
    FM_b.F_I_vel = F_vel_b;
    FM_b.M_I_vel = M_vel_b;
    FM_b.F_I_vel_L = F_vel_b_L;
    FM_b.F_I_vel_R = F_vel_b_R;
    FM_b.M_I_vel_L = M_vel_b_L;
    FM_b.M_I_vel_R = M_vel_b_R;
    FM_b.M_mat_b = [M11 M12; M21 M22];
    FM_b.LM      = M_matrix.Lin_momentum;
    FM_b.AM      = M_matrix.Ang_momentum;
    FM_b.KE      = 0.5*[vb; wb; wL; wR]'*M*[vb; wb; wL; wR];
    FM_b.KE_lin  = 0.5*vb'*M11*vb;
    FM_b.KE_ang  = 0.5*[wb; wL; wR]'*[ M22 M23 M24; M32 M33 M34; M42 M43 M44]*[wb; wL; wR];
    
    FM_strkpln.F_I_acc = R_strk*F_acc_b;
    FM_strkpln.M_I_acc = R_strk*M_acc_b;
    FM_strkpln.F_I_vel = R_strk*F_vel_b;
    FM_strkpln.M_I_vel = R_strk*M_vel_b;

end

