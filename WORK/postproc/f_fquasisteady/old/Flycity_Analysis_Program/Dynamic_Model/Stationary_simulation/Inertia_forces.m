function [ FM_strkpln ] = Inertia_forces( body_model, wing_model, kine )


    % Compute the forces and moments on the body due to body and wing
    % inertia:
    
    u_strk          = kine.u_strk;          % [ mm/s ]
    w_strk          = kine.w_strk;          % [ rad/s ]
    a_strk          = kine.
    wL              = kine.wL;              % [ rad/s ]
    wR              = kine.wR;              % [ rad/s ]
    RL              = kine.RL;
    RR              = kine.RR;
    R_strk          = kine.R_strk;
    
    % Convert from strokeplane reference frame to body reference frame:
    
    % 
    
    [ M_matrix ] = Mass_matrix( body_model, wing_model, body_kin, wing_kin )
    
    F_acc_b(:,i) = [M11 M12]*[a_b(:,i); w_dot_b(:,i)]; %[kg*mm/s^2]
    F_vel_b(:,i) = [w_b_dL_dv_b+M_dot_12*w_b(:,i)+M_dot_13*w_bL(:,i)+M13*w_dot_bL(:,i)+M_dot_14*w_bR(:,i)+M14*w_dot_bR(:,i)]; %[kg*mm/s^2]
    M_acc_b(:,i) = [M21 M22]*[a_b(:,i); w_dot_b(:,i)]; %[kg*mm^2/s^2]
    M_vel_b(:,i) = [v_b_dL_dv_b+w_b_dL_dw_b+M_dot_21*v_b(:,i)+M_dot_22*w_b(:,i)+M_dot_23*w_bL(:,i)+M_23*w_dot_bL(:,i)+M_dot_24*w_bR(:,i)+M_24*w_dot_bR(:,i)]; %[kg*mm^2/s^2]
    F_acc_st(:,i) = R_beta*F_acc_b(:,i);
    F_vel_st(:,i) = R_beta*F_vel_b(:,i);
    M_acc_st(:,i) = R_beta*M_acc_b(:,i);
    M_vel_st(:,i) = R_beta*M_vel_b(:,i);
                
    F_I_b(:,i) = F_acc_b(:,i)+F_vel_b(:,i);
    M_I_b(:,i) = M_acc_b(:,i)+M_vel_b(:,i);
                
    F_I_st(:,i) = R_beta*[F_acc_b(:,i)+F_vel_b(:,i)];         
    M_I_st(:,i) = R_beta*[M_acc_b(:,i)+M_vel_b(:,i)];

    F_I_st(:,i) = R_beta*[F_acc_b(:,i)];
    M_I_st(:,i) = R_beta*[M_acc_b(:,i)];



end

