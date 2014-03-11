function [ Inertia_f ] = Inertia_forces( wing_kin, Body, n, v_st, a_st, w_st, w_dot_st)

    % Import body parameters and wing kinematics:
    
    m_fly = Body.m_wl; % kg    
    m_vw = Body.m_v_w; % kg    
    wing_l = Body.wing_l; % mm    
    JL = Body.joint_L; % mm    
    JR = Body.joint_R; % mm    
    cg_body = Body.cg_body; % mm
    cg_L = Body.cg_L; % mm    
    cg_R = Body.cg_R; % mm    
    I_body = Body.I_body; % kg/mm^2    
    I_v_wing = Body.I_v_wing; % kg/mm^2    
    
    w_L = wing_kin.w_L; % rad/s
    w_R = wing_kin.w_R; % rad/s
    w_dot_L = wing_kin.w_dot_L; % rad/s^2
    w_dot_R = wing_kin.w_dot_R; % rad/s^2
    R_L = wing_kin.Rot_L;
    R_R = wing_kin.Rot_R;
    
   
    % Rotation matrix strokeplane reference frame:
    
    beta = (55/180)*pi;
    
    R_beta = [cos(-beta) 0 -sin(-beta); ...
              0 1 0; ...
              sin(-beta) 0 cos(-beta)]; 
    
    
    % Define joint arms and wing cg arms and wing angular velocity/acceleration in
    % strokeplane reference frame:
    
    arm_JL = JL-cg_body;
    arm_JR = JR-cg_body;
    
    arm_cg_L = zeros(3,n);
    arm_cg_R = zeros(3,n);
    w_bL = zeros(3,n);
    w_bR = zeros(3,n);
    w_dot_bL = zeros(3,n);
    w_dot_bR = zeros(3,n);
    
    for i = 1:n
    
        arm_cg_L(:,i) = R_L(:,:,i)'*cg_L;
        arm_cg_R(:,i) = R_R(:,:,i)'*cg_R;
        w_bL(:,i) = R_L(:,:,i)'*w_L(:,i);
        w_bR(:,i) = R_R(:,:,i)'*w_R(:,i);
        w_dot_bL(:,i) = R_L(:,:,i)'*w_dot_L(:,i);
        w_dot_bR(:,i) = R_R(:,:,i)'*w_dot_R(:,i);
    
    end
    
    % Compute the inertia for the current wingbeat sequence:
            
    F_acc_b = zeros(3,n);
    F_vel_b = zeros(3,n);      
    M_acc_b = zeros(3,n);
    M_vel_b = zeros(3,n);
    
    F_acc_st = zeros(3,n);
    F_vel_st = zeros(3,n);      
    M_acc_st = zeros(3,n);
    M_vel_st = zeros(3,n);
            
    F_I_b = zeros(3,n);
    M_I_b = zeros(3,n);
    
    F_I_st = zeros(3,n);
    M_I_st = zeros(3,n);
    
    Lin_momentum_b = zeros(3,n);
    
    Ang_momentum_b = zeros(3,n);
    
    Lin_momentum_st = zeros(3,n);
    
    Ang_momentum_st = zeros(3,n);
    
    % Body kinematics in body reference frame:
            
    v_b = zeros(3,n);
    a_b = zeros(3,n);
    w_b = zeros(3,n);
    w_dot_b = zeros(3,n);
            
    for i = 1:n
            
        v_b(:,i) = R_beta'*v_st(:,i);
        a_b(:,i) = R_beta'*a_st(:,i);
        w_b(:,i) = R_beta'*w_st(:,i);
        w_dot_b(:,i) = R_beta'*w_dot_st(:,i);
             
    end
    
    
    Inertia_f = {};
            
    for i = 1:n

                
        M_matrix = Mass_matrix( m_fly, m_vw , R_L(:,:,i), R_R(:,:,i), arm_JL, arm_JR, arm_cg_L(:,i), arm_cg_R(:,i), I_body, I_v_wing, v_b(:,i), w_b(:,i), w_bL(:,i), w_bR(:,i));
                
        M11 = M_matrix.M11; %[kg]
        M12 = M_matrix.M12; %[kg*mm]
        M_dot_12 = M_matrix.M_dot_12; %[kg*mm/s]
        M13 = M_matrix.M13; %[kg*mm]
        M_dot_13 = M_matrix.M_dot_13; %[kg*mm/s]
        M14 = M_matrix.M14; %[kg*mm]
        M_dot_14 = M_matrix.M_dot_14; %[kg*mm/s]
                
        M21 = M_matrix.M21; %[kg*mm]
        M22 = M_matrix.M22; %[kg*mm^2]
        M_dot_21 = M_matrix.M_dot_21; %[kg*mm/s]
        M_dot_22 = M_matrix.M_dot_22; %[kg*mm^2/s]
        M_23 = M_matrix.M23; %[kg*mm^2]
        M_dot_23 = M_matrix.M_dot_23; %[kg*mm^2/s]
        M_24 = M_matrix.M24; %[kg*mm^2]
        M_dot_24 = M_matrix.M_dot_24; %[kg*mm^2/s]
                
        w_b_dL_dv_b = M_matrix.wb_dL_dv_b; %[kg*mm/s^2]
        v_b_dL_dv_b = M_matrix.v_b_dL_dv_b; %[kg*mm^2/s^2]
        w_b_dL_dw_b = M_matrix.w_b_dL_dw_b; %[kg*mm^2/s^2]
                
        F_acc_b(:,i) = [M11 M12]*[a_b(:,i); w_dot_b(:,i)]; %[kg*mm/s^2]
       
        F_vel_b(:,i) = [w_b_dL_dv_b+M_dot_12*w_b(:,i)+M_dot_13*w_bL(:,i)+M13*w_dot_bL(:,i)+M_dot_14*w_bR(:,i)+M14*w_dot_bR(:,i)]; %[kg*mm/s^2]
                               
        M_acc_b(:,i) = [M21 M22]*[a_b(:,i); w_dot_b(:,i)]; %[kg*mm^2/s^2]
                
        M_vel_b(:,i) = [v_b_dL_dv_b+w_b_dL_dw_b+M_dot_21*v_b(:,i)+M_dot_22*w_b(:,i)+M_dot_23*w_bL(:,i)+M_23*w_dot_bL(:,i)+M_dot_24*w_bR(:,i)+M_24*w_dot_bR(:,i)]; %[kg*mm^2/s^2]
        
        F_acc_st(:,i) = R_beta*F_acc_b(:,i);
        
        F_vel_st(:,i) = R_beta*F_vel_b(:,i);
        
        M_acc_st(:,i) = R_beta*M_acc_b(:,i);
        
        M_vel_st(:,i) = R_beta*M_vel_b(:,i);
                
%         F_I_b(:,i) = F_acc_b(:,i)+F_vel_b(:,i);
%                  
%         M_I_b(:,i) = M_acc_b(:,i)+M_vel_b(:,i);

        F_I_b(:,i) = F_acc_b(:,i);
                 
        M_I_b(:,i) = M_acc_b(:,i);
                
%         F_I_st(:,i) = R_beta*[F_acc_b(:,i)+F_vel_b(:,i)];
%                 
%         M_I_st(:,i) = R_beta*[M_acc_b(:,i)+M_vel_b(:,i)];

        F_I_st(:,i) = R_beta*[F_acc_b(:,i)];
                
        M_I_st(:,i) = R_beta*[M_acc_b(:,i)];
        
        Lin_momentum_b(:,i) = m_fly*v_b(:,i)+m_vw*(cross(w_b(:,i),arm_JL)+cross(w_b(:,i),arm_JR)+cross(w_bL(:,i),arm_cg_L(:,i))+cross(w_bR(:,i),arm_cg_R(:,i)));
        
        Ang_momentum_b(:,i) = I_body*w_b(:,i)+I_v_wing*(w_b(:,i)+w_bL(:,i))+I_v_wing*(w_b(:,i)+w_bR(:,i))+cross(arm_JL,(v_b(:,i)+cross(w_b(:,i),arm_JL)+cross(w_bL(:,i),arm_cg_L(:,i))))*m_vw+cross(arm_JR,(v_b(:,i)+cross(w_b(:,i),arm_JR)+cross(w_bR(:,i),arm_cg_R(:,i))))*m_vw;
        
        Lin_momentum_st(:,i) = R_beta*Lin_momentum_b(:,i);
        
        Ang_momentum_st(:,i) = R_beta*Ang_momentum_b(:,i);
                
    end
    
    Inertia_f.F_acc_b = F_acc_b;
    Inertia_f.F_vel_b = F_vel_b;
    Inertia_f.M_acc_b = M_acc_b;
    Inertia_f.M_vel_b = M_vel_b;
    Inertia_f.F_acc_st = F_acc_st;
    Inertia_f.F_vel_st = F_vel_st;
    Inertia_f.M_acc_st = M_acc_st;
    Inertia_f.M_vel_st = M_vel_st;
    Inertia_f.F_I_b = F_I_b;
    Inertia_f.M_I_b = M_I_b;
    Inertia_f.F_I_st = F_I_st;
    Inertia_f.M_I_st = M_I_st;
    Inertia_f.Lin_momentum_b = Lin_momentum_b;
    Inertia_f.Ang_momentum_b = Ang_momentum_b;
    Inertia_f.Lin_momentum_st = Lin_momentum_st;
    Inertia_f.Ang_momentum_st = Ang_momentum_st;
   
end

