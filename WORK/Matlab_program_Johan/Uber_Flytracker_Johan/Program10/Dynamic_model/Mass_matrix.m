function [ M_matrix ] = Mass_matrix( m_fly, m_vw , R_L, R_R, arm_JL, arm_JR, arm_cg_L, arm_cg_R, I_body, I_vw, v_b, w_b, w_bL, w_bR)


    % Return the kinetic energy matrix and some of the derivatives of the
    % kinetic energy matrix.
    
    
    c_JL = [ 0 -arm_JL(3) arm_JL(2); ...
             arm_JL(3) 0 -arm_JL(1); ...
             -arm_JL(2) arm_JL(1) 0];
    
    c_JR = [ 0 -arm_JR(3) arm_JR(2); ...
             arm_JR(3) 0 -arm_JR(1); ...
             -arm_JR(2) arm_JR(1) 0];
    
    c_cg_L = [ 0 -arm_cg_L(3) arm_cg_L(2); ...
               arm_cg_L(3) 0 -arm_cg_L(1); ...
               -arm_cg_L(2) arm_cg_L(1) 0]; 
    
    c_cg_R = [ 0 -arm_cg_R(3) arm_cg_R(2); ...
               arm_cg_R(3) 0 -arm_cg_R(1); ...
               -arm_cg_R(2) arm_cg_R(1) 0]; 
    
           
    R_cg_L = R_L'*c_cg_L*R_L; %[mm]
    
    R_cg_R = R_R'*c_cg_R*R_R; %[mm]
    
    I_vw_L = R_L'*I_vw*R_L; %[kg*mm^2]
    
    I_vw_R = R_R'*I_vw*R_R; %[kg*mm^2]
    
         
    M11 = eye(3)*m_fly; %[kg]
    
    M12 = m_vw*(-c_JL-c_JR-R_cg_L-R_cg_R); %[kg*mm]
    
    M13 = m_vw*(-R_cg_L); %[kg*mm]
    
    M14 = m_vw*(-R_cg_R); %[kg*mm]
    
    M21 = M12'; %[kg*mm]
    
    M22 = I_body+I_vw_L+I_vw_R+m_vw*(-c_JL-R_cg_L)'*(-c_JL-R_cg_L)+m_vw*(-c_JR-R_cg_R)'*(-c_JR-R_cg_R); %[kg*mm^2]
    
    M23 = I_vw_L+m_vw*(-c_JL-R_cg_L)'*(-R_cg_L); %[kg*mm^2]
    
    M24 = I_vw_R+m_vw*(-c_JR-R_cg_R)'*(-R_cg_R); %[kg*mm^2]
    
    M31 = M13'; %[kg*mm]
    
    M32 = M23'; %[kg*mm^2]
    
    M33 = I_vw_L+m_vw*(-R_cg_L)'*(-R_cg_L); %[kg*mm^2]
    
    M34 = zeros(3); %[kg*mm^2]
    
    M41 = M14'; %[kg*mm]
    
    M42 = M24'; %[kg*mm^2]
    
    M43 = M34'; %[kg*mm^2]
    
    M44 = I_vw_R+m_vw*(-R_cg_R)'*(-R_cg_R); %[kg*mm^2]
    
    % Construct the Mass matrix M:
    
    M = [ M11 M12 M13 M14; ...
          M21 M22 M23 M24; ...
          M31 M32 M33 M34; ...
          M41 M42 M43 M44];

      
    % Determine the derivatives of M:
    
    v_b_cross = [0 -v_b(3) v_b(2); ...
                 v_b(3) 0 -v_b(1); ...
                 -v_b(2) v_b(1) 0]; %[mm/s]
    
    w_b_cross = [0 -w_b(3) w_b(2); ...
                 w_b(3) 0 -w_b(1); ...
                 -w_b(2) w_b(1) 0]; %[rad/s]
    
    w_bL_cross = [0 -w_bL(3) w_bL(2); ...
                 w_bL(3) 0 -w_bL(1); ...
                 -w_bL(2) w_bL(1) 0]; %[rad/s]
    
    w_bR_cross = [0 -w_bR(3) w_bR(2); ...
                  w_bR(3) 0 -w_bR(1); ...
                  -w_bR(2) w_bR(1) 0]; %[rad/s]
    
    M_dot_12 = w_b_cross*M12; %[kg*mm/s]
    
    M_dot_13 = w_bL_cross*M13; %[kg*mm/s]
    
    M_dot_14 = w_bR_cross*M14; %[kg*mm/s]
    
    M_dot_21 = w_b_cross*M21; %[kg*mm/s]
    
    M_dot_22 = w_b_cross*M22; %[kg*mm^2/s]
    
    M_dot_23 = w_bL_cross*M23; %[kg*mm^2/s]
    
    M_dot_24 = w_bR_cross*M24; %[kg*mm^2/s]
    
    
    % Construct the Lagrangian derivative terms: dL/dvb, dl/dwb:
    
%     dL_dv_b = 0.5*(v_b'*M11+w_b'*M21+w_bL'*M31+w_bR'*M41)'+0.5*(M11*v_b+M12*w_b+M13*w_bL+M14*w_bR); %[kg*mm/s]
%     
%     dL_dw_b = 0.5*(v_b'*M12+w_b'*M22+w_bL'*M32+w_bR'*M42)'+0.5*(M12*v_b+M22*w_b+M32*w_bL+M42*w_bR); %[kg*mm^2/s]

    dL_dv_b = M11*v_b+M12*w_b+M13*w_bL+M14*w_bR; %[kg*mm/s]
    
    dL_dw_b = M12*v_b+M22*w_b+M32*w_bL+M42*w_bR; %[kg*mm^2/s]

    wb_dL_dv_b = w_b_cross*dL_dv_b; %[kg*mm/s^2]
    
    v_b_dL_dv_b = v_b_cross*dL_dv_b; %[kg*mm^2/s^2]
    
    w_b_dL_dw_b = w_b_cross*dL_dw_b; %[kg*mm^2/s^2]
    
    
    % Return structure M_matrix:
    
    M_matrix = {};
    
    M_matrix.M11 = M11;
    
    M_matrix.M12 = M12;
    
    M_matrix.M13 = M13;
    
    M_matrix.M14 = M14;
    
    M_matrix.M21 = M21;
    
    M_matrix.M22 = M22;
    
    M_matrix.M23 = M23;
    
    M_matrix.M24 = M24;
    
    M_matrix.M31 = M31;
    
    M_matrix.M32 = M32;
    
    M_matrix.M33 = M33;
    
    M_matrix.M34 = M34;
    
    M_matrix.M41 = M14;
    
    M_matrix.M42 = M42;
    
    M_matrix.M43 = M43;
    
    M_matrix.M44 = M44;
    
    
    M_matrix.M = M;
    
    
    M_matrix.M_dot_12 = M_dot_12;
    
    M_matrix.M_dot_13 = M_dot_13;
    
    M_matrix.M_dot_14 = M_dot_14;
    
    M_matrix.M_dot_21 = M_dot_21;
    
    M_matrix.M_dot_22 = M_dot_22;
    
    M_matrix.M_dot_23 = M_dot_23;
    
    M_matrix.M_dot_24 = M_dot_24;
    
    
    M_matrix.dL_dv_b = dL_dv_b;
    
    M_matrix.dL_dw_b = dL_dw_b;
    
    M_matrix.wb_dL_dv_b = wb_dL_dv_b;
    
    M_matrix.v_b_dL_dv_b = v_b_dL_dv_b;
    
    M_matrix.w_b_dL_dw_b = w_b_dL_dw_b;

end

