function [ State_update, Force_param ] = State_space(Pattern,Body,State,t,dt)


    % State-space system of a fruit fly:

    
    % Extract body parameters from the structure body:
    
    m_wl = Body.m_wl;
    
    m_w = Body.m_w;
    
    m_v_w = Body.m_v_w;
    
    wing_l = Body.wing_l;
    
    joint_L = Body.joint_L;
    
    joint_R = Body.joint_R;
    
    cg_body = Body.cg_body;
    
    cg_L = Body.cg_L;
    
    cg_R = Body.cg_R;
    
    I_body = Body.I_body; % [ kg * mm^2 ]
    
    I_wing = Body.I_wing; % [ kg * mm^2 ]
    
    I_v_wing = Body.I_v_wing; % [ kg * mm^2 ]
    
    y_sect = Body.y_sect;
    
    chords = Body.chords;
    
    
    % Get current body translation and rotation parameters:
    
    u_b = [State(4); State(5); State(6)];
    
    a_b = [State(7); State(8); State(9)];
    
    q_b = [State(10); State(11); State(12); State(13)];
    
    w_b = [State(14); State(15); State(16)];
    
    w_dot_b = [State(17); State(18); State(19)];
    
    
    % Get current wingbeat kinematics:
    
    [q_L,w_L,w_dot_L,q_R,w_R,w_dot_R,up,down] = Wing_Pattern_Generator(Pattern,t,dt);
    
        
    % Rotation matrices:
        
    R_b = quat2matNEW(q_b);
        
    R_L = quat2matNEW(q_L);
        
    R_R = quat2matNEW(q_R);

        
    % Wing angular velocity in body frame:
    
    w_bL = w_b + R_L'*w_L;
        
    w_bR = w_b + R_R'*w_R;
        
        
    % Wing angular acceleration in body frame:
        
    w_dot_bL = w_dot_b + R_L'*w_dot_L;
        
    w_dot_bR = w_dot_b + R_R'*w_dot_R;
    
    
    % State vector:
        
    V = [u_b; w_b; w_bL; w_bR];
    
    
    
    c_bjL = joint_L - cg_body;
    
    c_bjR = joint_R - cg_body;
    
    c_jwL = cg_L;
    
    c_jwR = cg_R;  
    
    
    % Cross matrices:
        
    C_cross_bjL = [0 -c_bjL(3) c_bjL(2); ...
                   c_bjL(3) 0 -c_bjL(1); ...
                   -c_bjL(1) c_bjL(2) 0];
    
    C_cross_bjR = [0 -c_bjR(3) c_bjR(2); ...
                   c_bjR(3) 0 -c_bjR(1); ...
                   -c_bjR(1) c_bjR(2) 0];
    
    C_cross_jwL = [0 -c_jwL(3) c_jwL(2); ...
                   c_jwL(3) 0 -c_jwL(1); ...
                   -c_jwL(1) c_jwL(2) 0];
    
    C_cross_jwR = [0 -c_jwR(3) c_jwR(2); ...
                   c_jwR(3) 0 -c_jwR(1); ...
                   -c_jwR(1) c_jwR(2) 0];
            
    
    % Spatial inertia matrix:

    M11 = m_wl*eye(3);  % [kg]

    M12 = m_v_w*(-C_cross_bjL - (R_L'*C_cross_jwL*R_L) - C_cross_bjR - (R_R'*C_cross_jwR*R_R));  % [kg*mm]

    M13 = m_v_w*(-R_L'*C_cross_jwL*R_L); % [kg*mm]

    M14 = m_v_w*(-R_R'*C_cross_jwR*R_R); % [kg*mm]

    M21 = M12'; % [kg*mm]

    M22 = I_body + (R_L'*I_v_wing*R_L) + (R_R'*I_v_wing*R_R) + m_v_w*(-C_cross_bjL - (R_L'*C_cross_jwL*R_L))'*(-C_cross_bjL - (R_L'*C_cross_jwL*R_L)) + ...
          m_v_w*(-C_cross_bjR - (R_R'*C_cross_jwR*R_R))'*(-C_cross_bjR - (R_R'*C_cross_jwR*R_R)); % [kg*mm^2]

    M23 = (R_L'*I_v_wing*R_L) + m_v_w*(-C_cross_bjL - (R_L'*C_cross_jwL*R_L))'*(-(R_L'*C_cross_jwL*R_L)); % [kg*mm^2]

    M24 = (R_R'*I_v_wing*R_R) + m_v_w*(-C_cross_bjR - (R_R'*C_cross_jwR*R_R))'*(-(R_R'*C_cross_jwR*R_R)); % [kg*mm^2]

    M31 = M13'; % [kg*mm]

    M32 = M23'; % [kg*mm^2]

    M33 = (R_L'*I_v_wing*R_L) + m_v_w*(R_L'*C_cross_jwL*R_L)'*(R_L'*C_cross_jwL*R_L); % [kg*mm^2]

    M34 = zeros(3); % [kg*mm^2]

    M41 = M14'; % [kg*mm^2]

    M42 = M24'; % [kg*mm^2]

    M43 = zeros(3); % [kg*mm^2]

    M44 = R_R'*I_v_wing*R_R + m_v_w*(-R_R'*C_cross_jwR*R_R)'*(-R_R'*C_cross_jwR*R_R); % [kg*mm^2]
        
        
    M = [ M11 M12 M13 M14; ...
          M21 M22 M23 M24; ...
          M31 M32 M33 M34; ...
          M41 M42 M43 M44];
          
        % [ kg    kg*mm   kg*mm   kg*mm   ]
        % [ kg*mm kg*mm^2 kg*mm^2 kg*mm^2 ]
        % [ kg*mm kg*mm^2 kg*mm^2 kg*mm^2 ] 
        % [ kg*mm kg*mm^2 kg*mm^2 kg*mm^2 ]
        
        
        
    % Kinetic Energy:
        
    KE = 0.5*V'*M*V;
        
    KE_trans = 0.5*V(1:3)'*[M11 M12 M13 M14]*V;
        
    KE_rot = 0.5*V(4:12)'*[M21 M22 M23 M24; M31 M32 M33 M34; M41 M42 M43 M44]*V;    
    
    
    % Derivatives rotation matrices:

    w_L_cross = [0 -w_L(3) w_L(2); ...
                 w_L(3) 0 -w_L(1); ...
                 -w_L(2) w_L(1) 0];
        
    w_R_cross = [0 -w_R(3) w_R(2); ...
                 w_R(3) 0 -w_R(1); ...
                 -w_R(2) w_R(1) 0];
        
%     w_bL_cross = [0 -w_bL(3) w_bL(2); ...
%                   w_bL(3) 0 -w_bL(1); ...
%                   -w_bL(2) w_bL(1) 0];
%          
%     w_bR_cross = [0 -w_bR(3) w_bR(2); ...
%                   w_bR(3) 0 -w_bR(1); ...
%                   -w_bR(2) w_bR(1) 0];         
                 
    R_L_dot = w_L_cross*R_L;

    R_R_dot = w_R_cross*R_R;
    
    
    
    % Derivatives spatial inertia submatrices:

    M12_dot = m_v_w*(-R_L_dot'*C_cross_jwL*R_L) + m_v_w*(-R_L'*C_cross_jwL*R_L_dot) + m_v_w*(-R_R_dot'*C_cross_jwR*R_R) + m_v_w*(-R_R'*C_cross_jwR*R_R_dot);  % [kg*mm/s]
        
    M13_dot = m_v_w*(-R_L_dot'*C_cross_jwL*R_L) + m_v_w*(-R_L'*C_cross_jwL*R_L_dot); % [kg*mm/s]
        
    M14_dot = m_v_w*(-R_R_dot'*C_cross_jwR*R_R) + m_v_w*(-R_R'*C_cross_jwR*R_R_dot); % [kg*mm/s]
        
    M21_dot = M12_dot'; % [kg*mm/s]
        
    M22_dot = R_L_dot'*I_wing*R_L + R_L'*I_wing*R_L_dot + R_R_dot'*I_wing*R_R + R_R'*I_wing*R_R_dot + m_v_w*(-R_L_dot'*C_cross_jwL*R_L)'*(-C_cross_bjL - R_L'*C_cross_jwL*R_L) ...
                   + m_v_w*(-R_L'*C_cross_jwL*R_L_dot)'*(-C_cross_bjL - R_L'*C_cross_jwL*R_L)  + m_v_w*(-R_L'*C_cross_jwL*R_L)'*(-R_L_dot'*C_cross_jwL*R_L) + ...
                   m_v_w*(-R_L'*C_cross_jwL*R_L)'*(-R_L'*C_cross_jwL*R_L_dot) + m_v_w*(-R_R_dot'*C_cross_jwR*R_R)'*(-C_cross_bjR - R_R'*C_cross_jwR*R_R) ...
                   + m_v_w*(-R_R'*C_cross_jwR*R_R_dot)'*(-C_cross_bjR - R_R'*C_cross_jwR*R_R)  + m_v_w*(-R_R'*C_cross_jwR*R_R)'*(-R_R_dot'*C_cross_jwR*R_R) + ...
                   m_v_w*(-R_R'*C_cross_jwR*R_R)'*(-R_R'*C_cross_jwR*R_R_dot); % [kg*mm^2/s]
       
    M23_dot = R_L_dot'*I_wing*R_L + R_L'*I_wing*R_L_dot + m_v_w*(-R_L_dot'*C_cross_jwL*R_L)'*(-R_L'*C_cross_jwL*R_L) + m_v_w*(-R_L'*C_cross_jwL*R_L_dot)'*(-R_L'*C_cross_jwL*R_L) ...
                  + m_v_w*(-C_cross_bjL - R_L'*C_cross_jwL*R_L)'*(-R_L_dot'*C_cross_jwL*R_L) + m_v_w*(-C_cross_bjL - R_L'*C_cross_jwL*R_L)'*(-R_L'*C_cross_jwL*R_L_dot); % [kg*mm^2/s]
        
    M24_dot = R_R_dot'*I_wing*R_R + R_R'*I_wing*R_R_dot + m_v_w*(-R_R_dot'*C_cross_jwR*R_R)'*(-R_R'*C_cross_jwR*R_R) + m_v_w*(-R_R'*C_cross_jwR*R_R_dot)'*(-R_R'*C_cross_jwR*R_R) ...
                  + m_v_w*(-C_cross_bjR - R_R'*C_cross_jwR*R_R)'*(-R_R_dot'*C_cross_jwR*R_R) + m_v_w*(-C_cross_bjR - R_R'*C_cross_jwR*R_R)'*(-R_R'*C_cross_jwR*R_R_dot); % [kg*mm^2/s]
        
    M31_dot = M13_dot'; % [kg*mm/s]
        
    M32_dot = M23_dot'; % [kg*mm^2/s]
        
    M41_dot = M14_dot'; % [kg*mm/s]
        
    M42_dot = M24_dot'; % [kg*mm^2/s]    
    
    
    % Lagrange's equation derivatives:

    dL_dub = (M11*u_b+M12*w_b+M13*w_bL+M14*w_bR); % [kg*mm/s]
        
    dL_wub = (M21*u_b+M22*w_b+M23*w_bL+M24*w_bR); % [kg*mm^2/s]   
    
    
    % Inertia forces and moments:
        
    F1 = M11*a_b; %[kg*mm/s^2]
        
    F2 = M12*w_dot_b; %[kg*mm/s^2]
        
    F3 = M12_dot*w_b; %[kg*mm/s^2]
        
    F4 = M13_dot*w_bL; %[kg*mm/s^2]
        
    F5 = M14_dot*w_bR; %[kg*mm/s^2]
        
    F6 = M13*w_dot_bL; %[kg*mm/s^2]
        
    F7 = M14*w_dot_bR; %[kg*mm/s^2]
        
    F8 = cross(w_b,dL_dub); %[kg*mm/s^2]
        
    %F_I = 1e-3*(M11*a_b+M12*w_dot_b+M12_dot*w_b+M13_dot*w_bL+M14_dot*w_bR+M13*w_dot_bL+M14*w_dot_bR+cross(w_b,dL_dub)); %[kg*m/s^2] [N]
    
    F_I = (F1+F2+F3+F4+F5+F6+F7+F8); %[kg*m/s^2] [N]
        
    M1 = M21*a_b; %[kg*mm/s^2*mm]
        
    M2 = M22*w_dot_b; %[kg*mm/s^2*mm]
        
    M3 = cross(u_b,dL_dub); %[kg*mm/s^2*mm]
        
    M4 = cross(w_b,dL_wub); %[kg*mm/s^2*mm]
        
    M5 = M21_dot*u_b; % [kg*mm/s^2 *mm]
        
    M6 = M22_dot*w_b; % [kg*mm/s^2 *mm]
        
    M7 = M23_dot*w_bL; % [kg*mm/s^2 *mm]
        
    M8 = M24_dot*w_bR; % [kg*mm/s^2 *mm]
        
    M9 = M23*w_dot_bL; %[kg*mm/s^2 *mm]
        
    M10 = M24*w_dot_bR; %[kg*mm/s^2 *mm]
        
    %M_I = 1e-3*(M21*a_b+M22*w_dot_b+cross(u_b,dL_dub)+cross(w_b,dL_wub)+M21_dot*u_b+M22_dot*w_b+M23_dot*w_bL+M24_dot*w_bR+M23*w_dot_bL+M24*w_dot_bR); %[kg*m/s^2*mm] [N*mm]
    
    M_I = (M1+M2+M3+M4+M5+M6+M7+M8+M9+M10); %[kg*m/s^2*mm] [N*mm]
    
%     % Inertia forces and moments:
%         
%     F1 = M11*a_b; %[kg*mm/s^2]
%         
%     F2 = M12*w_dot_b; %[kg*mm/s^2]
%         
%     F3 = M12_dot*w_b; %[kg*mm/s^2]
%         
%     F4 = M13_dot*w_bL; %[kg*mm/s^2]
%         
%     F5 = M14_dot*w_bR; %[kg*mm/s^2]
%         
%     F6 = M13*w_dot_bL; %[kg*mm/s^2]
%         
%     F7 = M14*w_dot_bR; %[kg*mm/s^2]
%         
%     F8 = cross(w_b,dL_dub); %[kg*mm/s^2]
%         
%     %F_I = 1e-3*(M11*a_b+M12*w_dot_b+M12_dot*w_b+M13_dot*w_bL+M14_dot*w_bR+M13*w_dot_bL+M14*w_dot_bR+cross(w_b,dL_dub)); %[kg*m/s^2] [N]
%     
%     F_I = F1+F2+F3+F4+F5+F6+F7+F8; %[kg*m/s^2] [N]
%         
%     M1 = M21*a_b; %[kg*mm/s^2*mm]
%         
%     M2 = M22*w_dot_b; %[kg*mm/s^2*mm]
%         
%     M3 = cross(u_b,dL_dub); %[kg*mm/s^2*mm]
%         
%     M4 = cross(w_b,dL_wub); %[kg*mm/s^2*mm]
%         
%     M5 = M21_dot*u_b; % [kg*mm/s^2 *mm]
%         
%     M6 = M22_dot*w_b; % [kg*mm/s^2 *mm]
%         
%     M7 = M23_dot*w_bL; % [kg*mm/s^2 *mm]
%         
%     M8 = M24_dot*w_bR; % [kg*mm/s^2 *mm]
%         
%     M9 = M23*w_dot_bL; %[kg*mm/s^2 *mm]
%         
%     M10 = M24*w_dot_bR; %[kg*mm/s^2 *mm]
%         
%     %M_I = 1e-3*(M21*a_b+M22*w_dot_b+cross(u_b,dL_dub)+cross(w_b,dL_wub)+M21_dot*u_b+M22_dot*w_b+M23_dot*w_bL+M24_dot*w_bR+M23*w_dot_bL+M24*w_dot_bR); %[kg*m/s^2*mm] [N*mm]
%     
%     M_I = M1+M2+M3+M4+M5+M6+M7+M8+M9+M10; %[kg*m/s^2*mm] [N*mm]    
    
    % Get Aerodynamic forces and center of pressures:
      
    [ F_aero_L, F_aero_R, c_pres_L, c_pres_R ] = Aerodynamic_Model(wing_l,y_sect,chords,joint_L,joint_R,u_b,w_b,w_L,w_R,R_b,R_L,R_R,up,down);
    
        
    F_aero_L_arm = c_bjL+R_L*c_pres_L;     
        
    F_aero_R_arm = c_bjR+R_R*c_pres_R;  
        
        
    F_A = R_L*F_aero_L+R_R*F_aero_R; % [N]
        
    M_A = cross(F_aero_L_arm,R_L*F_aero_L)+cross(F_aero_R_arm,R_R*F_aero_R); % [N*mm]
    
    
    
    % Gravitational forces and moments:
        
    F_G = R_b'*[0; 0; -m_wl*9810]; %[N]
        
    M_G = [0; 0; 0]; %[N*mm]
    
    

    % Resultant forces:
        
%     F_R = F_A+F_G;
%         
%     M_R = M_A+M_G;    
    
    F_R = F_G;
        
    M_R = M_G;     
    
    % Calculate acceleration and angular acceleration
        
    M_inv = inv([M11 M12; M21 M22]);
        
%     a_b_update = M_inv(1:3,:)*[F_R-F_I; M_R-M_I];
%         
%     w_dot_b_update = M_inv(4:6,:)*[F_R-F_I; M_R-M_I];
    
    a_b_update = M_inv(1:3,:)*[F_R-(F3+F4+F5+F6+F7+F8); M_R-(M3+M4+M5+M6+M7+M8+M9+M10)];
        
    w_dot_b_update = M_inv(4:6,:)*[F_R-(F3+F4+F5+F6+F7+F8); M_R-(M3+M4+M5+M6+M7+M8+M9+M10)];

    
%     % Calculate the radius of gyration:
%     
%     I_tot = 
%     
%     r_gyration = [sqrt(
    
    
    
    
    
    % Update the state vector and save the forces and moments and other
    % parameters
    
    Force_param = {};
    
    
    Force_param.F1 = F1;
    
    Force_param.F2 = F2;
    
    Force_param.F3 = F3;
    
    Force_param.F4 = F4;
    
    Force_param.F5 = F5;
    
    Force_param.F6 = F6;
    
    Force_param.F7 = F7;
    
    Force_param.F8 = F8;
    
    Force_param.F_I = F_I;
    
    
    Force_param.M1 = M1;
    
    Force_param.M2 = M2;
    
    Force_param.M3 = M3;
    
    Force_param.M4 = M4;
    
    Force_param.M5 = M5;
    
    Force_param.M6 = M6;
    
    Force_param.M7 = M7;
    
    Force_param.M8 = M8;
    
    Force_param.M9 = M9;
    
    Force_param.M10 = M10;
    
    Force_param.M_I = M_I;   
    
    
    Force_param.F_aero_L = F_aero_L;
    
    Force_param.F_aero_R = F_aero_R;
    
    Force_param.F_aero_L_arm = F_aero_L_arm;
    
    Force_param.F_aero_R_arm = F_aero_R_arm;
    
    Force_param.F_A = F_A;
    
    Force_param.M_A = M_A;
    
    
    Force_param.F_g = F_G;
    
    Force_param.M_g = M_G;
    
    
    Force_param.F_R = F_R;
    
    Force_param.M_R = M_R;
    
    
    Force_param.M_inv = M_inv;
    
    
    Force_param.KE = KE;
    
    Force_param.KE_trans = KE_trans;
    
    Force_param.KE_rot = KE_rot;
    
    
    Force_param.R_b = R_b;
    
    Force_param.R_L = R_L;
    
    Force_param.R_R = R_R;
    
    
    
    Force_param.q_L = q_L;
    
    Force_param.q_R = q_R;
    
    Force_param.w_L = w_L;
    
    Force_param.w_R = w_R;
    
    Force_param.w_dot_L = w_dot_L;
    
    Force_param.w_dot_R = w_dot_R;
    
    
    
    State_update = zeros(length(State),1);
    
    State_update = State;
    
    State_update(7:9) = a_b_update;
    
    State_update(17:19) = w_dot_b_update;
    

    
    
end

