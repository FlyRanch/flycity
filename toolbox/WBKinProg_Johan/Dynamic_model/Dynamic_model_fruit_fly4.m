function Dynamic_model_fruit_fly4(settings, pathDB, seq_nr)

    % Program contains the dynamic model of a fruitfly
    
    
    % Time:
    
    dt = pathDB.t(2)-pathDB.t(1);
    
    t = pathDB.t;
    
    start = find(isnan(pathDB.x(:,seq_nr))==0, 1 );
    stop = find(isnan(pathDB.x(:,seq_nr))==0, 1, 'last' );
    
    
    % Body parameters:
    
    wing_l = pathDB.wing_l(seq_nr);
    
    m_wl = 1e-6*0.1078*wing_l^3.008;
    
    [cg_body, cg_L, cg_R, I_body, I_wing, I_v_wing, m_b, m_w, m_v_w ] = cg_plus_Inertia(settings, pathDB, m_wl, wing_l ,seq_nr );
    
    c_bjL = pathDB.joint_pos_L(:,seq_nr) - cg_body;
    
    c_bjR = pathDB.joint_pos_R(:,seq_nr) - cg_body;
    
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
    
               
               
    KE = zeros(stop-start+1,1);
    
    KE_trans = zeros(stop-start+1,1);
    
    KE_rot = zeros(stop-start+1,1);    
    
    R_L_dot_test = zeros(3,3,stop-start+1);
    
    R_R_dot_test = zeros(3,3,stop-start+1);
    
    
    
    
    % Test matrices inertia:
    
    F1 = zeros(3,stop-start+1);
    
    F2 = zeros(3,stop-start+1);
    
    F3 = zeros(3,stop-start+1);
    
    F4 = zeros(3,stop-start+1);
    
    F5 = zeros(3,stop-start+1);
    
    F6 = zeros(3,stop-start+1);
    
    F7 = zeros(3,stop-start+1);
    
    F8 = zeros(3,stop-start+1);
    
    F_I = zeros(3,stop-start+1);
    
    
    M1 = zeros(3,stop-start+1);
    
    M2 = zeros(3,stop-start+1);
    
    M3 = zeros(3,stop-start+1);
    
    M4 = zeros(3,stop-start+1);
    
    M5 = zeros(3,stop-start+1);
    
    M6 = zeros(3,stop-start+1);
    
    M7 = zeros(3,stop-start+1);
    
    M8 = zeros(3,stop-start+1);
    
    M9 = zeros(3,stop-start+1);
    
    M10 = zeros(3,stop-start+1);
    
    M_I = zeros(3,stop-start+1);
    
    A_b = zeros(3,stop-start+1);
    
    W_dot_b = zeros(3,stop-start+1);
    
    
    w_dot_bL_test = zeros(3,stop-start+1);
    
    w_dot_bR_test = zeros(3,stop-start+1);
    
    a_b_calc = zeros(3,stop-start+1);
    
    w_b_calc = zeros(3,stop-start+1);
    
    F_R = zeros(3,stop-start+1);
    
    M_R = zeros(3,stop-start+1);
    
    
    M_I_strk = zeros(3,stop-start+1);
    
    I_strk = zeros(3,stop-start+1);
    

    for i = start:stop
        
        
        % Position and orientation:
        
        %u_b = [pathDB.u_filt(i,seq_nr); pathDB.v_filt(i,seq_nr); pathDB.w_filt(i,seq_nr)];

        u_b = [pathDB.u_body(i,seq_nr); pathDB.v_body(i,seq_nr); pathDB.w_body(i,seq_nr)];
        
        %a_b = [pathDB.ax_filt(i,seq_nr); pathDB.ay_filt(i,seq_nr); pathDB.az_filt(i,seq_nr)];
        
        a_b = [pathDB.ax_body(i,seq_nr); pathDB.ay_body(i,seq_nr); pathDB.az_body(i,seq_nr)];

        A_b(:,i-start+1) = a_b;
        
        q_b = [pathDB.qb1_filt(i,seq_nr); pathDB.qb2_filt(i,seq_nr); pathDB.qb3_filt(i,seq_nr); pathDB.qb4_filt(i,seq_nr)];
        
        q_L = [pathDB.qL1_filt2(i,seq_nr); pathDB.qL2_filt2(i,seq_nr); pathDB.qL3_filt2(i,seq_nr); pathDB.qL4_filt2(i,seq_nr)];
        
        q_R = [pathDB.qR1_filt2(i,seq_nr); pathDB.qR2_filt2(i,seq_nr); pathDB.qR3_filt2(i,seq_nr); pathDB.qR4_filt2(i,seq_nr)];
        
        w_b = [pathDB.b_omega1(i,seq_nr); pathDB.b_omega2(i,seq_nr); pathDB.b_omega3(i,seq_nr)];
        
        w_L = [pathDB.omega1_L(i,seq_nr); pathDB.omega2_L(i,seq_nr); pathDB.omega3_L(i,seq_nr)];
        
        w_R = [pathDB.omega1_R(i,seq_nr); pathDB.omega2_R(i,seq_nr); pathDB.omega3_R(i,seq_nr)];
        
                
        % Angular acceleration:
        
        if i > start
        
        w_dot_b = (w_b-[pathDB.b_omega1(i-1,seq_nr); pathDB.b_omega2(i-1,seq_nr); pathDB.b_omega3(i-1,seq_nr)])./dt;
    
        w_dot_L = (w_L-[pathDB.omega1_L(i-1,seq_nr); pathDB.omega2_L(i-1,seq_nr); pathDB.omega3_L(i-1,seq_nr)])./dt;
    
        w_dot_R = (w_R-[pathDB.omega1_R(i-1,seq_nr); pathDB.omega2_R(i-1,seq_nr); pathDB.omega3_R(i-1,seq_nr)])./dt;   
        
        W_dot_b(:,i-start+1) = w_dot_b;
        
        else
            
        w_dot_b = [0; 0; 0];
        
        w_dot_L = [0; 0; 0];
        
        w_dot_R = [0; 0; 0];
        
        end

        
%         % Vacuum conditions:
%         
%         u_b = [0; 0; 0];
%         
%         a_b = [0; 0; 0];
% 
%         w_b = [0; 0; 0];
% %         
%         w_dot_b = [0; 0; 0];
%         
%         q_b = [0; sin((55/180)*pi/2); 0; cos((55/180)*pi/2)];
%   
%         W_dot_b(:,i-start+1) = w_dot_b;
%         
%         w_dot_L = [0; 0; 0];
%         
%         w_dot_R = [0; 0; 0];

        
               
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
        
        
        w_dot_bL_test(:,i-start+1) = w_dot_bL;
        
        w_dot_bR_test(:,i-start+1) = w_dot_bR;
        
        
        
        % State vector:
        
        V = [u_b; w_b; w_bL; w_bR];
        
        
        % Derivative state vector:
        
        V_dot = [a_b; w_dot_b; w_dot_bL; w_dot_bR];
        
        
        % Spatial inertia matrix:

        M11 = m_wl*eye(3);  % [kg]

        M12 = m_w*(-C_cross_bjL - (R_L'*C_cross_jwL*R_L) - C_cross_bjR - (R_R'*C_cross_jwR*R_R));  % [kg*mm]

        M13 = m_w*(-R_L'*C_cross_jwL*R_L); % [kg*mm]

        M14 = m_w*(-R_R'*C_cross_jwR*R_R); % [kg*mm]

        M21 = M12'; % [kg*mm]

        M22 = I_body + (R_L'*I_v_wing*R_L) + (R_R'*I_v_wing*R_R) + m_w*(-C_cross_bjL - (R_L'*C_cross_jwL*R_L))'*(-C_cross_bjL - (R_L'*C_cross_jwL*R_L)) + ...
              m_w*(-C_cross_bjR - (R_R'*C_cross_jwR*R_R))'*(-C_cross_bjR - (R_R'*C_cross_jwR*R_R)); % [kg*mm^2]

        M23 = (R_L'*I_v_wing*R_L) + m_w*(-C_cross_bjL - (R_L'*C_cross_jwL*R_L))'*(-(R_L'*C_cross_jwL*R_L)); % [kg*mm^2]

        M24 = (R_R'*I_v_wing*R_R) + m_w*(-C_cross_bjR - (R_R'*C_cross_jwR*R_R))'*(-(R_R'*C_cross_jwR*R_R)); % [kg*mm^2]

        M31 = M13'; % [kg*mm]

        M32 = M23'; % [kg*mm^2]

        M33 = (R_L'*I_v_wing*R_L) + m_w*(R_L'*C_cross_jwL*R_L)'*(R_L'*C_cross_jwL*R_L); % [kg*mm^2]

        M34 = zeros(3); % [kg*mm^2]

        M41 = M14'; % [kg*mm^2]

        M42 = M24'; % [kg*mm^2]

        M43 = zeros(3); % [kg*mm^2]

        M44 = R_R'*I_v_wing*R_R + m_w*(-R_R'*C_cross_jwR*R_R)'*(-R_R'*C_cross_jwR*R_R); % [kg*mm^2]
        
        
        M = [ M11 M12 M13 M14; ...
              M21 M22 M23 M24; ...
              M31 M32 M33 M34; ...
              M41 M42 M43 M44];
          
        % [ kg    kg*mm   kg*mm   kg*mm   ]
        % [ kg*mm kg*mm^2 kg*mm^2 kg*mm^2 ]
        % [ kg*mm kg*mm^2 kg*mm^2 kg*mm^2 ] 
        % [ kg*mm kg*mm^2 kg*mm^2 kg*mm^2 ]        
        
        
        % Kinetic Energy:
        
        KE(i-start+1) = 0.5*V'*M*V;
        
        KE_trans(i-start+1) = 0.5*V(1:3)'*[M11 M12 M13 M14]*V;
        
        KE_rot(i-start+1) = 0.5*V(4:12)'*[M21 M22 M23 M24; M31 M32 M33 M34; M41 M42 M43 M44]*V;
        
        
        
        % Derivatives rotation matrices:

        w_L_cross = [0 -w_L(3) w_L(2); ...
                     w_L(3) 0 -w_L(1); ...
                     -w_L(2) w_L(1) 0];
        
        w_R_cross = [0 -w_R(3) w_R(2); ...
                     w_R(3) 0 -w_R(1); ...
                     -w_R(2) w_R(1) 0];
        
%         w_bL_cross = [0 -w_bL(3) w_bL(2); ...
%                      w_bL(3) 0 -w_bL(1); ...
%                      -w_bL(2) w_bL(1) 0];
%         
%         w_bR_cross = [0 -w_bR(3) w_bR(2); ...
%                      w_bR(3) 0 -w_bR(1); ...
%                      -w_bR(2) w_bR(1) 0];         
                 
%         R_L_dot = w_L_cross*R_L;
% 
%         R_R_dot = w_R_cross*R_R;
        
        R_L_dot = w_L_cross*R_L;

        R_R_dot = w_R_cross*R_R;
        
        
        
        % Derivatives spatial inertia submatrices:

        M12_dot = m_w*(-R_L_dot'*C_cross_jwL*R_L) + m_w*(-R_L'*C_cross_jwL*R_L_dot) + m_w*(-R_R_dot'*C_cross_jwR*R_R) + m_w*(-R_R'*C_cross_jwR*R_R_dot);  % [kg*mm]
        
        M13_dot = m_w*(-R_L_dot'*C_cross_jwL*R_L) + m_w*(-R_L'*C_cross_jwL*R_L_dot); % [kg*mm]
        
        M14_dot = m_w*(-R_R_dot'*C_cross_jwR*R_R) + m_w*(-R_R'*C_cross_jwR*R_R_dot); % [kg*mm]
        
        M21_dot = M12_dot'; % [kg*mm]
        
        M22_dot = R_L_dot'*I_wing*R_L + R_L'*I_wing*R_L_dot + R_R_dot'*I_wing*R_R + R_R'*I_wing*R_R_dot + m_w*(-R_L_dot'*C_cross_jwL*R_L)'*(-C_cross_bjL - R_L'*C_cross_jwL*R_L) ...
                   + m_w*(-R_L'*C_cross_jwL*R_L_dot)'*(-C_cross_bjL - R_L'*C_cross_jwL*R_L)  + m_w*(-R_L'*C_cross_jwL*R_L)'*(-R_L_dot'*C_cross_jwL*R_L) + ...
                   m_w*(-R_L'*C_cross_jwL*R_L)'*(-R_L'*C_cross_jwL*R_L_dot) + m_w*(-R_R_dot'*C_cross_jwR*R_R)'*(-C_cross_bjR - R_R'*C_cross_jwR*R_R) ...
                   + m_w*(-R_R'*C_cross_jwR*R_R_dot)'*(-C_cross_bjR - R_R'*C_cross_jwR*R_R)  + m_w*(-R_R'*C_cross_jwR*R_R)'*(-R_R_dot'*C_cross_jwR*R_R) + ...
                   m_w*(-R_R'*C_cross_jwR*R_R)'*(-R_R'*C_cross_jwR*R_R_dot); % [kg*mm^2]
       
        M23_dot = R_L_dot'*I_wing*R_L + R_L'*I_wing*R_L_dot + m_w*(-R_L_dot'*C_cross_jwL*R_L)'*(-R_L'*C_cross_jwL*R_L) + m_w*(-R_L'*C_cross_jwL*R_L_dot)'*(-R_L'*C_cross_jwL*R_L) ...
                  + m_w*(-C_cross_bjL - R_L'*C_cross_jwL*R_L)'*(-R_L_dot'*C_cross_jwL*R_L) + m_w*(-C_cross_bjL - R_L'*C_cross_jwL*R_L)'*(-R_L'*C_cross_jwL*R_L_dot); % [kg*mm^2]
        
        M24_dot = R_R_dot'*I_wing*R_R + R_R'*I_wing*R_R_dot + m_w*(-R_R_dot'*C_cross_jwR*R_R)'*(-R_R'*C_cross_jwR*R_R) + m_w*(-R_R'*C_cross_jwR*R_R_dot)'*(-R_R'*C_cross_jwR*R_R) ...
                  + m_w*(-C_cross_bjR - R_R'*C_cross_jwR*R_R)'*(-R_R_dot'*C_cross_jwR*R_R) + m_w*(-C_cross_bjR - R_R'*C_cross_jwR*R_R)'*(-R_R'*C_cross_jwR*R_R_dot); % [kg*mm^2]
        
        M31_dot = M13_dot'; % [kg*mm]
        
        M32_dot = M23_dot'; % [kg*mm^2]
        
        M41_dot = M14_dot'; % [kg*mm]
        
        M42_dot = M24_dot'; % [kg*mm^2]
        
        
        % Lagrange's equation derivatives:

        dL_dub = (M11*u_b+M12*w_b+M13*w_bL+M14*w_bR); % [kg*mm/s]
        
        dL_wub = (M21*u_b+M22*w_b+M23*w_bL+M24*w_bR); % [kg*mm^2/s]
%         
%         d_dt_dL_dub = 
%         
%         d_dt_dL_wub = 


        % Inertia forces and moments:
        
        F1(:,i-start+1) = M11*a_b; %[kg*mm/s^2]
        
        F2(:,i-start+1) = M12*w_dot_b; %[kg*mm/s^2]
        
        F3(:,i-start+1) = M12_dot*w_b; %[kg*mm/s^2]
        
        F4(:,i-start+1) = M13_dot*w_bL; %[kg*mm/s^2]
        
        F5(:,i-start+1) = M14_dot*w_bR; %[kg*mm/s^2]
        
        F6(:,i-start+1) = M13*w_dot_bL; %[kg*mm/s^2]
        
        F7(:,i-start+1) = M14*w_dot_bR; %[kg*mm/s^2]
        
        F8(:,i-start+1) = cross(w_b,dL_dub); %[kg*mm/s^2]
        
        F_I(:,i-start+1) = 1e-3*(M11*a_b+M12*w_dot_b+M12_dot*w_b+M13_dot*w_bL+M14_dot*w_bR+M13*w_dot_bL+M14*w_dot_bR+cross(w_b,dL_dub)); %[kg*m/s^2] [N]
        
        M1(:,i-start+1) = M21*a_b; %[kg*mm/s^2*mm]
        
        M2(:,i-start+1) = M22*w_dot_b; %[kg*mm/s^2*mm]
        
        M3(:,i-start+1) = cross(u_b,dL_dub); %[kg*mm/s^2*mm]
        
        M4(:,i-start+1) = cross(w_b,dL_wub); %[kg*mm/s^2*mm]
        
        M5(:,i-start+1) = M21_dot*u_b; % [kg*mm/s^2 *mm]
        
        M6(:,i-start+1) = M22_dot*w_b; % [kg*mm/s^2 *mm]
        
        M7(:,i-start+1) = M23_dot*w_bL; % [kg*mm/s^2 *mm]
        
        M8(:,i-start+1) = M24_dot*w_bR; % [kg*mm/s^2 *mm]
        
        M9(:,i-start+1) = M23*w_dot_bL; %[kg*mm/s^2 *mm]
        
        M10(:,i-start+1) = M24*w_dot_bR; %[kg*mm/s^2 *mm]
        
        M_I(:,i-start+1) = 1e-3*(M21*a_b+M22*w_dot_b+cross(u_b,dL_dub)+cross(w_b,dL_wub)+M21_dot*u_b+M22_dot*w_b+M23_dot*w_bL+M24_dot*w_bR+M23*w_dot_bL+M24*w_dot_bR); %[kg*m/s^2*mm] [N*mm]

        
        % Aerodynamic forces and moments
        

        F_aero_L = pathDB.F_L(i,:,seq_nr)';
        
        F_aero_L_arm = c_bjL+R_L*[0; pathDB.F_L_arm(i,seq_nr); 0];      
        
        F_aero_R = pathDB.F_R(i,:,seq_nr)';
        
        F_aero_R_arm = c_bjR+R_R*[0; pathDB.F_R_arm(i,seq_nr); 0];  
        
        
        F_A = R_L*F_aero_L+R_R*F_aero_R; % [N]
        
        M_A = cross(F_aero_L_arm,R_L*F_aero_L)+cross(F_aero_R_arm,R_R*F_aero_R); % [N*mm]
        
        
        % Gravitational forces and moments
        
        F_G = R_b'*[0; 0; -m_wl*9.81];
        
        M_G = [0; 0; 0];
        
        
        % Calculate acceleration and angular acceleration
        
        M_inv = inv([M11 M12; M21 M22]);
        
        a_b_calc(:,i-start+1) = M_inv(1:3,:)*[F_A+F_G-(F6(:,i-start+1)+F7(:,i-start+1)+F8(:,i-start+1)); M_A+M_G-(M9(:,i-start+1)+M10(:,i-start+1))]; 
        
        w_b_calc(:,i-start+1) = M_inv(4:6,:)*[F_A+F_G-(F6(:,i-start+1)+F7(:,i-start+1)+F8(:,i-start+1)); M_A+M_G-(M9(:,i-start+1)+M10(:,i-start+1))]; 
        

        % Resultant forces:
        
        F_R(:,i-start+1) = F_A+F_G;
        
        M_R(:,i-start+1) = M_A+M_G;
        
        
%         % Calculate the inertia matrix in the strokeplane reference frame.
%         % The strokeplane reference frame axes are regarded as the
%         % principal axes:
% 
%         beta = -(55/180)*pi;
% 
%         R_strk = [cos(beta) 0 -sin(beta); ...
%                   0 1 0; ...
%                   sin(beta) 0 cos(beta)];        
%         
%         M_I_strk(:,i-start+1) = R_strk*M_I(:,i-start+1);
%         
%         w_b_strk = R_strk*w_b;
%         
%         w_dot_b_strk = R_strk*w_dot_b;
%         
%         W_strk = [ w_dot_b_strk(1) -w_b_strk(2)*w_b_strk(3) w_b_strk(2)*w_b_strk(3); ...
%                    w_b_strk(1)*w_b_strk(3) w_dot_b_strk(2) -w_b_strk(1)*w_b_strk(3); ...
%                    -w_b_strk(1)*w_b_strk(2) w_b_strk(1)*w_b_strk(2) w_dot_b_strk(3)];
%         
%         I_strk(:,i-start+1) = W_strk\M_I_strk(:,i-start+1);
        
        
    end
    
    
    
    
    figure()
    plot(KE,'b')
    hold on
    plot(KE_trans,'r')
    plot(KE_rot,'g')
    hold off
    
    figure()
    plot(t(start:stop),1000*sqrt(F_I(1,:).^2+F_I(2,:).^2+F_I(3,:).^2), ... 
    t(start:stop),sqrt(F1(1,:).^2+F1(2,:).^2+F1(3,:).^2), ...
    t(start:stop),sqrt(F2(1,:).^2+F2(2,:).^2+F2(3,:).^2), ...
    t(start:stop),sqrt(F3(1,:).^2+F3(2,:).^2+F3(3,:).^2), ...
    t(start:stop),sqrt(F4(1,:).^2+F4(2,:).^2+F4(3,:).^2), ...
    t(start:stop),sqrt(F5(1,:).^2+F5(2,:).^2+F5(3,:).^2), ...
    t(start:stop),sqrt(F6(1,:).^2+F6(2,:).^2+F6(3,:).^2), ...
    t(start:stop),sqrt(F7(1,:).^2+F7(2,:).^2+F7(3,:).^2), ...
    t(start:stop),sqrt(F8(1,:).^2+F8(2,:).^2+F8(3,:).^2))
    legend('F_I','F1','F2','F3','F4','F5','F6','F7','F8')
    
    
    
    figure()
    plot(t(start:stop),1000*sqrt(M_I(1,:).^2+M_I(2,:).^2+M_I(3,:).^2), ...
    t(start:stop),sqrt(M1(1,:).^2+M1(2,:).^2+M1(3,:).^2), ...
    t(start:stop),sqrt(M2(1,:).^2+M2(2,:).^2+M2(3,:).^2), ...
    t(start:stop),sqrt(M3(1,:).^2+M3(2,:).^2+M3(3,:).^2), ...
    t(start:stop),sqrt(M4(1,:).^2+M4(2,:).^2+M4(3,:).^2), ...
    t(start:stop),sqrt(M5(1,:).^2+M5(2,:).^2+M5(3,:).^2), ...
    t(start:stop),sqrt(M6(1,:).^2+M6(2,:).^2+M6(3,:).^2), ...
    t(start:stop),sqrt(M7(1,:).^2+M7(2,:).^2+M7(3,:).^2), ...
    t(start:stop),sqrt(M8(1,:).^2+M8(2,:).^2+M8(3,:).^2), ...
    t(start:stop),sqrt(M9(1,:).^2+M9(2,:).^2+M9(3,:).^2), ...
    t(start:stop),sqrt(M10(1,:).^2+M10(2,:).^2+M10(3,:).^2))
    legend('M_I','M1','M2','M3','M4','M5','M6','M7','M8','M9','M10')

    F_I2 = 1e-3*(F1+F6+F7+F8);
    
    M_I2 = 1e-3*(M9+M10);
    
    figure()
    plot(t(start:stop),sqrt(F_I(1,:).^2+F_I(2,:).^2+F_I(3,:).^2), ...
         t(start:stop),sqrt(F_I2(1,:).^2+F_I2(2,:).^2+F_I2(3,:).^2))
    legend('F_I','F_I2')
    
    figure()
    plot(t(start:stop),sqrt(M_I(1,:).^2+M_I(2,:).^2+M_I(3,:).^2), ...
         t(start:stop),sqrt(M_I2(1,:).^2+M_I2(2,:).^2+M_I2(3,:).^2))
    legend('M_I','M_I2')
    
    figure()
    subplot(3,1,1); plot(w_dot_bL_test(1,:))
    title('w dot bL')
    subplot(3,1,2); plot(w_dot_bL_test(2,:))
    subplot(3,1,3); plot(w_dot_bL_test(3,:))
    
    figure()
    subplot(3,1,1); plot(w_dot_bR_test(1,:))
    title('w dot bR')
    subplot(3,1,2); plot(w_dot_bR_test(2,:))
    subplot(3,1,3); plot(w_dot_bR_test(3,:))
         
    
    
    figure()
    subplot(3,1,1); plot(M_I2(1,:))
    subplot(3,1,2); plot(M_I2(2,:))
    subplot(3,1,3); plot(M_I2(3,:))
    
    
    figure()
    subplot(3,1,1); plot(t(start:stop),A_b(1,:),t(start:stop),a_b_calc(1,:))
    subplot(3,1,2); plot(t(start:stop),A_b(2,:),t(start:stop),a_b_calc(2,:))
    subplot(3,1,3); plot(t(start:stop),A_b(3,:),t(start:stop),a_b_calc(3,:))
    
    figure()
    subplot(3,1,1); plot(t(start:stop),W_dot_b(1,:),t(start:stop),w_b_calc(1,:))
    subplot(3,1,2); plot(t(start:stop),W_dot_b(2,:),t(start:stop),w_b_calc(2,:))
    subplot(3,1,3); plot(t(start:stop),W_dot_b(3,:),t(start:stop),w_b_calc(3,:))
       
    figure()
    subplot(3,1,1); plot(t(start:stop),F_I(1,:),t(start:stop),F_R(1,:))
    subplot(3,1,2); plot(t(start:stop),F_I(2,:),t(start:stop),F_R(2,:))
    subplot(3,1,3); plot(t(start:stop),F_I(3,:),t(start:stop),F_R(3,:))
    
    figure()
    subplot(3,1,1); plot(t(start:stop),M_I(1,:),t(start:stop),M_R(1,:))
    subplot(3,1,2); plot(t(start:stop),M_I(2,:),t(start:stop),M_R(2,:))
    subplot(3,1,3); plot(t(start:stop),M_I(3,:),t(start:stop),M_R(3,:))
    
    
    % Averaging over the wingbeat:
    
    nr_wb = find(isnan(pathDB.wingbeat_time(:,1,seq_nr))==0, 1, 'last' );
    
    F_R_mean = zeros(3,nr_wb);
    
    M_R_mean = zeros(3,nr_wb);
    
    F_I_mean = zeros(3,nr_wb);
    
    M_I_mean = zeros(3,nr_wb);
    
    t_mean = zeros(1,nr_wb);
    
    t_wb_start = zeros(1,nr_wb);
    
    for j = 1:nr_wb
        
        wb_end = find(isnan(pathDB.wingbeat_time(j,:,seq_nr))==0, 1, 'last' );
        
        wb_t = pathDB.wingbeat_time(j,1:wb_end,seq_nr);
        
        F_R_mean(:,j) = [mean(F_R(1,wb_t)); mean(F_R(2,wb_t)); mean(F_R(3,wb_t))];
        
        M_R_mean(:,j) = [mean(M_R(1,wb_t)); mean(M_R(2,wb_t)); mean(M_R(3,wb_t))];
        
        F_I_mean(:,j) = [mean(F_I(1,wb_t)); mean(F_I(2,wb_t)); mean(F_I(3,wb_t))];
        
        M_I_mean(:,j) = [mean(M_I(1,wb_t)); mean(M_I(2,wb_t)); mean(M_I(3,wb_t))];
        
        t_mean(j) = t(start+round(mean(wb_t)));
        
        t_wb_start(j) = pathDB.wingbeat_time(j,1,seq_nr);
        
    end
        
    figure()
    subplot(3,1,1); plot(t_mean,F_I_mean(1,:),t_mean,F_R_mean(1,:),t_mean,F_I_mean(1,:)+F_R_mean(1,:))
    subplot(3,1,2); plot(t_mean,F_I_mean(2,:),t_mean,F_R_mean(2,:),t_mean,F_I_mean(2,:)+F_R_mean(2,:))
    subplot(3,1,3); plot(t_mean,F_I_mean(3,:),t_mean,F_R_mean(3,:),t_mean,F_I_mean(3,:)+F_R_mean(3,:))
    legend('Inertia','Resultant')
    
    figure()
    subplot(3,1,1); plot(t_mean,M_I_mean(1,:),t_mean,M_R_mean(1,:),t_mean,M_I_mean(1,:)+M_R_mean(1,:))
    subplot(3,1,2); plot(t_mean,M_I_mean(2,:),t_mean,M_R_mean(2,:),t_mean,M_I_mean(2,:)+M_R_mean(2,:))
    subplot(3,1,3); plot(t_mean,M_I_mean(3,:),t_mean,M_R_mean(3,:),t_mean,M_I_mean(3,:)+M_R_mean(3,:))
    legend('Inertia','Resultant')
    
    
    
%     figure()
%     subplot(3,1,1); plot(t(start:stop),M_I_strk(1,:),t(start-1+t_wb_start),zeros(1,nr_wb),'o')
%     subplot(3,1,2); plot(t(start:stop),M_I_strk(2,:),t(start-1+t_wb_start),zeros(1,nr_wb),'o')
%     subplot(3,1,3); plot(t(start:stop),M_I_strk(3,:),t(start-1+t_wb_start),zeros(1,nr_wb),'o')
    
    
    
    
    
%     R_wt_L = [0; -1; 0];
%     
%     R_wt_R = [0; 1; 0];
%     
%     Wt_L = zeros(3,stop-start+1);
%     
%     Wt_R = zeros(3,stop-start+1);
%     
%     RL1 = quat2matNEW([pathDB.qL1_filt2(start,seq_nr); pathDB.qL2_filt2(start,seq_nr); pathDB.qL3_filt2(start,seq_nr); pathDB.qL4_filt2(start,seq_nr)]);
%     
%     RR1 = quat2matNEW([pathDB.qR1_filt2(start,seq_nr); pathDB.qR2_filt2(start,seq_nr); pathDB.qR3_filt2(start,seq_nr); pathDB.qR4_filt2(start,seq_nr)]);
%     
%     Wt_L(:,1) = RL1'*R_wt_L;
%     
%     Wt_R(:,1) = RR1'*R_wt_R;
%     
%     for j = 2:(stop-start+1)
%         
%         Wt_L(:,j) = R_L_dot_test(:,:,j)'*Wt_L(:,j-1);
%         
%         
%         Wt_R(:,j) = R_R_dot_test(:,:,j)'*Wt_R(:,j-1);
%         
%         
%     end
%     
%     
%     figure()
%     plot3(Wt_L(1,:), Wt_L(2,:), Wt_L(3,:), 'r')
%     hold on
%     plot3(Wt_R(1,:), Wt_R(2,:), Wt_R(3,:), 'g')
%     axis equal
%     hold off
end
