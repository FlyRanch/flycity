function [F_r, M_r, T_b_trans, T_b_rot, T_L_trans, T_L_rot, T_R_trans, T_R_rot, F_r_avg, M_r_avg, t_avg, M_r_meas, M_r_meas_avg, F_r_global_avg] = Dynamic_model_fruit_fly(settings, pathDB, seq_nr)


%[a_body_calc, w_body_calc, a_L_calc, w_L_calc, a_R_calc, w_R_calc, T_trans_body, T_rot_body, V_body, T_trans_L, T_rot_L, V_L, T_trans_R, T_rot_R, V_R ] = Dynamic_model_fruit_fly(settings, pathDB, i, seq_nr)


    % Function that calculates the forces and moments on a fruit fly at
    % time t(i) and returns the accelerations and angular velocities at the
    % fruit fly body and wings for those forces and moments. The actual
    % measured accelerations and angular velocities are compared with the
    % calculated ones to derive to what extend the model for generated
    % inertial and aerodynamic forces is capable of modelling the reality.
    
    
    dt = pathDB.t(2)-pathDB.t(1);
    
    
    % Fruit fly mass based on wing length:
    
    wing_l = pathDB.wing_l(seq_nr);
    
    m_wl = 1e-6*0.1078*wing_l^3.008;
    
    
    % Calculate body density (3D):
    
    [cg_body, cg_L, cg_R, I_body, I_wing, I_v_wing, m_b, m_w, m_v_w ] = cg_plus_Inertia(settings, pathDB, m_wl, wing_l ,seq_nr );

    
    start = find(isnan(pathDB.x(:,seq_nr))==0, 1 );
    stop = find(isnan(pathDB.x(:,seq_nr))==0, 1, 'last' );
    
    
        
    F_r = zeros(3,stop-start+1);
    
    M_r = zeros(3,stop-start+1);
    
    
    T_b_trans = zeros(stop-start+1,1);
    
    T_b_rot = zeros(stop-start+1,1);
    
    T_L_trans = zeros(stop-start+1,1);
    
    T_L_rot = zeros(stop-start+1,1);
    
    T_R_trans = zeros(stop-start+1,1);
    
    T_R_rot = zeros(stop-start+1,1);
    
%     
%     I_r = zeros(3,3,stop-start+1);
    

    % Give wingbeat-averaged forces
    
    wb_time_end = find(isnan(pathDB.wingbeat_time(:,1,seq_nr))==0, 1, 'last' );
    
    wb_time = pathDB.wingbeat_time(1:wb_time_end,:,seq_nr);


    
    M_r_meas = zeros(3,stop-start+1);
    
    F_r_global = zeros(3,stop-start+1);

    for i = (start+1):stop
   
    
    
    % Load rotation matrices body and wings:
    
    
    qB = [pathDB.qb1_filt(i,seq_nr); pathDB.qb2_filt(i,seq_nr); pathDB.qb3_filt(i,seq_nr); pathDB.qb4_filt(i,seq_nr)];
    
    qL = [pathDB.qL1_filt2(i,seq_nr); pathDB.qL2_filt2(i,seq_nr); pathDB.qL3_filt2(i,seq_nr); pathDB.qL4_filt2(i,seq_nr)];
    
    qR = [pathDB.qR1_filt2(i,seq_nr); pathDB.qR2_filt2(i,seq_nr); pathDB.qR3_filt2(i,seq_nr); pathDB.qR4_filt2(i,seq_nr)];
    
    R_b = quat2matNEW(qB);
    
    R_L = quat2matNEW(qL);
    
    R_R = quat2matNEW(qR);
    
    qB_m1 = [pathDB.qb1_filt(i-1,seq_nr); pathDB.qb2_filt(i-1,seq_nr); pathDB.qb3_filt(i-1,seq_nr); pathDB.qb4_filt(i-1,seq_nr)];
    
    qL_m1 = [pathDB.qL1_filt2(i-1,seq_nr); pathDB.qL2_filt2(i-1,seq_nr); pathDB.qL3_filt2(i-1,seq_nr); pathDB.qL4_filt2(i-1,seq_nr)];
    
    qR_m1 = [pathDB.qR1_filt2(i-1,seq_nr); pathDB.qR2_filt2(i-1,seq_nr); pathDB.qR3_filt2(i-1,seq_nr); pathDB.qR4_filt2(i-1,seq_nr)];
    
%     R_b_m1 = quat2matNEW(qB_m1);
%     
%     R_L_m1 = quat2matNEW(qL_m1);
%     
%     R_R_m1 = quat2matNEW(qR_m1);

        
    
    % Load body velocity, acceleration and angular velocity:
    
    u_b = [pathDB.u_filt(i,seq_nr); pathDB.v_filt(i,seq_nr); pathDB.w_filt(i,seq_nr)];
    
    a_b = [pathDB.ax_filt(i,seq_nr); pathDB.ay_filt(i,seq_nr); pathDB.az_filt(i,seq_nr)];
    
    w_b = [pathDB.b_omega1(i,seq_nr); pathDB.b_omega2(i,seq_nr); pathDB.b_omega3(i,seq_nr)];
    
    
    % Load wing (left and right) angular velocities and joint locations:
    
    w_L = [pathDB.omega1_L(i,seq_nr); pathDB.omega2_L(i,seq_nr); pathDB.omega3_L(i,seq_nr)];
    
    w_R = [pathDB.omega1_R(i,seq_nr); pathDB.omega2_R(i,seq_nr); pathDB.omega3_R(i,seq_nr)];
    
    J_L = pathDB.joint_pos_L(:,seq_nr);
    
    J_R = pathDB.joint_pos_R(:,seq_nr);
    
    
    % Load aerodynamic forces of the wings:
    
    F_aero_L = pathDB.F_L(i,:,seq_nr)';
    F_aero_R = pathDB.F_R(i,:,seq_nr)';
    
    F_aero_L_arm = [0; pathDB.F_L_arm(i,seq_nr); 0];
    F_aero_R_arm = [0; pathDB.F_R_arm(i,seq_nr); 0];
    
    
    % Calculate angular accelerations body and wings based on two
    % subsequent frames:
    
    w_dot_b = [(w_b-[pathDB.b_omega1(i-1,seq_nr); pathDB.b_omega2(i-1,seq_nr); pathDB.b_omega3(i-1,seq_nr)])./dt];
    
    w_dot_L = R_L'*w_dot_b+[(w_L-[pathDB.omega1_L(i-1,seq_nr); pathDB.omega2_L(i-1,seq_nr); pathDB.omega3_L(i-1,seq_nr)])./dt];
    
    w_dot_R = R_R'*w_dot_b+[(w_R-[pathDB.omega1_R(i-1,seq_nr); pathDB.omega2_R(i-1,seq_nr); pathDB.omega3_R(i-1,seq_nr)])./dt];

% 
%     w_dot_b = [0; 0; 0];
%     
%     w_dot_L = [0; 0; 0];
%     
%     w_dot_R = [0; 0; 0];
    
    % Calculate (rotational) inertial moments on the body and wings:
    
    
    M_I_b = 1e-3*(I_body*w_dot_b+cross(w_b,I_body*w_b)); % [N*mm]
    
    M_I_L = 1e-3*(I_v_wing*w_dot_L+cross([w_L+R_L'*w_b],I_v_wing*[w_L+R_L'*w_b])); % [N*mm]
    
    M_I_R = 1e-3*(I_v_wing*w_dot_R+cross([w_R+R_R'*w_b],I_v_wing*[w_R+R_R'*w_b])); % [N*mm]

%     M_I_b = 1e-3*(I_body*w_dot_b+cross(w_b,I_body*w_b)); % [N*mm]
%     
%     M_I_L = 1e-3*(I_v_wing*w_dot_L+cross([w_L],I_v_wing*[w_L])); % [N*mm]
%     
%     M_I_R = 1e-3*(I_v_wing*w_dot_R+cross([w_R],I_v_wing*[w_R])); % [N*mm]

    
    
    % Calculate gravitational forces body and wings:
    
    g = 9.81; % [m/s^2]
    
    F_g_b = R_b'*[0; 0; -m_b*g]; %[N]
    
    F_g_L = R_L'*R_b'*[0; 0; -m_w*g]; %[N]
    
    F_g_R = R_R'*R_b'*[0; 0; -m_w*g]; %[N]
    
    
    % Calculate (translational) acceleration forces on body and wings:
    
%     F_a_b = 1e-3*m_b*(a_b+cross(w_dot_b,cg_body)+cross(w_b,cross(w_b,cg_body))); %[N]

%     F_a_L = 1e-3*m_v_w*(R_L'*(a_b+cross(w_dot_b,J_L)+cross(w_b,cross(w_b,J_L)))+cross(w_dot_L,cg_L)+cross(w_L,cross(w_L,cg_L))); %[N]
%     
%     F_a_R = 1e-3*m_v_w*(R_R'*(a_b+cross(w_dot_b,J_R)+cross(w_b,cross(w_b,J_R)))+cross(w_dot_R,cg_R)+cross(w_R,cross(w_R,cg_R))); %[N]


    F_a_b = 1e-3*m_b*(cross(w_dot_b,cg_body)+cross(w_b,cross(w_b,cg_body))); %[N]

    F_a_L = 1e-3*m_v_w*(R_L'*(cross(w_dot_b,J_L)+cross(w_b,cross(w_b,J_L)))+cross(w_dot_L,cg_L)+cross(w_L,cross(w_L,cg_L))); %[N]
    
    F_a_R = 1e-3*m_v_w*(R_R'*(cross(w_dot_b,J_R)+cross(w_b,cross(w_b,J_R)))+cross(w_dot_R,cg_R)+cross(w_R,cross(w_R,cg_R))); %[N]
    
    
%     F_r(:,i-start+1) = (R_L*(F_aero_L + F_g_L)+R_R*(F_aero_R + F_g_R) + F_g_b);
%     
%     M_r(:,i-start+1) = cross(R_L*(F_aero_L + F_g_L + F_a_L + cross(cg_L,M_I_L)),[J_L-cg_body])+cross(R_R*(F_aero_R + F_g_R + F_a_R + cross(cg_R,M_I_R)),[J_R-cg_body]);
     
    F_r(:,i-start+1) = (R_L*(F_aero_L + F_g_L + F_a_L + cross(cg_L,M_I_L))+R_R*(F_aero_R + F_g_R + F_a_R + cross(cg_R,M_I_R)) + F_g_b + F_a_b);

%      M_r(:,i-start+1) = cross(R_L*(F_aero_L + F_g_L + F_a_L + cross(cg_L,M_I_L)),[J_L-cg_body])+cross(R_R*(F_aero_R + F_g_R + F_a_R + cross(cg_R,M_I_R)),[J_R-cg_body])+M_I_b;

     M_r(:,i-start+1) = (cross([J_L-cg_body],R_L*(F_aero_L + F_g_L + F_a_L + cross(cg_L,M_I_L)))+cross([J_R-cg_body],R_R*(F_aero_R + F_g_R + F_a_R + cross(cg_R,M_I_R))));

%     M_r(:,i-start+1) = -(cross(R_L*(F_aero_L + F_g_L + F_a_L),[J_L-cg_body])+cross(R_R*(F_aero_R + F_g_R + F_a_R),[J_R-cg_body]));  
    

    F_r_global(:,i-start+1) = R_b*F_r(:,i-start+1);

    
    T_b_trans(i-start+1) = 0.5*m_b*norm(u_b)^2;
    
    T_b_rot(i-start+1) = 0.5*(R_b*w_b)'*I_body*(R_b*w_b);
    
    T_L_trans(i-start+1) = 0.5*m_w*norm(u_b+R_b*R_L*cross(R_L'*w_b+w_L,cg_L))^2;
    
    T_L_rot(i-start+1) = 0.5*(R_b*w_b+R_b*R_L*w_L)'*I_v_wing*(R_b*w_b+R_b*R_L*w_L);
    
    T_R_trans(i-start+1) = 0.5*m_w*norm(u_b+R_b*R_R*cross(R_R'*w_b+w_R,cg_R))^2;
    
    T_R_rot(i-start+1) = 0.5*(R_b*w_b+R_b*R_R*w_R)'*I_v_wing*(R_b*w_b+R_b*R_R*w_R);
    

    % Total moment of inertia of the fruitfly is equal to body moment of
    % inertia as the body moment of inertia is 2 orders of magnitude bigger
    % than the moment of inertia of the two wings combined.

    
    % Calculate resultant moments of the fruitfly as a function of its
    % moments of inertia and its angular velocities and angular
    % accelerations:
    
%     M_r_meas(:,i-start+1) = 1e-3*(I_body*w_dot_b + cross(w_b,I_body*w_b) + I_wing*w_dot_L + cross(w_L,I_wing*w_L) + I_wing*w_dot_R + cross(w_R,I_wing*w_R) + I_wing*R_L*w_dot_b ); % + cross(R_L*w_b,I_wing*R_L*w_b) + I_wing*R_R*w_dot_b + cross(R_R*w_b,I_wing*R_R*w_b)
   
    M_r_meas(:,i-start+1) = 1e-3*(I_body*w_dot_b + cross(w_b,I_body*w_b) + R_L*(I_wing*w_dot_L + cross(w_L,I_wing*w_L)) + R_R*(I_wing*w_dot_R + cross(w_R,I_wing*w_R)) );
     
%     % Calculate moment of inertia of the whole fruitfly, the center of
%     % gravity can be assumed to be the body center of gravity:
%     
%     I_r(:,:,i-start+1) = m_v_w*(J_L-cg_body+R_L*cg_L).^2+m_v_w*(J_R-cg_body+R_R*cg_R).^2+I_body

     
     
     
     
    
%     % Plot forces in body model:
%     
%     
%     % Load body model:
% 
%     [x, y, z] = load_body_model(settings, seq_nr, [0 0 0 qB' qL' qR']);
%     
%     % Gravity:
%     
%     quiv_scale = 1e5;
%     
%     F_g_b_quiv = quiv_scale.*R_b*F_g_b;
%     
%     F_g_b_pos = R_b*cg_body;
%     
%     % Aero
%     
%     F_aero_L_quiv = quiv_scale.*R_b*R_L*F_aero_L;
%     
%     F_aero_R_quiv = quiv_scale.*R_b*R_R*F_aero_R;
%     
%     F_aero_L_pos = R_b*(J_L+R_L*F_aero_L_arm);
%     
%     F_aero_R_pos = R_b*(J_R+R_R*F_aero_R_arm);
%     
%     % Acceleration
%     
%     F_a_b_quiv = quiv_scale*R_b*F_a_b;
%     
%     F_a_L_quiv = quiv_scale*R_b*R_L*F_a_L;
%     
%     F_a_R_quiv = quiv_scale*R_b*R_L*F_a_R;
%     
%     F_g_L_pos = R_b*(J_L+R_L*cg_L);
%     
%     F_g_R_pos = R_b*(J_R+R_R*cg_R);
% 
%     
% %     figure()
% %     hold on
% %     for k = 2:3
% %     surf(x{k},y{k},z{k},'facecolor',[0.3 0.3 0.3],'edgecolor','k','facelighting','phong')
% %     end
% %     quiver3(F_g_b_pos(1),F_g_b_pos(2),F_g_b_pos(3),F_g_b_quiv(1),F_g_b_quiv(2),F_g_b_quiv(3),'Color','k')
% % %     quiver3(F_g_L_pos(1),F_g_L_pos(2),F_g_L_pos(3),F_g_L_quiv(1),F_g_L_quiv(2),F_g_L_quiv(3),'Color','k')
% % %     quiver3(F_g_R_pos(1),F_g_R_pos(2),F_g_R_pos(3),F_g_R_quiv(1),F_g_R_quiv(2),F_g_R_quiv(3),'Color','k')
% %     quiver3(F_aero_L_pos(1),F_aero_L_pos(2),F_aero_L_pos(3),F_aero_L_quiv(1),F_aero_L_quiv(2),F_aero_L_quiv(3),'Color','b')
% %     quiver3(F_aero_R_pos(1),F_aero_R_pos(2),F_aero_R_pos(3),F_aero_R_quiv(1),F_aero_R_quiv(2),F_aero_R_quiv(3),'Color','b')
% %     quiver3(F_g_b_pos(1),F_g_b_pos(2),F_g_b_pos(3),F_a_b_quiv(1),F_a_b_quiv(2),F_a_b_quiv(3),'Color','r')
% %     quiver3(F_g_L_pos(1),F_g_L_pos(2),F_g_L_pos(3),F_a_L_quiv(1),F_a_L_quiv(2),F_a_L_quiv(3),'Color','r')
% %     quiver3(F_g_R_pos(1),F_g_R_pos(2),F_g_R_pos(3),F_a_R_quiv(1),F_a_R_quiv(2),F_a_R_quiv(3),'Color','r')
% %     xlabel('x')
% %     ylabel('y')
% %     zlabel('z')
% %     axis equal
% %     hold off
% 
% 
%     
%     % Calculate forces generated on wing joints:
% 
%     
%     F_Jw_L = quiv_scale*R_b*R_L*(F_aero_L + F_g_L + F_a_L + cross(cg_L,M_I_L));
%     
%     F_Jw_R = quiv_scale*R_b*R_R*(F_aero_R + F_g_R + F_a_R + cross(cg_R,M_I_R));    
%     
%     J_L_pos = R_b*J_L;
%     
%     J_R_pos = R_b*J_R;
%     
%     
%     F_r = quiv_scale*R_b*(R_L*(F_aero_L + F_g_L + F_a_L + cross(cg_L,M_I_L))+R_R*(F_aero_R + F_g_R + F_a_R + cross(cg_R,M_I_R))+F_g_b+F_a_b);
% 
%     
%     figure()
%     hold on
%     for k = 1:3
%     surf(x{k},y{k},z{k},'facecolor',[0.3 0.3 0.3],'edgecolor','k','facelighting','phong')
%     end
%     quiver3(F_g_b_pos(1),F_g_b_pos(2),F_g_b_pos(3),F_g_b_quiv(1),F_g_b_quiv(2),F_g_b_quiv(3),'Color','k')
%     quiver3(F_g_b_pos(1),F_g_b_pos(2),F_g_b_pos(3),F_a_b_quiv(1),F_a_b_quiv(2),F_a_b_quiv(3),'Color','r')
%     quiver3(F_aero_L_pos(1),F_aero_L_pos(2),F_aero_L_pos(3),F_aero_L_quiv(1),F_aero_L_quiv(2),F_aero_L_quiv(3),'Color','b')
%     quiver3(F_aero_R_pos(1),F_aero_R_pos(2),F_aero_R_pos(3),F_aero_R_quiv(1),F_aero_R_quiv(2),F_aero_R_quiv(3),'Color','b')
%     quiver3(J_L_pos(1),J_L_pos(2),J_L_pos(3),F_Jw_L(1),F_Jw_L(2),F_Jw_L(3),'Color','b')
%     quiver3(J_R_pos(1),J_R_pos(2),J_R_pos(3),F_Jw_R(1),F_Jw_R(2),F_Jw_R(3),'Color','b')
%     quiver3(F_g_b_pos(1),F_g_b_pos(2),F_g_b_pos(3),F_r(1),F_r(2),F_r(3),'Color','g')
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     view(45,45)
%     axis equal
%     hold off
%     
%     pause
    
%     1e3*-R_b*(R_L*(F_aero_L) +R_R*(F_aero_R )+F_g_b)./(m_wl)
%     
%     a_b
    
    end
    
    
    % Calculate wingbeat averaged forces and moments:
    
    F_r_avg = zeros(3,size(wb_time,1));
    
    M_r_avg = zeros(3,size(wb_time,1));
    
    t_avg = zeros(1,size(wb_time,1));
    
    M_r_meas_avg = zeros(3,size(wb_time,1));
    
    F_r_global_avg = zeros(3,size(wb_time,1));
    
    for k = 1:size(wb_time,1)
        
        if k < size(wb_time,1)
        
        F_r_avg(1,k) = mean(F_r(1,wb_time(k,1):(wb_time(k+1,1)-1)));
        
        F_r_avg(2,k) = mean(F_r(2,wb_time(k,1):(wb_time(k+1,1)-1)));
        
        F_r_avg(3,k) = mean(F_r(3,wb_time(k,1):(wb_time(k+1,1)-1)));
        
        M_r_avg(1,k) = mean(M_r(1,wb_time(k,1):(wb_time(k+1,1)-1)));
        
        M_r_avg(2,k) = mean(M_r(2,wb_time(k,1):(wb_time(k+1,1)-1)));
        
        M_r_avg(3,k) = mean(M_r(3,wb_time(k,1):(wb_time(k+1,1)-1)));
        
        F_r_global_avg(1,k) = mean(F_r_global(1,wb_time(k,1):(wb_time(k+1,1)-1)));
        
        F_r_global_avg(2,k) = mean(F_r_global(2,wb_time(k,1):(wb_time(k+1,1)-1)));
        
        F_r_global_avg(3,k) = mean(F_r_global(3,wb_time(k,1):(wb_time(k+1,1)-1)));        
        
        M_r_meas_avg(1,k) = mean(M_r_meas(1,wb_time(k,1):(wb_time(k+1,1)-1)));
        
        M_r_meas_avg(2,k) = mean(M_r_meas(2,wb_time(k,1):(wb_time(k+1,1)-1)));
        
        M_r_meas_avg(3,k) = mean(M_r_meas(3,wb_time(k,1):(wb_time(k+1,1)-1)));
        
        t_avg(k) = round((wb_time(k,1)+(wb_time(k+1,1)-1))/2);
        
        elseif k == size(wb_time,1)
            
        wb_t_end = find(isnan(wb_time(k,:))==0, 1, 'last' );
            
        F_r_avg(1,k) = mean(F_r(1,wb_time(k,1):wb_t_end));
        
        F_r_avg(2,k) = mean(F_r(2,wb_time(k,1):wb_t_end));
        
        F_r_avg(3,k) = mean(F_r(3,wb_time(k,1):wb_t_end));
        
        M_r_avg(1,k) = mean(M_r(1,wb_time(k,1):wb_t_end));
        
        M_r_avg(2,k) = mean(M_r(2,wb_time(k,1):wb_t_end));
        
        M_r_avg(3,k) = mean(M_r(3,wb_time(k,1):wb_t_end));
        
        F_r_global_avg(1,k) = mean(F_r_global(1,wb_time(k,1):wb_t_end));
        
        F_r_global_avg(2,k) = mean(F_r_global(2,wb_time(k,1):wb_t_end));
        
        F_r_global_avg(3,k) = mean(F_r_global(3,wb_time(k,1):wb_t_end));
        
        M_r_meas_avg(1,k) = mean(M_r_meas(1,wb_time(k,1):wb_t_end));
        
        M_r_meas_avg(2,k) = mean(M_r_meas(2,wb_time(k,1):wb_t_end));
        
        M_r_meas_avg(3,k) = mean(M_r_meas(3,wb_time(k,1):wb_t_end));        
        
        t_avg(k) = round((wb_time(k,1)+(wb_time(k,wb_t_end)))/2);
            
        end
        
    end
    
    clear k

    
  
end

