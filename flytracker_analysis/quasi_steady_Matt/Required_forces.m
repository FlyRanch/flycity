function [ Req_f ] = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up_vect, f_vect, dt, seq_nr, q_b, v_st, a_st, w_st, w_dot_st )


    % Compute the required force generation for a given sets of wingbeats
    % with a given time-course of body and wing kinematics:

    Req_f = {};

    nr_wb = length(nr_points_vect);
    
   
    % Set up coefficient matrix for left and right wings:
    
    a_theta_L = a_sim.theta_L;
    a_eta_L = a_sim.eta_L;
    a_phi_L = a_sim.phi_L;
    
    a_theta_R = a_sim.theta_R;
    a_eta_R = a_sim.eta_R;
    a_phi_R = a_sim.phi_R;

    
    n_pol_theta = (size(a_theta_L,1)-2)/2;
    n_pol_eta = (size(a_eta_L,1)-2)/2;
    n_pol_phi = (size(a_phi_L,1)-2)/2;
    
   
    theta_L = zeros(sum(nr_points_vect)-nr_wb+1,1);
    eta_L = zeros(sum(nr_points_vect)-nr_wb+1,1);
    phi_L = zeros(sum(nr_points_vect)-nr_wb+1,1);
    
    theta_R = zeros(sum(nr_points_vect)-nr_wb+1,1);
    eta_R = zeros(sum(nr_points_vect)-nr_wb+1,1);
    phi_R = zeros(sum(nr_points_vect)-nr_wb+1,1);

    theta_dot_L = zeros(sum(nr_points_vect)-nr_wb+1,1);
    eta_dot_L = zeros(sum(nr_points_vect)-nr_wb+1,1);
    phi_dot_L = zeros(sum(nr_points_vect)-nr_wb+1,1);
    
    theta_dot_R = zeros(sum(nr_points_vect)-nr_wb+1,1);
    eta_dot_R = zeros(sum(nr_points_vect)-nr_wb+1,1);
    phi_dot_R = zeros(sum(nr_points_vect)-nr_wb+1,1);
    
    theta_ddot_L = zeros(sum(nr_points_vect)-nr_wb+1,1);
    eta_ddot_L = zeros(sum(nr_points_vect)-nr_wb+1,1);
    phi_ddot_L = zeros(sum(nr_points_vect)-nr_wb+1,1);
    
    theta_ddot_R = zeros(sum(nr_points_vect)-nr_wb+1,1);
    eta_ddot_R = zeros(sum(nr_points_vect)-nr_wb+1,1);
    phi_ddot_R = zeros(sum(nr_points_vect)-nr_wb+1,1);
    
%     size(theta_L)
    
    time = zeros(sum(nr_points_vect)-nr_wb+1,1);
    
    down = zeros(sum(nr_points_vect)-nr_wb+1,1);
    up = zeros(sum(nr_points_vect)-nr_wb+1,1);
    
    f_vect
    
    dt
    
    for i = 1:nr_wb
        
        [t_temp, X_theta] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_vect(i), nr_points_vect(i), t_vect(i), t_vect(i+1),0 );
        [~, X_eta] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_vect(i), nr_points_vect(i), t_vect(i), t_vect(i+1),0 );
        [~, X_phi] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_vect(i), nr_points_vect(i), t_vect(i), t_vect(i+1),0 );
        
        [~, X_dot_theta] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_vect(i), nr_points_vect(i), t_vect(i), t_vect(i+1),1 );
        [~, X_dot_eta] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_vect(i), nr_points_vect(i), t_vect(i), t_vect(i+1),1 );
        [~, X_dot_phi] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_vect(i), nr_points_vect(i), t_vect(i), t_vect(i+1),1 );
        
        [~, X_ddot_theta] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_vect(i), nr_points_vect(i), t_vect(i), t_vect(i+1),2 );
        [~, X_ddot_eta] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_vect(i), nr_points_vect(i), t_vect(i), t_vect(i+1),2 );
        [~, X_ddot_phi] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_vect(i), nr_points_vect(i), t_vect(i), t_vect(i+1),2 );
        
        i_int = (sum(nr_points_vect(1:(i-1)))-i+2):(sum(nr_points_vect(1:i))-i+1);
        
        down_int = (sum(nr_points_vect(1:(i-1)))-i+2):((sum(nr_points_vect(1:(i-1)))-i+2)+round(length(i_int)*down_up_vect)-1);
        up_int = ((sum(nr_points_vect(1:(i-1)))-i+2)+round(length(i_int)*down_up_vect)):(sum(nr_points_vect(1:i))-i);
        
        
        down(down_int) = ones(length(down_int),1);
        up(up_int) = ones(length(up_int),1);
        up(end) = 1;
        
        theta_L(i_int) = X_theta*a_theta_L(:,i);
        eta_L(i_int) = X_eta*a_eta_L(:,i);
        phi_L(i_int) = X_phi*a_phi_L(:,i);

        theta_R(i_int) = X_theta*a_theta_R(:,i);
        eta_R(i_int) = X_eta*a_eta_R(:,i);
        phi_R(i_int) = X_phi*a_phi_R(:,i);

%         theta_dot_L(i_int) = (nr_points_vect(i)-1)*X_dot_theta*a_theta_L(:,i);
%         eta_dot_L(i_int) = (nr_points_vect(i)-1)*X_dot_eta*a_eta_L(:,i);
%         phi_dot_L(i_int) = (nr_points_vect(i)-1)*X_dot_phi*a_phi_L(:,i);
% 
%         theta_dot_R(i_int) = (nr_points_vect(i)-1)*X_dot_theta*a_theta_R(:,i);
%         eta_dot_R(i_int) = (nr_points_vect(i)-1)*X_dot_eta*a_eta_R(:,i);
%         phi_dot_R(i_int) = (nr_points_vect(i)-1)*X_dot_phi*a_phi_R(:,i);
% 
%         theta_ddot_L(i_int) = (nr_points_vect(i)-1)^2*X_ddot_theta*a_theta_L(:,i);
%         eta_ddot_L(i_int) = (nr_points_vect(i)-1)^2*X_ddot_eta*a_eta_L(:,i);
%         phi_ddot_L(i_int) = (nr_points_vect(i)-1)^2*X_ddot_phi*a_phi_L(:,i);
% 
%         theta_ddot_R(i_int) = (nr_points_vect(i)-1)^2*X_ddot_theta*a_theta_R(:,i);
%         eta_ddot_R(i_int) = (nr_points_vect(i)-1)^2*X_ddot_eta*a_eta_R(:,i);
%         phi_ddot_R(i_int) = (nr_points_vect(i)-1)^2*X_ddot_phi*a_phi_R(:,i);

        theta_dot_L(i_int) = pi*f_vect(i)*X_dot_theta*a_theta_L(:,i);
        eta_dot_L(i_int) = pi*f_vect(i)*X_dot_eta*a_eta_L(:,i);
        phi_dot_L(i_int) = pi*f_vect(i)*X_dot_phi*a_phi_L(:,i);

        theta_dot_R(i_int) = pi*f_vect(i)*X_dot_theta*a_theta_R(:,i);
        eta_dot_R(i_int) = pi*f_vect(i)*X_dot_eta*a_eta_R(:,i);
        phi_dot_R(i_int) = pi*f_vect(i)*X_dot_phi*a_phi_R(:,i);

        theta_ddot_L(i_int) = (pi*f_vect(i))^2*X_ddot_theta*a_theta_L(:,i);
        eta_ddot_L(i_int) = (pi*f_vect(i))^2*X_ddot_eta*a_eta_L(:,i);
        phi_ddot_L(i_int) = (pi*f_vect(i))^2*X_ddot_phi*a_phi_L(:,i);

        theta_ddot_R(i_int) = (pi*f_vect(i))^2*X_ddot_theta*a_theta_R(:,i);
        eta_ddot_R(i_int) = (pi*f_vect(i))^2*X_ddot_eta*a_eta_R(:,i);
        phi_ddot_R(i_int) = (pi*f_vect(i))^2*X_ddot_phi*a_phi_R(:,i);
% 
%         theta_dot_L(i_int) = (1/(dt))*X_dot_theta*a_theta_L(:,i);
%         eta_dot_L(i_int) = (1/(dt))*X_dot_eta*a_eta_L(:,i);
%         phi_dot_L(i_int) = (1/(dt))*X_dot_phi*a_phi_L(:,i);
% 
%         theta_dot_R(i_int) = (1/(dt))*X_dot_theta*a_theta_R(:,i);
%         eta_dot_R(i_int) = (1/(dt))*X_dot_eta*a_eta_R(:,i);
%         phi_dot_R(i_int) = (1/(dt))*X_dot_phi*a_phi_R(:,i);
% 
%         theta_ddot_L(i_int) = (1/(dt))^2*X_ddot_theta*a_theta_L(:,i);
%         eta_ddot_L(i_int) = (1/(dt))^2*X_ddot_eta*a_eta_L(:,i);
%         phi_ddot_L(i_int) = (1/(dt))^2*X_ddot_phi*a_phi_L(:,i);
% 
%         theta_ddot_R(i_int) = (1/(dt))^2*X_ddot_theta*a_theta_R(:,i);
%         eta_ddot_R(i_int) = (1/(dt))^2*X_ddot_eta*a_eta_R(:,i);
%         phi_ddot_R(i_int) = (1/(dt))^2*X_ddot_phi*a_phi_R(:,i);
% 
%         theta_dot_L(i_int) = X_dot_theta*a_theta_L(:,i);
%         eta_dot_L(i_int) = X_dot_eta*a_eta_L(:,i);
%         phi_dot_L(i_int) = X_dot_phi*a_phi_L(:,i);
% 
%         theta_dot_R(i_int) = X_dot_theta*a_theta_R(:,i);
%         eta_dot_R(i_int) = X_dot_eta*a_eta_R(:,i);
%         phi_dot_R(i_int) = X_dot_phi*a_phi_R(:,i);
% 
%         theta_ddot_L(i_int) = X_ddot_theta*a_theta_L(:,i);
%         eta_ddot_L(i_int) = X_ddot_eta*a_eta_L(:,i);
%         phi_ddot_L(i_int) = X_ddot_phi*a_phi_L(:,i);
% 
%         theta_ddot_R(i_int) = X_ddot_theta*a_theta_R(:,i);
%         eta_ddot_R(i_int) = X_ddot_eta*a_eta_R(:,i);
%         phi_ddot_R(i_int) = X_ddot_phi*a_phi_R(:,i);

        time(i_int) = t_temp;
        
    end

    
    % Insert wing kinematics into the structure wing_kin:
    
    wing_kin = {};
    
    wing_kin.theta_L = theta_L; % rad
    wing_kin.eta_L = eta_L; % rad
    wing_kin.phi_L = phi_L; % rad
    
    wing_kin.theta_R = theta_R; % rad
    wing_kin.eta_R = eta_R; % rad
    wing_kin.phi_R = phi_R; % rad
    
    wing_kin.theta_dot_L = theta_dot_L; % rad/s
    wing_kin.eta_dot_L = eta_dot_L; % rad/s
    wing_kin.phi_dot_L = phi_dot_L; % rad/s
    
    wing_kin.theta_dot_R = theta_dot_R; % rad/s
    wing_kin.eta_dot_R = eta_dot_R; % rad/s
    wing_kin.phi_dot_R = phi_dot_R; % rad/s
    
    wing_kin.theta_ddot_L = theta_ddot_L; % rad/s^2
    wing_kin.eta_ddot_L = eta_ddot_L; % rad/s^2
    wing_kin.phi_ddot_L = phi_ddot_L; % rad/s^2
    
    wing_kin.theta_ddot_R = theta_ddot_R; % rad/s^2
    wing_kin.eta_ddot_R = eta_ddot_R; % rad/s^2
    wing_kin.phi_ddot_R = phi_ddot_R;  % rad/s^2
    
    % Transfer the wing kinematic angles into quaternions, angular
    % velocities and angular accelerations:
    
    [q_L,w_L,w_dot_L,q_R,w_R,w_dot_R,Rot_L,Rot_R] = Wingkin_2_quaternion( wing_kin );
    
    wing_kin.q_L = q_L;
    wing_kin.q_R = q_R;
    wing_kin.w_L = w_L; % rad/s
    wing_kin.w_R = w_R; % rad/s
    wing_kin.w_dot_L = w_dot_L; % rad/s^2
    wing_kin.w_dot_R = w_dot_R; % rad/s^2
    wing_kin.Rot_L = Rot_L;
    wing_kin.Rot_R = Rot_R;
    
    
    % Obtain body data for current fruit fly:
    
    wing_l = pathDB.wing_l(seq_nr); % mm
    
    mass_fly = 1.85e-6*(wing_l/3)^3; % kg
    
    nr_sect = 20;
    
    [cg_body, cg_L, cg_R, I_body, I_wing, I_v_wing, m_b, m_w, m_v_w, y_sect, chords] = cg_plus_Inertia(settings, pathDB, mass_fly, wing_l ,seq_nr ,nr_sect);

    
    Body = {};   
    Body.m_wl = mass_fly; % kg    
    Body.m_w = m_w; % kg    
    Body.m_v_w = m_v_w; % kg    
    Body.wing_l = wing_l; % mm    
    Body.joint_L = pathDB.joint_pos_L(:,seq_nr); % mm    
    Body.joint_R = pathDB.joint_pos_R(:,seq_nr); % mm    
    Body.cg_body = cg_body; % mm
    Body.cg_L = cg_L; % mm    
    Body.cg_R = cg_R; % mm    
    Body.I_body = I_body; % kg/mm^2    
    Body.I_wing = I_wing; % kg/mm^2    
    Body.I_v_wing = I_v_wing; % kg/mm^2    
    Body.y_sect = y_sect; % mm    
    Body.chords = chords; % mm
    
    
    n = sum(nr_points_vect)-nr_wb+1;
    
    Inertia_f = Inertia_forces( wing_kin, Body, n, v_st, a_st, w_st, w_dot_st);
    Aero_f = Aerodynamic_forces( wing_kin, Body, n, v_st, a_st, w_st, w_dot_st, down, up );
    gravity_f = Gravity_forces( q_b , mass_fly, n );
    
    Req_f.time = time;
    
    Req_f.F_I_acc = Inertia_f.F_acc_st;
    Req_f.F_I_vel = Inertia_f.F_vel_st;
    Req_f.M_I_acc = Inertia_f.M_acc_st;
    Req_f.M_I_vel = Inertia_f.M_vel_st;
    Req_f.F_I = Inertia_f.F_I_st;
    Req_f.M_I = Inertia_f.M_I_st;
    
    Req_f.c_pres_L = Aero_f.c_pres_L;
    Req_f.c_pres_R = Aero_f.c_pres_R;
    Req_f.F_a_wL = Aero_f.F_L;
    Req_f.F_a_wR = Aero_f.F_R;
    Req_f.F_a_L = Aero_f.F_L_st;
    Req_f.F_a_R = Aero_f.F_R_st;
    Req_f.M_a_L = Aero_f.M_L_st;
    Req_f.M_a_R = Aero_f.M_R_st;
    
    Req_f.F_g = gravity_f.Fg_st;
    
    Req_f.R_b = gravity_f.R_b;
    Req_f.R_L = wing_kin.Rot_L;
    Req_f.R_R = wing_kin.Rot_R;
    
    
    
    
    % Compute the mean forces and moments per wingbeat:
    
    F_I_acc_mean = zeros(3,nr_wb_man);  % Inertia_f.F_acc_st;
    F_I_vel_mean = zeros(3,nr_wb_man);  % Inertia_f.F_vel_st;
    M_I_acc_mean = zeros(3,nr_wb_man);  % Inertia_f.M_acc_st;
    M_I_vel_mean = zeros(3,nr_wb_man);  % Inertia_f.M_vel_st;
    F_I_mean = zeros(3,nr_wb_man);      % Inertia_f.F_I_st;
    M_I_mean = zeros(3,nr_wb_man);      % Inertia_f.M_I_st;
    F_a_L_mean = zeros(3,nr_wb_man);    % Aero_f.F_L_st;
    F_a_R_mean = zeros(3,nr_wb_man);    % Aero_f.F_R_st;
    M_a_L_mean = zeros(3,nr_wb_man);    % Aero_f.M_L_st;
    M_a_R_mean = zeros(3,nr_wb_man);    % Aero_f.M_R_st;
    F_g_mean = zeros(3,nr_wb_man);      % gravity_f.Fg_st;   

    
       for k = 1:nr_wb_man
           if k == 1
               b = 1:nr_points_vect(1);
           else           
               b = (sum(nr_points_vect(1:(k-1)))-(k-1)+1):(sum(nr_points_vect(1:k))-(k-1));
           end 
           
            F_I_acc_mean(:,k) = mean(Inertia_f.F_acc_st(:,b),2);
            F_I_vel_mean(:,k) = mean(Inertia_f.F_vel_st(:,b),2);
            M_I_acc_mean(:,k) = mean(Inertia_f.M_acc_st(:,b),2);
            M_I_vel_mean(:,k) = mean(Inertia_f.M_vel_st(:,b),2);
            F_I_mean(:,k) =     mean(Inertia_f.F_I_st(:,b),2);
            M_I_mean(:,k) =     mean(Inertia_f.M_I_st(:,b),2);
            F_a_L_mean(:,k) =   mean(Aero_f.F_L_st(:,b),2);
            F_a_R_mean(:,k) =   mean(Aero_f.F_R_st(:,b),2);
            M_a_L_mean(:,k) =   mean(Aero_f.M_L_st(:,b),2);
            M_a_R_mean(:,k) =   mean(Aero_f.M_R_st(:,b),2);
            F_g_mean(:,k) =     mean(gravity_f.Fg_st(:,b),2);   

       end    
       
    Req_f.F_I_acc_mean = F_I_acc_mean;
    Req_f.F_I_vel_mean = F_I_vel_mean;
    Req_f.M_I_acc_mean = M_I_acc_mean;
    Req_f.M_I_vel_mean = M_I_vel_mean;
    Req_f.F_I_mean = F_I_mean;
    Req_f.M_I_mean = M_I_mean;
    Req_f.F_a_L_mean = F_a_L_mean;
    Req_f.F_a_R_mean = F_a_R_mean;
    Req_f.M_a_L_mean = M_a_L_mean;
    Req_f.M_a_R_mean = M_a_R_mean;
    Req_f.F_g_mean = F_g_mean;

end

