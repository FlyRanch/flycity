function Dynamic_model( settings, pathDB )


%--------------------------------------------------------------------------

    % Use a state-space model of a fruit fly to compute the trajectory of
    % the fruit fly based on an inertial and aerodynamic model of the fruit
    % fly for a given set of wing kinematics.

%--------------------------------------------------------------------------

%%

%--------------------------------------------------------------------------

    %Input: Set up the wing kinematic pattern for a specified time with
    %specified time steps. Load the body and inertia parameters of the
    %specific fruit fly in the structures Pattern and Body.
    
%--------------------------------------------------------------------------

    

    dt = 1e-5; %[sec]
    
    seq_nr = 22; % sequence number
    
    nr_wb = 10; % number of modeled wingbeats
    
    trigger_wb = nr_wb+1; % wingbeat number when maneuver starts
    
    % Load the wing length and average wingbeat coefficients:
    
    a_avg = pathDB.a_avg.(char(['a_avg_' int2str(seq_nr)]));
    
    down_up_avg = pathDB.down_up_avg.(char(['down_up_avg_' int2str(seq_nr)]));
    
    f_avg = pathDB.f_avg.(char(['f_avg_' int2str(seq_nr)]));
    
    % Set up a vector with the wingbeat frequency for each wingbeat:
    
    f_vect = ones(nr_wb,1)*f_avg;
    
    % Set up a vector with the downstroke/upstroke ratio for each wingbeat:
    
    down_up_vect = ones(nr_wb,1)*down_up_avg;
    
    % Set up the time-line of the modeled sequence of wingbeats:
    
    t_vect = zeros(nr_wb+1,1);
    
    n_pol_theta = length(a_avg.theta_L1)-1;
    n_pol_eta = length(a_avg.eta_L1)-1;
    n_pol_phi = length(a_avg.phi_L1)-1;
    
    % Set up coefficient matrix for left and right wings:
    
    a_theta_L = zeros(2*(n_pol_theta+1),nr_wb);
    a_eta_L = zeros(2*(n_pol_eta+1),nr_wb);
    a_phi_L = zeros(2*(n_pol_phi+1),nr_wb);
    
    a_theta_R = zeros(2*(n_pol_theta+1),nr_wb);
    a_eta_R = zeros(2*(n_pol_eta+1),nr_wb);
    a_phi_R = zeros(2*(n_pol_phi+1),nr_wb);
    
    nr_points_vect = zeros(nr_wb,1);
    
    for i = 1:nr_wb
        
       nr_of_points = round((1/f_vect(i))/dt)+1;
       
       nr_points_vect(i) = nr_of_points;
        
       t_vect(i+1) = t_vect(i)+(nr_of_points-1)*dt;
    
       a_theta_L(:,i) = [a_avg.theta_LR1; a_avg.theta_LR2];
       a_eta_L(:,i) = [a_avg.eta_LR1; a_avg.eta_LR2];
       a_phi_L(:,i) = [a_avg.phi_LR1; a_avg.phi_LR2];
       
       a_theta_R(:,i) = [a_avg.theta_LR1; a_avg.theta_LR2];
       a_eta_R(:,i) = [a_avg.eta_LR1; a_avg.eta_LR2];
       a_phi_R(:,i) = [a_avg.phi_LR1; a_avg.phi_LR2];
        
    end
    
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
    
    size(theta_L)
    
    time = zeros(sum(nr_points_vect)-nr_wb+1,1);
    
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
        
        theta_L(i_int) = X_theta*a_theta_L(:,i);
        eta_L(i_int) = X_eta*a_eta_L(:,i);
        phi_L(i_int) = X_phi*a_phi_L(:,i);

        theta_R(i_int) = X_theta*a_theta_R(:,i);
        eta_R(i_int) = X_eta*a_eta_R(:,i);
        phi_R(i_int) = X_phi*a_phi_R(:,i);

        theta_dot_L(i_int) = (nr_points_vect(i)-1)*X_dot_theta*a_theta_L(:,i);
        eta_dot_L(i_int) = (nr_points_vect(i)-1)*X_dot_eta*a_eta_L(:,i);
        phi_dot_L(i_int) = (nr_points_vect(i)-1)*X_dot_phi*a_phi_L(:,i);

        theta_dot_R(i_int) = (nr_points_vect(i)-1)*X_dot_theta*a_theta_R(:,i);
        eta_dot_R(i_int) = (nr_points_vect(i)-1)*X_dot_eta*a_eta_R(:,i);
        phi_dot_R(i_int) = (nr_points_vect(i)-1)*X_dot_phi*a_phi_R(:,i);

        theta_ddot_L(i_int) = (nr_points_vect(i)-1)^2*X_ddot_theta*a_theta_L(:,i);
        eta_ddot_L(i_int) = (nr_points_vect(i)-1)^2*X_ddot_eta*a_eta_L(:,i);
        phi_ddot_L(i_int) = (nr_points_vect(i)-1)^2*X_ddot_phi*a_phi_L(:,i);

        theta_ddot_R(i_int) = (nr_points_vect(i)-1)^2*X_ddot_theta*a_theta_R(:,i);
        eta_ddot_R(i_int) = (nr_points_vect(i)-1)^2*X_ddot_eta*a_eta_R(:,i);
        phi_ddot_R(i_int) = (nr_points_vect(i)-1)^2*X_ddot_phi*a_phi_R(:,i);
        
        time(i_int) = t_temp;
        
        if i >= trigger_wb
            
            'maneuver'
            
        end
        
    end
    
    % Plot wing kinematics:
    
    figure()
    subplot(2,1,1); plot(time,theta_L,time,eta_L,time,phi_L)
    title('left wing kinematics')
    xlabel('time [s]')
    ylabel('angles rad')
    legend('\theta','\eta','\phi')
    subplot(2,1,2); plot(time,theta_R,time,eta_R,time,phi_R)
    title('right wing kinematics')
    xlabel('time [s]')
    ylabel('angles rad')
    legend('\theta','\eta','\phi')
    
    figure()
    subplot(2,1,1); plot(time,theta_dot_L,time,eta_dot_L,time,phi_dot_L)
    title('left wing kinematics')
    xlabel('time [s]')
    ylabel('angular velocity [rad/s]')
    legend('\theta^*','\eta^*','\phi^*')
    subplot(2,1,2); plot(time,theta_dot_R,time,eta_dot_R,time,phi_dot_R)
    title('right wing kinematics')
    xlabel('time [s]')
    ylabel('angular velocity [rad/s]')
    legend('\theta^*','\eta^*','\phi^*')
    
    figure()
    subplot(2,1,1); plot(time,theta_ddot_L,time,eta_ddot_L,time,phi_ddot_L)
    title('left wing kinematics')
    xlabel('time [s]')
    ylabel('angular acceleration [rad/s^2]')
    legend('\theta^{**}','\eta^{**}','\phi^{**}')
    subplot(2,1,2); plot(time,theta_ddot_R,time,eta_ddot_R,time,phi_ddot_R)
    title('right wing kinematics')
    xlabel('time [s]')
    ylabel('angular acceleration [rad/s^2]')
    legend('\theta^{**}','\eta^{**}','\phi^{**}')
    
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
    
    
    % Plot quaternion wing kinematics:
    
    figure()
    subplot(2,1,1); plot(time,q_L(1,:),time,q_L(2,:),time,q_L(3,:),time,q_L(4,:))
    title('left wing kinematics')
    xlabel('time [s]')
    legend('q_1','q_2','q_3','q_4')
    subplot(2,1,2); plot(time,q_R(1,:),time,q_R(2,:),time,q_R(3,:),time,q_R(4,:))
    title('right wing kinematics')
    xlabel('time [s]')
    legend('q_1','q_2','q_3','q_4')
    
    figure()
    subplot(2,1,1); plot(time,w_L(1,:),time,w_L(2,:),time,w_L(3,:))
    title('left wing kinematics')
    xlabel('time [s]')
    ylabel('angular velocity [rad/s]')
    legend('w_x','w_y','w_z')
    subplot(2,1,2); plot(time,w_R(1,:),time,w_R(2,:),time,w_R(3,:))
    title('right wing kinematics')
    xlabel('time [s]')
    ylabel('angular velocity [rad/s]')
    legend('w_x','w_y','w_z')
    
    figure()
    subplot(2,1,1); plot(time,w_dot_L(1,:),time,w_dot_L(2,:),time,w_dot_L(3,:))
    title('left wing kinematics')
    xlabel('time [s]')
    ylabel('angular acceleration [rad/s^2]')
    legend('w_x^*','w_y^*','w_z^*')
    subplot(2,1,2); plot(time,w_dot_R(1,:),time,w_dot_R(2,:),time,w_dot_R(3,:))
    title('right wing kinematics')
    xlabel('time [s]')
    ylabel('angular acceleration [rad/s^2]')
    legend('w_x^*','w_y^*','w_z^*')
    
    figure()
    hold on
    for j = 1:5:length(q_L(1,:))
        xyz_L = Rot_L(:,:,j)'*[0; -1; 0];
        xyz_R = Rot_R(:,:,j)'*[0; 1; 0];
        plot3(xyz_L(1), -xyz_L(2), -xyz_L(3),'r')
        plot3(xyz_R(1), -xyz_R(2), -xyz_R(3),'g')
    end
    axis equal
    hold off

    
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
    
    
%     mass_fly
%     m_w
%     m_v_w
%     wing_l
%     pathDB.joint_pos_L(:,seq_nr)
%     pathDB.joint_pos_R(:,seq_nr)
%     cg_body
%     cg_L
%     cg_R
%     I_body
%     I_wing
%     I_v_wing
%     y_sect
%     chords


    % Test whether inerita model works:

    Inertia_test( wing_kin, Body, time )
    
    
    %----------------------------------------------------------------------
    
    
    %%
    
    
    
    %----------------------------------------------------------------------
    
    % Simulation: Use the wing kinematics and body parameters in
    % combination with an aerodynamic model to compute the forces and
    % moments at the c.g. of the fruit fly. Use a state-space system and a
    % Runge-Kutta scheme to predict the trajectory based on these forces
    % and moments.
    
    %----------------------------------------------------------------------
    
%     t = 0;
%     
%     dt2 = 2*dt;
%     
%     t_end = floor(length(theta_L)/2)*dt2;
%     
%     while t <= t_end
%         
%         t = t+dt2;
%         
%     end
    
    %----------------------------------------------------------------------
    
    %%
    
end

