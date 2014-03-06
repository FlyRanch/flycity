function Dynamic_test(settings, pathDB)

    % Program that tests the dynamics of a fruit fly model for a given
    % sequence of wing kinematics.
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Create the wingbeat pattern:
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Get an average wingbeat:
    
    seq_nr = 1;
    
%     n_pol_theta = 13; % Order of used polynomials
%     n_pol_eta = 16; % Order of used polynomials
%     n_pol_phi = 8; % Order of used polynomials

    n_pol_theta = 20; % Order of used polynomials
    n_pol_eta = 20; % Order of used polynomials
    n_pol_phi = 10; % Order of used polynomials
    
    [a_fit_L,a_fit_R,a_avg_L,a_avg_R,a_avg_LR,f_avg,down_up,trigger_wb] = Standard_wingbeat( settings, pathDB, seq_nr , n_pol_theta, n_pol_eta, n_pol_phi);
    
    wb_time = 1/f_avg;
    
    wb_sample_nr = 101; % Always take an odd number
    
    dt = wb_time/((wb_sample_nr-1)/2);
    
    sample_loc = -1:(2/(wb_sample_nr-1)):1; % Legendre sample locations
    
    wb_nr = 10; % Number of simulated wingbeats
    
    nr_samples = wb_sample_nr*wb_nr-(wb_nr-1); % Total number of time steps
    
        
    PN_theta = Legendre_polynomial(n_pol_theta,3,sample_loc);
    PN_eta = Legendre_polynomial(n_pol_eta,3,sample_loc);
    PN_phi = Legendre_polynomial(n_pol_phi,3,sample_loc);  
    
    theta = zeros(nr_samples,2);
    eta = zeros(nr_samples,2);
    phi = zeros(nr_samples,2);
    
    theta_dot = zeros(nr_samples,2);
    eta_dot = zeros(nr_samples,2);
    phi_dot = zeros(nr_samples,2);

    theta_ddot = zeros(nr_samples,2);
    eta_ddot = zeros(nr_samples,2);
    phi_ddot = zeros(nr_samples,2);
    
    down = nan(nr_samples);
    
    up = nan(nr_samples);

    
    for i = 1:wb_nr
        
        if i < wb_nr
            
            n = (((i-1)*wb_sample_nr-(i-1))+1):(i*wb_sample_nr-i);
            
            down_str_l = round(down_up*(wb_sample_nr-1)/(1+down_up));
            
            m = (((i-1)*wb_sample_nr-(i-1))+1):(((i-1)*wb_sample_nr-(i-1))+down_str_l);
            
            l = (((i-1)*wb_sample_nr-(i-1))+down_str_l+1):(i*wb_sample_nr-i);
            
            theta(n,1) = PN_theta(:,1:(end-1),1)'*a_avg_LR.theta;
            eta(n,1) = PN_eta(:,1:(end-1),1)'*a_avg_LR.eta;
            phi(n,1) = PN_phi(:,1:(end-1),1)'*a_avg_LR.phi;
    
            theta_dot(n,1) = PN_theta(:,1:(end-1),2)'*a_avg_LR.theta;
            eta_dot(n,1) = PN_eta(:,1:(end-1),2)'*a_avg_LR.eta;
            phi_dot(n,1) = PN_phi(:,1:(end-1),2)'*a_avg_LR.phi;

            theta_ddot(n,1) = PN_theta(:,1:(end-1),3)'*a_avg_LR.theta;
            eta_ddot(n,1) = PN_eta(:,1:(end-1),3)'*a_avg_LR.eta;
            phi_ddot(n,1) = PN_phi(:,1:(end-1),3)'*a_avg_LR.phi;
            
            theta(n,2) = PN_theta(:,1:(end-1),1)'*a_avg_LR.theta;
            eta(n,2) = PN_eta(:,1:(end-1),1)'*a_avg_LR.eta;
            phi(n,2) = PN_phi(:,1:(end-1),1)'*a_avg_LR.phi;
    
            theta_dot(n,2) = PN_theta(:,1:(end-1),2)'*a_avg_LR.theta;
            eta_dot(n,2) = PN_eta(:,1:(end-1),2)'*a_avg_LR.eta;
            phi_dot(n,2) = PN_phi(:,1:(end-1),2)'*a_avg_LR.phi;

            theta_ddot(n,2) = PN_theta(:,1:(end-1),3)'*a_avg_LR.theta;
            eta_ddot(n,2) = PN_eta(:,1:(end-1),3)'*a_avg_LR.eta;
            phi_ddot(n,2) = PN_phi(:,1:(end-1),3)'*a_avg_LR.phi;
            
            down(m) = 1;
            
            down(l) = 0;
            
            up(l) = 1;
            
            up(m) = 0;
            
        elseif i == wb_nr
            
            n = (((i-1)*wb_sample_nr-(i-1))+1):(i*wb_sample_nr-i+1);
            
            m = (((i-1)*wb_sample_nr-(i-1))+1):(((i-1)*wb_sample_nr-(i-1))+down_str_l);
            
            l = (((i-1)*wb_sample_nr-(i-1))+down_str_l+1):(i*wb_sample_nr-i+1);
            
            theta(n,1) = PN_theta(:,:,1)'*a_avg_LR.theta;
            eta(n,1) = PN_eta(:,:,1)'*a_avg_LR.eta;
            phi(n,1) = PN_phi(:,:,1)'*a_avg_LR.phi;
    
            theta_dot(n,1) = PN_theta(:,:,2)'*a_avg_LR.theta;
            eta_dot(n,1) = PN_eta(:,:,2)'*a_avg_LR.eta;
            phi_dot(n,1) = PN_phi(:,:,2)'*a_avg_LR.phi;

            theta_ddot(n,1) = PN_theta(:,:,3)'*a_avg_LR.theta;
            eta_ddot(n,1) = PN_eta(:,:,3)'*a_avg_LR.eta;
            phi_ddot(n,1) = PN_phi(:,:,3)'*a_avg_LR.phi;
            
            theta(n,2) = PN_theta(:,:,1)'*a_avg_LR.theta;
            eta(n,2) = PN_eta(:,:,1)'*a_avg_LR.eta;
            phi(n,2) = PN_phi(:,:,1)'*a_avg_LR.phi;
    
            theta_dot(n,2) = PN_theta(:,:,2)'*a_avg_LR.theta;
            eta_dot(n,2) = PN_eta(:,:,2)'*a_avg_LR.eta;
            phi_dot(n,2) = PN_phi(:,:,2)'*a_avg_LR.phi;

            theta_ddot(n,2) = PN_theta(:,:,3)'*a_avg_LR.theta;
            eta_ddot(n,2) = PN_eta(:,:,3)'*a_avg_LR.eta;
            phi_ddot(n,2) = PN_phi(:,:,3)'*a_avg_LR.phi;            

            down(m) = 1;
            
            down(l) = 0;
            
            up(l) = 1;
            
            up(m) = 0;
            
        end
        
        
        
    end
    
    
    Pattern = {};
    
    Pattern.theta = theta;
    Pattern.eta = eta;
    Pattern.phi = phi;
    
    Pattern.theta_dot = theta_dot;
    Pattern.eta_dot = eta_dot;
    Pattern.phi_dot = phi_dot;
    
    Pattern.theta_ddot = theta_ddot;
    Pattern.eta_ddot = eta_ddot;
    Pattern.phi_ddot = phi_ddot;
    
    Pattern.down = down;
    Pattern.up = up;
    
    
%     figure()
%     plot(theta(:,1))
%     
%     figure()
%     plot(eta(:,1))
%     
%     figure()
%     plot(phi(:,1))
%     
%     figure()
%     plot(theta_dot(:,1))
%     
%     figure()
%     plot(eta_dot(:,1))
%     
%     figure()
%     plot(phi_dot(:,1))
%     
%     figure()
%     plot(theta_ddot(:,1))
%     
%     figure()
%     plot(eta_ddot(:,1))
%     
%     figure()
%     plot(phi_ddot(:,1))
%     
%     pause

    
    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Compute the body parameters:
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    wing_l = pathDB.wing_l(seq_nr);
    
    m_wl = 1e-6*0.1078*wing_l^3.008;
    
    nr_sect = 20;
    
    [cg_body, cg_L, cg_R, I_body, I_wing, I_v_wing, m_b, m_w, m_v_w, y_sect, chords] = cg_plus_Inertia(settings, pathDB, m_wl, wing_l ,seq_nr ,nr_sect);
    
    
    Body = {};
    
    Body.m_wl = m_wl;
    
    Body.m_w = m_w;
    
    Body.m_v_w = m_v_w;
    
    Body.wing_l = wing_l;
    
    Body.joint_L = pathDB.joint_pos_L(:,seq_nr);
    
    Body.joint_R = pathDB.joint_pos_R(:,seq_nr);
    
    Body.cg_body = cg_body;
    
    Body.cg_L = cg_L;
    
    Body.cg_R = cg_R;
    
    Body.I_body = I_body;
    
    Body.I_wing = I_wing;
    
    Body.I_v_wing = I_v_wing;
    
    Body.y_sect = y_sect;
    
    Body.chords = chords;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     q_L_test = zeros(3,3,nr_samples);
%     
%     q_R_test = zeros(3,3,nr_samples);
% 
%     
%     t = 0;
%     
%     t_end = ((nr_samples-1)/2)*dt;
%     
%     k = 1;
%     
%     while t <= t_end
%     
%     [q_L,w_L,w_dot_L,q_R,w_R,w_dot_R,up,down] = Wing_Pattern_Generator(Pattern,t, dt);
%     
%     q_L_test(:,:,k) = q_L;
%     
%     q_R_test(:,:,k) = q_R;
%     
%     t = t + 0.5*dt
%     
%     k = k+1
%     
%     end    
%     
%     L_wt_test = zeros(3,length(q_L_test(1,1,:)));
%     
%     R_wt_test = zeros(3,length(q_R_test(1,1,:)));
%     
%     L_wt_test2 = zeros(3,length(q_L_test(1,1,:)));
%     
%     R_wt_test2 = zeros(3,length(q_R_test(1,1,:)));
%     
%     for i = 1:length(q_L_test(1,1,:))
%         
%         L_wt_test(:,i) = q_L_test(:,:,i)'*[0; -1; 0];
%         
%         R_wt_test(:,i) = q_R_test(:,:,i)'*[0; 1; 0];
% 
%         L_wt_test2(:,i) = q_L_test(:,:,i)'*[0.05; -sqrt(1-0.05^2); 0];
%         
%         R_wt_test2(:,i) = q_R_test(:,:,i)'*[0.05; sqrt(1-0.05^2); 0];
%         
%     end
%     
%     b = 5;
%     c = 2^b-1;
%     [x,y,z] = sphere(c);
%     
% %     figure()
% %     plot3(L_wt_test(1,:), L_wt_test(2,:), L_wt_test(3,:))
% %     hold on
% %     plot3([0 L_wt_test(1,1)], [0 L_wt_test(2,1)], [0 L_wt_test(3,1)])
% %     grid on
% %     hold off
%     
%     fig_nr1 = 1;
%     
%     figure(fig_nr1)
%     surf(x,y,z,'FaceColor','black','EdgeColor','none');
%     alpha(0.2)
%     axis equal
%     hold on
%     
%     figure(fig_nr1)
%     plot3(L_wt_test(1,:), L_wt_test(2,:), L_wt_test(3,:),'r')
%     axis equal
%     hold on
%     
%     figure(fig_nr1)
%     plot3(R_wt_test(1,:), R_wt_test(2,:), R_wt_test(3,:),'g')
%     axis equal   
%     
%     figure(fig_nr1)
%     hold on
%     for l = 1:length(q_L_test(1,1,:))
%         plot3([L_wt_test(1,l) L_wt_test2(1,l)], [L_wt_test(2,l) L_wt_test2(2,l)], [L_wt_test(3,l) L_wt_test2(3,l)],'k')
%     end
%     
%     figure(fig_nr1)
%     hold on
%     for l = 1:length(q_L_test(1,1,:))
%         plot3([R_wt_test(1,l) R_wt_test2(1,l)], [R_wt_test(2,l) R_wt_test2(2,l)], [R_wt_test(3,l) R_wt_test2(3,l)],'k')
%     end
%     hold off

    
    % Test Pattern generator:
    
    q_L_test = zeros(4,nr_samples);
    
    q_R_test = zeros(4,nr_samples);
    
    w_L_test = zeros(3,nr_samples);
    
    w_R_test = zeros(3,nr_samples);
    
    w_dot_L_test = zeros(3,nr_samples);
    
    w_dot_R_test = zeros(3,nr_samples);
    
    t = 0;
    
    t_end = ((nr_samples-1)/2)*dt;
    
    k = 1;
    
    while t <= t_end
    
    [q_L,w_L,w_dot_L,q_R,w_R,w_dot_R,up,down] = Wing_Pattern_Generator(Pattern,t, dt);
    
    q_L_test(:,k) = q_L;
    
    q_R_test(:,k) = q_R;
    
    w_L_test(:,k) = w_L;
    
    w_R_test(:,k) = w_R;
    
    w_dot_L_test(:,k) = w_dot_L;
    
    w_dot_R_test(:,k) = w_dot_R;
    
    t = t + 0.5*dt
    
    k = k+1
    
    end
    

    

    
    figure()
    subplot(4,1,1); plot(q_L_test(1,:));
    subplot(4,1,2); plot(q_L_test(2,:));
    subplot(4,1,3); plot(q_L_test(3,:));
    subplot(4,1,4); plot(q_L_test(4,:));
    
    figure()
    subplot(4,1,1); plot(q_R_test(1,:));
    subplot(4,1,2); plot(q_R_test(2,:));
    subplot(4,1,3); plot(q_R_test(3,:));
    subplot(4,1,4); plot(q_R_test(4,:));
    
    figure()
    subplot(3,1,1); plot(w_L_test(1,:));
    subplot(3,1,2); plot(w_L_test(2,:));
    subplot(3,1,3); plot(w_L_test(3,:));
    
    figure()
    subplot(3,1,1); plot(w_R_test(1,:));
    subplot(3,1,2); plot(w_R_test(2,:));
    subplot(3,1,3); plot(w_R_test(3,:));
    
    figure()
    subplot(3,1,1); plot(w_dot_L_test(1,:));
    subplot(3,1,2); plot(w_dot_L_test(2,:));
    subplot(3,1,3); plot(w_dot_L_test(3,:));
    
    figure()
    subplot(3,1,1); plot(w_dot_R_test(1,:));
    subplot(3,1,2); plot(w_dot_R_test(2,:));
    subplot(3,1,3); plot(w_dot_R_test(3,:));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Run the state-space system in combination with a Runge-Kutta update
    % for the duration of the wingbeat pattern:
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Initial conditions;
    
%     t = 0;
%     
%     t_end = ((nr_samples-1)/2)*dt;
%     
%     sin((pi*(55/180))/2)
%     
%     cos((pi*(55/180))/2)
%     
%     State = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; sin((pi*(55/180))/2); 0; cos((pi*(55/180))/2); 0; 0; 0; 0; 0; 0];
%     
%     
%     while t <= t_end
%         
%         State
%         
%         [ State, Force_param ] = State_space(Pattern,Body,State,t,dt);
%         
% 
%         
%         [ State ] = Runge_Kutta_update(Pattern,Body,State,t,dt);
%         
% 
%         
%         t = t + dt;
%         
%     end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Plot the output:
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

