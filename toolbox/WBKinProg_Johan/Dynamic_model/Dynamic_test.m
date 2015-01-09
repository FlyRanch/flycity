function Dynamic_test(settings, pathDB)

    % Program that tests the dynamics of a fruit fly model for a given
    % sequence of wing kinematics.
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Create the wingbeat pattern:
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    strokeplane_WBkin = settings.strokeplane_WBkin;
    beta = (strokeplane_WBkin/180)*pi;

    % Get an average wingbeat:
    
    seq_nr = 4;
    
%     n_pol_theta = 13; % Order of used polynomials
%     n_pol_eta = 16; % Order of used polynomials
%     n_pol_phi = 8; % Order of used polynomials
% 
    n_pol_theta = 12; % Order of used polynomials
    n_pol_eta = 14; % Order of used polynomials
    n_pol_phi = 8; % Order of used polynomials

% 
%     n_pol_theta = 10; % Order of used polynomials
%     n_pol_eta = 12; % Order of used polynomials
%     n_pol_phi = 6; % Order of used polynomials
    
    
    [a_fit,a_avg,f_avg,down_up,trigger_wb,ratio_1_2,ratio_1_2_avg] = Standard_wingbeat( settings, pathDB, seq_nr, n_pol_theta, n_pol_eta, n_pol_phi);

%     a_fit
%     
%     a_avg
%     
%     f_avg
%     
%     down_up
%     
%     trigger_wb
%     
%     ratio_1_2
%     
%     ratio_1_2_avg
    
    wb_time = 1/f_avg;
    
    wb_sample_nr = 51; % Always take an odd number
    
    dt = wb_time/((wb_sample_nr-1)/2);
    
    sample_loc1 = -1:(2/round((wb_sample_nr-1)*ratio_1_2_avg/2)):1;
    
    sample_loc2 = -1:(2/(round((wb_sample_nr-1)*(2-ratio_1_2_avg)/2))):1;
    
    sample_nr1 = length(sample_loc1);
    
    sample_nr2 = length(sample_loc2);
    
    wb_nr = 40; % Number of simulated wingbeats
    
    nr_samples = (wb_sample_nr*wb_nr-(wb_nr-1)); % Total number of time steps
        
    PN_theta_1 = Legendre_polynomial(n_pol_theta,3,sample_loc1);
    PN_eta_1 = Legendre_polynomial(n_pol_eta,3,sample_loc1);
    PN_phi_1 = Legendre_polynomial(n_pol_phi,3,sample_loc1);  
    
    PN_theta_2 = Legendre_polynomial(n_pol_theta,3,sample_loc2);
    PN_eta_2 = Legendre_polynomial(n_pol_eta,3,sample_loc2);
    PN_phi_2 = Legendre_polynomial(n_pol_phi,3,sample_loc2);  
    
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
            
            n1 = (((i-1)*wb_sample_nr-(i-1))+1):(((i-1)*wb_sample_nr-(i-1))+sample_nr1-1);
            
            n2 = (((i-1)*wb_sample_nr-(i-1))+sample_nr1):(i*wb_sample_nr-i);
            
            down_str_l = round(down_up*(wb_sample_nr-1)/(1+down_up));
            
            m = (((i-1)*wb_sample_nr-(i-1))+1):(((i-1)*wb_sample_nr-(i-1))+down_str_l);
            
            l = (((i-1)*wb_sample_nr-(i-1))+down_str_l+1):(i*wb_sample_nr-i);
            
            theta(n1,1) = PN_theta_1(:,1:(end-1),1)'*a_avg.theta_LR1;
            eta(n1,1) = PN_eta_1(:,1:(end-1),1)'*a_avg.eta_LR1;
            phi(n1,1) = PN_phi_1(:,1:(end-1),1)'*a_avg.phi_LR1;
    
            theta_dot(n1,1) = PN_theta_1(:,1:(end-1),2)'*a_avg.theta_LR1;
            eta_dot(n1,1) = PN_eta_1(:,1:(end-1),2)'*a_avg.eta_LR1;
            phi_dot(n1,1) = PN_phi_1(:,1:(end-1),2)'*a_avg.phi_LR1;

            theta_ddot(n1,1) = PN_theta_1(:,1:(end-1),3)'*a_avg.theta_LR1;
            eta_ddot(n1,1) = PN_eta_1(:,1:(end-1),3)'*a_avg.eta_LR1;
            phi_ddot(n1,1) = PN_phi_1(:,1:(end-1),3)'*a_avg.phi_LR1;
            
            theta(n1,2) = PN_theta_1(:,1:(end-1),1)'*a_avg.theta_LR1;
            eta(n1,2) = PN_eta_1(:,1:(end-1),1)'*a_avg.eta_LR1;
            phi(n1,2) = PN_phi_1(:,1:(end-1),1)'*a_avg.phi_LR1;
    
            theta_dot(n1,2) = PN_theta_1(:,1:(end-1),2)'*a_avg.theta_LR1;
            eta_dot(n1,2) = PN_eta_1(:,1:(end-1),2)'*a_avg.eta_LR1;
            phi_dot(n1,2) = PN_phi_1(:,1:(end-1),2)'*a_avg.phi_LR1;

            theta_ddot(n1,2) = PN_theta_1(:,1:(end-1),3)'*a_avg.theta_LR1;
            eta_ddot(n1,2) = PN_eta_1(:,1:(end-1),3)'*a_avg.eta_LR1;
            phi_ddot(n1,2) = PN_phi_1(:,1:(end-1),3)'*a_avg.phi_LR1;
            
            
            
            
            theta(n2,1) = PN_theta_2(:,1:(end-1),1)'*a_avg.theta_LR2;
            eta(n2,1) = PN_eta_2(:,1:(end-1),1)'*a_avg.eta_LR2;
            phi(n2,1) = PN_phi_2(:,1:(end-1),1)'*a_avg.phi_LR2;
    
            theta_dot(n2,1) = PN_theta_2(:,1:(end-1),2)'*a_avg.theta_LR2;
            eta_dot(n2,1) = PN_eta_2(:,1:(end-1),2)'*a_avg.eta_LR2;
            phi_dot(n2,1) = PN_phi_2(:,1:(end-1),2)'*a_avg.phi_LR2;

            theta_ddot(n2,1) = PN_theta_2(:,1:(end-1),3)'*a_avg.theta_LR2;
            eta_ddot(n2,1) = PN_eta_2(:,1:(end-1),3)'*a_avg.eta_LR2;
            phi_ddot(n2,1) = PN_phi_2(:,1:(end-1),3)'*a_avg.phi_LR2;
            
            theta(n2,2) = PN_theta_2(:,1:(end-1),1)'*a_avg.theta_LR2;
            eta(n2,2) = PN_eta_2(:,1:(end-1),1)'*a_avg.eta_LR2;
            phi(n2,2) = PN_phi_2(:,1:(end-1),1)'*a_avg.phi_LR2;
    
            theta_dot(n2,2) = PN_theta_2(:,1:(end-1),2)'*a_avg.theta_LR2;
            eta_dot(n2,2) = PN_eta_2(:,1:(end-1),2)'*a_avg.eta_LR2;
            phi_dot(n2,2) = PN_phi_2(:,1:(end-1),2)'*a_avg.phi_LR2;

            theta_ddot(n2,2) = PN_theta_2(:,1:(end-1),3)'*a_avg.theta_LR2;
            eta_ddot(n2,2) = PN_eta_2(:,1:(end-1),3)'*a_avg.eta_LR2;
            phi_ddot(n2,2) = PN_phi_2(:,1:(end-1),3)'*a_avg.phi_LR2;            
            
            
            down(m) = 1;
            
            down(l) = 0;
            
            up(l) = 1;
            
            up(m) = 0;
            
        elseif i == wb_nr
            
            n1 = (((i-1)*wb_sample_nr-(i-1))+1):(((i-1)*wb_sample_nr-(i-1))+sample_nr1-1);
            
            n2 = (((i-1)*wb_sample_nr-(i-1))+sample_nr1):(i*wb_sample_nr-i+1);
            
            m = (((i-1)*wb_sample_nr-(i-1))+1):(((i-1)*wb_sample_nr-(i-1))+down_str_l);
            
            l = (((i-1)*wb_sample_nr-(i-1))+down_str_l+1):(i*wb_sample_nr-i+1);
            
            theta(n1,1) = PN_theta_1(:,1:(end-1),1)'*a_avg.theta_LR1;
            eta(n1,1) = PN_eta_1(:,1:(end-1),1)'*a_avg.eta_LR1;
            phi(n1,1) = PN_phi_1(:,1:(end-1),1)'*a_avg.phi_LR1;
    
            theta_dot(n1,1) = PN_theta_1(:,1:(end-1),2)'*a_avg.theta_LR1;
            eta_dot(n1,1) = PN_eta_1(:,1:(end-1),2)'*a_avg.eta_LR1;
            phi_dot(n1,1) = PN_phi_1(:,1:(end-1),2)'*a_avg.phi_LR1;

            theta_ddot(n1,1) = PN_theta_1(:,1:(end-1),3)'*a_avg.theta_LR1;
            eta_ddot(n1,1) = PN_eta_1(:,1:(end-1),3)'*a_avg.eta_LR1;
            phi_ddot(n1,1) = PN_phi_1(:,1:(end-1),3)'*a_avg.phi_LR1;
            
            theta(n1,2) = PN_theta_1(:,1:(end-1),1)'*a_avg.theta_LR1;
            eta(n1,2) = PN_eta_1(:,1:(end-1),1)'*a_avg.eta_LR1;
            phi(n1,2) = PN_phi_1(:,1:(end-1),1)'*a_avg.phi_LR1;
    
            theta_dot(n1,2) = PN_theta_1(:,1:(end-1),2)'*a_avg.theta_LR1;
            eta_dot(n1,2) = PN_eta_1(:,1:(end-1),2)'*a_avg.eta_LR1;
            phi_dot(n1,2) = PN_phi_1(:,1:(end-1),2)'*a_avg.phi_LR1;

            theta_ddot(n1,2) = PN_theta_1(:,1:(end-1),3)'*a_avg.theta_LR1;
            eta_ddot(n1,2) = PN_eta_1(:,1:(end-1),3)'*a_avg.eta_LR1;
            phi_ddot(n1,2) = PN_phi_1(:,1:(end-1),3)'*a_avg.phi_LR1;
            
            
            
            
            theta(n2,1) = PN_theta_2(:,:,1)'*a_avg.theta_LR2;
            eta(n2,1) = PN_eta_2(:,:,1)'*a_avg.eta_LR2;
            phi(n2,1) = PN_phi_2(:,:,1)'*a_avg.phi_LR2;
    
            theta_dot(n2,1) = PN_theta_2(:,:,2)'*a_avg.theta_LR2;
            eta_dot(n2,1) = PN_eta_2(:,:,2)'*a_avg.eta_LR2;
            phi_dot(n2,1) = PN_phi_2(:,:,2)'*a_avg.phi_LR2;

            theta_ddot(n2,1) = PN_theta_2(:,:,3)'*a_avg.theta_LR2;
            eta_ddot(n2,1) = PN_eta_2(:,:,3)'*a_avg.eta_LR2;
            phi_ddot(n2,1) = PN_phi_2(:,:,3)'*a_avg.phi_LR2;
            
            theta(n2,2) = PN_theta_2(:,:,1)'*a_avg.theta_LR2;
            eta(n2,2) = PN_eta_2(:,:,1)'*a_avg.eta_LR2;
            phi(n2,2) = PN_phi_2(:,:,1)'*a_avg.phi_LR2;
    
            theta_dot(n2,2) = PN_theta_2(:,:,2)'*a_avg.theta_LR2;
            eta_dot(n2,2) = PN_eta_2(:,:,2)'*a_avg.eta_LR2;
            phi_dot(n2,2) = PN_phi_2(:,:,2)'*a_avg.phi_LR2;

            theta_ddot(n2,2) = PN_theta_2(:,:,3)'*a_avg.theta_LR2;
            eta_ddot(n2,2) = PN_eta_2(:,:,3)'*a_avg.eta_LR2;
            phi_ddot(n2,2) = PN_phi_2(:,:,3)'*a_avg.phi_LR2;          

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
    
    
    figure()
    plot(theta(:,1),'r')
    hold on
    plot(theta(:,2),'b')
    hold off
    
    figure()
    plot(eta(:,1),'r')
    hold on
    plot(eta(:,2),'b')
    hold off
    
    figure()
    plot(phi(:,1),'r')
    hold on
    plot(phi(:,2),'b')
    hold off
    
    figure()
    plot(theta_dot(:,1),'r')
    hold on
    plot(theta_dot(:,2),'b')
    hold off
    
    
    figure()
    plot(eta_dot(:,1),'r')
    hold on
    plot(eta_dot(:,2),'b')
    hold off
    
    figure()
    plot(phi_dot(:,1),'r')
    hold on
    plot(phi_dot(:,2),'b')
    hold off
    
    figure()
    plot(theta_ddot(:,1),'r')
    hold on
    plot(theta_ddot(:,2),'b')
    hold off
    
    figure()
    plot(eta_ddot(:,1),'r')
    hold on
    plot(eta_ddot(:,2),'b')
    hold off
    
    figure()
    plot(phi_ddot(:,1),'r')
    hold on
    plot(phi_ddot(:,2),'b')
    hold off
    


    
    
        
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

    
%     % Test Pattern generator:
%     
%     q_L_test = zeros(4,nr_samples);
%     
%     q_R_test = zeros(4,nr_samples);
%     
%     w_L_test = zeros(3,nr_samples);
%     
%     w_R_test = zeros(3,nr_samples);
%     
%     w_dot_L_test = zeros(3,nr_samples);
%     
%     w_dot_R_test = zeros(3,nr_samples);
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
%     q_L_test(:,k) = q_L;
%     
%     q_R_test(:,k) = q_R;
% 
%    
%     w_L_test(:,k) = w_L;
%     
%     w_R_test(:,k) = w_R;
%     
%     w_dot_L_test(:,k) = w_dot_L;
%     
%     w_dot_R_test(:,k) = w_dot_R;
%     
%     t = t + 0.5*dt
%     
%     k = k+1
%     
%     end
%     
% %     q_L_test = [pathDB.qL1_filt2(2000:2200); pathDB.qL2_filt2(2000:2200); pathDB.qL3_filt2(2000:2200); pathDB.qL4_filt2(2000:2200)];
% %     
% %     q_R_test = [pathDB.qR1_filt2(2000:2200); pathDB.qR2_filt2(2000:2200); pathDB.qR3_filt2(2000:2200); pathDB.qR4_filt2(2000:2200)];
%     
%     L_wt_test = zeros(3,length(q_L_test(1,:)));
%     
%     R_wt_test = zeros(3,length(q_R_test(1,:)));
%     
%     L_wt_test2 = zeros(3,length(q_L_test(1,:)));
%     
%     R_wt_test2 = zeros(3,length(q_R_test(1,:)));
%     
%     for i = 1:length(q_L_test(1,:))
%         
%         L_wt_test(:,i) = quat2matNEW(q_L_test(:,i))*[0; -1; 0];
%         
%         R_wt_test(:,i) = quat2matNEW(q_R_test(:,i))*[0; 1; 0];
% 
%         L_wt_test2(:,i) = quat2matNEW(q_L_test(:,i))*[0.05; -sqrt(1-0.05^2); 0];
%         
%         R_wt_test2(:,i) = quat2matNEW(q_R_test(:,i))*[0.05; sqrt(1-0.05^2); 0];
%         
%     end
%     
%     
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
%     for l = 1:length(q_L_test(1,:))
%         plot3([L_wt_test(1,l) L_wt_test2(1,l)], [L_wt_test(2,l) L_wt_test2(2,l)], [L_wt_test(3,l) L_wt_test2(3,l)],'k')
%     end
%     
%     figure(fig_nr1)
%     hold on
%     for l = 1:length(q_L_test(1,:))
%         plot3([R_wt_test(1,l) R_wt_test2(1,l)], [R_wt_test(2,l) R_wt_test2(2,l)], [R_wt_test(3,l) R_wt_test2(3,l)],'k')
%     end
%     hold off
%     
%     
%     
%     
% %     
% % 
% %     
% % 
% %     
%     figure()
%     subplot(4,1,1); plot(q_L_test(1,:));
%     subplot(4,1,2); plot(q_L_test(2,:));
%     subplot(4,1,3); plot(q_L_test(3,:));
%     subplot(4,1,4); plot(q_L_test(4,:));
%     
%     figure()
%     subplot(4,1,1); plot(q_R_test(1,:));
%     subplot(4,1,2); plot(q_R_test(2,:));
%     subplot(4,1,3); plot(q_R_test(3,:));
%     subplot(4,1,4); plot(q_R_test(4,:));
%     
%     figure()
%     subplot(3,1,1); plot(w_L_test(1,:));
%     subplot(3,1,2); plot(w_L_test(2,:));
%     subplot(3,1,3); plot(w_L_test(3,:));
%     
%     figure()
%     subplot(3,1,1); plot(w_R_test(1,:));
%     subplot(3,1,2); plot(w_R_test(2,:));
%     subplot(3,1,3); plot(w_R_test(3,:));
%     
%     figure()
%     subplot(3,1,1); plot(w_dot_L_test(1,:));
%     subplot(3,1,2); plot(w_dot_L_test(2,:));
%     subplot(3,1,3); plot(w_dot_L_test(3,:));
%     
%     figure()
%     subplot(3,1,1); plot(w_dot_R_test(1,:));
%     subplot(3,1,2); plot(w_dot_R_test(2,:));
%     subplot(3,1,3); plot(w_dot_R_test(3,:));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Run the state-space system in combination with a Runge-Kutta update
    % for the duration of the wingbeat pattern:
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Initial conditions;
    
    t = 0;
    
    time = 0;
    
    t_end = ((nr_samples-3)/2)*dt;
    
%     R_1 = [cos(pi*(55/180)) 0 -sin(pi*(55/180)); ...
%            0 1 0; ...
%            sin(pi*(55/180)) 0 cos(pi*(55/180))];
           
    R_1 = [cos(beta) 0 sin(beta); ...
           0 1 0; ...
           -sin(beta) 0 cos(beta)];
           
    R_2 = [1 0 0; ...
           0 cos(pi) sin(pi); ...
           0 -sin(pi) cos(pi)];
    
    q_0 = quat2matNEW((R_1*R_2)');
             
    State = [0; 0; 0; 0; 0; 0; 0; 0; 0; q_0(1); q_0(2); q_0(3); q_0(4); 0; 0; 0; 0; 0; 0];
    
    State_vector = State;
    
    Force_vector = {};
    
    Force_vector.F_I = zeros(3,1);
        
    Force_vector.M_I = zeros(3,1);
        
    Force_vector.F_A = zeros(3,1);
        
    Force_vector.M_A = zeros(3,1);
        
    Force_vector.F_g = zeros(3,1);
        
    Force_vector.M_g = zeros(3,1);
        
    Force_vector.F_r = zeros(3,1);
        
    Force_vector.M_r = zeros(3,1);
        
    Force_vector.F_aero_L = zeros(3,1);
        
    Force_vector.F_aero_R = zeros(3,1);
        
    Force_vector.F_aero_L_arm = zeros(3,1);
        
    Force_vector.F_aero_R_arm = zeros(3,1);
    
    Force_vector.KE = 0;
    
    Force_vector.KE_trans = 0;
    
    Force_vector.KE_rot = 0;
    
    
    
    Wing_kin = {};
    
    Wing_kin.q_L = zeros(4,1);
    
    Wing_kin.q_R = zeros(4,1);
    
    Wing_kin.w_L = zeros(3,1);
    
    Wing_kin.w_R = zeros(3,1);
    
    Wing_kin.w_dot_L = zeros(3,1);
    
    Wing_kin.w_dot_R = zeros(3,1);
    
    
    
    while t <= t_end
        
        %State
        
        [ State, Force_param ] = State_space(Pattern,Body,State,t,dt);
        
        Force_vector.F_I = [Force_vector.F_I Force_param.F_I];
        
        Force_vector.M_I = [Force_vector.M_I Force_param.M_I];
        
        Force_vector.F_A = [Force_vector.F_A Force_param.F_A];
        
        Force_vector.M_A = [Force_vector.M_A Force_param.M_A];
        
        Force_vector.F_g = [Force_vector.F_g Force_param.F_g];
        
        Force_vector.M_g = [Force_vector.M_g Force_param.M_g];
        
        Force_vector.F_r = [Force_vector.F_r Force_param.F_R];
        
        Force_vector.M_r = [Force_vector.M_r Force_param.M_R];
        
        Force_vector.F_aero_L = [Force_vector.F_aero_L Force_param.F_aero_L];
        
        Force_vector.F_aero_R = [Force_vector.F_aero_R Force_param.F_aero_R];
        
        Force_vector.F_aero_L_arm = [Force_vector.F_aero_L_arm Force_param.F_aero_L_arm];
        
        Force_vector.F_aero_R_arm = [Force_vector.F_aero_R_arm Force_param.F_aero_R_arm];
        
        Force_vector.KE = [Force_vector.KE Force_param.KE];
        
        Force_vector.KE_trans = [Force_vector.KE_trans Force_param.KE_trans];
        
        Force_vector.KE_rot = [Force_vector.KE_rot Force_param.KE_rot];
        
        Wing_kin.q_L = [Wing_kin.q_L Force_param.q_L];

        Wing_kin.q_R = [Wing_kin.q_R Force_param.q_R];

        Wing_kin.w_L = [Wing_kin.w_L Force_param.w_L];

        Wing_kin.w_R = [Wing_kin.w_R Force_param.w_R];

        Wing_kin.w_dot_L = [Wing_kin.w_dot_L Force_param.w_dot_L];

        Wing_kin.w_dot_R = [Wing_kin.w_dot_R Force_param.w_dot_R];
        
        [ State ] = Runge_Kutta_update(Pattern,Body,State,t,dt);
        
        State_vector = [State_vector State];
        
        time = [time t+dt];
        
        t
        
        t = t + dt;
        
        
        
        
    end
    
    State_vector
    
    Force_vector
    
    Wing_kin
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Plot the output:
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure()
    plot3(State_vector(1,:), State_vector(2,:), State_vector(3,:))
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    axis('equal')
    title('Trajectory')
    
    figure()
    subplot(3,1,1); plot(time,State_vector(4,:))
    xlabel('time [s]')
    ylabel('u [mm/s]')
    title('Velocity')
    hold on
    subplot(3,1,2); plot(time,State_vector(5,:))
    xlabel('time [s]')
    ylabel('v [mm/s]')
    subplot(3,1,3); plot(time,State_vector(6,:))
    xlabel('time [s]')
    ylabel('w [mm/s]')
    hold off
    
    figure()
    subplot(3,1,1); plot(time,State_vector(7,:))
    xlabel('time [s]')
    ylabel('ax [mm/s^2]')
    title('Acceleration')
    hold on
    subplot(3,1,2); plot(time,State_vector(8,:))
    xlabel('time [s]')
    ylabel('ay [mm/s^2]')
    subplot(3,1,3); plot(time,State_vector(9,:))
    xlabel('time [s]')
    ylabel('az [mm/s^2]')
    hold off    
    
    figure()
    subplot(4,1,1); plot(time,State_vector(10,:))
    xlabel('time [s]')
    ylabel('q1')
    title('Body quaternion')
    hold on
    subplot(4,1,2); plot(time,State_vector(11,:))
    xlabel('time [s]')
    ylabel('q2')
    subplot(4,1,3); plot(time,State_vector(12,:))
    xlabel('time [s]')
    ylabel('q3')
    subplot(4,1,4); plot(time,State_vector(13,:))
    xlabel('time [s]')
    ylabel('q4')
    hold off    

    figure()
    subplot(3,1,1); plot(time,State_vector(14,:))
    xlabel('time [s]')
    ylabel('wx [rad/s]')
    title('Angular velocity body')
    hold on
    subplot(3,1,2); plot(time,State_vector(15,:))
    xlabel('time [s]')
    ylabel('wy [rad/s]')
    subplot(3,1,3); plot(time,State_vector(16,:))
    xlabel('time [s]')
    ylabel('wz [rad/s]')
    hold off    
    
    figure()
    subplot(3,1,1); plot(time,State_vector(17,:))
    xlabel('time [s]')
    ylabel('w dot x [rad/s^2]')
    title('Angular acceleration body')
    hold on
    subplot(3,1,2); plot(time,State_vector(18,:))
    xlabel('time [s]')
    ylabel('w dot y [rad/s^2]')
    subplot(3,1,3); plot(time,State_vector(19,:))
    xlabel('time [s]')
    ylabel('w dot z [rad/s^2]')
    hold off      
    
    figure()
    subplot(3,1,1); plot(time,Force_vector.F_I(1,:),'b',time,mean(Force_vector.F_I(1,:)),'r')
    xlabel('time [s]')
    ylabel('Fx [N]')
    title('Inertial forces body')
    hold on
    subplot(3,1,2); plot(time,Force_vector.F_I(2,:),'b',time,mean(Force_vector.F_I(2,:)),'r')
    xlabel('time [s]')
    ylabel('Fy [N]')
    subplot(3,1,3); plot(time,Force_vector.F_I(3,:),'b',time,mean(Force_vector.F_I(3,:)),'r')
    xlabel('time [s]')
    ylabel('Fz [N]')
    hold off        

    
    figure()
    subplot(3,1,1); plot(time,Force_vector.M_I(1,:),'b',time,mean(Force_vector.M_I(1,:)),'r')
    xlabel('time [s]')
    ylabel('Mx [N*mm]')
    title('Inertial moments body')
    hold on
    subplot(3,1,2); plot(time,Force_vector.M_I(2,:),'b',time,mean(Force_vector.M_I(2,:)),'r')
    xlabel('time [s]')
    ylabel('My [N*mm]')
    subplot(3,1,3); plot(time,Force_vector.M_I(3,:),'b',time,mean(Force_vector.M_I(3,:)),'r')
    xlabel('time [s]')
    ylabel('Mz [N*mm]')
    hold off        
    

    figure()
    subplot(3,1,1); plot(time,Force_vector.F_A(1,:),'b',time,mean(Force_vector.F_A(1,:)),'r')
    xlabel('time [s]')
    ylabel('Fx [N]')
    title('Aerodynamic forces body')
    hold on
    subplot(3,1,2); plot(time,Force_vector.F_A(2,:),'b',time,mean(Force_vector.F_A(2,:)),'r')
    xlabel('time [s]')
    ylabel('Fy [N]')
    subplot(3,1,3); plot(time,Force_vector.F_A(3,:),'b',time,mean(Force_vector.F_A(3,:)),'r')
    xlabel('time [s]')
    ylabel('Fz [N]')
    hold off        

    
    figure()
    subplot(3,1,1); plot(time,Force_vector.M_A(1,:),'b',time,mean(Force_vector.M_A(1,:)),'r')
    xlabel('time [s]')
    ylabel('Mx [N*mm]')
    title('Aerodynamic moments body')
    hold on
    subplot(3,1,2); plot(time,Force_vector.M_A(2,:),'b',time,mean(Force_vector.M_A(2,:)),'r')
    xlabel('time [s]')
    ylabel('My [N*mm]')
    subplot(3,1,3); plot(time,Force_vector.M_A(3,:),'b',time,mean(Force_vector.M_A(3,:)),'r')
    xlabel('time [s]')
    ylabel('Mz [N*mm]')
    hold off 
    
    figure()
    subplot(3,1,1); plot(time,Force_vector.F_aero_L(1,:),'b',time,mean(Force_vector.F_aero_L(1,:)),'r')
    xlabel('time [s]')
    ylabel('Fx [N]')
    title('Aerodynamic forces left wing')
    hold on
    subplot(3,1,2); plot(time,Force_vector.F_aero_L(2,:),'b',time,mean(Force_vector.F_aero_L(2,:)),'r')
    xlabel('time [s]')
    ylabel('Fy [N]')
    subplot(3,1,3); plot(time,Force_vector.F_aero_L(3,:),'b',time,mean(Force_vector.F_aero_L(3,:)),'r')
    xlabel('time [s]')
    ylabel('Fz [N]')
    hold off        

  
    figure()
    subplot(3,1,1); plot(time,Force_vector.F_aero_R(1,:),'b',time,mean(Force_vector.F_aero_R(1,:)),'r')
    xlabel('time [s]')
    ylabel('Fx [N]')
    title('Aerodynamic forces right wing')
    hold on
    subplot(3,1,2); plot(time,Force_vector.F_aero_R(2,:),'b',time,mean(Force_vector.F_aero_R(2,:)),'r')
    xlabel('time [s]')
    ylabel('Fy [N]')
    subplot(3,1,3); plot(time,Force_vector.F_aero_R(3,:),'b',time,mean(Force_vector.F_aero_R(3,:)),'r')
    xlabel('time [s]')
    ylabel('Fz [N]')
    hold off        

    
%     figure()
%     subplot(3,1,1); plot(time,Force_vector.F_aero_L_arm(1,:))
%     xlabel('time [s]')
%     ylabel('x [mm]')
%     title('Aerodynamic moment arm left wing')
%     hold on
%     subplot(3,1,2); plot(time,Force_vector.F_aero_L_arm(2,:))
%     xlabel('time [s]')
%     ylabel('y [mm]')
%     subplot(3,1,3); plot(time,Force_vector.F_aero_L_arm(3,:))
%     xlabel('time [s]')
%     ylabel('z [mm]')
%     hold off        
% 
%   
%     figure()
%     subplot(3,1,1); plot(time,Force_vector.F_aero_R_arm(1,:))
%     xlabel('time [s]')
%     ylabel('x [mm]')
%     title('Aerodynamic moment arm right wing')
%     hold on
%     subplot(3,1,2); plot(time,Force_vector.F_aero_R_arm(2,:))
%     xlabel('time [s]')
%     ylabel('y [mm]')
%     subplot(3,1,3); plot(time,Force_vector.F_aero_R_arm(3,:))
%     xlabel('time [s]')
%     ylabel('z [mm]')
%     hold off   
    
    
    figure()
    subplot(3,1,1); plot(time,Force_vector.F_g(1,:),'b',time,mean(Force_vector.F_g(1,:)),'r')
    xlabel('time [s]')
    ylabel('Fx [N]')
    title('Gravitational forces body')
    hold on
    subplot(3,1,2); plot(time,Force_vector.F_g(2,:),'b',time,mean(Force_vector.F_g(2,:)),'r')
    xlabel('time [s]')
    ylabel('Fy [N]')
    subplot(3,1,3); plot(time,Force_vector.F_g(3,:),'b',time,mean(Force_vector.F_g(3,:)),'r')
    xlabel('time [s]')
    ylabel('Fz [N]')
    hold off        

    
    figure()
    subplot(3,1,1); plot(time,Force_vector.M_g(1,:))
    xlabel('time [s]')
    ylabel('Mx [N*mm]')
    title('Gravitational moments body')
    hold on
    subplot(3,1,2); plot(time,Force_vector.M_g(2,:))
    xlabel('time [s]')
    ylabel('My [N*mm]')
    subplot(3,1,3); plot(time,Force_vector.M_g(3,:))
    xlabel('time [s]')
    ylabel('Mz [N*mm]')
    hold off 
    
    
    figure()
    subplot(3,1,1); plot(time,Force_vector.F_r(1,:),'b',time,mean(Force_vector.F_r(1,:)),'r')
    xlabel('time [s]')
    ylabel('Fx [N]')
    title('Resultant forces body')
    hold on
    subplot(3,1,2); plot(time,Force_vector.F_r(2,:),'b',time,mean(Force_vector.F_r(2,:)),'r')
    xlabel('time [s]')
    ylabel('Fy [N]')
    subplot(3,1,3); plot(time,Force_vector.F_r(3,:),'b',time,mean(Force_vector.F_r(3,:)),'r')
    xlabel('time [s]')
    ylabel('Fz [N]')
    hold off        

    
    figure()
    subplot(3,1,1); plot(time,Force_vector.M_r(1,:),'b',time,mean(Force_vector.M_r(1,:)),'r')
    xlabel('time [s]')
    ylabel('Mx [N*mm]')
    title('Resultant moments body')
    hold on
    subplot(3,1,2); plot(time,Force_vector.M_r(2,:),'b',time,mean(Force_vector.M_r(2,:)),'r')
    xlabel('time [s]')
    ylabel('My [N*mm]')
    subplot(3,1,3); plot(time,Force_vector.M_r(3,:),'b',time,mean(Force_vector.M_r(3,:)),'r')
    xlabel('time [s]')
    ylabel('Mz [N*mm]')
    hold off    
    
    
    figure()
    plot(time,Force_vector.KE)
    xlabel('time [s]')
    ylabel('Kinetic energy [N*mm]')    
    title('Kinetic energy body and wings')
    
    figure()
    plot(time,Force_vector.KE_trans)
    xlabel('time [s]')
    ylabel('Kinetic energy [N*mm]')    
    title('Translational Kinetic energy body and wings')    
    
    figure()
    plot(time,Force_vector.KE_rot)
    xlabel('time [s]')
    ylabel('Kinetic energy [N*mm]')    
    title('Rotational Kinetic energy body and wings')  
    
    
    figure()
    subplot(4,1,1); plot(time,Wing_kin.q_L(1,:))
    xlabel('time [s]')
    ylabel('q(1)')
    title('Left wing quaternion')
    hold on
    subplot(4,1,2); plot(time,Wing_kin.q_L(2,:))
    xlabel('time [s]')
    ylabel('q(2)')
    subplot(4,1,3); plot(time,Wing_kin.q_L(3,:))
    xlabel('time [s]')
    ylabel('q(3)')
    subplot(4,1,4); plot(time,Wing_kin.q_L(4,:))
    xlabel('time [s]')
    ylabel('q(4)')
    hold off

    figure()
    subplot(4,1,1); plot(time,Wing_kin.q_R(1,:))
    xlabel('time [s]')
    ylabel('q(1)')
    title('Right wing quaternion')
    hold on
    subplot(4,1,2); plot(time,Wing_kin.q_R(2,:))
    xlabel('time [s]')
    ylabel('q(2)')
    subplot(4,1,3); plot(time,Wing_kin.q_R(3,:))
    xlabel('time [s]')
    ylabel('q(3)')
    subplot(4,1,4); plot(time,Wing_kin.q_R(4,:))
    xlabel('time [s]')
    ylabel('q(4)')
    hold off
    
    figure()
    subplot(3,1,1); plot(time,Wing_kin.w_L(1,:))
    xlabel('time [s]')
    ylabel('\omega_x')
    title('Angular velocity left wing')
    hold on
    subplot(3,1,2); plot(time,Wing_kin.w_L(2,:))
    xlabel('time [s]')
    ylabel('\omega_y')
    subplot(3,1,3); plot(time,Wing_kin.w_L(3,:))
    xlabel('time [s]')
    ylabel('\omega_z')
    hold off
    
    figure()
    subplot(3,1,1); plot(time,Wing_kin.w_R(1,:))
    xlabel('time [s]')
    ylabel('\omega_x')
    title('Angular velocity right wing')
    hold on
    subplot(3,1,2); plot(time,Wing_kin.w_R(2,:))
    xlabel('time [s]')
    ylabel('\omega_y')
    subplot(3,1,3); plot(time,Wing_kin.w_R(3,:))
    xlabel('time [s]')
    ylabel('\omega_z')
    hold off
    
    figure()
    subplot(3,1,1); plot(time,Wing_kin.w_dot_L(1,:))
    xlabel('time [s]')
    ylabel('w_dot_x')
    title('Angular acceleration left wing')
    hold on
    subplot(3,1,2); plot(time,Wing_kin.w_dot_L(2,:))
    xlabel('time [s]')
    ylabel('w_dot_y')
    subplot(3,1,3); plot(time,Wing_kin.w_dot_L(3,:))
    xlabel('time [s]')
    ylabel('w_dot_z')
    hold off
    
    figure()
    subplot(3,1,1); plot(time,Wing_kin.w_dot_R(1,:))
    xlabel('time [s]')
    ylabel('w_dot_x')
    title('Angular acceleration right wing')
    hold on
    subplot(3,1,2); plot(time,Wing_kin.w_dot_R(2,:))
    xlabel('time [s]')
    ylabel('w_dot_y')
    subplot(3,1,3); plot(time,Wing_kin.w_dot_R(3,:))
    xlabel('time [s]')
    ylabel('w_dot_z')
    hold off   
    
    % Plot body euler angles:

    phi_body = zeros(1,length(State_vector(10,:)));
    eta_body = zeros(1,length(State_vector(10,:)));
    xsi_body = zeros(1,length(State_vector(10,:)));
    
    for k = 1:length(State_vector(10,:))
        
        q_t1 = State_vector(10,k);
        q_t2 = State_vector(11,k);
        q_t3 = State_vector(12,k);
        q_t4 = State_vector(13,k);
        
        phi_body(k) = atan2(2*(q_t4*q_t1+q_t3*q_t2),1-2*(q_t1^2+q_t2^2));
        eta_body(k) = asin(2*(q_t4*q_t2-q_t3*q_t1));
        xsi_body(k) = atan2(2*(q_t4*q_t3+q_t1*q_t2),1-2*(q_t2^2+q_t3^2));
        
    end
    
    figure()
    subplot(3,1,1); plot(time,radtodeg(phi_body)-180)
    xlabel('time [s]')
    ylabel('phi body')
    title('Euler angles body')
    hold on
    subplot(3,1,2); plot(time,radtodeg(eta_body))
    xlabel('time [s]')
    ylabel('eta body')
    subplot(3,1,3); plot(time,radtodeg(xsi_body))
    xlabel('time [s]')
    ylabel('xsi body')
    hold off
    
    
    wing_l
    
    %trajectory_plus_strokeplane_plot(settings,pathDB,seq_nr,save_on_off,fig_nr)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

