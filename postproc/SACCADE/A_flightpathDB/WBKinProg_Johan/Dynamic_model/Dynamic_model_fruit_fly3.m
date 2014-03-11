function Dynamic_model_fruit_fly3(settings, pathDB, seq_nr)

    % Dynamical model of a fruit fly:
    
    strokeplane_WBkin = settings.strokeplane_WBkin;
    beta = (strokeplane_WBkin/180)*pi;

    dt = pathDB.t(2)-pathDB.t(1);
    
    t = pathDB.t;
    
    
    % Fruit fly mass based on wing length:
    
    wing_l = pathDB.wing_l(seq_nr);
    
    m_wl = 1e-6*0.1078*wing_l^3.008;
    
    
    % Calculate body density (3D):
    
    [cg_body, cg_L, cg_R, I_body, I_wing, I_v_wing, m_b, m_w, m_v_w ] = cg_plus_Inertia(settings, pathDB, m_wl, wing_l ,seq_nr );

    
    start = find(isnan(pathDB.x(:,seq_nr))==0, 1 );
    stop = find(isnan(pathDB.x(:,seq_nr))==0, 1, 'last' );
    
    
    % Load joint locations:
    
    J_L = pathDB.joint_pos_L(:,seq_nr) - cg_body;

    J_R = pathDB.joint_pos_R(:,seq_nr) - cg_body;
    
    
    
    % Wrenches:
    
    W_aero_L = zeros(6,stop-start+1);
    
    W_aero_R = zeros(6,stop-start+1);
    
    W_Inertia_L = zeros(6,stop-start+1);
    
    W_Inertia_R = zeros(6,stop-start+1);
    
    W_aero = zeros(6,stop-start+1);
    
    W_Inertia_w = zeros(6,stop-start+1);
    
    W_Inertia_b = zeros(6,stop-start+1);
    
    W_gravity_b = zeros(6,stop-start+1);
    
    W_Stokes_b = zeros(6,stop-start+1);
    
    
    W_Inertia_L_1 = zeros(6,stop-start+1);
    
    W_Inertia_L_2 = zeros(6,stop-start+1);
    
    Tau_L_b = zeros(3,stop-start+1);
    
    Tau_R_b = zeros(3,stop-start+1);
    
    Tau_w_b = zeros(3,stop-start+1);
    
    Tau_L_omega_dot = zeros(6,stop-start+1);
    
    Tau_L_omega = zeros(6,stop-start+1);
    
    % Inertia:
    
    T_fruitfly = zeros(3,stop-start+1);
    
    T_L = zeros(3,3,stop-start+1);
    
    T_R = zeros(3,stop-start+1);
    
    T_w1 = zeros(3,stop-start+1);
    
    T_w2 = zeros(3,stop-start+1);
    
    
    % Calculate
    
    for i = (start+1):stop
        
        
        j = i-start+1;
        
        
        % Load rotation matrices body and wings:
        
        qB = [pathDB.qb1_filt(i,seq_nr); pathDB.qb2_filt(i,seq_nr); pathDB.qb3_filt(i,seq_nr); pathDB.qb4_filt(i,seq_nr)];

        qL = [pathDB.qL1_filt2(i,seq_nr); pathDB.qL2_filt2(i,seq_nr); pathDB.qL3_filt2(i,seq_nr); pathDB.qL4_filt2(i,seq_nr)];

        qR = [pathDB.qR1_filt2(i,seq_nr); pathDB.qR2_filt2(i,seq_nr); pathDB.qR3_filt2(i,seq_nr); pathDB.qR4_filt2(i,seq_nr)];

        R_b = quat2matNEW(qB);

        R_L = quat2matNEW(qL);

        R_R = quat2matNEW(qR);
        
        
        % Load body velocity, acceleration and angular velocity:
    
        u_b = [pathDB.u_filt(i,seq_nr); pathDB.v_filt(i,seq_nr); pathDB.w_filt(i,seq_nr)];

        a_b = [pathDB.ax_filt(i,seq_nr); pathDB.ay_filt(i,seq_nr); pathDB.az_filt(i,seq_nr)];

        w_b = [pathDB.b_omega1(i,seq_nr); pathDB.b_omega2(i,seq_nr); pathDB.b_omega3(i,seq_nr)];
        
        
        % Load wing (left and right) angular velocities:
    
        w_L = [pathDB.omega1_L(i,seq_nr); pathDB.omega2_L(i,seq_nr); pathDB.omega3_L(i,seq_nr)];

        w_R = [pathDB.omega1_R(i,seq_nr); pathDB.omega2_R(i,seq_nr); pathDB.omega3_R(i,seq_nr)];
        
        
        % Compute the angular accelerations of the wings and the body:
        
        w_dot_b = [(w_b-[pathDB.b_omega1(i-1,seq_nr); pathDB.b_omega2(i-1,seq_nr); pathDB.b_omega3(i-1,seq_nr)])./dt];
    
        w_dot_L = [(w_L-[pathDB.omega1_L(i-1,seq_nr); pathDB.omega2_L(i-1,seq_nr); pathDB.omega3_L(i-1,seq_nr)])./dt];
    
        w_dot_R = [(w_R-[pathDB.omega1_R(i-1,seq_nr); pathDB.omega2_R(i-1,seq_nr); pathDB.omega3_R(i-1,seq_nr)])./dt];   
      
       
        % Calculate left wing aerodynamic wrench:

        F_aero_L = pathDB.F_L(i,:,seq_nr)';
        
        F_aero_L_arm = [0; pathDB.F_L_arm(i,seq_nr); 0];      
        
        W_aero_L(:,j) = [F_aero_L; cross(F_aero_L_arm,F_aero_L)]; % [ N ; N*mm ]
        
        
        % Calculate left wing inertial wrench:
        
            
        a_cg_L = R_L'*(1e-3*a_b+cross(w_dot_b,1e-3*J_L))+cross(w_dot_L,1e-3*cg_L);
        
        w_cg_L = R_L'*w_b + w_L;
        
        w_dot_cg_L = R_L'*w_dot_b + w_dot_L;
        
        C_w_L = 1e-3*[ 0 -cg_L(3) cg_L(2); ...
                       cg_L(3) 0 -cg_L(1); ...
                       -cg_L(2) cg_L(1) 0];
        
        I_w_L = [m_v_w*eye(3) -m_v_w*C_w_L; ...
                 m_v_w*C_w_L 1e-3*I_v_wing-m_v_w*C_w_L*C_w_L];
             
        W_Inertia_L(:,j) = I_w_L*[a_cg_L; w_dot_cg_L] + [ m_v_w*cross(w_cg_L,cross(w_cg_L,1e-3*cg_L)); cross(w_cg_L,(1e-3*I_v_wing-m_v_w*C_w_L*C_w_L)*w_cg_L) ];
        
%         beta = (55/180)*pi;
%         
%         R_stroke = [cos(beta) 0 -sin(beta); ...
%                     0 1 0; ...
%                     sin(beta) 0 cos(beta)];
        
        R_stroke = [cos(beta) 0 sin(beta); ...
                    0 1 0; ...
                    -sin(beta) 0 cos(beta)];
        
        Tau_L_b(:,j) = R_stroke*R_L*W_Inertia_L(4:6,j);
        
        Tau_L_omega_dot(:,j) = I_w_L*[a_cg_L; w_dot_cg_L];
        
        Tau_L_omega(:,j) = [ m_v_w*cross(w_cg_L,cross(w_cg_L,1e-3*cg_L)); cross(w_cg_L,(1e-3*I_v_wing-m_v_w*C_w_L*C_w_L)*w_cg_L) ];
        
        %W_Inertia_L(:,j) = I_w_L*[a_cg_L; w_dot_cg_L];
        
        
        
        
        % Calculate right wing aerodynamic wrench:
        
        F_aero_R = pathDB.F_R(i,:,seq_nr)';
        
        F_aero_R_arm = [0; pathDB.F_R_arm(i,seq_nr); 0];  
        
        W_aero_R(:,j) = [F_aero_R; cross(F_aero_R_arm,F_aero_R)];  % [ N ; N*mm ]
        

        
        
        % Calculate right wing inertial wrench:
        
       
        a_cg_R = R_R'*(1e-3*a_b+cross(w_dot_b,1e-3*J_R))+cross(w_dot_R,1e-3*cg_R);
        
        w_cg_R = R_R'*w_b + w_R;
        
        w_dot_cg_R = R_R'*w_dot_b + w_dot_R;
        
        C_w_R = 1e-3*[ 0 -cg_R(3) cg_R(2); ...
                       cg_R(3) 0 -cg_R(1); ...
                       -cg_R(2) cg_R(1) 0];
        
        I_w_R = [m_v_w*eye(3) -m_v_w*C_w_R; ...
                 m_v_w*C_w_R 1e-3*I_v_wing-m_v_w*C_w_R*C_w_R];
             
        W_Inertia_R(:,j) = I_w_R*[a_cg_R; w_dot_cg_R] + [ m_v_w*cross(w_cg_R,cross(w_cg_R,1e-3*cg_R)); cross(w_cg_R,(1e-3*I_v_wing-m_v_w*C_w_R*C_w_R)*w_cg_R) ];  % [ N ; N*mm ]
        
        Tau_R_b(:,j) = R_stroke*R_R*W_Inertia_R(4:6,j);

        %W_Inertia_R(:,j) = I_w_R*[a_cg_R; w_dot_cg_R];

        % Calculate body_wing aerodynamic wrench
        
        W_aero(:,j) = [ R_L*W_aero_L(1:3,j) + cross(R_L*W_aero_L(4:6,j),J_L); R_L*W_aero_L(4:6,j) + cross(J_L,R_L*W_aero_L(1:3,j))] + ...
                      [ R_R*W_aero_R(1:3,j) + cross(R_R*W_aero_R(4:6,j),J_R); R_R*W_aero_R(4:6,j) + cross(J_R,R_R*W_aero_R(1:3,j))]; % [ N ; N*mm ]
        
        
        % Calculate body_wing inertial wrench
        
        W_Inertia_w(:,j) = [ R_L*W_Inertia_L(1:3,j) + cross(R_L*W_Inertia_L(4:6,j),J_L); R_L*W_Inertia_L(4:6,j) + cross(J_L,R_L*W_Inertia_L(1:3,j))] + ...
                           [ R_R*W_Inertia_R(1:3,j) + cross(R_R*W_Inertia_R(4:6,j),J_R); R_R*W_Inertia_R(4:6,j) + cross(J_R,R_R*W_Inertia_R(1:3,j))]; % [ N ; N*mm ]

        
        %Tau_w_b(:,j) = Tau_L_b(:,j) + Tau_R_b(:,j) + cross(R_stroke*R_L*W_Inertia_L(1:3,j),[0; J_L(2); 0])    + cross(R_stroke*R_R*W_Inertia_R(1:3,j),[0; J_R(2); 0]);
        
        Tau_w_b(:,j) = Tau_L_b(:,j) - Tau_R_b(:,j);
        
        % Calculate virtual strokeplane inertia tensor:
        
%         beta = (55/180)*pi;
%         
%         R_stroke = [cos(beta) 0 -sin(beta); ...
%                     0 1 0; ...
%                     sin(beta) 0 cos(beta)];
                
        tau_w_L = [ R_stroke*R_L*W_Inertia_L(4:6,j) + cross([0; J_L(2); 0],R_stroke*R_L*W_Inertia_L(1:3,j))];
            
        tau_w_R = [ R_stroke*R_R*W_Inertia_R(4:6,j) + cross([0; J_R(2); 0],R_stroke*R_R*W_Inertia_R(1:3,j))];
        
        w_b_stroke = R_stroke*w_b;
        
        w_dot_b_stroke = R_stroke*w_dot_b;
        
        T_temp1 = [ 1 -w_b_stroke(2)*w_b_stroke(3)/w_dot_b_stroke(1) w_b_stroke(2)*w_b_stroke(3)/w_dot_b_stroke(1); ...
                    w_b_stroke(1)*w_b_stroke(3)/w_dot_b_stroke(2) 1 -w_b_stroke(1)*w_b_stroke(3)/w_dot_b_stroke(2); ...
                    -w_b_stroke(1)*w_b_stroke(2)/w_dot_b_stroke(3) w_b_stroke(1)*w_b_stroke(2)/w_dot_b_stroke(3) 1];
                
        T_temp2 = [tau_w_L(1)/w_dot_b_stroke(1); tau_w_L(2)/w_dot_b_stroke(2); tau_w_L(3)/w_dot_b_stroke(3)];
        
        T_temp3 = [tau_w_R(1)/w_dot_b_stroke(1); tau_w_R(2)/w_dot_b_stroke(2); tau_w_R(3)/w_dot_b_stroke(3)];


        T_temp4 = inv(T_temp1)*T_temp2;
        
        T_temp5 = inv(T_temp1)*T_temp3;
        
                             
        Tw1(:,j) = T_temp4;       
        
        Tw2(:,j) = T_temp5;      
                       
                       
        % Calculate body inertial wrench
        
        W_Inertia_b(:,j) = [ m_wl*eye(3) zeros(3); zeros(3) 1e-3*I_body ] * [1e-3*a_b; w_dot_b] + [ zeros(3,1) ; cross(w_b,1e-3*I_body*w_b) ]; % [ N ; N*mm ]

        
        % Calculate body gravitational wrench
        
        W_gravity_b(:,j) = [ R_b'*[0; 0; -m_wl*9.81]; zeros(3,1) ]; % [ N ; N*mm ]
        
        
%         % Calculate fruit fly inertia tensor
%         
%         T_temp1 = [ 1 -w_b(2)*w_b(3)/w_dot_b(1) w_b(2)*w_b(3)/w_dot_b(1); ...
%                     w_b(1)*w_b(3)/w_dot_b(2) 1 -w_b(1)*w_b(3)/w_dot_b(2); ...
%                     -w_b(1)*w_b(2)/w_dot_b(3) w_b(1)*w_b(2)/w_dot_b(3) 1];
%         
%         T_temp2 = [ (W_Inertia_b(4,j)+W_Inertia_w(4,j))/w_dot_b(1); (W_Inertia_b(5,j)+W_Inertia_w(5,j))/w_dot_b(2); (W_Inertia_b(6,j)+W_Inertia_w(6,j))/w_dot_b(3)];
%      
%         %T_temp2 = [ (W_aero(4,j)+W_gravity_b(4,j))/w_dot_b(1); (W_aero(5,j)+W_gravity_b(5,j))/w_dot_b(2); (W_aero(6,j)+W_gravity_b(6,j))/w_dot_b(3)];
% 
%         T_fruitfly(:,j) = inv(T_temp1)*T_temp2;
        
        
    end
    
    % Set first frame to the second frame:
    
    W_aero_L(:,1) = W_aero_L(:,2);
    
    W_aero_R(:,1) = W_aero_R(:,2);
    
    W_Inertia_L(:,1) = W_Inertia_L(:,2);
    
    W_Inertia_R(:,1) = W_Inertia_R(:,2);
    
    W_aero(:,1) = W_aero(:,2);
    
    W_Inertia_w(:,1) = W_Inertia_w(:,2);
    
    W_Inertia_b(:,1) = W_Inertia_b(:,2);
    
    W_gravity_b(:,1) = W_gravity_b(:,2);
    
    
    
%     figure()
%     plot(W_Inertia_1(1,:),'r')
%     
%     figure()
%     plot(W_Inertia_2(1,:),'b')
%     
%     figure()
%     plot(W_Inertia_1(2,:),'r')
%     
%     figure()
%     plot(W_Inertia_2(2,:),'b')
%     
%     figure()
%     plot(W_Inertia_1(3,:),'r')
%     
%     figure()
%     plot(W_Inertia_2(3,:),'b')
%     
%     figure()
%     plot(W_Inertia_1(4,:),'r')
%     
%     figure()
%     plot(W_Inertia_2(4,:),'b')
%     
%     figure()
%     plot(W_Inertia_1(5,:),'r')
%     
%     figure()
%     plot(W_Inertia_2(5,:),'b')
%     
%     figure()
%     plot(W_Inertia_1(6,:),'r')
%     
%     figure()
%     plot(W_Inertia_2(6,:),'b')
    
    
    
    
    % Plot the body inertial wrench:
    
%     figure()
%     subplot(3,1,1); plot(W_Inertia_b(1,:))
%     subplot(3,1,2); plot(W_Inertia_b(2,:))
%     subplot(3,1,3); plot(W_Inertia_b(3,:))
%     
%     figure()
%     subplot(3,1,1); plot(W_Inertia_b(4,:))
%     subplot(3,1,2); plot(W_Inertia_b(5,:))
%     subplot(3,1,3); plot(W_Inertia_b(6,:))
%     
%     
    figure()
    subplot(3,1,1); plot(W_Inertia_w(1,:))
    subplot(3,1,2); plot(W_Inertia_w(2,:))
    subplot(3,1,3); plot(W_Inertia_w(3,:))
    
    figure()
    subplot(3,1,1); plot(W_Inertia_w(4,:))
    subplot(3,1,2); plot(W_Inertia_w(5,:))
    subplot(3,1,3); plot(W_Inertia_w(6,:))
% 
% 
    figure()
    plot(W_Inertia_b(1,:)-W_Inertia_w(1,:),'r')
    hold on
    plot(W_Inertia_b(2,:)-W_Inertia_w(2,:),'b')
    plot(W_Inertia_b(3,:)-W_Inertia_w(3,:),'g')
    hold off
    
    figure()
    plot(W_Inertia_b(4,:)-W_Inertia_w(4,:),'r')
    hold on
    plot(W_Inertia_b(5,:)-W_Inertia_w(5,:),'b')
    plot(W_Inertia_b(6,:)-W_Inertia_w(6,:),'g')
    hold off
    
    
    nr_of_wb = find(isnan(pathDB.wingbeat_time(:,1,seq_nr))==0, 1, 'last' );
    
    x_inertia_b_wb = zeros(nr_of_wb,100);
    
    y_inertia_b_wb = zeros(nr_of_wb,100);
    
    z_inertia_b_wb = zeros(nr_of_wb,100);
    
    for k = 1:nr_of_wb
       
        wb_end = find(isnan(pathDB.wingbeat_time(k,:,seq_nr))==0, 1, 'last' );
        
        sect_wb = pathDB.wingbeat_time(k,1,seq_nr):((pathDB.wingbeat_time(k,wb_end,seq_nr)-pathDB.wingbeat_time(k,1,seq_nr))/99):pathDB.wingbeat_time(k,wb_end,seq_nr);
        
        a = pathDB.wingbeat_time(k,1,seq_nr):pathDB.wingbeat_time(k,wb_end,seq_nr);
        
        x_b_inertia_wb(k,:) = interp1(a,W_Inertia_b(1,a)-W_Inertia_w(1,a),sect_wb);
        
        y_b_inertia_wb(k,:) = interp1(a,W_Inertia_b(2,a)-W_Inertia_w(2,a),sect_wb);
        
        z_b_inertia_wb(k,:) = interp1(a,W_Inertia_b(3,a)-W_Inertia_w(3,a),sect_wb);
    end
    
    x_inertia_b_wb_mean = zeros(1,100);
    
    y_inertia_b_wb_mean = zeros(1,100);
    
    z_inertia_b_wb_mean = zeros(1,100);
    
    for m = 1:100
        
        x_inertia_b_wb_mean(m) = mean(x_b_inertia_wb(:,m));
        
        y_inertia_b_wb_mean(m) = mean(y_b_inertia_wb(:,m));
        
        z_inertia_b_wb_mean(m) = mean(z_b_inertia_wb(:,m));
        
    end
  
    figure()
    hold on
    for k = 1:nr_of_wb
        
        plot(x_b_inertia_wb(k,:))
        
    end
    plot(zeros(1,100),'r')
    plot(x_inertia_b_wb_mean,'g')
    hold off
    
    figure()
    hold on
    for k = 1:nr_of_wb
        
        plot(y_b_inertia_wb(k,:))
        
    end
    plot(zeros(1,100),'r')
    plot(y_inertia_b_wb_mean,'g')
    hold off
    
    figure()
    hold on
    for k = 1:nr_of_wb
        
        plot(z_b_inertia_wb(k,:))
        
    end
    plot(zeros(1,100),'r')
    plot(z_inertia_b_wb_mean,'g')
    hold off
    
    
    
%     
%     figure()
%     plot(W_aero(1,:),'r')
%     hold on
%     plot(W_aero(2,:),'b')
%     plot(W_aero(3,:),'g')
%     hold off
%     
%     figure()
%     plot(W_aero(4,:),'r')
%     hold on
%     plot(W_aero(5,:),'b')
%     plot(W_aero(6,:),'g')
%     hold off
%     
%     figure()
%     plot(W_gravity_b(1,:),'r')
%     hold on
%     plot(W_gravity_b(2,:),'b')
%     plot(W_gravity_b(3,:),'g')
%     hold off
%     
%     figure()
%     plot(W_gravity_b(4,:),'r')
%     hold on
%     plot(W_gravity_b(5,:),'b')
%     plot(W_gravity_b(6,:),'g')
%     hold off
% 
% 
%     figure()
%     subplot(3,1,1); plot(t(start:stop),-W_Inertia_b(1,:)-W_Inertia_w(1,:)+W_aero(1,:)+W_gravity_b(1,:),t(start:stop),W_aero(1,:)+W_gravity_b(1,:))
%     subplot(3,1,2); plot(t(start:stop),-W_Inertia_b(2,:)-W_Inertia_w(2,:)+W_aero(2,:)+W_gravity_b(2,:),t(start:stop),W_aero(2,:)+W_gravity_b(2,:))
%     subplot(3,1,3); plot(t(start:stop),-W_Inertia_b(3,:)-W_Inertia_w(3,:)+W_aero(3,:)+W_gravity_b(3,:),t(start:stop),W_aero(3,:)+W_gravity_b(3,:))
%     legend('resultant','aero and gravity')
%     
%     figure()
%     subplot(3,1,1); plot(t(start:stop),-W_Inertia_b(4,:)-W_Inertia_w(4,:)+W_aero(4,:)+W_gravity_b(4,:),t(start:stop),W_aero(4,:)+W_gravity_b(4,:))
%     subplot(3,1,2); plot(t(start:stop),-W_Inertia_b(5,:)-W_Inertia_w(5,:)+W_aero(5,:)+W_gravity_b(5,:),t(start:stop),W_aero(5,:)+W_gravity_b(5,:))
%     subplot(3,1,3); plot(t(start:stop),-W_Inertia_b(6,:)-W_Inertia_w(6,:)+W_aero(6,:)+W_gravity_b(6,:),t(start:stop),W_aero(6,:)+W_gravity_b(6,:))
%     legend('resultant','aero and gravity')
%     
%     figure()
%     subplot(3,1,1); plot(t(start:stop),W_Inertia_b(1,:)+W_Inertia_w(1,:),t(start:stop),W_aero(1,:)+W_gravity_b(1,:),t(start:stop),W_Inertia_b(1,:)+W_Inertia_w(1,:)+W_aero(1,:)+W_gravity_b(1,:))
%     subplot(3,1,2); plot(t(start:stop),W_Inertia_b(2,:)+W_Inertia_w(2,:),t(start:stop),W_aero(2,:)+W_gravity_b(2,:),t(start:stop),W_Inertia_b(2,:)+W_Inertia_w(2,:)+W_aero(2,:)+W_gravity_b(2,:))
%     subplot(3,1,3); plot(t(start:stop),W_Inertia_b(3,:)+W_Inertia_w(3,:),t(start:stop),W_aero(3,:)+W_gravity_b(3,:),t(start:stop),W_Inertia_b(3,:)+W_Inertia_w(3,:)+W_aero(3,:)+W_gravity_b(3,:))
%     legend('inertia','aero and gravity')
%     
%     figure()
%     subplot(3,1,1); plot(t(start:stop),W_Inertia_b(4,:)+W_Inertia_w(4,:),t(start:stop),W_aero(4,:)+W_gravity_b(4,:),t(start:stop),W_Inertia_b(4,:)+W_Inertia_w(4,:)+W_aero(4,:)+W_gravity_b(4,:))
%     subplot(3,1,2); plot(t(start:stop),W_Inertia_b(5,:)+W_Inertia_w(5,:),t(start:stop),W_aero(5,:)+W_gravity_b(5,:),t(start:stop),W_Inertia_b(5,:)+W_Inertia_w(5,:)+W_aero(5,:)+W_gravity_b(5,:))
%     subplot(3,1,3); plot(t(start:stop),W_Inertia_b(6,:)+W_Inertia_w(6,:),t(start:stop),W_aero(6,:)+W_gravity_b(6,:),t(start:stop),W_Inertia_b(6,:)+W_Inertia_w(6,:)+W_aero(6,:)+W_gravity_b(6,:))
%     legend('inertia','aero and gravity')
%     
%     figure()
%     subplot(3,1,1); plot(t(start:stop),W_Inertia_L(1,:),'r')
%     subplot(3,1,2); plot(t(start:stop),W_Inertia_L(2,:),'b')
%     subplot(3,1,3); plot(t(start:stop),W_Inertia_L(3,:),'g')
%     
%     figure()
%     subplot(3,1,1); plot(t(start:stop),W_Inertia_L(4,:),'r')
%     subplot(3,1,2); plot(t(start:stop),W_Inertia_L(5,:),'b')
%     subplot(3,1,3); plot(t(start:stop),W_Inertia_L(6,:),'g')
    
%     T_m1 = nanmean(T_L(1,1,:));
%     
%     T_m2 = nanmean(T_L(2,2,:));
%     
%     T_m3 = nanmean(T_L(3,3,:));
%     
%     T_m1*1e3
%    
%     T_m2*1e3
%     
%     T_m3*1e3
%     
%     I_v_wing
%     
%     
%     figure()
%     subplot(3,1,1); plot(t(start:stop),T_fruitfly(1,:))
%     subplot(3,1,2); plot(t(start:stop),T_fruitfly(2,:))
%     subplot(3,1,3); plot(t(start:stop),T_fruitfly(3,:))

    figure()
    subplot(3,1,1); plot(t(start:stop),Tau_L_b(1,:),t(start:stop),Tau_R_b(1,:),t(start:stop),Tau_w_b(1,:))
    subplot(3,1,2); plot(t(start:stop),Tau_L_b(2,:),t(start:stop),Tau_R_b(2,:),t(start:stop),Tau_w_b(2,:))
    subplot(3,1,3); plot(t(start:stop),Tau_L_b(3,:),t(start:stop),Tau_R_b(3,:),t(start:stop),Tau_w_b(3,:))      
    
    figure()
    subplot(3,1,1); plot(t(start:stop),Tau_L_omega_dot(4,:),t(start:stop),Tau_L_omega(4,:))
    subplot(3,1,2); plot(t(start:stop),Tau_L_omega_dot(5,:),t(start:stop),Tau_L_omega(5,:))
    subplot(3,1,3); plot(t(start:stop),Tau_L_omega_dot(6,:),t(start:stop),Tau_L_omega(6,:)) 
    
    
%     nr_of_wb = find(isnan(pathDB.wingbeat_time(:,1,seq_nr))==0, 1, 'last' );
    
    x_inertia_wb = zeros(nr_of_wb,100);

    y_inertia_wb = zeros(nr_of_wb,100);
    
    z_inertia_wb = zeros(nr_of_wb,100);
    
    for k = 1:nr_of_wb
       
        wb_end = find(isnan(pathDB.wingbeat_time(k,:,seq_nr))==0, 1, 'last' );
        
        sect_wb = pathDB.wingbeat_time(k,1,seq_nr):((pathDB.wingbeat_time(k,wb_end,seq_nr)-pathDB.wingbeat_time(k,1,seq_nr))/99):pathDB.wingbeat_time(k,wb_end,seq_nr);
        
        a = pathDB.wingbeat_time(k,1,seq_nr):pathDB.wingbeat_time(k,wb_end,seq_nr);
        
        x_inertia_wb(k,:) = interp1(a,Tau_w_b(1,a),sect_wb);
        
        y_inertia_wb(k,:) = interp1(a,Tau_w_b(2,a),sect_wb);
        
        z_inertia_wb(k,:) = interp1(a,Tau_w_b(3,a),sect_wb);
    end
  
    figure()
    hold on
    for k = 1:nr_of_wb
        
        plot(x_inertia_wb(k,:))
        
    end
    plot(zeros(1,100),'r')
    hold off
    
    figure()
    hold on
    for k = 1:nr_of_wb
        
        plot(y_inertia_wb(k,:))
        
    end
    plot(zeros(1,100),'r')
    hold off
    
    figure()
    hold on
    for k = 1:nr_of_wb
        
        plot(z_inertia_wb(k,:))
        
    end
    plot(zeros(1,100),'r')
    hold off


%     Tw1
%     
%     Tw2
% 
%     figure()
%     subplot(3,1,1); plot(t(start:stop),Tw1(1,:))
%     subplot(3,1,2); plot(t(start:stop),Tw1(2,:))
%     subplot(3,1,3); plot(t(start:stop),Tw1(3,:))    
%     
%     figure()
%     subplot(3,1,1); plot(t(start:stop),Tw2(1,:))
%     subplot(3,1,2); plot(t(start:stop),Tw2(2,:))
%     subplot(3,1,3); plot(t(start:stop),Tw2(3,:))    
    
    
%     nanmean(T_fruitfly(1,:))*1e6
%     
%     nanmean(T_fruitfly(2,:))*1e6
%     
%     nanmean(T_fruitfly(3,:))*1e6
%     
%     I_body

end
