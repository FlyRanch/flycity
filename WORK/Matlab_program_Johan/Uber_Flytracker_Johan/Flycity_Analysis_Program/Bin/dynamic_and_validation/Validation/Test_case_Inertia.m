function Test_case_Inertia( settings, pathDB )


    % Test the Inertia model for conservation of linear and angular
    % momentum:
    
    seq_nr      = 1;
    nr_points   = 1000;
    nr_wb       = 5;
    
    % Setup the body_model and wing_model:
    
    body_model.mass_fly         = pathDB.body_model.mass_fly(seq_nr);
    body_model.mass_body        = pathDB.body_model.mass_body(seq_nr);
    body_model.Joint_left       = pathDB.body_model.Joint_left(seq_nr,:)';
    body_model.Joint_right      = pathDB.body_model.Joint_right(seq_nr,:)';
    body_model.cg_b             = pathDB.body_model.cg(seq_nr,:)';
    body_model.Inertia          = pathDB.body_model.Inertia(:,:,seq_nr);
    body_model.g                = 1e-3*settings.g;
    body_model.x_mod            = pathDB.body_model.x_mod(:,:,seq_nr);
    body_model.y_mod            = pathDB.body_model.y_mod(:,:,seq_nr);
    body_model.z_mod            = pathDB.body_model.z_mod(:,:,seq_nr);
    
    wing_model.virtual_mass     = pathDB.wing_model.virtual_mass(seq_nr);
    wing_model.wing_cg_L        = pathDB.wing_model.wing_cg_L(seq_nr,:)';
    wing_model.wing_cg_R        = pathDB.wing_model.wing_cg_R(seq_nr,:)';
    wing_model.virtual_Inertia  = pathDB.wing_model.virtual_Inertia(:,:,seq_nr);
    wing_model.y_sect_L         = pathDB.wing_model.y_sect_L(:,:,seq_nr)';
    wing_model.chords_L         = pathDB.wing_model.chords_L(seq_nr,:)';
    wing_model.y_sect_R         = pathDB.wing_model.y_sect_R(:,:,seq_nr)';
    wing_model.chords_R         = pathDB.wing_model.chords_R(seq_nr,:)';
    wing_model.rho              = settings.rho_air;
    wing_model.length           = pathDB.wing_model.length(seq_nr);
    
    f                           = pathDB.poly_fit.a_avg.f(seq_nr);
    dt                          = (1/f)/(nr_points-1);
    R_strk                      = pathDB.rot_mat.Rstr;
    
    % Create the wing kinematics:
    
    t_func          = nan(1,nr_points*nr_wb+1);
    theta_L         = nan(1,nr_points*nr_wb+1);
    eta_L           = nan(1,nr_points*nr_wb+1);
    phi_L           = nan(1,nr_points*nr_wb+1);
    theta_dot_L     = nan(1,nr_points*nr_wb+1);
    eta_dot_L       = nan(1,nr_points*nr_wb+1);
    phi_dot_L       = nan(1,nr_points*nr_wb+1);
    theta_ddot_L    = nan(1,nr_points*nr_wb+1);
    eta_ddot_L      = nan(1,nr_points*nr_wb+1);
    phi_ddot_L      = nan(1,nr_points*nr_wb+1);
    theta_R         = nan(1,nr_points*nr_wb+1);
    eta_R           = nan(1,nr_points*nr_wb+1);
    phi_R           = nan(1,nr_points*nr_wb+1);
    theta_dot_R     = nan(1,nr_points*nr_wb+1);
    eta_dot_R       = nan(1,nr_points*nr_wb+1);
    phi_dot_R       = nan(1,nr_points*nr_wb+1);
    theta_ddot_R    = nan(1,nr_points*nr_wb+1);
    eta_ddot_R      = nan(1,nr_points*nr_wb+1);
    phi_ddot_R      = nan(1,nr_points*nr_wb+1);
    qL              = nan(4,nr_points*nr_wb+1);
    qR              = nan(4,nr_points*nr_wb+1);
    RL              = nan(3,3,nr_points*nr_wb+1);
    RR              = nan(3,3,nr_points*nr_wb+1);
    wL              = nan(3,nr_points*nr_wb+1);
    wR              = nan(3,nr_points*nr_wb+1);
    w_dot_L         = nan(3,nr_points*nr_wb+1);
    w_dot_R         = nan(3,nr_points*nr_wb+1);
    wL_b            = nan(3,nr_points*nr_wb+1);
    wR_b            = nan(3,nr_points*nr_wb+1);
    w_dot_L_b       = nan(3,nr_points*nr_wb+1);
    w_dot_R_b       = nan(3,nr_points*nr_wb+1);
    wL_strk         = nan(3,nr_points*nr_wb+1);
    wR_strk         = nan(3,nr_points*nr_wb+1);
    w_dot_L_strk    = nan(3,nr_points*nr_wb+1);
    w_dot_R_strk    = nan(3,nr_points*nr_wb+1);
    
    theta_L0      = (5*pi)/180;
    eta_L0        = pi/2;
    phi_L0        = (10*pi)/180;
    theta_R0      = (5*pi)/180;
    eta_R0        = pi/2;
    phi_R0        = (10*pi)/180;
    A_theta_L     = (10*pi)/180;
    A_eta_L       = -(45*pi)/180;
    A_phi_L       = (60*pi)/180;
    A_theta_R     = (10*pi)/180;
    A_eta_R       = -(45*pi)/180;
    A_phi_R       = (60*pi)/180;
     
    for i = 1:nr_wb
        
        a_fit.theta_L0      = theta_L0;
        a_fit.eta_L0        = eta_L0;
        a_fit.phi_L0        = phi_L0;
        a_fit.theta_R0      = theta_R0;
        a_fit.eta_R0        = eta_R0;
        a_fit.phi_R0        = phi_R0;
        a_fit.A_theta_L     = A_theta_L;
        a_fit.A_eta_L       = A_eta_L;
        a_fit.A_phi_L       = A_phi_L;
        a_fit.A_theta_R     = A_theta_R;
        a_fit.A_eta_R       = A_eta_R;
        a_fit.A_phi_R       = A_phi_R;
        a_fit.f             = f;
        a_fit.nr_points     = nr_points+1;
        a_fit.R_strk        = R_strk;
        
        [ kine ] = angular_velocities_sine( a_fit );
        
        if i == 1
            
            range_i         = ((i-1)*nr_points+1):(i*nr_points+1);
            range_temp      = 1:(nr_points+1);
            
            t_func(range_i)         = kine.t(range_temp);
            theta_L(range_i)        = kine.theta_L(range_temp);
            eta_L(range_i)          = kine.eta_L(range_temp);
            phi_L(range_i)          = kine.phi_L(range_temp);
            theta_dot_L(range_i)    = kine.theta_dot_L(range_temp);
            eta_dot_L(range_i)      = kine.eta_dot_L(range_temp);
            phi_dot_L(range_i)      = kine.phi_dot_L(range_temp);
            theta_ddot_L(range_i)   = kine.theta_ddot_L(range_temp);
            eta_ddot_L(range_i)     = kine.eta_ddot_L(range_temp);
            phi_ddot_L(range_i)     = kine.phi_ddot_L(range_temp);
            theta_R(range_i)        = kine.theta_R(range_temp);
            eta_R(range_i)          = kine.eta_R(range_temp);
            phi_R(range_i)          = kine.phi_R(range_temp);
            theta_dot_R(range_i)    = kine.theta_dot_R(range_temp);
            eta_dot_R(range_i)      = kine.eta_dot_R(range_temp);
            phi_dot_R(range_i)      = kine.phi_dot_R(range_temp);
            theta_ddot_R(range_i)   = kine.theta_ddot_R(range_temp);
            eta_ddot_R(range_i)     = kine.eta_ddot_R(range_temp);
            phi_ddot_R(range_i)     = kine.phi_ddot_R(range_temp);
            RL(:,:,range_i)         = kine.RL(:,:,range_temp);
            RR(:,:,range_i)         = kine.RR(:,:,range_temp);
            qL(:,range_i)           = kine.qL(:,range_temp);
            qR(:,range_i)           = kine.qR(:,range_temp);
            wL(:,range_i)           = kine.wL(:,range_temp);
            wR(:,range_i)           = kine.wR(:,range_temp);
            w_dot_L(:,range_i)      = kine.w_dot_L(:,range_temp);
            w_dot_R(:,range_i)      = kine.w_dot_R(:,range_temp);
            wL_b(:,range_i)         = kine.wL_b(:,range_temp);
            wR_b(:,range_i)         = kine.wR_b(:,range_temp);
            w_dot_L_b(:,range_i)    = kine.w_dot_L_b(:,range_temp);
            w_dot_R_b(:,range_i)    = kine.w_dot_R_b(:,range_temp);
            wL_strk(:,range_i)      = kine.wL_strk(:,range_temp);
            wR_strk(:,range_i)      = kine.wR_strk(:,range_temp);
            w_dot_L_strk(:,range_i) = kine.w_dot_L_strk(:,range_temp);
            w_dot_R_strk(:,range_i) = kine.w_dot_R_strk(:,range_temp);
            
        elseif i == nr_wb && i > 1
            
            range_i         = ((i-1)*nr_points+1):(i*nr_points+1);
            range_temp      = 1:(nr_points+1);
            
            t_func(range_i)         = t_func((i-1)*nr_points)+kine.dt+kine.t(range_temp);
            theta_L(range_i)        = kine.theta_L(range_temp);
            eta_L(range_i)          = kine.eta_L(range_temp);
            phi_L(range_i)          = kine.phi_L(range_temp);
            theta_dot_L(range_i)    = kine.theta_dot_L(range_temp);
            eta_dot_L(range_i)      = kine.eta_dot_L(range_temp);
            phi_dot_L(range_i)      = kine.phi_dot_L(range_temp);
            theta_ddot_L(range_i)   = kine.theta_ddot_L(range_temp);
            eta_ddot_L(range_i)     = kine.eta_ddot_L(range_temp);
            phi_ddot_L(range_i)     = kine.phi_ddot_L(range_temp);
            theta_R(range_i)        = kine.theta_R(range_temp);
            eta_R(range_i)          = kine.eta_R(range_temp);
            phi_R(range_i)          = kine.phi_R(range_temp);
            theta_dot_R(range_i)    = kine.theta_dot_R(range_temp);
            eta_dot_R(range_i)      = kine.eta_dot_R(range_temp);
            phi_dot_R(range_i)      = kine.phi_dot_R(range_temp);
            theta_ddot_R(range_i)   = kine.theta_ddot_R(range_temp);
            eta_ddot_R(range_i)     = kine.eta_ddot_R(range_temp);
            phi_ddot_R(range_i)     = kine.phi_ddot_R(range_temp);
            RL(:,:,range_i)         = kine.RL(:,:,range_temp);
            RR(:,:,range_i)         = kine.RR(:,:,range_temp);
            qL(:,range_i)           = kine.qL(:,range_temp);
            qR(:,range_i)           = kine.qR(:,range_temp);
            wL(:,range_i)           = kine.wL(:,range_temp);
            wR(:,range_i)           = kine.wR(:,range_temp);
            w_dot_L(:,range_i)      = kine.w_dot_L(:,range_temp);
            w_dot_R(:,range_i)      = kine.w_dot_R(:,range_temp);
            wL_b(:,range_i)         = kine.wL_b(:,range_temp);
            wR_b(:,range_i)         = kine.wR_b(:,range_temp);
            w_dot_L_b(:,range_i)    = kine.w_dot_L_b(:,range_temp);
            w_dot_R_b(:,range_i)    = kine.w_dot_R_b(:,range_temp);
            wL_strk(:,range_i)      = kine.wL_strk(:,range_temp);
            wR_strk(:,range_i)      = kine.wR_strk(:,range_temp);
            w_dot_L_strk(:,range_i) = kine.w_dot_L_strk(:,range_temp);
            w_dot_R_strk(:,range_i) = kine.w_dot_R_strk(:,range_temp);
           
        else
            
            range_i         = ((i-1)*nr_points+1):(i*nr_points+1);
            
            t_func(range_i)         = t_func((i-1)*nr_points)+kine.dt+kine.t;
            theta_L(range_i)        = kine.theta_L;
            eta_L(range_i)          = kine.eta_L;
            phi_L(range_i)          = kine.phi_L;
            theta_dot_L(range_i)    = kine.theta_dot_L;
            eta_dot_L(range_i)      = kine.eta_dot_L;
            phi_dot_L(range_i)      = kine.phi_dot_L;
            theta_ddot_L(range_i)   = kine.theta_ddot_L;
            eta_ddot_L(range_i)     = kine.eta_ddot_L;
            phi_ddot_L(range_i)     = kine.phi_ddot_L;
            theta_R(range_i)        = kine.theta_R;
            eta_R(range_i)          = kine.eta_R;
            phi_R(range_i)          = kine.phi_R;
            theta_dot_R(range_i)    = kine.theta_dot_R;
            eta_dot_R(range_i)      = kine.eta_dot_R;
            phi_dot_R(range_i)      = kine.phi_dot_R;
            theta_ddot_R(range_i)   = kine.theta_ddot_R;
            eta_ddot_R(range_i)     = kine.eta_ddot_R;
            phi_ddot_R(range_i)     = kine.phi_ddot_R;
            RL(:,:,range_i)         = kine.RL;
            RR(:,:,range_i)         = kine.RR;
            qL(:,range_i)           = kine.qL(:,range_temp);
            qR(:,range_i)           = kine.qR(:,range_temp);
            wL(:,range_i)           = kine.wL;
            wR(:,range_i)           = kine.wR;
            w_dot_L(:,range_i)      = kine.w_dot_L;
            w_dot_R(:,range_i)      = kine.w_dot_R;
            wL_b(:,range_i)         = kine.wL_b;
            wR_b(:,range_i)         = kine.wR_b;
            w_dot_L_b(:,range_i)    = kine.w_dot_L_b;
            w_dot_R_b(:,range_i)    = kine.w_dot_R_b;
            wL_strk(:,range_i)      = kine.wL_strk;
            wR_strk(:,range_i)      = kine.wR_strk;
            w_dot_L_strk(:,range_i) = kine.w_dot_L_strk;
            w_dot_R_strk(:,range_i) = kine.w_dot_R_strk;
            
        end
        
        clear kine a_fit
        
    end
    
    
    Fly_plot_3D_wingtip( [0; 0; 0], R_strk'*[1 0 0; 0 -1 0; 0 0 -1], RL(:,:,1:20:end), RR(:,:,1:20:end), body_model, wing_model )
    
    pause
    
    % Wing kinematics:
    
    wing_kin.R_strk     = R_strk;
    wing_kin.RL         = RL;
    wing_kin.RR         = RR;
    wing_kin.qL         = qL;
    wing_kin.qR         = qR;
    wing_kin.wL         = wL;
    wing_kin.wR         = wR;
    wing_kin.wL_b       = wL_b;
    wing_kin.wR_b       = wR_b;
    wing_kin.w_dot_L    = w_dot_L;
    wing_kin.w_dot_R    = w_dot_R;
    wing_kin.w_dot_L_b  = w_dot_L_b;
    wing_kin.w_dot_R_b  = w_dot_R_b;
    
    % Solve the ODE:
    
    t_span = [ t_func(1) t_func(end) ];
    IC = zeros(6,1);
    
    options = odeset('RelTol',1e-4,'AbsTol',[1e-10 1e-10 1e-10 1e-10 1e-10 1e-10]);
    [ T X ] = ode45(@(t,x) state_space_func(t,x,t_func,wing_kin,body_model,wing_model),t_span,IC,options);
    
    X_dot = zeros(6,length(T));
    
    for j = 1:length(T)
    
        X_dot(:,j) = state_space_func(T(j),X(j,:)',t_func,wing_kin,body_model,wing_model);
        
    end
    
    
    
    X_mean_1 = mean(interp1(T,X(:,1),linspace(t_func(1),t_func(end),100)))
    X_mean_2 = mean(interp1(T,X(:,2),linspace(t_func(1),t_func(end),100)))
    X_mean_3 = mean(interp1(T,X(:,3),linspace(t_func(1),t_func(end),100)))
    X_mean_4 = mean(interp1(T,X(:,4),linspace(t_func(1),t_func(end),100)))
    X_mean_5 = mean(interp1(T,X(:,5),linspace(t_func(1),t_func(end),100)))
    X_mean_6 = mean(interp1(T,X(:,6),linspace(t_func(1),t_func(end),100)))
    
    X_dot_mean_1 = mean(interp1(T,X_dot(1,:),linspace(t_func(1),t_func(end),100)))
    X_dot_mean_2 = mean(interp1(T,X_dot(2,:),linspace(t_func(1),t_func(end),100)))
    X_dot_mean_3 = mean(interp1(T,X_dot(3,:),linspace(t_func(1),t_func(end),100)))
    X_dot_mean_4 = mean(interp1(T,X_dot(4,:),linspace(t_func(1),t_func(end),100)))
    X_dot_mean_5 = mean(interp1(T,X_dot(5,:),linspace(t_func(1),t_func(end),100)))
    X_dot_mean_6 = mean(interp1(T,X_dot(6,:),linspace(t_func(1),t_func(end),100)))
    
    body_kin.vb         = X(:,1:3)';
    body_kin.wb         = X(:,4:6)';
    body_kin.ab         = X_dot(1:3,:);
    body_kin.w_dot_b    = X_dot(4:6,:);
    
    [ FMI_b, FMI_strk ] = Inertia_FM( body_kin, wing_kin, body_model, wing_model, T, t_func );
    
    vb_0 = FMI_b.vb_0(:,1)
    wb_0 = FMI_b.wb_0(:,1)
    ab_0 = X_dot(1:3,1)
    w_dot_b_0 = X_dot(4:6,1)
    
    % Solve the ODE:
    
    t_span = [ t_func(1) t_func(end) ];
    IC = [vb_0; wb_0];
    
    options = odeset('RelTol',1e-7,'AbsTol',[1e-15 1e-15 1e-15 1e-15 1e-15 1e-15]);
    [ T X ] = ode45(@(t,x) state_space_func(t,x,t_func,wing_kin,body_model,wing_model),t_span,IC,options);
    
    X_dot = zeros(6,length(T));
    
    for j = 1:length(T)
    
        X_dot(:,j) = state_space_func(T(j),X(j,:)',t_func,wing_kin,body_model,wing_model);
        
    end
    
    
    
    X_mean_1 = mean(interp1(T,X(:,1),linspace(t_func(1),t_func(end),100)))
    X_mean_2 = mean(interp1(T,X(:,2),linspace(t_func(1),t_func(end),100)))
    X_mean_3 = mean(interp1(T,X(:,3),linspace(t_func(1),t_func(end),100)))
    X_mean_4 = mean(interp1(T,X(:,4),linspace(t_func(1),t_func(end),100)))
    X_mean_5 = mean(interp1(T,X(:,5),linspace(t_func(1),t_func(end),100)))
    X_mean_6 = mean(interp1(T,X(:,6),linspace(t_func(1),t_func(end),100)))
    
    X_dot_mean_1 = mean(interp1(T,X_dot(1,:),linspace(t_func(1),t_func(end),100)))
    X_dot_mean_2 = mean(interp1(T,X_dot(2,:),linspace(t_func(1),t_func(end),100)))
    X_dot_mean_3 = mean(interp1(T,X_dot(3,:),linspace(t_func(1),t_func(end),100)))
    X_dot_mean_4 = mean(interp1(T,X_dot(4,:),linspace(t_func(1),t_func(end),100)))
    X_dot_mean_5 = mean(interp1(T,X_dot(5,:),linspace(t_func(1),t_func(end),100)))
    X_dot_mean_6 = mean(interp1(T,X_dot(6,:),linspace(t_func(1),t_func(end),100)))
    
    body_kin.vb         = X(:,1:3)';
    body_kin.wb         = X(:,4:6)';
    body_kin.ab         = X_dot(1:3,:);
    body_kin.w_dot_b    = X_dot(4:6,:);
    
    [ FMI_b, FMI_strk ] = Inertia_FM( body_kin, wing_kin, body_model, wing_model, T, t_func );
       
    FI_mean = mean(FMI_b.FI_vel,2)
    
    MI_mean = mean(FMI_b.MI_vel,2)
    
    LM_mean = mean(FMI_b.LM,2)
    
    AM_mean = mean(FMI_b.AM,2)
    
    figure()
    hold on
    subplot(3,1,1); plot(T,X(:,1))
    subplot(3,1,2); plot(T,X(:,2))
    subplot(3,1,3); plot(T,X(:,3))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(T,X(:,4))
    subplot(3,1,2); plot(T,X(:,5))
    subplot(3,1,3); plot(T,X(:,6))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(T,X_dot(1,:))
    subplot(3,1,2); plot(T,X_dot(2,:))
    subplot(3,1,3); plot(T,X_dot(3,:))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(T,X_dot(4,:))
    subplot(3,1,2); plot(T,X_dot(5,:))
    subplot(3,1,3); plot(T,X_dot(6,:))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t_func,FMI_b.FI_vel(1,:),'b',t_func,FMI_b.FI_acc(1,:),'r')
    subplot(3,1,2); plot(t_func,FMI_b.FI_vel(2,:),'b',t_func,FMI_b.FI_acc(2,:),'r')
    subplot(3,1,3); plot(t_func,FMI_b.FI_vel(3,:),'b',t_func,FMI_b.FI_acc(3,:),'r')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t_func,FMI_b.MI_vel(1,:),'b',t_func,FMI_b.MI_acc(1,:),'r')
    subplot(3,1,2); plot(t_func,FMI_b.MI_vel(2,:),'b',t_func,FMI_b.MI_acc(2,:),'r')
    subplot(3,1,3); plot(t_func,FMI_b.MI_vel(3,:),'b',t_func,FMI_b.MI_acc(3,:),'r')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t_func,FMI_strk.FI_vel(1,:),'b',t_func,FMI_strk.FI_acc(1,:),'r')
    subplot(3,1,2); plot(t_func,FMI_strk.FI_vel(2,:),'b',t_func,FMI_strk.FI_acc(2,:),'r')
    subplot(3,1,3); plot(t_func,FMI_strk.FI_vel(3,:),'b',t_func,FMI_strk.FI_acc(3,:),'r')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t_func,FMI_strk.MI_vel(1,:),'b',t_func,FMI_strk.MI_acc(1,:),'r')
    subplot(3,1,2); plot(t_func,FMI_strk.MI_vel(2,:),'b',t_func,FMI_strk.MI_acc(2,:),'r')
    subplot(3,1,3); plot(t_func,FMI_strk.MI_vel(3,:),'b',t_func,FMI_strk.MI_acc(3,:),'r')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t_func,FMI_b.LM(1,:))
    subplot(3,1,2); plot(t_func,FMI_b.LM(2,:))
    subplot(3,1,3); plot(t_func,FMI_b.LM(3,:))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t_func,FMI_b.AM(1,:))
    subplot(3,1,2); plot(t_func,FMI_b.AM(2,:))
    subplot(3,1,3); plot(t_func,FMI_b.AM(3,:))
    hold off
    
end
    
function x_dot = state_space_func(t,x,t_func,wing_kin,body_model,wing_model)

    qL_temp          = [interp1(t_func,wing_kin.qL(1,:),t,'spline'); ...
                        interp1(t_func,wing_kin.qL(2,:),t,'spline'); ...
                        interp1(t_func,wing_kin.qL(3,:),t,'spline'); ...
                        interp1(t_func,wing_kin.qL(4,:),t,'spline')];
                    
    qR_temp          = [interp1(t_func,wing_kin.qR(1,:),t,'spline'); ...
                        interp1(t_func,wing_kin.qR(2,:),t,'spline'); ...
                        interp1(t_func,wing_kin.qR(3,:),t,'spline'); ...
                        interp1(t_func,wing_kin.qR(4,:),t,'spline')];
                    
    wL_temp          = [interp1(t_func,wing_kin.wL(1,:),t,'spline'); ...
                        interp1(t_func,wing_kin.wL(2,:),t,'spline'); ...
                        interp1(t_func,wing_kin.wL(3,:),t,'spline')];
    
    wR_temp          = [interp1(t_func,wing_kin.wR(1,:),t,'spline'); ...
                        interp1(t_func,wing_kin.wR(2,:),t,'spline'); ...
                        interp1(t_func,wing_kin.wR(3,:),t,'spline')];
    
    w_dot_L_temp     = [interp1(t_func,wing_kin.w_dot_L(1,:),t,'spline'); ...
                        interp1(t_func,wing_kin.w_dot_L(2,:),t,'spline'); ...
                        interp1(t_func,wing_kin.w_dot_L(3,:),t,'spline')];
    
    w_dot_R_temp     = [interp1(t_func,wing_kin.w_dot_R(1,:),t,'spline'); ...
                        interp1(t_func,wing_kin.w_dot_R(2,:),t,'spline'); ...
                        interp1(t_func,wing_kin.w_dot_R(3,:),t,'spline')];
    
    wL_b_temp        = [interp1(t_func,wing_kin.wL_b(1,:),t,'spline'); ...
                        interp1(t_func,wing_kin.wL_b(2,:),t,'spline'); ...
                        interp1(t_func,wing_kin.wL_b(3,:),t,'spline')];
    
    wR_b_temp        = [interp1(t_func,wing_kin.wR_b(1,:),t,'spline'); ...
                        interp1(t_func,wing_kin.wR_b(2,:),t,'spline'); ...
                        interp1(t_func,wing_kin.wR_b(3,:),t,'spline')];
    
    w_dot_L_b_temp   = [interp1(t_func,wing_kin.w_dot_L_b(1,:),t,'spline'); ...
                        interp1(t_func,wing_kin.w_dot_L_b(2,:),t,'spline'); ...
                        interp1(t_func,wing_kin.w_dot_L_b(3,:),t,'spline')];
    
    w_dot_R_b_temp   = [interp1(t_func,wing_kin.w_dot_R_b(1,:),t,'spline'); ...
                        interp1(t_func,wing_kin.w_dot_R_b(2,:),t,'spline'); ...
                        interp1(t_func,wing_kin.w_dot_R_b(3,:),t,'spline')];
                        

    kine.qL          = qL_temp/norm(qL_temp);
    kine.qR          = qR_temp/norm(qR_temp);
    kine.RL          = quat2mat(kine.qL);
    kine.RR          = quat2mat(kine.qR);
    kine.vb          = x(1:3);
    kine.wb          = x(4:6);
    kine.wL          = wL_temp;
    kine.wR          = wR_temp;
    kine.wL_b        = wL_b_temp;
    kine.wR_b        = wR_b_temp;
    kine.w_dot_L     = w_dot_L_temp;
    kine.w_dot_R     = w_dot_R_temp;
    kine.w_dot_L_b   = w_dot_L_b_temp;
    kine.w_dot_R_b   = w_dot_R_b_temp;
    kine.R_strk      = wing_kin.R_strk;

    [ FMI_b, ~ ] = Inertia_instantaneous( kine, body_model, wing_model );
    
    M_mat = FMI_b.M_mat_b;
    FI    = FMI_b.F_I_vel;
    MI    = FMI_b.M_I_vel;
    
    x_dot = zeros(6,1);
    
    x_dot            = 1000*(M_mat\[ FI; MI]);

    
end

function [ FMI_b, FMI_strk ] = Inertia_FM( body_kin, wing_kin, body_model, wing_model, t, t_func )

    N = length(wing_kin.qL(1,:));

    vb          = body_kin.vb;
    wb          = body_kin.wb;
    ab          = body_kin.ab;
    w_dot_b     = body_kin.w_dot_b;
    
    FI_vel_b    = zeros(3,N);
    MI_vel_b    = zeros(3,N);
    FI_acc_b    = zeros(3,N);
    MI_acc_b    = zeros(3,N);
    
    FI_vel_strk = zeros(3,N);
    MI_vel_strk = zeros(3,N);
    FI_acc_strk = zeros(3,N);
    MI_acc_strk = zeros(3,N);
    
    LM          = zeros(3,N);
    AM          = zeros(3,N);
    KE          = zeros(1,N);
    KE_lin      = zeros(1,N);
    KE_ang      = zeros(1,N);
    
    vb_0        = zeros(3,N);
    wb_0        = zeros(3,N);
    
    
    for i = 1:N
        
        vb_temp         = [ interp1(t,vb(1,:),t_func(i),'spline'); ...
                            interp1(t,vb(2,:),t_func(i),'spline'); ...
                            interp1(t,vb(3,:),t_func(i),'spline') ];
                        
        wb_temp         = [ interp1(t,wb(1,:),t_func(i),'spline'); ...
                            interp1(t,wb(2,:),t_func(i),'spline'); ...
                            interp1(t,wb(3,:),t_func(i),'spline') ];
                        
        ab_temp         = [ interp1(t,ab(1,:),t_func(i),'spline'); ...
                            interp1(t,ab(2,:),t_func(i),'spline'); ...
                            interp1(t,ab(3,:),t_func(i),'spline') ];
                        
        w_dot_b_temp    = [ interp1(t,w_dot_b(1,:),t_func(i),'spline'); ...
                            interp1(t,w_dot_b(2,:),t_func(i),'spline'); ...
                            interp1(t,w_dot_b(3,:),t_func(i),'spline') ];

        kine.qL          = wing_kin.qL(:,i);
        kine.qR          = wing_kin.qR(:,i);
        kine.RL          = wing_kin.RL(:,:,i);
        kine.RR          = wing_kin.RR(:,:,i);
        kine.vb          = vb_temp;
        kine.wb          = wb_temp;
        kine.wL          = wing_kin.wL(:,i);
        kine.wR          = wing_kin.wR(:,i);
        kine.wL_b        = wing_kin.wL_b(:,i);
        kine.wR_b        = wing_kin.wR_b(:,i);
        kine.w_dot_L     = wing_kin.w_dot_L(:,i);
        kine.w_dot_R     = wing_kin.w_dot_R(:,i);
        kine.w_dot_L_b   = wing_kin.w_dot_L_b(:,i);
        kine.w_dot_R_b   = wing_kin.w_dot_R_b(:,i);
        kine.R_strk      = wing_kin.R_strk;    
        
        [ FM_b, FM_strk ] = Inertia_instantaneous( kine, body_model, wing_model );
        
        M_mat = FM_b.M_mat_b;
        
        FM_acc = 1e-3*M_mat*[ ab_temp; w_dot_b_temp ];
        
        FI_vel_b(:,i)       = FM_b.F_I_vel;
        MI_vel_b(:,i)       = FM_b.M_I_vel;
        FI_acc_b(:,i)       = FM_acc(1:3);
        MI_acc_b(:,i)       = FM_acc(4:6);

        FI_vel_strk(:,i)    = FM_strk.F_I_vel;
        MI_vel_strk(:,i)    = FM_strk.M_I_vel;
        FI_acc_strk(:,i)    = wing_kin.R_strk*FM_acc(1:3);
        MI_acc_strk(:,i)    = wing_kin.R_strk*FM_acc(4:6);

        LM(:,i)             = FM_b.LM;
        AM(:,i)             = FM_b.AM;
        KE(i)               = FM_b.KE;
        KE_lin(i)           = FM_b.KE_lin;
        KE_ang(i)           = FM_b.KE_ang;

        vb_0(:,i)           = FM_b.vb_0;
        wb_0(:,i)           = FM_b.wb_0;
        
    end
    
    FMI_b.FI_vel        = FI_vel_b;
    FMI_b.MI_vel        = MI_vel_b;
    FMI_b.FI_acc        = FI_acc_b;
    FMI_b.MI_acc        = MI_acc_b;
    FMI_b.LM            = LM;
    FMI_b.AM            = AM;
    FMI_b.KE            = KE;
    FMI_b.KE_lin        = KE_lin;
    FMI_b.KE_ang        = KE_ang;
    FMI_b.vb_0          = vb_0;
    FMI_b.wb_0          = wb_0;
    
    FMI_strk.FI_vel     = FI_vel_strk;
    FMI_strk.MI_vel     = MI_vel_strk;
    FMI_strk.FI_acc     = FI_acc_strk;
    FMI_strk.MI_acc     = MI_acc_strk;


end

