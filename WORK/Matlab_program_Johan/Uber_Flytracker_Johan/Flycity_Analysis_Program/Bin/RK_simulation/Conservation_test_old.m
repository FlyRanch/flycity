function Conservation_test_old( settings, pathDB )


    % Test the Inertia model for conservation of linear and angular
    % momentum:
    
    seq_nr      = 1;
    nr_points   = 100;
    nr_wb       = 3;
    
    % Setup the body_model and wing_model:
    
    body_model.mass_fly         = pathDB.body_model.mass_fly(seq_nr);
    body_model.mass_body        = pathDB.body_model.mass_body(seq_nr);
    body_model.Joint_left       = pathDB.body_model.Joint_left(seq_nr,:)';
    body_model.Joint_right      = pathDB.body_model.Joint_right(seq_nr,:)';
    body_model.cg_b             = pathDB.body_model.cg(seq_nr,:)';
    body_model.Inertia          = pathDB.body_model.Inertia(:,:,seq_nr);
    body_model.g                = 1e-3*settings.g;
    
    wing_model.virtual_mass     = pathDB.wing_model.virtual_mass(seq_nr);
    wing_model.wing_cg_L        = pathDB.wing_model.wing_cg_L(seq_nr,:)';
    wing_model.wing_cg_R        = pathDB.wing_model.wing_cg_R(seq_nr,:)';
    wing_model.virtual_Inertia  = pathDB.wing_model.virtual_Inertia(:,:,seq_nr);
    wing_model.y_sect_L         = pathDB.wing_model.y_sect_L(:,:,seq_nr)';
    wing_model.chords_L         = pathDB.wing_model.chords_L(seq_nr,:)';
    wing_model.y_sect_R         = pathDB.wing_model.y_sect_R(:,:,seq_nr)';
    wing_model.chords_R         = pathDB.wing_model.chords_R(seq_nr,:)';
    wing_model.rho              = settings.rho_air;
    
    f                           = pathDB.poly_fit.a_avg.f(seq_nr);
    dt                          = (1/f)/(nr_points-1);
    R_strk                      = pathDB.rot_mat.Rstr;
    
    % Create the wing kinematics:
    
    t               = nan(1,(nr_points*2)*nr_wb+1);
    theta_L         = nan(1,(nr_points*2)*nr_wb+1);
    eta_L           = nan(1,(nr_points*2)*nr_wb+1);
    phi_L           = nan(1,(nr_points*2)*nr_wb+1);
    theta_dot_L     = nan(1,(nr_points*2)*nr_wb+1);
    eta_dot_L       = nan(1,(nr_points*2)*nr_wb+1);
    phi_dot_L       = nan(1,(nr_points*2)*nr_wb+1);
    theta_ddot_L    = nan(1,(nr_points*2)*nr_wb+1);
    eta_ddot_L      = nan(1,(nr_points*2)*nr_wb+1);
    phi_ddot_L      = nan(1,(nr_points*2)*nr_wb+1);
    theta_R         = nan(1,(nr_points*2)*nr_wb+1);
    eta_R           = nan(1,(nr_points*2)*nr_wb+1);
    phi_R           = nan(1,(nr_points*2)*nr_wb+1);
    theta_dot_R     = nan(1,(nr_points*2)*nr_wb+1);
    eta_dot_R       = nan(1,(nr_points*2)*nr_wb+1);
    phi_dot_R       = nan(1,(nr_points*2)*nr_wb+1);
    theta_ddot_R    = nan(1,(nr_points*2)*nr_wb+1);
    eta_ddot_R      = nan(1,(nr_points*2)*nr_wb+1);
    phi_ddot_R      = nan(1,(nr_points*2)*nr_wb+1);
    RL              = nan(3,3,(nr_points*2)*nr_wb+1);
    RR              = nan(3,3,(nr_points*2)*nr_wb+1);
    wL              = nan(3,(nr_points*2)*nr_wb+1);
    wR              = nan(3,(nr_points*2)*nr_wb+1);
    w_dot_L         = nan(3,(nr_points*2)*nr_wb+1);
    w_dot_R         = nan(3,(nr_points*2)*nr_wb+1);
    wL_b            = nan(3,(nr_points*2)*nr_wb+1);
    wR_b            = nan(3,(nr_points*2)*nr_wb+1);
    w_dot_L_b       = nan(3,(nr_points*2)*nr_wb+1);
    w_dot_R_b       = nan(3,(nr_points*2)*nr_wb+1);
    wL_strk         = nan(3,(nr_points*2)*nr_wb+1);
    wR_strk         = nan(3,(nr_points*2)*nr_wb+1);
    w_dot_L_strk    = nan(3,(nr_points*2)*nr_wb+1);
    w_dot_R_strk    = nan(3,(nr_points*2)*nr_wb+1);
    
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

%     theta_L0      = (5*pi)/180;
%     eta_L0        = 0;
%     phi_L0        = (10*pi)/180;
%     theta_R0      = (5*pi)/180;
%     eta_R0        = 0;
%     phi_R0        = (10*pi)/180;
%     A_theta_L     = (10*pi)/180;
%     A_eta_L       = 0;
%     A_phi_L       = (60*pi)/180;
%     A_theta_R     = (10*pi)/180;
%     A_eta_R       = 0;
%     A_phi_R       = (60*pi)/180;

%     theta_L0      = 0;
%     eta_L0        = 0;
%     phi_L0        = (10*pi)/180;
%     theta_R0      = 0;
%     eta_R0        = 0;
%     phi_R0        = (10*pi)/180;
%     A_theta_L     = (10*pi)/180;
%     A_eta_L       = 0;
%     A_phi_L       = (60*pi)/180;
%     A_theta_R     = (10*pi)/180;
%     A_eta_R       = 0;
%     A_phi_R       = (60*pi)/180;
     
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
        a_fit.nr_points     = nr_points*2+1;
        a_fit.R_strk        = R_strk;

        [ kine ] = angular_velocities2( a_fit );
        
        if i == 1
            
            range_i         = ((i-1)*(2*nr_points)+1):(i*(2*nr_points)+1);
            range_temp      = 1:(2*nr_points+1);
            
            t(range_i)              = kine.t(range_temp);
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
            
            range_i         = ((i-1)*(2*nr_points)+1):(i*(2*nr_points)+1);
            range_temp      = 1:(2*nr_points+1);
            
            t(range_i)              = t((i-1)*(2*nr_points))+kine.dt+kine.t(range_temp);
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
            
            range_i         = ((i-1)*(2*nr_points)+1):(i*(2*nr_points)+1);
            
            t(range_i)              = t((i-1)*(2*nr_points))+kine.dt+kine.t;
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
    
    phi_b_0         = pi;
    theta_b_0       = -settings.beta_strk;
    
    qb_0_t          = [ sin(phi_b_0)*cos(theta_b_0); ...
                        cos(phi_b_0)*sin(theta_b_0); ...
                       -sin(phi_b_0)*cos(theta_b_0); ...
                        cos(phi_b_0)*cos(theta_b_0)];
                    
    qb_0            = qb_0_t/norm(qb_0_t);
    

    vb      = zeros(3,2*nr_points);
    wb      = zeros(3,2*nr_points);
    ab      = zeros(3,2*nr_points);
    w_dot_b = zeros(3,2*nr_points);
    FI_acc  = zeros(3,2*nr_points);
    FI_vel  = zeros(3,2*nr_points);
    MI_acc  = zeros(3,2*nr_points);
    MI_vel  = zeros(3,2*nr_points);
    LM_t    = zeros(3,2*nr_points);
    AM_t    = zeros(3,2*nr_points);
    
       
    for j = 1:(2*nr_points)
        
        kine.RL          = RL(:,:,j);
        kine.RR          = RR(:,:,j);
        kine.Rb          = quat2mat(qb_0);
        kine.vb          = vb(:,j);
        kine.wb          = wb(:,j);
        kine.ab          = ab(:,j);
        kine.w_dot_b     = w_dot_b(:,j);
        kine.wL          = wL(:,j);
        kine.wR          = wR(:,j);
        kine.wL_b        = wL_b(:,j);
        kine.wR_b        = wR_b(:,j);
        kine.w_dot_L     = w_dot_L(:,j);
        kine.w_dot_R     = w_dot_R(:,j);
        kine.w_dot_L_b   = w_dot_L_b(:,j);
        kine.w_dot_R_b   = w_dot_R_b(:,j);
        kine.R_strk      = R_strk;

        [ FM_b, ~ ] = Inertia_instantaneous( kine, body_model, wing_model );

        FI_acc(:,j) = FM_b.F_I_acc;
        MI_acc(:,j) = FM_b.M_I_acc;
        FI_vel(:,j) = FM_b.F_I_vel;
        MI_vel(:,j) = FM_b.M_I_vel;
        LM_t(:,j)   = FM_b.LM;
        AM_t(:,j)   = FM_b.AM;
        
    end
    
    FI_mean = mean(FI_vel,2)
    MI_mean = mean(MI_vel,2)

    
    % Wing kinematics:
    
    wing_kin.R_strk     = R_strk;
    wing_kin.RL         = RL;
    wing_kin.RR         = RR;
    wing_kin.wL         = wL;
    wing_kin.wR         = wR;
    wing_kin.wL_b       = wL_b;
    wing_kin.wR_b       = wR_b;
    wing_kin.w_dot_L_b  = w_dot_L_b;
    wing_kin.w_dot_R_b  = w_dot_R_b;
    

    % Initial conditions:
    
    
    state.id        = 1;
    state.dt        = dt;
    state.xyz       = [0; 0; 0];
    state.qb        = qb_0;
    state.Rb        = quat2mat(qb_0);
    state.v_system  = [0; 0; 0];
    state.vb        = [0; 0; 0];
    state.wb        = [0; 0; 0];
    state.ab        = [0; 0; 0];
    state.w_dot_b   = [0; 0; 0];
    
%     [ state_new, ~ ] = RK_updater3( state, body_model, wing_model, wing_kin );
%     
%     state.id        = 1;
%     state.dt        = dt;
%     state.vb        = state_new.vb_0;
%     state.wb        = state_new.wb_0;
%     state.ab        = state_new.ab_old;
%     state.w_dot_b   = state_new.w_dot_b_old;
    
    state.xyz
    state.qb
    state.Rb
    state.v_system
    state.vb
    state.wb
    state.ab
    state.w_dot_b

    xyz     = zeros(3,nr_points*nr_wb);
    qb      = zeros(4,nr_points*nr_wb);
    Rb      = zeros(3,3,nr_points*nr_wb);
    v_system= zeros(3,nr_points*nr_wb);
    vb      = zeros(3,nr_points*nr_wb);
    wb      = zeros(3,nr_points*nr_wb);
    ab      = zeros(3,nr_points*nr_wb);
    w_dot_b = zeros(3,nr_points*nr_wb);
    vb_0    = zeros(3,nr_points*nr_wb);
    wb_0    = zeros(3,nr_points*nr_wb);
    FI_acc  = zeros(3,nr_points*nr_wb);
    FI_vel  = zeros(3,nr_points*nr_wb);
    MI_acc  = zeros(3,nr_points*nr_wb);
    MI_vel  = zeros(3,nr_points*nr_wb);
    LM      = zeros(3,nr_points*nr_wb);
    AM      = zeros(3,nr_points*nr_wb);
    KE      = zeros(1,nr_points*nr_wb);
    
    for j = 1:(nr_points*nr_wb)
    
        xyz(:,j)        = state.xyz;
        qb(:,j)         = state.qb;
        Rb(:,:,j)       = state.Rb;
        v_system(:,j)   = state.v_system;
        vb(:,j)         = state.vb;
        wb(:,j)         = state.wb;
        ab(:,j)         = state.ab;
        w_dot_b(:,j)    = state.w_dot_b;
        
        [ state_new, FM_old ] = RK_updater3( state, body_model, wing_model, wing_kin );
        
        state.id        = state_new.id;
        state.dt        = state_new.dt;
        state.xyz       = state_new.xyz;
        state.qb        = state_new.qb;
        state.Rb        = state_new.Rb;
        state.v_system  = state_new.v_system;
        state.vb        = state_new.vb;
        state.wb        = state_new.wb;
        state.ab        = state_new.ab;
        state.w_dot_b   = state_new.w_dot_b;
                
        FI_acc(:,j)     = FM_old.FI_acc_b;
        FI_vel(:,j)     = FM_old.FI_vel_b;
        MI_acc(:,j)     = FM_old.MI_acc_b;
        MI_vel(:,j)     = FM_old.MI_vel_b;
        
        LM(:,j)         = FM_old.LM;
        AM(:,j)         = FM_old.AM;
        
        vb_0(:,j)       = state_new.vb_0;
        wb_0(:,j)       = state_new.wb_0;
        
        KE(j)           = FM_old.KE;

    end
        
    figure()
    plot(KE)
    
    figure()
    hold on
    subplot(3,1,1); plot(xyz(1,:))
    subplot(3,1,2); plot(xyz(2,:))
    subplot(3,1,3); plot(xyz(3,:))
    hold off
    
    figure()
    hold on
    subplot(4,1,1); plot(qb(1,:))
    subplot(4,1,2); plot(qb(2,:))
    subplot(4,1,3); plot(qb(3,:))
    subplot(4,1,4); plot(qb(4,:))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(vb(1,:))
    title('vb')
    subplot(3,1,2); plot(vb(2,:))
    subplot(3,1,3); plot(vb(3,:))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(wb(1,:))
    title('wb')
    subplot(3,1,2); plot(wb(2,:))
    subplot(3,1,3); plot(wb(3,:))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(ab(1,:))
    title('ab')
    subplot(3,1,2); plot(ab(2,:))
    subplot(3,1,3); plot(ab(3,:))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(w_dot_b(1,:))
    title('w dot b')
    subplot(3,1,2); plot(w_dot_b(2,:))
    subplot(3,1,3); plot(w_dot_b(3,:))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(1:(nr_points*nr_wb),FI_acc(1,:),'r',1:(nr_points*nr_wb),FI_vel(1,:),'b',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(FI_vel(1,1:nr_points)))
    title('FI')
    subplot(3,1,2); plot(1:(nr_points*nr_wb),FI_acc(2,:),'r',1:(nr_points*nr_wb),FI_vel(2,:),'b',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(FI_vel(2,1:nr_points)))
    subplot(3,1,3); plot(1:(nr_points*nr_wb),FI_acc(3,:),'r',1:(nr_points*nr_wb),FI_vel(3,:),'b',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(FI_vel(3,1:nr_points)))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(1:(nr_points*nr_wb),MI_acc(1,:),'r',1:(nr_points*nr_wb),MI_vel(1,:),'b',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(MI_vel(1,1:nr_points)))
    title('MI')
    subplot(3,1,2); plot(1:(nr_points*nr_wb),MI_acc(2,:),'r',1:(nr_points*nr_wb),MI_vel(2,:),'b',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(MI_vel(2,1:nr_points)))
    subplot(3,1,3); plot(1:(nr_points*nr_wb),MI_acc(3,:),'r',1:(nr_points*nr_wb),MI_vel(3,:),'b',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(MI_vel(3,1:nr_points)))
    hold off
    
    
    figure()
    hold on
    subplot(3,1,1); plot(1:(nr_points*nr_wb),LM(1,:),'r',1:(nr_points*nr_wb),LM(1,:),'b',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(LM(1,1:nr_points)))
    title('LM')
    subplot(3,1,2); plot(1:(nr_points*nr_wb),LM(2,:),'r',1:(nr_points*nr_wb),LM(2,:),'b',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(LM(2,1:nr_points)))
    subplot(3,1,3); plot(1:(nr_points*nr_wb),LM(3,:),'r',1:(nr_points*nr_wb),LM(3,:),'b',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(LM(3,1:nr_points)))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(1:(nr_points*nr_wb),AM(1,:),'r',1:(nr_points*nr_wb),AM(1,:),'b',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(AM(1,1:nr_points)))
    title('AM')
    subplot(3,1,2); plot(1:(nr_points*nr_wb),AM(2,:),'r',1:(nr_points*nr_wb),AM(2,:),'b',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(AM(2,1:nr_points)))
    subplot(3,1,3); plot(1:(nr_points*nr_wb),AM(3,:),'r',1:(nr_points*nr_wb),AM(3,:),'b',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(AM(3,1:nr_points)))
    hold off

    vb_mean = mean(vb(:,1:(nr_points)),2)
    wb_mean = mean(wb(:,1:(nr_points)),2)
    ab_mean = mean(ab(:,1:(nr_points)),2)
    w_dot_b_mean = mean(w_dot_b(:,1:(nr_points)),2)
    FI_vel_mean = mean(FI_vel(:,1:(nr_points)),2)
    MI_vel_mean = mean(MI_vel(:,1:(nr_points)),2)
    LM_mean = mean(LM(:,1:(nr_points)),2)
    AM_mean = mean(AM(:,1:(nr_points)),2)

    
end

