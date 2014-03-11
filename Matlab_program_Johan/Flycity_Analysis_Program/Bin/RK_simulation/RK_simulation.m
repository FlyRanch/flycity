function RK_simulation(settings,pathDB)

    % Create the wing kinematics:
    
    seq_nr      = 11;
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
    
    nr_sect                     = length(wing_model.chords_L);
    
    f                           = pathDB.poly_fit.a_avg.f(seq_nr);
    down_up_avg                 = pathDB.poly_fit.a_avg.down_up(seq_nr);
    dt                          = (1/f)/(nr_points-1);
    R_strk                      = pathDB.rot_mat.Rstr;
    
    % Create the wing kinematics:
    
    n_pol_theta     = settings.n_pol_theta;
    n_pol_eta       = settings.n_pol_eta;
    n_pol_phi       = settings.n_pol_phi;
    
    a_theta_L       = nan(n_pol_theta*2+2,nr_wb);
    a_eta_L         = nan(n_pol_eta*2+2,nr_wb);
    a_phi_L         = nan(n_pol_phi*2+2,nr_wb);
    a_theta_R       = nan(n_pol_theta*2+2,nr_wb);
    a_eta_R         = nan(n_pol_eta*2+2,nr_wb);
    a_phi_R         = nan(n_pol_phi*2+2,nr_wb);
    
%     a_glob_theta    = pathDB.poly_fit.a_glob.theta;
%     a_glob_eta      = pathDB.poly_fit.a_glob.eta;
%     a_glob_phi      = pathDB.poly_fit.a_glob.phi;
%     down_up_glob    = pathDB.poly_fit.a_glob.down_up;

    a_avg_theta    = pathDB.poly_fit.a_avg.theta_LR(:,seq_nr);
    a_avg_eta      = pathDB.poly_fit.a_avg.eta_LR(:,seq_nr);
    a_avg_phi      = pathDB.poly_fit.a_avg.phi_LR(:,seq_nr);
    down_up_avg    = pathDB.poly_fit.a_avg.down_up(seq_nr);
    
    down_up         = nan(nr_wb,1);
    
%     for k = 1:nr_wb
%         
%         a_theta_L(:,k)  = a_glob_theta;
%         a_eta_L(:,k)    = a_glob_eta;
%         a_phi_L(:,k)    = a_glob_phi;
%         a_theta_R(:,k)  = a_glob_theta;
%         a_eta_R(:,k)    = a_glob_eta;
%         a_phi_R(:,k)    = a_glob_phi;
%         down_up(k)      = down_up_glob;
%         
%     end

    for k = 1:nr_wb
        
        a_theta_L(:,k)  = a_avg_theta;
        a_eta_L(:,k)    = a_avg_eta;
        a_phi_L(:,k)    = a_avg_phi;
        a_theta_R(:,k)  = a_avg_theta;
        a_eta_R(:,k)    = a_avg_eta;
        a_phi_R(:,k)    = a_avg_phi;
        down_up(k)      = down_up_avg;
        
    end
    
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
    
    
    for i = 1:nr_wb
        
        a_fit.a_theta_L     = a_theta_L(:,i);
        a_fit.a_eta_L       = a_eta_L(:,i);
        a_fit.a_phi_L       = a_phi_L(:,i);
        a_fit.a_theta_R     = a_theta_R(:,i);
        a_fit.a_eta_R       = a_eta_R(:,i);
        a_fit.a_phi_R       = a_phi_R(:,i);
        a_fit.f             = f;
        a_fit.down_up       = down_up(i);
        a_fit.nr_points     = nr_points*2+1;
        a_fit.R_strk        = R_strk;

        [ kine ] = angular_velocities( a_fit );
        
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

        [ FM_b, ~ ] = Inertia_instantaneous( kine, body_model, wing_model);

        FI_acc(:,j) = FM_b.F_I_acc;
        MI_acc(:,j) = FM_b.M_I_acc;
        FI_vel(:,j) = FM_b.F_I_vel;
        MI_vel(:,j) = FM_b.M_I_vel;
        
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
    
    wb_loc      = zeros(nr_wb,2);
    down_loc    = zeros(nr_wb,2);
    up_loc      = zeros(nr_wb,2);
    
    for k = 1:nr_wb
            
        wb_loc(k,:)      = [ ((k-1)*2*nr_points)+1 ((k-1)*2*nr_points)+2*nr_points];
        down_loc(k,:)    = [ ((k-1)*2*nr_points)+1 ((k-1)*2*nr_points)+round(down_up(k)*2*nr_points)];
        up_loc(k,:)      = [ ((k-1)*2*nr_points)+round(down_up(k)*2*nr_points)+1 ((k-1)*2*nr_points)+2*nr_points];

    end
    
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
    state.alfa_L_old = zeros(nr_sect,1);
    state.alfa_R_old = zeros(nr_sect,1);
    
%     [ state_new, ~ ] = RK_updater( state, body_model, wing_model, wing_kin );
%     
%     state.id        = 1;
%     state.dt        = dt;
%     state.vb        = state_new.vb_0;
%     state.wb        = state_new.wb_0;
%     state.ab        = state_new.ab_old;
%     state.w_dot_b   = state_new.w_dot_b_old;
%     
%     state.xyz
%     state.qb
%     state.Rb
%     state.v_system
%     state.vb
%     state.wb
%     state.ab
%     state.w_dot_b

    state.id
    state.dt
    state.xyz
    state.qb
    state.Rb
    state.v_system
    state.vb
    state.wb
    state.ab
    state.w_dot_b
    state.alfa_L_old
    state.alfa_R_old

    xyz         = zeros(3,nr_points*nr_wb);
    qb          = zeros(4,nr_points*nr_wb);
    Rb          = zeros(3,3,nr_points*nr_wb);
    v_system    = zeros(3,nr_points*nr_wb);
    vb          = zeros(3,nr_points*nr_wb);
    wb          = zeros(3,nr_points*nr_wb);
    ab          = zeros(3,nr_points*nr_wb);
    w_dot_b     = zeros(3,nr_points*nr_wb);
    vb_0        = zeros(3,nr_points*nr_wb);
    wb_0        = zeros(3,nr_points*nr_wb);
    FI_acc_b    = zeros(3,nr_points*nr_wb);
    FI_vel_b    = zeros(3,nr_points*nr_wb);
    MI_acc_b    = zeros(3,nr_points*nr_wb);
    MI_vel_b    = zeros(3,nr_points*nr_wb);
    FA_b        = zeros(3,nr_points*nr_wb);
    MA_b        = zeros(3,nr_points*nr_wb);
    Fg_b        = zeros(3,nr_points*nr_wb);
    FI_acc_strk = zeros(3,nr_points*nr_wb);
    FI_vel_strk = zeros(3,nr_points*nr_wb);
    MI_acc_strk = zeros(3,nr_points*nr_wb);
    MI_vel_strk = zeros(3,nr_points*nr_wb);
    FA_strk     = zeros(3,nr_points*nr_wb);
    MA_strk     = zeros(3,nr_points*nr_wb);
    Fg_strk     = zeros(3,nr_points*nr_wb);
    LM          = zeros(3,nr_points*nr_wb);
    AM          = zeros(3,nr_points*nr_wb);
    KE          = zeros(1,nr_points*nr_wb);
    
    wb_count = 1;
    
    for j = 1:(nr_points*nr_wb)
    
        xyz(:,j)        = state.xyz;
        qb(:,j)         = state.qb;
        Rb(:,:,j)       = state.Rb;
        v_system(:,j)   = state.v_system;
        vb(:,j)         = state.vb;
        wb(:,j)         = state.wb;
        ab(:,j)         = state.ab;
        w_dot_b(:,j)    = state.w_dot_b;
        
                   
        wingbeat.wb_loc     = wb_loc(wb_count,:);
        wingbeat.down_loc   = down_loc(wb_count,:);
        wingbeat.up_loc     = up_loc(wb_count,:);

        if mod(j,nr_points) == 0
            
            wb_count = wb_count+1;
            
        end
        
        [ state_new, FM_old ] = RK_updater( state, body_model, wing_model, wing_kin, wingbeat );
        
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
        state.alfa_L_old = state_new.alfa_L_old;
        state.alfa_R_old = state_new.alfa_R_old;
        
        FI_acc_b(:,j)    = FM_old.FI_acc_b;
        FI_vel_b(:,j)    = FM_old.FI_vel_b;
        MI_acc_b(:,j)    = FM_old.MI_acc_b;
        MI_vel_b(:,j)    = FM_old.MI_vel_b;
        FA_b(:,j)        = FM_old.FA_b;
        MA_b(:,j)        = FM_old.MA_b;
        Fg_b(:,j)        = FM_old.Fg_b;
        FI_acc_strk(:,j) = FM_old.FI_acc_strk;
        FI_vel_strk(:,j) = FM_old.FI_vel_strk;
        MI_acc_strk(:,j) = FM_old.MI_acc_strk;
        MI_vel_strk(:,j) = FM_old.MI_vel_strk;
        FA_strk(:,j)     = FM_old.FA_strk;
        MA_strk(:,j)     = FM_old.MA_strk;
        Fg_strk(:,j)     = FM_old.Fg_strk;
        
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
    
    % Body reference frame:
    
    figure()
    hold on
    subplot(3,1,1); plot(1:(nr_points*nr_wb),FI_acc_b(1,:),'r',1:(nr_points*nr_wb),FI_vel_b(1,:),'b',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(FI_vel_b(1,1:nr_points)))
    title('FI')
    subplot(3,1,2); plot(1:(nr_points*nr_wb),FI_acc_b(2,:),'r',1:(nr_points*nr_wb),FI_vel_b(2,:),'b',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(FI_vel_b(2,1:nr_points)))
    subplot(3,1,3); plot(1:(nr_points*nr_wb),FI_acc_b(3,:),'r',1:(nr_points*nr_wb),FI_vel_b(3,:),'b',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(FI_vel_b(3,1:nr_points)))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(1:(nr_points*nr_wb),MI_acc_b(1,:),'r',1:(nr_points*nr_wb),MI_vel_b(1,:),'b',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(MI_vel_b(1,1:nr_points)))
    title('MI')
    subplot(3,1,2); plot(1:(nr_points*nr_wb),MI_acc_b(2,:),'r',1:(nr_points*nr_wb),MI_vel_b(2,:),'b',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(MI_vel_b(2,1:nr_points)))
    subplot(3,1,3); plot(1:(nr_points*nr_wb),MI_acc_b(3,:),'r',1:(nr_points*nr_wb),MI_vel_b(3,:),'b',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(MI_vel_b(3,1:nr_points)))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(1:(nr_points*nr_wb),FA_b(1,:),'r',1:(nr_points*nr_wb),FI_vel_b(1,:),'b',1:(nr_points*nr_wb),Fg_b(1,:),'g',1:(nr_points*nr_wb),FI_vel_b(1,:)+FA_b(1,:)+Fg_b(1,:),'k',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(FI_vel_b(1,1:nr_points)+FA_b(1,1:nr_points)+Fg_b(1,1:nr_points)))
    title('FA')
    subplot(3,1,2); plot(1:(nr_points*nr_wb),FA_b(2,:),'r',1:(nr_points*nr_wb),FI_vel_b(2,:),'b',1:(nr_points*nr_wb),Fg_b(2,:),'g',1:(nr_points*nr_wb),FI_vel_b(2,:)+FA_b(2,:)+Fg_b(2,:),'k',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(FI_vel_b(2,1:nr_points)+FA_b(2,1:nr_points)+Fg_b(2,1:nr_points)))
    subplot(3,1,3); plot(1:(nr_points*nr_wb),FA_b(3,:),'r',1:(nr_points*nr_wb),FI_vel_b(3,:),'b',1:(nr_points*nr_wb),Fg_b(3,:),'g',1:(nr_points*nr_wb),FI_vel_b(3,:)+FA_b(3,:)+Fg_b(3,:),'k',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(FI_vel_b(3,1:nr_points)+FA_b(3,1:nr_points)+Fg_b(3,1:nr_points)))
    
    figure()
    hold on
    subplot(3,1,1); plot(1:(nr_points*nr_wb),MA_b(1,:),'r',1:(nr_points*nr_wb),MI_vel_b(1,:),'b',1:(nr_points*nr_wb),MI_vel_b(1,:)+MA_b(1,:),'k',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(MI_vel_b(1,1:nr_points)+MA_b(1,1:nr_points)))
    title('MA')
    subplot(3,1,2); plot(1:(nr_points*nr_wb),MA_b(2,:),'r',1:(nr_points*nr_wb),MI_vel_b(2,:),'b',1:(nr_points*nr_wb),MI_vel_b(2,:)+MA_b(2,:),'k',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(MI_vel_b(2,1:nr_points)+MA_b(2,1:nr_points)))
    subplot(3,1,3); plot(1:(nr_points*nr_wb),MA_b(3,:),'r',1:(nr_points*nr_wb),MI_vel_b(3,:),'b',1:(nr_points*nr_wb),MI_vel_b(3,:)+MA_b(3,:),'k',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(MI_vel_b(3,1:nr_points)+MA_b(3,1:nr_points)))
    hold off
    
    % Strokeplane reference frame:
    
    figure()
    hold on
    subplot(3,1,1); plot(1:(nr_points*nr_wb),FI_acc_strk(1,:),'r',1:(nr_points*nr_wb),FI_vel_strk(1,:),'b',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(FI_vel_strk(1,1:nr_points)))
    title('FI')
    subplot(3,1,2); plot(1:(nr_points*nr_wb),FI_acc_strk(2,:),'r',1:(nr_points*nr_wb),FI_vel_strk(2,:),'b',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(FI_vel_strk(2,1:nr_points)))
    subplot(3,1,3); plot(1:(nr_points*nr_wb),FI_acc_strk(3,:),'r',1:(nr_points*nr_wb),FI_vel_strk(3,:),'b',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(FI_vel_strk(3,1:nr_points)))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(1:(nr_points*nr_wb),MI_acc_strk(1,:),'r',1:(nr_points*nr_wb),MI_vel_strk(1,:),'b',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(MI_vel_strk(1,1:nr_points)))
    title('MI')
    subplot(3,1,2); plot(1:(nr_points*nr_wb),MI_acc_strk(2,:),'r',1:(nr_points*nr_wb),MI_vel_strk(2,:),'b',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(MI_vel_strk(2,1:nr_points)))
    subplot(3,1,3); plot(1:(nr_points*nr_wb),MI_acc_strk(3,:),'r',1:(nr_points*nr_wb),MI_vel_strk(3,:),'b',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(MI_vel_strk(3,1:nr_points)))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(1:(nr_points*nr_wb),FA_strk(1,:),'r',1:(nr_points*nr_wb),FI_vel_strk(1,:),'b',1:(nr_points*nr_wb),Fg_strk(1,:),'g',1:(nr_points*nr_wb),FI_vel_strk(1,:)+FA_strk(1,:)+Fg_strk(1,:),'k',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(FI_vel_strk(1,1:nr_points)+FA_strk(1,1:nr_points)+Fg_strk(1,1:nr_points)))
    title('FA')
    subplot(3,1,2); plot(1:(nr_points*nr_wb),FA_strk(2,:),'r',1:(nr_points*nr_wb),FI_vel_strk(2,:),'b',1:(nr_points*nr_wb),Fg_strk(2,:),'g',1:(nr_points*nr_wb),FI_vel_strk(2,:)+FA_strk(2,:)+Fg_strk(2,:),'k',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(FI_vel_strk(2,1:nr_points)+FA_strk(2,1:nr_points)+Fg_strk(2,1:nr_points)))
    subplot(3,1,3); plot(1:(nr_points*nr_wb),FA_strk(3,:),'r',1:(nr_points*nr_wb),FI_vel_strk(3,:),'b',1:(nr_points*nr_wb),Fg_strk(3,:),'g',1:(nr_points*nr_wb),FI_vel_strk(3,:)+FA_strk(3,:)+Fg_strk(3,:),'k',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(FI_vel_strk(3,1:nr_points)+FA_strk(3,1:nr_points)+Fg_strk(3,1:nr_points)))
    
    figure()
    hold on
    subplot(3,1,1); plot(1:(nr_points*nr_wb),MA_strk(1,:),'r',1:(nr_points*nr_wb),MI_vel_strk(1,:),'b',1:(nr_points*nr_wb),MI_vel_strk(1,:)+MA_strk(1,:),'k',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(MI_vel_strk(1,1:nr_points)+MA_strk(1,1:nr_points)))
    title('MA')
    subplot(3,1,2); plot(1:(nr_points*nr_wb),MA_strk(2,:),'r',1:(nr_points*nr_wb),MI_vel_strk(2,:),'b',1:(nr_points*nr_wb),MI_vel_strk(2,:)+MA_strk(2,:),'k',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(MI_vel_strk(2,1:nr_points)+MA_strk(2,1:nr_points)))
    subplot(3,1,3); plot(1:(nr_points*nr_wb),MA_strk(3,:),'r',1:(nr_points*nr_wb),MI_vel_strk(3,:),'b',1:(nr_points*nr_wb),MI_vel_strk(3,:)+MA_strk(3,:),'k',1:(nr_points*nr_wb),ones(nr_points*nr_wb,1)*mean(MI_vel_strk(3,1:nr_points)+MA_strk(3,1:nr_points)))
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

end

