function Conservation_test( settings, pathDB )


    % Test the Inertia model for conservation of linear and angular
    % momentum:
    
    seq_nr      = 1;
    nr_points   = 200;
    nr_wb       = 10;
    
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
    
    nr_sect                     = length(wing_model.y_sect_L(1,:));
    
    a_avg_theta_L   = pathDB.poly_fit.a_avg.theta_LR(:,seq_nr);
    a_avg_eta_L     = pathDB.poly_fit.a_avg.eta_LR(:,seq_nr);
    a_avg_phi_L     = pathDB.poly_fit.a_avg.phi_LR(:,seq_nr);  
    
    a_avg_theta_R   = pathDB.poly_fit.a_avg.theta_LR(:,seq_nr);
    a_avg_eta_R     = pathDB.poly_fit.a_avg.eta_LR(:,seq_nr);
    a_avg_phi_R     = pathDB.poly_fit.a_avg.phi_LR(:,seq_nr);  
   
    down_up         = pathDB.poly_fit.a_avg.down_up(seq_nr);
    
    f               = pathDB.poly_fit.a_avg.f(seq_nr);
    
    dt              = (1/f)/(nr_points-1);
    
    R_strk = pathDB.rot_mat.Rstr;
    
    n_pol_theta     = (length(a_avg_theta_L)-2)/2;
    n_pol_eta       = (length(a_avg_eta_L)-2)/2;
    n_pol_phi       = (length(a_avg_phi_L)-2)/2;
    
    % Set up the polynomial coefficients per wingbeat:
    
    a_theta_L       = nan(2*(n_pol_theta+1),nr_wb);
    a_eta_L         = nan(2*(n_pol_eta+1),nr_wb);
    a_phi_L         = nan(2*(n_pol_phi+1),nr_wb);
    a_theta_R       = nan(2*(n_pol_theta+1),nr_wb);
    a_eta_R         = nan(2*(n_pol_eta+1),nr_wb);
    a_phi_R         = nan(2*(n_pol_phi+1),nr_wb);
    
    for i = 1:nr_wb
        
        a_theta_L(:,i)  = a_avg_theta_L;
        a_eta_L(:,i)    = a_avg_eta_L;
        a_phi_L(:,i)    = a_avg_phi_L;
        a_theta_R(:,i)  = a_avg_theta_R;
        a_eta_R(:,i)    = a_avg_eta_R;
        a_phi_R(:,i)    = a_avg_phi_R;
        
    end
    
        
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
    
    wb_loc          = nan(nr_wb,2);
    down_loc        = nan(nr_wb,2);
    up_loc          = nan(nr_wb,2);
    
        
    for i = 1:nr_wb
        
        a_fit.a_theta_L     = a_theta_L(:,i);
        a_fit.a_eta_L       = a_eta_L(:,i);
        a_fit.a_phi_L       = a_phi_L(:,i);
        a_fit.a_theta_R     = a_theta_R(:,i);
        a_fit.a_eta_R       = a_eta_R(:,i);
        a_fit.a_phi_R       = a_phi_R(:,i);
        a_fit.f             = f;
        a_fit.down_up       = down_up;
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
        
        wb_loc(i,:)         = [ ((i-1)*nr_points+1) (i*nr_points) ];
        down_loc(i,:)       = [ ((i-1)*nr_points+1) (((i-1)*nr_points)+round(down_up*nr_points)) ];
        up_loc(i,:)         = [ (((i-1)*nr_points)+round(down_up*nr_points)+1) (i*nr_points) ];
            
        clear kine a_fit
        
    end
    
    figure()
    hold on
    subplot(3,1,1); plot(t,wL(1,:),'r',t,wR(1,:),'b',t,ones(nr_points*2*nr_wb+1,1)*mean(wL(1,1:(end-1))),'r',t,ones(nr_points*2*nr_wb+1,1)*mean(wR(1,1:(end-1))),'b')
    subplot(3,1,2); plot(t,wL(2,:),'r',t,wR(2,:),'b',t,ones(nr_points*2*nr_wb+1,1)*mean(wL(2,1:(end-1))),'r',t,ones(nr_points*2*nr_wb+1,1)*mean(wR(2,1:(end-1))),'b')
    subplot(3,1,3); plot(t,wL(3,:),'r',t,wR(3,:),'b',t,ones(nr_points*2*nr_wb+1,1)*mean(wL(3,1:(end-1))),'r',t,ones(nr_points*2*nr_wb+1,1)*mean(wR(3,1:(end-1))),'b')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t,wL_b(1,:),'r',t,wR_b(1,:),'b',t,ones(nr_points*2*nr_wb+1,1)*mean(wL_b(1,1:(end-1))),'r',t,ones(nr_points*2*nr_wb+1,1)*mean(wR_b(1,1:(end-1))),'b')
    subplot(3,1,2); plot(t,wL_b(2,:),'r',t,wR_b(2,:),'b',t,ones(nr_points*2*nr_wb+1,1)*mean(wL_b(2,1:(end-1))),'r',t,ones(nr_points*2*nr_wb+1,1)*mean(wR_b(2,1:(end-1))),'b')
    subplot(3,1,3); plot(t,wL_b(3,:),'r',t,wR_b(3,:),'b',t,ones(nr_points*2*nr_wb+1,1)*mean(wL_b(3,1:(end-1))),'r',t,ones(nr_points*2*nr_wb+1,1)*mean(wR_b(3,1:(end-1))),'b')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t,wL_strk(1,:),'r',t,wR_strk(1,:),'b',t,ones(nr_points*2*nr_wb+1,1)*mean(wL_strk(1,1:(end-1))),'r',t,ones(nr_points*2*nr_wb+1,1)*mean(wR_strk(1,1:(end-1))),'b')
    subplot(3,1,2); plot(t,wL_strk(2,:),'r',t,wR_strk(2,:),'b',t,ones(nr_points*2*nr_wb+1,1)*mean(wL_strk(2,1:(end-1))),'r',t,ones(nr_points*2*nr_wb+1,1)*mean(wR_strk(2,1:(end-1))),'b')
    subplot(3,1,3); plot(t,wL_strk(3,:),'r',t,wR_strk(3,:),'b',t,ones(nr_points*2*nr_wb+1,1)*mean(wL_strk(3,1:(end-1))),'r',t,ones(nr_points*2*nr_wb+1,1)*mean(wR_strk(3,1:(end-1))),'b')
    hold off
    
    FI_acc = nan(3,nr_points*nr_wb);
    FI_vel = nan(3,nr_points*nr_wb);
    MI_acc = nan(3,nr_points*nr_wb);
    MI_vel = nan(3,nr_points*nr_wb);
    
    vb      = nan(3,nr_points*nr_wb);
    wb      = nan(3,nr_points*nr_wb);
    ab      = nan(3,nr_points*nr_wb);
    w_dot_b = nan(3,nr_points*nr_wb);
    
    LM      = nan(3,nr_points*nr_wb);
    AM      = nan(3,nr_points*nr_wb);
    
    vb_0    = nan(3,nr_points*nr_wb);
    wb_0    = nan(3,nr_points*nr_wb);
    
    RcgLwL  = nan(3,nr_points*nr_wb);
    RcgRwR  = nan(3,nr_points*nr_wb);
    
    % Wing kinematics:
    
    wing_kin.R_strk     = R_strk;
    wing_kin.RL         = RL;
    wing_kin.RR         = RR;
    wing_kin.wL_b       = wL_b;
    wing_kin.wR_b       = wR_b;
    wing_kin.w_dot_L_b  = w_dot_L_b;
    wing_kin.w_dot_R_b  = w_dot_R_b;
    
    % Initial conditions:
    
    state.id        = 1;
    state.dt        = dt;
    state.vb        = [0; 0; 0];
    state.wb        = [0; 0; 0];
    state.ab        = [0; 0; 0];
    state.w_dot_b   = [0; 0; 0];
    
    [ state_new, ~ ] = RK_updater( state, body_model, wing_model, wing_kin );
    
    state.id        = 1;
    state.dt        = dt;
    state.vb        = state_new.vb_0;
    state.wb        = state_new.wb_0;
    state.ab        = state_new.ab;
    state.w_dot_b   = state_new.w_dot_b;

    
    % RK4-simulation:
    
    for j = 1:(nr_points*nr_wb)
        
        
        vb(:,j)         = state.vb;
        wb(:,j)         = state.wb;
        ab(:,j)         = state.ab;
        w_dot_b(:,j)    = state.w_dot_b;
        
        [ state_new, FM_old ] = RK_updater( state, body_model, wing_model, wing_kin );
        
        state.id        = state_new.id;
        state.dt        = state_new.dt;
        state.vb        = state_new.vb;
        state.wb        = state_new.wb;
        state.ab        = state_new.ab;
        state.w_dot_b   = state_new.w_dot_b;
                
        FI_acc(:,j) = FM_old.FI_acc_b;
        FI_vel(:,j) = FM_old.FI_vel_b;
        MI_acc(:,j) = FM_old.MI_acc_b;
        MI_vel(:,j) = FM_old.MI_vel_b;
        
        LM(:,j)       = FM_old.LM;
        AM(:,j)       = FM_old.AM;
        
        vb_0(:,j)     = state_new.vb_0;
        wb_0(:,j)     = state_new.wb_0;
        
        RcgLwL(:,j)   = state_new.RcgLwL;
        RcgRwR(:,j)   = state_new.RcgRwR;
        
    end
    
    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),FI_acc(1,:),'y',t(1:2:(end-1)),FI_vel(1,:),'r')
    subplot(3,1,2); plot(t(1:2:(end-1)),FI_acc(2,:),'y',t(1:2:(end-1)),FI_vel(2,:),'r')
    subplot(3,1,3); plot(t(1:2:(end-1)),FI_acc(3,:),'y',t(1:2:(end-1)),FI_vel(3,:),'r')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),MI_acc(1,:),'y',t(1:2:(end-1)),MI_vel(1,:),'r')
    subplot(3,1,2); plot(t(1:2:(end-1)),MI_acc(2,:),'y',t(1:2:(end-1)),MI_vel(2,:),'r')
    subplot(3,1,3); plot(t(1:2:(end-1)),MI_acc(3,:),'y',t(1:2:(end-1)),MI_vel(3,:),'r')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),vb(1,:))
    subplot(3,1,2); plot(t(1:2:(end-1)),vb(2,:))
    subplot(3,1,3); plot(t(1:2:(end-1)),vb(3,:))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),wb(1,:))
    subplot(3,1,2); plot(t(1:2:(end-1)),wb(2,:))
    subplot(3,1,3); plot(t(1:2:(end-1)),wb(3,:))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),ab(1,:))
    subplot(3,1,2); plot(t(1:2:(end-1)),ab(2,:))
    subplot(3,1,3); plot(t(1:2:(end-1)),ab(3,:))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),w_dot_b(1,:))
    subplot(3,1,2); plot(t(1:2:(end-1)),w_dot_b(2,:))
    subplot(3,1,3); plot(t(1:2:(end-1)),w_dot_b(3,:))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),LM(1,:))
    subplot(3,1,2); plot(t(1:2:(end-1)),LM(2,:))
    subplot(3,1,3); plot(t(1:2:(end-1)),LM(3,:))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),AM(1,:))
    subplot(3,1,2); plot(t(1:2:(end-1)),AM(2,:))
    subplot(3,1,3); plot(t(1:2:(end-1)),AM(3,:))
    hold off 
    
    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),vb_0(1,:),'b',t(1:2:(end-1)),ones(nr_points*nr_wb)*mean(vb_0(1,:)),'r')
    subplot(3,1,2); plot(t(1:2:(end-1)),vb_0(2,:),'b',t(1:2:(end-1)),ones(nr_points*nr_wb)*mean(vb_0(2,:)),'r')
    subplot(3,1,3); plot(t(1:2:(end-1)),vb_0(3,:),'b',t(1:2:(end-1)),ones(nr_points*nr_wb)*mean(vb_0(3,:)),'r')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),wb_0(1,:),'b',t(1:2:(end-1)),ones(nr_points*nr_wb)*mean(wb_0(1,:)),'r')
    subplot(3,1,2); plot(t(1:2:(end-1)),wb_0(2,:),'b',t(1:2:(end-1)),ones(nr_points*nr_wb)*mean(wb_0(2,:)),'r')
    subplot(3,1,3); plot(t(1:2:(end-1)),wb_0(3,:),'b',t(1:2:(end-1)),ones(nr_points*nr_wb)*mean(wb_0(3,:)),'r')
    hold off 
    
    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),RcgLwL(1,:),'b',t(1:2:(end-1)),ones(nr_points*nr_wb)*mean(RcgLwL(1,:)),'r')
    subplot(3,1,2); plot(t(1:2:(end-1)),RcgLwL(2,:),'b',t(1:2:(end-1)),ones(nr_points*nr_wb)*mean(RcgLwL(2,:)),'r')
    subplot(3,1,3); plot(t(1:2:(end-1)),RcgLwL(3,:),'b',t(1:2:(end-1)),ones(nr_points*nr_wb)*mean(RcgLwL(3,:)),'r')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),RcgRwR(1,:),'b',t(1:2:(end-1)),ones(nr_points*nr_wb)*mean(RcgRwR(1,:)),'r')
    subplot(3,1,2); plot(t(1:2:(end-1)),RcgRwR(2,:),'b',t(1:2:(end-1)),ones(nr_points*nr_wb)*mean(RcgRwR(2,:)),'r')
    subplot(3,1,3); plot(t(1:2:(end-1)),RcgRwR(3,:),'b',t(1:2:(end-1)),ones(nr_points*nr_wb)*mean(RcgRwR(3,:)),'r')
    hold off 

end

