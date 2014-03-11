function Conservation_test_old( settings, pathDB )


    % Test for conservation of linear and angular momentum:
    
    seq_nr      = 5;
    nr_points   = 200;
    
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
    
    a_theta_L = pathDB.poly_fit.a_avg.theta_LR(:,seq_nr);
    a_eta_L   = pathDB.poly_fit.a_avg.eta_LR(:,seq_nr);
    a_phi_L   = pathDB.poly_fit.a_avg.phi_LR(:,seq_nr);  
    
    a_theta_R = pathDB.poly_fit.a_avg.theta_LR(:,seq_nr);
    a_eta_R   = pathDB.poly_fit.a_avg.eta_LR(:,seq_nr);
    a_phi_R   = pathDB.poly_fit.a_avg.phi_LR(:,seq_nr);  
   
    down_up = pathDB.poly_fit.a_avg.down_up(seq_nr);
    
    f       = pathDB.poly_fit.a_avg.f(seq_nr);
    
    dt = (1/f)/(nr_points-1);
    
    a_fit_L.a_theta     = a_theta_L;
    a_fit_L.a_eta       = a_eta_L;
    a_fit_L.a_phi       = a_phi_L;
    a_fit_L.f           = f;
    a_fit_L.down_up     = down_up;
    
    a_fit_R.a_theta     = a_theta_R;
    a_fit_R.a_eta       = a_eta_R;
    a_fit_R.a_phi       = a_phi_R;
    a_fit_R.f           = f;
    a_fit_R.down_up     = down_up;
    
    nr_wb_sim       = 10;
    
    t               = nan(1,(nr_points*2)*nr_wb_sim+1);
    theta_L         = nan(1,(nr_points*2)*nr_wb_sim+1);
    eta_L           = nan(1,(nr_points*2)*nr_wb_sim+1);
    phi_L           = nan(1,(nr_points*2)*nr_wb_sim+1);
    theta_dot_L     = nan(1,(nr_points*2)*nr_wb_sim+1);
    eta_dot_L       = nan(1,(nr_points*2)*nr_wb_sim+1);
    phi_dot_L       = nan(1,(nr_points*2)*nr_wb_sim+1);
    theta_ddot_L    = nan(1,(nr_points*2)*nr_wb_sim+1);
    eta_ddot_L      = nan(1,(nr_points*2)*nr_wb_sim+1);
    phi_ddot_L      = nan(1,(nr_points*2)*nr_wb_sim+1);
    theta_R         = nan(1,(nr_points*2)*nr_wb_sim+1);
    eta_R           = nan(1,(nr_points*2)*nr_wb_sim+1);
    phi_R           = nan(1,(nr_points*2)*nr_wb_sim+1);
    theta_dot_R     = nan(1,(nr_points*2)*nr_wb_sim+1);
    eta_dot_R       = nan(1,(nr_points*2)*nr_wb_sim+1);
    phi_dot_R       = nan(1,(nr_points*2)*nr_wb_sim+1);
    theta_ddot_R    = nan(1,(nr_points*2)*nr_wb_sim+1);
    eta_ddot_R      = nan(1,(nr_points*2)*nr_wb_sim+1);
    phi_ddot_R      = nan(1,(nr_points*2)*nr_wb_sim+1);
    
    wb_loc      = nan(nr_wb_sim,2);
    down_loc    = nan(nr_wb_sim,2);
    up_loc      = nan(nr_wb_sim,2);
    
    for k = 1:nr_wb_sim
    
    [ kine_L ] = poly_derivatives( a_fit_L, nr_points*2+1 );
    
    [ kine_R ] = poly_derivatives( a_fit_R, nr_points*2+1 );
    
        t               = kine_L.t;    
        theta_L         = kine_L.theta;
        eta_L           = kine_L.eta;
        phi_L           = kine_L.phi;
        theta_dot_L     = kine_L.theta_dot;
        eta_dot_L       = -kine_L.eta_dot;
        phi_dot_L       = -kine_L.phi_dot;
        theta_ddot_L    = kine_L.theta_ddot;
        eta_ddot_L      = -kine_L.eta_ddot;
        phi_ddot_L      = -kine_L.phi_ddot;
        theta_R         = kine_R.theta;
        eta_R           = kine_R.eta;
        phi_R           = kine_R.phi;
        theta_dot_R     = -kine_R.theta_dot;
        eta_dot_R       = -kine_R.eta_dot;
        phi_dot_R       = kine_R.phi_dot;
        theta_ddot_R    = -kine_R.theta_ddot;
        eta_ddot_R      = -kine_R.eta_ddot;
        phi_ddot_R      = kine_R.phi_ddot;

        wb_loc          = [ 1 nr_points];
        down_loc        = [ 1 round(down_up*nr_points)];
        up_loc          = [ round(down_up*nr_points)+1 nr_points];
    
    end
    
    wingbeat.wb_loc     = wb_loc;
    wingbeat.down_loc   = down_loc;
    wingbeat.up_loc     = up_loc;
    
    R_strk = pathDB.rot_mat.Rstr;
    
    wingbeat.R_strk = R_strk;
    
    kine.R_strk  = R_strk;
    kine.theta_L = theta_L;
    kine.eta_L   = eta_L;
    kine.phi_L   = phi_L;
    kine.theta_R = theta_R;
    kine.eta_R   = eta_R;
    kine.phi_R   = phi_R;
    kine.theta_dot_L = theta_dot_L;
    kine.eta_dot_L   = eta_dot_L;
    kine.phi_dot_L   = phi_dot_L;
    kine.theta_dot_R = theta_dot_R;
    kine.eta_dot_R   = eta_dot_R;
    kine.phi_dot_R   = phi_dot_R;
    kine.theta_ddot_L = theta_ddot_L;
    kine.eta_ddot_L   = eta_ddot_L;
    kine.phi_ddot_L   = phi_ddot_L;
    kine.theta_ddot_R = theta_ddot_R;
    kine.eta_ddot_R   = eta_ddot_R;
    kine.phi_ddot_R   = phi_ddot_R;
    
    [ rot_mat ] = wingkin_rotation_mat( kine, nr_points*2+1 );
    
    wing_kin.RL          = rot_mat.RL;
    wing_kin.RR          = rot_mat.RR;
    wing_kin.wL          = rot_mat.wL;
    wing_kin.wR          = rot_mat.wR;
    wing_kin.w_dot_L     = rot_mat.w_dot_L;
    wing_kin.w_dot_R     = rot_mat.w_dot_R;

    
    
    
    % Retrieve initial conditions:

    
    state.dt            = dt;
    state.id            = 1;
    state.xyz           = zeros(3,1);
    state.qb            = quat2mat(R_strk')';
    state.Rb            = R_strk';
    state.vb            = zeros(3,1);
    state.wb            = zeros(3,1);
    state.ab            = zeros(3,1);
    state.w_dot_b       = zeros(3,1);
    
    state.dt
    state.id
    state.xyz
    state.qb
    state.Rb
    state.vb
    state.wb
    state.ab
    state.w_dot_b
    
    % Run the simulation:
    
    xyz_sim             = nan(3,nr_points);
    qb_sim              = nan(4,nr_points);
    Rb_sim              = nan(3,3,nr_points);
    vb_sim              = nan(3,nr_points);
    wb_sim              = nan(3,nr_points);
    ab_sim              = nan(3,nr_points);
    w_dot_b_sim         = nan(3,nr_points);
    alfa_L_sim          = nan(nr_sect,nr_points);
    alfa_R_sim          = nan(nr_sect,nr_points);
    alfa_dot_L_sim      = nan(nr_sect,nr_points);
    alfa_dot_R_sim      = nan(nr_sect,nr_points);
    FMA_b_sim           = nan(6,nr_points);
    FMA_strkpln_sim     = nan(6,nr_points);
    FI_acc_b_sim        = nan(3,nr_points);
    FI_vel_b_sim        = nan(3,nr_points);
    FI_acc_strkpln_sim  = nan(3,nr_points);
    FI_vel_strkpln_sim  = nan(3,nr_points);
    MI_acc_b_sim        = nan(3,nr_points);
    MI_vel_b_sim        = nan(3,nr_points);
    MI_acc_strkpln_sim  = nan(3,nr_points);
    MI_vel_strkpln_sim  = nan(3,nr_points);
    Fg_b_sim            = nan(3,nr_points);
    Fg_strkpln_sim      = nan(3,nr_points);
    Lin_momentum        = nan(3,nr_points);
    Ang_momentum        = nan(3,nr_points);
    KE_sim              = nan(nr_points,1);
    KE_lin_sim          = nan(nr_points,1);
    KE_ang_sim          = nan(nr_points,1);
        
    
    
%     [ vb_0_t, wb_0_t ] = vb_0_and_wb_0( body_model, wing_model, wing_kin );
%     
%     vb_0 = mean(vb_0_t,2);
%     
%     wb_0 = mean(wb_0_t,2);
%     
%     state.vb            = vb_0_t(:,1);
%     state.wb            = wb_0_t(:,1);
    
    for i = 1:nr_points
        
        [ state_new ] = RK_updater_inertia_only( state, wing_kin, body_model, wing_model, wingbeat );
        
        xyz_sim(:,i)            = state_new.xyz;
        qb_sim(:,i)             = state_new.qb;
        Rb_sim(:,:,i)           = state_new.Rb;
        vb_sim(:,i)             = state_new.vb;
        wb_sim(:,i)             = state_new.wb;
        ab_sim(:,i)             = state_new.ab;
        w_dot_b_sim(:,i)        = state_new.w_dot_b;
        FI_acc_b_sim(:,i)       = state_new.FMI_b.F_I_acc;
        FI_vel_b_sim(:,i)       = state_new.FMI_b.F_I_vel;
        FI_acc_strkpln_sim(:,i) = state_new.FMI_strkpln.F_I_acc;
        FI_vel_strkpln_sim(:,i) = state_new.FMI_strkpln.F_I_vel;
        MI_acc_b_sim(:,i)       = state_new.FMI_b.M_I_acc;
        MI_vel_b_sim(:,i)       = state_new.FMI_b.M_I_vel;
        MI_acc_strkpln_sim(:,i) = state_new.FMI_strkpln.M_I_acc;
        MI_vel_strkpln_sim(:,i) = state_new.FMI_strkpln.M_I_vel;
        Lin_momentum(:,i)       = state_new.LM;
        Ang_momentum(:,i)       = state_new.AM;
        KE_sim(i)               = state_new.KE;
        KE_lin_sim(i)           = state_new.KE_lin;
        KE_ang_sim(i)           = state_new.KE_ang;
        
        state.id            = state_new.id;
        state.xyz           = state_new.xyz;
        state.qb            = state_new.qb;
        state.Rb            = state_new.Rb;
        state.vb            = state_new.vb;
        state.wb            = state_new.wb;
        state.ab            = state_new.ab;
        state.w_dot_b       = state_new.w_dot_b;
        
        i
    
    end

    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),xyz_sim(1,:))
    title('position cg body')
    subplot(3,1,2); plot(t(1:2:(end-1)),xyz_sim(2,:))
    subplot(3,1,3); plot(t(1:2:(end-1)),xyz_sim(3,:))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),vb_sim(1,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(vb_sim(1,:)),'k')
    title('Velocity body')
    subplot(3,1,2); plot(t(1:2:(end-1)),vb_sim(2,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(vb_sim(2,:)),'k')
    subplot(3,1,3); plot(t(1:2:(end-1)),vb_sim(3,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(vb_sim(3,:)),'k')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),wb_sim(1,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(wb_sim(1,:)),'k')
    title('Angular velocity body')
    subplot(3,1,2); plot(t(1:2:(end-1)),wb_sim(2,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(wb_sim(2,:)),'k')
    subplot(3,1,3); plot(t(1:2:(end-1)),wb_sim(3,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(wb_sim(3,:)),'k')
    hold off
        
    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),ab_sim(1,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(ab_sim(1,:)),'k')
    title('Acceleration body')
    subplot(3,1,2); plot(t(1:2:(end-1)),ab_sim(2,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(ab_sim(2,:)),'k')
    subplot(3,1,3); plot(t(1:2:(end-1)),ab_sim(3,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(ab_sim(3,:)),'k')
    hold off
        
    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),w_dot_b_sim(1,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(w_dot_b_sim(1,:)),'k')
    title('Angular acceleration body')
    subplot(3,1,2); plot(t(1:2:(end-1)),w_dot_b_sim(2,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(w_dot_b_sim(2,:)),'k')
    subplot(3,1,3); plot(t(1:2:(end-1)),w_dot_b_sim(3,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(w_dot_b_sim(3,:)),'k')
    hold off

    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),FI_vel_strkpln_sim(1,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(FI_vel_strkpln_sim(1,:)),'k')
    title('FI')
    subplot(3,1,2); plot(t(1:2:(end-1)),FI_vel_strkpln_sim(2,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(FI_vel_strkpln_sim(2,:)),'k')
    subplot(3,1,3); plot(t(1:2:(end-1)),FI_vel_strkpln_sim(3,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(FI_vel_strkpln_sim(3,:)),'k')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),MI_vel_strkpln_sim(1,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(MI_vel_strkpln_sim(1,:)),'k')
    title('MI')
    subplot(3,1,2); plot(t(1:2:(end-1)),MI_vel_strkpln_sim(2,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(MI_vel_strkpln_sim(2,:)),'k')
    subplot(3,1,3); plot(t(1:2:(end-1)),MI_vel_strkpln_sim(3,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(MI_vel_strkpln_sim(3,:)),'k')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),Lin_momentum(1,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(Lin_momentum(1,:)),'k')
    title('Linear Momentum')
    subplot(3,1,2); plot(t(1:2:(end-1)),Lin_momentum(2,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(Lin_momentum(2,:)),'k')
    subplot(3,1,3); plot(t(1:2:(end-1)),Lin_momentum(3,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(Lin_momentum(3,:)),'k')
    hold off

    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),Ang_momentum(1,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(Ang_momentum(1,:)),'k')
    title('Angular Momentum')
    subplot(3,1,2); plot(t(1:2:(end-1)),Ang_momentum(2,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(Ang_momentum(2,:)),'k')
    subplot(3,1,3); plot(t(1:2:(end-1)),Ang_momentum(3,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(Ang_momentum(3,:)),'k')
    hold off
    
    figure()
    plot(t(1:2:(end-1)),KE_sim,t(1:2:(end-1)),KE_lin_sim,t(1:2:(end-1)),KE_ang_sim)
    
    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),-FI_acc_strkpln_sim(1,:),'y',t(1:2:(end-1)),FI_vel_strkpln_sim(1,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(-FI_acc_strkpln_sim(1,:)+FI_vel_strkpln_sim(1,:),2),'k')
    subplot(3,1,2); plot(t(1:2:(end-1)),-FI_acc_strkpln_sim(2,:),'y',t(1:2:(end-1)),FI_vel_strkpln_sim(2,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(-FI_acc_strkpln_sim(2,:)+FI_vel_strkpln_sim(2,:),2),'k')
    subplot(3,1,3); plot(t(1:2:(end-1)),-FI_acc_strkpln_sim(3,:),'y',t(1:2:(end-1)),FI_vel_strkpln_sim(3,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(-FI_acc_strkpln_sim(3,:)+FI_vel_strkpln_sim(3,:),2),'k')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),-MI_acc_strkpln_sim(1,:),'y',t(1:2:(end-1)),MI_vel_strkpln_sim(1,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(-MI_acc_strkpln_sim(1,:)+MI_vel_strkpln_sim(1,:),2),'k')
    subplot(3,1,2); plot(t(1:2:(end-1)),-MI_acc_strkpln_sim(2,:),'y',t(1:2:(end-1)),MI_vel_strkpln_sim(2,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(-MI_acc_strkpln_sim(2,:)+MI_vel_strkpln_sim(2,:),2),'k')
    subplot(3,1,3); plot(t(1:2:(end-1)),-MI_acc_strkpln_sim(3,:),'y',t(1:2:(end-1)),MI_vel_strkpln_sim(3,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1))))*mean(-MI_acc_strkpln_sim(3,:)+MI_vel_strkpln_sim(3,:),2),'k')
    hold off
    
end

