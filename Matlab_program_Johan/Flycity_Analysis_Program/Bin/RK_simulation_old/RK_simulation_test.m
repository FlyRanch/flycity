function RK_simulation_test( settings, pathDB )

    % Run a Runge Kutta simulation for the length of a wingbeat and
    % simulate the aerodynamic, inertial and gravtitational forces on the
    % body.
    
    seq_nr      = 5;
    wb_id       = 35;
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
    
    % Setup wing kinematics during the wingbeat:
    
%     a_theta_L   = pathDB.poly_fit.a_fit.theta_L(:,wb_id,seq_nr);
%     a_eta_L     = pathDB.poly_fit.a_fit.eta_L(:,wb_id,seq_nr);
%     a_phi_L     = pathDB.poly_fit.a_fit.phi_L(:,wb_id,seq_nr);  
%     
%     a_theta_R   = pathDB.poly_fit.a_fit.theta_R(:,wb_id,seq_nr);
%     a_eta_R     = pathDB.poly_fit.a_fit.eta_R(:,wb_id,seq_nr);
%     a_phi_R     = pathDB.poly_fit.a_fit.phi_R(:,wb_id,seq_nr);  
%    
%     down_up     = pathDB.poly_fit.a_fit.down_up(wb_id,seq_nr);
%     
%     f           = pathDB.poly_fit.a_fit.f(wb_id,seq_nr);
    
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
    
    wb_loc      = [ 1 nr_points];
    down_loc    = [ 1 round(down_up*nr_points)];
    up_loc      = [ round(down_up*nr_points)+1 nr_points];
    
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
    
%     % Compute the baseline forces and moments for the average wingbeat:
%     
%     kine.u_strk                 = zeros(3,nr_points*2+1);
%     kine.w_strk                 = zeros(3,nr_points*2+1);
%     kine.wL                     = wing_kin.wL;
%     kine.wR                     = wing_kin.wR;
%     kine.RL                     = wing_kin.RL;
%     kine.RR                     = wing_kin.RR;
%     kine.R_strk                 = R_strk;
%     
%     wb.wb_loc     = [1 nr_points*2+1];
%     wb.down_loc   = [1 round((nr_points*2+1)*down_up)];
%     wb.up_loc     = [round((nr_points*2+1)*down_up)+1 nr_points*2+1];
%     wb.dt         = dt/2;
%     
%     [ FM_strkpln_avg, ~, ~ , ~, ~, ~, ~, ~, ~] = Aerodynamic_forces( kine, body_model, wing_model, wb, 1);
%     
%     FA_avg_strkpln = mean(FM_strkpln_avg(1:3,:),2);
%     MA_avg_strkpln = mean(FM_strkpln_avg(4:6,:),2);
%     
%     FA_avg_b        = R_strk'*FA_avg_strkpln;
%     MA_avg_b        = R_strk'*MA_avg_strkpln;
%        
%     % Compute location center of gravity:
%     
%     r_cg_z = body_model.Joint_left(3);
%     r_cg_x = (-MA_avg_b(2)-r_cg_z*FA_avg_b(1))/FA_avg_b(3);
%     
%     body_model.cg_b = body_model.cg_b+[ body_model.Joint_left(1)+r_cg_x; 0; 0];
%     
%     body_model.Inertia = body_model.Inertia+[0 0 0; 0 body_model.mass_fly*body_model.cg_b(2)^2 0; 0 0 body_model.mass_fly*body_model.cg_b(2)^2];
%     
%     body_model.Inertia
%     
%     MA_avg_strkpln

    figure()
    hold on
    subplot(3,1,1); plot(t,kine.theta_L)
    subplot(3,1,2); plot(t,kine.eta_L)
    subplot(3,1,3); plot(t,kine.phi_L)
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t,kine.theta_R)
    subplot(3,1,2); plot(t,kine.eta_R)
    subplot(3,1,3); plot(t,kine.phi_R)
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t,wing_kin.wL(1,:))
    subplot(3,1,2); plot(t,wing_kin.wL(2,:))
    subplot(3,1,3); plot(t,wing_kin.wL(3,:))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t,wing_kin.wR(1,:))
    subplot(3,1,2); plot(t,wing_kin.wR(2,:))
    subplot(3,1,3); plot(t,wing_kin.wR(3,:))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t,wing_kin.w_dot_L(1,:))
    subplot(3,1,2); plot(t,wing_kin.w_dot_L(2,:))
    subplot(3,1,3); plot(t,wing_kin.w_dot_L(3,:))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t,wing_kin.w_dot_R(1,:))
    subplot(3,1,2); plot(t,wing_kin.w_dot_R(2,:))
    subplot(3,1,3); plot(t,wing_kin.w_dot_R(3,:))
    hold off
    
    
    
    % Retrieve initial conditions:
    
%     wb_loc_pathDB       = pathDB.wingbeats.wingbeat_loc(wb_id,1,seq_nr);
%     
%     state.dt            = dt;
%     state.id            = 1;
%     state.xyz           = pathDB.filt.xyz(wb_loc_pathDB,:,seq_nr)';
%     state.qb            = pathDB.filt.qB(wb_loc_pathDB,:,seq_nr)';
%     state.Rb            = pathDB.rot_mat.RB(:,:,wb_loc_pathDB,seq_nr);
%     state.vb            = R_strk*pathDB.strkpln_kin.uvw(wb_loc_pathDB,:,seq_nr)';
%     state.wb            = pathDB.filt.wB(wb_loc_pathDB,:,seq_nr)';
%     state.ab            = R_strk*pathDB.strkpln_kin.a_xyz(wb_loc_pathDB,:,seq_nr)';
%     state.w_dot_b       = zeros(3,1);
%     state.alfa_L_old    = nan(nr_sect,1);
%     state.alfa_R_old    = nan(nr_sect,1);
    
    state.dt            = dt;
    state.id            = 1;
    state.xyz           = zeros(3,1);
    state.qb            = quat2mat(R_strk')';
    state.Rb            = R_strk';
    state.vb            = zeros(3,1);
    state.wb            = zeros(3,1);
    state.ab            = zeros(3,1);
    state.w_dot_b       = zeros(3,1);
    state.alfa_L_old    = nan(nr_sect,1);
    state.alfa_R_old    = nan(nr_sect,1);
    
    state.dt
    state.id
    state.xyz
    state.qb
    state.Rb
    state.vb
    state.wb
    state.ab
    state.w_dot_b
    state.alfa_L_old
    state.alfa_R_old
    
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
    
        
        
    for i = 1:nr_points
    
        [ state_new ] = RK_updater( state, wing_kin, body_model, wing_model, wingbeat );
        
        xyz_sim(:,i)            = state_new.xyz;
        qb_sim(:,i)             = state_new.qb;
        Rb_sim(:,:,i)           = state_new.Rb;
        vb_sim(:,i)             = state_new.vb;
        wb_sim(:,i)             = state_new.wb;
        ab_sim(:,i)             = state_new.ab;
        w_dot_b_sim(:,i)        = state_new.w_dot_b;
        alfa_L_sim(:,i)         = state_new.alfa_L_old;
        alfa_R_sim(:,i)         = state_new.alfa_R_old;
        alfa_dot_L_sim(:,i)     = state_new.alfa_dot_L;
        alfa_dot_R_sim(:,i)     = state_new.alfa_dot_R;
        FMA_b_sim(:,i)          = state_new.FMA_b;
        FMA_strkpln_sim(:,i)    = state_new.FMA_strkpln;
        FI_acc_b_sim(:,i)       = state_new.FMI_b.F_I_acc;
        FI_vel_b_sim(:,i)       = state_new.FMI_b.F_I_vel;
        FI_acc_strkpln_sim(:,i) = state_new.FMI_strkpln.F_I_acc;
        FI_vel_strkpln_sim(:,i) = state_new.FMI_strkpln.F_I_vel;
        MI_acc_b_sim(:,i)       = state_new.FMI_b.M_I_acc;
        MI_vel_b_sim(:,i)       = state_new.FMI_b.M_I_vel;
        MI_acc_strkpln_sim(:,i) = state_new.FMI_strkpln.M_I_acc;
        MI_vel_strkpln_sim(:,i) = state_new.FMI_strkpln.M_I_vel;
        Fg_b_sim(:,i)           = state_new.Fg_b;
        Fg_strkpln_sim(:,i)     = state_new.Fg_strkpln;
        
        state.id            = state_new.id;
        state.xyz           = state_new.xyz;
        state.qb            = state_new.qb;
        state.Rb            = state_new.Rb;
        state.vb            = state_new.vb;
        state.wb            = state_new.wb;
        state.ab            = state_new.ab;
        state.w_dot_b       = state_new.w_dot_b;
        state.alfa_L_old    = state_new.alfa_L_old;
        state.alfa_R_old    = state_new.alfa_R_old;
        
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
    subplot(3,1,1); plot(t(1:2:(end-1)),vb_sim(1,:))
    title('Velocity body')
    subplot(3,1,2); plot(t(1:2:(end-1)),vb_sim(2,:))
    subplot(3,1,3); plot(t(1:2:(end-1)),vb_sim(3,:))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),wb_sim(1,:))
    title('Angular velocity body')
    subplot(3,1,2); plot(t(1:2:(end-1)),wb_sim(2,:))
    subplot(3,1,3); plot(t(1:2:(end-1)),wb_sim(3,:))
    hold off
        
    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),ab_sim(1,:))
    title('Acceleration body')
    subplot(3,1,2); plot(t(1:2:(end-1)),ab_sim(2,:))
    subplot(3,1,3); plot(t(1:2:(end-1)),ab_sim(3,:))
    hold off
        
    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),w_dot_b_sim(1,:))
    title('Angular acceleration body')
    subplot(3,1,2); plot(t(1:2:(end-1)),w_dot_b_sim(2,:))
    subplot(3,1,3); plot(t(1:2:(end-1)),w_dot_b_sim(3,:))
    hold off
        
    figure()
    hold on
    for i = 1:nr_sect
    plot(t(1:2:(end-1)),alfa_L_sim(i,:))
    end
    title('alfa L')
    hold off
    
    figure()
    hold on
    for i = 1:nr_sect
    plot(t(1:2:(end-1)),alfa_R_sim(i,:))
    end
    title('alfa R')
    hold off
    
    figure()
    hold on
    for i = 1:nr_sect
    plot(t(1:2:(end-1)),alfa_dot_L_sim(i,:))
    end
    title('alfa dot L')
    hold off
    
    figure()
    hold on
    for i = 1:nr_sect
    plot(t(1:2:(end-1)),alfa_dot_R_sim(i,:))
    end
    title('alfa dot R')
    hold off
    
    figure()
    hold on
    subplot(2,1,1); plot(t(1:2:(end-1)),FMA_b_sim(1,:),'r',t(1:2:(end-1)),FMA_b_sim(2,:),'g',t(1:2:(end-1)),FMA_b_sim(3,:),'b')
    title('Aerodynamic forces body frame')
    subplot(2,1,2); plot(t(1:2:(end-1)),FMA_b_sim(4,:),'r',t(1:2:(end-1)),FMA_b_sim(5,:),'g',t(1:2:(end-1)),FMA_b_sim(6,:),'b')
    title('Aerodynamic moments body frame')
    hold off
    
    figure()
    hold on
    subplot(2,1,1); plot(t(1:2:(end-1)),FMA_strkpln_sim(1,:),'r',t(1:2:(end-1)),FMA_strkpln_sim(2,:),'g',t(1:2:(end-1)),FMA_strkpln_sim(3,:),'b')
    title('Aerodynamic forces strokeplane')
    subplot(2,1,2); plot(t(1:2:(end-1)),FMA_strkpln_sim(4,:),'r',t(1:2:(end-1)),FMA_strkpln_sim(5,:),'g',t(1:2:(end-1)),FMA_strkpln_sim(6,:),'b')
    title('Aerodynamic moments strokeplane')
    hold off
    
    figure()
    hold on
    subplot(2,1,1); plot(t(1:2:(end-1)),FI_acc_b_sim(1,:),'r',t(1:2:(end-1)),FI_acc_b_sim(2,:),'g',t(1:2:(end-1)),FI_acc_b_sim(3,:),'b')
    title('Inertial forces due to accelerations body frame')
    subplot(2,1,2); plot(t(1:2:(end-1)),FI_vel_b_sim(1,:),'r',t(1:2:(end-1)),FI_vel_b_sim(2,:),'g',t(1:2:(end-1)),FI_vel_b_sim(3,:),'b')
    title('Inertial forces due to velocities body frame')
    hold off
    
    figure()
    hold on
    subplot(2,1,1); plot(t(1:2:(end-1)),FI_acc_strkpln_sim(1,:),'r',t(1:2:(end-1)),FI_acc_strkpln_sim(2,:),'g',t(1:2:(end-1)),FI_acc_strkpln_sim(3,:),'b')
    title('Inertial forces due to accelerations strokeplane frame')
    subplot(2,1,2); plot(t(1:2:(end-1)),FI_vel_strkpln_sim(1,:),'r',t(1:2:(end-1)),FI_vel_strkpln_sim(2,:),'g',t(1:2:(end-1)),FI_vel_strkpln_sim(3,:),'b')
    title('Inertial forces due to velocities strokeplane frame')
    hold off
    
    figure()
    hold on
    subplot(2,1,1); plot(t(1:2:(end-1)),MI_acc_b_sim(1,:),'r',t(1:2:(end-1)),MI_acc_b_sim(2,:),'g',t(1:2:(end-1)),MI_acc_b_sim(3,:),'b')
    title('Inertial moments due to accelerations body frame')
    subplot(2,1,2); plot(t(1:2:(end-1)),MI_vel_b_sim(1,:),'r',t(1:2:(end-1)),MI_vel_b_sim(2,:),'g',t(1:2:(end-1)),MI_vel_b_sim(3,:),'b')
    title('Inertial moments due to velocities body frame')
    hold off
    
    figure()
    hold on
    subplot(2,1,1); plot(t(1:2:(end-1)),MI_acc_strkpln_sim(1,:),'r',t(1:2:(end-1)),MI_acc_strkpln_sim(2,:),'g',t(1:2:(end-1)),MI_acc_strkpln_sim(3,:),'b')
    title('Inertial moments due to accelerations strokeplane frame')
    subplot(2,1,2); plot(t(1:2:(end-1)),MI_vel_strkpln_sim(1,:),'r',t(1:2:(end-1)),MI_vel_strkpln_sim(2,:),'g',t(1:2:(end-1)),MI_vel_strkpln_sim(3,:),'b')
    title('Inertial moments due to velocities strokeplane frame')
    hold off
    
    figure()
    plot(t(1:2:(end-1)),Fg_b_sim(1,:),'r',t(1:2:(end-1)),Fg_b_sim(2,:),'g',t(1:2:(end-1)),Fg_b_sim(3,:),'b')
    title('Gravity force body frame')
    
    figure()
    plot(t(1:2:(end-1)),Fg_strkpln_sim(1,:),'r',t(1:2:(end-1)),Fg_strkpln_sim(2,:),'g',t(1:2:(end-1)),Fg_strkpln_sim(3,:),'b')
    title('Gravity force strokeplane frame')
    
    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),FMA_strkpln_sim(1,:),'b',t(1:2:(end-1)),-FI_acc_strkpln_sim(1,:),'y',t(1:2:(end-1)),FI_vel_strkpln_sim(1,:),'r',t(1:2:(end-1)),Fg_strkpln_sim(1,:),'g',t(1:2:(end-1)),FMA_strkpln_sim(1,:)-FI_acc_strkpln_sim(1,:)+FI_vel_strkpln_sim(1,:)+Fg_strkpln_sim(1,:),'k')
    subplot(3,1,2); plot(t(1:2:(end-1)),FMA_strkpln_sim(2,:),'b',t(1:2:(end-1)),-FI_acc_strkpln_sim(2,:),'y',t(1:2:(end-1)),FI_vel_strkpln_sim(2,:),'r',t(1:2:(end-1)),Fg_strkpln_sim(2,:),'g',t(1:2:(end-1)),FMA_strkpln_sim(2,:)-FI_acc_strkpln_sim(2,:)+FI_vel_strkpln_sim(2,:)+Fg_strkpln_sim(2,:),'k')
    subplot(3,1,3); plot(t(1:2:(end-1)),FMA_strkpln_sim(3,:),'b',t(1:2:(end-1)),-FI_acc_strkpln_sim(3,:),'y',t(1:2:(end-1)),FI_vel_strkpln_sim(3,:),'r',t(1:2:(end-1)),Fg_strkpln_sim(3,:),'g',t(1:2:(end-1)),FMA_strkpln_sim(3,:)-FI_acc_strkpln_sim(3,:)+FI_vel_strkpln_sim(3,:)+Fg_strkpln_sim(3,:),'k')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t(1:2:(end-1)),FMA_strkpln_sim(4,:),'b',t(1:2:(end-1)),-MI_acc_strkpln_sim(1,:),'y',t(1:2:(end-1)),MI_vel_strkpln_sim(1,:),'r',t(1:2:(end-1)),FMA_strkpln_sim(4,:)-MI_acc_strkpln_sim(1,:)+MI_vel_strkpln_sim(1,:),'k')
    subplot(3,1,2); plot(t(1:2:(end-1)),FMA_strkpln_sim(5,:),'b',t(1:2:(end-1)),-MI_acc_strkpln_sim(2,:),'y',t(1:2:(end-1)),MI_vel_strkpln_sim(2,:),'r',t(1:2:(end-1)),FMA_strkpln_sim(5,:)-MI_acc_strkpln_sim(2,:)+MI_vel_strkpln_sim(2,:),'k')
    subplot(3,1,3); plot(t(1:2:(end-1)),FMA_strkpln_sim(6,:),'b',t(1:2:(end-1)),-MI_acc_strkpln_sim(3,:),'y',t(1:2:(end-1)),MI_vel_strkpln_sim(3,:),'r',t(1:2:(end-1)),FMA_strkpln_sim(6,:)-MI_acc_strkpln_sim(3,:)+MI_vel_strkpln_sim(3,:),'k')
    hold off
    
    mean(FMA_strkpln_sim(5,:))
    mean(MI_vel_strkpln_sim(2,:))
    mean(FMA_strkpln_sim(5,:)+MI_vel_strkpln_sim(2,:))
    
    figure()
    plot(t(1:2:(end-1)),FMA_strkpln_sim(5,:),'b',t(1:2:(end-1)),MI_vel_strkpln_sim(2,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1)))).*mean(FMA_strkpln_sim(5,:)+MI_vel_strkpln_sim(2,:)),'k')
    
    figure()
    plot(t(1:2:(end-1)),MI_vel_strkpln_sim(2,:),'r',t(1:2:(end-1)),ones(1,length(t(1:2:(end-1)))).*mean(MI_vel_strkpln_sim(2,:)),'k')
    
    x_circ = cos((2*pi/(99))*(1:100));
    y_circ = sin((2*pi/(99))*(1:100));
    
    nr_plot = length(1:10:nr_points);
    
    circ_strkpln = nan(3,100,nr_plot);
    
    for i = 1:10:nr_points
        circ_strkpln(:,:,i) = R_strk*Rb_sim(:,:,i)*[x_circ; y_circ; zeros(1,100)];
    end
    
    j = 1;
    
    figure()
    hold on
    for i = 1:10:nr_points
        plot3(xyz_sim(1,i)+circ_strkpln(1,:,i),xyz_sim(2,i)+circ_strkpln(2,:,i),xyz_sim(3,i)+circ_strkpln(3,:,i),'Color',[0.5*((nr_plot+j)/(nr_plot+1)) 0.5*((nr_plot-j)/(nr_plot+1)) 0.5*((nr_plot-j)/(nr_plot+1))])  
        j = j+1;
    end
    axis equal
    hold off
    
    
end

