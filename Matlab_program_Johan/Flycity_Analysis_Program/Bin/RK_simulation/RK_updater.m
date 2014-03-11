function [ state_new, FM_old ] = RK_updater( state, body_model, wing_model, wing_kin, wingbeat )

    i   = state.id;
    dt  = state.dt;

    % RK4 scheme:
    
    % Step 1:
    
    kine.R_strk         = wing_kin.R_strk;
    kine.dt             = dt;
    kine.id             = i*2-1;
    kine.vb             = state.vb;
    kine.wb             = state.wb;
    kine.ab             = state.ab;
    kine.w_dot_b        = state.w_dot_b;
    kine.RL             = wing_kin.RL(:,:,(i*2-1));
    kine.RR             = wing_kin.RR(:,:,(i*2-1));
    kine.wL             = wing_kin.wL(:,(i*2-1));
    kine.wR             = wing_kin.wR(:,(i*2-1));
    kine.wL_b           = wing_kin.wL_b(:,(i*2-1));
    kine.wR_b           = wing_kin.wR_b(:,(i*2-1));
    kine.w_dot_L_b      = wing_kin.w_dot_L_b(:,(i*2-1));
    kine.w_dot_R_b      = wing_kin.w_dot_R_b(:,(i*2-1));
    kine.xyz            = state.xyz;
    kine.qb             = state.qb;
    kine.Rb             = state.Rb;
    kine.v_system       = state.v_system;
    kine.alfa_L_old     = state.alfa_L_old;
    kine.alfa_R_old     = state.alfa_R_old;
    kine.g              = body_model.g;
    
    [ FMI_b, FMI_strkpln ] = Inertia_instantaneous( kine, body_model, wing_model );
    [ FMA_b, FMA_strkpln, alfa_L_old, alfa_R_old, ~, ~ ] = Aero_instantaneous( kine, body_model, wing_model, wingbeat, 1 );
    [ Fg_b, Fg_strkpln ] = Gravity_instantaneous( body_model, kine );
    
    M_mat = FMI_b.M_mat_b;
    FI    = FMI_b.F_I_vel;
    MI    = FMI_b.M_I_vel;
    FA    = FMA_b(1:3);
    MA    = FMA_b(4:6);
    Fg    = Fg_b;
    
    k1                  = 1000*(M_mat\[ FA+Fg+FI; MA+MI]);
    
    FM_old.FI_acc_b     = FMI_b.F_I_acc;
    FM_old.FI_vel_b     = FMI_b.F_I_vel;
    FM_old.MI_acc_b     = FMI_b.M_I_acc;
    FM_old.MI_vel_b     = FMI_b.M_I_vel;
    FM_old.LM           = FMI_b.LM;
    FM_old.AM           = FMI_b.AM;
    FM_old.KE           = FMI_b.KE;
    FM_old.KE_lin       = FMI_b.KE_lin;
    FM_old.KE_ang       = FMI_b.KE_ang;
    FM_old.FI_acc_strk  = FMI_strkpln.F_I_acc;
    FM_old.FI_vel_strk  = FMI_strkpln.F_I_vel;
    FM_old.MI_acc_strk  = FMI_strkpln.M_I_acc;
    FM_old.MI_vel_strk  = FMI_strkpln.M_I_vel;
    FM_old.FA_b         = FA;
    FM_old.MA_b         = MA;
    FM_old.Fg_b         = Fg;
    FM_old.FA_strk      = FMA_strkpln(1:3);
    FM_old.MA_strk      = FMA_strkpln(4:6);
    FM_old.Fg_strk      = Fg_strkpln;
    ab_old              = k1(1:3);
    w_dot_b_old         = k1(4:6);
    
    vb_0 = FMI_b.vb_0;
    wb_0 = FMI_b.wb_0;
    
    clear FM_b FM_strkpln M_mat FI MI kine xyz q_new Rb
    
    % Step 2:
        
    kine.R_strk         = wing_kin.R_strk;
    kine.dt             = dt/2;
    kine.id             = i*2;
    kine.vb             = state.vb+kine.dt*k1(1:3);
    kine.wb             = state.wb+kine.dt*k1(4:6);
    kine.ab             = k1(1:3);
    kine.w_dot_b        = k1(4:6);
    kine.RL             = wing_kin.RL(:,:,(i*2));
    kine.RR             = wing_kin.RR(:,:,(i*2));
    kine.wL             = wing_kin.wL(:,(i*2));
    kine.wR             = wing_kin.wR(:,(i*2));
    kine.wL_b           = wing_kin.wL_b(:,(i*2));
    kine.wR_b           = wing_kin.wR_b(:,(i*2));
    kine.w_dot_L_b      = wing_kin.w_dot_L_b(:,(i*2));
    kine.w_dot_R_b      = wing_kin.w_dot_R_b(:,(i*2));
    [ q_new, Rb ]       = quat_update( state.qb, kine.wb, kine.dt );
    kine.xyz            = state.xyz+kine.dt*Rb*kine.vb;
    kine.qb             = q_new;
    kine.Rb             = Rb;
    kine.v_system       = Rb*kine.vb;
    kine.alfa_L_old     = alfa_L_old;
    kine.alfa_R_old     = alfa_R_old;
    kine.g              = body_model.g;
    
    [ FMI_b, ~ ] = Inertia_instantaneous( kine, body_model, wing_model );
    [ FMA_b, ~, alfa_L_old, alfa_R_old, ~, ~ ] = Aero_instantaneous( kine, body_model, wing_model, wingbeat, 1 );
    [ Fg_b, ~ ] = Gravity_instantaneous( body_model, kine );
    
    M_mat = FMI_b.M_mat_b;
    FI    = FMI_b.F_I_vel;
    MI    = FMI_b.M_I_vel;
    FA    = FMA_b(1:3);
    MA    = FMA_b(4:6);
    Fg    = Fg_b;
    
    k2                  = 1000*(M_mat\[ FA+Fg+FI; MA+MI]);
    
    clear FM_b FM_strkpln M_mat FI MI kine xyz q_new Rb
        
    % Step 3:
    
    kine.R_strk         = wing_kin.R_strk;
    kine.dt             = dt/2;
    kine.id             = i*2;
    kine.vb             = state.vb+kine.dt*k2(1:3);
    kine.wb             = state.wb+kine.dt*k2(4:6);
    kine.ab             = k2(1:3);
    kine.w_dot_b        = k2(4:6);
    kine.RL             = wing_kin.RL(:,:,(i*2));
    kine.RR             = wing_kin.RR(:,:,(i*2));
    kine.wL             = wing_kin.wL(:,(i*2));
    kine.wR             = wing_kin.wR(:,(i*2));
    kine.wL_b           = wing_kin.wL_b(:,(i*2));
    kine.wR_b           = wing_kin.wR_b(:,(i*2));
    kine.w_dot_L_b      = wing_kin.w_dot_L_b(:,(i*2));
    kine.w_dot_R_b      = wing_kin.w_dot_R_b(:,(i*2));
    [ q_new, Rb ]       = quat_update( state.qb, kine.wb, kine.dt );
    kine.xyz            = state.xyz+kine.dt*Rb*kine.vb;
    kine.qb             = q_new;
    kine.Rb             = Rb;
    kine.v_system       = Rb*kine.vb;
    kine.alfa_L_old     = alfa_L_old;
    kine.alfa_R_old     = alfa_R_old;
    kine.g              = body_model.g;
    
    [ FMI_b, ~ ] = Inertia_instantaneous( kine, body_model, wing_model );
    [ FMA_b, ~, alfa_L_old, alfa_R_old, ~, ~ ] = Aero_instantaneous( kine, body_model, wing_model, wingbeat, 1 );
    [ Fg_b, ~ ] = Gravity_instantaneous( body_model, kine );
    
    M_mat = FMI_b.M_mat_b;
    FI    = FMI_b.F_I_vel;
    MI    = FMI_b.M_I_vel;
    FA    = FMA_b(1:3);
    MA    = FMA_b(4:6);
    Fg    = Fg_b;
    
    k3                  = 1000*(M_mat\[ FA+Fg+FI; MA+MI]);
    
    clear FM_b FM_strkpln M_mat FI MI kine xyz q_new Rb
        
    % Step 4:
    
    kine.R_strk         = wing_kin.R_strk;
    kine.dt             = dt;
    kine.id             = i*2+1;
    kine.vb             = state.vb+kine.dt*k3(1:3);
    kine.wb             = state.wb+kine.dt*k3(4:6);
    kine.ab             = k3(1:3);
    kine.w_dot_b        = k3(4:6);
    kine.RL             = wing_kin.RL(:,:,(i*2+1));
    kine.RR             = wing_kin.RR(:,:,(i*2+1));
    kine.wL             = wing_kin.wL(:,(i*2+1));
    kine.wR             = wing_kin.wR(:,(i*2+1));
    kine.wL_b           = wing_kin.wL_b(:,(i*2+1));
    kine.wR_b           = wing_kin.wR_b(:,(i*2+1));
    kine.w_dot_L_b      = wing_kin.w_dot_L_b(:,(i*2+1));
    kine.w_dot_R_b      = wing_kin.w_dot_R_b(:,(i*2+1));
    [ q_new, Rb ]       = quat_update( state.qb, kine.wb, kine.dt );
    kine.xyz            = state.xyz+kine.dt*Rb*kine.vb;
    kine.qb             = q_new;
    kine.Rb             = Rb;
    kine.v_system       = Rb*kine.vb;
    kine.alfa_L_old     = alfa_L_old;
    kine.alfa_R_old     = alfa_R_old;
    kine.g              = body_model.g;
    
    [ FMI_b, ~ ] = Inertia_instantaneous( kine, body_model, wing_model );
    [ FMA_b, ~, alfa_L_old, alfa_R_old, ~, ~ ] = Aero_instantaneous( kine, body_model, wing_model, wingbeat, 1 );
    [ Fg_b, ~ ] = Gravity_instantaneous( body_model, kine );
    
    M_mat = FMI_b.M_mat_b;
    FI    = FMI_b.F_I_vel;
    MI    = FMI_b.M_I_vel;
    FA    = FMA_b(1:3);
    MA    = FMA_b(4:6);
    Fg    = Fg_b;
    
    k4                  = 1000*(M_mat\[ FA+Fg+FI; MA+MI]);
    
    clear FM_b FM_strkpln M_mat FI MI kine xyz q_new Rb
    
    % Compute ab_new and w_dot_b_new:
    
    kine.R_strk         = wing_kin.R_strk;
    kine.dt             = dt;
    kine.id             = i*2+1;
    kine.vb             = state.vb+(dt/6)*(k1(1:3)+2*k2(1:3)+2*k3(1:3)+k4(1:3));
    kine.wb             = state.wb+(dt/6)*(k1(4:6)+2*k2(4:6)+2*k3(4:6)+k4(4:6));
    kine.ab             = k4(1:3);
    kine.w_dot_b        = k4(4:6);
    kine.RL             = wing_kin.RL(:,:,(i*2+1));
    kine.RR             = wing_kin.RR(:,:,(i*2+1));
    kine.wL             = wing_kin.wL(:,(i*2+1));
    kine.wR             = wing_kin.wR(:,(i*2+1));
    kine.wL_b           = wing_kin.wL_b(:,(i*2+1));
    kine.wR_b           = wing_kin.wR_b(:,(i*2+1));
    kine.w_dot_L_b      = wing_kin.w_dot_L_b(:,(i*2+1));
    kine.w_dot_R_b      = wing_kin.w_dot_R_b(:,(i*2+1));
    [ q_new, Rb ]       = quat_update( state.qb, kine.wb, kine.dt );
    kine.xyz            = state.xyz+kine.dt*Rb*kine.vb;
    kine.qb             = q_new;
    kine.Rb             = Rb;
    kine.v_system       = Rb*kine.vb;
    kine.alfa_L_old     = alfa_L_old;
    kine.alfa_R_old     = alfa_R_old;
    kine.g              = body_model.g;
    
    [ FMI_b, ~ ] = Inertia_instantaneous( kine, body_model, wing_model );
    [ FMA_b, ~, alfa_L_old, alfa_R_old, ~, ~ ] = Aero_instantaneous( kine, body_model, wing_model, wingbeat, 1 );
    [ Fg_b, ~ ] = Gravity_instantaneous( body_model, kine );
    
    M_mat = FMI_b.M_mat_b;
    FI    = FMI_b.F_I_vel;
    MI    = FMI_b.M_I_vel;
    FA    = FMA_b(1:3);
    MA    = FMA_b(4:6);
    Fg    = Fg_b;
    
    k5                  = 1000*(M_mat\[ FA+Fg+FI; MA+MI]);
    
    clear FM_b FM_strkpln M_mat FI MI kine xyz q_new Rb
    
    % Update state:
    
    state_new.id        = i+1;
    state_new.dt        = dt;
    state_new.vb        = state.vb+(dt/6)*(k1(1:3)+2*k2(1:3)+2*k3(1:3)+k4(1:3));
    state_new.wb        = state.wb+(dt/6)*(k1(4:6)+2*k2(4:6)+2*k3(4:6)+k4(4:6));
    state_new.ab        = k5(1:3);
    state_new.w_dot_b   = k5(4:6);
    state_new.vb_0      = vb_0;
    state_new.wb_0      = wb_0;
    state_new.ab_old            = ab_old;
    state_new.w_dot_b_old       = w_dot_b_old;
    [ q_new, Rb ]       = quat_update( state.qb, state_new.wb, state_new.dt );
    state_new.xyz       = state.xyz+state_new.dt*Rb*state_new.vb;
    state_new.qb        = q_new;
    state_new.Rb        = Rb;
    state_new.v_system  = Rb*state_new.vb;
    state_new.alfa_L_old = alfa_L_old;
    state_new.alfa_R_old = alfa_R_old;
    
    

    
end