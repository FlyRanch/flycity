function [ state_new ] = RK2_inertia_only( state, wing_kin, body_model, wing_model, wingbeat )

    if state.id == 1
        
        kine.R_strk         = wingbeat.R_strk;
        kine.dt             = state.dt;

        kine.id             = state.id;
        kine.xyz            = state.xyz;
        kine.qb             = state.qb;
        kine.Rb             = state.Rb;
        kine.vb             = state.vb;
        kine.wb             = state.wb;
        kine.ab             = state.ab;
        kine.w_dot_b        = state.w_dot_b;
        kine.RL             = wing_kin.RL(:,:,(kine.id*2-1));
        kine.RR             = wing_kin.RR(:,:,(kine.id*2-1));
        kine.wL             = wing_kin.wL(:,(kine.id*2-1));
        kine.wR             = wing_kin.wR(:,(kine.id*2-1));
        kine.w_dot_L        = wing_kin.w_dot_L(:,(kine.id*2-1));
        kine.w_dot_R        = wing_kin.w_dot_R(:,(kine.id*2-1));

        [ FMI_b, ~ ] = Inertia_instantaneous( kine, body_model, wing_model );

        M_mat = FMI_b.M_mat_b;
        FI    = FMI_b.F_I_vel;
        MI    = FMI_b.M_I_vel;

        k3 = 1000*(M_mat\[ FI; MI]);

        clear FI MI M_mat
        
        ab_0      = k3(1:3);
        w_dot_b_0 = k3(4:6);
        
        state_new.ab_0 = ab_0;
        state_new.w_dot_b_0 = w_dot_b_0;
        
    end


    % Step 1: -------------------------------------------------------------
    
    if state.id > 1
        
        ab_0 = state.ab_0;
        w_dot_b_0 = state.w_dot_b_0;
        
    end
    
    kine.R_strk         = wingbeat.R_strk;
    kine.dt             = state.dt;
    
    kine.id             = state.id;
    kine.xyz            = state.xyz;
    kine.qb             = state.qb;
    kine.Rb             = state.Rb;
    kine.vb             = state.vb;
    kine.wb             = state.wb;
    kine.ab             = state.ab;
    kine.w_dot_b        = state.w_dot_b;
    kine.RL             = wing_kin.RL(:,:,(kine.id*2-1));
    kine.RR             = wing_kin.RR(:,:,(kine.id*2-1));
    kine.wL             = wing_kin.wL(:,(kine.id*2-1));
    kine.wR             = wing_kin.wR(:,(kine.id*2-1));
    kine.w_dot_L        = wing_kin.w_dot_L(:,(kine.id*2-1));
    kine.w_dot_R        = wing_kin.w_dot_R(:,(kine.id*2-1));
        
    [ FMI_b, ~ ] = Inertia_instantaneous( kine, body_model, wing_model );
    
    M_mat = FMI_b.M_mat_b;
    FI    = FMI_b.F_I_vel;
    MI    = FMI_b.M_I_vel;
    
    k1 = 1000*(M_mat\[ FI; MI]);
    
    clear FI MI M_mat
    
    
    
    % Step 2: -------------------------------------------------------------
    
    kine.R_strk         = wingbeat.R_strk;
    kine.dt             = state.dt/2;
    
    kine.id             = state.id;
    kine.xyz            = state.xyz;
    kine.qb             = state.qb;
    kine.Rb             = state.Rb;
    kine.vb             = state.vb+0.5*k1(1:3)*state.dt;
    kine.wb             = state.wb+0.5*k1(4:6)*state.dt;
    kine.ab             = k1(1:3);
    kine.w_dot_b        = k1(4:6);
    kine.RL             = wing_kin.RL(:,:,(kine.id*2));
    kine.RR             = wing_kin.RR(:,:,(kine.id*2));
    kine.wL             = wing_kin.wL(:,(kine.id*2));
    kine.wR             = wing_kin.wR(:,(kine.id*2));
    kine.w_dot_L        = wing_kin.w_dot_L(:,(kine.id*2));
    kine.w_dot_R        = wing_kin.w_dot_R(:,(kine.id*2));
        
    [ FMI_b, ~ ] = Inertia_instantaneous( kine, body_model, wing_model );
    
    M_mat = FMI_b.M_mat_b;
    FI    = FMI_b.F_I_vel;
    MI    = FMI_b.M_I_vel;
    
    k2 = 1000*(M_mat\[ FI; MI]);
    
    clear FI MI M_mat
    
    % Step 3, compute new state:
    
    kine.R_strk         = wingbeat.R_strk;
    kine.dt             = state.dt;
    
    kine.id             = state.id;
    kine.xyz            = state.xyz;
    kine.qb             = state.qb;
    kine.Rb             = state.Rb;
    kine.vb             = state.vb+0.5*k2(1:3)*state.dt;
    kine.wb             = state.wb+0.5*k2(4:6)*state.dt;
    kine.ab             = k2(1:3);
    kine.w_dot_b        = k2(4:6);
    kine.RL             = wing_kin.RL(:,:,(kine.id*2));
    kine.RR             = wing_kin.RR(:,:,(kine.id*2));
    kine.wL             = wing_kin.wL(:,(kine.id*2));
    kine.wR             = wing_kin.wR(:,(kine.id*2));
    kine.w_dot_L        = wing_kin.w_dot_L(:,(kine.id*2));
    kine.w_dot_R        = wing_kin.w_dot_R(:,(kine.id*2));
        
    [ FMI_b, ~ ] = Inertia_instantaneous( kine, body_model, wing_model );
    
    M_mat = FMI_b.M_mat_b;
    FI    = FMI_b.F_I_vel;
    MI    = FMI_b.M_I_vel;
    
    k3 = 1000*(M_mat\[ FI; MI]);
    
    clear FI MI M_mat

    % Update the state:
    
    state_new.id         = state.id+1;
    state_new.xyz        = kine.xyz;
    state_new.qb         = kine.qb;
    state_new.Rb         = kine.Rb;
%     state_new.vb         = state.vb+0.5*k2(1:3)*state.dt;
%     state_new.wb         = state.wb+0.5*k2(4:6)*state.dt;
    state_new.ab         = k3(1:3);
    state_new.w_dot_b    = k3(4:6);
    state_new.vb         = kine.vb;
    state_new.wb         = kine.wb;
    
    
    


end

