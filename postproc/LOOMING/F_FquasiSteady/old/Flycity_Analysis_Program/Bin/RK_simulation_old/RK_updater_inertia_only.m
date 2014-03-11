function [ state_new ] = RK_updater_inertia_only( state, wing_kin, body_model, wing_model, wingbeat )

    
    % Use 4th order Runge-Kutta scheme to update the state vector from i to
    % i+1:
    
   
    % Step 1: -------------------------------------------------------------
    
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
    
    k1 = nan(12,1);
    
    k1(1:3) = kine.vb;
    k1(4:6) = kine.wb;
    k1(7:12) = 1000*(M_mat\[ FI; MI]);
    
    clear FI MI M_mat FMI_b
    
    
    % Step 2: -------------------------------------------------------------
    
    
    kine.dt             = state.dt/2;
    kine.id             = state.id;   
    
    dt_t = kine.dt;
    qb_t = state.qb;
    wb_t = k1(4:6);
    
    PHI = [ qb_t(4)*(dt_t/2)   -qb_t(3)*(dt_t/2)    qb_t(2)*(dt_t/2)    1   wb_t(3)*(dt_t/2)   -wb_t(2)*(dt_t/2)    wb_t(1)*(dt_t/2) ; ...
            qb_t(3)*(dt_t/2)    qb_t(4)*(dt_t/2)   -qb_t(1)*(dt_t/2)   -wb_t(3)*(dt_t/2)    1   wb_t(1)*(dt_t/2)    wb_t(2)*(dt_t/2) ; ...
           -qb_t(2)*(dt_t/2)    qb_t(1)*(dt_t/2)    qb_t(4)*(dt_t/2)    wb_t(2)*(dt_t/2)   -wb_t(1)*(dt_t/2)    1   wb_t(3)*(dt_t/2) ; ...
           -qb_t(1)*(dt_t/2)   -qb_t(2)*(dt_t/2)   -qb_t(3)*(dt_t/2)   -wb_t(1)*(dt_t/2)   -wb_t(2)*(dt_t/2)   -wb_t(3)*(dt_t/2)   1 ];
       
    qb_n = PHI*[wb_t(1); wb_t(2); wb_t(3); qb_t(1); qb_t(2); qb_t(3); qb_t(4)];
    
    kine.qb             = qb_n./norm(qb_n);
    kine.Rb             = quat2mat(kine.qb);
    kine.xyz            = state.xyz+kine.Rb'*k1(1:3)*kine.dt;
    kine.vb             = state.vb+kine.dt*k1(7:9);
    kine.wb             = state.wb+kine.dt*k1(10:12);
    kine.ab             = k1(7:9);
    kine.w_dot_b        = k1(10:12);
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
    
    k2 = nan(12,1);
    
    k2(1:3) = kine.vb;
    k2(4:6) = kine.wb;
    k2(7:12) = 1000*(M_mat\[ FI; MI]);
    
    clear FI MI M_mat FMI_b
    
    % Step 3: -------------------------------------------------------------
    
    kine.dt             = state.dt/2;
    kine.id             = state.id;
    
    dt_t = kine.dt;
    qb_t = state.qb;
    wb_t = k2(4:6);
    
    PHI = [ qb_t(4)*(dt_t/2)   -qb_t(3)*(dt_t/2)    qb_t(2)*(dt_t/2)    1   wb_t(3)*(dt_t/2)   -wb_t(2)*(dt_t/2)    wb_t(1)*(dt_t/2) ; ...
            qb_t(3)*(dt_t/2)    qb_t(4)*(dt_t/2)   -qb_t(1)*(dt_t/2)   -wb_t(3)*(dt_t/2)    1   wb_t(1)*(dt_t/2)    wb_t(2)*(dt_t/2) ; ...
           -qb_t(2)*(dt_t/2)    qb_t(1)*(dt_t/2)    qb_t(4)*(dt_t/2)    wb_t(2)*(dt_t/2)   -wb_t(1)*(dt_t/2)    1   wb_t(3)*(dt_t/2) ; ...
           -qb_t(1)*(dt_t/2)   -qb_t(2)*(dt_t/2)   -qb_t(3)*(dt_t/2)   -wb_t(1)*(dt_t/2)   -wb_t(2)*(dt_t/2)   -wb_t(3)*(dt_t/2)   1 ];
       
    qb_n = PHI*[wb_t(1); wb_t(2); wb_t(3); qb_t(1); qb_t(2); qb_t(3); qb_t(4)];
    
    kine.qb             = qb_n./norm(qb_n);
    kine.Rb             = quat2mat(kine.qb);
    kine.xyz            = state.xyz+kine.Rb'*k2(1:3)*kine.dt;
    kine.vb             = state.vb+kine.dt*k2(7:9);
    kine.wb             = state.wb+kine.dt*k2(10:12);
    kine.ab             = k2(7:9);
    kine.w_dot_b        = k2(10:12);
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
    
    k3 = nan(12,1);
    
    k3(1:3) = kine.vb;
    k3(4:6) = kine.wb;
    k3(7:12) = 1000*(M_mat\[ FI; MI]);
    
    clear FI MI M_mat FMI_b
    
    
    
    % Step 4: -------------------------------------------------------------
    
    kine.dt             = state.dt;
    kine.id             = state.id+1;
    
    dt_t = kine.dt;
    qb_t = state.qb;
    wb_t = k3(4:6);
    
    PHI = [ qb_t(4)*(dt_t/2)   -qb_t(3)*(dt_t/2)    qb_t(2)*(dt_t/2)    1   wb_t(3)*(dt_t/2)   -wb_t(2)*(dt_t/2)    wb_t(1)*(dt_t/2) ; ...
            qb_t(3)*(dt_t/2)    qb_t(4)*(dt_t/2)   -qb_t(1)*(dt_t/2)   -wb_t(3)*(dt_t/2)    1   wb_t(1)*(dt_t/2)    wb_t(2)*(dt_t/2) ; ...
           -qb_t(2)*(dt_t/2)    qb_t(1)*(dt_t/2)    qb_t(4)*(dt_t/2)    wb_t(2)*(dt_t/2)   -wb_t(1)*(dt_t/2)    1   wb_t(3)*(dt_t/2) ; ...
           -qb_t(1)*(dt_t/2)   -qb_t(2)*(dt_t/2)   -qb_t(3)*(dt_t/2)   -wb_t(1)*(dt_t/2)   -wb_t(2)*(dt_t/2)   -wb_t(3)*(dt_t/2)   1 ];
       
    qb_n = PHI*[wb_t(1); wb_t(2); wb_t(3); qb_t(1); qb_t(2); qb_t(3); qb_t(4)];
    
    kine.qb             = qb_n./norm(qb_n);
    kine.Rb             = quat2mat(kine.qb);
    kine.xyz            = state.xyz+kine.Rb'*k3(1:3)*kine.dt;
    kine.vb             = state.vb+kine.dt*k3(7:9);
    kine.wb             = state.wb+kine.dt*k3(10:12);
    kine.ab             = k3(7:9);
    kine.w_dot_b        = k3(10:12);
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
    
    k4 = nan(12,1);
    
    k4(1:3) = kine.vb;
    k4(4:6) = kine.wb;
    k4(7:12) = 1000*(M_mat\[ FI; MI]);
    
    clear FI MI M_mat FMI_b
    
    
    
    % Compute new state with coefficients k1, k2, k3 and k4: --------------
    
    % Compute ab and w_dot_b;
    
    kine.dt             = state.dt;
    kine.id             = state.id+1;
    
    dt_t = kine.dt;
    qb_t = state.qb;
    wb_t = (1/6)*(k1(4:6)+2*k2(4:6)+2*k3(4:6)+k4(4:6));
    
    PHI = [ qb_t(4)*(dt_t/2)   -qb_t(3)*(dt_t/2)    qb_t(2)*(dt_t/2)    1   wb_t(3)*(dt_t/2)   -wb_t(2)*(dt_t/2)    wb_t(1)*(dt_t/2) ; ...
            qb_t(3)*(dt_t/2)    qb_t(4)*(dt_t/2)   -qb_t(1)*(dt_t/2)   -wb_t(3)*(dt_t/2)    1   wb_t(1)*(dt_t/2)    wb_t(2)*(dt_t/2) ; ...
           -qb_t(2)*(dt_t/2)    qb_t(1)*(dt_t/2)    qb_t(4)*(dt_t/2)    wb_t(2)*(dt_t/2)   -wb_t(1)*(dt_t/2)    1   wb_t(3)*(dt_t/2) ; ...
           -qb_t(1)*(dt_t/2)   -qb_t(2)*(dt_t/2)   -qb_t(3)*(dt_t/2)   -wb_t(1)*(dt_t/2)   -wb_t(2)*(dt_t/2)   -wb_t(3)*(dt_t/2)   1 ];
       
    qb_n = PHI*[wb_t(1); wb_t(2); wb_t(3); qb_t(1); qb_t(2); qb_t(3); qb_t(4)];
    
    kine.qb             = qb_n./norm(qb_n);
    kine.Rb             = quat2mat(kine.qb);
    kine.xyz            = state.xyz+kine.Rb'*(1/6)*(k1(1:3)+2*k2(1:3)+2*k3(1:3)+k4(1:3))*state.dt;
    kine.vb             = state.vb+(1/6)*(k1(7:9)+2*k2(7:9)+2*k3(7:9)+k4(7:9))*state.dt;
    kine.wb             = state.wb+(1/6)*(k1(10:12)+2*k2(10:12)+2*k3(10:12)+k4(10:12))*state.dt;
    kine.ab             = k4(7:9);
    kine.w_dot_b        = k4(10:12);
    kine.RL             = wing_kin.RL(:,:,(kine.id*2-1));
    kine.RR             = wing_kin.RR(:,:,(kine.id*2-1));
    kine.wL             = wing_kin.wL(:,(kine.id*2-1));
    kine.wR             = wing_kin.wR(:,(kine.id*2-1));
    kine.w_dot_L        = wing_kin.w_dot_L(:,(kine.id*2-1));
    kine.w_dot_R        = wing_kin.w_dot_R(:,(kine.id*2-1));
    
    [ FMI_b, FMI_strkpln ] = Inertia_instantaneous( kine, body_model, wing_model );
    
    M_mat = FMI_b.M_mat_b;
    FI    = FMI_b.F_I_vel;
    MI    = FMI_b.M_I_vel;
    
    k5 = 1000*(M_mat\[ FI; MI]);
    
    % Update state:
    
    state_new.id         = state.id+1;
    state_new.xyz        = kine.xyz;
    state_new.qb         = kine.qb;
    state_new.Rb         = kine.Rb;
    state_new.vb         = kine.vb;
    state_new.wb         = kine.wb;      
    state_new.ab         = k5(1:3);
    state_new.w_dot_b    = k5(4:6);
    state_new.FMI_b      = FMI_b;
    state_new.FMI_strkpln= FMI_strkpln;
%     state_new.LM         = FMI_b.LM;
%     state_new.AM         = FMI_b.AM;
%     state_new.KE         = FMI_b.KE;
%     state_new.KE_lin     = FMI_b.KE_lin;
%     state_new.KE_ang     = FMI_b.KE_ang;
    
    state_new.FI_vel     = FMI_b.F_I_vel;
    state_new.MI_vel     = FMI_b.M_I_vel;
    
    state_new.wL_sim      = kine.wL;
    state_new.wR_sim      = kine.wR;
    state_new.w_dot_L_sim = kine.w_dot_L;
    state_new.w_dot_R_sim = kine.w_dot_R;
    
end

