function [ state_new ] = RK_updater( state, wing_kin, body_model, wing_model, wingbeat )

    
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
    kine.alfa_L_old     = state.alfa_L_old;
    kine.alfa_R_old     = state.alfa_R_old;
    kine.g              = body_model.g;
        
    
    [ FMA_b, FMA_strkpln, alfa_L, alfa_R, alfa_dot_L, alfa_dot_R ] = Aero_instantaneous( kine, body_model, wing_model, wingbeat, 1 );
    [ FMI_b, FMI_strkpln ] = Inertia_instantaneous( kine, body_model, wing_model );
    [ Fg_b, Fg_strkpln ] = Gravity_instantaneous( body_model, kine );
    
    FMA_b_t         = FMA_b;
    FMA_strkpln_t   = FMA_strkpln;
    FMI_b_t         = FMI_b;
    FMI_strkpln_t   = FMI_strkpln;
    Fg_b_t          = Fg_b;
    Fg_strkpln_t    = Fg_strkpln;
    
    M_mat = FMI_b.M_mat_b;
    FA    = FMA_b(1:3);
    FI    = FMI_b.F_I_vel;
    Fg    = Fg_b;
    MA    = FMA_b(4:6);
    MI    = FMI_b.M_I_vel;
    
    k1 = nan(12,1);
    
    k1(1:3) = kine.vb;
    k1(4:6) = kine.wb;
    k1(7:12) = 1000*(M_mat\[ FA+FI+Fg; MA+MI]);
    
    clear FMA_b FMI_b Fg_b M_mat FA FI Fg MA MI
    
    
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
    kine.alfa_L_old     = alfa_L;
    kine.alfa_R_old     = alfa_R;
    
    [ FMA_b, ~, ~, ~, ~, ~ ] = Aero_instantaneous( kine, body_model, wing_model, wingbeat, 1 );
    [ FMI_b, ~ ] = Inertia_instantaneous( kine, body_model, wing_model );
    [ Fg_b, ~ ] = Gravity_instantaneous( body_model, kine );
    
    M_mat = FMI_b.M_mat_b;
    FA    = FMA_b(1:3);
    FI    = FMI_b.F_I_vel;
    Fg    = Fg_b;
    MA    = FMA_b(4:6);
    MI    = FMI_b.M_I_vel;
    
    k2 = nan(12,1);
    
    k2(1:3) = kine.vb;
    k2(4:6) = kine.wb;
    k2(7:12) = 1000*(M_mat\[ FA+FI+Fg; MA+MI]);
    
    clear FMA_b FMI_b Fg_b M_mat FA FI Fg MA MI PHI

    
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
    kine.alfa_L_old     = alfa_L;
    kine.alfa_R_old     = alfa_R;
    
    [ FMA_b, ~, ~, ~, ~, ~ ] = Aero_instantaneous( kine, body_model, wing_model, wingbeat, 1 );
    [ FMI_b, ~ ] = Inertia_instantaneous( kine, body_model, wing_model );
    [ Fg_b, ~ ] = Gravity_instantaneous( body_model, kine );
    
    M_mat = FMI_b.M_mat_b;
    FA    = FMA_b(1:3);
    FI    = FMI_b.F_I_vel;
    Fg    = Fg_b;
    MA    = FMA_b(4:6);
    MI    = FMI_b.M_I_vel;
    
    k3 = nan(12,1);
    
    k3(1:3) = kine.vb;
    k3(4:6) = kine.wb;
    k3(7:12) = 1000*(M_mat\[ FA+FI+Fg; MA+MI]);
    
    clear FMA_b FMI_b Fg_b M_mat FA FI Fg MA MI PHI
    
    
    
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
    kine.alfa_L_old     = alfa_L;
    kine.alfa_R_old     = alfa_R;
    
    [ FMA_b, ~, ~, ~, ~, ~ ] = Aero_instantaneous( kine, body_model, wing_model, wingbeat, 1 );
    [ FMI_b, ~ ] = Inertia_instantaneous( kine, body_model, wing_model );
    [ Fg_b, ~ ] = Gravity_instantaneous( body_model, kine );
    
    M_mat = FMI_b.M_mat_b;
    FA    = FMA_b(1:3);
    FI    = FMI_b.F_I_vel;
    Fg    = Fg_b;
    MA    = FMA_b(4:6);
    MI    = FMI_b.M_I_vel;
    
    k4 = nan(12,1);
    
    k4(1:3) = kine.vb;
    k4(4:6) = kine.wb;
    k4(7:12) = 1000*(M_mat\[ FA+FI+Fg; MA+MI]);
    
    clear FMA_b FMI_b Fg_b M_mat FA FI Fg MA MI PHI
    
    
    
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
    kine.alfa_L_old     = alfa_L;
    kine.alfa_R_old     = alfa_R;
    
    [ FMA_b, ~, ~, ~, ~, ~ ] = Aero_instantaneous( kine, body_model, wing_model, wingbeat, 1 );
    [ FMI_b, ~ ] = Inertia_instantaneous( kine, body_model, wing_model );
    [ Fg_b, ~ ] = Gravity_instantaneous( body_model, kine );
    
    M_mat = FMI_b.M_mat_b;
    FA    = FMA_b(1:3);
    FI    = FMI_b.F_I_vel;
    Fg    = Fg_b;
    MA    = FMA_b(4:6);
    MI    = FMI_b.M_I_vel;
    
        
    k5 = 1000*(M_mat\[ FA+FI+Fg; MA+MI]);
    
    
    % Update state:
    
    state_new.id         = state.id+1;
    state_new.xyz        = kine.xyz;
    state_new.qb         = kine.qb;
    state_new.Rb         = kine.Rb;
    state_new.vb         = kine.vb;
    state_new.wb         = kine.wb;      
    state_new.ab         = k5(1:3);
    state_new.w_dot_b    = k5(4:6);
    state_new.alfa_L_old = alfa_L;
    state_new.alfa_R_old = alfa_R;
    state_new.alfa_dot_L = alfa_dot_L;
    state_new.alfa_dot_R = alfa_dot_R;
    state_new.FMA_b      = FMA_b_t;
    state_new.FMA_strkpln= FMA_strkpln_t;
    state_new.FMI_b      = FMI_b_t;
    state_new.FMI_strkpln= FMI_strkpln_t;
    state_new.Fg_b       = Fg_b_t;
    state_new.Fg_strkpln = Fg_strkpln_t;
    
end

