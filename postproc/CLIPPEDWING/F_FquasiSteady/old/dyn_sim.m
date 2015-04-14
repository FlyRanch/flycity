function sim_data = dyn_sim(t_func,wing_kin,body_model,wing_model,body_IC,options,case_nr)


    nr_sect = wing_model.nr_sect;
    
    % Create a temporary .mat file to save data:
    
    temp_file_avg_wb    = 'temp_file_avg_wb.mat';
    
    if exist((char(temp_file_avg_wb))) == 2
                    
        delete((char(temp_file_avg_wb)));
        
    end
        
    % Solve the ODE:
    
    t_span = [ t_func(1) t_func(end) ];
    IC = [body_IC.vb_0; body_IC.wb_0];
    
    [ T X ] = ode45(@(t,x) state_space_func(t,x,t_func,wing_kin,body_model,wing_model,body_IC,temp_file_avg_wb,case_nr) ,t_span,IC,options);
    
    N = length(T);
    
    X_dot = zeros(6,N);
    
    for j = 1:N
    
        X_dot(:,j) = state_space_func(T(j),X(j,:)',t_func,wing_kin,body_model,wing_model,body_IC,temp_file_avg_wb,case_nr);
        
    end

    % Recover forces and moments:
    
    
    FI_vel_b        = zeros(3,N);
    MI_vel_b        = zeros(3,N);
    FI_vel_strk     = zeros(3,N);
    MI_vel_strk     = zeros(3,N);
    FI_acc_b        = zeros(3,N);
    MI_acc_b        = zeros(3,N);
    FI_acc_strk     = zeros(3,N);
    MI_acc_strk     = zeros(3,N);
    LM              = zeros(3,N);
    AM              = zeros(3,N);
    KE              = zeros(1,N);
    KE_lin          = zeros(1,N);
    KE_ang          = zeros(1,N);
    vb_0            = zeros(3,N);
    wb_0            = zeros(3,N);
    
    FA_b            = zeros(3,N);
    MA_b            = zeros(3,N);
    FA_strk         = zeros(3,N);
    MA_strk         = zeros(3,N);
    alfa_L          = zeros(nr_sect,N);
    alfa_R          = zeros(nr_sect,N);
    alfa_dot_L      = zeros(nr_sect,N);
    alfa_dot_R      = zeros(nr_sect,N);
    
    Fg_b            = zeros(3,N);
    Fg_strk         = zeros(3,N);
    
    xyz             = zeros(3,N);
    q_b             = zeros(4,N);
    
    
    for j = 1:length(T)
    
        [ pos, FMI, FMA, FMg ] = FM_recovery(T(j),X(j,:)',X_dot(:,j),t_func,wing_kin,body_model,wing_model,body_IC,temp_file_avg_wb,case_nr);
        
        FI_vel_b(:,j)        = FMI.F_vel_b;
        MI_vel_b(:,j)        = FMI.M_vel_b;
        FI_vel_strk(:,j)     = FMI.F_vel_strk;
        MI_vel_strk(:,j)     = FMI.M_vel_strk;
        FI_acc_b(:,j)        = FMI.F_acc_b;
        MI_acc_b(:,j)        = FMI.M_acc_b;
        FI_acc_strk(:,j)     = FMI.F_acc_strk;
        MI_acc_strk(:,j)     = FMI.M_acc_strk;
        LM(:,j)              = FMI.LM;
        AM(:,j)              = FMI.AM;
        KE(:,j)              = FMI.KE;
        KE_lin(:,j)          = FMI.KE_lin;
        KE_ang(:,j)          = FMI.KE_ang;
        vb_0(:,j)            = FMI.vb_0;
        wb_0(:,j)            = FMI.wb_0;

        FA_b(:,j)            = FMA.F_b;
        MA_b(:,j)            = FMA.M_b;
        FA_strk(:,j)         = FMA.F_strk;
        MA_strk(:,j)         = FMA.M_strk;
        alfa_L(:,j)          = FMA.alfa_L;
        alfa_R(:,j)          = FMA.alfa_R;
        alfa_dot_L(:,j)      = FMA.alfa_dot_L;
        alfa_dot_R(:,j)      = FMA.alfa_dot_R;

        Fg_b(:,j)            = FMg.Fg_b;
        Fg_strk(:,j)         = FMg.Fg_strk;

        xyz(:,j)             = pos.xyz;
        q_b(:,j)             = pos.q_b;
    end
    
    
    sim_data.T                  = T';
    sim_data.xyz                = xyz;
    sim_data.qb                 = q_b;
    sim_data.vb                 = X(:,1:3)';
    sim_data.wb                 = X(:,4:6)';
    sim_data.ab                 = X_dot(1:3,:);
    sim_data.w_dot_b            = X_dot(4:6,:);
    sim_data.vb_mean            = FM_mean_func(T,t_func,sim_data.vb);
    sim_data.wb_mean            = FM_mean_func(T,t_func,sim_data.wb);
    sim_data.ab_mean            = FM_mean_func(T,t_func,sim_data.ab);
    sim_data.w_dot_b_mean       = FM_mean_func(T,t_func,sim_data.w_dot_b);
    sim_data.FI_vel_b           = FI_vel_b;
    sim_data.MI_vel_b           = MI_vel_b;
    sim_data.FI_vel_strk        = FI_vel_strk;
    sim_data.MI_vel_strk        = MI_vel_strk;
    sim_data.FI_acc_b           = FI_acc_b;
    sim_data.MI_acc_b           = MI_acc_b;
    sim_data.FI_acc_strk     	= FI_acc_strk;
    sim_data.MI_acc_strk        = MI_acc_strk;
    sim_data.LM                 = LM;
    sim_data.AM                 = AM;
    sim_data.FI_vel_b_mean      = FM_mean_func(T,t_func,sim_data.FI_vel_b);
    sim_data.MI_vel_b_mean      = FM_mean_func(T,t_func,sim_data.MI_vel_b);
    sim_data.FI_vel_strk_mean   = FM_mean_func(T,t_func,sim_data.FI_vel_strk);
    sim_data.MI_vel_strk_mean   = FM_mean_func(T,t_func,sim_data.MI_vel_strk);
    sim_data.FI_acc_b_mean      = FM_mean_func(T,t_func,sim_data.FI_acc_b);
    sim_data.MI_acc_b_mean      = FM_mean_func(T,t_func,sim_data.MI_acc_b);
    sim_data.FI_acc_strk_mean   = FM_mean_func(T,t_func,sim_data.FI_acc_strk);
    sim_data.MI_acc_strk_mean   = FM_mean_func(T,t_func,sim_data.MI_acc_strk);
    sim_data.LM_mean            = FM_mean_func(T,t_func,sim_data.LM);
    sim_data.AM_mean            = FM_mean_func(T,t_func,sim_data.AM);
    sim_data.KE                 = KE;
    sim_data.KE_lin             = KE_lin;
    sim_data.KE_ang             = KE_ang;
    sim_data.vb_0               = vb_0;
    sim_data.wb_0               = wb_0;
    sim_data.FA_b               = FA_b;
    sim_data.MA_b               = MA_b;
    sim_data.FA_strk            = FA_strk;
    sim_data.MA_strk            = MA_strk;
    sim_data.alfa_L             = alfa_L;
    sim_data.alfa_R             = alfa_R;
    sim_data.alfa_dot_L         = alfa_dot_L;
    sim_data.alfa_dot_R         = alfa_dot_R;
    sim_data.FA_b_mean          = FM_mean_func(T,t_func,sim_data.FA_b);
    sim_data.MA_b_mean          = FM_mean_func(T,t_func,sim_data.MA_b);
    sim_data.FA_strk_mean       = FM_mean_func(T,t_func,sim_data.FA_strk);
    sim_data.MA_strk_mean       = FM_mean_func(T,t_func,sim_data.MA_strk);
    sim_data.Fg_b               = Fg_b;
    sim_data.Fg_strk            = Fg_strk;
    sim_data.Fg_b_mean          = FM_mean_func(T,t_func,sim_data.Fg_b);
    sim_data.Fg_strk_mean       = FM_mean_func(T,t_func,sim_data.Fg_strk);
    sim_data.F_0                = body_IC.F_0;
    sim_data.M_0                = body_IC.M_0;
    
end