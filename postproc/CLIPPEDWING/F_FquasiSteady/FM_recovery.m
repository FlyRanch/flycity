function [ pos, FMI, FMA, FMg ] = FM_recovery(t,X,X_dot,t_func,wing_kin,body_model,wing_model,body_IC,temp_file_name,case_nr)

    qL_temp          = [interp1(t_func,wing_kin.qL(1,:),t,'spline'); ...
                        interp1(t_func,wing_kin.qL(2,:),t,'spline'); ...
                        interp1(t_func,wing_kin.qL(3,:),t,'spline'); ...
                        interp1(t_func,wing_kin.qL(4,:),t,'spline')];
                    
    qR_temp          = [interp1(t_func,wing_kin.qR(1,:),t,'spline'); ...
                        interp1(t_func,wing_kin.qR(2,:),t,'spline'); ...
                        interp1(t_func,wing_kin.qR(3,:),t,'spline'); ...
                        interp1(t_func,wing_kin.qR(4,:),t,'spline')];
                    
    wL_temp          = [interp1(t_func,wing_kin.wL(1,:),t,'spline'); ...
                        interp1(t_func,wing_kin.wL(2,:),t,'spline'); ...
                        interp1(t_func,wing_kin.wL(3,:),t,'spline')];
    
    wR_temp          = [interp1(t_func,wing_kin.wR(1,:),t,'spline'); ...
                        interp1(t_func,wing_kin.wR(2,:),t,'spline'); ...
                        interp1(t_func,wing_kin.wR(3,:),t,'spline')];
    
    w_dot_L_temp     = [interp1(t_func,wing_kin.w_dot_L(1,:),t,'spline'); ...
                        interp1(t_func,wing_kin.w_dot_L(2,:),t,'spline'); ...
                        interp1(t_func,wing_kin.w_dot_L(3,:),t,'spline')];
    
    w_dot_R_temp     = [interp1(t_func,wing_kin.w_dot_R(1,:),t,'spline'); ...
                        interp1(t_func,wing_kin.w_dot_R(2,:),t,'spline'); ...
                        interp1(t_func,wing_kin.w_dot_R(3,:),t,'spline')];
    
    wL_b_temp        = [interp1(t_func,wing_kin.wL_b(1,:),t,'spline'); ...
                        interp1(t_func,wing_kin.wL_b(2,:),t,'spline'); ...
                        interp1(t_func,wing_kin.wL_b(3,:),t,'spline')];
    
    wR_b_temp        = [interp1(t_func,wing_kin.wR_b(1,:),t,'spline'); ...
                        interp1(t_func,wing_kin.wR_b(2,:),t,'spline'); ...
                        interp1(t_func,wing_kin.wR_b(3,:),t,'spline')];
    
    w_dot_L_b_temp   = [interp1(t_func,wing_kin.w_dot_L_b(1,:),t,'spline'); ...
                        interp1(t_func,wing_kin.w_dot_L_b(2,:),t,'spline'); ...
                        interp1(t_func,wing_kin.w_dot_L_b(3,:),t,'spline')];
    
    w_dot_R_b_temp   = [interp1(t_func,wing_kin.w_dot_R_b(1,:),t,'spline'); ...
                        interp1(t_func,wing_kin.w_dot_R_b(2,:),t,'spline'); ...
                        interp1(t_func,wing_kin.w_dot_R_b(3,:),t,'spline')];
                    
    nr_sect = wing_model.nr_sect;
                    
    if exist((char(temp_file_name))) == 2
        
        temp = load((char(temp_file_name)));

        t_old           = temp.t_old;
        xyz_old         = temp.xyz_old;
        qb_old          = temp.qb_old;
        vb_old          = temp.vb_old;
        wb_old          = temp.wb_old;
        
        dt              = t-t_old(end);
        
        if abs(dt) < 1e-12
            qb = qb_old(:,end);
            Rb = quat2mat(qb);
            xyz = xyz_old(:,end)+0.5*vb_old(:,end)*dt;
        elseif abs(dt) > (t_func(end)/2)
            xyz             = body_IC.xyz_0;
            qb              = body_IC.qb_0;
            Rb              = quat2mat(qb);
        else
            [ qb, Rb ]      = quat_update( qb_old(:,end), wb_old(:,end), dt );
            xyz             = xyz_old(:,end)+0.5*vb_old(:,end)*dt;
        end
        
        clear temp
        
        delete((char(temp_file_name)));
        
    else
        
       
        t_old           = [];
        xyz_old         = [];
        qb_old          = [];
        vb_old          = [];
        wb_old          = [];
        
        xyz             = body_IC.xyz_0;
        qb              = body_IC.qb_0;
        Rb              = quat2mat(qb);
        
    end
    
    
    % kine:
                       
    kine.qL          = qL_temp/norm(qL_temp);
    kine.qR          = qR_temp/norm(qR_temp);
    kine.RL          = quat2mat(kine.qL);
    kine.RR          = quat2mat(kine.qR);
    kine.Rb          = Rb;
    kine.vb          = X(1:3);
    kine.wb          = X(4:6);
    kine.wL          = wL_temp;
    kine.wR          = wR_temp;
    kine.wL_b        = wL_b_temp;
    kine.wR_b        = wR_b_temp;
    kine.w_dot_L     = w_dot_L_temp;
    kine.w_dot_R     = w_dot_R_temp;
    kine.w_dot_L_b   = w_dot_L_b_temp;
    kine.w_dot_R_b   = w_dot_R_b_temp;
    kine.R_strk      = wing_kin.R_strk;
    
    
    if case_nr == 0
        
        [ FMI_b, FMI_strk ] = Inertia_instantaneous( kine, body_model, wing_model );
        
        M_mat = FMI_b.M_mat_b;

        FMI_acc = 1e-3*M_mat*X_dot;

        FMI.F_vel_b         = FMI_b.F_I_vel;
        FMI.M_vel_b         = FMI_b.M_I_vel;
        FMI.F_vel_strk      = FMI_strk.F_I_vel;
        FMI.M_vel_strk      = FMI_strk.M_I_vel;
        FMI.F_acc_b         = FMI_acc(1:3);
        FMI.M_acc_b         = FMI_acc(4:6);
        FMI.F_acc_strk      = kine.R_strk*FMI_acc(1:3);
        FMI.M_acc_strk      = kine.R_strk*FMI_acc(4:6);
        FMI.LM              = FMI_b.LM;
        FMI.AM              = FMI_b.AM;
        FMI.KE              = FMI_b.KE;
        FMI.KE_lin          = FMI_b.KE_lin;
        FMI.KE_ang          = FMI_b.KE_ang;
        FMI.vb_0            = FMI_b.vb_0;
        FMI.wb_0            = FMI_b.wb_0;

        FMA.F_b             = zeros(3,1);
        FMA.M_b             = zeros(3,1);
        FMA.F_strk          = zeros(3,1);
        FMA.M_strk          = zeros(3,1);
        FMA.alfa_L          = zeros(nr_sect,1);
        FMA.alfa_R          = zeros(nr_sect,1);
        FMA.alfa_dot_L      = zeros(nr_sect,1);
        FMA.alfa_dot_R      = zeros(nr_sect,1);

        FMg.Fg_b            = zeros(3,1);
        FMg.Fg_strk         = zeros(3,1);
                        
    elseif case_nr == 1
        
        [ FMI_b, FMI_strk ] = Inertia_instantaneous( kine, body_model, wing_model );
        [ FMA_b, FMA_strk, alfa_L, alfa_R, alfa_dot_L, alfa_dot_R ] = Aero_instantaneous( kine, body_model, wing_model, 0 );
        [ Fg_b,  Fg_strk ] = Gravity_instantaneous( body_model, kine );
        
        M_mat = FMI_b.M_mat_b;

        FMI_acc = 1e-3*M_mat*X_dot;

        FMI.F_vel_b         = FMI_b.F_I_vel;
        FMI.M_vel_b         = FMI_b.M_I_vel;
        FMI.F_vel_strk      = FMI_strk.F_I_vel;
        FMI.M_vel_strk      = FMI_strk.M_I_vel;
        FMI.F_acc_b         = FMI_acc(1:3);
        FMI.M_acc_b         = FMI_acc(4:6);
        FMI.F_acc_strk      = kine.R_strk*FMI_acc(1:3);
        FMI.M_acc_strk      = kine.R_strk*FMI_acc(4:6);
        FMI.LM              = FMI_b.LM;
        FMI.AM              = FMI_b.AM;
        FMI.KE              = FMI_b.KE;
        FMI.KE_lin          = FMI_b.KE_lin;
        FMI.KE_ang          = FMI_b.KE_ang;
        FMI.vb_0            = FMI_b.vb_0;
        FMI.wb_0            = FMI_b.wb_0;

        FMA.F_b             = FMA_b(1:3);
        FMA.M_b             = FMA_b(4:6);
        FMA.F_strk          = FMA_strk(1:3);
        FMA.M_strk          = FMA_strk(4:6);
        FMA.alfa_L          = alfa_L;
        FMA.alfa_R          = alfa_R;
        FMA.alfa_dot_L      = alfa_dot_L;
        FMA.alfa_dot_R      = alfa_dot_R;

        FMg.Fg_b            = Fg_b;
        FMg.Fg_strk         = Fg_strk;
        
    elseif case_nr == 2
        
        [ FMI_b, FMI_strk ] = Inertia_instantaneous( kine, body_model, wing_model );
        [ FMA_b, FMA_strk, alfa_L, alfa_R, alfa_dot_L, alfa_dot_R ] = Aero_instantaneous( kine, body_model, wing_model, 1 );
        [ Fg_b,  Fg_strk ] = Gravity_instantaneous( body_model, kine );
        
        M_mat = FMI_b.M_mat_b;

        FMI_acc = 1e-3*M_mat*X_dot;

        FMI.F_vel_b         = FMI_b.F_I_vel;
        FMI.M_vel_b         = FMI_b.M_I_vel;
        FMI.F_vel_strk      = FMI_strk.F_I_vel;
        FMI.M_vel_strk      = FMI_strk.M_I_vel;
        FMI.F_acc_b         = FMI_acc(1:3);
        FMI.M_acc_b         = FMI_acc(4:6);
        FMI.F_acc_strk      = kine.R_strk*FMI_acc(1:3);
        FMI.M_acc_strk      = kine.R_strk*FMI_acc(4:6);
        FMI.LM              = FMI_b.LM;
        FMI.AM              = FMI_b.AM;
        FMI.KE              = FMI_b.KE;
        FMI.KE_lin          = FMI_b.KE_lin;
        FMI.KE_ang          = FMI_b.KE_ang;
        FMI.vb_0            = FMI_b.vb_0;
        FMI.wb_0            = FMI_b.wb_0;

        FMA.F_b             = FMA_b(1:3);
        FMA.M_b             = FMA_b(4:6);
        FMA.F_strk          = FMA_strk(1:3);
        FMA.M_strk          = FMA_strk(4:6);
        FMA.alfa_L          = alfa_L;
        FMA.alfa_R          = alfa_R;
        FMA.alfa_dot_L      = alfa_dot_L;
        FMA.alfa_dot_R      = alfa_dot_R;

        FMg.Fg_b            = Fg_b;
        FMg.Fg_strk         = Fg_strk;
        
    end
    
    pos.xyz             = xyz;
    pos.q_b             = qb;
    
    vb              = X(1:3);
    wb              = X(4:6);
    
    t_old           = [t_old                t];
    xyz_old         = [xyz_old              xyz];
    qb_old          = [qb_old               qb];
    vb_old          = [vb_old               vb];
    wb_old          = [wb_old               wb];
    
    save((char(temp_file_name)),'t_old','xyz_old','qb_old','vb_old','wb_old');

end