function x_dot = state_space_func(t,x,t_func,wing_kin,body_model,wing_model,body_IC,temp_file_name,case_nr)

    % case_nr == 0: inertia forces only
    % case_nr == 1; inertia, translational aerodynamics and gravity on
    % case_nr == 2; inertia, translational and rotational aerodynamics and gravity on

    % Wing kinematics:

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
                    
    % Propagation position and orientation body:
                    
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
            xyz = xyz_old(:,end)+Rb'*(0.5*vb_old(:,end)*dt);
        elseif abs(dt) > (t_func(end)/2)
            xyz             = body_IC.xyz_0;
            qb              = body_IC.qb_0;
            Rb              = quat2mat(qb);
        else
            [ qb, Rb ]      = quat_update( qb_old(:,end), wb_old(:,end), dt );
            xyz             = xyz_old(:,end)+Rb'*(0.5*vb_old(:,end)*dt);
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
    kine.vb          = x(1:3);
    kine.wb          = x(4:6);
    kine.wL          = wL_temp;
    kine.wR          = wR_temp;
    kine.wL_b        = wL_b_temp;
    kine.wR_b        = wR_b_temp;
    kine.w_dot_L     = w_dot_L_temp;
    kine.w_dot_R     = w_dot_R_temp;
    kine.w_dot_L_b   = w_dot_L_b_temp;
    kine.w_dot_R_b   = w_dot_R_b_temp;
    kine.R_strk      = wing_kin.R_strk;
    
    % Forces and moments:
    
    if case_nr == 0
        
        [ FMI_b, ~ ] = Inertia_instantaneous( kine, body_model, wing_model );
        
        M_mat = FMI_b.M_mat_b;
        FI    = FMI_b.F_I_vel;
        MI    = FMI_b.M_I_vel;
        
        x_dot            = 1000*(M_mat\[ FI; MI]);
                
    elseif case_nr == 1
        
        [ FMI_b, ~ ] = Inertia_instantaneous( kine, body_model, wing_model );
        [ FMA_b, ~, ~, ~, ~, ~ ] = Aero_instantaneous( kine, body_model, wing_model, 0 );
        [ Fg_b,  ~ ] = Gravity_instantaneous( body_model, kine );
        
        M_mat = FMI_b.M_mat_b;
        FI    = FMI_b.F_I_vel;
        MI    = FMI_b.M_I_vel;

        FA    = FMA_b(1:3);
        MA    = FMA_b(4:6);

        Fg    = Fg_b;

        F_0   = body_IC.F_0;
        M_0   = body_IC.M_0;

        x_dot = 1000*(M_mat\[ FI+FA+Fg-F_0; MI+MA-M_0]);
        
    elseif case_nr == 2
        
        [ FMI_b, ~ ] = Inertia_instantaneous( kine, body_model, wing_model );
        [ FMA_b, ~, ~, ~, ~, ~ ] = Aero_instantaneous( kine, body_model, wing_model, 1 );
        [ Fg_b,  ~ ] = Gravity_instantaneous( body_model, kine );
        
        M_mat = FMI_b.M_mat_b;
        FI    = FMI_b.F_I_vel;
        MI    = FMI_b.M_I_vel;
        
        FA    = FMA_b(1:3);
        MA    = FMA_b(4:6);   
        
        Fg    = Fg_b;

        F_0   = body_IC.F_0;
        M_0   = body_IC.M_0;
        
        x_dot = 1000*(M_mat\[ FI+FA+Fg-F_0; MI+MA-M_0]);
        
    end

    
    vb              = x(1:3);
    wb              = x(4:6);
    
    t_old           = [t_old                t];
    xyz_old         = [xyz_old              xyz];
    qb_old          = [qb_old               qb];
    vb_old          = [vb_old               vb];
    wb_old          = [wb_old               wb];
    
    save((char(temp_file_name)),'t_old','xyz_old','qb_old','vb_old','wb_old');
    
    
end