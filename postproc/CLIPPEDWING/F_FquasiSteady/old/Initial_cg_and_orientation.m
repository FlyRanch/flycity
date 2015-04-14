function [ F_0, M_0 ] = Initial_cg_and_orientation( body_model, wing_model, a_avg, vb_0, wb_0)
    
    R_strk = a_avg.R_strk;
    
    [ kine ] = angular_velocities_polynomial( a_avg );

    t_func              = kine.t;
    kine.R_strk         = R_strk;
    kine.down_up_t      = t_func(end)*a_avg.down_up;
    
    body_IC.xyz_0   = [0; 0; 0];
    body_IC.qb_0    = quat2mat((R_strk'*[1 0 0; 0 -1 0; 0 0 -1]))';
    body_IC.vb_0    = vb_0;
    body_IC.wb_0    = wb_0;
    body_IC.F_0     = [0; 0; 0];
    body_IC.M_0     = [0; 0; 0];

    clear sim_data

    % Tolerance settings ode45:

    options = odeset('RelTol',1e-5,'AbsTol',[1e-15 1e-15 1e-15 1e-15 1e-15 1e-15]);

    % case_nr == 0: Inertia only
    % case_nr == 1: Inertia, translational aerodynamics, gravity
    % case_nr == 2: Inertia, translational and rotational aerodynamics,
    % gravity

    case_nr = 2;

    sim_data = dyn_sim(t_func,kine,body_model,wing_model,body_IC,options,case_nr);
    
    F_tot1  = sim_data.FI_acc_b_mean;
    M_tot1  = sim_data.MI_acc_b_mean;
    FA     = sim_data.FA_b_mean;
    MA     = sim_data.MA_b_mean;
    FI     = sim_data.FI_vel_b_mean;
    MI     = sim_data.MI_vel_b_mean;
    Fg     = sim_data.Fg_b_mean;
    
    clear sim_data

    
    body_IC.xyz_0   = [0; 0; 0];
    body_IC.qb_0    = quat2mat((R_strk'*[1 0 0; 0 -1 0; 0 0 -1]))';
    body_IC.vb_0    = vb_0;
    body_IC.wb_0    = wb_0;
    body_IC.F_0     = F_tot1;
    body_IC.M_0     = M_tot1;

    sim_data = dyn_sim(t_func,kine,body_model,wing_model,body_IC,options,case_nr);
    
    F_tot2  = sim_data.FI_acc_b_mean;
    M_tot2  = sim_data.MI_acc_b_mean;
    FA     = sim_data.FA_b_mean;
    MA     = sim_data.MA_b_mean;
    FI     = sim_data.FI_vel_b_mean;
    MI     = sim_data.MI_vel_b_mean;
    Fg     = sim_data.Fg_b_mean;
    
    F_0 = F_tot1+F_tot2;
    M_0 = M_tot1+M_tot2;


end


%     x_cg = (M_tot(2)+(FA(1)+FI(1))*JL(3))/(FA(3)+FI(3));
%     
%     sup_ratio = sqrt(FA(1)^2+FA(3)^2)/sqrt(Fg(1)^2+Fg(3)^2);
%     
%     m_fly = body_model.mass_fly;
%     
%     m_fly_new = m_fly;
%     
%     m_body_new = m_fly_new-2*wing_model.mass;
%     
%     theta = asin(cross([FA(1)+FI(1); 0; FA(3)+FI(3)],[Fg(1); 0; Fg(3)])/(norm([FA(1)+FI(1); 0; FA(3)+FI(3)])*norm([Fg(1); 0; Fg(3)])))
%     
%     settings.beta_strk
%     
%     beta = settings.beta_strk+theta(2)
%     
%     R_0 = [cos(beta) 0 -sin(beta); ...
%            0         1 0; ...
%            sin(beta) 0 cos(beta)];
%        
%     cg_b = [x_cg-JL(1); 0; 0];
%        
%     D_body = [ (cg_b(2)^2+cg_b(3)^2) -cg_b(1)*cg_b(2) -cg_b(1)*cg_b(3); ...
%                -cg_b(1)*cg_b(2) (cg_b(1)^2+cg_b(3)^2) -cg_b(2)*cg_b(3); ...
%                -cg_b(1)*cg_b(3) -cg_b(2)*cg_b(3) (cg_b(1)^2+cg_b(2)^2)];
%        
%     Inertia_new = body_model.Inertia; %+m_body_new*D_body;
%        
%     body_model_new.mass_fly         = m_fly_new;
%     body_model_new.mass_body        = body_model.mass_body;
%     body_model_new.Joint_left       = [x_cg; JL(2); JL(3)];
%     body_model_new.Joint_right      = [x_cg; -JL(2); JL(3)];
%     body_model_new.cg_b             = body_model.cg_b;
%     body_model_new.Inertia          = Inertia_new;
%     body_model_new.x_mod            = body_model.x_mod;
%     body_model_new.y_mod            = body_model.y_mod;
%     body_model_new.z_mod            = body_model.z_mod;
%     body_model_new.g                = 1e-3*settings.g;
