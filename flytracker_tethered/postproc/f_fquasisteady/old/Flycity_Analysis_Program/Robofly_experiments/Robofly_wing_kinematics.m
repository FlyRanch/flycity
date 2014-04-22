function Robofly_wing_kinematics( settings, pathDB )

    savefile1 = 'pathDB7.mat';
    
    savefile2 = 'Robofly_maneuvering_tests.mat';

    % Determine the range of non-dimensional forces and moments of the
    % data-set:
    
    nr_steps = 20;
    nr_wb_exp = 7;
    nr_points = 101;
    
    theta_max_range     = 30;
    phi_max_range       = 90;
    
    Fx_forward_factor   = 1.2;
    Fx_back_factor      = 1.2;
    Fy_factor           = 1.2;
    Fz_down_factor      = 1.2;
    Fz_up_factor        = 1.2;
    Mx_factor           = 1.2;
    My_up_factor        = 1.2;
    My_down_factor      = 1.2;
    Mz_factor           = 1.2;
    
    FM_ax =  pathDB.maneuver.c_fit_ax.FM.FM_1;
    FM_ay =  pathDB.maneuver.c_fit_ay.FM.FM_2;
    FM_az =  pathDB.maneuver.c_fit_az.FM.FM_3;
    FM_wx =  pathDB.maneuver.c_fit_wx.FM.FM_4;
    FM_wy =  pathDB.maneuver.c_fit_wy.FM.FM_5;
    FM_wz =  pathDB.maneuver.c_fit_wz.FM.FM_6;
    
    N_Fx_forward    = pathDB.maneuver.c_fit_ax.N_p;
    N_Fx_back       = pathDB.maneuver.c_fit_ax.N_n;
    N_Fy            = pathDB.maneuver.c_fit_ay.N;
    N_Fz_down       = pathDB.maneuver.c_fit_az.N_n;
    N_Fz_up         = pathDB.maneuver.c_fit_az.N_p;
    N_Mx            = pathDB.maneuver.c_fit_wx.N;
    N_My_up         = pathDB.maneuver.c_fit_wy.N_n;
    N_My_down       = pathDB.maneuver.c_fit_wy.N_p;
    N_Mz            = pathDB.maneuver.c_fit_wz.N;
    
    
    
    N_Fx_forward
    N_Fx_back
    N_Fy
    N_Fz_down
    N_Fz_up
    N_Mx
    N_My_up
    N_My_down
    N_Mz
    
    Fx_forward_max  = 0;
    Fx_back_max     = 0;
    Fy_max          = 0;
    Fz_down_max     = 0;
    Fz_up_max       = 0;
    Mx_max          = 0;
    My_up_max       = 0;
    My_down_max     = 0;
    Mz_max          = 0;
    
    if max(FM_ax) > 0    
        Fx_forward_max  = max(FM_ax);
    end
    if min(FM_ax) < 0
        Fx_back_max     = min(FM_ax);
    end
    Fy_max              = max(abs(FM_ay));
    if max(FM_az) > 0
        Fz_down_max     = max(FM_az);
    end
    if min(FM_az) < 0
        Fz_up_max       = min(FM_az);
    end
    Mx_max              = max(abs(FM_wx));
    if max(FM_wy) > 0
        My_up_max       = max(FM_wy);
    end
    if min(FM_wy) < 0
        My_down_max     = min(FM_wy);
    end
    Mz_max              = max(abs(FM_wz));
    
    
    
    if abs(Fx_forward_max) > 0
        Fx_forward_range    = 0:((Fx_forward_max*Fx_forward_factor)/(nr_steps-1)):(Fx_forward_max*Fx_forward_factor);
    else
        Fx_forward_range    = zeros(1,nr_steps);
    end
    if abs(Fx_back_max) > 0
        Fx_back_range       = 0:((Fx_back_max*Fx_back_factor)/(nr_steps-1)):(Fx_back_max*Fx_back_factor);
    else
        Fx_back_range       = zeros(1,nr_steps);
    end
    if abs(Fy_max) > 0
        Fy_range            = 0:((Fy_max*Fy_factor)/(nr_steps-1)):(Fy_max*Fy_factor);
    else
        Fy_range            = zeros(1,nr_steps);
    end
    if abs(Fz_up_max) > 0
        Fz_up_range         = 0:((Fz_up_max*Fz_up_factor)/(nr_steps-1)):(Fz_up_max*Fz_up_factor);
    else
        Fz_up_range         = zeros(1,nr_steps);
    end
    if abs(Fz_down_max) > 0
        Fz_down_range       = 0:((Fz_down_max*Fz_down_factor)/(nr_steps-1)):(Fz_down_max*Fz_down_factor);
    else
        Fz_down_range       = zeros(1,nr_steps);
    end
    if abs(Mx_max) > 0
        Mx_range            = 0:((Mx_max*Mx_factor)/(nr_steps-1)):(Mx_max*Mx_factor);
    else
        Mx_range            = zeros(1,nr_steps);
    end
    if abs(My_up_max) > 0
        My_up_range         = 0:((My_up_max*My_up_factor)/(nr_steps-1)):(My_up_max*My_up_factor);
    else
        My_up_range         = zeros(1,nr_steps);
    end
    if abs(My_down_max) > 0
        My_down_range       = 0:((My_down_max*My_down_factor)/(nr_steps-1)):(My_down_max*My_down_factor);
    else
        My_down_range       = zeros(1,nr_steps);
    end
    if abs(Mz_max) > 0
        Mz_range            = 0:((Mz_max*Mz_factor)/(nr_steps-1)):(Mz_max*Mz_factor);
    else
        Mz_range            = zeros(1,nr_steps);
    end
    
    
    % Set up the maneuvering wing kinematics:
    
       
    n_pol_theta     = settings.n_pol_theta;
    n_pol_eta       = settings.n_pol_eta;
    n_pol_phi       = settings.n_pol_phi;
    
    a_theta_glob    = pathDB.poly_fit.a_glob.theta;
    a_eta_glob      = pathDB.poly_fit.a_glob.eta;
    a_phi_glob      = pathDB.poly_fit.a_glob.phi;
    
    m_fly_glob      = pathDB.poly_fit.a_glob.m_fly;
    wing_l_glob     = pathDB.poly_fit.a_glob.wing_l;
    down_up_glob    = pathDB.poly_fit.a_glob.down_up;
    f_glob          = pathDB.poly_fit.a_glob.f;
    
    a_theta_Fx_forward_L    = zeros(n_pol_theta*2+2,nr_steps);
    a_eta_Fx_forward_L      = zeros(n_pol_eta*2+2,nr_steps);
    a_phi_Fx_forward_L      = zeros(n_pol_phi*2+2,nr_steps);
    
    a_theta_Fx_forward_R    = zeros(n_pol_theta*2+2,nr_steps);
    a_eta_Fx_forward_R      = zeros(n_pol_eta*2+2,nr_steps);
    a_phi_Fx_forward_R      = zeros(n_pol_phi*2+2,nr_steps);
    
    a_theta_Fx_back_L       = zeros(n_pol_theta*2+2,nr_steps);
    a_eta_Fx_back_L         = zeros(n_pol_eta*2+2,nr_steps);
    a_phi_Fx_back_L         = zeros(n_pol_phi*2+2,nr_steps);
    
    a_theta_Fx_back_R       = zeros(n_pol_theta*2+2,nr_steps);
    a_eta_Fx_back_R         = zeros(n_pol_eta*2+2,nr_steps);
    a_phi_Fx_back_R         = zeros(n_pol_phi*2+2,nr_steps);
    
    a_theta_Fy_L            = zeros(n_pol_theta*2+2,nr_steps);
    a_eta_Fy_L              = zeros(n_pol_eta*2+2,nr_steps);
    a_phi_Fy_L              = zeros(n_pol_phi*2+2,nr_steps);
    
    a_theta_Fy_R            = zeros(n_pol_theta*2+2,nr_steps);
    a_eta_Fy_R              = zeros(n_pol_eta*2+2,nr_steps);
    a_phi_Fy_R              = zeros(n_pol_phi*2+2,nr_steps);
    
    a_theta_Fz_down_L       = zeros(n_pol_theta*2+2,nr_steps);
    a_eta_Fz_down_L         = zeros(n_pol_eta*2+2,nr_steps);
    a_phi_Fz_down_L         = zeros(n_pol_phi*2+2,nr_steps);
    
    a_theta_Fz_down_R       = zeros(n_pol_theta*2+2,nr_steps);
    a_eta_Fz_down_R         = zeros(n_pol_eta*2+2,nr_steps);
    a_phi_Fz_down_R         = zeros(n_pol_phi*2+2,nr_steps);
    
    a_theta_Fz_up_L         = zeros(n_pol_theta*2+2,nr_steps);
    a_eta_Fz_up_L           = zeros(n_pol_eta*2+2,nr_steps);
    a_phi_Fz_up_L           = zeros(n_pol_phi*2+2,nr_steps);
    
    a_theta_Fz_up_R         = zeros(n_pol_theta*2+2,nr_steps);
    a_eta_Fz_up_R           = zeros(n_pol_eta*2+2,nr_steps);
    a_phi_Fz_up_R           = zeros(n_pol_phi*2+2,nr_steps);
    
    a_theta_Mx_L            = zeros(n_pol_theta*2+2,nr_steps);
    a_eta_Mx_L              = zeros(n_pol_eta*2+2,nr_steps);
    a_phi_Mx_L              = zeros(n_pol_phi*2+2,nr_steps);
    
    a_theta_Mx_R            = zeros(n_pol_theta*2+2,nr_steps);
    a_eta_Mx_R              = zeros(n_pol_eta*2+2,nr_steps);
    a_phi_Mx_R              = zeros(n_pol_phi*2+2,nr_steps);
    
    a_theta_My_up_L         = zeros(n_pol_theta*2+2,nr_steps);
    a_eta_My_up_L           = zeros(n_pol_eta*2+2,nr_steps);
    a_phi_My_up_L           = zeros(n_pol_phi*2+2,nr_steps);
    
    a_theta_My_up_R         = zeros(n_pol_theta*2+2,nr_steps);
    a_eta_My_up_R           = zeros(n_pol_eta*2+2,nr_steps);
    a_phi_My_up_R           = zeros(n_pol_phi*2+2,nr_steps);

    a_theta_My_down_L       = zeros(n_pol_theta*2+2,nr_steps);
    a_eta_My_down_L         = zeros(n_pol_eta*2+2,nr_steps);
    a_phi_My_down_L         = zeros(n_pol_phi*2+2,nr_steps);
    
    a_theta_My_down_R       = zeros(n_pol_theta*2+2,nr_steps);
    a_eta_My_down_R         = zeros(n_pol_eta*2+2,nr_steps);
    a_phi_My_down_R         = zeros(n_pol_phi*2+2,nr_steps);
    
    a_theta_Mz_L            = zeros(n_pol_theta*2+2,nr_steps);
    a_eta_Mz_L              = zeros(n_pol_eta*2+2,nr_steps);
    a_phi_Mz_L              = zeros(n_pol_phi*2+2,nr_steps);
    
    a_theta_Mz_R            = zeros(n_pol_theta*2+2,nr_steps);
    a_eta_Mz_R              = zeros(n_pol_eta*2+2,nr_steps);
    a_phi_Mz_R              = zeros(n_pol_phi*2+2,nr_steps);
    
    b_theta_Fx_forward      = pathDB.maneuver.c_fit_ax.b_theta_p(:,1,1);
    b_eta_Fx_forward        = pathDB.maneuver.c_fit_ax.b_eta_p(:,1,1);
    b_phi_Fx_forward        = pathDB.maneuver.c_fit_ax.b_phi_p(:,1,1);
    
    b_theta_Fx_back         = pathDB.maneuver.c_fit_ax.b_theta_n(:,1,1);
    b_eta_Fx_back           = pathDB.maneuver.c_fit_ax.b_eta_n(:,1,1);
    b_phi_Fx_back           = pathDB.maneuver.c_fit_ax.b_phi_n(:,1,1);
    
    b_theta_Fy_L            = pathDB.maneuver.c_fit_ay.b_theta_L(:,1,2);
    b_eta_Fy_L              = pathDB.maneuver.c_fit_ay.b_eta_L(:,1,2);
    b_phi_Fy_L              = pathDB.maneuver.c_fit_ay.b_phi_L(:,1,2);
    
    b_theta_Fy_R            = pathDB.maneuver.c_fit_ay.b_theta_R(:,1,2);
    b_eta_Fy_R              = pathDB.maneuver.c_fit_ay.b_eta_R(:,1,2);
    b_phi_Fy_R              = pathDB.maneuver.c_fit_ay.b_phi_R(:,1,2);
    
    b_theta_Fz_down         = pathDB.maneuver.c_fit_az.b_theta_p(:,1,3);
    b_eta_Fz_down           = pathDB.maneuver.c_fit_az.b_eta_p(:,1,3);
    b_phi_Fz_down           = pathDB.maneuver.c_fit_az.b_phi_p(:,1,3);
    
    b_theta_Fz_up           = pathDB.maneuver.c_fit_az.b_theta_n(:,1,3);
    b_eta_Fz_up             = pathDB.maneuver.c_fit_az.b_eta_n(:,1,3);
    b_phi_Fz_up             = pathDB.maneuver.c_fit_az.b_phi_n(:,1,3);
    
    b_theta_Mx_L            = pathDB.maneuver.c_fit_wx.b_theta_L(:,1,4);
    b_eta_Mx_L              = pathDB.maneuver.c_fit_wx.b_eta_L(:,1,4);
    b_phi_Mx_L              = pathDB.maneuver.c_fit_wx.b_phi_L(:,1,4);
    
    b_theta_Mx_R            = pathDB.maneuver.c_fit_wx.b_theta_R(:,1,4);
    b_eta_Mx_R              = pathDB.maneuver.c_fit_wx.b_eta_R(:,1,4);
    b_phi_Mx_R              = pathDB.maneuver.c_fit_wx.b_phi_R(:,1,4);
    
    b_theta_My_up           = pathDB.maneuver.c_fit_wy.b_theta_p(:,1,5);
    b_eta_My_up             = pathDB.maneuver.c_fit_wy.b_eta_p(:,1,5);
    b_phi_My_up             = pathDB.maneuver.c_fit_wy.b_phi_p(:,1,5);

    b_theta_My_down         = pathDB.maneuver.c_fit_wy.b_theta_n(:,1,5);
    b_eta_My_down           = pathDB.maneuver.c_fit_wy.b_eta_n(:,1,5);
    b_phi_My_down           = pathDB.maneuver.c_fit_wy.b_phi_n(:,1,5);
    
    b_theta_Mz_L            = pathDB.maneuver.c_fit_wz.b_theta_L(:,1,6);
    b_eta_Mz_L              = pathDB.maneuver.c_fit_wz.b_eta_L(:,1,6);
    b_phi_Mz_L              = pathDB.maneuver.c_fit_wz.b_phi_L(:,1,6);
    
    b_theta_Mz_R            = pathDB.maneuver.c_fit_wz.b_theta_R(:,1,6);
    b_eta_Mz_R              = pathDB.maneuver.c_fit_wz.b_eta_R(:,1,6);
    b_phi_Mz_R              = pathDB.maneuver.c_fit_wz.b_phi_R(:,1,6);

    
    for i = 1:nr_steps
        
        a_theta_Fx_forward_L(:,i)    = a_theta_glob+b_theta_Fx_forward*Fx_forward_range(i);
        a_eta_Fx_forward_L(:,i)      = a_eta_glob+b_eta_Fx_forward*Fx_forward_range(i);
        a_phi_Fx_forward_L(:,i)      = a_phi_glob+b_phi_Fx_forward*Fx_forward_range(i);

        a_theta_Fx_forward_R(:,i)    = a_theta_glob+b_theta_Fx_forward*Fx_forward_range(i);
        a_eta_Fx_forward_R(:,i)      = a_eta_glob+b_eta_Fx_forward*Fx_forward_range(i);
        a_phi_Fx_forward_R(:,i)      = a_phi_glob+b_phi_Fx_forward*Fx_forward_range(i);

        a_theta_Fx_back_L(:,i)       = a_theta_glob+b_theta_Fx_back*Fx_back_range(i);
        a_eta_Fx_back_L(:,i)         = a_eta_glob+b_eta_Fx_back*Fx_back_range(i);
        a_phi_Fx_back_L(:,i)         = a_phi_glob+b_phi_Fx_back*Fx_back_range(i);

        a_theta_Fx_back_R(:,i)       = a_theta_glob+b_theta_Fx_back*Fx_back_range(i);
        a_eta_Fx_back_R(:,i)         = a_eta_glob+b_eta_Fx_back*Fx_back_range(i);
        a_phi_Fx_back_R(:,i)         = a_phi_glob+b_phi_Fx_back*Fx_back_range(i);

        a_theta_Fy_L(:,i)            = a_theta_glob+b_theta_Fy_L*Fy_range(i);
        a_eta_Fy_L(:,i)              = a_eta_glob+b_eta_Fy_L*Fy_range(i);
        a_phi_Fy_L(:,i)              = a_phi_glob+b_phi_Fy_L*Fy_range(i);

        a_theta_Fy_R(:,i)            = a_theta_glob+b_theta_Fy_R*Fy_range(i);
        a_eta_Fy_R(:,i)              = a_eta_glob+b_eta_Fy_R*Fy_range(i);
        a_phi_Fy_R(:,i)              = a_phi_glob+b_phi_Fy_R*Fy_range(i);

        a_theta_Fz_down_L(:,i)       = a_theta_glob+b_theta_Fz_down*Fz_down_range(i);
        a_eta_Fz_down_L(:,i)         = a_eta_glob+b_eta_Fz_down*Fz_down_range(i);
        a_phi_Fz_down_L(:,i)         = a_phi_glob+b_phi_Fz_down*Fz_down_range(i);

        a_theta_Fz_down_R(:,i)       = a_theta_glob+b_theta_Fz_down*Fz_down_range(i);
        a_eta_Fz_down_R(:,i)         = a_eta_glob+b_eta_Fz_down*Fz_down_range(i);
        a_phi_Fz_down_R(:,i)         = a_phi_glob+b_phi_Fz_down*Fz_down_range(i);

        a_theta_Fz_up_L(:,i)         = a_theta_glob+b_theta_Fz_up*Fz_up_range(i);
        a_eta_Fz_up_L(:,i)           = a_eta_glob+b_eta_Fz_up*Fz_up_range(i);
        a_phi_Fz_up_L(:,i)           = a_phi_glob+b_phi_Fz_up*Fz_up_range(i);

        a_theta_Fz_up_R(:,i)         = a_theta_glob+b_theta_Fz_up*Fz_up_range(i);
        a_eta_Fz_up_R(:,i)           = a_eta_glob+b_eta_Fz_up*Fz_up_range(i);
        a_phi_Fz_up_R(:,i)           = a_phi_glob+b_phi_Fz_up*Fz_up_range(i);

        a_theta_Mx_L(:,i)            = a_theta_glob+b_theta_Mx_L*Mx_range(i);
        a_eta_Mx_L(:,i)              = a_eta_glob+b_eta_Mx_L*Mx_range(i);
        a_phi_Mx_L(:,i)              = a_phi_glob+b_phi_Mx_L*Mx_range(i);

        a_theta_Mx_R(:,i)            = a_theta_glob+b_theta_Mx_R*Mx_range(i);
        a_eta_Mx_R(:,i)              = a_eta_glob+b_eta_Mx_R*Mx_range(i);
        a_phi_Mx_R(:,i)              = a_phi_glob+b_phi_Mx_R*Mx_range(i);

        a_theta_My_up_L(:,i)         = a_theta_glob+b_theta_My_up*My_up_range(i);
        a_eta_My_up_L(:,i)           = a_eta_glob+b_eta_My_up*My_up_range(i);
        a_phi_My_up_L(:,i)           = a_phi_glob+b_phi_My_up*My_up_range(i);
        
        a_theta_My_up_R(:,i)         = a_theta_glob+b_theta_My_up*My_up_range(i);
        a_eta_My_up_R(:,i)           = a_eta_glob+b_eta_My_up*My_up_range(i);
        a_phi_My_up_R(:,i)           = a_phi_glob+b_phi_My_up*My_up_range(i);

        a_theta_My_down_L(:,i)       = a_theta_glob+b_theta_My_down*My_down_range(i);
        a_eta_My_down_L(:,i)         = a_eta_glob+b_eta_My_down*My_down_range(i);
        a_phi_My_down_L(:,i)         = a_phi_glob+b_phi_My_down*My_down_range(i);

        a_theta_My_down_R(:,i)       = a_theta_glob+b_theta_My_down*My_down_range(i);
        a_eta_My_down_R(:,i)         = a_eta_glob+b_eta_My_down*My_down_range(i);
        a_phi_My_down_R(:,i)         = a_phi_glob+b_phi_My_down*My_down_range(i);

        a_theta_Mz_L(:,i)            = a_theta_glob+b_theta_Mz_L*Mz_range(i);
        a_eta_Mz_L(:,i)              = a_eta_glob+b_eta_Mz_L*Mz_range(i);
        a_phi_Mz_L(:,i)              = a_phi_glob+b_phi_Mz_L*Mz_range(i);

        a_theta_Mz_R(:,i)            = a_theta_glob+b_theta_Mz_R*Mz_range(i);
        a_eta_Mz_R(:,i)              = a_eta_glob+b_eta_Mz_R*Mz_range(i);
        a_phi_Mz_R(:,i)              = a_phi_glob+b_phi_Mz_R*Mz_range(i);
        
    end
    
    
    % Plot the maneuvering wing kinematics:
    
    t = nan(1,nr_wb_exp*(nr_points-1)+1);
    
    theta_Fx_forward_L  = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    eta_Fx_forward_L    = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    phi_Fx_forward_L    = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);

    theta_Fx_forward_R  = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    eta_Fx_forward_R    = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    phi_Fx_forward_R    = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);

    theta_Fy_L          = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    eta_Fy_L            = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    phi_Fy_L            = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);

    theta_Fy_R          = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    eta_Fy_R            = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    phi_Fy_R            = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);

    theta_Fx_back_L     = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    eta_Fx_back_L       = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    phi_Fx_back_L       = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);

    theta_Fx_back_R     = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    eta_Fx_back_R       = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    phi_Fx_back_R       = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);

    theta_Fz_down_L     = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    eta_Fz_down_L       = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    phi_Fz_down_L       = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);

    theta_Fz_down_R     = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    eta_Fz_down_R       = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    phi_Fz_down_R       = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);

    theta_Fz_up_L       = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    eta_Fz_up_L         = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    phi_Fz_up_L         = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);

    theta_Fz_up_R       = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    eta_Fz_up_R         = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    phi_Fz_up_R         = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);

    theta_Mx_L          = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    eta_Mx_L            = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    phi_Mx_L            = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);

    theta_Mx_R          = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    eta_Mx_R            = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    phi_Mx_R            = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);

    theta_My_up_L       = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    eta_My_up_L         = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    phi_My_up_L         = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);

    theta_My_up_R       = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    eta_My_up_R         = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    phi_My_up_R         = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);

    theta_My_down_L     = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    eta_My_down_L       = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    phi_My_down_L       = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);

    theta_My_down_R     = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    eta_My_down_R       = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    phi_My_down_R       = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);

    theta_Mz_L          = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    eta_Mz_L            = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    phi_Mz_L            = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);

    theta_Mz_R          = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    eta_Mz_R            = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    phi_Mz_R            = nan(nr_steps,nr_wb_exp*(nr_points-1)+1);
    
    
    
    theta_Fx_forward_L_max  = nan(nr_steps,1);
    phi_Fx_forward_L_max    = nan(nr_steps,1);

    theta_Fx_forward_R_max  = nan(nr_steps,1);
    phi_Fx_forward_R_max    = nan(nr_steps,1);

    theta_Fy_L_max          = nan(nr_steps,1);
    phi_Fy_L_max            = nan(nr_steps,1);

    theta_Fy_R_max          = nan(nr_steps,1);
    phi_Fy_R_max            = nan(nr_steps,1);

    theta_Fx_back_L_max     = nan(nr_steps,1);
    phi_Fx_back_L_max       = nan(nr_steps,1);

    theta_Fx_back_R_max     = nan(nr_steps,1);
    phi_Fx_back_R_max       = nan(nr_steps,1);
    
    theta_Fz_down_L_max     = nan(nr_steps,1);
    phi_Fz_down_L_max       = nan(nr_steps,1);

    theta_Fz_down_R_max     = nan(nr_steps,1);
    phi_Fz_down_R_max       = nan(nr_steps,1);

    theta_Fz_up_L_max       = nan(nr_steps,1);
    phi_Fz_up_L_max         = nan(nr_steps,1);

    theta_Fz_up_R_max       = nan(nr_steps,1);
    phi_Fz_up_R_max         = nan(nr_steps,1);

    theta_Mx_L_max          = nan(nr_steps,1);
    phi_Mx_L_max            = nan(nr_steps,1);

    theta_Mx_R_max          = nan(nr_steps,1);
    phi_Mx_R_max            = nan(nr_steps,1);

    theta_My_up_L_max       = nan(nr_steps,1);
    phi_My_up_L_max         = nan(nr_steps,1);

    theta_My_up_R_max       = nan(nr_steps,1);
    phi_My_up_R_max         = nan(nr_steps,1);

    theta_My_down_L_max     = nan(nr_steps,1);
    phi_My_down_L_max       = nan(nr_steps,1);

    theta_My_down_R_max     = nan(nr_steps,1);
    phi_My_down_R_max       = nan(nr_steps,1);
    
    theta_Mz_L_max          = nan(nr_steps,1);
    phi_Mz_L_max            = nan(nr_steps,1);

    theta_Mz_R_max          = nan(nr_steps,1);
    phi_Mz_R_max            = nan(nr_steps,1);
    
    
    
    for k = 1:nr_wb_exp

        [ time, X_theta ]  = Wingbeat_Legendre_matrix( n_pol_theta, down_up_glob, nr_points, 0, 1/f_glob, 0 );
        [ ~, X_eta ]    = Wingbeat_Legendre_matrix( n_pol_eta, down_up_glob, nr_points, 0, 1/f_glob, 0 );
        [ ~, X_phi ]    = Wingbeat_Legendre_matrix( n_pol_phi, down_up_glob, nr_points, 0, 1/f_glob, 0 );

        k_range = ((k-1)*(nr_points-1)+1):((k-1)*(nr_points-1)+nr_points);
        
        if k == 1
            
            t(k_range) = time;
            
        else

            t(k_range) = t((k-1)*(nr_points-1)+1)+time;
        
        end
            
        for j = 1:nr_steps

            theta_Fx_forward_L(j,k_range)  = X_theta*a_theta_Fx_forward_L(:,j);
            eta_Fx_forward_L(j,k_range)    = X_eta*a_eta_Fx_forward_L(:,j);
            phi_Fx_forward_L(j,k_range)    = X_phi*a_phi_Fx_forward_L(:,j);

            theta_Fx_forward_R(j,k_range)  = X_theta*a_theta_Fx_forward_R(:,j);
            eta_Fx_forward_R(j,k_range)    = X_eta*a_eta_Fx_forward_R(:,j);
            phi_Fx_forward_R(j,k_range)    = X_phi*a_phi_Fx_forward_R(:,j);

            theta_Fx_back_L(j,k_range)     = X_theta*a_theta_Fx_back_L(:,j);
            eta_Fx_back_L(j,k_range)       = X_eta*a_eta_Fx_back_L(:,j);
            phi_Fx_back_L(j,k_range)       = X_phi*a_phi_Fx_back_L(:,j);

            theta_Fx_back_R(j,k_range)     = X_theta*a_theta_Fx_back_R(:,j);
            eta_Fx_back_R(j,k_range)       = X_eta*a_eta_Fx_back_R(:,j);
            phi_Fx_back_R(j,k_range)       = X_phi*a_phi_Fx_back_R(:,j);

            theta_Fy_L(j,k_range)          = X_theta*a_theta_Fy_L(:,j);
            eta_Fy_L(j,k_range)            = X_eta*a_eta_Fy_L(:,j);
            phi_Fy_L(j,k_range)            = X_phi*a_phi_Fy_L(:,j);

            theta_Fy_R(j,k_range)          = X_theta*a_theta_Fy_R(:,j);
            eta_Fy_R(j,k_range)            = X_eta*a_eta_Fy_R(:,j);
            phi_Fy_R(j,k_range)            = X_phi*a_phi_Fy_R(:,j);

            theta_Fz_down_L(j,k_range)     = X_theta*a_theta_Fz_down_L(:,j);
            eta_Fz_down_L(j,k_range)       = X_eta*a_eta_Fz_down_L(:,j);
            phi_Fz_down_L(j,k_range)       = X_phi*a_phi_Fz_down_L(:,j);

            theta_Fz_down_R(j,k_range)     = X_theta*a_theta_Fz_down_R(:,j);
            eta_Fz_down_R(j,k_range)       = X_eta*a_eta_Fz_down_R(:,j);
            phi_Fz_down_R(j,k_range)       = X_phi*a_phi_Fz_down_R(:,j);

            theta_Fz_up_L(j,k_range)       = X_theta*a_theta_Fz_up_L(:,j);
            eta_Fz_up_L(j,k_range)         = X_eta*a_eta_Fz_up_L(:,j);
            phi_Fz_up_L(j,k_range)         = X_phi*a_phi_Fz_up_L(:,j);

            theta_Fz_up_R(j,k_range)       = X_theta*a_theta_Fz_up_R(:,j);
            eta_Fz_up_R(j,k_range)         = X_eta*a_eta_Fz_up_R(:,j);
            phi_Fz_up_R(j,k_range)         = X_phi*a_phi_Fz_up_R(:,j);

            theta_Mx_L(j,k_range)          = X_theta*a_theta_Mx_L(:,j);
            eta_Mx_L(j,k_range)            = X_eta*a_eta_Mx_L(:,j);
            phi_Mx_L(j,k_range)            = X_phi*a_phi_Mx_L(:,j);

            theta_Mx_R(j,k_range)          = X_theta*a_theta_Mx_R(:,j);
            eta_Mx_R(j,k_range)            = X_eta*a_eta_Mx_R(:,j);
            phi_Mx_R(j,k_range)            = X_phi*a_phi_Mx_R(:,j);

            theta_My_up_L(j,k_range)       = X_theta*a_theta_My_up_L(:,j);
            eta_My_up_L(j,k_range)         = X_eta*a_eta_My_up_L(:,j);
            phi_My_up_L(j,k_range)         = X_phi*a_phi_My_up_L(:,j);

            theta_My_up_R(j,k_range)       = X_theta*a_theta_My_up_R(:,j);
            eta_My_up_R(j,k_range)         = X_eta*a_eta_My_up_R(:,j);
            phi_My_up_R(j,k_range)         = X_phi*a_phi_My_up_R(:,j);

            theta_My_down_L(j,k_range)     = X_theta*a_theta_My_down_L(:,j);
            eta_My_down_L(j,k_range)       = X_eta*a_eta_My_down_L(:,j);
            phi_My_down_L(j,k_range)       = X_phi*a_phi_My_down_L(:,j);

            theta_My_down_R(j,k_range)      = X_theta*a_theta_My_down_R(:,j);
            eta_My_down_R(j,k_range)       = X_eta*a_eta_My_down_R(:,j);
            phi_My_down_R(j,k_range)       = X_phi*a_phi_My_down_R(:,j);

            theta_Mz_L(j,k_range)          = X_theta*a_theta_Mz_L(:,j);
            eta_Mz_L(j,k_range)            = X_eta*a_eta_Mz_L(:,j);
            phi_Mz_L(j,k_range)            = X_phi*a_phi_Mz_L(:,j);

            theta_Mz_R(j,k_range)          = X_theta*a_theta_Mz_R(:,j);
            eta_Mz_R(j,k_range)            = X_eta*a_eta_Mz_R(:,j);
            phi_Mz_R(j,k_range)            = X_phi*a_phi_Mz_R(:,j);
            
            if k == 1
                
                theta_Fx_forward_L_max(j)  = max(abs(theta_Fx_forward_L(j,k_range)));
                phi_Fx_forward_L_max(j)    = max(abs(phi_Fx_forward_L(j,k_range)));
                
                theta_Fx_forward_R_max(j)  = max(abs(theta_Fx_forward_R(j,k_range)));
                phi_Fx_forward_R_max(j)    = max(abs(phi_Fx_forward_R(j,k_range)));

                theta_Fy_L_max(j)          = max(abs(theta_Fy_L(j,k_range)));
                phi_Fy_L_max(j)            = max(abs(phi_Fy_L(j,k_range)));

                theta_Fy_R_max(j)          = max(abs(theta_Fy_R(j,k_range)));
                phi_Fy_R_max(j)            = max(abs(phi_Fy_R(j,k_range)));

                theta_Fx_back_L_max(j)     = max(abs(theta_Fx_back_L(j,k_range)));
                phi_Fx_back_L_max(j)       = max(abs(phi_Fx_back_L(j,k_range)));

                theta_Fx_back_R_max(j)     = max(abs(theta_Fx_back_R(j,k_range)));
                phi_Fx_back_R_max(j)       = max(abs(phi_Fx_back_R(j,k_range)));

                theta_Fz_down_L_max(j)     = max(abs(theta_Fz_down_L(j,k_range)));
                phi_Fz_down_L_max(j)       = max(abs(phi_Fz_down_L(j,k_range)));

                theta_Fz_down_R_max(j)     = max(abs(theta_Fz_down_R(j,k_range)));
                phi_Fz_down_R_max(j)       = max(abs(phi_Fz_down_R(j,k_range)));

                theta_Fz_up_L_max(j)       = max(abs(theta_Fz_up_L(j,k_range)));
                phi_Fz_up_L_max(j)         = max(abs(phi_Fz_up_L(j,k_range)));

                theta_Fz_up_R_max(j)       = max(abs(theta_Fz_up_R(j,k_range)));
                phi_Fz_up_R_max(j)         = max(abs(phi_Fz_up_R(j,k_range)));

                theta_Mx_L_max(j)          = max(abs(theta_Mx_L(j,k_range)));
                phi_Mx_L_max(j)            = max(abs(phi_Mx_L(j,k_range)));

                theta_Mx_R_max(j)          = max(abs(theta_Mx_R(j,k_range)));
                phi_Mx_R_max(j)            = max(abs(phi_Mx_R(j,k_range)));

                theta_My_up_L_max(j)       = max(abs(theta_My_up_L(j,k_range)));
                phi_My_up_L_max(j)         = max(abs(phi_My_up_L(j,k_range)));

                theta_My_up_R_max(j)       = max(abs(theta_My_up_R(j,k_range)));
                phi_My_up_R_max(j)         = max(abs(phi_My_up_R(j,k_range)));

                theta_My_down_L_max(j)     = max(abs(theta_My_down_L(j,k_range)));
                phi_My_down_L_max(j)       = max(abs(phi_My_down_L(j,k_range)));

                theta_My_down_R_max(j)     = max(abs(theta_My_down_R(j,k_range)));
                phi_My_down_R_max(j)       = max(abs(phi_My_down_R(j,k_range)));

                theta_Mz_L_max(j)          = max(abs(theta_Mz_L(j,k_range)));
                phi_Mz_L_max(j)            = max(abs(phi_Mz_L(j,k_range)));

                theta_Mz_R_max(j)          = max(abs(theta_Mz_R(j,k_range)));
                phi_Mz_R_max(j)            = max(abs(phi_Mz_R(j,k_range)));
                
            end
            
        end
        
    end
    
    Fx_forward_factor_t   = Fx_forward_factor;
    Fx_back_factor_t      = Fx_back_factor;
    Fy_factor_t           = Fy_factor;
    Fz_down_factor_t      = Fz_down_factor;
    Fz_up_factor_t        = Fz_up_factor;
    Mx_factor_t           = Mx_factor;
    My_up_factor_t        = My_up_factor;
    My_down_factor_t      = My_down_factor;
    Mz_factor_t           = Mz_factor;
    
    for j = 2:nr_steps
        
        % Fx forward
    
        if abs(radtodeg(theta_Fx_forward_L_max(j-1))) < theta_max_range && abs(radtodeg(theta_Fx_forward_L_max(j))) >= theta_max_range
            if Fx_forward_factor_t > (((j-2)/(nr_steps-1))*Fx_forward_factor)
                Fx_forward_factor_t = (((j-2)/(nr_steps-1))*Fx_forward_factor);
            end
        end
        if abs(radtodeg(phi_Fx_forward_L_max(j-1))) < phi_max_range && abs(radtodeg(phi_Fx_forward_L_max(j))) >= phi_max_range
            if Fx_forward_factor_t > (((j-2)/(nr_steps-1))*Fx_forward_factor)
                Fx_forward_factor_t = (((j-2)/(nr_steps-1))*Fx_forward_factor);
            end
        end
        if abs(radtodeg(theta_Fx_forward_R_max(j-1))) < theta_max_range && abs(radtodeg(theta_Fx_forward_R_max(j))) >= theta_max_range
            if Fx_forward_factor_t > (((j-2)/(nr_steps-1))*Fx_forward_factor)
                Fx_forward_factor_t = (((j-2)/(nr_steps-1))*Fx_forward_factor);
            end
        end
        if abs(radtodeg(phi_Fx_forward_R_max(j-1))) < phi_max_range && abs(radtodeg(phi_Fx_forward_R_max(j))) >= phi_max_range
            if Fx_forward_factor_t > (((j-2)/(nr_steps-1))*Fx_forward_factor)
                Fx_forward_factor_t = (((j-2)/(nr_steps-1))*Fx_forward_factor);
            end
        end
        
        % Fx back
        
        if abs(radtodeg(theta_Fx_back_L_max(j-1))) < theta_max_range && abs(radtodeg(theta_Fx_back_L_max(j))) >= theta_max_range
            if Fx_back_factor_t > (((j-2)/(nr_steps-1))*Fx_back_factor)
                Fx_back_factor_t = (((j-2)/(nr_steps-1))*Fx_back_factor);
            end
        end
        if abs(radtodeg(phi_Fx_back_L_max(j-1))) < phi_max_range && abs(radtodeg(phi_Fx_back_L_max(j))) >= phi_max_range
            if Fx_back_factor_t > (((j-2)/(nr_steps-1))*Fx_back_factor)
                Fx_back_factor_t = (((j-2)/(nr_steps-1))*Fx_back_factor);
            end
        end
        if abs(radtodeg(theta_Fx_back_R_max(j-1))) < theta_max_range && abs(radtodeg(theta_Fx_back_R_max(j))) >= theta_max_range
            if Fx_back_factor_t > (((j-2)/(nr_steps-1))*Fx_back_factor)
                Fx_back_factor_t = (((j-2)/(nr_steps-1))*Fx_back_factor);
            end
        end
        if abs(radtodeg(phi_Fx_back_R_max(j-1))) < phi_max_range && abs(radtodeg(phi_Fx_back_R_max(j))) >= phi_max_range
            if Fx_back_factor_t > (((j-2)/(nr_steps-1))*Fx_back_factor)
                Fx_back_factor_t = (((j-2)/(nr_steps-1))*Fx_back_factor);
            end
        end
        
        % Fy
        
        if abs(radtodeg(theta_Fy_L_max(j-1))) < theta_max_range && abs(radtodeg(theta_Fy_L_max(j))) >= theta_max_range
            if Fy_factor_t > (((j-2)/(nr_steps-1))*Fy_factor)
                Fy_factor_t = (((j-2)/(nr_steps-1))*Fy_factor);
            end
        end
        if abs(radtodeg(phi_Fy_L_max(j-1))) < phi_max_range && abs(radtodeg(phi_Fy_L_max(j))) >= phi_max_range
            if Fy_factor_t > (((j-2)/(nr_steps-1))*Fy_factor)
                Fy_factor_t = (((j-2)/(nr_steps-1))*Fy_factor);
            end
        end
        if abs(radtodeg(theta_Fy_R_max(j-1))) < theta_max_range && abs(radtodeg(theta_Fy_R_max(j))) >= theta_max_range
            if Fy_factor_t > (((j-2)/(nr_steps-1))*Fy_factor)
                Fy_factor_t = (((j-2)/(nr_steps-1))*Fy_factor);
            end
        end
        if abs(radtodeg(phi_Fy_R_max(j-1))) < phi_max_range && abs(radtodeg(phi_Fy_R_max(j))) >= phi_max_range
            if Fy_factor_t > (((j-2)/(nr_steps-1))*Fy_factor)
                Fy_factor_t = (((j-2)/(nr_steps-1))*Fy_factor);
            end
        end
        
        % Fz down
        
        if abs(radtodeg(theta_Fz_down_L_max(j-1))) < theta_max_range && abs(radtodeg(theta_Fz_down_L_max(j))) >= theta_max_range
            if Fz_down_factor_t > (((j-2)/(nr_steps-1))*Fz_down_factor)
                Fz_down_factor_t = (((j-2)/(nr_steps-1))*Fz_down_factor);
            end
        end
        if abs(radtodeg(phi_Fz_down_L_max(j-1))) < phi_max_range && abs(radtodeg(phi_Fz_down_L_max(j))) >= phi_max_range
            if Fz_down_factor_t > (((j-2)/(nr_steps-1))*Fz_down_factor)
                Fz_down_factor_t = (((j-2)/(nr_steps-1))*Fz_down_factor);
            end
        end
        if abs(radtodeg(theta_Fz_down_R_max(j-1))) < theta_max_range && abs(radtodeg(theta_Fz_down_R_max(j))) >= theta_max_range
            if Fz_down_factor_t > (((j-2)/(nr_steps-1))*Fz_down_factor)
                Fz_down_factor_t = (((j-2)/(nr_steps-1))*Fz_down_factor);
            end
        end
        if abs(radtodeg(phi_Fz_down_R_max(j-1))) < phi_max_range && abs(radtodeg(phi_Fz_down_R_max(j))) >= phi_max_range
            if Fz_down_factor_t > (((j-2)/(nr_steps-1))*Fz_down_factor)
                Fz_down_factor_t = (((j-2)/(nr_steps-1))*Fz_down_factor);
            end
        end
        
        % Fz up
        
        if abs(radtodeg(theta_Fz_up_L_max(j-1))) < theta_max_range && abs(radtodeg(theta_Fz_up_L_max(j))) >= theta_max_range
            if Fz_up_factor_t > (((j-2)/(nr_steps-1))*Fz_up_factor)
                Fz_up_factor_t = (((j-2)/(nr_steps-1))*Fz_up_factor);
            end
        end
        if abs(radtodeg(phi_Fz_up_L_max(j-1))) < phi_max_range && abs(radtodeg(phi_Fz_up_L_max(j))) >= phi_max_range
            if Fz_up_factor_t > (((j-2)/(nr_steps-1))*Fz_up_factor)
                Fz_up_factor_t = (((j-2)/(nr_steps-1))*Fz_up_factor);
            end
        end
        if abs(radtodeg(theta_Fz_up_R_max(j-1))) < theta_max_range && abs(radtodeg(theta_Fz_up_R_max(j))) >= theta_max_range
            if Fz_up_factor_t > (((j-2)/(nr_steps-1))*Fz_up_factor)
                Fz_up_factor_t = (((j-2)/(nr_steps-1))*Fz_up_factor);
            end
        end
        if abs(radtodeg(phi_Fz_up_R_max(j-1))) < phi_max_range && abs(radtodeg(phi_Fz_up_R_max(j))) >= phi_max_range
            if Fz_up_factor_t > (((j-2)/(nr_steps-1))*Fz_up_factor)
                Fz_up_factor_t = (((j-2)/(nr_steps-1))*Fz_up_factor);
            end
        end
        
        % Mx
        
        if abs(radtodeg(theta_Mx_L_max(j-1))) < theta_max_range && abs(radtodeg(theta_Mx_L_max(j))) >= theta_max_range
            if Mx_factor_t > (((j-2)/(nr_steps-1))*Mx_factor)
                Mx_factor_t = (((j-2)/(nr_steps-1))*Mx_factor);
            end
        end
        if abs(radtodeg(phi_Mx_L_max(j-1))) < phi_max_range && abs(radtodeg(phi_Mx_L_max(j))) >= phi_max_range
            if Mx_factor_t > (((j-2)/(nr_steps-1))*Mx_factor)
                Mx_factor_t = (((j-2)/(nr_steps-1))*Mx_factor);
            end
        end
        if abs(radtodeg(theta_Mx_R_max(j-1))) < theta_max_range && abs(radtodeg(theta_Mx_R_max(j))) >= theta_max_range
            if Mx_factor_t > (((j-2)/(nr_steps-1))*Mx_factor)
                Mx_factor_t = (((j-2)/(nr_steps-1))*Mx_factor);
            end
        end
        if abs(radtodeg(phi_Mx_R_max(j-1))) < phi_max_range && abs(radtodeg(phi_Mx_R_max(j))) >= phi_max_range
            if Mx_factor_t > (((j-2)/(nr_steps-1))*Mx_factor)
                Mx_factor_t = (((j-2)/(nr_steps-1))*Mx_factor);
            end
        end
        
        % My_up
        
        if abs(radtodeg(theta_My_up_L_max(j-1))) < theta_max_range && abs(radtodeg(theta_My_up_L_max(j))) >= theta_max_range
            if My_up_factor_t > (((j-2)/(nr_steps-1))*My_up_factor)
                My_up_factor_t = (((j-2)/(nr_steps-1))*My_up_factor);
            end
        end
        if abs(radtodeg(phi_My_up_L_max(j-1))) < phi_max_range && abs(radtodeg(phi_My_up_L_max(j))) >= phi_max_range
            if My_up_factor_t > (((j-2)/(nr_steps-1))*My_up_factor)
                My_up_factor_t = (((j-2)/(nr_steps-1))*My_up_factor);
            end
        end
        if abs(radtodeg(theta_My_up_R_max(j-1))) < theta_max_range && abs(radtodeg(theta_My_up_R_max(j))) >= theta_max_range
            if My_up_factor_t > (((j-2)/(nr_steps-1))*My_up_factor)
                My_up_factor_t = (((j-2)/(nr_steps-1))*My_up_factor);
            end
        end
        if abs(radtodeg(phi_My_up_R_max(j-1))) < phi_max_range && abs(radtodeg(phi_My_up_R_max(j))) >= phi_max_range
            if My_up_factor_t > (((j-2)/(nr_steps-1))*My_up_factor)
                My_up_factor_t = (((j-2)/(nr_steps-1))*My_up_factor);
            end
        end
        
        
        % My_down
        
        if abs(radtodeg(theta_My_down_L_max(j-1))) < theta_max_range && abs(radtodeg(theta_My_down_L_max(j))) >= theta_max_range
            if My_down_factor_t > (((j-2)/(nr_steps-1))*My_down_factor)
                My_down_factor_t = (((j-2)/(nr_steps-1))*My_down_factor);
            end
        end
        if abs(radtodeg(phi_My_down_L_max(j-1))) < phi_max_range && abs(radtodeg(phi_My_down_L_max(j))) >= phi_max_range
            if My_down_factor_t > (((j-2)/(nr_steps-1))*My_down_factor)
                My_down_factor_t = (((j-2)/(nr_steps-1))*My_down_factor);
            end
        end
        if abs(radtodeg(theta_My_down_R_max(j-1))) < theta_max_range && abs(radtodeg(theta_My_down_R_max(j))) >= theta_max_range
            if My_down_factor_t > (((j-2)/(nr_steps-1))*My_down_factor)
                My_down_factor_t = (((j-2)/(nr_steps-1))*My_down_factor);
            end
        end
        if abs(radtodeg(phi_My_down_R_max(j-1))) < phi_max_range && abs(radtodeg(phi_My_down_R_max(j))) >= phi_max_range
            if My_down_factor_t > (((j-2)/(nr_steps-1))*My_down_factor)
                My_down_factor_t = (((j-2)/(nr_steps-1))*My_down_factor);
            end
        end
        
        % Mz
        
        if abs(radtodeg(theta_Mz_L_max(j-1))) < theta_max_range && abs(radtodeg(theta_Mz_L_max(j))) >= theta_max_range
            if Mz_factor_t > (((j-2)/(nr_steps-1))*Mz_factor)
                Mz_factor_t = (((j-2)/(nr_steps-1))*Mz_factor);
            end
        end
        if abs(radtodeg(phi_Mz_L_max(j-1))) < phi_max_range && abs(radtodeg(phi_Mz_L_max(j))) >= phi_max_range
            if Mz_factor_t > (((j-2)/(nr_steps-1))*Mz_factor)
                Mz_factor_t = (((j-2)/(nr_steps-1))*Mz_factor);
            end
        end
        if abs(radtodeg(theta_Mz_R_max(j-1))) < theta_max_range && abs(radtodeg(theta_Mz_R_max(j))) >= theta_max_range
            if Mz_factor_t > (((j-2)/(nr_steps-1))*Mz_factor)
                Mz_factor_t = (((j-2)/(nr_steps-1))*Mz_factor);
            end
        end
        if abs(radtodeg(phi_Mz_R_max(j-1))) < phi_max_range && abs(radtodeg(phi_Mz_R_max(j))) >= phi_max_range
            if Mz_factor_t > (((j-2)/(nr_steps-1))*Mz_factor)
                Mz_factor_t = (((j-2)/(nr_steps-1))*Mz_factor);
            end
        end
        
    end
    
    Fx_forward_factor   = Fx_forward_factor_t;
    Fx_back_factor      = Fx_back_factor_t;
    Fy_factor           = Fy_factor_t;
    Fz_down_factor      = Fz_down_factor_t;
    Fz_up_factor        = Fz_up_factor_t;
    Mx_factor           = Mx_factor_t;
    My_up_factor        = My_up_factor_t;
    My_down_factor      = My_down_factor_t;
    Mz_factor           = Mz_factor_t;
    
    Fx_forward_factor
    Fx_back_factor
    Fy_factor
    Fz_down_factor
    Fz_up_factor
    Mx_factor
    My_up_factor
    My_down_factor
    Mz_factor
    
    if abs(Fx_forward_max) > 0
        Fx_forward_range    = 0:((Fx_forward_max*Fx_forward_factor)/(nr_steps-1)):(Fx_forward_max*Fx_forward_factor);
    else
        Fx_forward_range    = zeros(1,nr_steps);
    end
    if abs(Fx_back_max) > 0
        Fx_back_range       = 0:((Fx_back_max*Fx_back_factor)/(nr_steps-1)):(Fx_back_max*Fx_back_factor);
    else
        Fx_back_range       = zeros(1,nr_steps);
    end
    if abs(Fy_max) > 0
        Fy_range            = 0:((Fy_max*Fy_factor)/(nr_steps-1)):(Fy_max*Fy_factor);
    else
        Fy_range            = zeros(1,nr_steps);
    end
    if abs(Fz_up_max) > 0
        Fz_up_range         = 0:((Fz_up_max*Fz_up_factor)/(nr_steps-1)):(Fz_up_max*Fz_up_factor);
    else
        Fz_up_range         = zeros(1,nr_steps);
    end
    if abs(Fz_down_max) > 0
        Fz_down_range       = 0:((Fz_down_max*Fz_down_factor)/(nr_steps-1)):(Fz_down_max*Fz_down_factor);
    else
        Fz_down_range       = zeros(1,nr_steps);
    end
    if abs(Mx_max) > 0
        Mx_range            = 0:((Mx_max*Mx_factor)/(nr_steps-1)):(Mx_max*Mx_factor);
    else
        Mx_range            = zeros(1,nr_steps);
    end
    if abs(My_up_max) > 0
        My_up_range         = 0:((My_up_max*My_up_factor)/(nr_steps-1)):(My_up_max*My_up_factor);
    else
        My_up_range         = zeros(1,nr_steps);
    end
    if abs(My_down_max) > 0
        My_down_range       = 0:((My_down_max*My_down_factor)/(nr_steps-1)):(My_down_max*My_down_factor);
    else
        My_down_range       = zeros(1,nr_steps);
    end
    if abs(Mz_max) > 0
        Mz_range            = 0:((Mz_max*Mz_factor)/(nr_steps-1)):(Mz_max*Mz_factor);
    else
        Mz_range            = zeros(1,nr_steps);
    end
    
    Fx_forward_range
    Fx_back_range
    Fy_range
    Fz_up_range
    Fz_down_range
    Mx_range
    My_up_range
    My_down_range
    Mz_range
    
    
    for i = 1:nr_steps
        
        a_theta_Fx_forward_L(:,i)    = a_theta_glob+b_theta_Fx_forward*Fx_forward_range(i);
        a_eta_Fx_forward_L(:,i)      = a_eta_glob+b_eta_Fx_forward*Fx_forward_range(i);
        a_phi_Fx_forward_L(:,i)      = a_phi_glob+b_phi_Fx_forward*Fx_forward_range(i);

        a_theta_Fx_forward_R(:,i)    = a_theta_glob+b_theta_Fx_forward*Fx_forward_range(i);
        a_eta_Fx_forward_R(:,i)      = a_eta_glob+b_eta_Fx_forward*Fx_forward_range(i);
        a_phi_Fx_forward_R(:,i)      = a_phi_glob+b_phi_Fx_forward*Fx_forward_range(i);

        a_theta_Fx_back_L(:,i)       = a_theta_glob+b_theta_Fx_back*Fx_back_range(i);
        a_eta_Fx_back_L(:,i)         = a_eta_glob+b_eta_Fx_back*Fx_back_range(i);
        a_phi_Fx_back_L(:,i)         = a_phi_glob+b_phi_Fx_back*Fx_back_range(i);

        a_theta_Fx_back_R(:,i)       = a_theta_glob+b_theta_Fx_back*Fx_back_range(i);
        a_eta_Fx_back_R(:,i)         = a_eta_glob+b_eta_Fx_back*Fx_back_range(i);
        a_phi_Fx_back_R(:,i)         = a_phi_glob+b_phi_Fx_back*Fx_back_range(i);

        a_theta_Fy_L(:,i)            = a_theta_glob+b_theta_Fy_L*Fy_range(i);
        a_eta_Fy_L(:,i)              = a_eta_glob+b_eta_Fy_L*Fy_range(i);
        a_phi_Fy_L(:,i)              = a_phi_glob+b_phi_Fy_L*Fy_range(i);

        a_theta_Fy_R(:,i)            = a_theta_glob+b_theta_Fy_R*Fy_range(i);
        a_eta_Fy_R(:,i)              = a_eta_glob+b_eta_Fy_R*Fy_range(i);
        a_phi_Fy_R(:,i)              = a_phi_glob+b_phi_Fy_R*Fy_range(i);

        a_theta_Fz_down_L(:,i)       = a_theta_glob+b_theta_Fz_down*Fz_down_range(i);
        a_eta_Fz_down_L(:,i)         = a_eta_glob+b_eta_Fz_down*Fz_down_range(i);
        a_phi_Fz_down_L(:,i)         = a_phi_glob+b_phi_Fz_down*Fz_down_range(i);

        a_theta_Fz_down_R(:,i)       = a_theta_glob+b_theta_Fz_down*Fz_down_range(i);
        a_eta_Fz_down_R(:,i)         = a_eta_glob+b_eta_Fz_down*Fz_down_range(i);
        a_phi_Fz_down_R(:,i)         = a_phi_glob+b_phi_Fz_down*Fz_down_range(i);

        a_theta_Fz_up_L(:,i)         = a_theta_glob+b_theta_Fz_up*Fz_up_range(i);
        a_eta_Fz_up_L(:,i)           = a_eta_glob+b_eta_Fz_up*Fz_up_range(i);
        a_phi_Fz_up_L(:,i)           = a_phi_glob+b_phi_Fz_up*Fz_up_range(i);

        a_theta_Fz_up_R(:,i)         = a_theta_glob+b_theta_Fz_up*Fz_up_range(i);
        a_eta_Fz_up_R(:,i)           = a_eta_glob+b_eta_Fz_up*Fz_up_range(i);
        a_phi_Fz_up_R(:,i)           = a_phi_glob+b_phi_Fz_up*Fz_up_range(i);

        a_theta_Mx_L(:,i)            = a_theta_glob+b_theta_Mx_L*Mx_range(i);
        a_eta_Mx_L(:,i)              = a_eta_glob+b_eta_Mx_L*Mx_range(i);
        a_phi_Mx_L(:,i)              = a_phi_glob+b_phi_Mx_L*Mx_range(i);

        a_theta_Mx_R(:,i)            = a_theta_glob+b_theta_Mx_R*Mx_range(i);
        a_eta_Mx_R(:,i)              = a_eta_glob+b_eta_Mx_R*Mx_range(i);
        a_phi_Mx_R(:,i)              = a_phi_glob+b_phi_Mx_R*Mx_range(i);

        a_theta_My_up_L(:,i)         = a_theta_glob+b_theta_My_up*My_up_range(i);
        a_eta_My_up_L(:,i)           = a_eta_glob+b_eta_My_up*My_up_range(i);
        a_phi_My_up_L(:,i)           = a_phi_glob+b_phi_My_up*My_up_range(i);
        
        a_theta_My_up_R(:,i)         = a_theta_glob+b_theta_My_up*My_up_range(i);
        a_eta_My_up_R(:,i)           = a_eta_glob+b_eta_My_up*My_up_range(i);
        a_phi_My_up_R(:,i)           = a_phi_glob+b_phi_My_up*My_up_range(i);

        a_theta_My_down_L(:,i)       = a_theta_glob+b_theta_My_down*My_down_range(i);
        a_eta_My_down_L(:,i)         = a_eta_glob+b_eta_My_down*My_down_range(i);
        a_phi_My_down_L(:,i)         = a_phi_glob+b_phi_My_down*My_down_range(i);

        a_theta_My_down_R(:,i)       = a_theta_glob+b_theta_My_down*My_down_range(i);
        a_eta_My_down_R(:,i)         = a_eta_glob+b_eta_My_down*My_down_range(i);
        a_phi_My_down_R(:,i)         = a_phi_glob+b_phi_My_down*My_down_range(i);

        a_theta_Mz_L(:,i)            = a_theta_glob+b_theta_Mz_L*Mz_range(i);
        a_eta_Mz_L(:,i)              = a_eta_glob+b_eta_Mz_L*Mz_range(i);
        a_phi_Mz_L(:,i)              = a_phi_glob+b_phi_Mz_L*Mz_range(i);

        a_theta_Mz_R(:,i)            = a_theta_glob+b_theta_Mz_R*Mz_range(i);
        a_eta_Mz_R(:,i)              = a_eta_glob+b_eta_Mz_R*Mz_range(i);
        a_phi_Mz_R(:,i)              = a_phi_glob+b_phi_Mz_R*Mz_range(i);
        
    end
    
    for k = 1:nr_wb_exp

        [ ~, X_theta ]  = Wingbeat_Legendre_matrix( n_pol_theta, down_up_glob, nr_points, 0, 1/f_glob, 0 );
        [ ~, X_eta ]    = Wingbeat_Legendre_matrix( n_pol_eta, down_up_glob, nr_points, 0, 1/f_glob, 0 );
        [ ~, X_phi ]    = Wingbeat_Legendre_matrix( n_pol_phi, down_up_glob, nr_points, 0, 1/f_glob, 0 );

        k_range = ((k-1)*(nr_points-1)+1):((k-1)*(nr_points-1)+nr_points);
            
        for j = 1:nr_steps

            theta_Fx_forward_L(j,k_range)  = X_theta*a_theta_Fx_forward_L(:,j);
            eta_Fx_forward_L(j,k_range)    = X_eta*a_eta_Fx_forward_L(:,j);
            phi_Fx_forward_L(j,k_range)    = X_phi*a_phi_Fx_forward_L(:,j);

            theta_Fx_forward_R(j,k_range)  = X_theta*a_theta_Fx_forward_R(:,j);
            eta_Fx_forward_R(j,k_range)    = X_eta*a_eta_Fx_forward_R(:,j);
            phi_Fx_forward_R(j,k_range)    = X_phi*a_phi_Fx_forward_R(:,j);

            theta_Fx_back_L(j,k_range)     = X_theta*a_theta_Fx_back_L(:,j);
            eta_Fx_back_L(j,k_range)       = X_eta*a_eta_Fx_back_L(:,j);
            phi_Fx_back_L(j,k_range)       = X_phi*a_phi_Fx_back_L(:,j);

            theta_Fx_back_R(j,k_range)     = X_theta*a_theta_Fx_back_R(:,j);
            eta_Fx_back_R(j,k_range)       = X_eta*a_eta_Fx_back_R(:,j);
            phi_Fx_back_R(j,k_range)       = X_phi*a_phi_Fx_back_R(:,j);

            theta_Fy_L(j,k_range)          = X_theta*a_theta_Fy_L(:,j);
            eta_Fy_L(j,k_range)            = X_eta*a_eta_Fy_L(:,j);
            phi_Fy_L(j,k_range)            = X_phi*a_phi_Fy_L(:,j);

            theta_Fy_R(j,k_range)          = X_theta*a_theta_Fy_R(:,j);
            eta_Fy_R(j,k_range)            = X_eta*a_eta_Fy_R(:,j);
            phi_Fy_R(j,k_range)            = X_phi*a_phi_Fy_R(:,j);

            theta_Fz_down_L(j,k_range)     = X_theta*a_theta_Fz_down_L(:,j);
            eta_Fz_down_L(j,k_range)       = X_eta*a_eta_Fz_down_L(:,j);
            phi_Fz_down_L(j,k_range)       = X_phi*a_phi_Fz_down_L(:,j);

            theta_Fz_down_R(j,k_range)     = X_theta*a_theta_Fz_down_R(:,j);
            eta_Fz_down_R(j,k_range)       = X_eta*a_eta_Fz_down_R(:,j);
            phi_Fz_down_R(j,k_range)       = X_phi*a_phi_Fz_down_R(:,j);

            theta_Fz_up_L(j,k_range)       = X_theta*a_theta_Fz_up_L(:,j);
            eta_Fz_up_L(j,k_range)         = X_eta*a_eta_Fz_up_L(:,j);
            phi_Fz_up_L(j,k_range)         = X_phi*a_phi_Fz_up_L(:,j);

            theta_Fz_up_R(j,k_range)       = X_theta*a_theta_Fz_up_R(:,j);
            eta_Fz_up_R(j,k_range)         = X_eta*a_eta_Fz_up_R(:,j);
            phi_Fz_up_R(j,k_range)         = X_phi*a_phi_Fz_up_R(:,j);

            theta_Mx_L(j,k_range)          = X_theta*a_theta_Mx_L(:,j);
            eta_Mx_L(j,k_range)            = X_eta*a_eta_Mx_L(:,j);
            phi_Mx_L(j,k_range)            = X_phi*a_phi_Mx_L(:,j);

            theta_Mx_R(j,k_range)          = X_theta*a_theta_Mx_R(:,j);
            eta_Mx_R(j,k_range)            = X_eta*a_eta_Mx_R(:,j);
            phi_Mx_R(j,k_range)            = X_phi*a_phi_Mx_R(:,j);

            theta_My_up_L(j,k_range)       = X_theta*a_theta_My_up_L(:,j);
            eta_My_up_L(j,k_range)         = X_eta*a_eta_My_up_L(:,j);
            phi_My_up_L(j,k_range)         = X_phi*a_phi_My_up_L(:,j);

            theta_My_up_R(j,k_range)       = X_theta*a_theta_My_up_R(:,j);
            eta_My_up_R(j,k_range)         = X_eta*a_eta_My_up_R(:,j);
            phi_My_up_R(j,k_range)         = X_phi*a_phi_My_up_R(:,j);

            theta_My_down_L(j,k_range)     = X_theta*a_theta_My_down_L(:,j);
            eta_My_down_L(j,k_range)       = X_eta*a_eta_My_down_L(:,j);
            phi_My_down_L(j,k_range)       = X_phi*a_phi_My_down_L(:,j);

            theta_My_down_R(j,k_range)     = X_theta*a_theta_My_down_R(:,j);
            eta_My_down_R(j,k_range)       = X_eta*a_eta_My_down_R(:,j);
            phi_My_down_R(j,k_range)       = X_phi*a_phi_My_down_R(:,j);

            theta_Mz_L(j,k_range)          = X_theta*a_theta_Mz_L(:,j);
            eta_Mz_L(j,k_range)            = X_eta*a_eta_Mz_L(:,j);
            phi_Mz_L(j,k_range)            = X_phi*a_phi_Mz_L(:,j);

            theta_Mz_R(j,k_range)          = X_theta*a_theta_Mz_R(:,j);
            eta_Mz_R(j,k_range)            = X_eta*a_eta_Mz_R(:,j);
            phi_Mz_R(j,k_range)            = X_phi*a_phi_Mz_R(:,j);
            
        end
        
    end

    
%     %----------------------------------------------------------------------
%     
%     figure(1)
%     hFig = figure(1);
%     set(gcf,'PaperPositionMode','auto');
%     set(hFig,'Position',[0 0 600 600]);
%     hold on
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,1); plot(t,radtodeg(theta_Fx_forward_L(j,:)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('theta L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,3); plot(t,radtodeg(eta_Fx_forward_L(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('eta L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,5); plot(t,radtodeg(phi_Fx_forward_L(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('phi L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,2); plot(t,radtodeg(theta_Fx_forward_R(j,:)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('theta R [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,4); plot(t,radtodeg(eta_Fx_forward_R(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('eta R [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,6); plot(t,radtodeg(phi_Fx_forward_R(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('phi R [deg]')
%     hold off
%     [~,h1] = suplabel('Maneuvering wing kinematics Fx forward', 't');
%     set(h1,'FontSize',10)
%     print ([char(settings.plot_loc) '/Maneuvering_wing_kin/FX_forward.eps'] ,'-depsc2');
%     
%     %----------------------------------------------------------------------
%     
%     figure(2)
%     hFig = figure(2);
%     set(gcf,'PaperPositionMode','auto');
%     set(hFig,'Position',[0 0 600 600]);
%     hold on
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,1); plot(t,radtodeg(theta_Fx_back_L(j,:)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('theta L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,3); plot(t,radtodeg(eta_Fx_back_L(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('eta L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,5); plot(t,radtodeg(phi_Fx_back_L(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('phi L [deg]')
%     for j = 1:nr_steps
%         hold on
%     	subplot(3,2,2); plot(t,radtodeg(theta_Fx_back_R(j,:)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('theta R [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,4); plot(t,radtodeg(eta_Fx_back_R(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('eta R [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,6); plot(t,radtodeg(phi_Fx_back_R(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('phi R [deg]')
%     hold off
%     [~,h1] = suplabel('Maneuvering wing kinematics Fx back', 't');
%     set(h1,'FontSize',10)
%     print ([char(settings.plot_loc) '/Maneuvering_wing_kin/FX_back.eps'] ,'-depsc2');
%     
%     
%     %----------------------------------------------------------------------
%     
%     figure(3)
%     hFig = figure(3);
%     set(gcf,'PaperPositionMode','auto');
%     set(hFig,'Position',[0 0 600 600]);
%     hold on
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,1); plot(t,radtodeg(theta_Fy_L(j,:)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('theta L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,3); plot(t,radtodeg(eta_Fy_L(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('eta L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,5); plot(t,radtodeg(phi_Fy_L(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('phi L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,2); plot(t,radtodeg(theta_Fy_R(j,:)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('theta R [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,4); plot(t,radtodeg(eta_Fy_R(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('eta R [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,6); plot(t,radtodeg(phi_Fy_R(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('phi R [deg]')
%     hold off
%     [~,h1] = suplabel('Maneuvering wing kinematics Fy', 't');
%     set(h1,'FontSize',10)
%     print ([char(settings.plot_loc) '/Maneuvering_wing_kin/FY.eps'] ,'-depsc2');
%     
%     
%     %----------------------------------------------------------------------
%     
%     figure(4)
%     hFig = figure(4);
%     set(gcf,'PaperPositionMode','auto');
%     set(hFig,'Position',[0 0 600 600]);
%     hold on
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,1); plot(t,radtodeg(theta_Fz_down_L(j,:)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('theta L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,3); plot(t,radtodeg(eta_Fz_down_L(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('eta L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,5); plot(t,radtodeg(phi_Fz_down_L(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('phi L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,2); plot(t,radtodeg(theta_Fz_down_R(j,:)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('theta R [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,4); plot(t,radtodeg(eta_Fz_down_R(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('eta R [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,6); plot(t,radtodeg(phi_Fz_down_R(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('phi R [deg]')
%     hold off
%     [~,h1] = suplabel('Maneuvering wing kinematics Fz down', 't');
%     set(h1,'FontSize',10)
%     print ([char(settings.plot_loc) '/Maneuvering_wing_kin/FZ_down.eps'] ,'-depsc2');
%     
%     
%     %----------------------------------------------------------------------
%     
%     figure(5)
%     hFig = figure(5);
%     set(gcf,'PaperPositionMode','auto');
%     set(hFig,'Position',[0 0 600 600]);
%     hold on
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,1); plot(t,radtodeg(theta_Fz_up_L(j,:)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('theta L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,3); plot(t,radtodeg(eta_Fz_up_L(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('eta L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,5); plot(t,radtodeg(phi_Fz_up_L(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('phi L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,2); plot(t,radtodeg(theta_Fz_up_R(j,:)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('theta R [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,4); plot(t,radtodeg(eta_Fz_up_R(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('eta R [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,6); plot(t,radtodeg(phi_Fz_up_R(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('phi R [deg]')
%     hold off
%     [~,h1] = suplabel('Maneuvering wing kinematics Fz up', 't');
%     set(h1,'FontSize',10)
%     print ([char(settings.plot_loc) '/Maneuvering_wing_kin/FZ_up.eps'] ,'-depsc2');
%     
%     
%     %----------------------------------------------------------------------
%     
%     figure(6)
%     hFig = figure(6);
%     set(gcf,'PaperPositionMode','auto');
%     set(hFig,'Position',[0 0 600 600]);
%     hold on
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,1); plot(t,radtodeg(theta_Mx_L(j,:)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('theta L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,3); plot(t,radtodeg(eta_Mx_L(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('eta L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,5); plot(t,radtodeg(phi_Mx_L(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('phi L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,2); plot(t,radtodeg(theta_Mx_R(j,:)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('theta R [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,4); plot(t,radtodeg(eta_Mx_R(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('eta R [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,6); plot(t,radtodeg(phi_Mx_R(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('phi R [deg]')
%     hold off
%     [~,h1] = suplabel('Maneuvering wing kinematics Mx', 't');
%     set(h1,'FontSize',10)
%     print ([char(settings.plot_loc) '/Maneuvering_wing_kin/MX.eps'] ,'-depsc2');
%     
%     
%     %----------------------------------------------------------------------
%     
%     figure(7)
%     hFig = figure(7);
%     set(gcf,'PaperPositionMode','auto');
%     set(hFig,'Position',[0 0 600 600]);
%     hold on
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,1); plot(t,radtodeg(theta_My_up_L(j,:)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('theta L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,3); plot(t,radtodeg(eta_My_up_L(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('eta L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,5); plot(t,radtodeg(phi_My_up_L(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('phi L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,2); plot(t,radtodeg(theta_My_up_R(j,:)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('theta R [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,4); plot(t,radtodeg(eta_My_up_R(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('eta R [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,6); plot(t,radtodeg(phi_My_up_R(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('phi R [deg]')
%     hold off
%     [~,h1] = suplabel('Maneuvering wing kinematics My up', 't');
%     set(h1,'FontSize',10)
%     print ([char(settings.plot_loc) '/Maneuvering_wing_kin/MY_up.eps'] ,'-depsc2');
%     
%     
%     %----------------------------------------------------------------------
%     
%     figure(8)
%     hFig = figure(8);
%     set(gcf,'PaperPositionMode','auto');
%     set(hFig,'Position',[0 0 600 600]);
%     hold on
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,1); plot(t,radtodeg(theta_My_down_L(j,:)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('theta L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,3); plot(t,radtodeg(eta_My_down_L(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('eta L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,5); plot(t,radtodeg(phi_My_down_L(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%     	hold off
%     end
%     xlabel('t [s]')
%     ylabel('phi L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,2); plot(t,radtodeg(theta_My_down_R(j,:)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('theta R [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,4); plot(t,radtodeg(eta_My_down_R(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('eta R [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,6); plot(t,radtodeg(phi_My_down_R(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('phi R [deg]')
%     hold off
%     [~,h1] = suplabel('Maneuvering wing kinematics My down', 't');
%     set(h1,'FontSize',10)
%     print ([char(settings.plot_loc) '/Maneuvering_wing_kin/MY_down.eps'] ,'-depsc2');
%     
%     
%     %----------------------------------------------------------------------
%     
%     figure(9)
%     hFig = figure(9);
%     set(gcf,'PaperPositionMode','auto');
%     set(hFig,'Position',[0 0 600 600]);
%     hold on
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,1); plot(t,radtodeg(theta_Mz_L(j,:)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('theta L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,3); plot(t,radtodeg(eta_Mz_L(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('eta L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,5); plot(t,radtodeg(phi_Mz_L(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('phi L [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,2); plot(t,radtodeg(theta_Mz_R(j,:)),'Color',[ (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('theta R [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,4); plot(t,radtodeg(eta_Mz_R(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('eta R [deg]')
%     for j = 1:nr_steps
%         hold on
%         subplot(3,2,6); plot(t,radtodeg(phi_Mz_R(j,:)),'Color',[ (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5-(0.5*(j-1))/(nr_steps-1)) (0.5+(0.5*(j-1))/(nr_steps-1)) ])
%         hold off
%     end
%     xlabel('t [s]')
%     ylabel('phi R [deg]')
%     hold off
%     [~,h1] = suplabel('Maneuvering wing kinematics Mz', 't');
%     set(h1,'FontSize',10)
%     print ([char(settings.plot_loc) '/Maneuvering_wing_kin/MZ.eps'] ,'-depsc2');
    
    
    %----------------------------------------------------------------------

    
    
    % Save results in structure:
    
    Fx_forward = {};
    
    Fx_forward.nr_steps     = nr_steps;
    Fx_forward.nr_wb_exp    = nr_wb_exp;
    Fx_forward.nr_points    = nr_points;
    Fx_forward.factor       = Fx_forward_factor;
    Fx_forward.max          = Fx_forward_max;
    Fx_forward.range        = Fx_forward_range;
    Fx_forward.b_theta      = b_theta_Fx_forward;
    Fx_forward.b_eta        = b_eta_Fx_forward;
    Fx_forward.b_phi        = b_phi_Fx_forward;
    Fx_forward.a_theta_L    = a_theta_Fx_forward_L;
    Fx_forward.a_eta_L      = a_eta_Fx_forward_L;
    Fx_forward.a_phi_L      = a_phi_Fx_forward_L;
    Fx_forward.a_theta_R    = a_theta_Fx_forward_R;
    Fx_forward.a_eta_R      = a_eta_Fx_forward_R;
    Fx_forward.a_phi_R      = a_phi_Fx_forward_R;
    Fx_forward.theta_L      = theta_Fx_forward_L;
    Fx_forward.eta_L        = eta_Fx_forward_L;
    Fx_forward.phi_L        = phi_Fx_forward_L;
    Fx_forward.theta_R      = theta_Fx_forward_R;
    Fx_forward.eta_R        = eta_Fx_forward_R;
    Fx_forward.phi_R        = phi_Fx_forward_R;
    Fx_forward.t            = t;
    Fx_forward.f            = f_glob;
    Fx_forward.down_up      = down_up_glob;
    Fx_forward.m_fly        = m_fly_glob;
    Fx_forward.wing_l       = wing_l_glob;
    
    Fx_back = {};
    
    Fx_back.nr_steps     = nr_steps;
    Fx_back.nr_wb_exp    = nr_wb_exp;
    Fx_back.nr_points    = nr_points;
    Fx_back.factor       = Fx_back_factor;
    Fx_back.max          = Fx_back_max;
    Fx_back.range        = Fx_back_range;
    Fx_back.b_theta      = b_theta_Fx_back;
    Fx_back.b_eta        = b_eta_Fx_back;
    Fx_back.b_phi        = b_phi_Fx_back;
    Fx_back.a_theta_L    = a_theta_Fx_back_L;
    Fx_back.a_eta_L      = a_eta_Fx_back_L;
    Fx_back.a_phi_L      = a_phi_Fx_back_L;
    Fx_back.a_theta_R    = a_theta_Fx_back_R;
    Fx_back.a_eta_R      = a_eta_Fx_back_R;
    Fx_back.a_phi_R      = a_phi_Fx_back_R;
    Fx_back.theta_L      = theta_Fx_back_L;
    Fx_back.eta_L        = eta_Fx_back_L;
    Fx_back.phi_L        = phi_Fx_back_L;
    Fx_back.theta_R      = theta_Fx_back_R;
    Fx_back.eta_R        = eta_Fx_back_R;
    Fx_back.phi_R        = phi_Fx_back_R;
    Fx_back.t            = t;
    Fx_back.f            = f_glob;
    Fx_back.down_up      = down_up_glob;
    Fx_back.m_fly        = m_fly_glob;
    Fx_back.wing_l       = wing_l_glob;
    
    Fy = {};
    
    Fy.nr_steps     = nr_steps;
    Fy.nr_wb_exp    = nr_wb_exp;
    Fy.nr_points    = nr_points;
    Fy.factor       = Fy_factor;
    Fy.max          = Fy_max;
    Fy.range        = Fy_range;
    Fy.b_theta_L    = b_theta_Fy_L;
    Fy.b_eta_L      = b_eta_Fy_L;
    Fy.b_phi_L      = b_phi_Fy_L;
    Fy.b_theta_R    = b_theta_Fy_R;
    Fy.b_eta_R      = b_eta_Fy_R;
    Fy.b_phi_R      = b_phi_Fy_R;
    Fy.a_theta_L    = a_theta_Fy_L;
    Fy.a_eta_L      = a_eta_Fy_L;
    Fy.a_phi_L      = a_phi_Fy_L;
    Fy.a_theta_R    = a_theta_Fy_R;
    Fy.a_eta_R      = a_eta_Fy_R;
    Fy.a_phi_R      = a_phi_Fy_R;
    Fy.theta_L      = theta_Fy_L;
    Fy.eta_L        = eta_Fy_L;
    Fy.phi_L        = phi_Fy_L;
    Fy.theta_R      = theta_Fy_R;
    Fy.eta_R        = eta_Fy_R;
    Fy.phi_R        = phi_Fy_R;
    Fy.t            = t;
    Fy.f            = f_glob;
    Fy.down_up      = down_up_glob;
    Fy.m_fly        = m_fly_glob;
    Fy.wing_l       = wing_l_glob;    
    
    Fz_down = {};
    
    Fz_down.nr_steps     = nr_steps;
    Fz_down.nr_wb_exp    = nr_wb_exp;
    Fz_down.nr_points    = nr_points;
    Fz_down.factor       = Fz_down_factor;
    Fz_down.max          = Fz_down_max;
    Fz_down.range        = Fz_down_range;
    Fz_down.b_theta      = b_theta_Fz_down;
    Fz_down.b_eta        = b_eta_Fz_down;
    Fz_down.b_phi        = b_phi_Fz_down;
    Fz_down.a_theta_L    = a_theta_Fz_down_L;
    Fz_down.a_eta_L      = a_eta_Fz_down_L;
    Fz_down.a_phi_L      = a_phi_Fz_down_L;
    Fz_down.a_theta_R    = a_theta_Fz_down_R;
    Fz_down.a_eta_R      = a_eta_Fz_down_R;
    Fz_down.a_phi_R      = a_phi_Fz_down_R;
    Fz_down.theta_L      = theta_Fz_down_L;
    Fz_down.eta_L        = eta_Fz_down_L;
    Fz_down.phi_L        = phi_Fz_down_L;
    Fz_down.theta_R      = theta_Fz_down_R;
    Fz_down.eta_R        = eta_Fz_down_R;
    Fz_down.phi_R        = phi_Fz_down_R;
    Fz_down.t            = t;
    Fz_down.f            = f_glob;
    Fz_down.down_up      = down_up_glob;
    Fz_down.m_fly        = m_fly_glob;
    Fz_down.wing_l       = wing_l_glob;
    
    Fz_up = {};
    
    Fz_up.nr_steps     = nr_steps;
    Fz_up.nr_wb_exp    = nr_wb_exp;
    Fz_up.nr_points    = nr_points;
    Fz_up.factor       = Fz_up_factor;
    Fz_up.max          = Fz_up_max;
    Fz_up.range        = Fz_up_range;
    Fz_up.b_theta      = b_theta_Fz_up;
    Fz_up.b_eta        = b_eta_Fz_up;
    Fz_up.b_phi        = b_phi_Fz_up;
    Fz_up.a_theta_L    = a_theta_Fz_up_L;
    Fz_up.a_eta_L      = a_eta_Fz_up_L;
    Fz_up.a_phi_L      = a_phi_Fz_up_L;
    Fz_up.a_theta_R    = a_theta_Fz_up_R;
    Fz_up.a_eta_R      = a_eta_Fz_up_R;
    Fz_up.a_phi_R      = a_phi_Fz_up_R;
    Fz_up.theta_L      = theta_Fz_up_L;
    Fz_up.eta_L        = eta_Fz_up_L;
    Fz_up.phi_L        = phi_Fz_up_L;
    Fz_up.theta_R      = theta_Fz_up_R;
    Fz_up.eta_R        = eta_Fz_up_R;
    Fz_up.phi_R        = phi_Fz_up_R;
    Fz_up.t            = t;
    Fz_up.f            = f_glob;
    Fz_up.down_up      = down_up_glob;
    Fz_up.m_fly        = m_fly_glob;
    Fz_up.wing_l       = wing_l_glob;
    
    Mx = {};
    
    Mx.nr_steps     = nr_steps;
    Mx.nr_wb_exp    = nr_wb_exp;
    Mx.nr_points    = nr_points;
    Mx.factor       = Mx_factor;
    Mx.max          = Mx_max;
    Mx.range        = Mx_range;
    Mx.b_theta_L    = b_theta_Mx_L;
    Mx.b_eta_L      = b_eta_Mx_L;
    Mx.b_phi_L      = b_phi_Mx_L;
    Mx.b_theta_R    = b_theta_Mx_R;
    Mx.b_eta_R      = b_eta_Mx_R;
    Mx.b_phi_R      = b_phi_Mx_R;
    Mx.a_theta_L    = a_theta_Mx_L;
    Mx.a_eta_L      = a_eta_Mx_L;
    Mx.a_phi_L      = a_phi_Mx_L;
    Mx.a_theta_R    = a_theta_Mx_R;
    Mx.a_eta_R      = a_eta_Mx_R;
    Mx.a_phi_R      = a_phi_Mx_R;
    Mx.theta_L      = theta_Mx_L;
    Mx.eta_L        = eta_Mx_L;
    Mx.phi_L        = phi_Mx_L;
    Mx.theta_R      = theta_Mx_R;
    Mx.eta_R        = eta_Mx_R;
    Mx.phi_R        = phi_Mx_R;
    Mx.t            = t;
    Mx.f            = f_glob;
    Mx.down_up      = down_up_glob;
    Mx.m_fly        = m_fly_glob;
    Mx.wing_l       = wing_l_glob; 
    
    My_up = {};
    
    My_up.nr_steps     = nr_steps;
    My_up.nr_wb_exp    = nr_wb_exp;
    My_up.nr_points    = nr_points;
    My_up.factor       = My_up_factor;
    My_up.max          = My_up_max;
    My_up.range        = My_up_range;
    My_up.b_theta      = b_theta_My_up;
    My_up.b_eta        = b_eta_My_up;
    My_up.b_phi        = b_phi_My_up;
    My_up.a_theta_L    = a_theta_My_up_L;
    My_up.a_eta_L      = a_eta_My_up_L;
    My_up.a_phi_L      = a_phi_My_up_L;
    My_up.a_theta_R    = a_theta_My_up_R;
    My_up.a_eta_R      = a_eta_My_up_R;
    My_up.a_phi_R      = a_phi_My_up_R;
    My_up.theta_L      = theta_My_up_L;
    My_up.eta_L        = eta_My_up_L;
    My_up.phi_L        = phi_My_up_L;
    My_up.theta_R      = theta_My_up_R;
    My_up.eta_R        = eta_My_up_R;
    My_up.phi_R        = phi_My_up_R;
    My_up.t            = t;
    My_up.f            = f_glob;
    My_up.down_up      = down_up_glob;
    My_up.m_fly        = m_fly_glob;
    My_up.wing_l       = wing_l_glob;
    
    My_down = {};
    
    My_down.nr_steps     = nr_steps;
    My_down.nr_wb_exp    = nr_wb_exp;
    My_down.nr_points    = nr_points;
    My_down.factor       = My_down_factor;
    My_down.max          = My_down_max;
    My_down.range        = My_down_range;
    My_down.b_theta      = b_theta_My_down;
    My_down.b_eta        = b_eta_My_down;
    My_down.b_phi        = b_phi_My_down;
    My_down.a_theta_L    = a_theta_My_down_L;
    My_down.a_eta_L      = a_eta_My_down_L;
    My_down.a_phi_L      = a_phi_My_down_L;
    My_down.a_theta_R    = a_theta_My_down_R;
    My_down.a_eta_R      = a_eta_My_down_R;
    My_down.a_phi_R      = a_phi_My_down_R;
    My_down.theta_L      = theta_My_down_L;
    My_down.eta_L        = eta_My_down_L;
    My_down.phi_L        = phi_My_down_L;
    My_down.theta_R      = theta_My_down_R;
    My_down.eta_R        = eta_My_down_R;
    My_down.phi_R        = phi_My_down_R;
    My_down.t            = t;
    My_down.f            = f_glob;
    My_down.down_up      = down_up_glob;
    My_down.m_fly        = m_fly_glob;
    My_down.wing_l       = wing_l_glob;
    
    Mz = {};
    
    Mz.nr_steps     = nr_steps;
    Mz.nr_wb_exp    = nr_wb_exp;
    Mz.nr_points    = nr_points;
    Mz.factor       = Mz_factor;
    Mz.max          = Mz_max;
    Mz.range        = Mz_range;
    Mz.b_theta_L    = b_theta_Mz_L;
    Mz.b_eta_L      = b_eta_Mz_L;
    Mz.b_phi_L      = b_phi_Mz_L;
    Mz.b_theta_R    = b_theta_Mz_R;
    Mz.b_eta_R      = b_eta_Mz_R;
    Mz.b_phi_R      = b_phi_Mz_R;
    Mz.a_theta_L    = a_theta_Mz_L;
    Mz.a_eta_L      = a_eta_Mz_L;
    Mz.a_phi_L      = a_phi_Mz_L;
    Mz.a_theta_R    = a_theta_Mz_R;
    Mz.a_eta_R      = a_eta_Mz_R;
    Mz.a_phi_R      = a_phi_Mz_R;
    Mz.theta_L      = theta_Mz_L;
    Mz.eta_L        = eta_Mz_L;
    Mz.phi_L        = phi_Mz_L;
    Mz.theta_R      = theta_Mz_R;
    Mz.eta_R        = eta_Mz_R;
    Mz.phi_R        = phi_Mz_R;
    Mz.t            = t;
    Mz.f            = f_glob;
    Mz.down_up      = down_up_glob;
    Mz.m_fly        = m_fly_glob;
    Mz.wing_l       = wing_l_glob; 
    
    save(savefile1,'Fx_forward','Fx_back','Fy','Fz_down','Fz_up','Mx','My_up','My_down','Mz')
    
    
    % Generate .mat file for Robofly experiments:
    
    a_glob.theta    = a_theta_glob;
    a_glob.eta      = a_eta_glob;
    a_glob.phi      = a_phi_glob;
    
    FX_forward_test     = Robofly_tests( Fx_forward, a_glob );
    FX_back_test        = Robofly_tests( Fx_back, a_glob );
    FY_test             = Robofly_tests( Fy, a_glob );
    FZ_up_test          = Robofly_tests( Fz_up, a_glob );
    MX_test             = Robofly_tests( Mx, a_glob );
    MY_up_test          = Robofly_tests( My_up, a_glob );
    MY_down_test        = Robofly_tests( My_down, a_glob );
    MZ_test             = Robofly_tests( Mz, a_glob );
    
    save(savefile2,'FX_forward_test','FX_back_test','FY_test','FZ_up_test','MX_test','MY_up_test','MY_down_test','MZ_test')
    
end

