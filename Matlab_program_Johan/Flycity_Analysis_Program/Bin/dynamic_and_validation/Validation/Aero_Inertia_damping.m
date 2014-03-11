function Aero_Inertia_damping( settings, pathDB )

    close all

    savefile = 'Damp_temp.mat';

    nr_steps = 4;
    
    % Step 1: load and create the wing kinematics

    % Load a_glob and maneuvering wing kinematic functions:
    
    a_avg_theta     = pathDB.poly_fit.a_glob.theta;
    a_avg_eta       = pathDB.poly_fit.a_glob.eta;
    a_avg_phi       = pathDB.poly_fit.a_glob.phi;
    
    f               = pathDB.poly_fit.a_glob.f;
    down_up         = pathDB.poly_fit.a_glob.down_up;
    
    R_strk          = pathDB.rot_mat.Rstr;
    
    n_pol_theta     = settings.n_pol_theta;
    n_pol_eta       = settings.n_pol_eta;
    n_pol_phi       = settings.n_pol_phi;
    
    nr_points       = 1000;
    
        
    % Construct maneuvering wing kinematic polynomial coefficients:
    
    b_ax.forward.theta_L    = pathDB.maneuver.c_fit_ax.b_theta_p(:,1,1);
    b_ax.forward.eta_L      = pathDB.maneuver.c_fit_ax.b_eta_p(:,1,1);
    b_ax.forward.phi_L      = pathDB.maneuver.c_fit_ax.b_phi_p(:,1,1);
    b_ax.forward.theta_R    = pathDB.maneuver.c_fit_ax.b_theta_p(:,1,1);
    b_ax.forward.eta_R      = pathDB.maneuver.c_fit_ax.b_eta_p(:,1,1);
    b_ax.forward.phi_R      = pathDB.maneuver.c_fit_ax.b_phi_p(:,1,1);
    
    b_ax.back.theta_L       = pathDB.maneuver.c_fit_ax.b_theta_n(:,1,1);
    b_ax.back.eta_L         = pathDB.maneuver.c_fit_ax.b_eta_n(:,1,1);
    b_ax.back.phi_L         = pathDB.maneuver.c_fit_ax.b_phi_n(:,1,1);
    b_ax.back.theta_R       = pathDB.maneuver.c_fit_ax.b_theta_n(:,1,1);
    b_ax.back.eta_R         = pathDB.maneuver.c_fit_ax.b_eta_n(:,1,1);
    b_ax.back.phi_R         = pathDB.maneuver.c_fit_ax.b_phi_n(:,1,1);
    
    b_ay.theta_L            = pathDB.maneuver.c_fit_ay.b_theta_L(:,1,2);
    b_ay.eta_L              = pathDB.maneuver.c_fit_ay.b_eta_L(:,1,2);
    b_ay.phi_L              = pathDB.maneuver.c_fit_ay.b_phi_L(:,1,2);
    b_ay.theta_R            = pathDB.maneuver.c_fit_ay.b_theta_R(:,1,2);
    b_ay.eta_R              = pathDB.maneuver.c_fit_ay.b_eta_R(:,1,2);
    b_ay.phi_R              = pathDB.maneuver.c_fit_ay.b_phi_R(:,1,2);
    
    b_az.up.theta_L         = pathDB.maneuver.c_fit_az.b_theta_n(:,1,3);
    b_az.up.eta_L           = pathDB.maneuver.c_fit_az.b_eta_n(:,1,3);
    b_az.up.phi_L           = pathDB.maneuver.c_fit_az.b_phi_n(:,1,3);
    b_az.up.theta_R         = pathDB.maneuver.c_fit_az.b_theta_n(:,1,3);
    b_az.up.eta_R           = pathDB.maneuver.c_fit_az.b_eta_n(:,1,3);
    b_az.up.phi_R           = pathDB.maneuver.c_fit_az.b_phi_n(:,1,3);
    
    b_wx.theta_L            = pathDB.maneuver.c_fit_wx.b_theta_L(:,1,4);
    b_wx.eta_L              = pathDB.maneuver.c_fit_wx.b_eta_L(:,1,4);
    b_wx.phi_L              = pathDB.maneuver.c_fit_wx.b_phi_L(:,1,4);
    b_wx.theta_R            = pathDB.maneuver.c_fit_wx.b_theta_R(:,1,4);
    b_wx.eta_R              = pathDB.maneuver.c_fit_wx.b_eta_R(:,1,4);
    b_wx.phi_R              = pathDB.maneuver.c_fit_wx.b_phi_R(:,1,4);
    
    b_wy.down.theta_L       = pathDB.maneuver.c_fit_wy.b_theta_n(:,1,5);
    b_wy.down.eta_L         = pathDB.maneuver.c_fit_wy.b_eta_n(:,1,5);
    b_wy.down.phi_L         = pathDB.maneuver.c_fit_wy.b_phi_n(:,1,5);
    b_wy.down.theta_R       = pathDB.maneuver.c_fit_wy.b_theta_n(:,1,5);
    b_wy.down.eta_R         = pathDB.maneuver.c_fit_wy.b_eta_n(:,1,5);
    b_wy.down.phi_R         = pathDB.maneuver.c_fit_wy.b_phi_n(:,1,5);
    
    b_wy.up.theta_L         = pathDB.maneuver.c_fit_wy.b_theta_p(:,1,5);
    b_wy.up.eta_L           = pathDB.maneuver.c_fit_wy.b_eta_p(:,1,5);
    b_wy.up.phi_L           = pathDB.maneuver.c_fit_wy.b_phi_p(:,1,5);
    b_wy.up.theta_R         = pathDB.maneuver.c_fit_wy.b_theta_p(:,1,5);
    b_wy.up.eta_R           = pathDB.maneuver.c_fit_wy.b_eta_p(:,1,5);
    b_wy.up.phi_R           = pathDB.maneuver.c_fit_wy.b_phi_p(:,1,5);
    
    b_wz.theta_L            = pathDB.maneuver.c_fit_wz.b_theta_L(:,1,6);
    b_wz.eta_L              = pathDB.maneuver.c_fit_wz.b_eta_L(:,1,6);
    b_wz.phi_L              = pathDB.maneuver.c_fit_wz.b_phi_L(:,1,6);
    b_wz.theta_R            = pathDB.maneuver.c_fit_wz.b_theta_R(:,1,6);
    b_wz.eta_R              = pathDB.maneuver.c_fit_wz.b_eta_R(:,1,6);
    b_wz.phi_R              = pathDB.maneuver.c_fit_wz.b_phi_R(:,1,6);
    
    
    % Non-dimensional forces:
    
    F_ax_forward_max        = pathDB.maneuver.c_fit_ax.FM.FM_1(end);
    F_ax_back_max           = pathDB.maneuver.c_fit_ax.FM.FM_1(1);
    F_ay_max                = pathDB.maneuver.c_fit_ay.FM.FM_2(end);
    F_az_up_max             = pathDB.maneuver.c_fit_az.FM.FM_3(1);
    M_wx_max                = pathDB.maneuver.c_fit_wx.FM.FM_4(end);
    M_wy_down_max           = pathDB.maneuver.c_fit_wy.FM.FM_5(1);
    M_wy_up_max             = pathDB.maneuver.c_fit_wy.FM.FM_5(end);
    M_wz_max                = pathDB.maneuver.c_fit_wz.FM.FM_6(end);
    
        
    F_ax_forward_range      = 0:(F_ax_forward_max/nr_steps):F_ax_forward_max;
    F_ax_back_range         = 0:(F_ax_back_max/nr_steps):F_ax_back_max;
    F_ay_range              = 0:(F_ay_max/nr_steps):F_ay_max;
    F_az_up_range           = 0:(F_az_up_max/nr_steps):F_az_up_max;
    M_wx_range              = 0:(M_wx_max/nr_steps):M_wx_max;
    M_wy_down_range         = 0:(M_wy_down_max/nr_steps):M_wy_down_max;
    M_wy_up_range           = 0:(M_wy_up_max/nr_steps):M_wy_up_max;
    M_wz_range              = 0:(M_wz_max/nr_steps):M_wz_max;
    
    vx_forward_max          = 500;
    vx_back_max             = -180;
    vy_max                  = 550;
    vz_up_max               = -550;
    wx_max                  = 120;
    wy_down_max             = -30;
    wy_up_max               = 90;
    wz_max                  = 70;
    
    vx_forward_range        = 0:(vx_forward_max/nr_steps):vx_forward_max;
    vx_back_range           = 0:(vx_back_max/nr_steps):vx_back_max;
    vy_range                = 0:(vy_max/nr_steps):vy_max;
    vz_up_range             = 0:(vz_up_max/nr_steps):vz_up_max;
    wx_range                = 0:(wx_max/nr_steps):wx_max;
    wy_down_range           = 0:(wy_down_max/nr_steps):wy_down_max;
    wy_up_range             = 0:(wy_up_max/nr_steps):wy_up_max;
    wz_range                = 0:(wz_max/nr_steps):wz_max;
    
    a_x_forward.theta_L     = zeros(n_pol_theta*2+2,nr_steps+1);
    a_x_forward.eta_L       = zeros(n_pol_eta*2+2,nr_steps+1);
    a_x_forward.phi_L       = zeros(n_pol_phi*2+2,nr_steps+1);
    a_x_forward.theta_R     = zeros(n_pol_theta*2+2,nr_steps+1);
    a_x_forward.eta_R       = zeros(n_pol_eta*2+2,nr_steps+1);
    a_x_forward.phi_R       = zeros(n_pol_phi*2+2,nr_steps+1);
    
    a_x_back.theta_L        = zeros(n_pol_theta*2+2,nr_steps+1);
    a_x_back.eta_L          = zeros(n_pol_eta*2+2,nr_steps+1);
    a_x_back.phi_L          = zeros(n_pol_phi*2+2,nr_steps+1);
    a_x_back.theta_R        = zeros(n_pol_theta*2+2,nr_steps+1);
    a_x_back.eta_R          = zeros(n_pol_eta*2+2,nr_steps+1);
    a_x_back.phi_R       	= zeros(n_pol_phi*2+2,nr_steps+1);   
    
    a_y.theta_L             = zeros(n_pol_theta*2+2,nr_steps+1);
    a_y.eta_L               = zeros(n_pol_eta*2+2,nr_steps+1);
    a_y.phi_L               = zeros(n_pol_phi*2+2,nr_steps+1);
    a_y.theta_R             = zeros(n_pol_theta*2+2,nr_steps+1);
    a_y.eta_R               = zeros(n_pol_eta*2+2,nr_steps+1);
    a_y.phi_R               = zeros(n_pol_phi*2+2,nr_steps+1);
    
    a_z_up.theta_L          = zeros(n_pol_theta*2+2,nr_steps+1);
    a_z_up.eta_L            = zeros(n_pol_eta*2+2,nr_steps+1);
    a_z_up.phi_L            = zeros(n_pol_phi*2+2,nr_steps+1);
    a_z_up.theta_R          = zeros(n_pol_theta*2+2,nr_steps+1);
    a_z_up.eta_R            = zeros(n_pol_eta*2+2,nr_steps+1);
    a_z_up.phi_R            = zeros(n_pol_phi*2+2,nr_steps+1); 
    
    w_x.theta_L             = zeros(n_pol_theta*2+2,nr_steps+1);
    w_x.eta_L               = zeros(n_pol_eta*2+2,nr_steps+1);
    w_x.phi_L               = zeros(n_pol_phi*2+2,nr_steps+1);
    w_x.theta_R             = zeros(n_pol_theta*2+2,nr_steps+1);
    w_x.eta_R               = zeros(n_pol_eta*2+2,nr_steps+1);
    w_x.phi_R               = zeros(n_pol_phi*2+2,nr_steps+1);
    
    w_y_down.theta_L        = zeros(n_pol_theta*2+2,nr_steps+1);
    w_y_down.eta_L          = zeros(n_pol_eta*2+2,nr_steps+1);
    w_y_down.phi_L          = zeros(n_pol_phi*2+2,nr_steps+1);
    w_y_down.theta_R        = zeros(n_pol_theta*2+2,nr_steps+1);
    w_y_down.eta_R          = zeros(n_pol_eta*2+2,nr_steps+1);
    w_y_down.phi_R          = zeros(n_pol_phi*2+2,nr_steps+1);
    
    w_y_up.theta_L          = zeros(n_pol_theta*2+2,nr_steps+1);
    w_y_up.eta_L            = zeros(n_pol_eta*2+2,nr_steps+1);
    w_y_up.phi_L            = zeros(n_pol_phi*2+2,nr_steps+1);
    w_y_up.theta_R          = zeros(n_pol_theta*2+2,nr_steps+1);
    w_y_up.eta_R            = zeros(n_pol_eta*2+2,nr_steps+1);
    w_y_up.phi_R            = zeros(n_pol_phi*2+2,nr_steps+1);
    
    w_z.theta_L             = zeros(n_pol_theta*2+2,nr_steps+1);
    w_z.eta_L               = zeros(n_pol_eta*2+2,nr_steps+1);
    w_z.phi_L               = zeros(n_pol_phi*2+2,nr_steps+1);
    w_z.theta_R             = zeros(n_pol_theta*2+2,nr_steps+1);
    w_z.eta_R               = zeros(n_pol_eta*2+2,nr_steps+1);
    w_z.phi_R               = zeros(n_pol_phi*2+2,nr_steps+1);
        
    for i = 1:(nr_steps+1)
        
        a_x_forward.theta_L(:,i)     = a_avg_theta+b_ax.forward.theta_L*F_ax_forward_range(i);
        a_x_forward.eta_L(:,i)       = a_avg_eta+b_ax.forward.eta_L*F_ax_forward_range(i);
        a_x_forward.phi_L(:,i)       = a_avg_phi+b_ax.forward.phi_L*F_ax_forward_range(i);
        a_x_forward.theta_R(:,i)     = a_avg_theta+b_ax.forward.theta_R*F_ax_forward_range(i);
        a_x_forward.eta_R(:,i)       = a_avg_eta+b_ax.forward.eta_R*F_ax_forward_range(i);
        a_x_forward.phi_R(:,i)       = a_avg_phi+b_ax.forward.phi_R*F_ax_forward_range(i);

        a_x_back.theta_L(:,i)        = a_avg_theta+b_ax.back.theta_L*F_ax_back_range(i);
        a_x_back.eta_L(:,i)          = a_avg_eta+b_ax.back.eta_L*F_ax_back_range(i);
        a_x_back.phi_L(:,i)          = a_avg_phi+b_ax.back.phi_L*F_ax_back_range(i);
        a_x_back.theta_R(:,i)        = a_avg_theta+b_ax.back.theta_R*F_ax_back_range(i);
        a_x_back.eta_R(:,i)          = a_avg_eta+b_ax.back.eta_R*F_ax_back_range(i);
        a_x_back.phi_R(:,i)       	 = a_avg_phi+b_ax.back.phi_R*F_ax_back_range(i);

        a_y.theta_L(:,i)             = a_avg_theta+b_ay.theta_L*F_ay_range(i);
        a_y.eta_L(:,i)               = a_avg_eta+b_ay.eta_L*F_ay_range(i);
        a_y.phi_L(:,i)               = a_avg_phi+b_ay.phi_L*F_ay_range(i);
        a_y.theta_R(:,i)             = a_avg_theta+b_ay.theta_R*F_ay_range(i);
        a_y.eta_R(:,i)               = a_avg_eta+b_ay.eta_R*F_ay_range(i);
        a_y.phi_R(:,i)               = a_avg_phi+b_ay.phi_R*F_ay_range(i);

        a_z_up.theta_L(:,i)          = a_avg_theta+b_az.up.theta_L*F_az_up_range(i);
        a_z_up.eta_L(:,i)            = a_avg_eta+b_az.up.eta_L*F_az_up_range(i);
        a_z_up.phi_L(:,i)            = a_avg_phi+b_az.up.phi_L*F_az_up_range(i);
        a_z_up.theta_R(:,i)          = a_avg_theta+b_az.up.theta_R*F_az_up_range(i);
        a_z_up.eta_R(:,i)            = a_avg_eta+b_az.up.eta_R*F_az_up_range(i);
        a_z_up.phi_R(:,i)            = a_avg_phi+b_az.up.phi_R*F_az_up_range(i);

        w_x.theta_L(:,i)             = a_avg_theta+b_wx.theta_L*M_wx_range(i);
        w_x.eta_L(:,i)               = a_avg_eta+b_wx.eta_L*M_wx_range(i);
        w_x.phi_L(:,i)               = a_avg_phi+b_wx.phi_L*M_wx_range(i);
        w_x.theta_R(:,i)             = a_avg_theta+b_wx.theta_R*M_wx_range(i);
        w_x.eta_R(:,i)               = a_avg_eta+b_wx.eta_R*M_wx_range(i);
        w_x.phi_R(:,i)               = a_avg_phi+b_wx.phi_R*M_wx_range(i);

        w_y_down.theta_L(:,i)        = a_avg_theta+b_wy.down.theta_L*M_wy_down_range(i);
        w_y_down.eta_L(:,i)          = a_avg_eta+b_wy.down.eta_L*M_wy_down_range(i);
        w_y_down.phi_L(:,i)          = a_avg_phi+b_wy.down.phi_L*M_wy_down_range(i);
        w_y_down.theta_R(:,i)        = a_avg_theta+b_wy.down.theta_R*M_wy_down_range(i);
        w_y_down.eta_R(:,i)          = a_avg_eta+b_wy.down.eta_R*M_wy_down_range(i);
        w_y_down.phi_R(:,i)          = a_avg_phi+b_wy.down.phi_R*M_wy_down_range(i);

        w_y_up.theta_L(:,i)          = a_avg_theta+b_wy.up.theta_L*M_wy_up_range(i);
        w_y_up.eta_L(:,i)            = a_avg_eta+b_wy.up.eta_L*M_wy_up_range(i);
        w_y_up.phi_L(:,i)            = a_avg_phi+b_wy.up.phi_L*M_wy_up_range(i);
        w_y_up.theta_R(:,i)          = a_avg_theta+b_wy.up.theta_R*M_wy_up_range(i);
        w_y_up.eta_R(:,i)            = a_avg_eta+b_wy.up.eta_R*M_wy_up_range(i);
        w_y_up.phi_R(:,i)            = a_avg_phi+b_wy.up.phi_R*M_wy_up_range(i);

        w_z.theta_L(:,i)             = a_avg_theta+b_wz.theta_L*M_wz_range(i);
        w_z.eta_L(:,i)               = a_avg_eta+b_wz.eta_L*M_wz_range(i);
        w_z.phi_L(:,i)               = a_avg_phi+b_wz.phi_L*M_wz_range(i);
        w_z.theta_R(:,i)             = a_avg_theta+b_wz.theta_R*M_wz_range(i);
        w_z.eta_R(:,i)               = a_avg_eta+b_wz.eta_R*M_wz_range(i);
        w_z.phi_R(:,i)               = a_avg_phi+b_wz.phi_R*M_wz_range(i);
        
    end
    
    
    % Seq_nr 98 is cooresponding the best to the global fly conditions:
    
    seq_nr = 98;
    
    % Step 2: Setup the body_model and wing_model:
    
    body_model.mass_fly         = pathDB.body_model.mass_fly(seq_nr);
    body_model.mass_body        = pathDB.body_model.mass_body(seq_nr);
    body_model.Joint_left       = pathDB.body_model.Joint_left(seq_nr,:)';
    body_model.Joint_right      = pathDB.body_model.Joint_right(seq_nr,:)';
    body_model.cg_b             = pathDB.body_model.cg(seq_nr,:)';
    body_model.Inertia          = pathDB.body_model.Inertia(:,:,seq_nr);
    body_model.x_mod            = pathDB.body_model.x_mod(:,:,seq_nr);
    body_model.y_mod            = pathDB.body_model.y_mod(:,:,seq_nr);
    body_model.z_mod            = pathDB.body_model.z_mod(:,:,seq_nr);
    body_model.g                = 1e-3*settings.g;
    
    wing_model.virtual_mass     = pathDB.wing_model.virtual_mass(seq_nr);
    wing_model.mass             = pathDB.wing_model.mass(seq_nr);
    wing_model.wing_cg_L        = pathDB.wing_model.wing_cg_L(seq_nr,:)';
    wing_model.wing_cg_R        = pathDB.wing_model.wing_cg_R(seq_nr,:)';
    wing_model.virtual_Inertia  = pathDB.wing_model.virtual_Inertia(:,:,seq_nr);
    wing_model.y_sect_L         = pathDB.wing_model.y_sect_L(:,:,seq_nr)';
    wing_model.chords_L         = pathDB.wing_model.chords_L(seq_nr,:)';
    wing_model.y_sect_R         = pathDB.wing_model.y_sect_R(:,:,seq_nr)';
    wing_model.chords_R         = pathDB.wing_model.chords_R(seq_nr,:)';
    wing_model.x_mod_L          = pathDB.wing_model.x_mod_L(:,:,seq_nr);
    wing_model.y_mod_L          = pathDB.wing_model.y_mod_L(:,:,seq_nr);
    wing_model.z_mod_L          = pathDB.wing_model.z_mod_L(:,:,seq_nr);
    wing_model.x_mod_R          = pathDB.wing_model.x_mod_R(:,:,seq_nr);
    wing_model.y_mod_R          = pathDB.wing_model.y_mod_R(:,:,seq_nr);
    wing_model.z_mod_R          = pathDB.wing_model.z_mod_R(:,:,seq_nr);
    wing_model.rho              = settings.rho_air;
    wing_model.nr_sect          = size(wing_model.y_sect_L,2);
    
    a_avg.a_theta_L     = a_x_forward.theta_L(:,1);
    a_avg.a_eta_L       = a_x_forward.eta_L(:,1);
    a_avg.a_phi_L       = a_x_forward.phi_L(:,1);
    a_avg.a_theta_R     = a_x_forward.theta_R(:,1);
    a_avg.a_eta_R       = a_x_forward.eta_R(:,1);
    a_avg.a_phi_R       = a_x_forward.phi_R(:,1);
    a_avg.f             = f;
    a_avg.down_up       = down_up;
    a_avg.nr_points     = nr_points;
    a_avg.R_strk        = R_strk; 
    
    [ F_0, M_0 ] = Initial_cg_and_orientation( body_model, wing_model, a_avg);
    
    
    if exist(char(savefile))==2
        
        temp = load(char(savefile));
        
        F_tot   = temp.F_tot;
        M_tot   = temp.M_tot;
        FA      = temp.FA;
        MA      = temp.MA;
        FI      = temp.FI;
        MI      = temp.MI;
        Fg      = temp.Fg;
        
    else
    
    % Step 3: Simulate the test grid:
    
    F_tot   = zeros(3,nr_steps+1,nr_steps+1,8);
    M_tot   = zeros(3,nr_steps+1,nr_steps+1,8);
    FA      = zeros(3,nr_steps+1,nr_steps+1,8);
    MA      = zeros(3,nr_steps+1,nr_steps+1,8);
    FI      = zeros(3,nr_steps+1,nr_steps+1,8);
    MI      = zeros(3,nr_steps+1,nr_steps+1,8);
    Fg      = zeros(3,nr_steps+1,nr_steps+1,8);
    
    for i = 1:8
        
        for j = 1:(nr_steps+1)
            
            % Create wing kinematics:
            
            if i == 1
        
                a_fit.a_theta_L     = a_x_forward.theta_L(:,j);
                a_fit.a_eta_L       = a_x_forward.eta_L(:,j);
                a_fit.a_phi_L       = a_x_forward.phi_L(:,j);
                a_fit.a_theta_R     = a_x_forward.theta_R(:,j);
                a_fit.a_eta_R       = a_x_forward.eta_R(:,j);
                a_fit.a_phi_R       = a_x_forward.phi_R(:,j);
                
                vb_0_range = [vx_forward_range; zeros(1,nr_steps+1); zeros(1,nr_steps+1)];
                wb_0_range = zeros(3,nr_steps+1);
            
            elseif i == 2
                
                a_fit.a_theta_L     = a_x_back.theta_L(:,j);
                a_fit.a_eta_L       = a_x_back.eta_L(:,j);
                a_fit.a_phi_L       = a_x_back.phi_L(:,j);
                a_fit.a_theta_R     = a_x_back.theta_R(:,j);
                a_fit.a_eta_R       = a_x_back.eta_R(:,j);
                a_fit.a_phi_R       = a_x_back.phi_R(:,j);
                
                vb_0_range = [vx_back_range; zeros(1,nr_steps+1); zeros(1,nr_steps+1)];
                wb_0_range = zeros(3,nr_steps+1);
                
            elseif i == 3
                
                a_fit.a_theta_L     = a_y.theta_L(:,j);
                a_fit.a_eta_L       = a_y.eta_L(:,j);
                a_fit.a_phi_L       = a_y.phi_L(:,j);
                a_fit.a_theta_R     = a_y.theta_R(:,j);
                a_fit.a_eta_R       = a_y.eta_R(:,j);
                a_fit.a_phi_R       = a_y.phi_R(:,j);
                
                vb_0_range = [ zeros(1,nr_steps+1); vy_range; zeros(1,nr_steps+1)];
                wb_0_range = zeros(3,nr_steps+1);
                
            elseif i == 4
                
                a_fit.a_theta_L     = a_z_up.theta_L(:,j);
                a_fit.a_eta_L       = a_z_up.eta_L(:,j);
                a_fit.a_phi_L       = a_z_up.phi_L(:,j);
                a_fit.a_theta_R     = a_z_up.theta_R(:,j);
                a_fit.a_eta_R       = a_z_up.eta_R(:,j);
                a_fit.a_phi_R       = a_z_up.phi_R(:,j);
                
                vb_0_range = [ zeros(1,nr_steps+1); zeros(1,nr_steps+1); vz_up_range];
                wb_0_range = zeros(3,nr_steps+1);
                
            elseif i == 5
                
                a_fit.a_theta_L     = w_x.theta_L(:,j);
                a_fit.a_eta_L       = w_x.eta_L(:,j);
                a_fit.a_phi_L       = w_x.phi_L(:,j);
                a_fit.a_theta_R     = w_x.theta_R(:,j);
                a_fit.a_eta_R       = w_x.eta_R(:,j);
                a_fit.a_phi_R       = w_x.phi_R(:,j);
                
                vb_0_range = zeros(3,nr_steps+1);
                wb_0_range = [ wx_range; zeros(1,nr_steps+1); zeros(1,nr_steps+1)];
                
            elseif i == 6
                
                a_fit.a_theta_L     = w_y_down.theta_L(:,j);
                a_fit.a_eta_L       = w_y_down.eta_L(:,j);
                a_fit.a_phi_L       = w_y_down.phi_L(:,j);
                a_fit.a_theta_R     = w_y_down.theta_R(:,j);
                a_fit.a_eta_R       = w_y_down.eta_R(:,j);
                a_fit.a_phi_R       = w_y_down.phi_R(:,j);
                
                vb_0_range = zeros(3,nr_steps+1);
                wb_0_range = [ zeros(1,nr_steps+1); wy_down_range; zeros(1,nr_steps+1)];
                
            elseif i == 7
                
                a_fit.a_theta_L     = w_y_up.theta_L(:,j);
                a_fit.a_eta_L       = w_y_up.eta_L(:,j);
                a_fit.a_phi_L       = w_y_up.phi_L(:,j);
                a_fit.a_theta_R     = w_y_up.theta_R(:,j);
                a_fit.a_eta_R       = w_y_up.eta_R(:,j);
                a_fit.a_phi_R       = w_y_up.phi_R(:,j);
                
                vb_0_range = zeros(3,nr_steps+1);
                wb_0_range = [ zeros(1,nr_steps+1); wy_up_range; zeros(1,nr_steps+1)];
                
            elseif i == 8
                
                a_fit.a_theta_L     = w_z.theta_L(:,j);
                a_fit.a_eta_L       = w_z.eta_L(:,j);
                a_fit.a_phi_L       = w_z.phi_L(:,j);
                a_fit.a_theta_R     = w_z.theta_R(:,j);
                a_fit.a_eta_R       = w_z.eta_R(:,j);
                a_fit.a_phi_R       = w_z.phi_R(:,j);
                
                vb_0_range = zeros(3,nr_steps+1);
                wb_0_range = [ zeros(1,nr_steps+1); zeros(1,nr_steps+1); wz_range];
                
            end
            
            a_fit.f             = f;
            a_fit.down_up       = down_up;
            a_fit.nr_points     = nr_points;
            a_fit.R_strk        = R_strk;   

            [ kine ] = angular_velocities_polynomial( a_fit );

            t_func              = kine.t;
            kine.R_strk         = R_strk;
            kine.down_up_t      = t_func(end)*down_up;  
            
            for k = 1:(nr_steps+1)
                        
                body_IC.xyz_0   = [0; 0; 0];
                body_IC.qb_0    = quat2mat((R_strk'*[1 0 0; 0 -1 0; 0 0 -1]))';
                body_IC.vb_0    = vb_0_range(:,k);
                body_IC.wb_0    = wb_0_range(:,k);
                body_IC.F_0     = F_0;
                body_IC.M_0     = M_0;

                clear sim_data

                % Tolerance settings ode45:

                options = odeset('RelTol',1e-5,'AbsTol',[1e-15 1e-15 1e-15 1e-15 1e-15 1e-15]);

                % case_nr == 0: Inertia only
                % case_nr == 1: Inertia, translational aerodynamics, gravity
                % case_nr == 2: Inertia, translational and rotational aerodynamics,
                % gravity

                case_nr = 2;

                tic

                sim_data = dyn_sim(t_func,kine,body_model,wing_model,body_IC,options,case_nr);

                toc

%                 n = 1;
% 
%                 [n] = plot_dyn_sim( sim_data, case_nr, n+1 );
% 
%                 clear body_IC case_nr options
% 
%                 sim_data
                
                FI_acc_mean  = sim_data.FI_acc_strk_mean;
                MI_acc_mean  = sim_data.MI_acc_strk_mean;
                FA_mean      = sim_data.FA_strk_mean;
                MA_mean      = sim_data.MA_strk_mean;
                FI_vel_mean  = sim_data.FI_vel_strk_mean;
                MI_vel_mean  = sim_data.MI_vel_strk_mean;
                Fg_mean      = sim_data.Fg_strk_mean;
                
                F_tot(:,j,k,i)    = FI_acc_mean;
                M_tot(:,j,k,i)    = MI_acc_mean;
                FA(:,j,k,i)       = FA_mean;
                MA(:,j,k,i)       = MA_mean;
                FI(:,j,k,i)       = FI_vel_mean;
                MI(:,j,k,i)       = MI_vel_mean;
                Fg(:,j,k,i)       = Fg_mean;
            
            end
            
        end
        
    end
    
    save(savefile,'F_tot','M_tot','FA','MA','FI','MI','Fg')
    
    end
    
    Fx_forward_tot      = zeros(3,nr_steps+1,nr_steps+1);
    Mx_forward_tot      = zeros(3,nr_steps+1,nr_steps+1);
    
    Fx_back_tot         = zeros(3,nr_steps+1,nr_steps+1);
    Mx_back_tot         = zeros(3,nr_steps+1,nr_steps+1);
    
    Fy_tot              = zeros(3,nr_steps+1,nr_steps+1);
    My_tot              = zeros(3,nr_steps+1,nr_steps+1);
    
    Fz_up_tot           = zeros(3,nr_steps+1,nr_steps+1);
    Mz_up_tot           = zeros(3,nr_steps+1,nr_steps+1);
    
    Fx_tot              = zeros(3,nr_steps+1,nr_steps+1);
    Mx_tot              = zeros(3,nr_steps+1,nr_steps+1);
    
    Fy_down_tot         = zeros(3,nr_steps+1,nr_steps+1);
    My_down_tot         = zeros(3,nr_steps+1,nr_steps+1);
    
    Fy_up_tot           = zeros(3,nr_steps+1,nr_steps+1);
    My_up_tot           = zeros(3,nr_steps+1,nr_steps+1);
    
    Fz_tot              = zeros(3,nr_steps+1,nr_steps+1);
    Mz_tot              = zeros(3,nr_steps+1,nr_steps+1);

    Fx_forward_A        = zeros(3,nr_steps+1,nr_steps+1);
    Mx_forward_A        = zeros(3,nr_steps+1,nr_steps+1);
    
    Fx_back_A         = zeros(3,nr_steps+1,nr_steps+1);
    Mx_back_A         = zeros(3,nr_steps+1,nr_steps+1);
    
    Fy_A              = zeros(3,nr_steps+1,nr_steps+1);
    My_A              = zeros(3,nr_steps+1,nr_steps+1);
    
    Fz_up_A           = zeros(3,nr_steps+1,nr_steps+1);
    Mz_up_A           = zeros(3,nr_steps+1,nr_steps+1);
    
    Fx_A              = zeros(3,nr_steps+1,nr_steps+1);
    Mx_A              = zeros(3,nr_steps+1,nr_steps+1);
    
    Fy_down_A         = zeros(3,nr_steps+1,nr_steps+1);
    My_down_A         = zeros(3,nr_steps+1,nr_steps+1);
    
    Fy_up_A           = zeros(3,nr_steps+1,nr_steps+1);
    My_up_A           = zeros(3,nr_steps+1,nr_steps+1);
    
    Fz_A              = zeros(3,nr_steps+1,nr_steps+1);
    Mz_A              = zeros(3,nr_steps+1,nr_steps+1);
    
    Fx_forward_I        = zeros(3,nr_steps+1,nr_steps+1);
    Mx_forward_I        = zeros(3,nr_steps+1,nr_steps+1);
    
    Fx_back_I         = zeros(3,nr_steps+1,nr_steps+1);
    Mx_back_I         = zeros(3,nr_steps+1,nr_steps+1);
    
    Fy_I              = zeros(3,nr_steps+1,nr_steps+1);
    My_I              = zeros(3,nr_steps+1,nr_steps+1);
    
    Fz_up_I           = zeros(3,nr_steps+1,nr_steps+1);
    Mz_up_I           = zeros(3,nr_steps+1,nr_steps+1);
    
    Fx_I              = zeros(3,nr_steps+1,nr_steps+1);
    Mx_I              = zeros(3,nr_steps+1,nr_steps+1);
    
    Fy_down_I         = zeros(3,nr_steps+1,nr_steps+1);
    My_down_I         = zeros(3,nr_steps+1,nr_steps+1);
    
    Fy_up_I           = zeros(3,nr_steps+1,nr_steps+1);
    My_up_I           = zeros(3,nr_steps+1,nr_steps+1);
    
    Fz_I              = zeros(3,nr_steps+1,nr_steps+1);
    Mz_I              = zeros(3,nr_steps+1,nr_steps+1);
    
    for i = 1:(nr_steps+1)
        
        for j = 1:(nr_steps+1)
            
%             Fx_forward_tot(:,j,i) = F_tot(:,i,j,1);
%             Mx_forward_tot(:,j,i) = M_tot(:,i,j,1);
% 
%             Fx_back_tot(:,j,i)    = F_tot(:,i,j,2);
%             Mx_back_tot(:,j,i)    = M_tot(:,i,j,2);
% 
%             Fy_tot(:,j,i)         = F_tot(:,i,j,3);
%             My_tot(:,j,i)         = M_tot(:,i,j,3);
% 
%             Fz_up_tot(:,j,i)      = F_tot(:,i,j,4);
%             Mz_up_tot(:,j,i)      = M_tot(:,i,j,4);
% 
%             Fx_tot(:,j,i)         = F_tot(:,i,j,5);
%             Mx_tot(:,j,i)         = M_tot(:,i,j,5);
% 
%             Fy_down_tot(:,j,i)    = F_tot(:,i,j,6);
%             My_down_tot(:,j,i)    = M_tot(:,i,j,6);
% 
%             Fy_up_tot(:,j,i)      = F_tot(:,i,j,7);
%             My_up_tot(:,j,i)      = M_tot(:,i,j,7);
% 
%             Fz_tot(:,j,i)         = F_tot(:,i,j,8);
%             Mz_tot(:,j,i)         = M_tot(:,i,j,8);

            Fx_forward_tot(:,j,i) = F_tot(:,i,j,1)-FI(:,i,j,1);
            Mx_forward_tot(:,j,i) = M_tot(:,i,j,1)-MI(:,i,j,1);

            Fx_back_tot(:,j,i)    = F_tot(:,i,j,2)-FI(:,i,j,2);
            Mx_back_tot(:,j,i)    = M_tot(:,i,j,2)-MI(:,i,j,2);

            Fy_tot(:,j,i)         = F_tot(:,i,j,3)-FI(:,i,j,3);
            My_tot(:,j,i)         = M_tot(:,i,j,3)-MI(:,i,j,3);

            Fz_up_tot(:,j,i)      = F_tot(:,i,j,4)-FI(:,i,j,4);
            Mz_up_tot(:,j,i)      = M_tot(:,i,j,4)-MI(:,i,j,4);

            Fx_tot(:,j,i)         = F_tot(:,i,j,5)-FI(:,i,j,5);
            Mx_tot(:,j,i)         = M_tot(:,i,j,5)-MI(:,i,j,5);

            Fy_down_tot(:,j,i)    = F_tot(:,i,j,6)-FI(:,i,j,6);
            My_down_tot(:,j,i)    = M_tot(:,i,j,6)-MI(:,i,j,6);

            Fy_up_tot(:,j,i)      = F_tot(:,i,j,7)-FI(:,i,j,7);
            My_up_tot(:,j,i)      = M_tot(:,i,j,7)-MI(:,i,j,7);

            Fz_tot(:,j,i)         = F_tot(:,i,j,8)-FI(:,i,j,8);
            Mz_tot(:,j,i)         = M_tot(:,i,j,8)-MI(:,i,j,8);
            
            Fx_forward_A(:,j,i)   = FA(:,i,j,1);
            Mx_forward_A(:,j,i)   = MA(:,i,j,1);

            Fx_back_A(:,j,i)      = FA(:,i,j,2);
            Mx_back_A(:,j,i)      = MA(:,i,j,2);

            Fy_A(:,j,i)           = FA(:,i,j,3);
            My_A(:,j,i)           = MA(:,i,j,3);

            Fz_up_A(:,j,i)        = FA(:,i,j,4);
            Mz_up_A(:,j,i)        = MA(:,i,j,4);

            Fx_A(:,j,i)           = FA(:,i,j,5);
            Mx_A(:,j,i)           = MA(:,i,j,5);

            Fy_down_A(:,j,i)      = FA(:,i,j,6);
            My_down_A(:,j,i)      = MA(:,i,j,6);

            Fy_up_A(:,j,i)        = FA(:,i,j,7);
            My_up_A(:,j,i)        = MA(:,i,j,7);

            Fz_A(:,j,i)           = FA(:,i,j,8);
            Mz_A(:,j,i)           = MA(:,i,j,8);
            
%             Fx_forward_I(:,j,i)   = FI(:,i,j,1);
%             Mx_forward_I(:,j,i)   = MI(:,i,j,1);
% 
%             Fx_back_I(:,j,i)      = FI(:,i,j,2);
%             Mx_back_I(:,j,i)      = MI(:,i,j,2);
% 
%             Fy_I(:,j,i)           = FI(:,i,j,3);
%             My_I(:,j,i)           = MI(:,i,j,3);
% 
%             Fz_up_I(:,j,i)        = FI(:,i,j,4);
%             Mz_up_I(:,j,i)        = MI(:,i,j,4);
% 
%             Fx_I(:,j,i)           = FI(:,i,j,5);
%             Mx_I(:,j,i)           = MI(:,i,j,5);
% 
%             Fy_down_I(:,j,i)      = FI(:,i,j,6);
%             My_down_I(:,j,i)      = MI(:,i,j,6);
% 
%             Fy_up_I(:,j,i)        = FI(:,i,j,7);
%             My_up_I(:,j,i)        = MI(:,i,j,7);
% 
%             Fz_I(:,j,i)           = FI(:,i,j,8);
%             Mz_I(:,j,i)           = MI(:,i,j,8);

            Fx_forward_I(:,j,i)   = -FI(:,i,j,1);
            Mx_forward_I(:,j,i)   = -MI(:,i,j,1);

            Fx_back_I(:,j,i)      = -FI(:,i,j,2);
            Mx_back_I(:,j,i)      = -MI(:,i,j,2);

            Fy_I(:,j,i)           = -FI(:,i,j,3);
            My_I(:,j,i)           = -MI(:,i,j,3);

            Fz_up_I(:,j,i)        = -FI(:,i,j,4);
            Mz_up_I(:,j,i)        = -MI(:,i,j,4);

            Fx_I(:,j,i)           = -FI(:,i,j,5);
            Mx_I(:,j,i)           = -MI(:,i,j,5);

            Fy_down_I(:,j,i)      = -FI(:,i,j,6);
            My_down_I(:,j,i)      = -MI(:,i,j,6);

            Fy_up_I(:,j,i)        = -FI(:,i,j,7);
            My_up_I(:,j,i)        = -MI(:,i,j,7);

            Fz_I(:,j,i)           = -FI(:,i,j,8);
            Mz_I(:,j,i)           = -MI(:,i,j,8);
            
        end
        
    end
    
    % Plot total damping:
       
    figure(1)
    hFig = figure(1);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,1); plot(vx_forward_range,Fx_forward_tot(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('F_x [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,3); plot(vx_forward_range,Fx_forward_tot(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('F_y [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,5); plot(vx_forward_range,Fx_forward_tot(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('F_z [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,2); plot(vx_forward_range,Mx_forward_tot(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('M_x [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,4); plot(vx_forward_range,Mx_forward_tot(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('M_y [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,6); plot(vx_forward_range,Mx_forward_tot(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('M_z [N*mm]')
    hold off
    [~,h1] = suplabel('Total damping v_x forward', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Damping/Total_vx_forward.eps'] ,'-depsc2');
    
    figure(2)
    hFig = figure(2);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,1); plot(vx_back_range,Fx_back_tot(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('F_x [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,3); plot(vx_back_range,Fx_back_tot(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('F_y [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,5); plot(vx_back_range,Fx_back_tot(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('F_z [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,2); plot(vx_back_range,Mx_back_tot(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('M_x [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,4); plot(vx_back_range,Mx_back_tot(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('M_y [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,6); plot(vx_back_range,Mx_back_tot(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('M_z [N*mm]')
    hold off
    [~,h1] = suplabel('Total damping v_x back', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Damping/Total_vx_back.eps'] ,'-depsc2');

    figure(3)
    hFig = figure(3);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,1); plot(vy_range,Fy_tot(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_y [mm/s]')
    ylabel('F_x [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,3); plot(vy_range,Fy_tot(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_y [mm/s]')
    ylabel('F_y [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,5); plot(vy_range,Fy_tot(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_y [mm/s]')
    ylabel('F_z [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,2); plot(vy_range,My_tot(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_y [mm/s]')
    ylabel('M_x [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,4); plot(vy_range,My_tot(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_y [mm/s]')
    ylabel('M_y [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,6); plot(vy_range,My_tot(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_y [mm/s]')
    ylabel('M_z [N*mm]')
    hold off
    [~,h1] = suplabel('Total damping v_y', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Damping/Total_vy.eps'] ,'-depsc2');
    
    figure(4)
    hFig = figure(4);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,1); plot(vz_up_range,Fz_up_tot(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_z [mm/s]')
    ylabel('F_x [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,3); plot(vz_up_range,Fz_up_tot(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_z [mm/s]')
    ylabel('F_y [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,5); plot(vz_up_range,Fz_up_tot(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_z [mm/s]')
    ylabel('F_z [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,2); plot(vz_up_range,Mz_up_tot(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_z [mm/s]')
    ylabel('M_x [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,4); plot(vz_up_range,Mz_up_tot(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_z [mm/s]')
    ylabel('M_y [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,6); plot(vz_up_range,Mz_up_tot(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_z [mm/s]')
    ylabel('M_z [N*mm]')
    hold off
    [~,h1] = suplabel('Total damping v_z up', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Damping/Total_vz_up.eps'] ,'-depsc2');

    
    figure(5)
    hFig = figure(5);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,1); plot(wx_range,Fx_tot(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_x [rad/s]')
    ylabel('F_x [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,3); plot(wx_range,Fx_tot(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_x [rad/s]')
    ylabel('F_y [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,5); plot(wx_range,Fx_tot(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_x [rad/s]')
    ylabel('F_z [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,2); plot(wx_range,Mx_tot(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_x [rad/s]')
    ylabel('M_x [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,4); plot(wx_range,Mx_tot(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_x [rad/s]')
    ylabel('M_y [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,6); plot(wx_range,Mx_tot(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_x [rad/s]')
    ylabel('M_z [N*mm]')
    hold off
    [~,h1] = suplabel('Total damping \omega_x', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Damping/Total_wx.eps'] ,'-depsc2');
    

    figure(6)
    hFig = figure(6);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,1); plot(wy_down_range,Fy_down_tot(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('F_x [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,3); plot(wy_down_range,Fy_down_tot(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('F_y [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,5); plot(wy_down_range,Fy_down_tot(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('F_z [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,2); plot(wy_down_range,My_down_tot(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('M_x [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,4); plot(wy_down_range,My_down_tot(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('M_y [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,6); plot(wy_down_range,My_down_tot(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('M_z [N*mm]')
    hold off
    [~,h1] = suplabel('Total damping \omega_y down', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Damping/Total_wy_down.eps'] ,'-depsc2');
    
    
    
    figure(7)
    hFig = figure(7);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,1); plot(wy_up_range,Fy_up_tot(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('F_x [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,3); plot(wy_up_range,Fy_up_tot(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('F_y [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,5); plot(wy_up_range,Fy_up_tot(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('F_z [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,2); plot(wy_up_range,My_up_tot(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('M_x [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,4); plot(wy_up_range,My_up_tot(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('M_y [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,6); plot(wy_up_range,My_up_tot(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('M_z [N*mm]')
    hold off
    [~,h1] = suplabel('Total damping \omega_y up', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Damping/Total_wy_up.eps'] ,'-depsc2');
    
    
    figure(8)
    hFig = figure(8);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,1); plot(wz_range,Fz_tot(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_z [rad/s]')
    ylabel('F_x [N]')
    for j = 1:nr_steps
        hold on
        subplot(3,2,3); plot(wz_range,Fz_tot(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_z [rad/s]')
    ylabel('F_y [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,5); plot(wz_range,Fz_tot(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_z [rad/s]')
    ylabel('F_z [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,2); plot(wz_range,Mz_tot(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_z [rad/s]')
    ylabel('M_x [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,4); plot(wz_range,Mz_tot(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_z [rad/s]')
    ylabel('M_y [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,6); plot(wz_range,Mz_tot(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_x [rad/s]')
    ylabel('M_z [N*mm]')
    hold off
    [~,h1] = suplabel('Total damping \omega_z', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Damping/Total_wz.eps'] ,'-depsc2');    

    
    % Plot aerodynamic damping:
    
    figure(9)
    hFig = figure(9);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,1); plot(vx_forward_range,Fx_forward_A(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('F_x [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,3); plot(vx_forward_range,Fx_forward_A(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('F_y [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,5); plot(vx_forward_range,Fx_forward_A(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('F_z [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,2); plot(vx_forward_range,Mx_forward_A(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('M_x [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,4); plot(vx_forward_range,Mx_forward_A(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('M_y [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,6); plot(vx_forward_range,Mx_forward_A(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('M_z [N*mm]')
    hold off
    [~,h1] = suplabel('Aerodynamic damping v_x forward', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Damping/Aero_vx_forward.eps'] ,'-depsc2');
    
    figure(10)
    hFig = figure(10);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,1); plot(vx_back_range,Fx_back_A(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('F_x [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,3); plot(vx_back_range,Fx_back_A(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('F_y [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,5); plot(vx_back_range,Fx_back_A(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('F_z [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,2); plot(vx_back_range,Mx_back_A(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('M_x [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,4); plot(vx_back_range,Mx_back_A(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('M_y [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,6); plot(vx_back_range,Mx_back_A(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('M_z [N*mm]')
    hold off
    [~,h1] = suplabel('Aerodynamic damping v_x back', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Damping/Aero_vx_back.eps'] ,'-depsc2');

    figure(11)
    hFig = figure(11);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,1); plot(vy_range,Fy_A(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_y [mm/s]')
    ylabel('F_x [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,3); plot(vy_range,Fy_A(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_y [mm/s]')
    ylabel('F_y [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,5); plot(vy_range,Fy_A(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_y [mm/s]')
    ylabel('F_z [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,2); plot(vy_range,My_A(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_y [mm/s]')
    ylabel('M_x [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,4); plot(vy_range,My_A(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_y [mm/s]')
    ylabel('M_y [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,6); plot(vy_range,My_A(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_y [mm/s]')
    ylabel('M_z [N*mm]')
    hold off
    [~,h1] = suplabel('Aerodynamic damping v_y', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Damping/Aero_vy.eps'] ,'-depsc2');
    
    figure(12)
    hFig = figure(12);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,1); plot(vz_up_range,Fz_up_A(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_z [mm/s]')
    ylabel('F_x [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,3); plot(vz_up_range,Fz_up_A(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_z [mm/s]')
    ylabel('F_y [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,5); plot(vz_up_range,Fz_up_A(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_z [mm/s]')
    ylabel('F_z [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,2); plot(vz_up_range,Mz_up_A(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_z [mm/s]')
    ylabel('M_x [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,4); plot(vz_up_range,Mz_up_A(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_z [mm/s]')
    ylabel('M_y [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,6); plot(vz_up_range,Mz_up_A(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_z [mm/s]')
    ylabel('M_z [N*mm]')
    hold off
    [~,h1] = suplabel('Aerodynamic damping v_z up', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Damping/Aero_vz_up.eps'] ,'-depsc2');

    
    figure(13)
    hFig = figure(13);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,1); plot(wx_range,Fx_A(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_x [rad/s]')
    ylabel('F_x [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,3); plot(wx_range,Fx_A(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_x [rad/s]')
    ylabel('F_y [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,5); plot(wx_range,Fx_A(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_x [rad/s]')
    ylabel('F_z [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,2); plot(wx_range,Mx_A(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_x [rad/s]')
    ylabel('M_x [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,4); plot(wx_range,Mx_A(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_x [rad/s]')
    ylabel('M_y [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,6); plot(wx_range,Mx_A(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_x [rad/s]')
    ylabel('M_z [N*mm]')
    hold off
    [~,h1] = suplabel('Aerodynamic damping \omega_x', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Damping/Aero_wx.eps'] ,'-depsc2');
    

    figure(14)
    hFig = figure(14);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,1); plot(wy_down_range,Fy_down_A(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('F_x [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,3); plot(wy_down_range,Fy_down_A(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('F_y [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,5); plot(wy_down_range,Fy_down_A(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('F_z [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,2); plot(wy_down_range,My_down_A(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('M_x [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,4); plot(wy_down_range,My_down_A(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('M_y [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,6); plot(wy_down_range,My_down_A(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('M_z [N*mm]')
    hold off
    [~,h1] = suplabel('Aerodynamic damping \omega_y down', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Damping/Aero_wy_down.eps'] ,'-depsc2');
    
    
    
    figure(15)
    hFig = figure(15);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,1); plot(wy_up_range,Fy_up_A(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('F_x [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,3); plot(wy_up_range,Fy_up_A(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('F_y [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,5); plot(wy_up_range,Fy_up_A(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('F_z [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,2); plot(wy_up_range,My_up_A(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('M_x [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,4); plot(wy_up_range,My_up_A(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('M_y [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,6); plot(wy_up_range,My_up_A(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('M_z [N*mm]')
    hold off
    [~,h1] = suplabel('Aerodynamic damping \omega_y up', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Damping/Aero_wy_up.eps'] ,'-depsc2');
    
    
    figure(16)
    hFig = figure(16);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,1); plot(wz_range,Fz_A(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_z [rad/s]')
    ylabel('F_x [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,3); plot(wz_range,Fz_A(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_z [rad/s]')
    ylabel('F_y [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,5); plot(wz_range,Fz_A(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_z [rad/s]')
    ylabel('F_z [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,2); plot(wz_range,Mz_A(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_z [rad/s]')
    ylabel('M_x [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,4); plot(wz_range,Mz_A(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_z [rad/s]')
    ylabel('M_y [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,6); plot(wz_range,Mz_A(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_x [rad/s]')
    ylabel('M_z [N*mm]')
    hold off
    [~,h1] = suplabel('Aerodynamic damping \omega_z', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Damping/Aero_wz.eps'] ,'-depsc2');   
    
    
    % Plot inertial damping:
    
    figure(17)
    hFig = figure(17);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,1); plot(vx_forward_range,Fx_forward_I(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('F_x [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,3); plot(vx_forward_range,Fx_forward_I(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('F_y [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,5); plot(vx_forward_range,Fx_forward_I(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('F_z [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,2); plot(vx_forward_range,Mx_forward_I(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('M_x [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,4); plot(vx_forward_range,Mx_forward_I(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('M_y [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,6); plot(vx_forward_range,Mx_forward_I(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('M_z [N*mm]')
    hold off
    [~,h1] = suplabel('Inertial damping v_x forward', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Damping/Inertial_vx_forward.eps'] ,'-depsc2');
    
    figure(18)
    hFig = figure(18);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,1); plot(vx_back_range,Fx_back_I(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('F_x [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,3); plot(vx_back_range,Fx_back_I(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('F_y [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,5); plot(vx_back_range,Fx_back_I(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('F_z [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,2); plot(vx_back_range,Mx_back_I(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('M_x [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,4); plot(vx_back_range,Mx_back_I(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('M_y [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,6); plot(vx_back_range,Mx_back_I(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_x [mm/s]')
    ylabel('M_z [N*mm]')
    hold off
    [~,h1] = suplabel('Inertial damping v_x back', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Damping/Inertial_vx_back.eps'] ,'-depsc2');

    figure(19)
    hFig = figure(19);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,1); plot(vy_range,Fy_I(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_y [mm/s]')
    ylabel('F_x [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,3); plot(vy_range,Fy_I(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_y [mm/s]')
    ylabel('F_y [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,5); plot(vy_range,Fy_I(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_y [mm/s]')
    ylabel('F_z [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,2); plot(vy_range,My_I(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_y [mm/s]')
    ylabel('M_x [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,4); plot(vy_range,My_I(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_y [mm/s]')
    ylabel('M_y [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,6); plot(vy_range,My_I(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_y [mm/s]')
    ylabel('M_z [N*mm]')
    hold off
    [~,h1] = suplabel('Inertial damping v_y', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Damping/Inertial_vy.eps'] ,'-depsc2');
    
    figure(20)
    hFig = figure(20);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,1); plot(vz_up_range,Fz_up_I(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_z [mm/s]')
    ylabel('F_x [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,3); plot(vz_up_range,Fz_up_I(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_z [mm/s]')
    ylabel('F_y [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,5); plot(vz_up_range,Fz_up_I(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_z [mm/s]')
    ylabel('F_z [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,2); plot(vz_up_range,Mz_up_I(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_z [mm/s]')
    ylabel('M_x [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,4); plot(vz_up_range,Mz_up_I(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_z [mm/s]')
    ylabel('M_y [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,6); plot(vz_up_range,Mz_up_I(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('v_z [mm/s]')
    ylabel('M_z [N*mm]')
    hold off
    [~,h1] = suplabel('Inertial damping v_z up', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Damping/Inertial_vz_up.eps'] ,'-depsc2');

    
    figure(21)
    hFig = figure(21);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,1); plot(wx_range,Fx_I(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_x [rad/s]')
    ylabel('F_x [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,3); plot(wx_range,Fx_I(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_x [rad/s]')
    ylabel('F_y [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,5); plot(wx_range,Fx_I(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_x [rad/s]')
    ylabel('F_z [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,2); plot(wx_range,Mx_I(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_x [rad/s]')
    ylabel('M_x [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,4); plot(wx_range,Mx_I(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_x [rad/s]')
    ylabel('M_y [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,6); plot(wx_range,Mx_I(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_x [rad/s]')
    ylabel('M_z [N*mm]')
    hold off
    [~,h1] = suplabel('Inertial damping \omega_x', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Damping/Inertial_wx.eps'] ,'-depsc2');
    

    figure(22)
    hFig = figure(22);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,1); plot(wy_down_range,Fy_down_I(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('F_x [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,3); plot(wy_down_range,Fy_down_I(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('F_y [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,5); plot(wy_down_range,Fy_down_I(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('F_z [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,2); plot(wy_down_range,My_down_I(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('M_x [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,4); plot(wy_down_range,My_down_I(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('M_y [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,6); plot(wy_down_range,My_down_I(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('M_z [N*mm]')
    hold off
    [~,h1] = suplabel('Inertial damping \omega_y down', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Damping/Inertial_wy_down.eps'] ,'-depsc2');
    
    
    
    figure(23)
    hFig = figure(23);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,1); plot(wy_up_range,Fy_up_I(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('F_x [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,3); plot(wy_up_range,Fy_up_I(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('F_y [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,5); plot(wy_up_range,Fy_up_I(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('F_z [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,2); plot(wy_up_range,My_up_I(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('M_x [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,4); plot(wy_up_range,My_up_I(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('M_y [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,6); plot(wy_up_range,My_up_I(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_y [rad/s]')
    ylabel('M_z [N*mm]')
    hold off
    [~,h1] = suplabel('Inertial damping \omega_y up', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Damping/Inertial_wy_up.eps'] ,'-depsc2');
    
    
    figure(24)
    hFig = figure(24);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,1); plot(wz_range,Fz_I(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_z [rad/s]')
    ylabel('F_x [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,3); plot(wz_range,Fz_I(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_z [rad/s]')
    ylabel('F_y [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,5); plot(wz_range,Fz_I(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_z [rad/s]')
    ylabel('F_z [N]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,2); plot(wz_range,Mz_I(1,:,j),'Color',[ (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_z [rad/s]')
    ylabel('M_x [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,4); plot(wz_range,Mz_I(2,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_z [rad/s]')
    ylabel('M_y [N*mm]')
    for j = 1:(nr_steps+1)
        hold on
        subplot(3,2,6); plot(wz_range,Mz_I(3,:,j),'Color',[ (0.5-(0.5*(j-1))/nr_steps) (0.5-(0.5*(j-1))/nr_steps) (0.5+(0.5*(j-1))/nr_steps) ])
        hold off
    end
    xlabel('\omega_x [rad/s]')
    ylabel('M_z [N*mm]')
    hold off
    [~,h1] = suplabel('Inertial damping \omega_z', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Damping/Inertial_wz.eps'] ,'-depsc2'); 

end

