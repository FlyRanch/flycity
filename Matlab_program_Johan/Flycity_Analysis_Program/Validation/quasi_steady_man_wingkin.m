function [ FM ] = quasi_steady_man_wingkin( settings, pathDB, rot_on )


    nr_steps = 19;

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
    
    nr_points       = 200;

        
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
    
%     F_ax_forward_max        = 0.5038e-4;
%     F_ax_back_max           = -0.4015e-4;
%     F_ay_max                = 0.4685e-4;
%     F_az_up_max             = -0.7458e-4;
%     M_wx_max                = 0.2305e-4;
%     M_wy_down_max           = -0.2675e-4;
%     M_wy_up_max             = 0.2692e-4;
%     M_wz_max                = 0.4288e-4;
    
    F_ax_forward_max        = 0.5038e-4;
    F_ax_back_max           = -0.4015e-4;
    F_ay_max                = 0.4685e-4;
    F_az_up_max             = -0.6561e-4;
    M_wx_max                = 0.2756e-4;
    M_wy_down_max           = -0.3203e-4;
    M_wy_up_max             = 0.2146e-4;
    M_wz_max                = 0.4255e-4;
        
    F_ax_forward_range      = 0:(F_ax_forward_max/nr_steps):F_ax_forward_max;
    F_ax_back_range         = 0:(F_ax_back_max/nr_steps):F_ax_back_max;
    F_ay_range              = 0:(F_ay_max/nr_steps):F_ay_max;
    F_az_up_range           = 0:(F_az_up_max/nr_steps):F_az_up_max;
    M_wx_range              = 0:(M_wx_max/nr_steps):M_wx_max;
    M_wy_down_range         = 0:(M_wy_down_max/nr_steps):M_wy_down_max;
    M_wy_up_range           = 0:(M_wy_up_max/nr_steps):M_wy_up_max;
    M_wz_range              = 0:(M_wz_max/nr_steps):M_wz_max;
        
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
    wing_model.x_LE_L           = pathDB.wing_model.x_LE_L(seq_nr,:)';
    wing_model.y_sect_R         = pathDB.wing_model.y_sect_R(:,:,seq_nr)';
    wing_model.chords_R         = pathDB.wing_model.chords_R(seq_nr,:)';
    wing_model.x_LE_R           = pathDB.wing_model.x_LE_R(seq_nr,:)';
    wing_model.x_mod_L          = pathDB.wing_model.x_mod_L(:,:,seq_nr);
    wing_model.y_mod_L          = pathDB.wing_model.y_mod_L(:,:,seq_nr);
    wing_model.z_mod_L          = pathDB.wing_model.z_mod_L(:,:,seq_nr);
    wing_model.x_mod_R          = pathDB.wing_model.x_mod_R(:,:,seq_nr);
    wing_model.y_mod_R          = pathDB.wing_model.y_mod_R(:,:,seq_nr);
    wing_model.z_mod_R          = pathDB.wing_model.z_mod_R(:,:,seq_nr);
    wing_model.rho              = settings.rho_air;
    wing_model.nr_sect          = size(wing_model.y_sect_L,2);
    
    
    % Compute the global quasi-steady forces and moments:
    
    a_fit.a_theta_L     = a_avg_theta;
    a_fit.a_eta_L       = a_avg_eta;
    a_fit.a_phi_L       = a_avg_phi;
    a_fit.a_theta_R     = a_avg_theta;
    a_fit.a_eta_R       = a_avg_eta;
    a_fit.a_phi_R       = a_avg_phi;
    a_fit.f             = f;
    a_fit.down_up       = down_up;
    a_fit.nr_points     = 2000;
    a_fit.R_strk        = R_strk;

    [ wingkin ] = angular_velocities_polynomial( a_fit );

    kine.R_strk          = R_strk;
    kine.u_strk          = zeros(3,a_fit.nr_points);
    kine.w_strk          = zeros(3,a_fit.nr_points);
    kine.wL              = wingkin.wL;
    kine.wR              = wingkin.wR;
    kine.w_dot_L         = wingkin.w_dot_L;
    kine.w_dot_R         = wingkin.w_dot_R;
    kine.RL              = wingkin.RL;
    kine.RR              = wingkin.RR;

    [ FM_strkpln, ~, ~ , ~, ~, ~, ~, ~, ~ ] = Aerodynamic_forces( kine, body_model, wing_model, rot_on);

    FM_glob = FM_strkpln;
    
    
    clear a_fit FM_strkpln wingkin kine
    
    % Store Forces and moments:
    
    FM_ax_forward       = zeros(6,nr_steps+1);
    
    'FM_ax_forward'
   
    for i = 1:(nr_steps+1)

        a_fit.a_theta_L     = a_x_forward.theta_L(:,i);
        a_fit.a_eta_L       = a_x_forward.eta_L(:,i);
        a_fit.a_phi_L       = a_x_forward.phi_L(:,i);
        a_fit.a_theta_R     = a_x_forward.theta_R(:,i);
        a_fit.a_eta_R       = a_x_forward.eta_R(:,i);
        a_fit.a_phi_R       = a_x_forward.phi_R(:,i);
        a_fit.f             = f;
        a_fit.down_up       = down_up;
        a_fit.nr_points     = nr_points;
        a_fit.R_strk        = R_strk;

        [ wingkin ] = angular_velocities_polynomial( a_fit );

        kine.R_strk          = R_strk;
        kine.u_strk          = zeros(3,nr_points);
        kine.w_strk          = zeros(3,nr_points);
        kine.wL              = wingkin.wL;
        kine.wR              = wingkin.wR;
        kine.w_dot_L         = wingkin.w_dot_L;
        kine.w_dot_R         = wingkin.w_dot_R;
        kine.RL              = wingkin.RL;
        kine.RR              = wingkin.RR;

        [ FM_strkpln, ~, ~ , ~, ~, ~, ~, ~, ~ ] = Aerodynamic_forces( kine, body_model, wing_model, rot_on);

        FM_ax_forward(:,i) = mean(FM_strkpln,2);
    
    end
    
    FM_ax_back          = zeros(6,nr_steps+1);
    
    'FM_ax_back'
    
    for i = 1:(nr_steps+1)

        a_fit.a_theta_L     = a_x_back.theta_L(:,i);
        a_fit.a_eta_L       = a_x_back.eta_L(:,i);
        a_fit.a_phi_L       = a_x_back.phi_L(:,i);
        a_fit.a_theta_R     = a_x_back.theta_R(:,i);
        a_fit.a_eta_R       = a_x_back.eta_R(:,i);
        a_fit.a_phi_R       = a_x_back.phi_R(:,i);
        a_fit.f             = f;
        a_fit.down_up       = down_up;
        a_fit.nr_points     = nr_points;
        a_fit.R_strk        = R_strk;

        [ wingkin ] = angular_velocities_polynomial( a_fit );

        kine.R_strk          = R_strk;
        kine.u_strk          = zeros(3,nr_points);
        kine.w_strk          = zeros(3,nr_points);
        kine.wL              = wingkin.wL;
        kine.wR              = wingkin.wR;
        kine.w_dot_L         = wingkin.w_dot_L;
        kine.w_dot_R         = wingkin.w_dot_R;
        kine.RL              = wingkin.RL;
        kine.RR              = wingkin.RR;

        [ FM_strkpln, ~, ~ , ~, ~, ~, ~, ~, ~ ] = Aerodynamic_forces( kine, body_model, wing_model, rot_on);

        FM_ax_back(:,i) = mean(FM_strkpln,2);
    
    end
    
    FM_ay               = zeros(6,nr_steps+1);
    
    'FM_ay'
    
    for i = 1:(nr_steps+1)

        a_fit.a_theta_L     = a_y.theta_L(:,i);
        a_fit.a_eta_L       = a_y.eta_L(:,i);
        a_fit.a_phi_L       = a_y.phi_L(:,i);
        a_fit.a_theta_R     = a_y.theta_R(:,i);
        a_fit.a_eta_R       = a_y.eta_R(:,i);
        a_fit.a_phi_R       = a_y.phi_R(:,i);
        a_fit.f             = f;
        a_fit.down_up       = down_up;
        a_fit.nr_points     = nr_points;
        a_fit.R_strk        = R_strk;

        [ wingkin ] = angular_velocities_polynomial( a_fit );

        kine.R_strk          = R_strk;
        kine.u_strk          = zeros(3,nr_points);
        kine.w_strk          = zeros(3,nr_points);
        kine.wL              = wingkin.wL;
        kine.wR              = wingkin.wR;
        kine.w_dot_L         = wingkin.w_dot_L;
        kine.w_dot_R         = wingkin.w_dot_R;
        kine.RL              = wingkin.RL;
        kine.RR              = wingkin.RR;

        [ FM_strkpln, ~, ~ , ~, ~, ~, ~, ~, ~ ] = Aerodynamic_forces( kine, body_model, wing_model, rot_on);

        FM_ay(:,i) = mean(FM_strkpln,2);
    
    end
    
    FM_az_up            = zeros(6,nr_steps+1);
    
    'FM_az_up'
    
	for i = 1:(nr_steps+1)

        a_fit.a_theta_L     = a_z_up.theta_L(:,i);
        a_fit.a_eta_L       = a_z_up.eta_L(:,i);
        a_fit.a_phi_L       = a_z_up.phi_L(:,i);
        a_fit.a_theta_R     = a_z_up.theta_R(:,i);
        a_fit.a_eta_R       = a_z_up.eta_R(:,i);
        a_fit.a_phi_R       = a_z_up.phi_R(:,i);
        a_fit.f             = f;
        a_fit.down_up       = down_up;
        a_fit.nr_points     = nr_points;
        a_fit.R_strk        = R_strk;

        [ wingkin ] = angular_velocities_polynomial( a_fit );

        kine.R_strk          = R_strk;
        kine.u_strk          = zeros(3,nr_points);
        kine.w_strk          = zeros(3,nr_points);
        kine.wL              = wingkin.wL;
        kine.wR              = wingkin.wR;
        kine.w_dot_L         = wingkin.w_dot_L;
        kine.w_dot_R         = wingkin.w_dot_R;
        kine.RL              = wingkin.RL;
        kine.RR              = wingkin.RR;

        [ FM_strkpln, ~, ~ , ~, ~, ~, ~, ~, ~ ] = Aerodynamic_forces( kine, body_model, wing_model, rot_on);

        FM_az_up(:,i) = mean(FM_strkpln,2);
    
	end
    
    FM_wx               = zeros(6,nr_steps+1);
    
    'FM_wx'
    
    for i = 1:(nr_steps+1)

        a_fit.a_theta_L     = w_x.theta_L(:,i);
        a_fit.a_eta_L       = w_x.eta_L(:,i);
        a_fit.a_phi_L       = w_x.phi_L(:,i);
        a_fit.a_theta_R     = w_x.theta_R(:,i);
        a_fit.a_eta_R       = w_x.eta_R(:,i);
        a_fit.a_phi_R       = w_x.phi_R(:,i);
        a_fit.f             = f;
        a_fit.down_up       = down_up;
        a_fit.nr_points     = nr_points;
        a_fit.R_strk        = R_strk;

        [ wingkin ] = angular_velocities_polynomial( a_fit );

        kine.R_strk          = R_strk;
        kine.u_strk          = zeros(3,nr_points);
        kine.w_strk          = zeros(3,nr_points);
        kine.wL              = wingkin.wL;
        kine.wR              = wingkin.wR;
        kine.w_dot_L         = wingkin.w_dot_L;
        kine.w_dot_R         = wingkin.w_dot_R;
        kine.RL              = wingkin.RL;
        kine.RR              = wingkin.RR;

        [ FM_strkpln, ~, ~ , ~, ~, ~, ~, ~, ~ ] = Aerodynamic_forces( kine, body_model, wing_model, rot_on);

        FM_wx(:,i) = mean(FM_strkpln,2);

    end
    
    FM_wy_down          = zeros(6,nr_steps+1);
    
    'FM_wy_down'
    
    for i = 1:(nr_steps+1)

        a_fit.a_theta_L     = w_y_down.theta_L(:,i);
        a_fit.a_eta_L       = w_y_down.eta_L(:,i);
        a_fit.a_phi_L       = w_y_down.phi_L(:,i);
        a_fit.a_theta_R     = w_y_down.theta_R(:,i);
        a_fit.a_eta_R       = w_y_down.eta_R(:,i);
        a_fit.a_phi_R       = w_y_down.phi_R(:,i);
        a_fit.f             = f;
        a_fit.down_up       = down_up;
        a_fit.nr_points     = nr_points;
        a_fit.R_strk        = R_strk;

        [ wingkin ] = angular_velocities_polynomial( a_fit );

        kine.R_strk          = R_strk;
        kine.u_strk          = zeros(3,nr_points);
        kine.w_strk          = zeros(3,nr_points);
        kine.wL              = wingkin.wL;
        kine.wR              = wingkin.wR;
        kine.w_dot_L         = wingkin.w_dot_L;
        kine.w_dot_R         = wingkin.w_dot_R;
        kine.RL              = wingkin.RL;
        kine.RR              = wingkin.RR;

        [ FM_strkpln, ~, ~ , ~, ~, ~, ~, ~, ~ ] = Aerodynamic_forces( kine, body_model, wing_model, rot_on);

        FM_wy_down(:,i) = mean(FM_strkpln,2);
    
    end
    
    FM_wy_up            = zeros(6,nr_steps+1);
    
    'FM_wy_up'
    
    for i = 1:(nr_steps+1)

        a_fit.a_theta_L     = w_y_up.theta_L(:,i);
        a_fit.a_eta_L       = w_y_up.eta_L(:,i);
        a_fit.a_phi_L       = w_y_up.phi_L(:,i);
        a_fit.a_theta_R     = w_y_up.theta_R(:,i);
        a_fit.a_eta_R       = w_y_up.eta_R(:,i);
        a_fit.a_phi_R       = w_y_up.phi_R(:,i);
        a_fit.f             = f;
        a_fit.down_up       = down_up;
        a_fit.nr_points     = nr_points;
        a_fit.R_strk        = R_strk;

        [ wingkin ] = angular_velocities_polynomial( a_fit );

        kine.R_strk          = R_strk;
        kine.u_strk          = zeros(3,nr_points);
        kine.w_strk          = zeros(3,nr_points);
        kine.wL              = wingkin.wL;
        kine.wR              = wingkin.wR;
        kine.w_dot_L         = wingkin.w_dot_L;
        kine.w_dot_R         = wingkin.w_dot_R;
        kine.RL              = wingkin.RL;
        kine.RR              = wingkin.RR;

        [ FM_strkpln, ~, ~ , ~, ~, ~, ~, ~, ~ ] = Aerodynamic_forces( kine, body_model, wing_model, rot_on);

        FM_wy_up(:,i) = mean(FM_strkpln,2);
    
    end
    
    FM_wz               = zeros(6,nr_steps+1);
    
    'FM_wz'
    
    for i = 1:(nr_steps+1)

        a_fit.a_theta_L     = w_z.theta_L(:,i);
        a_fit.a_eta_L       = w_z.eta_L(:,i);
        a_fit.a_phi_L       = w_z.phi_L(:,i);
        a_fit.a_theta_R     = w_z.theta_R(:,i);
        a_fit.a_eta_R       = w_z.eta_R(:,i);
        a_fit.a_phi_R       = w_z.phi_R(:,i);
        a_fit.f             = f;
        a_fit.down_up       = down_up;
        a_fit.nr_points     = nr_points;
        a_fit.R_strk        = R_strk;

        [ wingkin ] = angular_velocities_polynomial( a_fit );

        kine.R_strk          = R_strk;
        kine.u_strk          = zeros(3,nr_points);
        kine.w_strk          = zeros(3,nr_points);
        kine.wL              = wingkin.wL;
        kine.wR              = wingkin.wR;
        kine.w_dot_L         = wingkin.w_dot_L;
        kine.w_dot_R         = wingkin.w_dot_R;
        kine.RL              = wingkin.RL;
        kine.RR              = wingkin.RR;

        [ FM_strkpln, ~, ~ , ~, ~, ~, ~, ~, ~ ] = Aerodynamic_forces( kine, body_model, wing_model, rot_on);

        FM_wz(:,i) = mean(FM_strkpln,2);
    
    end
    
    FM.glob         = FM_glob;
    FM.ax_forward   = FM_ax_forward;
    FM.ax_back      = FM_ax_back;
    FM.ay           = FM_ay;
    FM.az_up        = FM_az_up;
    FM.wx           = FM_wx;
    FM.wy_down      = FM_wy_down;
    FM.wy_up        = FM_wy_up;
    FM.wz           = FM_wz;
    
end

