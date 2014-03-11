function Maneuvering_wingkin_plots( settings, pathDB )


    % Test the Inertia model for conservation of linear and angular
    % momentum:
    
    seq_nr      = 1;
    nr_points   = 40;
    
    % Setup the body_model and wing_model:
    
    body_model.mass_fly         = pathDB.body_model.mass_fly(seq_nr);
    body_model.mass_body        = pathDB.body_model.mass_body(seq_nr);
    body_model.Joint_left       = pathDB.body_model.Joint_left(seq_nr,:)';
    body_model.Joint_right      = pathDB.body_model.Joint_right(seq_nr,:)';
    body_model.cg_b             = pathDB.body_model.cg(seq_nr,:)';
    body_model.Inertia          = pathDB.body_model.Inertia(:,:,seq_nr);
    body_model.g                = 1e-3*settings.g;
    body_model.x_mod            = pathDB.body_model.x_mod(:,:,seq_nr);
    body_model.y_mod            = pathDB.body_model.y_mod(:,:,seq_nr);
    body_model.z_mod            = pathDB.body_model.z_mod(:,:,seq_nr);
    
    wing_model.virtual_mass     = pathDB.wing_model.virtual_mass(seq_nr);
    wing_model.wing_cg_L        = pathDB.wing_model.wing_cg_L(seq_nr,:)';
    wing_model.wing_cg_R        = pathDB.wing_model.wing_cg_R(seq_nr,:)';
    wing_model.virtual_Inertia  = pathDB.wing_model.virtual_Inertia(:,:,seq_nr);
    wing_model.y_sect_L         = pathDB.wing_model.y_sect_L(:,:,seq_nr)';
    wing_model.chords_L         = pathDB.wing_model.chords_L(seq_nr,:)';
    wing_model.y_sect_R         = pathDB.wing_model.y_sect_R(:,:,seq_nr)';
    wing_model.chords_R         = pathDB.wing_model.chords_R(seq_nr,:)';
    wing_model.rho              = settings.rho_air;
    wing_model.length           = pathDB.wing_model.length(seq_nr);
    
    f                           = pathDB.poly_fit.a_glob.f;
    dt                          = (1/f)/(nr_points-1);
    R_strk                      = pathDB.rot_mat.Rstr;
    
    % Create the wing kinematics:
    
    n_pol_theta     = settings.n_pol_theta;
    n_pol_eta       = settings.n_pol_eta;
    n_pol_phi       = settings.n_pol_phi;

    a_avg_theta    = pathDB.poly_fit.a_glob.theta;
    a_avg_eta      = pathDB.poly_fit.a_glob.eta;
    a_avg_phi      = pathDB.poly_fit.a_glob.phi;
    down_up_avg    = pathDB.poly_fit.a_glob.down_up;

    a_theta_L_avg  = a_avg_theta;
    a_eta_L_avg    = a_avg_eta;
    a_phi_L_avg    = a_avg_phi;
    a_theta_R_avg  = a_avg_theta;
    a_eta_R_avg    = a_avg_eta;
    a_phi_R_avg    = a_avg_phi;
       
    % Create the average wing kinematics:

    a_fit_avg.a_theta_L     = a_theta_L_avg;
    a_fit_avg.a_eta_L       = a_eta_L_avg;
    a_fit_avg.a_phi_L       = a_phi_L_avg;
    a_fit_avg.a_theta_R     = a_theta_R_avg;
    a_fit_avg.a_eta_R       = a_eta_R_avg;
    a_fit_avg.a_phi_R       = a_phi_R_avg;
    a_fit_avg.f             = f;
    a_fit_avg.down_up       = down_up_avg;
    a_fit_avg.nr_points     = nr_points;
    a_fit_avg.R_strk        = R_strk;
        
    [ kine_avg ] = angular_velocities_polynomial( a_fit_avg );
    
    
    F_ax_forward_max        = pathDB.maneuver.c_fit_ax.FM.FM_1(end);
    F_ax_back_max           = pathDB.maneuver.c_fit_ax.FM.FM_1(1);
    F_ay_max                = pathDB.maneuver.c_fit_ay.FM.FM_2(end);
    F_az_up_max             = pathDB.maneuver.c_fit_az.FM.FM_3(1);
    M_wx_max                = pathDB.maneuver.c_fit_wx.FM.FM_4(end);
    M_wy_down_max           = pathDB.maneuver.c_fit_wy.FM.FM_5(1);
    M_wy_up_max             = pathDB.maneuver.c_fit_wy.FM.FM_5(end);
    M_wz_max                = pathDB.maneuver.c_fit_wz.FM.FM_6(end);
        
%     F_ax_forward_range      = 0:(F_ax_forward_max/nr_steps):F_ax_forward_max;
%     F_ax_back_range         = 0:(F_ax_back_max/nr_steps):F_ax_back_max;
%     F_ay_range              = 0:(F_ay_max/nr_steps):F_ay_max;
%     F_az_up_range           = 0:(F_az_up_max/nr_steps):F_az_up_max;
%     M_wx_range              = 0:(M_wx_max/nr_steps):M_wx_max;
%     M_wy_down_range         = 0:(M_wy_down_max/nr_steps):M_wy_down_max;
%     M_wy_up_range           = 0:(M_wy_up_max/nr_steps):M_wy_up_max;
%     M_wz_range              = 0:(M_wz_max/nr_steps):M_wz_max;

    for i = 1:8
        
        if i == 1
            
            b_theta_L_man  = pathDB.maneuver.c_fit_ax.b_theta_p(:,1,1);
            b_eta_L_man    = pathDB.maneuver.c_fit_ax.b_eta_p(:,1,1);
            b_phi_L_man    = pathDB.maneuver.c_fit_ax.b_phi_p(:,1,1);
            b_theta_R_man  = pathDB.maneuver.c_fit_ax.b_theta_p(:,1,1);
            b_eta_R_man    = pathDB.maneuver.c_fit_ax.b_eta_p(:,1,1);
            b_phi_R_man    = pathDB.maneuver.c_fit_ax.b_phi_p(:,1,1);
    
            a_theta_L_man  = a_avg_theta+b_theta_L_man*F_ax_forward_max;
            a_eta_L_man    = a_avg_eta+b_eta_L_man*F_ax_forward_max;
            a_phi_L_man    = a_avg_phi+b_phi_L_man*F_ax_forward_max;
            a_theta_R_man  = a_avg_theta+b_theta_R_man*F_ax_forward_max;
            a_eta_R_man    = a_avg_eta+b_eta_R_man*F_ax_forward_max;
            a_phi_R_man    = a_avg_phi+b_phi_R_man*F_ax_forward_max;
            
            man_title_name      = 'Fx forward maneuver';
            man_save_name       = 'Fx_forward';
            
        elseif i == 2
            
            b_theta_L_man  = pathDB.maneuver.c_fit_ax.b_theta_n(:,1,1);
            b_eta_L_man    = pathDB.maneuver.c_fit_ax.b_eta_n(:,1,1);
            b_phi_L_man    = pathDB.maneuver.c_fit_ax.b_phi_n(:,1,1);
            b_theta_R_man  = pathDB.maneuver.c_fit_ax.b_theta_n(:,1,1);
            b_eta_R_man    = pathDB.maneuver.c_fit_ax.b_eta_n(:,1,1);
            b_phi_R_man    = pathDB.maneuver.c_fit_ax.b_phi_n(:,1,1);
    
            a_theta_L_man  = a_avg_theta+b_theta_L_man*F_ax_back_max;
            a_eta_L_man    = a_avg_eta+b_eta_L_man*F_ax_back_max;
            a_phi_L_man    = a_avg_phi+b_phi_L_man*F_ax_back_max;
            a_theta_R_man  = a_avg_theta+b_theta_R_man*F_ax_back_max;
            a_eta_R_man    = a_avg_eta+b_eta_R_man*F_ax_back_max;
            a_phi_R_man    = a_avg_phi+b_phi_R_man*F_ax_back_max;      
            
            man_title_name      = 'Fx back maneuver';
            man_save_name       = 'Fx_back';
            
        elseif i == 3
            
            b_theta_L_man  = pathDB.maneuver.c_fit_ay.b_theta_L(:,1,2);
            b_eta_L_man    = pathDB.maneuver.c_fit_ay.b_eta_L(:,1,2);
            b_phi_L_man    = pathDB.maneuver.c_fit_ay.b_phi_L(:,1,2);
            b_theta_R_man  = pathDB.maneuver.c_fit_ay.b_theta_R(:,1,2);
            b_eta_R_man    = pathDB.maneuver.c_fit_ay.b_eta_R(:,1,2);
            b_phi_R_man    = pathDB.maneuver.c_fit_ay.b_phi_R(:,1,2);
    
            a_theta_L_man  = a_avg_theta+b_theta_L_man*F_ay_max;
            a_eta_L_man    = a_avg_eta+b_eta_L_man*F_ay_max;
            a_phi_L_man    = a_avg_phi+b_phi_L_man*F_ay_max;
            a_theta_R_man  = a_avg_theta+b_theta_R_man*F_ay_max;
            a_eta_R_man    = a_avg_eta+b_eta_R_man*F_ay_max;
            a_phi_R_man    = a_avg_phi+b_phi_R_man*F_ay_max;      
            
            man_title_name      = 'Fy maneuver';
            man_save_name       = 'Fy';
            
        elseif i == 4
            
            b_theta_L_man  = pathDB.maneuver.c_fit_az.b_theta_n(:,1,3);
            b_eta_L_man    = pathDB.maneuver.c_fit_az.b_eta_n(:,1,3);
            b_phi_L_man    = pathDB.maneuver.c_fit_az.b_phi_n(:,1,3);
            b_theta_R_man  = pathDB.maneuver.c_fit_az.b_theta_n(:,1,3);
            b_eta_R_man    = pathDB.maneuver.c_fit_az.b_eta_n(:,1,3);
            b_phi_R_man    = pathDB.maneuver.c_fit_az.b_phi_n(:,1,3);
    
            a_theta_L_man  = a_avg_theta+b_theta_L_man*F_az_up_max;
            a_eta_L_man    = a_avg_eta+b_eta_L_man*F_az_up_max;
            a_phi_L_man    = a_avg_phi+b_phi_L_man*F_az_up_max;
            a_theta_R_man  = a_avg_theta+b_theta_R_man*F_az_up_max;
            a_eta_R_man    = a_avg_eta+b_eta_R_man*F_az_up_max;
            a_phi_R_man    = a_avg_phi+b_phi_R_man*F_az_up_max; 
            
            man_title_name      = 'Fz up maneuver';
            man_save_name       = 'Fz_up';
            
        elseif i == 5
            
            b_theta_L_man  = pathDB.maneuver.c_fit_wx.b_theta_L(:,1,4);
            b_eta_L_man    = pathDB.maneuver.c_fit_wx.b_eta_L(:,1,4);
            b_phi_L_man    = pathDB.maneuver.c_fit_wx.b_phi_L(:,1,4);
            b_theta_R_man  = pathDB.maneuver.c_fit_wx.b_theta_R(:,1,4);
            b_eta_R_man    = pathDB.maneuver.c_fit_wx.b_eta_R(:,1,4);
            b_phi_R_man    = pathDB.maneuver.c_fit_wx.b_phi_R(:,1,4);
    
            a_theta_L_man  = a_avg_theta+b_theta_L_man*M_wx_max;
            a_eta_L_man    = a_avg_eta+b_eta_L_man*M_wx_max;
            a_phi_L_man    = a_avg_phi+b_phi_L_man*M_wx_max;
            a_theta_R_man  = a_avg_theta+b_theta_R_man*M_wx_max;
            a_eta_R_man    = a_avg_eta+b_eta_R_man*M_wx_max;
            a_phi_R_man    = a_avg_phi+b_phi_R_man*M_wx_max;
            
            man_title_name      = 'Roll maneuver';
            man_save_name       = 'Mx';
            
        elseif i == 6
            
            b_theta_L_man  = pathDB.maneuver.c_fit_wy.b_theta_n(:,1,5);
            b_eta_L_man    = pathDB.maneuver.c_fit_wy.b_eta_n(:,1,5);
            b_phi_L_man    = pathDB.maneuver.c_fit_wy.b_phi_n(:,1,5);
            b_theta_R_man  = pathDB.maneuver.c_fit_wy.b_theta_n(:,1,5);
            b_eta_R_man    = pathDB.maneuver.c_fit_wy.b_eta_n(:,1,5);
            b_phi_R_man    = pathDB.maneuver.c_fit_wy.b_phi_n(:,1,5);
    
            a_theta_L_man  = a_avg_theta+b_theta_L_man*M_wy_down_max;
            a_eta_L_man    = a_avg_eta+b_eta_L_man*M_wy_down_max;
            a_phi_L_man    = a_avg_phi+b_phi_L_man*M_wy_down_max;
            a_theta_R_man  = a_avg_theta+b_theta_R_man*M_wy_down_max;
            a_eta_R_man    = a_avg_eta+b_eta_R_man*M_wy_down_max;
            a_phi_R_man    = a_avg_phi+b_phi_R_man*M_wy_down_max;
            
            man_title_name      = 'Pitch down maneuver';
            man_save_name       = 'My_down';
            
        elseif i == 7
            
            b_theta_L_man  = pathDB.maneuver.c_fit_wy.b_theta_p(:,1,5);
            b_eta_L_man    = pathDB.maneuver.c_fit_wy.b_eta_p(:,1,5);
            b_phi_L_man    = pathDB.maneuver.c_fit_wy.b_phi_p(:,1,5);
            b_theta_R_man  = pathDB.maneuver.c_fit_wy.b_theta_p(:,1,5);
            b_eta_R_man    = pathDB.maneuver.c_fit_wy.b_eta_p(:,1,5);
            b_phi_R_man    = pathDB.maneuver.c_fit_wy.b_phi_p(:,1,5);
    
            a_theta_L_man  = a_avg_theta+b_theta_L_man*M_wy_up_max;
            a_eta_L_man    = a_avg_eta+b_eta_L_man*M_wy_up_max;
            a_phi_L_man    = a_avg_phi+b_phi_L_man*M_wy_up_max;
            a_theta_R_man  = a_avg_theta+b_theta_R_man*M_wy_up_max;
            a_eta_R_man    = a_avg_eta+b_eta_R_man*M_wy_up_max;
            a_phi_R_man    = a_avg_phi+b_phi_R_man*M_wy_up_max;
            
            man_title_name      = 'Pitch up maneuver';
            man_save_name       = 'My_up';
            
        elseif i == 8
            
            b_theta_L_man  = pathDB.maneuver.c_fit_wz.b_theta_L(:,1,6);
            b_eta_L_man    = pathDB.maneuver.c_fit_wz.b_eta_L(:,1,6);
            b_phi_L_man    = pathDB.maneuver.c_fit_wz.b_phi_L(:,1,6);
            b_theta_R_man  = pathDB.maneuver.c_fit_wz.b_theta_R(:,1,6);
            b_eta_R_man    = pathDB.maneuver.c_fit_wz.b_eta_R(:,1,6);
            b_phi_R_man    = pathDB.maneuver.c_fit_wz.b_phi_R(:,1,6);
    
            a_theta_L_man  = a_avg_theta+b_theta_L_man*M_wz_max;
            a_eta_L_man    = a_avg_eta+b_eta_L_man*M_wz_max;
            a_phi_L_man    = a_avg_phi+b_phi_L_man*M_wz_max;
            a_theta_R_man  = a_avg_theta+b_theta_R_man*M_wz_max;
            a_eta_R_man    = a_avg_eta+b_eta_R_man*M_wz_max;
            a_phi_R_man    = a_avg_phi+b_phi_R_man*M_wz_max; 
            
            man_title_name      = 'Yaw maneuver';
            man_save_name       = 'Mz';
            
        end

        
    a_fit_man.a_theta_L     = a_theta_L_man;
    a_fit_man.a_eta_L       = a_eta_L_man;
    a_fit_man.a_phi_L       = a_phi_L_man;
    a_fit_man.a_theta_R     = a_theta_R_man;
    a_fit_man.a_eta_R       = a_eta_R_man;
    a_fit_man.a_phi_R       = a_phi_R_man;
    a_fit_man.f             = f;
    a_fit_man.down_up       = down_up_avg;
    a_fit_man.nr_points     = nr_points;
    a_fit_man.R_strk        = R_strk;
        
    [ kine_man ] = angular_velocities_polynomial( a_fit_man );
    
    
    RL_avg = kine_avg.RL;
    RR_avg = kine_avg.RR;
    
    RL_man = kine_man.RL;
    RR_man = kine_man.RR;
    
    minx = -3.5;
    maxx = 3.5;
    miny = -3.5;
    maxy = 3.5;
    minz = -3.5;
    maxz = 3.5;
    
    figure(2*(i-1)+1)
    g1 = figure(2*(i-1)+1);
    hold on
    Fly_plot_maneuvering_wing_kinematics( [0; 0; 0], R_strk'*[1 0 0; 0 -1 0; 0 0 -1], RL_avg, RR_avg, down_up_avg, RL_man, RR_man, down_up_avg, body_model, wing_model, 0)
    axis([minx maxx miny maxy minz maxz])
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    hold off
    saveas(g1,[char(settings.plot_loc) '/Manuevering_wing_kin2/' man_save_name '_3D' ],'fig')


    figure(2*i)
    hFig = figure(2*i);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    subplot(2,2,1); Fly_plot_maneuvering_wing_kinematics( [0; 0; 0], R_strk'*[1 0 0; 0 -1 0; 0 0 -1], RL_avg, RR_avg, down_up_avg, RL_man, RR_man, down_up_avg, body_model, wing_model, 1)
    axis([minx maxx miny maxy minz maxz])
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    view(180,0)
    subplot(2,2,2); Fly_plot_maneuvering_wing_kinematics( [0; 0; 0], R_strk'*[1 0 0; 0 -1 0; 0 0 -1], RL_avg, RR_avg, down_up_avg, RL_man, RR_man, down_up_avg, body_model, wing_model, 2)
    axis([minx maxx miny maxy minz maxz])
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    view(0,0)
    subplot(2,2,3); Fly_plot_maneuvering_wing_kinematics( [0; 0; 0], R_strk'*[1 0 0; 0 -1 0; 0 0 -1], RL_avg, RR_avg, down_up_avg, RL_man, RR_man, down_up_avg, body_model, wing_model, 0)
    axis([minx maxx miny maxy minz maxz])
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    view(90,0)
    subplot(2,2,4); Fly_plot_maneuvering_wing_kinematics( [0; 0; 0], R_strk'*[1 0 0; 0 -1 0; 0 0 -1], RL_avg, RR_avg, down_up_avg, RL_man, RR_man, down_up_avg, body_model, wing_model, 0)
    axis([minx maxx miny maxy minz maxz])
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    view(-90,90)
    hold off
    set(2*i,'Color','k')
    [~,h1] = suplabel(man_title_name, 't');
    set(h1,'FontSize',10)
    set(h1,'Color',[0.7 0.7 0.7])
    set(gcf, 'InvertHardCopy', 'off');
    print ([char(settings.plot_loc) '/Manuevering_wing_kin2/' man_save_name '_2D.eps'] ,'-depsc2');
    saveas(hFig,[char(settings.plot_loc) '/Manuevering_wing_kin2/' man_save_name '_2D'],'fig')

    
    end

end

