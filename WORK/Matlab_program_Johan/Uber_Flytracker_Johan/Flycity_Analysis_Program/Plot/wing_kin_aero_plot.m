function wing_kin_aero_plot( settings, pathDB, seq_nr, wb_nr)

    % Plot the 3D wing kinematics in combination with the local aerodynamic
    % forces:
    
    wb_domain = pathDB.wingbeats.wingbeat_loc(wb_nr,1,seq_nr):pathDB.wingbeats.wingbeat_loc(wb_nr,2,seq_nr);
    
    nr_points   = length(wb_domain);
    
    R_strk = pathDB.rot_mat.Rstr;
    
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
    
    % Create the wing kinematics

    a_fit.a_theta_L     = pathDB.poly_fit.a_fit.theta_L(:,wb_nr,seq_nr);
    a_fit.a_eta_L       = pathDB.poly_fit.a_fit.eta_L(:,wb_nr,seq_nr);
    a_fit.a_phi_L       = pathDB.poly_fit.a_fit.phi_L(:,wb_nr,seq_nr);
    a_fit.a_theta_R     = pathDB.poly_fit.a_fit.theta_R(:,wb_nr,seq_nr);
    a_fit.a_eta_R       = pathDB.poly_fit.a_fit.eta_R(:,wb_nr,seq_nr);
    a_fit.a_phi_R       = pathDB.poly_fit.a_fit.phi_R(:,wb_nr,seq_nr);
    a_fit.f             = pathDB.poly_fit.a_fit.f(wb_nr,seq_nr);
    a_fit.down_up       = pathDB.poly_fit.a_fit.down_up(wb_nr,seq_nr);
    a_fit.nr_points     = nr_points;
    a_fit.R_strk        = R_strk;
        
    [ kine ] = angular_velocities_polynomial( a_fit );
    
    
    kine.R_strk = R_strk;
    kine.u_strk = pathDB.strkpln_kin.uvw(wb_domain,:,seq_nr)';
    kine.w_strk = pathDB.strkpln_kin.w(wb_domain,:,seq_nr)';
    
    
    % Compute aerodynamic forces:
    
    [ ~, FM_L, FM_R ,~, ~, ~, ~, ~, ~ ] = Aerodynamic_forces( kine, body_model, wing_model, 1);

    FA_L = FM_L(1:3,:);
    FA_R = FM_R(1:3,:);
    
    % Plot the wingbeat:
    
%     qb = pathDB.filt.qB(wb_domain,:,seq_nr);
%     [ qb_avg ] = q_avg(qb(:,1), qb(:,2) ,qb(:,3), qb(:,4));
%     Rb_avg = quat2mat(qb_avg);

    wb_loc = pathDB.wingbeats.wingbeat_loc(wb_nr,:,seq_nr);

    r_avg = mean(pathDB.strkpln_kin.roll(seq_nr,wb_loc(1):wb_loc(2)));
    p_avg = mean(pathDB.strkpln_kin.pitch(seq_nr,wb_loc(1):wb_loc(2)));
    y_avg = mean(pathDB.strkpln_kin.yaw(seq_nr,wb_loc(1):wb_loc(2)));

    Rb_1 = R_strk'*[cos(p_avg) 0 -sin(p_avg); 0 1 0; sin(p_avg) 0 cos(p_avg)]*[1 0 0; 0 -1 0; 0 0 -1];
    Rb_2 = R_strk'*[cos(p_avg) 0 -sin(p_avg); 0 1 0; sin(p_avg) 0 cos(p_avg)]*[1 0 0; 0 -1 0; 0 0 -1];
    Rb_3 = R_strk'*[1 0 0; 0 cos(r_avg) sin(r_avg); 0 -sin(r_avg) cos(r_avg)]*[1 0 0; 0 -1 0; 0 0 -1];
    Rb_4 = R_strk'*[cos(y_avg) sin(y_avg) 0; -sin(y_avg) cos(y_avg) 0; 0 0 1]*[1 0 0; 0 -1 0; 0 0 -1];

    minx = -4;
    maxx = 4;
    miny = -4;
    maxy = 4;
    minz = -4;
    maxz = 4;
    
    hold on
    subplot(2,3,1); wingkin_plot( [0; 0; 0], Rb_1, kine.RL, kine.RR, a_fit.down_up, body_model, wing_model, 1)
    axis([minx maxx miny maxy minz maxz])
    xlabel('x [mm]','FontSize',10)
    set(gca, 'XTick', [-4,0,4],'FontSize',10);
    ylabel('y [mm]','FontSize',10)
    set(gca, 'YTick', [-4,0,4],'FontSize',10);
    zlabel('z [mm]','FontSize',10)
    set(gca, 'ZTick', [-4,0,4],'FontSize',10);
    title('Left view','FontSize',12)
    view([0,1,0])
    subplot(2,3,2); wingkin_plot( [0; 0; 0], Rb_2, kine.RL, kine.RR, a_fit.down_up, body_model, wing_model, 2)
    axis([minx maxx miny maxy minz maxz])
    xlabel('x [mm]','FontSize',10)
    set(gca, 'XTick', [-4,0,4],'FontSize',10);
    ylabel('y [mm]','FontSize',10)
    set(gca, 'YTick', [-4,0,4],'FontSize',10);
    zlabel('z [mm]','FontSize',10)
    set(gca, 'ZTick', [-4,0,4],'FontSize',10);
    title('Right view','FontSize',12)
    view([0,-1,0])
    subplot(2,3,4); wingkin_plot( [0; 0; 0], Rb_3, kine.RL, kine.RR, a_fit.down_up, body_model, wing_model, 0)
    axis([minx maxx miny maxy minz maxz])
    xlabel('x [mm]','FontSize',10)
    set(gca, 'XTick', [-4,0,4],'FontSize',10);
    ylabel('y [mm]','FontSize',10)
    set(gca, 'YTick', [-4,0,4],'FontSize',10);
    zlabel('z [mm]','FontSize',10)
    set(gca, 'ZTick', [-4,0,4],'FontSize',10);
    title('Front view','FontSize',12)
    view([1,0,0])
    subplot(2,3,5); wingkin_plot( [0; 0; 0], Rb_4, kine.RL, kine.RR, a_fit.down_up, body_model, wing_model, 0)
    axis([minx maxx miny maxy minz maxz])
    xlabel('x [mm]','FontSize',10)
    set(gca, 'XTick', [-4,0,4],'FontSize',10);
    ylabel('y [mm]','FontSize',10)
    set(gca, 'YTick', [-4,0,4],'FontSize',10);
    zlabel('z [mm]','FontSize',10)
    set(gca, 'ZTick', [-4,0,4],'FontSize',10);
    title('Top view','FontSize',12)
    view([0,0,1])
    subplot(2,3,3); Aero_force_plot( [0; 0; 0], Rb_1, kine.RL, kine.RR, a_fit.down_up, FA_L, FA_R, body_model, wing_model, 1)
    axis([minx maxx miny maxy minz maxz])
    xlabel('x [mm]','FontSize',10)
    set(gca, 'XTick', [-4,0,4],'FontSize',10);
    ylabel('y [mm]','FontSize',10)
    set(gca, 'YTick', [-4,0,4],'FontSize',10);
    zlabel('z [mm]','FontSize',10)
    set(gca, 'ZTick', [-4,0,4],'FontSize',10);
    title('Forces left','FontSize',12)
    view([0,1,0])
    subplot(2,3,6); Aero_force_plot( [0; 0; 0], Rb_2, kine.RL, kine.RR, a_fit.down_up, FA_L, FA_R, body_model, wing_model, 2)
    axis([minx maxx miny maxy minz maxz])
    xlabel('x [mm]','FontSize',10)
    set(gca, 'XTick', [-4,0,4],'FontSize',10);
    ylabel('y [mm]','FontSize',10)
    set(gca, 'YTick', [-4,0,4],'FontSize',10);
    zlabel('z [mm]','FontSize',10)
    set(gca, 'ZTick', [-4,0,4],'FontSize',10);
    title('Forces right','FontSize',12)
    view([0,-1,0])
    hold off
    
end