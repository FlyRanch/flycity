function DynSim_chapter_plots( settings, pathDB, seq_nr )

    % Create plots of the aerodynamic model:

    R_strk = pathDB.rot_mat.Rstr;
    
    nr_points = 100;
    
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

    % Create the wing kinematics

%     a_fit.a_theta_L     = pathDB.poly_fit.a_avg.theta_LR(:,seq_nr);
%     a_fit.a_eta_L       = pathDB.poly_fit.a_avg.eta_LR(:,seq_nr);
%     a_fit.a_phi_L       = pathDB.poly_fit.a_avg.phi_LR(:,seq_nr);
%     a_fit.a_theta_R     = pathDB.poly_fit.a_avg.theta_LR(:,seq_nr);
%     a_fit.a_eta_R       = pathDB.poly_fit.a_avg.eta_LR(:,seq_nr);
%     a_fit.a_phi_R       = pathDB.poly_fit.a_avg.phi_LR(:,seq_nr);
%     a_fit.f             = pathDB.poly_fit.a_avg.f(seq_nr);
%     a_fit.down_up       = pathDB.poly_fit.a_avg.down_up(seq_nr);
%     a_fit.nr_points     = nr_points;
%     a_fit.R_strk        = R_strk;
%         
%     [ kine ] = angular_velocities_polynomial( a_fit );

    a_fit.a_theta_L     = pathDB.poly_fit.a_glob.theta;
    a_fit.a_eta_L       = pathDB.poly_fit.a_glob.eta;
    a_fit.a_phi_L       = pathDB.poly_fit.a_glob.phi;
    a_fit.a_theta_R     = pathDB.poly_fit.a_glob.theta;
    a_fit.a_eta_R       = pathDB.poly_fit.a_glob.eta;
    a_fit.a_phi_R       = pathDB.poly_fit.a_glob.phi;
    a_fit.f             = pathDB.poly_fit.a_glob.f;
    a_fit.down_up       = pathDB.poly_fit.a_glob.down_up;
    a_fit.nr_points     = nr_points;
    a_fit.R_strk        = R_strk;
        
    [ kine ] = angular_velocities_polynomial( a_fit );
    
    
    kine.R_strk = R_strk;
    kine.u_strk = zeros(3,nr_points);
    kine.w_strk = zeros(3,nr_points);
    
    
    % Compute aerodynamic forces:
    
    [ FM_strkpln, FM_L, FM_R , U_L, U_R, alfa_L, alfa_R, alfa_dot_L, alfa_dot_R ] = Aerodynamic_forces( kine, body_model, wing_model, 1);

    t = kine.t(1:nr_points);
    
    FA_L = FM_L(1:3,:);
    FA_R = FM_R(1:3,:);
    
    FA_strk = FM_strkpln(1:3,:);
    MA_strk = FM_strkpln(4:6,:);
    
    figure(1)
    hFig = figure(1);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 1200 800]);
    hold on
    subplot(3,2,1); plot(t,FA_strk(1,:),'r',t,zeros(1,nr_points),'k')
    ylabel('F_x [N]')
    subplot(3,2,3); plot(t,-FA_strk(2,:),'r',t,zeros(1,nr_points),'k')
    ylabel('F_y [N]')
    subplot(3,2,5); plot(t,-FA_strk(3,:),'r',t,zeros(1,nr_points),'k')
    xlabel('t [s]')
    ylabel('F_z [N]')
    subplot(3,2,2); plot(t,MA_strk(1,:),'r',t,zeros(1,nr_points),'k')
    ylabel('M_x [N*mm]')
    subplot(3,2,4); plot(t,-MA_strk(2,:),'r',t,zeros(1,nr_points),'k')
    ylabel('M_y [N*mm]')
    subplot(3,2,6); plot(t,-MA_strk(3,:),'r',t,zeros(1,nr_points),'k')
    xlabel('t [s]')
    ylabel('M_z [N*mm]')
    hold off
    [~,h1] = suplabel('Quasi-steady aerodynamic model, forces and moments', 't');
    set(h1,'FontSize',10)
    saveas(hFig,[char(settings.plot_loc) char('/Dynamic_simulation_plots/FM_quasi_steady')],'fig')
    
    figure(2)
    hFig = figure(2);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 800 800]);
    hold on
    subplot(1,2,1); plot(t,(180/pi)*alfa_L(1,:),'b',t,(180/pi)*abs(alfa_L(1,:)),'r',t,zeros(1,nr_points),'k')
    title('\alpha left wing')
    xlabel('t [s]')
    ylabel('\alpha [deg]')
    subplot(1,2,2);plot(t,(180/pi)*alfa_dot_L(1,:),'b',t,zeros(1,nr_points),'k')
    title('d/dt(\alpha) left wing')
    xlabel('t [s]')
    ylabel('d/dt(\alpha) [deg/s]')
    hold off
    [~,h1] = suplabel('Angle-of-attack and derivative', 't');
    set(h1,'FontSize',10)
    saveas(hFig,[char(settings.plot_loc) char('/Dynamic_simulation_plots/alpha_and_alpha_dot')],'fig')
    
    Rb = R_strk'*[1 0 0; 0 -1 0; 0 0 -1];
    
    minx = -4;
    maxx = 4;
    miny = -4;
    maxy = 4;
    minz = -4;
    maxz = 4;
    
    figure(3)
    hFig = figure(3);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 800 800]);
    hold on
    subplot(2,2,1); Aero_force_plot( [0; 0; 0], Rb, kine.RL(:,:,1:4:end), kine.RR(:,:,1:4:end), a_fit.down_up, FA_L(:,1:4:end), FA_R(:,1:4:end), body_model, wing_model, 1)
    axis([minx maxx miny maxy minz maxz])
    xlabel('x [mm]','FontSize',10)
    set(gca, 'XTick', [-4,0,4],'FontSize',10);
    ylabel('y [mm]','FontSize',10)
    set(gca, 'YTick', [-4,0,4],'FontSize',10);
    zlabel('z [mm]','FontSize',10)
    set(gca, 'ZTick', [-4,0,4],'FontSize',10);
    title('Left view','FontSize',12)
    view([0,1,0])
    subplot(2,2,2); Aero_force_plot( [0; 0; 0], Rb, kine.RL(:,:,1:4:end), kine.RR(:,:,1:4:end), a_fit.down_up, FA_L(:,1:4:end), FA_R(:,1:4:end), body_model, wing_model, 2)
    axis([minx maxx miny maxy minz maxz])
    xlabel('x [mm]','FontSize',10)
    set(gca, 'XTick', [-4,0,4],'FontSize',10);
    ylabel('y [mm]','FontSize',10)
    set(gca, 'YTick', [-4,0,4],'FontSize',10);
    zlabel('z [mm]','FontSize',10)
    set(gca, 'ZTick', [-4,0,4],'FontSize',10);
    title('Right view','FontSize',12)
    view([0,-1,0])
    subplot(2,2,3); Aero_force_plot( [0; 0; 0], Rb, kine.RL(:,:,1:4:end), kine.RR(:,:,1:4:end), a_fit.down_up, FA_L(:,1:4:end), FA_R(:,1:4:end), body_model, wing_model, 0)
    axis([minx maxx miny maxy minz maxz])
    xlabel('x [mm]','FontSize',10)
    set(gca, 'XTick', [-4,0,4],'FontSize',10);
    ylabel('y [mm]','FontSize',10)
    set(gca, 'YTick', [-4,0,4],'FontSize',10);
    zlabel('z [mm]','FontSize',10)
    set(gca, 'ZTick', [-4,0,4],'FontSize',10);
    title('Front view','FontSize',12)
    view([1,0,0])
    subplot(2,2,4); Aero_force_plot( [0; 0; 0], Rb, kine.RL(:,:,1:4:end), kine.RR(:,:,1:4:end), a_fit.down_up, FA_L(:,1:4:end), FA_R(:,1:4:end), body_model, wing_model, 0)
    axis([minx maxx miny maxy minz maxz])
    xlabel('x [mm]','FontSize',10)
    set(gca, 'XTick', [-4,0,4],'FontSize',10);
    ylabel('y [mm]','FontSize',10)
    set(gca, 'YTick', [-4,0,4],'FontSize',10);
    zlabel('z [mm]','FontSize',10)
    set(gca, 'ZTick', [-4,0,4],'FontSize',10);
    title('Top view','FontSize',12)
    view([0,0,1])
    hold off
    [~,h1] = suplabel('Quasi-steady aerodynamic model, forces and wing kinematics', 't');
    set(h1,'FontSize',10)
    saveas(hFig,[char(settings.plot_loc) char('/Dynamic_simulation_plots/Force_wingkin_quasi_steady')],'fig')
    
end

