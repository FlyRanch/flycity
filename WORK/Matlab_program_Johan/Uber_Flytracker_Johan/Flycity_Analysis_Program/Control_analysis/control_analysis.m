function control_analysis(settings, pathDB)

    savefile = 'pathDB6.mat';
    
    savefile2 = 'FM_man_wb.mat';
    
    rot_lift_on = settings.rot_lift_on;

    % Substract the maneuvering wing kinematics from maneuvering wingbeats:
    
    wb_names_ax = fieldnames(pathDB.maneuver.ax);
    wb_names_ay = fieldnames(pathDB.maneuver.ay);
    wb_names_az = fieldnames(pathDB.maneuver.az);
    wb_names_wx = fieldnames(pathDB.maneuver.wx);
    wb_names_wy = fieldnames(pathDB.maneuver.wy);
    wb_names_wz = fieldnames(pathDB.maneuver.wz);
    
    nr_wb_ax = length(wb_names_ax)
    nr_wb_ay = length(wb_names_ay)
    nr_wb_az = length(wb_names_az)
    nr_wb_wx = length(wb_names_wx)
    nr_wb_wy = length(wb_names_wy)
    nr_wb_wz = length(wb_names_wz)
    
    n_pol_theta = settings.n_pol_theta;
    n_pol_eta   = settings.n_pol_eta;
    n_pol_phi   = settings.n_pol_phi;
    
    down_up_glob = pathDB.poly_fit.a_glob.down_up;
    
    R_strk      = pathDB.rot_mat.Rstr;
    
    % Set-up the maneuver structures:
    
    maneuver.ax = pathDB.maneuver.ax;
    maneuver.ay = pathDB.maneuver.ay;
    maneuver.az = pathDB.maneuver.az;
    maneuver.wx = pathDB.maneuver.wx;
    maneuver.wy = pathDB.maneuver.wy;
    maneuver.wz = pathDB.maneuver.wz;
    
    % Compute non-dimensional aerodynamic forces without body motion:
    
    'man_ax'
    man_ax = FM_aero_maneuver(maneuver.ax,R_strk,wb_names_ax,nr_wb_ax,n_pol_theta,n_pol_eta,n_pol_phi,rot_lift_on);
    'man_ay'
    man_ay = FM_aero_maneuver(maneuver.ay,R_strk,wb_names_ay,nr_wb_ay,n_pol_theta,n_pol_eta,n_pol_phi,rot_lift_on);
    'man_az'
    man_az = FM_aero_maneuver(maneuver.az,R_strk,wb_names_az,nr_wb_az,n_pol_theta,n_pol_eta,n_pol_phi,rot_lift_on);
    'man_wx'
    man_wx = FM_aero_maneuver(maneuver.wx,R_strk,wb_names_wx,nr_wb_wx,n_pol_theta,n_pol_eta,n_pol_phi,rot_lift_on);
    'man_wy'
    man_wy = FM_aero_maneuver(maneuver.wy,R_strk,wb_names_wy,nr_wb_wy,n_pol_theta,n_pol_eta,n_pol_phi,rot_lift_on);
    'man_wz'
    man_wz = FM_aero_maneuver(maneuver.wz,R_strk,wb_names_wz,nr_wb_wz,n_pol_theta,n_pol_eta,n_pol_phi,rot_lift_on);

    save(savefile2,'man_ax','man_ay','man_az','man_wx','man_wy','man_wz')
    
%     man_t = load('FM_man_wb.mat');
%     
%     man_ax = man_t.man_ax;
%     man_ay = man_t.man_ay;
%     man_az = man_t.man_az;
%     man_wx = man_t.man_wx;
%     man_wy = man_t.man_wy;
%     man_wz = man_t.man_wz;
        
    % Motion regression:
    
    c_fit_ax = symmetric_motion_regression( man_ax, down_up_glob );
    c_fit_ay = asymmetric_motion_regression( man_ay, down_up_glob );
    c_fit_az = symmetric_motion_regression( man_az, down_up_glob );
    c_fit_wx = asymmetric_motion_regression( man_wx, down_up_glob );
    c_fit_wy = symmetric_motion_regression( man_wy, down_up_glob );
    c_fit_wz = asymmetric_motion_regression( man_wz, down_up_glob );
    
    % Save the data:
    
    save(savefile,'c_fit_ax','c_fit_ay','c_fit_az','c_fit_wx','c_fit_wy','c_fit_wz')

end

