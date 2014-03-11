function flight_traject( settings, pathDB, seq_nr )

    % Plot the flight trajectory and force plots at critical points:
    
    nr_of_wb = pathDB.wingbeats.nr_of_wb(seq_nr);
    
    wb_vect = [ 1:1:nr_of_wb ];
    
    N = length(wb_vect);
    
    xyz_path = pathDB.filt.xyz(:,:,seq_nr);
    
    qb       = pathDB.filt.qB(:,:,seq_nr);
    
    xyz_wb   = zeros(3,N);
    Rb       = zeros(3,3,N);
    a_xyz    = zeros(3,N);
    
    for i = 1:N
        
        wb_range = pathDB.wingbeats.wingbeat_loc(wb_vect(i),1,seq_nr):pathDB.wingbeats.wingbeat_loc(wb_vect(i),2,seq_nr);
        
        xyz_wb(:,i) = mean(xyz_path(wb_range,:))';
        
        [ qb_avg ] = q_avg(qb(wb_range,1), qb(wb_range,2) ,qb(wb_range,3), qb(wb_range,4));
        Rb_avg = quat2mat(qb_avg);
        
        Rb(:,:,i) = Rb_avg;
        
        a_xyz(:,i) = mean(pathDB.filt.a_xyz(wb_range,:,seq_nr))';
        
    end
    
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
    
    
    x_min = min(xyz_path(:,1));
    x_max = max(xyz_path(:,1));
    y_min = min(xyz_path(:,2));
    y_max = max(xyz_path(:,2));
    z_min = min(xyz_path(:,3));
    z_max = max(xyz_path(:,3));

    hold on
    plot3(xyz_path(:,1),xyz_path(:,2),xyz_path(:,3),'Color','r')
    for i = 1:N
        Fly_helicopter_model( Rb(:,:,i), R_strk, xyz_wb(:,i), a_xyz(:,i), body_model, wing_model, 0.5 )
    end
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    axis equal
    grid on
    xlim([x_min-4 x_max+4])
    ylim([y_min-4 y_max+4])
    zlim([z_min-4 z_max+4])
    hold off


end

