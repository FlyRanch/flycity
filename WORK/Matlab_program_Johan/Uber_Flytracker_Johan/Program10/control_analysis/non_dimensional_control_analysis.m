function [all_maneuvers] = non_dimensional_control_analysis( settings, pathDB)


    nr_of_seq = size(pathDB.x,2);
    
    dt = pathDB.t(2)-pathDB.t(1);
    
    all_maneuvers = {};
    
    roll_maneuvers = {};
    pitch_maneuvers = {};
    yaw_maneuvers = {};
    ax_maneuvers = {};
    ay_maneuvers = {};
    az_maneuvers = {};
    pure_roll = {};
    pure_pitch = {};
    pure_yaw = {};
    pure_ax = {};
    pure_ay = {};
    pure_az = {};
    
    for i = 1:nr_of_seq
        
        maneuver = pathDB.maneuver.(char(['maneuver_' int2str(i)]));
        
        stroke_var = pathDB.stroke_var.(char(['stroke_var_' int2str(i)]));
        
        if isstruct(maneuver) == 1

            L_wing = pathDB.wing_l(i);

            mass_fly = 1.85e-6*(L_wing/3)^3;

%             [cg_body, cg_L, cg_R, I_body, I_wing, I_v_wing, m_b, m_w, m_v_w ,y_sect, chords] = cg_plus_Inertia(settings, pathDB, mass_fly, L_wing , i, 20 );
            
            nr_wb = length(maneuver.yaw_turns);
            
            roll_man_nr = 0;
            pitch_man_nr = 0;
            yaw_man_nr = 0;
            ax_man_nr = 0;
            ay_man_nr = 0;
            az_man_nr = 0;
            pure_roll_nr = 0; 
            pure_pitch_nr = 0;
            pure_yaw_nr = 0;
            pure_ax_nr = 0;
            pure_ay_nr = 0;
            pure_az_nr = 0;
            
            for j = 1:nr_wb
                
                if maneuver.yaw_turns(j) ~= 0
                    
                    if maneuver.yaw_turns(j-1) == 0
                        
                        yaw_man_nr = yaw_man_nr+1;
                        
                        temp_man = {};
                        
                        man_end = j+find(maneuver.yaw_turns(j:end)==0,1,'first')-2;
                        
                        movie_percentage = j/nr_wb;
                                                
                        temp_man.seq_nr = i;
                        temp_man.seq_name = settings.sequence_names(i);
                        temp_man.movie_percentage = movie_percentage;
                        temp_man.wing_l = L_wing;
                        temp_man.mass_fly = mass_fly;
                        temp_man.a_avg = pathDB.a_avg.(char(['a_avg_' int2str(i) ]));
                        temp_man.f_avg = pathDB.f_avg.(char(['f_avg_' int2str(i) ]));
                        temp_man.down_up_avg = pathDB.down_up_avg.(char(['down_up_avg_' int2str(i) ]));
                                                
                        man_a_dev = pathDB.a_dev.(char(['a_dev_' int2str(i) ]));
                        man_a_fit = pathDB.a_fit.(char(['a_fit_' int2str(i) ]));
                        man_down_up_ratio = pathDB.down_up_ratio.(char(['down_up_ratio_' int2str(i) ]));
                                                
                        temp_man.a_dev.theta_L1 = man_a_dev.theta_L1(:,j:man_end);
                        temp_man.a_dev.theta_L2 = man_a_dev.theta_L2(:,j:man_end);
                        temp_man.a_dev.eta_L1 = man_a_dev.eta_L1(:,j:man_end);
                        temp_man.a_dev.eta_L2 = man_a_dev.eta_L2(:,j:man_end);
                        temp_man.a_dev.phi_L1 = man_a_dev.phi_L1(:,j:man_end);
                        temp_man.a_dev.phi_L2 = man_a_dev.phi_L2(:,j:man_end);
                        
                        temp_man.a_dev.theta_R1 = man_a_dev.theta_R1(:,j:man_end);
                        temp_man.a_dev.theta_R2 = man_a_dev.theta_R2(:,j:man_end);
                        temp_man.a_dev.eta_R1 = man_a_dev.eta_R1(:,j:man_end);
                        temp_man.a_dev.eta_R2 = man_a_dev.eta_R2(:,j:man_end);
                        temp_man.a_dev.phi_R1 = man_a_dev.phi_R1(:,j:man_end);
                        temp_man.a_dev.phi_R2 = man_a_dev.phi_R2(:,j:man_end);
                        
                        temp_man.q_body_wb = stroke_var.q_body(:,j:man_end);
                        temp_man.Omega_strk = stroke_var.Omega_strk(:,j:man_end);
                        temp_man.Omega_dot_strk = stroke_var.Omega_dot_strk(:,j:man_end);
                        temp_man.V_strk = stroke_var.V_strk(:,j:man_end);
                        temp_man.A_strk = stroke_var.A_strk(:,j:man_end);
                            
                        temp_man.f = man_a_fit.f(j:man_end);
                        temp_man.down_up = man_down_up_ratio(j:man_end);
                        
                        yaw_maneuvers.(char(['man_' int2str(i) '_' int2str(yaw_man_nr) ])) = temp_man;
                        
                    end
                    
                end

                
                if maneuver.roll_turns(j) ~= 0
                    
                    if maneuver.roll_turns(j-1) == 0
                        
                        roll_man_nr = roll_man_nr+1;
                        
                        temp_man = {};
                        
                        man_end = j+find(maneuver.roll_turns(j:end)==0,1,'first')-2;
                        
                        movie_percentage = j/nr_wb;
                                                
                        temp_man.seq_nr = i;
                        temp_man.seq_name = settings.sequence_names(i);
                        temp_man.movie_percentage = movie_percentage;
                        temp_man.wing_l = L_wing;
                        temp_man.mass_fly = mass_fly;
                        temp_man.a_avg = pathDB.a_avg.(char(['a_avg_' int2str(i) ]));
                        temp_man.f_avg = pathDB.f_avg.(char(['f_avg_' int2str(i) ]));
                        temp_man.down_up_avg = pathDB.down_up_avg.(char(['down_up_avg_' int2str(i) ]));
                                                
                        man_a_dev = pathDB.a_dev.(char(['a_dev_' int2str(i) ]));
                        man_a_fit = pathDB.a_fit.(char(['a_fit_' int2str(i) ]));
                        man_down_up_ratio = pathDB.down_up_ratio.(char(['down_up_ratio_' int2str(i) ]));
                                                
                        temp_man.a_dev.theta_L1 = man_a_dev.theta_L1(:,j:man_end);
                        temp_man.a_dev.theta_L2 = man_a_dev.theta_L2(:,j:man_end);
                        temp_man.a_dev.eta_L1 = man_a_dev.eta_L1(:,j:man_end);
                        temp_man.a_dev.eta_L2 = man_a_dev.eta_L2(:,j:man_end);
                        temp_man.a_dev.phi_L1 = man_a_dev.phi_L1(:,j:man_end);
                        temp_man.a_dev.phi_L2 = man_a_dev.phi_L2(:,j:man_end);
                        
                        temp_man.a_dev.theta_R1 = man_a_dev.theta_R1(:,j:man_end);
                        temp_man.a_dev.theta_R2 = man_a_dev.theta_R2(:,j:man_end);
                        temp_man.a_dev.eta_R1 = man_a_dev.eta_R1(:,j:man_end);
                        temp_man.a_dev.eta_R2 = man_a_dev.eta_R2(:,j:man_end);
                        temp_man.a_dev.phi_R1 = man_a_dev.phi_R1(:,j:man_end);
                        temp_man.a_dev.phi_R2 = man_a_dev.phi_R2(:,j:man_end);
                        
                        temp_man.q_body_wb = stroke_var.q_body(:,j:man_end);
                        temp_man.Omega_strk = stroke_var.Omega_strk(:,j:man_end);
                        temp_man.Omega_dot_strk = stroke_var.Omega_dot_strk(:,j:man_end);
                        temp_man.V_strk = stroke_var.V_strk(:,j:man_end);
                        temp_man.A_strk = stroke_var.A_strk(:,j:man_end);
                            
                        temp_man.f = man_a_fit.f(j:man_end);
                        temp_man.down_up = man_down_up_ratio(j:man_end);
                        
                        roll_maneuvers.(char(['man_' int2str(i) '_' int2str(roll_man_nr) ])) = temp_man;
                        
                    end
                    
                end
                
                if maneuver.pitch_turns(j) ~= 0
                    
                    if maneuver.pitch_turns(j-1) == 0
                        
                        pitch_man_nr = pitch_man_nr+1;
                        
                        temp_man = {};
                        
                        man_end = j+find(maneuver.pitch_turns(j:end)==0,1,'first')-2;
                        
                        movie_percentage = j/nr_wb;
                                                
                        temp_man.seq_nr = i;
                        temp_man.seq_name = settings.sequence_names(i);
                        temp_man.movie_percentage = movie_percentage;
                        temp_man.wing_l = L_wing;
                        temp_man.mass_fly = mass_fly;
                        temp_man.a_avg = pathDB.a_avg.(char(['a_avg_' int2str(i) ]));
                        temp_man.f_avg = pathDB.f_avg.(char(['f_avg_' int2str(i) ]));
                        temp_man.down_up_avg = pathDB.down_up_avg.(char(['down_up_avg_' int2str(i) ]));
                                                
                        man_a_dev = pathDB.a_dev.(char(['a_dev_' int2str(i) ]));
                        man_a_fit = pathDB.a_fit.(char(['a_fit_' int2str(i) ]));
                        man_down_up_ratio = pathDB.down_up_ratio.(char(['down_up_ratio_' int2str(i) ]));
                                                
                        temp_man.a_dev.theta_L1 = man_a_dev.theta_L1(:,j:man_end);
                        temp_man.a_dev.theta_L2 = man_a_dev.theta_L2(:,j:man_end);
                        temp_man.a_dev.eta_L1 = man_a_dev.eta_L1(:,j:man_end);
                        temp_man.a_dev.eta_L2 = man_a_dev.eta_L2(:,j:man_end);
                        temp_man.a_dev.phi_L1 = man_a_dev.phi_L1(:,j:man_end);
                        temp_man.a_dev.phi_L2 = man_a_dev.phi_L2(:,j:man_end);
                        
                        temp_man.a_dev.theta_R1 = man_a_dev.theta_R1(:,j:man_end);
                        temp_man.a_dev.theta_R2 = man_a_dev.theta_R2(:,j:man_end);
                        temp_man.a_dev.eta_R1 = man_a_dev.eta_R1(:,j:man_end);
                        temp_man.a_dev.eta_R2 = man_a_dev.eta_R2(:,j:man_end);
                        temp_man.a_dev.phi_R1 = man_a_dev.phi_R1(:,j:man_end);
                        temp_man.a_dev.phi_R2 = man_a_dev.phi_R2(:,j:man_end);
                        
                        temp_man.q_body_wb = stroke_var.q_body(:,j:man_end);
                        temp_man.Omega_strk = stroke_var.Omega_strk(:,j:man_end);
                        temp_man.Omega_dot_strk = stroke_var.Omega_dot_strk(:,j:man_end);
                        temp_man.V_strk = stroke_var.V_strk(:,j:man_end);
                        temp_man.A_strk = stroke_var.A_strk(:,j:man_end);
                            
                        temp_man.f = man_a_fit.f(j:man_end);
                        temp_man.down_up = man_down_up_ratio(j:man_end);
                        
                        pitch_maneuvers.(char(['man_' int2str(i) '_' int2str(pitch_man_nr) ])) = temp_man;
                        
                    end
                    
                end                

                
                if maneuver.a_x(j) ~= 0
                    
                    if maneuver.a_x(j-1) == 0
                        
                        ax_man_nr = ax_man_nr+1;
                        
                        temp_man = {};
                        
                        man_end = j+find(maneuver.a_x(j:end)==0,1,'first')-2;
                        
                        movie_percentage = j/nr_wb;
                                                
                        temp_man.seq_nr = i;
                        temp_man.seq_name = settings.sequence_names(i);
                        temp_man.movie_percentage = movie_percentage;
                        temp_man.wing_l = L_wing;
                        temp_man.mass_fly = mass_fly;
                        temp_man.a_avg = pathDB.a_avg.(char(['a_avg_' int2str(i) ]));
                        temp_man.f_avg = pathDB.f_avg.(char(['f_avg_' int2str(i) ]));
                        temp_man.down_up_avg = pathDB.down_up_avg.(char(['down_up_avg_' int2str(i) ]));
                                                
                        man_a_dev = pathDB.a_dev.(char(['a_dev_' int2str(i) ]));
                        man_a_fit = pathDB.a_fit.(char(['a_fit_' int2str(i) ]));
                        man_down_up_ratio = pathDB.down_up_ratio.(char(['down_up_ratio_' int2str(i) ]));
                                                
                        temp_man.a_dev.theta_L1 = man_a_dev.theta_L1(:,j:man_end);
                        temp_man.a_dev.theta_L2 = man_a_dev.theta_L2(:,j:man_end);
                        temp_man.a_dev.eta_L1 = man_a_dev.eta_L1(:,j:man_end);
                        temp_man.a_dev.eta_L2 = man_a_dev.eta_L2(:,j:man_end);
                        temp_man.a_dev.phi_L1 = man_a_dev.phi_L1(:,j:man_end);
                        temp_man.a_dev.phi_L2 = man_a_dev.phi_L2(:,j:man_end);
                        
                        temp_man.a_dev.theta_R1 = man_a_dev.theta_R1(:,j:man_end);
                        temp_man.a_dev.theta_R2 = man_a_dev.theta_R2(:,j:man_end);
                        temp_man.a_dev.eta_R1 = man_a_dev.eta_R1(:,j:man_end);
                        temp_man.a_dev.eta_R2 = man_a_dev.eta_R2(:,j:man_end);
                        temp_man.a_dev.phi_R1 = man_a_dev.phi_R1(:,j:man_end);
                        temp_man.a_dev.phi_R2 = man_a_dev.phi_R2(:,j:man_end);
                        
                        temp_man.q_body_wb = stroke_var.q_body(:,j:man_end);
                        temp_man.Omega_strk = stroke_var.Omega_strk(:,j:man_end);
                        temp_man.Omega_dot_strk = stroke_var.Omega_dot_strk(:,j:man_end);
                        temp_man.V_strk = stroke_var.V_strk(:,j:man_end);
                        temp_man.A_strk = stroke_var.A_strk(:,j:man_end);
                            
                        temp_man.f = man_a_fit.f(j:man_end);
                        temp_man.down_up = man_down_up_ratio(j:man_end);
                        
                        ax_maneuvers.(char(['man_' int2str(i) '_' int2str(ax_man_nr) ])) = temp_man;
                        
                    end
                    
                end     

                
                if maneuver.a_y(j) ~= 0
                    
                    if maneuver.a_y(j-1) == 0
                        
                        ay_man_nr = ay_man_nr+1;
                        
                        temp_man = {};
                        
                        man_end = j+find(maneuver.a_y(j:end)==0,1,'first')-2;
                        
                        movie_percentage = j/nr_wb;
                                                
                        temp_man.seq_nr = i;
                        temp_man.seq_name = settings.sequence_names(i);
                        temp_man.movie_percentage = movie_percentage;
                        temp_man.wing_l = L_wing;
                        temp_man.mass_fly = mass_fly;
                        temp_man.a_avg = pathDB.a_avg.(char(['a_avg_' int2str(i) ]));
                        temp_man.f_avg = pathDB.f_avg.(char(['f_avg_' int2str(i) ]));
                        temp_man.down_up_avg = pathDB.down_up_avg.(char(['down_up_avg_' int2str(i) ]));
                                                
                        man_a_dev = pathDB.a_dev.(char(['a_dev_' int2str(i) ]));
                        man_a_fit = pathDB.a_fit.(char(['a_fit_' int2str(i) ]));
                        man_down_up_ratio = pathDB.down_up_ratio.(char(['down_up_ratio_' int2str(i) ]));
                                                
                        temp_man.a_dev.theta_L1 = man_a_dev.theta_L1(:,j:man_end);
                        temp_man.a_dev.theta_L2 = man_a_dev.theta_L2(:,j:man_end);
                        temp_man.a_dev.eta_L1 = man_a_dev.eta_L1(:,j:man_end);
                        temp_man.a_dev.eta_L2 = man_a_dev.eta_L2(:,j:man_end);
                        temp_man.a_dev.phi_L1 = man_a_dev.phi_L1(:,j:man_end);
                        temp_man.a_dev.phi_L2 = man_a_dev.phi_L2(:,j:man_end);
                        
                        temp_man.a_dev.theta_R1 = man_a_dev.theta_R1(:,j:man_end);
                        temp_man.a_dev.theta_R2 = man_a_dev.theta_R2(:,j:man_end);
                        temp_man.a_dev.eta_R1 = man_a_dev.eta_R1(:,j:man_end);
                        temp_man.a_dev.eta_R2 = man_a_dev.eta_R2(:,j:man_end);
                        temp_man.a_dev.phi_R1 = man_a_dev.phi_R1(:,j:man_end);
                        temp_man.a_dev.phi_R2 = man_a_dev.phi_R2(:,j:man_end);
                        
                        temp_man.q_body_wb = stroke_var.q_body(:,j:man_end);
                        temp_man.Omega_strk = stroke_var.Omega_strk(:,j:man_end);
                        temp_man.Omega_dot_strk = stroke_var.Omega_dot_strk(:,j:man_end);
                        temp_man.V_strk = stroke_var.V_strk(:,j:man_end);
                        temp_man.A_strk = stroke_var.A_strk(:,j:man_end);
                            
                        temp_man.f = man_a_fit.f(j:man_end);
                        temp_man.down_up = man_down_up_ratio(j:man_end);
                        
                        ay_maneuvers.(char(['man_' int2str(i) '_' int2str(ay_man_nr) ])) = temp_man;
                        
                    end
                    
                end 
                
                if maneuver.a_z(j) ~= 0
                    
                    if maneuver.a_z(j-1) == 0
                        
                        az_man_nr = az_man_nr+1;
                        
                        temp_man = {};
                        
                        man_end = j+find(maneuver.a_z(j:end)==0,1,'first')-2;
                        
                        movie_percentage = j/nr_wb;
                                                
                        temp_man.seq_nr = i;
                        temp_man.seq_name = settings.sequence_names(i);
                        temp_man.movie_percentage = movie_percentage;
                        temp_man.wing_l = L_wing;
                        temp_man.mass_fly = mass_fly;
                        temp_man.a_avg = pathDB.a_avg.(char(['a_avg_' int2str(i) ]));
                        temp_man.f_avg = pathDB.f_avg.(char(['f_avg_' int2str(i) ]));
                        temp_man.down_up_avg = pathDB.down_up_avg.(char(['down_up_avg_' int2str(i) ]));
                                                
                        man_a_dev = pathDB.a_dev.(char(['a_dev_' int2str(i) ]));
                        man_a_fit = pathDB.a_fit.(char(['a_fit_' int2str(i) ]));
                        man_down_up_ratio = pathDB.down_up_ratio.(char(['down_up_ratio_' int2str(i) ]));
                                                
                        temp_man.a_dev.theta_L1 = man_a_dev.theta_L1(:,j:man_end);
                        temp_man.a_dev.theta_L2 = man_a_dev.theta_L2(:,j:man_end);
                        temp_man.a_dev.eta_L1 = man_a_dev.eta_L1(:,j:man_end);
                        temp_man.a_dev.eta_L2 = man_a_dev.eta_L2(:,j:man_end);
                        temp_man.a_dev.phi_L1 = man_a_dev.phi_L1(:,j:man_end);
                        temp_man.a_dev.phi_L2 = man_a_dev.phi_L2(:,j:man_end);
                        
                        temp_man.a_dev.theta_R1 = man_a_dev.theta_R1(:,j:man_end);
                        temp_man.a_dev.theta_R2 = man_a_dev.theta_R2(:,j:man_end);
                        temp_man.a_dev.eta_R1 = man_a_dev.eta_R1(:,j:man_end);
                        temp_man.a_dev.eta_R2 = man_a_dev.eta_R2(:,j:man_end);
                        temp_man.a_dev.phi_R1 = man_a_dev.phi_R1(:,j:man_end);
                        temp_man.a_dev.phi_R2 = man_a_dev.phi_R2(:,j:man_end);
                        
                        temp_man.q_body_wb = stroke_var.q_body(:,j:man_end);
                        temp_man.Omega_strk = stroke_var.Omega_strk(:,j:man_end);
                        temp_man.Omega_dot_strk = stroke_var.Omega_dot_strk(:,j:man_end);
                        temp_man.V_strk = stroke_var.V_strk(:,j:man_end);
                        temp_man.A_strk = stroke_var.A_strk(:,j:man_end);
                            
                        temp_man.f = man_a_fit.f(j:man_end);
                        temp_man.down_up = man_down_up_ratio(j:man_end);
                        
                        az_maneuvers.(char(['man_' int2str(i) '_' int2str(az_man_nr) ])) = temp_man;
                        
                    end
                    
                end
                
                
                if maneuver.yaw_turns(j) ~= 0
                    
                    if maneuver.yaw_turns(j-1) == 0
                        
                        man_end = j+find(maneuver.yaw_turns(j:end)==0,1,'first')-2;
                        
                        pitch_man_on = sum(maneuver.pitch_turns(j:man_end));
                        
                        roll_man_on = sum(maneuver.roll_turns(j:man_end));
                        
                        if pitch_man_on == 0 && roll_man_on == 0
                        
                        pure_yaw_nr = pure_yaw_nr+1;
                        
                        temp_man = {};
                        
                        movie_percentage = j/nr_wb;
                        
                        temp_man.seq_nr = i;
                        temp_man.seq_name = settings.sequence_names(i);
                        temp_man.movie_percentage = movie_percentage;
                        temp_man.wing_l = L_wing;
                        temp_man.mass_fly = mass_fly;
                        temp_man.a_avg = pathDB.a_avg.(char(['a_avg_' int2str(i) ]));
                        temp_man.f_avg = pathDB.f_avg.(char(['f_avg_' int2str(i) ]));
                        temp_man.down_up_avg = pathDB.down_up_avg.(char(['down_up_avg_' int2str(i) ]));
                                                
                        man_a_dev = pathDB.a_dev.(char(['a_dev_' int2str(i) ]));
                        man_a_fit = pathDB.a_fit.(char(['a_fit_' int2str(i) ]));
                        man_down_up_ratio = pathDB.down_up_ratio.(char(['down_up_ratio_' int2str(i) ]));
                                                
                        temp_man.a_dev.theta_L1 = man_a_dev.theta_L1(:,j:man_end);
                        temp_man.a_dev.theta_L2 = man_a_dev.theta_L2(:,j:man_end);
                        temp_man.a_dev.eta_L1 = man_a_dev.eta_L1(:,j:man_end);
                        temp_man.a_dev.eta_L2 = man_a_dev.eta_L2(:,j:man_end);
                        temp_man.a_dev.phi_L1 = man_a_dev.phi_L1(:,j:man_end);
                        temp_man.a_dev.phi_L2 = man_a_dev.phi_L2(:,j:man_end);
                        
                        temp_man.a_dev.theta_R1 = man_a_dev.theta_R1(:,j:man_end);
                        temp_man.a_dev.theta_R2 = man_a_dev.theta_R2(:,j:man_end);
                        temp_man.a_dev.eta_R1 = man_a_dev.eta_R1(:,j:man_end);
                        temp_man.a_dev.eta_R2 = man_a_dev.eta_R2(:,j:man_end);
                        temp_man.a_dev.phi_R1 = man_a_dev.phi_R1(:,j:man_end);
                        temp_man.a_dev.phi_R2 = man_a_dev.phi_R2(:,j:man_end);
                        
                        temp_man.q_body_wb = stroke_var.q_body(:,j:man_end);
                        temp_man.Omega_strk = stroke_var.Omega_strk(:,j:man_end);
                        temp_man.Omega_dot_strk = stroke_var.Omega_dot_strk(:,j:man_end);
                        temp_man.V_strk = stroke_var.V_strk(:,j:man_end);
                        temp_man.A_strk = stroke_var.A_strk(:,j:man_end);
                            
                        temp_man.f = man_a_fit.f(j:man_end);
                        temp_man.down_up = man_down_up_ratio(j:man_end);
                        
                        pure_yaw.(char(['man_' int2str(i) '_' int2str(pure_yaw_nr) ])) = temp_man;
                        
                        end
                        
                    end
                    
                end

                if maneuver.pitch_turns(j) ~= 0
                    
                    if maneuver.pitch_turns(j-1) == 0
                        
                        man_end = j+find(maneuver.pitch_turns(j:end)==0,1,'first')-2;
                        
                        yaw_man_on = sum(maneuver.yaw_turns(j:man_end));
                        
                        roll_man_on = sum(maneuver.roll_turns(j:man_end));
                        
                        if yaw_man_on == 0 && roll_man_on == 0
                        
                        pure_pitch_nr = pure_pitch_nr+1;
                        
                        temp_man = {};
                        
                        movie_percentage = j/nr_wb;
                        
                        temp_man.seq_nr = i;
                        temp_man.seq_name = settings.sequence_names(i);
                        temp_man.movie_percentage = movie_percentage;
                        temp_man.wing_l = L_wing;
                        temp_man.mass_fly = mass_fly;
                        temp_man.a_avg = pathDB.a_avg.(char(['a_avg_' int2str(i) ]));
                        temp_man.f_avg = pathDB.f_avg.(char(['f_avg_' int2str(i) ]));
                        temp_man.down_up_avg = pathDB.down_up_avg.(char(['down_up_avg_' int2str(i) ]));
                                                
                        man_a_dev = pathDB.a_dev.(char(['a_dev_' int2str(i) ]));
                        man_a_fit = pathDB.a_fit.(char(['a_fit_' int2str(i) ]));
                        man_down_up_ratio = pathDB.down_up_ratio.(char(['down_up_ratio_' int2str(i) ]));
                                                
                        temp_man.a_dev.theta_L1 = man_a_dev.theta_L1(:,j:man_end);
                        temp_man.a_dev.theta_L2 = man_a_dev.theta_L2(:,j:man_end);
                        temp_man.a_dev.eta_L1 = man_a_dev.eta_L1(:,j:man_end);
                        temp_man.a_dev.eta_L2 = man_a_dev.eta_L2(:,j:man_end);
                        temp_man.a_dev.phi_L1 = man_a_dev.phi_L1(:,j:man_end);
                        temp_man.a_dev.phi_L2 = man_a_dev.phi_L2(:,j:man_end);
                        
                        temp_man.a_dev.theta_R1 = man_a_dev.theta_R1(:,j:man_end);
                        temp_man.a_dev.theta_R2 = man_a_dev.theta_R2(:,j:man_end);
                        temp_man.a_dev.eta_R1 = man_a_dev.eta_R1(:,j:man_end);
                        temp_man.a_dev.eta_R2 = man_a_dev.eta_R2(:,j:man_end);
                        temp_man.a_dev.phi_R1 = man_a_dev.phi_R1(:,j:man_end);
                        temp_man.a_dev.phi_R2 = man_a_dev.phi_R2(:,j:man_end);
                        
                        temp_man.q_body_wb = stroke_var.q_body(:,j:man_end);
                        temp_man.Omega_strk = stroke_var.Omega_strk(:,j:man_end);
                        temp_man.Omega_dot_strk = stroke_var.Omega_dot_strk(:,j:man_end);
                        temp_man.V_strk = stroke_var.V_strk(:,j:man_end);
                        temp_man.A_strk = stroke_var.A_strk(:,j:man_end);
                            
                        temp_man.f = man_a_fit.f(j:man_end);
                        temp_man.down_up = man_down_up_ratio(j:man_end);
                        
                        pure_pitch.(char(['man_' int2str(i) '_' int2str(pure_pitch_nr) ])) = temp_man;
                        
                        end
                        
                    end
                    
                end
                
                
                
                
                if maneuver.roll_turns(j) ~= 0
                   
                    if maneuver.roll_turns(j-1) == 0
                        
                        man_end = j+find(maneuver.roll_turns(j:end)==0,1,'first')-2;
                        
                        yaw_man_on = sum(maneuver.yaw_turns(j:man_end));
                        
                        pitch_man_on = sum(maneuver.pitch_turns(j:man_end));
                        
                        if yaw_man_on == 0 && pitch_man_on == 0
                        
                        pure_roll_nr = pure_roll_nr+1;
                        
                        temp_man = {};
                        
                        movie_percentage = j/nr_wb;
                        
                        temp_man.seq_nr = i;
                        temp_man.seq_name = settings.sequence_names(i);
                        temp_man.movie_percentage = movie_percentage;
                        temp_man.wing_l = L_wing;
                        temp_man.mass_fly = mass_fly;
                        temp_man.a_avg = pathDB.a_avg.(char(['a_avg_' int2str(i) ]));
                        temp_man.f_avg = pathDB.f_avg.(char(['f_avg_' int2str(i) ]));
                        temp_man.down_up_avg = pathDB.down_up_avg.(char(['down_up_avg_' int2str(i) ]));
                                                
                        man_a_dev = pathDB.a_dev.(char(['a_dev_' int2str(i) ]));
                        man_a_fit = pathDB.a_fit.(char(['a_fit_' int2str(i) ]));
                        man_down_up_ratio = pathDB.down_up_ratio.(char(['down_up_ratio_' int2str(i) ]));
                                                
                        temp_man.a_dev.theta_L1 = man_a_dev.theta_L1(:,j:man_end);
                        temp_man.a_dev.theta_L2 = man_a_dev.theta_L2(:,j:man_end);
                        temp_man.a_dev.eta_L1 = man_a_dev.eta_L1(:,j:man_end);
                        temp_man.a_dev.eta_L2 = man_a_dev.eta_L2(:,j:man_end);
                        temp_man.a_dev.phi_L1 = man_a_dev.phi_L1(:,j:man_end);
                        temp_man.a_dev.phi_L2 = man_a_dev.phi_L2(:,j:man_end);
                        
                        temp_man.a_dev.theta_R1 = man_a_dev.theta_R1(:,j:man_end);
                        temp_man.a_dev.theta_R2 = man_a_dev.theta_R2(:,j:man_end);
                        temp_man.a_dev.eta_R1 = man_a_dev.eta_R1(:,j:man_end);
                        temp_man.a_dev.eta_R2 = man_a_dev.eta_R2(:,j:man_end);
                        temp_man.a_dev.phi_R1 = man_a_dev.phi_R1(:,j:man_end);
                        temp_man.a_dev.phi_R2 = man_a_dev.phi_R2(:,j:man_end);
                        
                        temp_man.q_body_wb = stroke_var.q_body(:,j:man_end);
                        temp_man.Omega_strk = stroke_var.Omega_strk(:,j:man_end);
                        temp_man.Omega_dot_strk = stroke_var.Omega_dot_strk(:,j:man_end);
                        temp_man.V_strk = stroke_var.V_strk(:,j:man_end);
                        temp_man.A_strk = stroke_var.A_strk(:,j:man_end);
                            
                        temp_man.f = man_a_fit.f(j:man_end);
                        temp_man.down_up = man_down_up_ratio(j:man_end);
                        
                        pure_roll.(char(['man_' int2str(i) '_' int2str(pure_roll_nr) ])) = temp_man;
                        
                        end
                        
                    end
                    
                end
                
                
                if maneuver.a_x(j) ~= 0
                    
                    if maneuver.a_x(j-1) == 0
                        
                        man_end = j+find(maneuver.a_x(j:end)==0,1,'first')-2;
                        
                        ay_on = sum(maneuver.a_y(j:man_end));
                        
                        az_on = sum(maneuver.a_z(j:man_end));
                        
                        if ay_on == 0 && az_on == 0
                        
                        pure_ax_nr = pure_ax_nr+1;
                        
                        temp_man = {};
                        
                        movie_percentage = j/nr_wb;
                        
                        temp_man.seq_nr = i;
                        temp_man.seq_name = settings.sequence_names(i);
                        temp_man.movie_percentage = movie_percentage;
                        temp_man.wing_l = L_wing;
                        temp_man.mass_fly = mass_fly;
                        temp_man.a_avg = pathDB.a_avg.(char(['a_avg_' int2str(i) ]));
                        temp_man.f_avg = pathDB.f_avg.(char(['f_avg_' int2str(i) ]));
                        temp_man.down_up_avg = pathDB.down_up_avg.(char(['down_up_avg_' int2str(i) ]));
                                                
                        man_a_dev = pathDB.a_dev.(char(['a_dev_' int2str(i) ]));
                        man_a_fit = pathDB.a_fit.(char(['a_fit_' int2str(i) ]));
                        man_down_up_ratio = pathDB.down_up_ratio.(char(['down_up_ratio_' int2str(i) ]));
                                                
                        temp_man.a_dev.theta_L1 = man_a_dev.theta_L1(:,j:man_end);
                        temp_man.a_dev.theta_L2 = man_a_dev.theta_L2(:,j:man_end);
                        temp_man.a_dev.eta_L1 = man_a_dev.eta_L1(:,j:man_end);
                        temp_man.a_dev.eta_L2 = man_a_dev.eta_L2(:,j:man_end);
                        temp_man.a_dev.phi_L1 = man_a_dev.phi_L1(:,j:man_end);
                        temp_man.a_dev.phi_L2 = man_a_dev.phi_L2(:,j:man_end);
                        
                        temp_man.a_dev.theta_R1 = man_a_dev.theta_R1(:,j:man_end);
                        temp_man.a_dev.theta_R2 = man_a_dev.theta_R2(:,j:man_end);
                        temp_man.a_dev.eta_R1 = man_a_dev.eta_R1(:,j:man_end);
                        temp_man.a_dev.eta_R2 = man_a_dev.eta_R2(:,j:man_end);
                        temp_man.a_dev.phi_R1 = man_a_dev.phi_R1(:,j:man_end);
                        temp_man.a_dev.phi_R2 = man_a_dev.phi_R2(:,j:man_end);
                        
                        temp_man.q_body_wb = stroke_var.q_body(:,j:man_end);
                        temp_man.Omega_strk = stroke_var.Omega_strk(:,j:man_end);
                        temp_man.Omega_dot_strk = stroke_var.Omega_dot_strk(:,j:man_end);
                        temp_man.V_strk = stroke_var.V_strk(:,j:man_end);
                        temp_man.A_strk = stroke_var.A_strk(:,j:man_end);
                            
                        temp_man.f = man_a_fit.f(j:man_end);
                        temp_man.down_up = man_down_up_ratio(j:man_end);
                        
                        pure_ax.(char(['man_' int2str(i) '_' int2str(pure_ax_nr) ])) = temp_man;
                        
                        end
                        
                    end
                    
                end
                
                
                if maneuver.a_y(j) ~= 0
                    
                    if maneuver.a_y(j-1) == 0
                        
                        man_end = j+find(maneuver.a_y(j:end)==0,1,'first')-2;
                        
                        ax_on = sum(maneuver.a_x(j:man_end));
                        
                        az_on = sum(maneuver.a_z(j:man_end));
                        
                        if ax_on == 0 && az_on == 0
                        
                        pure_ay_nr = pure_ay_nr+1;
                        
                        temp_man = {};
                        
                        movie_percentage = j/nr_wb;
                        
                        temp_man.seq_nr = i;
                        temp_man.seq_name = settings.sequence_names(i);
                        temp_man.movie_percentage = movie_percentage;
                        temp_man.wing_l = L_wing;
                        temp_man.mass_fly = mass_fly;
                        temp_man.a_avg = pathDB.a_avg.(char(['a_avg_' int2str(i) ]));
                        temp_man.f_avg = pathDB.f_avg.(char(['f_avg_' int2str(i) ]));
                        temp_man.down_up_avg = pathDB.down_up_avg.(char(['down_up_avg_' int2str(i) ]));
                                                
                        man_a_dev = pathDB.a_dev.(char(['a_dev_' int2str(i) ]));
                        man_a_fit = pathDB.a_fit.(char(['a_fit_' int2str(i) ]));
                        man_down_up_ratio = pathDB.down_up_ratio.(char(['down_up_ratio_' int2str(i) ]));
                                                
                        temp_man.a_dev.theta_L1 = man_a_dev.theta_L1(:,j:man_end);
                        temp_man.a_dev.theta_L2 = man_a_dev.theta_L2(:,j:man_end);
                        temp_man.a_dev.eta_L1 = man_a_dev.eta_L1(:,j:man_end);
                        temp_man.a_dev.eta_L2 = man_a_dev.eta_L2(:,j:man_end);
                        temp_man.a_dev.phi_L1 = man_a_dev.phi_L1(:,j:man_end);
                        temp_man.a_dev.phi_L2 = man_a_dev.phi_L2(:,j:man_end);
                        
                        temp_man.a_dev.theta_R1 = man_a_dev.theta_R1(:,j:man_end);
                        temp_man.a_dev.theta_R2 = man_a_dev.theta_R2(:,j:man_end);
                        temp_man.a_dev.eta_R1 = man_a_dev.eta_R1(:,j:man_end);
                        temp_man.a_dev.eta_R2 = man_a_dev.eta_R2(:,j:man_end);
                        temp_man.a_dev.phi_R1 = man_a_dev.phi_R1(:,j:man_end);
                        temp_man.a_dev.phi_R2 = man_a_dev.phi_R2(:,j:man_end);
                        
                        temp_man.q_body_wb = stroke_var.q_body(:,j:man_end);
                        temp_man.Omega_strk = stroke_var.Omega_strk(:,j:man_end);
                        temp_man.Omega_dot_strk = stroke_var.Omega_dot_strk(:,j:man_end);
                        temp_man.V_strk = stroke_var.V_strk(:,j:man_end);
                        temp_man.A_strk = stroke_var.A_strk(:,j:man_end);
                            
                        temp_man.f = man_a_fit.f(j:man_end);
                        temp_man.down_up = man_down_up_ratio(j:man_end);
                        
                        pure_ay.(char(['man_' int2str(i) '_' int2str(pure_ay_nr) ])) = temp_man;
                        
                        end
                        
                    end
                    
                end
                
                if maneuver.a_z(j) ~= 0
                    
                    if maneuver.a_z(j-1) == 0
                        
                        man_end = j+find(maneuver.a_z(j:end)==0,1,'first')-2;
                        
                        ax_on = sum(maneuver.a_x(j:man_end));
                        
                        ay_on = sum(maneuver.a_y(j:man_end));
                        
                        if ax_on == 0 && ay_on == 0
                        
                        pure_az_nr = pure_az_nr+1;
                        
                        temp_man = {};
                        
                        movie_percentage = j/nr_wb;
                        
                        temp_man.seq_nr = i;
                        temp_man.seq_name = settings.sequence_names(i);
                        temp_man.movie_percentage = movie_percentage;
                        temp_man.wing_l = L_wing;
                        temp_man.mass_fly = mass_fly;
                        temp_man.a_avg = pathDB.a_avg.(char(['a_avg_' int2str(i) ]));
                        temp_man.f_avg = pathDB.f_avg.(char(['f_avg_' int2str(i) ]));
                        temp_man.down_up_avg = pathDB.down_up_avg.(char(['down_up_avg_' int2str(i) ]));
                                                
                        man_a_dev = pathDB.a_dev.(char(['a_dev_' int2str(i) ]));
                        man_a_fit = pathDB.a_fit.(char(['a_fit_' int2str(i) ]));
                        man_down_up_ratio = pathDB.down_up_ratio.(char(['down_up_ratio_' int2str(i) ]));
                                                
                        temp_man.a_dev.theta_L1 = man_a_dev.theta_L1(:,j:man_end);
                        temp_man.a_dev.theta_L2 = man_a_dev.theta_L2(:,j:man_end);
                        temp_man.a_dev.eta_L1 = man_a_dev.eta_L1(:,j:man_end);
                        temp_man.a_dev.eta_L2 = man_a_dev.eta_L2(:,j:man_end);
                        temp_man.a_dev.phi_L1 = man_a_dev.phi_L1(:,j:man_end);
                        temp_man.a_dev.phi_L2 = man_a_dev.phi_L2(:,j:man_end);
                        
                        temp_man.a_dev.theta_R1 = man_a_dev.theta_R1(:,j:man_end);
                        temp_man.a_dev.theta_R2 = man_a_dev.theta_R2(:,j:man_end);
                        temp_man.a_dev.eta_R1 = man_a_dev.eta_R1(:,j:man_end);
                        temp_man.a_dev.eta_R2 = man_a_dev.eta_R2(:,j:man_end);
                        temp_man.a_dev.phi_R1 = man_a_dev.phi_R1(:,j:man_end);
                        temp_man.a_dev.phi_R2 = man_a_dev.phi_R2(:,j:man_end);
                        
                        temp_man.q_body_wb = stroke_var.q_body(:,j:man_end);
                        temp_man.Omega_strk = stroke_var.Omega_strk(:,j:man_end);
                        temp_man.Omega_dot_strk = stroke_var.Omega_dot_strk(:,j:man_end);
                        temp_man.V_strk = stroke_var.V_strk(:,j:man_end);
                        temp_man.A_strk = stroke_var.A_strk(:,j:man_end);
                            
                        temp_man.f = man_a_fit.f(j:man_end);
                        temp_man.down_up = man_down_up_ratio(j:man_end);
                        
                        pure_az.(char(['man_' int2str(i) '_' int2str(pure_az_nr) ])) = temp_man;
                        
                        end
                        
                    end
                    
                end
                
                
            end


        end       

    end
    
    all_maneuvers.roll_maneuvers = roll_maneuvers; 
    all_maneuvers.pitch_maneuvers = pitch_maneuvers; 
    all_maneuvers.yaw_maneuvers = yaw_maneuvers; 
    all_maneuvers.ax_maneuvers = ax_maneuvers; 
    all_maneuvers.ay_maneuvers = ay_maneuvers; 
    all_maneuvers.az_maneuvers = az_maneuvers; 
    all_maneuvers.pure_roll = pure_roll;
    all_maneuvers.pure_pitch = pure_pitch;
    all_maneuvers.pure_yaw = pure_yaw;
    all_maneuvers.pure_ax = pure_ax;
    all_maneuvers.pure_ay = pure_ay;
    all_maneuvers.pure_az = pure_az;
    
    
end

