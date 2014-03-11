function pathDB9( settings, pathDB )


    savefile = 'pathDB9.mat';


        
    all_maneuvers1 = non_dimensional_control_analysis( settings, pathDB);
    
    [ F_M_r ] = F_M_req( settings, pathDB, all_maneuvers1 );
    
    
    
    roll_maneuvers = all_maneuvers1.roll_maneuvers;
    pitch_maneuvers = all_maneuvers1.pitch_maneuvers;
    yaw_maneuvers = all_maneuvers1.yaw_maneuvers;
    ax_maneuvers = all_maneuvers1.ax_maneuvers;
    ay_maneuvers = all_maneuvers1.ay_maneuvers;
    az_maneuvers = all_maneuvers1.az_maneuvers;
    pure_roll = all_maneuvers1.pure_roll;
    pure_pitch = all_maneuvers1.pure_pitch;
    pure_yaw = all_maneuvers1.pure_yaw;
    pure_ax = all_maneuvers1.pure_ax;
    pure_ay = all_maneuvers1.pure_ay;
    pure_az = all_maneuvers1.pure_az;
    
    roll_man_name = fieldnames(roll_maneuvers);
    pitch_man_name = fieldnames(pitch_maneuvers);
    yaw_man_name = fieldnames(yaw_maneuvers);
    ax_man_name = fieldnames(ax_maneuvers);
    ay_man_name = fieldnames(ay_maneuvers);
    az_man_name = fieldnames(az_maneuvers);
    pure_roll_name = fieldnames(pure_roll);
    pure_pitch_name = fieldnames(pure_pitch);
    pure_yaw_name = fieldnames(pure_yaw);
    pure_ax_name = fieldnames(pure_ax);
    pure_ay_name = fieldnames(pure_ay);
    pure_az_name = fieldnames(pure_az);
    
    nr_roll_man = size(roll_man_name,1);
    nr_pitch_man = size(pitch_man_name,1);
    nr_yaw_man = size(yaw_man_name,1);
    nr_ax_man = size(ax_man_name,1);
    nr_ay_man = size(ay_man_name,1);
    nr_az_man = size(az_man_name,1);
    nr_pure_roll = size(pure_roll_name,1);
    nr_pure_pitch = size(pure_pitch_name,1);
    nr_pure_yaw = size(pure_yaw_name,1);
    nr_pure_ax = size(pure_ax_name,1);
    nr_pure_ay = size(pure_ay_name,1);
    nr_pure_az = size(pure_az_name,1);
        
    
    all_maneuvers = {};
    
    for i = 1:nr_roll_man
        
        temp_man = roll_maneuvers.(char(roll_man_name(i)));
        
        temp_man.FM_req = F_M_r.roll_maneuver.(char(roll_man_name(i)));
        
        all_maneuvers.roll_maneuvers.(char(roll_man_name(i))) = temp_man;
        
    end
    
    
    
    for i = 1:nr_pitch_man
        
        temp_man = pitch_maneuvers.(char(pitch_man_name(i)));
        
        temp_man.FM_req = F_M_r.pitch_maneuver.(char(pitch_man_name(i)));
        
        all_maneuvers.pitch_maneuvers.(char(pitch_man_name(i))) = temp_man;
        
    end
    
    
    
    for i = 1:nr_yaw_man
        
        temp_man = yaw_maneuvers.(char(yaw_man_name(i)));
        
        temp_man.FM_req = F_M_r.yaw_maneuver.(char(yaw_man_name(i)));
        
        all_maneuvers.yaw_maneuvers.(char(yaw_man_name(i))) = temp_man;
        
    end
    
    
    
    for i = 1:nr_ax_man
        
        temp_man = ax_maneuvers.(char(ax_man_name(i)));
        
        temp_man.FM_req = F_M_r.ax_maneuver.(char(ax_man_name(i)));
        
        all_maneuvers.ax_maneuvers.(char(ax_man_name(i))) = temp_man;
        
    end
    
    
    
    for i = 1:nr_ay_man
        
        temp_man = ay_maneuvers.(char(ay_man_name(i)));
        
        temp_man.FM_req = F_M_r.ay_maneuver.(char(ay_man_name(i)));
        
        all_maneuvers.ay_maneuvers.(char(ay_man_name(i))) = temp_man;
        
    end
    
    
    
    for i = 1:nr_az_man
        
        temp_man = az_maneuvers.(char(az_man_name(i)));
        
        temp_man.FM_req = F_M_r.az_maneuver.(char(az_man_name(i)));
        
        all_maneuvers.az_maneuvers.(char(az_man_name(i))) = temp_man;
        
    end
    
    
    
    for i = 1:nr_pure_roll
        
        temp_man = pure_roll.(char(pure_roll_name(i)));
        
        temp_man.FM_req = F_M_r.pure_roll_maneuver.(char(pure_roll_name(i)));
        
        all_maneuvers.pure_roll.(char(pure_roll_name(i))) = temp_man;
        
    end
    
    
    for i = 1:nr_pure_pitch
        
        temp_man = pure_pitch.(char(pure_pitch_name(i)));
        
        temp_man.FM_req = F_M_r.pure_pitch_maneuver.(char(pure_pitch_name(i)));
        
        all_maneuvers.pure_pitch.(char(pure_pitch_name(i))) = temp_man;
        
    end
    
    
    
    for i = 1:nr_pure_yaw
        
        temp_man = pure_yaw.(char(pure_yaw_name(i)));
        
        temp_man.FM_req = F_M_r.pure_yaw_maneuver.(char(pure_yaw_name(i)));
        
        all_maneuvers.pure_yaw.(char(pure_yaw_name(i))) = temp_man;
        
    end
    
    
    
    for i = 1:nr_pure_ax
        
        temp_man = pure_ax.(char(pure_ax_name(i)));
        
        temp_man.FM_req = F_M_r.pure_ax_maneuver.(char(pure_ax_name(i)));
        
        all_maneuvers.pure_ax.(char(pure_ax_name(i))) = temp_man;
        
    end
    
    
    
    for i = 1:nr_pure_ay
        
        temp_man = pure_ay.(char(pure_ay_name(i)));
        
        temp_man.FM_req = F_M_r.pure_ay_maneuver.(char(pure_ay_name(i)));
        
        all_maneuvers.pure_ay.(char(pure_ay_name(i))) = temp_man;
        
    end
    
    
    
    for i = 1:nr_pure_az
        
        temp_man = pure_az.(char(pure_az_name(i)));
        
        temp_man.FM_req = F_M_r.pure_az_maneuver.(char(pure_az_name(i)));
        
        all_maneuvers.pure_az.(char(pure_az_name(i))) = temp_man;
        
    end
       
    save(savefile,'all_maneuvers')

end


