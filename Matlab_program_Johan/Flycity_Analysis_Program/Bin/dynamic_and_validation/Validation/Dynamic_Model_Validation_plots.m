function Dynamic_Model_Validation_plots( settings, pathDB, DynSim_avg, DynSim_man, DynSim_steady )

    % Validate the results of the dynamic model by comparison with the
    % actual movie data:
    
    R_strk = pathDB.rot_mat.Rstr;
    
    dt = pathDB.dt;
    
    steady_wb_names     = fieldnames(DynSim_steady);
    
    nr_wb_steady        = length(steady_wb_names);
    
    v_strk_sim_steady       = zeros(3,nr_wb_steady);
    w_strk_sim_steady       = zeros(3,nr_wb_steady);
    a_strk_sim_steady       = zeros(3,nr_wb_steady);
    w_dot_strk_sim_steady   = zeros(3,nr_wb_steady);
    
    v_strk_mov_steady       = zeros(3,nr_wb_steady);
    w_strk_mov_steady       = zeros(3,nr_wb_steady);
    a_strk_mov_steady       = zeros(3,nr_wb_steady);
    w_dot_strk_mov_steady   = zeros(3,nr_wb_steady);
    
    for k = 1:nr_wb_steady
        
        v_strk_sim_steady(:,k)      = R_strk*DynSim_steady.(char(steady_wb_names(k))).sim_data.vb_mean;
        w_strk_sim_steady(:,k)      = R_strk*DynSim_steady.(char(steady_wb_names(k))).sim_data.wb_mean;
        a_strk_sim_steady(:,k)      = R_strk*DynSim_steady.(char(steady_wb_names(k))).sim_data.ab_mean;
        w_dot_strk_sim_steady(:,k)  = R_strk*DynSim_steady.(char(steady_wb_names(k))).sim_data.w_dot_b_mean;
        
        qb_mean = pathDB.rand_wbs.(char(steady_wb_names(k))).body_kin.qb_mean;
        Rb_mean = quat2mat(qb_mean);
        
        v_strk_mov_steady(:,k)      = R_strk*Rb_mean*pathDB.rand_wbs.(char(steady_wb_names(k))).body_kin.uvw_mean';
        w_strk_mov_steady(:,k)      = R_strk*pathDB.rand_wbs.(char(steady_wb_names(k))).body_kin.wb_mean';
        a_strk_mov_steady(:,k)      = R_strk*Rb_mean*pathDB.rand_wbs.(char(steady_wb_names(k))).body_kin.a_xyz_mean';
        
        % Create w_dot_b_mean:
        
        w_dot_b_mean = mean(gradient(pathDB.rand_wbs.(char(steady_wb_names(k))).body_kin.wb')./dt,2);
        
        w_dot_strk_mov_steady(:,k)  = R_strk*w_dot_b_mean;
        
    end
    
    man_wb_names        = fieldnames(DynSim_man);
    
    nr_wb_man           = length(man_wb_names);
    
    v_strk_sim_man      = zeros(3,nr_wb_man);
    w_strk_sim_man      = zeros(3,nr_wb_man);
    a_strk_sim_man      = zeros(3,nr_wb_man);
    w_dot_strk_sim_man  = zeros(3,nr_wb_man);
    
    v_strk_mov_man      = zeros(3,nr_wb_man);
    w_strk_mov_man      = zeros(3,nr_wb_man);
    a_strk_mov_man      = zeros(3,nr_wb_man);
    w_dot_strk_mov_man  = zeros(3,nr_wb_man);
    
    seq_nr_list     = zeros(nr_wb_man,1);
    wb_nr_list      = zeros(nr_wb_man,1);
    man_type_list   = zeros(nr_wb_man,6);
    
    for k = 1:nr_wb_man
                
        if strncmpi(man_wb_names(k),'ax',2) == 1
                        
            t_string    = char(man_wb_names(k));
            man_type    = 'ax';
            man_type_list(k,1) = k;
            wb_nr       = str2num(t_string(7:end));
            seq_nr      = pathDB.maneuver.ax.(char(['wb_' num2str(wb_nr)])).seq_nr;
            
        end
                    
        if strncmpi(man_wb_names(k),'ay',2) == 1
                        
            t_string    = char(man_wb_names(k));
            man_type    = 'ay';
            man_type_list(k,2) = k;
            wb_nr       = str2num(t_string(7:end));
            seq_nr      = pathDB.maneuver.ay.(char(['wb_' num2str(wb_nr)])).seq_nr;
                        
        end
        
        if strncmpi(man_wb_names(k),'az',2) == 1
                        
            t_string    = char(man_wb_names(k));
            man_type    = 'az';
            man_type_list(k,3) = k;
            wb_nr       = str2num(t_string(7:end));
            seq_nr      = pathDB.maneuver.az.(char(['wb_' num2str(wb_nr)])).seq_nr;
            
        end
                    
        if strncmpi(man_wb_names(k),'wx',2) == 1
                        
            t_string    = char(man_wb_names(k));
            man_type    = 'wx';
            man_type_list(k,4) = k;
            wb_nr       = str2num(t_string(7:end));
            seq_nr      = pathDB.maneuver.wx.(char(['wb_' num2str(wb_nr)])).seq_nr;
                        
        end
        
        if strncmpi(man_wb_names(k),'wy',2) == 1
                        
            t_string    = char(man_wb_names(k));
            man_type    = 'wy';
            man_type_list(k,5) = k;
            wb_nr       = str2num(t_string(7:end));
            seq_nr      = pathDB.maneuver.wy.(char(['wb_' num2str(wb_nr)])).seq_nr;
                        
        end
        
        if strncmpi(man_wb_names(k),'wz',2) == 1
                        
            t_string    = char(man_wb_names(k));
            man_type    = 'wz';
            man_type_list(k,6) = k;
            wb_nr       = str2num(t_string(7:end));
            seq_nr      = pathDB.maneuver.wz.(char(['wb_' num2str(wb_nr)])).seq_nr;
                        
        end
        
        seq_nr_list(k)  = seq_nr;
        wb_nr_list(k)   = wb_nr;
                
        v_strk_sim_man(:,k) 	= R_strk*DynSim_man.(char([ char(man_type) '_wb_' num2str(wb_nr)])).sim_data.vb_mean;
        w_strk_sim_man(:,k)     = R_strk*DynSim_man.(char([ char(man_type) '_wb_' num2str(wb_nr)])).sim_data.wb_mean;
        a_strk_sim_man(:,k)     = R_strk*DynSim_man.(char([ char(man_type) '_wb_' num2str(wb_nr)])).sim_data.ab_mean;
        w_dot_strk_sim_man(:,k) = R_strk*DynSim_man.(char([ char(man_type) '_wb_' num2str(wb_nr)])).sim_data.w_dot_b_mean;
        
        qb_mean = pathDB.maneuver.(char(man_type)).(char(['wb_' num2str(wb_nr)])).body_kin.qb_mean;
        Rb_mean = quat2mat(qb_mean);
        
        v_strk_mov_man(:,k)     = R_strk*Rb_mean*pathDB.maneuver.(char(man_type)).(char(['wb_' num2str(wb_nr)])).body_kin.uvw_mean';
        w_strk_mov_man(:,k)     = R_strk*pathDB.maneuver.(char(man_type)).(char(['wb_' num2str(wb_nr)])).body_kin.wb_mean';
        a_strk_mov_man(:,k)     = R_strk*Rb_mean*pathDB.maneuver.(char(man_type)).(char(['wb_' num2str(wb_nr)])).body_kin.a_xyz_mean';
        
        % Create w_dot_b_mean:
        
        w_dot_b_mean = mean(gradient(pathDB.maneuver.(char(man_type)).(char(['wb_' num2str(wb_nr)])).body_kin.wb')./dt,2);
        
        w_dot_strk_mov_man(:,k)  = R_strk*w_dot_b_mean;
        
    end
    
    ax_man = find(man_type_list(:,1));
    ay_man = find(man_type_list(:,2));
    az_man = find(man_type_list(:,3));
    wx_man = find(man_type_list(:,4));
    wy_man = find(man_type_list(:,5));
    wz_man = find(man_type_list(:,6));
    
    % Plot velocity:
    
    error_plot( v_strk_sim_steady, v_strk_mov_steady, v_strk_sim_man, v_strk_mov_man, ax_man, 1);
    error_plot( v_strk_sim_steady, v_strk_mov_steady, v_strk_sim_man, v_strk_mov_man, ay_man, 2);
    error_plot( v_strk_sim_steady, v_strk_mov_steady, v_strk_sim_man, v_strk_mov_man, az_man, 3);
    
    % Plot angular velocity:
    
    error_plot( w_strk_sim_steady, w_strk_mov_steady, w_strk_sim_man, w_strk_mov_man, wx_man, 1);
    error_plot( w_strk_sim_steady, w_strk_mov_steady, w_strk_sim_man, w_strk_mov_man, wy_man, 2);
    error_plot( w_strk_sim_steady, w_strk_mov_steady, w_strk_sim_man, w_strk_mov_man, wz_man, 3);
    
    % Plot acceleration:
    
    error_plot( a_strk_sim_steady, a_strk_mov_steady, a_strk_sim_man, a_strk_mov_man, ax_man, 1);
    error_plot( a_strk_sim_steady, a_strk_mov_steady, a_strk_sim_man, a_strk_mov_man, ay_man, 2);
    error_plot( a_strk_sim_steady, a_strk_mov_steady, a_strk_sim_man, a_strk_mov_man, az_man, 3);
    
    % Plot angular acceleration:
    
    error_plot( w_dot_strk_sim_steady, w_dot_strk_mov_steady, w_dot_strk_sim_man, w_dot_strk_mov_man, wx_man, 1);
    error_plot( w_dot_strk_sim_steady, w_dot_strk_mov_steady, w_dot_strk_sim_man, w_dot_strk_mov_man, wy_man, 2);
    error_plot( w_dot_strk_sim_steady, w_dot_strk_mov_steady, w_dot_strk_sim_man, w_dot_strk_mov_man, wz_man, 3);
    

end