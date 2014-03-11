function Force_Moment_evaluation( settings, pathDB )

    % Evaluate the forces and moments of the dynamic simulation with the
    % forces and moments obtained from the movies:
    
    R_strk = pathDB.rot_mat.Rstr;
    
    DynSim_steady      = pathDB.DynSim_steady;
    steady_wb_names    = fieldnames(DynSim_steady);
    nr_wb_steady       = length(steady_wb_names);
    
       
    sim.FI_acc_steady  = zeros(3,nr_wb_steady);
    sim.MI_acc_steady  = zeros(3,nr_wb_steady);
    sim.FI_vel_steady  = zeros(3,nr_wb_steady);
    sim.MI_vel_steady  = zeros(3,nr_wb_steady);
    sim.FA_steady      = zeros(3,nr_wb_steady);
    sim.MA_steady      = zeros(3,nr_wb_steady);
    sim.Fg_steady      = zeros(3,nr_wb_steady);
    sim.F_0_steady     = zeros(3,nr_wb_steady);
    sim.M_0_steady     = zeros(3,nr_wb_steady);
    
    raw.qb_steady           = zeros(4,nr_wb_steady);
    raw.Rb_steady           = zeros(3,3,nr_wb_steady);
    raw.vb_steady           = zeros(3,nr_wb_steady);
    raw.wb_steady           = zeros(3,nr_wb_steady);
    raw.ab_steady           = zeros(3,nr_wb_steady);
    raw.w_dot_b_steady      = zeros(3,nr_wb_steady);
    raw.v_strk_steady       = zeros(3,nr_wb_steady);
    raw.w_strk_steady       = zeros(3,nr_wb_steady);
    raw.a_strk_steady       = zeros(3,nr_wb_steady);
    raw.w_dot_strk_steady   = zeros(3,nr_wb_steady);
    raw.FI_acc_steady       = zeros(3,nr_wb_steady);
    raw.MI_acc_steady       = zeros(3,nr_wb_steady);
    raw.FI_vel_steady       = zeros(3,nr_wb_steady);
    raw.MI_vel_steady       = zeros(3,nr_wb_steady);
    raw.FA_steady           = zeros(3,nr_wb_steady);
    raw.MA_steady           = zeros(3,nr_wb_steady);
    raw.Fg_steady           = zeros(3,nr_wb_steady);
    
    nr_wb_steady
    
    
    for k = 1:nr_wb_steady
        
        k
        
        seq_nr_steady   = pathDB.rand_wbs.(char(steady_wb_names(k))).seq_nr;
        wb_nr_steady    = pathDB.rand_wbs.(char(steady_wb_names(k))).wb_nr;
        
        FM_raw = FM_raw_movies(settings, pathDB, wb_nr_steady, seq_nr_steady);
        
        raw.qb_steady(:,k)           = FM_raw.qb;
        raw.Rb_steady(:,:,k)         = FM_raw.Rb;
        raw.vb_steady(:,k)           = FM_raw.vb;
        raw.wb_steady(:,k)           = FM_raw.wb;
        raw.ab_steady(:,k)           = FM_raw.ab;
        raw.w_dot_b_steady(:,k)      = FM_raw.w_dot_b;
        raw.v_strk_steady(:,k)       = R_strk*FM_raw.vb;
        raw.w_strk_steady(:,k)       = R_strk*FM_raw.wb;
        raw.a_strk_steady(:,k)       = R_strk*FM_raw.ab;
        raw.w_dot_strk_steady(:,k)   = R_strk*FM_raw.w_dot_b;
        raw.FI_acc_steady(:,k)       = FM_raw.FI_acc;
        raw.MI_acc_steady(:,k)       = FM_raw.MI_acc;
        raw.FI_vel_steady(:,k)       = FM_raw.FI_vel;
        raw.MI_vel_steady(:,k)       = FM_raw.MI_vel;
        raw.FA_steady(:,k)           = FM_raw.FA;
        raw.MA_steady(:,k)           = FM_raw.MA;
        raw.Fg_steady(:,k)           = FM_raw.Fg;
        
        sim_data = DynSim_steady.(char(steady_wb_names(k))).sim_data;
        
        sim.FI_acc_steady(:,k)       = sim_data.FI_acc_strk_mean;
        sim.MI_acc_steady(:,k)       = sim_data.MI_acc_strk_mean;
        sim.FI_vel_steady(:,k)       = sim_data.FI_vel_strk_mean;
        sim.MI_vel_steady(:,k)       = sim_data.MI_vel_strk_mean;
        sim.FA_steady(:,k)           = sim_data.FA_strk_mean;
        sim.MA_steady(:,k)           = sim_data.MA_strk_mean;
        sim.Fg_steady(:,k)           = sim_data.Fg_strk_mean;
        sim.F_0_steady(:,k)          = sim_data.F_0;
        sim.M_0_steady(:,k)          = sim_data.M_0;
        
    end
    
    DynSim_man          = pathDB.DynSim_man;
    man_wb_names        = fieldnames(DynSim_man);
    nr_wb_man           = length(man_wb_names);
    
    seq_nr_list     = zeros(nr_wb_man,1);
    wb_nr_list      = zeros(nr_wb_man,1);
    man_type_list   = zeros(nr_wb_man,6);
    
    sim.FI_acc_man  = zeros(3,nr_wb_man);
    sim.MI_acc_man  = zeros(3,nr_wb_man);
    sim.FI_vel_man  = zeros(3,nr_wb_man);
    sim.MI_vel_man  = zeros(3,nr_wb_man);
    sim.FA_man      = zeros(3,nr_wb_man);
    sim.MA_man      = zeros(3,nr_wb_man);
    sim.Fg_man      = zeros(3,nr_wb_man);
    sim.F_0_man     = zeros(3,nr_wb_man);
    sim.M_0_man     = zeros(3,nr_wb_man);
    
    raw.qb_man           = zeros(4,nr_wb_man);
    raw.Rb_man           = zeros(3,3,nr_wb_man);
    raw.vb_man           = zeros(3,nr_wb_man);
    raw.wb_man           = zeros(3,nr_wb_man);
    raw.ab_man           = zeros(3,nr_wb_man);
    raw.w_dot_b_man      = zeros(3,nr_wb_man);
    raw.v_strk_man       = zeros(3,nr_wb_man);
    raw.w_strk_man       = zeros(3,nr_wb_man);
    raw.a_strk_man       = zeros(3,nr_wb_man);
    raw.w_dot_strk_man   = zeros(3,nr_wb_man);
    raw.FI_acc_man       = zeros(3,nr_wb_man);
    raw.MI_acc_man       = zeros(3,nr_wb_man);
    raw.FI_vel_man       = zeros(3,nr_wb_man);
    raw.MI_vel_man       = zeros(3,nr_wb_man);
    raw.FA_man           = zeros(3,nr_wb_man);
    raw.MA_man           = zeros(3,nr_wb_man);
    raw.Fg_man           = zeros(3,nr_wb_man);   
    
    nr_wb_man
    
    for k = 1:nr_wb_man
        
        k
                
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
                
        seq_nr_man      = seq_nr;
        wb_nr_man       = pathDB.maneuver.(man_type).(['wb_' num2str(wb_nr)]).wb_nr;
        
        FM_raw = FM_raw_movies(settings, pathDB, wb_nr_man, seq_nr_man);
        
        raw.qb_man(:,k)                  = FM_raw.qb;
        raw.Rb_man(:,:,k)                = FM_raw.Rb;
        raw.vb_man(:,k)                  = FM_raw.vb;
        raw.wb_man(:,k)                  = FM_raw.wb;
        raw.ab_man(:,k)                  = FM_raw.ab;
        raw.w_dot_b_man(:,k)             = FM_raw.w_dot_b;
        raw.v_strk_man(:,k)              = R_strk*FM_raw.vb;
        raw.w_strk_man(:,k)              = R_strk*FM_raw.wb;
        raw.a_strk_man(:,k)              = R_strk*FM_raw.ab;
        raw.w_dot_strk_man(:,k)          = R_strk*FM_raw.w_dot_b;
        raw.FI_acc_man(:,k)              = FM_raw.FI_acc;
        raw.MI_acc_man(:,k)              = FM_raw.MI_acc;
        raw.FI_vel_man(:,k)              = FM_raw.FI_vel;
        raw.MI_vel_man(:,k)              = FM_raw.MI_vel;
        raw.FA_man(:,k)                  = FM_raw.FA;
        raw.MA_man(:,k)                  = FM_raw.MA;
        raw.Fg_man(:,k)                  = FM_raw.Fg;
        
        sim_data = DynSim_man.(char(man_wb_names(k))).sim_data;
        
        sim.FI_acc_man(:,k)              = sim_data.FI_acc_strk_mean;
        sim.MI_acc_man(:,k)              = sim_data.MI_acc_strk_mean;
        sim.FI_vel_man(:,k)              = sim_data.FI_vel_strk_mean;
        sim.MI_vel_man(:,k)              = sim_data.MI_vel_strk_mean;
        sim.FA_man(:,k)                  = sim_data.FA_strk_mean;
        sim.MA_man(:,k)                  = sim_data.MA_strk_mean;
        sim.Fg_man(:,k)                  = sim_data.Fg_strk_mean;
        sim.F_0_man(:,k)                 = sim_data.F_0;
        sim.M_0_man(:,k)                 = sim_data.M_0;
        
    end
    
    ax_man = find(man_type_list(:,1));
    ay_man = find(man_type_list(:,2));
    az_man = find(man_type_list(:,3));
    wx_man = find(man_type_list(:,4));
    wy_man = find(man_type_list(:,5));
    wz_man = find(man_type_list(:,6));
    
    FI_acc_avg_steady  	= mean(abs(raw.FI_acc_steady-sim.FI_acc_steady),2)
    FI_acc_avg_man      = mean(abs(raw.FI_acc_man-sim.FI_acc_man),2)
    MI_acc_avg_steady   = mean(abs(raw.MI_acc_steady-sim.MI_acc_steady),2)
    MI_acc_avg_man      = mean(abs(raw.MI_acc_man-sim.MI_acc_man),2)
    FI_vel_avg_steady   = mean(abs(raw.FI_vel_steady-sim.FI_vel_steady),2)
    FI_vel_avg_man      = mean(abs(raw.FI_vel_man-sim.FI_vel_man),2)
    MI_vel_avg_steady   = mean(abs(raw.MI_vel_steady-sim.MI_vel_steady),2)
    MI_vel_avg_man      = mean(abs(raw.MI_vel_man-sim.MI_vel_man),2)
    FA_avg_steady       = mean(abs(raw.FA_steady-sim.FA_steady),2)
    FA_avg_man          = mean(abs(raw.FA_man-sim.FA_man),2)
    MA_avg_steady       = mean(abs(raw.MA_steady-sim.MA_steady),2)
    MA_avg_man          = mean(abs(raw.MA_man-sim.MA_man),2)
    Fg_avg_steady       = mean(abs(raw.Fg_steady-sim.Fg_steady),2)
    Fg_avg_man          = mean(abs(raw.Fg_man-sim.Fg_man),2)
    
    Fg_ref = 1.85e-6*9.81;
    Mg_ref = 1.85e-6*9.81;
   
%     figure()
%     hold on
%     subplot(3,1,1); hold on
%     plot(raw.v_strk_steady(1,:),abs(raw.FI_acc_steady(1,:)-sim.FI_acc_steady(1,:))/Fg_ref,'o','Color','r')
%     plot(raw.v_strk_man(1,:),abs(raw.FI_acc_man(1,:)-sim.FI_acc_man(1,:))/Fg_ref,'o','Color','b')
%     hold off
%     subplot(3,1,2); hold on
%     plot(raw.v_strk_steady(2,:),abs(raw.FI_acc_steady(2,:)-sim.FI_acc_steady(2,:))/Fg_ref,'o','Color','r')
%     plot(raw.v_strk_man(2,:),abs(raw.FI_acc_man(2,:)-sim.FI_acc_man(2,:))/Fg_ref,'o','Color','b')
%     hold off
%     subplot(3,1,3); hold on
%     plot(raw.v_strk_steady(3,:),abs(raw.FI_acc_steady(3,:)-sim.FI_acc_steady(3,:))/Fg_ref,'o','Color','r')
%     plot(raw.v_strk_man(3,:),abs(raw.FI_acc_man(3,:)-sim.FI_acc_man(3,:))/Fg_ref,'o','Color','b')
%     hold off
%     title('v_b vs FI_{acc}')
%     hold off
    
    figure()
    hold on
    subplot(3,1,1); hold on
    plot([raw.FI_acc_steady(1,:) raw.FI_acc_man(1,:)],[raw.FI_acc_steady(1,:) raw.FI_acc_man(1,:)],'o','Color','r')
    plot([raw.FI_acc_steady(1,:) raw.FI_acc_man(1,:)],[sim.FI_acc_steady(1,:) sim.FI_acc_man(1,:)],'o','Color','b')
    hold off
    axis equal
    subplot(3,1,2); hold on
    plot([raw.FI_acc_steady(2,:) raw.FI_acc_man(2,:)],[raw.FI_acc_steady(2,:) raw.FI_acc_man(2,:)],'o','Color','r')
    plot([raw.FI_acc_steady(2,:) raw.FI_acc_man(2,:)],[sim.FI_acc_steady(2,:) sim.FI_acc_man(2,:)],'o','Color','b')
    hold off
    axis equal
    subplot(3,1,3); hold on
    plot([raw.FI_acc_steady(3,:) raw.FI_acc_man(3,:)],[raw.FI_acc_steady(3,:) raw.FI_acc_man(3,:)],'o','Color','r')
    plot([raw.FI_acc_steady(3,:) raw.FI_acc_man(3,:)],[sim.FI_acc_steady(3,:) sim.FI_acc_man(3,:)],'o','Color','b')
    hold off
    axis equal
    hold off
    
    figure()
    hold on
    subplot(3,1,1); hold on
    plot([raw.MI_acc_steady(1,:) raw.MI_acc_man(1,:)],[raw.MI_acc_steady(1,:) raw.MI_acc_man(1,:)],'o','Color','r')
    plot([raw.MI_acc_steady(1,:) raw.MI_acc_man(1,:)],[sim.MI_acc_steady(1,:) sim.MI_acc_man(1,:)],'o','Color','b')
    hold off
    axis equal
    subplot(3,1,2); hold on
    plot([raw.MI_acc_steady(2,:) raw.MI_acc_man(2,:)],[raw.MI_acc_steady(2,:) raw.MI_acc_man(2,:)],'o','Color','r')
    plot([raw.MI_acc_steady(2,:) raw.MI_acc_man(2,:)],[sim.MI_acc_steady(2,:) sim.MI_acc_man(2,:)],'o','Color','b')
    hold off
    axis equal
    subplot(3,1,3); hold on
    plot([raw.MI_acc_steady(3,:) raw.MI_acc_man(3,:)],[raw.MI_acc_steady(3,:) raw.MI_acc_man(3,:)],'o','Color','r')
    plot([raw.MI_acc_steady(3,:) raw.MI_acc_man(3,:)],[sim.MI_acc_steady(3,:) sim.MI_acc_man(3,:)],'o','Color','b')
    hold off
    axis equal
    hold off
    
    figure()
    hold on
    subplot(3,1,1); hold on
    plot([raw.Fg_steady(1,:) raw.Fg_man(1,:)],[raw.Fg_steady(1,:) raw.Fg_man(1,:)],'o','Color','r')
    plot([raw.Fg_steady(1,:) raw.Fg_man(1,:)],[sim.Fg_steady(1,:) sim.Fg_man(1,:)],'o','Color','b')
    hold off
    axis equal
    subplot(3,1,2); hold on
    plot([raw.Fg_steady(2,:) raw.Fg_man(2,:)],[raw.Fg_steady(2,:) raw.Fg_man(2,:)],'o','Color','r')
    plot([raw.Fg_steady(2,:) raw.Fg_man(2,:)],[sim.Fg_steady(2,:) sim.Fg_man(2,:)],'o','Color','b')
    hold off
    axis equal
    subplot(3,1,3); hold on
    plot([raw.Fg_steady(3,:) raw.Fg_man(3,:)],[raw.Fg_steady(3,:) raw.Fg_man(3,:)],'o','Color','r')
    plot([raw.Fg_steady(3,:) raw.Fg_man(3,:)],[sim.Fg_steady(3,:) sim.Fg_man(3,:)],'o','Color','b')
    hold off
    axis equal
    hold off   

    pause
    
%     figure()
%     hold on
%     subplot(3,1,1); hold on
%     plot(raw.v_strk_steady(1,:),abs(raw.FI_acc_steady(1,:)-(sim.FI_acc_steady(1,:)+sim.F_0_steady(1,:))),'o','Color','r')
%     plot(raw.v_strk_man(1,:),abs(raw.FI_acc_man(1,:)-(sim.FI_acc_man(1,:)+sim.F_0_man(1,:))),'o','Color','b')
%     hold off
%     subplot(3,1,2); hold on
%     plot(raw.v_strk_steady(2,:),abs(raw.FI_acc_steady(2,:)-(sim.FI_acc_steady(2,:)+sim.F_0_steady(2,:))),'o','Color','r')
%     plot(raw.v_strk_man(2,:),abs(raw.FI_acc_man(2,:)-(sim.FI_acc_man(2,:)+sim.F_0_man(2,:))),'o','Color','b')
%     hold off
%     subplot(3,1,3); hold on
%     plot(raw.v_strk_steady(3,:),abs(raw.FI_acc_steady(3,:)-(sim.FI_acc_steady(3,:)+sim.F_0_steady(3,:))),'o','Color','r')
%     plot(raw.v_strk_man(3,:),abs(raw.FI_acc_man(3,:)-(sim.FI_acc_man(3,:)+sim.F_0_man(3,:))),'o','Color','b')
%     hold off
%     title('v_b vs FI_{acc}')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); hold on
%     plot(raw.v_strk_steady(1,:),abs(raw.MI_acc_steady(1,:)-(sim.MI_acc_steady(1,:)+sim.M_0_steady(1,:))),'o','Color','r')
%     plot(raw.v_strk_man(1,:),abs(raw.MI_acc_man(1,:)-(sim.MI_acc_man(1,:)+sim.M_0_man(1,:))),'o','Color','b')
%     hold off
%     subplot(3,1,2); hold on
%     plot(raw.v_strk_steady(2,:),abs(raw.MI_acc_steady(2,:)-(sim.MI_acc_steady(2,:)+sim.M_0_steady(2,:))),'o','Color','r')
%     plot(raw.v_strk_man(2,:),abs(raw.MI_acc_man(2,:)-(sim.MI_acc_man(2,:)+sim.M_0_man(2,:))),'o','Color','b')
%     hold off
%     subplot(3,1,3); hold on
%     plot(raw.v_strk_steady(3,:),abs(raw.MI_acc_steady(3,:)-(sim.MI_acc_steady(3,:)+sim.M_0_steady(3,:))),'o','Color','r')
%     plot(raw.v_strk_man(3,:),abs(raw.MI_acc_man(3,:)-(sim.MI_acc_man(3,:)+sim.M_0_man(3,:))),'o','Color','b')
%     hold off
%     title('v_b vs FI_{acc}')
%     hold off
    
    figure()
    hold on
    subplot(3,1,1); hold on
    plot(raw.v_strk_steady(1,:),abs(raw.MI_acc_steady(1,:)-sim.MI_acc_steady(1,:))/Mg_ref,'o','Color','r')
    plot(raw.v_strk_man(1,:),abs(raw.MI_acc_man(1,:)-sim.MI_acc_man(1,:))/Mg_ref,'o','Color','b')
    hold off
    subplot(3,1,2); hold on
    plot(raw.v_strk_steady(2,:),abs(raw.MI_acc_steady(2,:)-sim.MI_acc_steady(2,:))/Mg_ref,'o','Color','r')
    plot(raw.v_strk_man(2,:),abs(raw.MI_acc_man(2,:)-sim.MI_acc_man(2,:))/Mg_ref,'o','Color','b')
    hold off
    subplot(3,1,3); hold on
    plot(raw.v_strk_steady(3,:),abs(raw.MI_acc_steady(3,:)-sim.MI_acc_steady(3,:))/Mg_ref,'o','Color','r')
    plot(raw.v_strk_man(3,:),abs(raw.MI_acc_man(3,:)-sim.MI_acc_man(3,:))/Mg_ref,'o','Color','b')
    hold off
    title('v_b vs MI_{acc}')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); hold on
    plot(raw.v_strk_steady(1,:),abs(raw.FI_vel_steady(1,:)-sim.FI_vel_steady(1,:))/Fg_ref,'o','Color','r')
    plot(raw.v_strk_man(1,:),abs(raw.FI_vel_man(1,:)-sim.FI_vel_man(1,:))/Fg_ref,'o','Color','b')
    hold off
    subplot(3,1,2); hold on
    plot(raw.v_strk_steady(2,:),abs(raw.FI_vel_steady(2,:)-sim.FI_vel_steady(2,:))/Fg_ref,'o','Color','r')
    plot(raw.v_strk_man(2,:),abs(raw.FI_vel_man(2,:)-sim.FI_vel_man(2,:))/Fg_ref,'o','Color','b')
    hold off
    subplot(3,1,3); hold on
    plot(raw.v_strk_steady(3,:),abs(raw.FI_vel_steady(3,:)-sim.FI_vel_steady(3,:))/Fg_ref,'o','Color','r')
    plot(raw.v_strk_man(3,:),abs(raw.FI_vel_man(3,:)-sim.FI_vel_man(3,:))/Fg_ref,'o','Color','b')
    hold off
    title('v_b vs FI_{vel}')
    hold off
        
    figure()
    hold on
    subplot(3,1,1); hold on
    plot(raw.v_strk_steady(1,:),abs(raw.MI_vel_steady(1,:)-sim.MI_vel_steady(1,:))/Mg_ref,'o','Color','r')
    plot(raw.v_strk_man(1,:),abs(raw.MI_vel_man(1,:)-sim.MI_vel_man(1,:))/Mg_ref,'o','Color','b')
    hold off
    subplot(3,1,2); hold on
    plot(raw.v_strk_steady(2,:),abs(raw.MI_vel_steady(2,:)-sim.MI_vel_steady(2,:))/Mg_ref,'o','Color','r')
    plot(raw.v_strk_man(2,:),abs(raw.MI_vel_man(2,:)-sim.MI_vel_man(2,:))/Mg_ref,'o','Color','b')
    hold off
    subplot(3,1,3); hold on
    plot(raw.v_strk_steady(3,:),abs(raw.MI_vel_steady(3,:)-sim.MI_vel_steady(3,:))/Mg_ref,'o','Color','r')
    plot(raw.v_strk_man(3,:),abs(raw.MI_vel_man(3,:)-sim.MI_vel_man(3,:))/Mg_ref,'o','Color','b')
    hold off
    title('v_b vs MI_{vel}')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); hold on
    plot(raw.v_strk_steady(1,:),abs(raw.FA_steady(1,:)-sim.FA_steady(1,:))/Fg_ref,'o','Color','r')
    plot(raw.v_strk_man(1,:),abs(raw.FA_man(1,:)-sim.FA_man(1,:))/Fg_ref,'o','Color','b')
    hold off
    subplot(3,1,2); hold on
    plot(raw.v_strk_steady(2,:),abs(raw.FA_steady(2,:)-sim.FA_steady(2,:))/Fg_ref,'o','Color','r')
    plot(raw.v_strk_man(2,:),abs(raw.FA_man(2,:)-sim.FA_man(2,:))/Fg_ref,'o','Color','b')
    hold off
    subplot(3,1,3); hold on
    plot(raw.v_strk_steady(3,:),abs(raw.FA_steady(3,:)-sim.FA_steady(3,:))/Fg_ref,'o','Color','r')
    plot(raw.v_strk_man(3,:),abs(raw.FA_man(3,:)-sim.FA_man(3,:))/Fg_ref,'o','Color','b')
    hold off
    title('v_b vs FA')
    hold off
        
    figure()
    hold on
    subplot(3,1,1); hold on
    plot(raw.v_strk_steady(1,:),abs(raw.MA_steady(1,:)-sim.MA_steady(1,:))/Mg_ref,'o','Color','r')
    plot(raw.v_strk_man(1,:),abs(raw.MA_man(1,:)-sim.MA_man(1,:))/Mg_ref,'o','Color','b')
    hold off
    subplot(3,1,2); hold on
    plot(raw.v_strk_steady(2,:),abs(raw.MA_steady(2,:)-sim.MA_steady(2,:))/Mg_ref,'o','Color','r')
    plot(raw.v_strk_man(2,:),abs(raw.MA_man(2,:)-sim.MA_man(2,:))/Mg_ref,'o','Color','b')
    hold off
    subplot(3,1,3); hold on
    plot(raw.v_strk_steady(3,:),abs(raw.MA_steady(3,:)-sim.MA_steady(3,:))/Mg_ref,'o','Color','r')
    plot(raw.v_strk_man(3,:),abs(raw.MA_man(3,:)-sim.MA_man(3,:))/Mg_ref,'o','Color','b')
    hold off
    title('v_b vs MA')
    hold off

    figure()
    hold on
    subplot(3,1,1); hold on
    plot(raw.v_strk_steady(1,:),abs(raw.Fg_steady(1,:))/Fg_ref,'o','Color','r')
    plot(raw.v_strk_man(1,:),abs(raw.Fg_man(1,:))/Fg_ref,'o','Color','b')
    plot(raw.v_strk_steady(1,:),abs(sim.Fg_steady(1,:))/Fg_ref,'o','Color','r')
    plot(raw.v_strk_man(1,:),abs(sim.Fg_man(1,:))/Fg_ref,'o','Color','b')
    hold off
    subplot(3,1,2); hold on
    plot(raw.v_strk_steady(2,:),abs(raw.Fg_steady(2,:))/Fg_ref,'o','Color','r')
    plot(raw.v_strk_man(2,:),abs(raw.Fg_man(2,:))/Fg_ref,'o','Color','b')
    plot(raw.v_strk_steady(2,:),abs(sim.Fg_steady(2,:))/Fg_ref,'o','Color','r')
    plot(raw.v_strk_man(2,:),abs(sim.Fg_man(2,:))/Fg_ref,'o','Color','b')
    hold off
    subplot(3,1,3); hold on
    plot(raw.v_strk_steady(3,:),abs(raw.Fg_steady(3,:))/Fg_ref,'o','Color','r')
    plot(raw.v_strk_man(3,:),abs(raw.Fg_man(3,:))/Fg_ref,'o','Color','b')
    plot(raw.v_strk_steady(3,:),abs(sim.Fg_steady(3,:))/Fg_ref,'o','Color','g')
    plot(raw.v_strk_man(3,:),abs(sim.Fg_man(3,:))/Fg_ref,'o','Color','b')
    hold off
    title('v_b vs F_g')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); hold on
    plot(raw.w_strk_steady(1,:),abs(raw.FI_acc_steady(1,:)-sim.FI_acc_steady(1,:))/Fg_ref,'o','Color','r')
    plot(raw.w_strk_man(1,:),abs(raw.FI_acc_man(1,:)-sim.FI_acc_man(1,:))/Fg_ref,'o','Color','b')
    hold off
    subplot(3,1,2); hold on
    plot(raw.w_strk_steady(2,:),abs(raw.FI_acc_steady(2,:)-sim.FI_acc_steady(2,:))/Fg_ref,'o','Color','r')
    plot(raw.w_strk_man(2,:),abs(raw.FI_acc_man(2,:)-sim.FI_acc_man(2,:))/Fg_ref,'o','Color','b')
    hold off
    subplot(3,1,3); hold on
    plot(raw.w_strk_steady(3,:),abs(raw.FI_acc_steady(3,:)-sim.FI_acc_steady(3,:))/Fg_ref,'o','Color','r')
    plot(raw.w_strk_man(3,:),abs(raw.FI_acc_man(3,:)-sim.FI_acc_man(3,:))/Fg_ref,'o','Color','b')
    hold off
    title('w_{strk} vs FI_{acc}')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); hold on
    plot(raw.w_strk_steady(1,:),abs(raw.MI_acc_steady(1,:)-sim.MI_acc_steady(1,:))/Mg_ref,'o','Color','r')
    plot(raw.w_strk_man(1,:),abs(raw.MI_acc_man(1,:)-sim.MI_acc_man(1,:))/Mg_ref,'o','Color','b')
    hold off
    subplot(3,1,2); hold on
    plot(raw.w_strk_steady(2,:),abs(raw.MI_acc_steady(2,:)-sim.MI_acc_steady(2,:))/Mg_ref,'o','Color','r')
    plot(raw.w_strk_man(2,:),abs(raw.MI_acc_man(2,:)-sim.MI_acc_man(2,:))/Mg_ref,'o','Color','b')
    hold off
    subplot(3,1,3); hold on
    plot(raw.w_strk_steady(3,:),abs(raw.MI_acc_steady(3,:)-sim.MI_acc_steady(3,:))/Mg_ref,'o','Color','r')
    plot(raw.w_strk_man(3,:),abs(raw.MI_acc_man(3,:)-sim.MI_acc_man(3,:))/Mg_ref,'o','Color','b')
    hold off
    title('w_{strk} vs MI_{acc}')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); hold on
    plot(raw.w_strk_steady(1,:),abs(raw.FI_vel_steady(1,:)-sim.FI_vel_steady(1,:))/Fg_ref,'o','Color','r')
    plot(raw.w_strk_man(1,:),abs(raw.FI_vel_man(1,:)-sim.FI_vel_man(1,:))/Fg_ref,'o','Color','b')
    hold off
    subplot(3,1,2); hold on
    plot(raw.w_strk_steady(2,:),abs(raw.FI_vel_steady(2,:)-sim.FI_vel_steady(2,:))/Fg_ref,'o','Color','r')
    plot(raw.w_strk_man(2,:),abs(raw.FI_vel_man(2,:)-sim.FI_vel_man(2,:))/Fg_ref,'o','Color','b')
    hold off
    subplot(3,1,3); hold on
    plot(raw.w_strk_steady(3,:),abs(raw.FI_vel_steady(3,:)-sim.FI_vel_steady(3,:))/Fg_ref,'o','Color','r')
    plot(raw.w_strk_man(3,:),abs(raw.FI_vel_man(3,:)-sim.FI_vel_man(3,:))/Fg_ref,'o','Color','b')
    hold off
    title('w_{strk} vs FI_{vel}')
    hold off
        
    figure()
    hold on
    subplot(3,1,1); hold on
    plot(raw.w_strk_steady(1,:),abs(raw.MI_vel_steady(1,:)-sim.MI_vel_steady(1,:))/Mg_ref,'o','Color','r')
    plot(raw.w_strk_man(1,:),abs(raw.MI_vel_man(1,:)-sim.MI_vel_man(1,:))/Mg_ref,'o','Color','b')
    hold off
    subplot(3,1,2); hold on
    plot(raw.w_strk_steady(2,:),abs(raw.MI_vel_steady(2,:)-sim.MI_vel_steady(2,:))/Mg_ref,'o','Color','r')
    plot(raw.w_strk_man(2,:),abs(raw.MI_vel_man(2,:)-sim.MI_vel_man(2,:))/Mg_ref,'o','Color','b')
    hold off
    subplot(3,1,3); hold on
    plot(raw.w_strk_steady(3,:),abs(raw.MI_vel_steady(3,:)-sim.MI_vel_steady(3,:))/Mg_ref,'o','Color','r')
    plot(raw.w_strk_man(3,:),abs(raw.MI_vel_man(3,:)-sim.MI_vel_man(3,:))/Mg_ref,'o','Color','b')
    hold off
    title('w_{strk} vs MI_{vel}')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); hold on
    plot(raw.w_strk_steady(1,:),abs(raw.FA_steady(1,:)-sim.FA_steady(1,:))/Fg_ref,'o','Color','r')
    plot(raw.w_strk_man(1,:),abs(raw.FA_man(1,:)-sim.FA_man(1,:))/Fg_ref,'o','Color','b')
    hold off
    subplot(3,1,2); hold on
    plot(raw.w_strk_steady(2,:),abs(raw.FA_steady(2,:)-sim.FA_steady(2,:))/Fg_ref,'o','Color','r')
    plot(raw.w_strk_man(2,:),abs(raw.FA_man(2,:)-sim.FA_man(2,:))/Fg_ref,'o','Color','b')
    hold off
    subplot(3,1,3); hold on
    plot(raw.w_strk_steady(3,:),abs(raw.FA_steady(3,:)-sim.FA_steady(3,:))/Fg_ref,'o','Color','r')
    plot(raw.w_strk_man(3,:),abs(raw.FA_man(3,:)-sim.FA_man(3,:))/Fg_ref,'o','Color','b')
    hold off
    title('w_{strk} vs FA')
    hold off
        
    figure()
    hold on
    subplot(3,1,1); hold on
    plot(raw.w_strk_steady(1,:),abs(raw.MA_steady(1,:)-sim.MA_steady(1,:))/Mg_ref,'o','Color','r')
    plot(raw.w_strk_man(1,:),abs(raw.MA_man(1,:)-sim.MA_man(1,:))/Mg_ref,'o','Color','b')
    hold off
    subplot(3,1,2); hold on
    plot(raw.w_strk_steady(2,:),abs(raw.MA_steady(2,:)-sim.MA_steady(2,:))/Mg_ref,'o','Color','r')
    plot(raw.w_strk_man(2,:),abs(raw.MA_man(2,:)-sim.MA_man(2,:))/Mg_ref,'o','Color','b')
    hold off
    subplot(3,1,3); hold on
    plot(raw.w_strk_steady(3,:),abs(raw.MA_steady(3,:)-sim.MA_steady(3,:))/Mg_ref,'o','Color','r')
    plot(raw.w_strk_man(3,:),abs(raw.MA_man(3,:)-sim.MA_man(3,:))/Mg_ref,'o','Color','b')
    hold off
    title('w_{strk} vs MA')
    hold off

    figure()
    hold on
    subplot(3,1,1); hold on
    plot(raw.w_strk_steady(1,:),abs(raw.Fg_steady(1,:))/Fg_ref,'o','Color','r')
    plot(raw.w_strk_man(1,:),abs(raw.Fg_man(1,:))/Fg_ref,'o','Color','b')
    plot(raw.w_strk_steady(1,:),abs(sim.Fg_steady(1,:))/Fg_ref,'o','Color','r')
    plot(raw.w_strk_man(1,:),abs(sim.Fg_man(1,:))/Fg_ref,'o','Color','b')
    hold off
    subplot(3,1,2); hold on
    plot(raw.w_strk_steady(2,:),abs(raw.Fg_steady(2,:))/Fg_ref,'o','Color','r')
    plot(raw.w_strk_man(2,:),abs(raw.Fg_man(2,:))/Fg_ref,'o','Color','b')
    plot(raw.w_strk_steady(2,:),abs(sim.Fg_steady(2,:))/Fg_ref,'o','Color','r')
    plot(raw.w_strk_man(2,:),abs(sim.Fg_man(2,:))/Fg_ref,'o','Color','b')
    hold off
    subplot(3,1,3); hold on
    plot(raw.w_strk_steady(3,:),abs(raw.Fg_steady(3,:))/Fg_ref,'o','Color','r')
    plot(raw.w_strk_man(3,:),abs(raw.Fg_man(3,:))/Fg_ref,'o','Color','b')
    plot(raw.w_strk_steady(3,:),abs(sim.Fg_steady(3,:))/Fg_ref,'o','Color','g')
    plot(raw.w_strk_man(3,:),abs(sim.Fg_man(3,:))/Fg_ref,'o','Color','b')
    hold off
    title('w_{strk} vs F_g')
    hold off
    
    
    figure()
    hold on
    subplot(3,1,1); hold on
    plot(raw.a_strk_steady(1,:),abs(raw.FI_acc_steady(1,:)-sim.FI_acc_steady(1,:))/Fg_ref,'o','Color','r')
    plot(raw.a_strk_man(1,:),abs(raw.FI_acc_man(1,:)-sim.FI_acc_man(1,:))/Fg_ref,'o','Color','b')
    hold off
    subplot(3,1,2); hold on
    plot(raw.a_strk_steady(2,:),abs(raw.FI_acc_steady(2,:)-sim.FI_acc_steady(2,:))/Fg_ref,'o','Color','r')
    plot(raw.a_strk_man(2,:),abs(raw.FI_acc_man(2,:)-sim.FI_acc_man(2,:))/Fg_ref,'o','Color','b')
    hold off
    subplot(3,1,3); hold on
    plot(raw.a_strk_steady(3,:),abs(raw.FI_acc_steady(3,:)-sim.FI_acc_steady(3,:))/Fg_ref,'o','Color','r')
    plot(raw.a_strk_man(3,:),abs(raw.FI_acc_man(3,:)-sim.FI_acc_man(3,:))/Fg_ref,'o','Color','b')
    hold off
    title('a_{strk} vs FI_{acc}')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); hold on
    plot(raw.a_strk_steady(1,:),abs(raw.MI_acc_steady(1,:)-sim.MI_acc_steady(1,:))/Mg_ref,'o','Color','r')
    plot(raw.a_strk_man(1,:),abs(raw.MI_acc_man(1,:)-sim.MI_acc_man(1,:))/Mg_ref,'o','Color','b')
    hold off
    subplot(3,1,2); hold on
    plot(raw.a_strk_steady(2,:),abs(raw.MI_acc_steady(2,:)-sim.MI_acc_steady(2,:))/Mg_ref,'o','Color','r')
    plot(raw.a_strk_man(2,:),abs(raw.MI_acc_man(2,:)-sim.MI_acc_man(2,:))/Mg_ref,'o','Color','b')
    hold off
    subplot(3,1,3); hold on
    plot(raw.a_strk_steady(3,:),abs(raw.MI_acc_steady(3,:)-sim.MI_acc_steady(3,:))/Mg_ref,'o','Color','r')
    plot(raw.a_strk_man(3,:),abs(raw.MI_acc_man(3,:)-sim.MI_acc_man(3,:))/Mg_ref,'o','Color','b')
    hold off
    title('a_{strk} vs MI_{acc}')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); hold on
    plot(raw.a_strk_steady(1,:),abs(raw.FI_vel_steady(1,:)-sim.FI_vel_steady(1,:))/Fg_ref,'o','Color','r')
    plot(raw.a_strk_man(1,:),abs(raw.FI_vel_man(1,:)-sim.FI_vel_man(1,:))/Fg_ref,'o','Color','b')
    hold off
    subplot(3,1,2); hold on
    plot(raw.a_strk_steady(2,:),abs(raw.FI_vel_steady(2,:)-sim.FI_vel_steady(2,:))/Fg_ref,'o','Color','r')
    plot(raw.a_strk_man(2,:),abs(raw.FI_vel_man(2,:)-sim.FI_vel_man(2,:))/Fg_ref,'o','Color','b')
    hold off
    subplot(3,1,3); hold on
    plot(raw.a_strk_steady(3,:),abs(raw.FI_vel_steady(3,:)-sim.FI_vel_steady(3,:))/Fg_ref,'o','Color','r')
    plot(raw.a_strk_man(3,:),abs(raw.FI_vel_man(3,:)-sim.FI_vel_man(3,:))/Fg_ref,'o','Color','b')
    hold off
    title('a_{strk} vs FI_{vel}')
    hold off
        
    figure()
    hold on
    subplot(3,1,1); hold on
    plot(raw.a_strk_steady(1,:),abs(raw.MI_vel_steady(1,:)-sim.MI_vel_steady(1,:))/Mg_ref,'o','Color','r')
    plot(raw.a_strk_man(1,:),abs(raw.MI_vel_man(1,:)-sim.MI_vel_man(1,:))/Mg_ref,'o','Color','b')
    hold off
    subplot(3,1,2); hold on
    plot(raw.a_strk_steady(2,:),abs(raw.MI_vel_steady(2,:)-sim.MI_vel_steady(2,:))/Mg_ref,'o','Color','r')
    plot(raw.a_strk_man(2,:),abs(raw.MI_vel_man(2,:)-sim.MI_vel_man(2,:))/Mg_ref,'o','Color','b')
    hold off
    subplot(3,1,3); hold on
    plot(raw.a_strk_steady(3,:),abs(raw.MI_vel_steady(3,:)-sim.MI_vel_steady(3,:))/Mg_ref,'o','Color','r')
    plot(raw.a_strk_man(3,:),abs(raw.MI_vel_man(3,:)-sim.MI_vel_man(3,:))/Mg_ref,'o','Color','b')
    hold off
    title('a_{strk} vs MI_{vel}')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); hold on
    plot(raw.a_strk_steady(1,:),abs(raw.FA_steady(1,:)-sim.FA_steady(1,:))/Fg_ref,'o','Color','r')
    plot(raw.a_strk_man(1,:),abs(raw.FA_man(1,:)-sim.FA_man(1,:))/Fg_ref,'o','Color','b')
    hold off
    subplot(3,1,2); hold on
    plot(raw.a_strk_steady(2,:),abs(raw.FA_steady(2,:)-sim.FA_steady(2,:))/Fg_ref,'o','Color','r')
    plot(raw.a_strk_man(2,:),abs(raw.FA_man(2,:)-sim.FA_man(2,:))/Fg_ref,'o','Color','b')
    hold off
    subplot(3,1,3); hold on
    plot(raw.a_strk_steady(3,:),abs(raw.FA_steady(3,:)-sim.FA_steady(3,:))/Fg_ref,'o','Color','r')
    plot(raw.a_strk_man(3,:),abs(raw.FA_man(3,:)-sim.FA_man(3,:))/Fg_ref,'o','Color','b')
    hold off
    title('a_{strk} vs FA')
    hold off
        
    figure()
    hold on
    subplot(3,1,1); hold on
    plot(raw.a_strk_steady(1,:),abs(raw.MA_steady(1,:)-sim.MA_steady(1,:))/Mg_ref,'o','Color','r')
    plot(raw.a_strk_man(1,:),abs(raw.MA_man(1,:)-sim.MA_man(1,:))/Mg_ref,'o','Color','b')
    hold off
    subplot(3,1,2); hold on
    plot(raw.a_strk_steady(2,:),abs(raw.MA_steady(2,:)-sim.MA_steady(2,:))/Mg_ref,'o','Color','r')
    plot(raw.a_strk_man(2,:),abs(raw.MA_man(2,:)-sim.MA_man(2,:))/Mg_ref,'o','Color','b')
    hold off
    subplot(3,1,3); hold on
    plot(raw.a_strk_steady(3,:),abs(raw.MA_steady(3,:)-sim.MA_steady(3,:))/Mg_ref,'o','Color','r')
    plot(raw.a_strk_man(3,:),abs(raw.MA_man(3,:)-sim.MA_man(3,:))/Mg_ref,'o','Color','b')
    hold off
    title('a_{strk} vs MA')
    hold off

    figure()
    hold on
    subplot(3,1,1); hold on
    plot(raw.a_strk_steady(1,:),abs(raw.Fg_steady(1,:))/Fg_ref,'o','Color','r')
    plot(raw.a_strk_man(1,:),abs(raw.Fg_man(1,:))/Fg_ref,'o','Color','b')
    plot(raw.a_strk_steady(1,:),abs(sim.Fg_steady(1,:))/Fg_ref,'o','Color','r')
    plot(raw.a_strk_man(1,:),abs(sim.Fg_man(1,:))/Fg_ref,'o','Color','b')
    hold off
    subplot(3,1,2); hold on
    plot(raw.a_strk_steady(2,:),abs(raw.Fg_steady(2,:))/Fg_ref,'o','Color','r')
    plot(raw.a_strk_man(2,:),abs(raw.Fg_man(2,:))/Fg_ref,'o','Color','b')
    plot(raw.a_strk_steady(2,:),abs(sim.Fg_steady(2,:))/Fg_ref,'o','Color','r')
    plot(raw.a_strk_man(2,:),abs(sim.Fg_man(2,:))/Fg_ref,'o','Color','b')
    hold off
    subplot(3,1,3); hold on
    plot(raw.a_strk_steady(3,:),abs(raw.Fg_steady(3,:))/Fg_ref,'o','Color','r')
    plot(raw.a_strk_man(3,:),abs(raw.Fg_man(3,:))/Fg_ref,'o','Color','b')
    plot(raw.a_strk_steady(3,:),abs(sim.Fg_steady(3,:))/Fg_ref,'o','Color','g')
    plot(raw.a_strk_man(3,:),abs(sim.Fg_man(3,:))/Fg_ref,'o','Color','b')
    hold off
    title('a_{strk} vs F_g')
    hold off
    
end

function [ FM_raw ] = FM_raw_movies(settings, pathDB, wb_nr, seq_nr)

    wb_range = pathDB.wingbeats.wingbeat_loc(wb_nr,1,seq_nr):pathDB.wingbeats.wingbeat_loc(wb_nr,2,seq_nr);
    dt       = pathDB.dt;
    
    R_strk   = pathDB.rot_mat.Rstr;
    
    N        = length(wb_range);
    
    qb       = pathDB.filt.qB(wb_range,:,seq_nr)';
    wb       = pathDB.filt.wB(wb_range,:,seq_nr)';
    
    Rb       = zeros(3,3,N);
    vb       = zeros(3,N);
    ab       = zeros(3,N);
    w_dot_b  = gradient(wb)./dt;
    
    u_strk   = zeros(3,N);
    w_strk   = zeros(3,N);
    
    for i = 1:N
        
        Rb(:,:,i)    = quat2mat(qb(:,i));
        vb(:,i)      = Rb(:,:,i)*pathDB.filt.uvw(wb_range(i),:,seq_nr)';
        ab(:,i)      = Rb(:,:,i)*pathDB.filt.a_xyz(wb_range(i),:,seq_nr)';
        u_strk(:,i)  = R_strk*Rb(:,:,i)*pathDB.filt.uvw(wb_range(i),:,seq_nr)';
        w_strk(:,i)  = R_strk*wb(:,i);
        
    end
    
    % Construct the body and wing model:
    
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
    
    % Create the wing kinematics:
    
%     a_theta_L           = pathDB.poly_fit.a_fit.theta_L(:,wb_nr,seq_nr);
%     a_eta_L             = pathDB.poly_fit.a_fit.eta_L(:,wb_nr,seq_nr);
%     a_phi_L             = pathDB.poly_fit.a_fit.phi_L(:,wb_nr,seq_nr);
%     a_theta_R           = pathDB.poly_fit.a_fit.theta_R(:,wb_nr,seq_nr);
%     a_eta_R             = pathDB.poly_fit.a_fit.eta_R(:,wb_nr,seq_nr);
%     a_phi_R             = pathDB.poly_fit.a_fit.phi_R(:,wb_nr,seq_nr);
%     f                   = pathDB.poly_fit.a_fit.f(wb_nr,seq_nr);
%     down_up             = pathDB.poly_fit.a_fit.down_up(wb_nr,seq_nr);

    a_theta_L           = pathDB.poly_fit.a_avg.theta_LR(:,seq_nr) + pathDB.poly_fit.a_dev.theta_L(:,wb_nr,seq_nr);
    a_eta_L             = pathDB.poly_fit.a_avg.eta_LR(:,seq_nr) + pathDB.poly_fit.a_dev.eta_L(:,wb_nr,seq_nr);
    a_phi_L             = pathDB.poly_fit.a_avg.phi_LR(:,seq_nr) + pathDB.poly_fit.a_dev.phi_L(:,wb_nr,seq_nr);
    a_theta_R           = pathDB.poly_fit.a_avg.theta_LR(:,seq_nr) + pathDB.poly_fit.a_dev.theta_R(:,wb_nr,seq_nr);
    a_eta_R             = pathDB.poly_fit.a_avg.eta_LR(:,seq_nr) + pathDB.poly_fit.a_dev.eta_R(:,wb_nr,seq_nr);
    a_phi_R             = pathDB.poly_fit.a_avg.phi_LR(:,seq_nr) + pathDB.poly_fit.a_dev.phi_R(:,wb_nr,seq_nr);
    f                   = pathDB.poly_fit.a_fit.f(wb_nr,seq_nr);
    down_up             = pathDB.poly_fit.a_fit.down_up(wb_nr,seq_nr);
    
    a_fit.a_theta_L     = a_theta_L;
    a_fit.a_eta_L       = a_eta_L;
    a_fit.a_phi_L       = a_phi_L;
    a_fit.a_theta_R     = a_theta_R;
    a_fit.a_eta_R       = a_eta_R;
    a_fit.a_phi_R       = a_phi_R;
    a_fit.f             = f;
    a_fit.down_up       = down_up;
    a_fit.nr_points     = N;
    a_fit.R_strk        = R_strk;
        
    [ kine ] = angular_velocities_polynomial( a_fit );
    
    t_func              = kine.t;
    kine.R_strk         = R_strk;
    kine.down_up_t      = t_func(end)*down_up;

    clear a_fit
    
    % Compute FI_acc & MI_acc:
    
    FI_acc = zeros(3,N);
    MI_acc = zeros(3,N);
    
    for i = 1:N
        
        body_kin.vb     = vb(:,i);
        body_kin.wb     = wb(:,i);
        
        wing_kin.RL     = kine.RL(:,:,i);
        wing_kin.RR     = kine.RR(:,:,i);

        wing_kin.wL     = kine.wL(:,i);
        wing_kin.wR     = kine.wR(:,i);
        wing_kin.wL_b   = kine.RL(:,:,i)'*kine.wL(:,i);
        wing_kin.wR_b   = kine.RR(:,:,i)'*kine.wR(:,i);
        
        M_matrix        = Mass_matrix( body_model, wing_model, body_kin, wing_kin );
        
        M = [M_matrix.M11 M_matrix.M12; M_matrix.M21 M_matrix.M22];

%         M = [M_matrix.M11 zeros(3); zeros(3) M_matrix.M22];
        
        FM_acc_t = 1e-3*M*[ab(:,i); w_dot_b(:,i)];
        
        FI_acc(:,i) = R_strk*FM_acc_t(1:3);
        MI_acc(:,i) = R_strk*FM_acc_t(4:6);
        
    end
    
    % Compute aerodynamic forces and moments:
    
    kine.u_strk = u_strk;
    kine.w_strk = w_strk;
    
    [ FM_strkpln, ~, ~ ,~, ~, ~, ~, ~, ~ ] = Aerodynamic_forces( kine, body_model, wing_model, 1);
    
    FA = FM_strkpln(1:3,:);
    MA = FM_strkpln(4:6,:);
    
    % Compute gravity vector:
    
    Fg = zeros(3,N);
    
    for i = 1:N
        
        body_kin.Rb     = Rb(:,:,i);
        body_kin.R_strk = R_strk;
    
        [ ~, Fg_t ] = Gravity_instantaneous( body_model, body_kin );
        
        Fg(:,i) = Fg_t;
    
    end
    
    % Compute the averages:
    
    FI_acc_avg      = mean(FI_acc,2);
    MI_acc_avg      = mean(MI_acc,2);
    
    FA_avg          = mean(FA,2);
    MA_avg          = mean(MA,2);
    
    Fg_avg          = mean(Fg,2);
    
    FI_vel_avg      = FI_acc_avg-FA_avg-Fg_avg;
    MI_vel_avg      = MI_acc_avg-MA_avg;
    
    FM_raw.qb       = mean(qb,2);
    FM_raw.Rb       = quat2mat(FM_raw.qb);
    FM_raw.vb       = mean(vb,2);
    FM_raw.wb       = mean(wb,2);
    FM_raw.ab       = mean(ab,2);
    FM_raw.w_dot_b  = mean(w_dot_b,2);
    FM_raw.FI_acc   = FI_acc_avg;
    FM_raw.MI_acc   = MI_acc_avg;
    FM_raw.FA       = FA_avg;
    FM_raw.MA       = MA_avg;
    FM_raw.Fg       = Fg_avg;
    FM_raw.FI_vel   = FI_vel_avg;
    FM_raw.MI_vel   = MI_vel_avg;
    
%     figure()
%     hold on
%     subplot(3,1,1); plot(FI_acc(1,:))
%     subplot(3,1,2); plot(FI_acc(2,:))
%     subplot(3,1,3); plot(FI_acc(3,:))
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(MI_acc(1,:))
%     subplot(3,1,2); plot(MI_acc(2,:))
%     subplot(3,1,3); plot(MI_acc(3,:))
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(FA(1,:))
%     subplot(3,1,2); plot(FA(2,:))
%     subplot(3,1,3); plot(FA(3,:))
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(MA(1,:))
%     subplot(3,1,2); plot(MA(2,:))
%     subplot(3,1,3); plot(MA(3,:))
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(Fg(1,:))
%     subplot(3,1,2); plot(Fg(2,:))
%     subplot(3,1,3); plot(Fg(3,:))
%     hold off
%     
%     pause

end


