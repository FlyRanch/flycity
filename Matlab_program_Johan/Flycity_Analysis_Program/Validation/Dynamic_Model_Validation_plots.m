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
    
    figure(1)
    hFig = figure(1);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 1200 800]);
    hold on
    subplot(3,2,1); error_plot( v_strk_sim_steady, v_strk_mov_steady, v_strk_sim_man, v_strk_mov_man, ax_man, 1, 0, 7)
    xlabel('v_x [mm/s]')
    ylabel('error v_x [mm/s]')
    ylim([0 50])
    subplot(3,2,2); error_plot( v_strk_sim_steady, v_strk_mov_steady, v_strk_sim_man, v_strk_mov_man, ax_man, 1, 1, 7)
    xlabel('v_x [mm/s]')
    ylabel('error angle v_x [deg]')
    ylim([0 180])
    subplot(3,2,3); error_plot( v_strk_sim_steady, v_strk_mov_steady, v_strk_sim_man, v_strk_mov_man, ay_man, 2, 0, 7)
    xlabel('v_y [mm/s]')
    ylabel('error v_y [mm/s]')
    ylim([0 50])
    subplot(3,2,4); error_plot( v_strk_sim_steady, v_strk_mov_steady, v_strk_sim_man, v_strk_mov_man, ay_man, 2, 1, 7)
    xlabel('v_y [mm/s]')
    ylabel('error angle v_y [deg]')
    ylim([0 180])
    subplot(3,2,5); error_plot( v_strk_sim_steady, v_strk_mov_steady, v_strk_sim_man, v_strk_mov_man, ay_man, 2, 0, 7)
    xlabel('v_z [mm/s]')
    ylabel('error v_z [mm/s]')
    ylim([0 50])
    subplot(3,2,6); error_plot( v_strk_sim_steady, v_strk_mov_steady, v_strk_sim_man, v_strk_mov_man, ay_man, 2, 1, 7)
    xlabel('v_z [mm/s]')
    ylabel('error angle v_z [deg]')
    ylim([0 180])
    hold off
    [~,h1] = suplabel('Error plot strokeplane velocity dynamic simulation', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Result_plots/Dyn_sim_validation_v.eps'] ,'-depsc2');
    saveas(hFig,[char(settings.plot_loc) '/Result_plots/Dyn_sim_validation_v'],'fig')

    % Plot angular velocity:
    
    figure(2)
    hFig = figure(2);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 1200 800]);
    hold on
    subplot(3,2,1); error_plot( w_strk_sim_steady, w_strk_mov_steady, w_strk_sim_man, w_strk_mov_man, wx_man, 1, 0, 7)
    xlabel('w_x [rad/s]')
    ylabel('error w_x [rad/s]')
    ylim([0 50])
    subplot(3,2,2); error_plot( w_strk_sim_steady, w_strk_mov_steady, w_strk_sim_man, w_strk_mov_man, wx_man, 1, 1, 7)
    xlabel('w_x [rad/s]')
    ylabel('error angle w_x [deg]')
    ylim([0 180])
    subplot(3,2,3); error_plot( w_strk_sim_steady, w_strk_mov_steady, w_strk_sim_man, w_strk_mov_man, wy_man, 2, 0, 7)
    xlabel('w_y [rad/s]')
    ylabel('error w_y [rad/s]')
    ylim([0 50])
    subplot(3,2,4); error_plot( w_strk_sim_steady, w_strk_mov_steady, w_strk_sim_man, w_strk_mov_man, wy_man, 2, 1, 7)
    xlabel('w_y [rad/s]')
    ylabel('error angle w_y [deg]')
    ylim([0 180])
    subplot(3,2,5); error_plot( w_strk_sim_steady, w_strk_mov_steady, w_strk_sim_man, w_strk_mov_man, wz_man, 3, 0, 7)
    xlabel('w_z [rad/s]')
    ylabel('error w_z [rad/s]')
    ylim([0 50])
    subplot(3,2,6); error_plot( w_strk_sim_steady, w_strk_mov_steady, w_strk_sim_man, w_strk_mov_man, wz_man, 3, 1, 7)
    xlabel('w_z [rad/s]')
    ylabel('error angle w_z [deg]')
    ylim([0 180])
    hold off
    [~,h1] = suplabel('Error plot strokeplane angular velocity dynamic simulation', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Result_plots/Dyn_sim_validation_w.eps'] ,'-depsc2');
    saveas(hFig,[char(settings.plot_loc) '/Result_plots/Dyn_sim_validation_w'],'fig')
    
    % Plot acceleration:

    figure(3)
    hFig = figure(3);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 1200 800]);
    hold on
    subplot(3,2,1); error_plot( a_strk_sim_steady, a_strk_mov_steady, a_strk_sim_man, a_strk_mov_man, ax_man, 1, 0, 9)
    xlabel('a_x [mm/s^2]')
    ylabel('error a_x [mm/s^2]')
    ylim([0 8000])
    subplot(3,2,2); error_plot( a_strk_sim_steady, a_strk_mov_steady, a_strk_sim_man, a_strk_mov_man, ax_man, 1, 1, 9)
    xlabel('a_x [mm/s^2]')
    ylabel('error angle a_x [deg]')
    ylim([0 180])
    subplot(3,2,3); error_plot( a_strk_sim_steady, a_strk_mov_steady, a_strk_sim_man, a_strk_mov_man, ay_man, 2, 0, 9)
    xlabel('a_y [mm/s^2]')
    ylabel('error a_y [mm/s^2]')
    ylim([0 8000])
    subplot(3,2,4); error_plot( a_strk_sim_steady, a_strk_mov_steady, a_strk_sim_man, a_strk_mov_man, ay_man, 2, 1, 9)
    xlabel('a_y [mm/s^2]')
    ylabel('error angle a_y [deg]')
    ylim([0 180])
    subplot(3,2,5); error_plot( a_strk_sim_steady, a_strk_mov_steady, a_strk_sim_man, a_strk_mov_man, az_man, 3, 0, 9)
    xlabel('a_z [mm/s^2]')
    ylabel('error a_z [mm/s^2]')
    ylim([0 8000])
    subplot(3,2,6); error_plot( a_strk_sim_steady, a_strk_mov_steady, a_strk_sim_man, a_strk_mov_man, az_man, 3, 1, 9)
    xlabel('a_z [mm/s^2]')
    ylabel('error angle a_z [deg]')
    ylim([0 180])
    hold off
    [~,h1] = suplabel('Error plot strokeplane acceleration dynamic simulation', 't');
    set(h1,'FontSize',10)
    print ([char(settings.plot_loc) '/Result_plots/Dyn_sim_validation_a.eps'] ,'-depsc2');
    saveas(hFig,[char(settings.plot_loc) '/Result_plots/Dyn_sim_validation_a'],'fig')
    
    % Plot angular acceleration:
    
    figure(4)
    hFig = figure(4);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 1200 800]);
    hold on
    subplot(3,2,1); error_plot( w_dot_strk_sim_steady, w_dot_strk_mov_steady, w_dot_strk_sim_man, w_dot_strk_mov_man, wx_man, 1, 0, 7)
    xlabel('w dot x [rad/s^2]')
    ylabel('error w dot x [rad/s^2]')
    ylim([0 1e4])
    subplot(3,2,2); error_plot( w_dot_strk_sim_steady, w_dot_strk_mov_steady, w_dot_strk_sim_man, w_dot_strk_mov_man, wx_man, 1, 1, 7)
    xlabel('w dot x [rad/s^2]')
    ylabel('error angle w dot x [deg]')
    ylim([0 180])
    subplot(3,2,3); error_plot( w_dot_strk_sim_steady, w_dot_strk_mov_steady, w_dot_strk_sim_man, w_dot_strk_mov_man, wy_man, 2, 0, 7)
    xlabel('w dot y [rad/s^2]')
    ylabel('error w dot y [rad/s^2]')
    ylim([0 1e4])
    subplot(3,2,4); error_plot( w_dot_strk_sim_steady, w_dot_strk_mov_steady, w_dot_strk_sim_man, w_dot_strk_mov_man, wy_man, 2, 1, 7)
    xlabel('w dot y [rad/s^2]')
    ylabel('error angle w dot y [deg]')
    ylim([0 180])
    subplot(3,2,5); error_plot( w_dot_strk_sim_steady, w_dot_strk_mov_steady, w_dot_strk_sim_man, w_dot_strk_mov_man, wz_man, 3, 0, 7)
    xlabel('w dot z [rad/s^2]')
    ylabel('error w dot z [rad/s^2]')
    ylim([0 1e4])
    subplot(3,2,6); error_plot( w_dot_strk_sim_steady, w_dot_strk_mov_steady, w_dot_strk_sim_man, w_dot_strk_mov_man, wz_man, 3, 1, 7)
    xlabel('w dot z [rad/s^2]')
    ylabel('error angle w dot z [deg]')
    ylim([0 180])
    hold off
    [~,h1] = suplabel('Error plot strokeplane angular acceleration dynamic simulation', 't');
    set(h1,'FontSize',10)    
    print ([char(settings.plot_loc) '/Result_plots/Dyn_sim_validation_w_dot.eps'] ,'-depsc2');
    saveas(hFig,[char(settings.plot_loc) '/Result_plots/Dyn_sim_validation_w_dot'],'fig')


end