function body_vel(settings, seq_nr, u_body, v_body, w_body, ax_body, ay_body, az_body, t , fig_nr1, fig_nr2, save_on_off)


    % Program that plots and saves the velocity and acceleration
    % experienced in the body frame of reference:
    
    U_body = sqrt(u_body.^2+v_body.^2+w_body.^2);
    
    a_body = sqrt(ax_body.^2+ay_body.^2+az_body.^2);
    
    figure(fig_nr1)
    subplot(4,1,1); plot(t,u_body)
    title('Velocity in the body frame of reference')
    xlabel('t [s]')
    ylabel('u [mm/s]')
    subplot(4,1,2); plot(t,v_body)
    xlabel('t [s]')
    ylabel('v [mm/s]')
    subplot(4,1,3); plot(t,w_body)
    xlabel('t [s]')
    ylabel('w [mm/s]')
    subplot(4,1,4); plot(t,U_body)
    xlabel('t [s]')
    ylabel('|U| [mm/s]')
    
    figure(fig_nr2)
    subplot(4,1,1); plot(t,ax_body)
    title('Acceleration in the body frame of reference')
    xlabel('t [s]')
    ylabel('a_x [mm/s^2]')
    subplot(4,1,2); plot(t,ay_body)
    xlabel('t [s]')
    ylabel('a_y [mm/s^2]')
    subplot(4,1,3); plot(t,az_body)
    xlabel('t [s]')
    ylabel('a_z [mm/s^2]')
    subplot(4,1,4); plot(t,a_body)
    xlabel('t [s]')
    ylabel('|a| [mm/s^2]')
    
    % Save plots
    
    if save_on_off == 1
    
    saveas(fig_nr1, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/body_U'], 'fig')
    
    saveas(fig_nr1, [char(settings.plot_folders(2)) '/body_vel/body_U_' int2str(seq_nr)], 'fig')
    
    saveas(fig_nr2, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/body_acc'], 'fig')
    
    saveas(fig_nr2, [char(settings.plot_folders(2)) '/body_vel/body_acc_' int2str(seq_nr)], 'fig')
    
    end




end