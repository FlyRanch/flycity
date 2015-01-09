function Plot_wing_vel_components(settings, seq_nr, u_trans_L, u_trans_R, u_body_rot_L, u_body_rot_R, u_wing_rot_L, u_wing_rot_R, t, fig_nr1, fig_nr2, fig_nr3, save_on_off)
    
    % Plot the 3 components which build up the velocity on the wing, being:
    % translational velocity, body rotation velocity and and wing rotation
    % velocity.
    
    
    N = length(u_trans_L(1,:,1));
    
    pos = [1 1-(0.5/(N-2)):-(1/(N-2)):(0.5/(N-2)) 0];

    
    for i = 2:N-1

        if pos(i-1) >= 0.75 >= pos(i)
            
            if (pos(i-1)-0.75) > (0.75-pos(i))
                
                sect_id = i;
                
            else
                
                sect_id = i-1;
                
            end
        
        end

    end
    
    
    U_trans_L = sqrt(u_trans_L(1:end,sect_id,1).^2+u_trans_L(1:end,sect_id,2).^2+u_trans_L(1:end,sect_id,3).^2);
    
    U_trans_R = sqrt(u_trans_R(1:end,sect_id,1).^2+u_trans_R(1:end,sect_id,2).^2+u_trans_R(1:end,sect_id,3).^2);
    
    U_body_rot_L = sqrt(u_body_rot_L(1:end,sect_id,1).^2+u_body_rot_L(1:end,sect_id,2).^2+u_body_rot_L(1:end,sect_id,3).^2);
    
    U_body_rot_R = sqrt(u_body_rot_R(1:end,sect_id,1).^2+u_body_rot_R(1:end,sect_id,2).^2+u_body_rot_R(1:end,sect_id,3).^2);
    
    U_wing_rot_L = sqrt(u_wing_rot_L(1:end,sect_id,1).^2+u_wing_rot_L(1:end,sect_id,2).^2+u_wing_rot_L(1:end,sect_id,3).^2);
    
    U_wing_rot_R = sqrt(u_wing_rot_R(1:end,sect_id,1).^2+u_wing_rot_R(1:end,sect_id,2).^2+u_wing_rot_R(1:end,sect_id,3).^2);

    % Plot velocities over the wing on different locations on the left and
    % right wing.
    
    figure(fig_nr1)
    plot(t,U_trans_L,t,U_trans_R)
    title('Translational air velocity at 75% span for left and right wing.')
    xlabel('t [s]')
    ylabel('|U| [mm/s]')
    legend('left','right')
    
    figure(fig_nr2)
    plot(t,U_body_rot_L,t,U_body_rot_R)
    title('Air velocity due to body rotation at 75% span for left and right wing.')
    xlabel('t [s]')
    ylabel('|U| [mm/s]')
    legend('left','right')
    
    figure(fig_nr3)
    plot(t,U_wing_rot_L,t,U_wing_rot_R)
    title('Air velocity due to wing rotation at 75% span for left and right wing.')
    xlabel('t [s]')
    ylabel('|U| [mm/s]')
    legend('left','right')
    
    % Save plots
    
    if save_on_off == 1
    
    saveas(fig_nr1, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/trans'], 'fig')
    
    saveas(fig_nr1, [char(settings.plot_folders(2)) '/trans_rot_wing/trans_' int2str(seq_nr)], 'fig')
    
    saveas(fig_nr2, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/body_rot'], 'fig')
    
    saveas(fig_nr2, [char(settings.plot_folders(2)) '/trans_rot_wing/body_rot_' int2str(seq_nr)], 'fig')
    
    saveas(fig_nr3, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/wing_rot'], 'fig')
    
    saveas(fig_nr3, [char(settings.plot_folders(2)) '/trans_rot_wing/wing_rot_' int2str(seq_nr)], 'fig')
    
    end


end

