function wing_vel(settings, seq_nr, u_wing_L, v_wing_L, w_wing_L, u_wing_R, v_wing_R, w_wing_R, t , fig_nr1, fig_nr2, fig_nr3, save_on_off)

    N = length(u_wing_L(1,:));
    
    pos = [1 1-(0.5/(N-2)):-(1/(N-2)):(0.5/(N-2)) 0];
    
    sect_id = ones(5,1);
    
    for i = 1:N-1
        
        if pos(i) == 1
            
            sect_id(1) = 1;
            
        elseif pos(i-1) >= 0.75 >= pos(i)
            
            if (pos(i-1)-0.75) > (0.75-pos(i))
                
                sect_id(2) = i;
                
            else
                
                sect_id(2) = i-1;
                
            end
            
        elseif pos(i-1) >= 0.5 >= pos(i)
            
            if (pos(i-1)-0.5) > (0.5-pos(i))
                
                sect_id(3) = i;
                
            else
                
                sect_id(3) = i-1;
                
            end
            
        elseif pos(i-1) >= 0.25 >= pos(i)
            
            if (pos(i-1)-0.25) > (0.25-pos(i))
                
                sect_id(4) = i;
                
            else
                
                sect_id(4) = i-1;
                
            end
            
            elseif pos(i-1) >= 0.1 >= pos(i)
            
            if (pos(i-1)-0.1) > (0.1-pos(i))
                
                sect_id(5) = i;
                
            else
                
                sect_id(5) = i-1;
                
            end
        end
    end





    U_wing_L = sqrt(u_wing_L(1:end,:).^2+v_wing_L(1:end,:).^2+w_wing_L(1:end,:).^2);
    
    U_wing_R = sqrt(u_wing_R(1:end,:).^2+v_wing_R(1:end,:).^2+w_wing_R(1:end,:).^2);

    % Plot velocities over the wing on different locations on the left and
    % right wing.
    
    figure(fig_nr1)
    subplot(4,1,1); plot(t,u_wing_L(:,sect_id(1)),t,u_wing_L(:,sect_id(2)),t,u_wing_L(:,sect_id(3)),t,u_wing_L(:,sect_id(4)),t,u_wing_L(:,sect_id(5)))
    title('Air velocities at the left wing')
    xlabel('t [s]')
    ylabel('u [mm/s]')
    legend('100%','75%','50%','25%','0%')
    subplot(4,1,2); plot(t,v_wing_L(:,sect_id(1)),t,v_wing_L(:,sect_id(2)),t,v_wing_L(:,sect_id(3)),t,v_wing_L(:,sect_id(4)),t,v_wing_L(:,sect_id(5)))
    xlabel('t [s]')
    ylabel('v [mm/s]')
    legend('100%','75%','50%','25%','0%')
    subplot(4,1,3); plot(t,w_wing_L(:,sect_id(1)),t,w_wing_L(:,sect_id(2)),t,w_wing_L(:,sect_id(3)),t,w_wing_L(:,sect_id(4)),t,w_wing_L(:,sect_id(5)))
    xlabel('t [s]')
    ylabel('w [mm/s]')
    legend('100%','75%','50%','25%','0%')
    subplot(4,1,4); plot(t,U_wing_L(:,sect_id(1))',t,U_wing_L(:,sect_id(2))',t,U_wing_L(:,sect_id(3))',t,U_wing_L(:,sect_id(4))',t,U_wing_L(:,sect_id(5))')
    xlabel('t [s]')
    ylabel('|U| [mm/s]')
    legend('100%','75%','50%','25%','0%')
    
    figure(fig_nr2)
    subplot(4,1,1); plot(t,u_wing_R(:,sect_id(1)),t,u_wing_R(:,sect_id(2)),t,u_wing_R(:,sect_id(3)),t,u_wing_R(:,sect_id(4)),t,u_wing_R(:,sect_id(5)))
    title('Air velocities at the right wing')
    xlabel('t [s]')
    ylabel('u [mm/s]')
    legend('100%','75%','50%','25%','0%')
    subplot(4,1,2); plot(t,v_wing_R(:,sect_id(1)),t,v_wing_R(:,sect_id(2)),t,v_wing_R(:,sect_id(3)),t,v_wing_R(:,sect_id(4)),t,v_wing_R(:,sect_id(5)))
    xlabel('t [s]')
    ylabel('v [mm/s]')
    legend('100%','75%','50%','25%','0%')
    subplot(4,1,3); plot(t,w_wing_R(:,sect_id(1)),t,w_wing_R(:,sect_id(2)),t,w_wing_R(:,sect_id(3)),t,w_wing_R(:,sect_id(4)),t,w_wing_R(:,sect_id(5)))
    xlabel('t [s]')
    ylabel('w [mm/s]')
    legend('100%','75%','50%','25%','0%')
    subplot(4,1,4); plot(t,U_wing_R(:,sect_id(1))',t,U_wing_R(:,sect_id(2))',t,U_wing_R(:,sect_id(3))',t,U_wing_R(:,sect_id(4))',t,U_wing_R(:,sect_id(5))')
    xlabel('t [s]')
    ylabel('|U| [mm/s]')
    legend('100%','75%','50%','25%','0%')
    
    figure(fig_nr3)
    subplot(4,1,1); plot(t,u_wing_L(:,sect_id(2)),t,u_wing_R(:,sect_id(2)))
    title('Air velocities at 75% to the tip for left and right wing.')
    xlabel('t [s]')
    ylabel('u [mm/s]')
    legend('left','right')
    subplot(4,1,2); plot(t,v_wing_L(:,sect_id(2)),t,v_wing_R(:,sect_id(2)))
    xlabel('t [s]')
    ylabel('v [mm/s]')
    legend('left','right')
    subplot(4,1,3); plot(t,w_wing_L(:,sect_id(2)),t,w_wing_R(:,sect_id(2)))
    xlabel('t [s]')
    ylabel('w [mm/s]')
    legend('left','right')
    subplot(4,1,4); plot(t,U_wing_L(:,sect_id(2))',t,U_wing_R(:,sect_id(2))')
    xlabel('t [s]')
    ylabel('|U| [mm/s]')
    legend('left','right')
    
    % Save plots
    
    if save_on_off == 1
    
    saveas(fig_nr1, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/wing_vel_L'], 'fig')
    
    saveas(fig_nr1, [char(settings.plot_folders(2)) '/wing_vel/wing_vel_L_' int2str(seq_nr)], 'fig')
    
    saveas(fig_nr2, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/wing_vel_R'], 'fig')
    
    saveas(fig_nr2, [char(settings.plot_folders(2)) '/wing_vel/wing_vel_R_' int2str(seq_nr)], 'fig')
    
    saveas(fig_nr3, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/wing_vel_75'], 'fig')
    
    saveas(fig_nr3, [char(settings.plot_folders(2)) '/wing_vel/wing_vel_75_' int2str(seq_nr)], 'fig')
    
    end


end

