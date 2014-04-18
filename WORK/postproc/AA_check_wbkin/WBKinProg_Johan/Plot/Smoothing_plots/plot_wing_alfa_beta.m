function plot_wing_alfa_beta(settings, seq_nr,alfa_L, beta_L, alfa_R, beta_R, t, fig_nr1, fig_nr2 , fig_nr3, save_on_off)

    % Plot the angles alfa & beta at several locations on the wing.
    
    % Plot at ~ 10%, 25%, 50%, 75% and 100%:
    
    N = length(alfa_L(1,:));
    
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
    
    figure(fig_nr1)
    subplot(2,1,1); plot(t,radtodeg(alfa_L(:,sect_id(1))),t,radtodeg(alfa_L(:,sect_id(2))),t,radtodeg(alfa_L(:,sect_id(3))),t,radtodeg(alfa_L(:,sect_id(4))),t,radtodeg(alfa_L(:,sect_id(5))))
    title('Angle of attack and side-slip angle at the left wing')
    xlabel('t [s]')
    ylabel('alfa [deg]')
    legend('100%','75%','50%','25%','10%')
    subplot(2,1,2); plot(t,radtodeg(beta_L(:,sect_id(1))),t,radtodeg(beta_L(:,sect_id(2))),t,radtodeg(beta_L(:,sect_id(3))),t,radtodeg(beta_L(:,sect_id(4))),t,radtodeg(beta_L(:,sect_id(5))))
    xlabel('t [s]')
    ylabel('beta [deg]')
    legend('100%','75%','50%','25%','10%')
    
    figure(fig_nr2)
    subplot(2,1,1); plot(t,radtodeg(alfa_R(:,sect_id(1))),t,radtodeg(alfa_R(:,sect_id(2))),t,radtodeg(alfa_R(:,sect_id(3))),t,radtodeg(alfa_R(:,sect_id(4))),t,radtodeg(alfa_R(:,sect_id(5))))
    title('Angle of attack and side-slip angle at the right wing')
    xlabel('t [s]')
    ylabel('alfa [deg]')
    legend('100%','75%','50%','25%','10%')
    subplot(2,1,2); plot(t,radtodeg(beta_R(:,sect_id(1))),t,radtodeg(beta_R(:,sect_id(2))),t,radtodeg(beta_R(:,sect_id(3))),t,radtodeg(beta_R(:,sect_id(4))),t,radtodeg(beta_R(:,sect_id(5))))
    xlabel('t [s]')
    ylabel('beta [deg]')
    legend('100%','75%','50%','25%','10%')
    
    figure(fig_nr3)
    subplot(2,1,1); plot(t,radtodeg(alfa_L(:,sect_id(2))),t,radtodeg(alfa_R(:,sect_id(2))))
    title('Angle of attack and side-slip angle at 75% wingspan at the left and right wing')
    xlabel('t [s]')
    ylabel('alfa [deg]')
    legend('left','right')
    subplot(2,1,2); plot(t,radtodeg(beta_L(:,sect_id(2))),t,radtodeg(beta_R(:,sect_id(2))))
    xlabel('t [s]')
    ylabel('beta [deg]')
    legend('left','right')
    
    
    if save_on_off == 1
    
    saveas(fig_nr1, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/alfa_left'], 'fig')
    
    saveas(fig_nr1, [char(settings.plot_folders(2)) '/alfa_beta_wing/alfa_left_' int2str(seq_nr)], 'fig')
    
    saveas(fig_nr2, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/alfa_right'], 'fig')
    
    saveas(fig_nr2, [char(settings.plot_folders(2)) '/alfa_beta_wing/alfa_right_' int2str(seq_nr)], 'fig')
    
    saveas(fig_nr3, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/alfa_lr_75'], 'fig')
    
    saveas(fig_nr3, [char(settings.plot_folders(2)) '/alfa_beta_wing/alfa_lr_75_' int2str(seq_nr)], 'fig')
    
    end

end

