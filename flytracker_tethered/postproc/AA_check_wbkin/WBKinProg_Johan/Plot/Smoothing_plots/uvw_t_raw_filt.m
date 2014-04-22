function uvw_t_raw_filt(settings, seq_nr, x_r, y_r, z_r,u_f, v_f, w_f, t, fig_nr1, fig_nr2, save_on_off)

    % Program that generates plots which compares the raw u,v,w data with
    % the filtered u,v,w data over time.
    
    N = length(x_r);
    dt = t(2) - t(1);
    
    % Calculate the velocities from time and raw position:
    u_r = zeros(N,1);
    v_r = zeros(N,1);
    w_r = zeros(N,1);
    
    for k = 1:(N-1)
    
        u_r(k) = (x_r(k+1)-x_r(k))/dt;
        v_r(k) = (y_r(k+1)-y_r(k))/dt;
        w_r(k) = (z_r(k+1)-z_r(k))/dt;
        
    end
    
    u_r(N) = u_r(N-1);
    v_r(N) = v_r(N-1);
    w_r(N) = w_r(N-1);
    
           
    figure(fig_nr1)
    subplot(3,1,1); plot(t,u_r,t,u_f)
    title('Filtered and unfiltered velocity over time')
    xlabel('t [s]')
    ylabel('u [mm/s]')
    subplot(3,1,2); plot(t,v_r,t,v_f)
    xlabel('t [s]')
    ylabel('v [mm/s]')
    subplot(3,1,3); plot(t,w_r,t,w_f)
    xlabel('t [s]')
    ylabel('w [mm/s]')

    figure(fig_nr2)
    subplot(3,1,1); plot(t,u_r-u_f,t,0)
    title('Deviation of unfiltered from filtered velocity')
    xlabel('t [s]')
    ylabel('u [mm/s]')
    subplot(3,1,2); plot(t,v_r-v_f,t,0)
    xlabel('t [s]')
    ylabel('v [mm/s]')
    subplot(3,1,3); plot(t,w_r-w_f,t,0)
    xlabel('t [s]')
    ylabel('w [mm/s]')
    
    % Save plots
    
    if save_on_off == 1
    
    saveas(fig_nr1, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/uvw1'], 'fig')
    
    saveas(fig_nr1, [char(settings.plot_folders(2)) '/uvw_body/uvw1_' int2str(seq_nr)], 'fig')
    
    saveas(fig_nr2, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/uvw2'], 'fig')
    
    saveas(fig_nr2, [char(settings.plot_folders(2)) '/uvw_body/uvw2_' int2str(seq_nr)], 'fig')
    
    end

end

