function axayaz_t_raw_filt(settings, seq_nr, u_f, v_f, w_f, ax_f, ay_f, az_f, t, fig_nr1, fig_nr2,save_on_off)

    % Program that generates plots which compares the raw ax,ay,az data with
    % the filtered ax,ay,az data over time.
    
    N = length(u_f);
    dt = t(2) - t(1);
       
    % Calculate acceleration from u_f, v_f and w_f
    ax_r = zeros(N,1);
    ay_r = zeros(N,1);
    az_r = zeros(N,1);
    
    for k = 1:(N-1)
    
        ax_r(k) = (u_f(k+1)-u_f(k))/(dt);
        ay_r(k) = (v_f(k+1)-v_f(k))/(dt);
        az_r(k) = (w_f(k+1)-w_f(k))/(dt);
        
    end
    
    ax_r(N) = ax_r(N-1);
    ay_r(N) = ay_r(N-1);
    az_r(N) = az_r(N-1);
    
    figure(fig_nr1)
    subplot(3,1,1); plot(t,ax_r,t,ax_f)
    title('Filtered and unfiltered acceleration over time')
    xlabel('t [s]')
    ylabel('a_x [mm/s^2]')
    subplot(3,1,2); plot(t,ay_r,t,ay_f)
    xlabel('t [s]')
    ylabel('a_y [mm/s^2]')
    subplot(3,1,3); plot(t,az_r,t,az_f)
    xlabel('t [s]')
    ylabel('a_z [mm/s^2]')

    figure(fig_nr2)
    subplot(3,1,1); plot(t,ax_r-ax_f,t,0)
    title('Deviation of unfiltered from filtered acceleration')
    xlabel('t [s]')
    ylabel('a_x [mm/s^2]')
    subplot(3,1,2); plot(t,ay_r-ay_f,t,0)
    xlabel('t [s]')
    ylabel('a_y [mm/s^2]')
    subplot(3,1,3); plot(t,az_r-az_f,t,0)
    xlabel('t [s]')
    ylabel('a_z [mm/s^2]')
    
    % Save plots
    
    if save_on_off == 1
    
    saveas(fig_nr1, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/axayaz1'], 'fig')
    
    saveas(fig_nr1, [char(settings.plot_folders(2)) '/axayaz_body/axayaz1_' int2str(seq_nr)], 'fig')
    
    saveas(fig_nr2, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/axayaz2'], 'fig')
    
    saveas(fig_nr2, [char(settings.plot_folders(2)) '/axayaz_body/axayaz2_' int2str(seq_nr)], 'fig')
    
    end

end

