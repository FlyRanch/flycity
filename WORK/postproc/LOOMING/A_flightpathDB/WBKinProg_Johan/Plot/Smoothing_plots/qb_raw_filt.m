function qb_raw_filt(settings, seq_nr, q1_r, q2_r, q3_r, q4_r, q1_f, q2_f, q3_f, q4_f, t, fig_nr1 , fig_nr2, save_on_off)

    % Program that plots the measured quaternions and the filtered
    % quaternions of the body.
    
    figure(fig_nr1)
    subplot(4,1,1); plot(t,q1_r,t,q1_f)
    title('Filtered and unfiltered quaternions over time')
    xlabel('t [s]')
    ylabel('Wx')
    subplot(4,1,2); plot(t,q2_r,t,q2_f)
    xlabel('t [s]')
    ylabel('Wy')
    subplot(4,1,3); plot(t,q3_r,t,q3_f)
    xlabel('t [s]')
    ylabel('Wz')
    subplot(4,1,4); plot(t,q4_r,t,q4_f)
    xlabel('t [s]')
    ylabel('S')
    
        
    figure(fig_nr2)
    subplot(4,1,1); plot(t,q1_r-q1_f,t,0)
    title('Deviation filtered and unfiltered quaternions over time')
    xlabel('t [s]')
    ylabel('Wx')
    subplot(4,1,2); plot(t,q2_r-q2_f,t,0)
    xlabel('t [s]')
    ylabel('Wy')
    subplot(4,1,3); plot(t,q3_r-q3_f,t,0)
    xlabel('t [s]')
    ylabel('Wz')
    subplot(4,1,4); plot(t,q4_r-q4_f,t,0)
    xlabel('t [s]')
    ylabel('S')
    
    % Save plots
    
    if save_on_off == 1
    
    saveas(fig_nr1, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/q_raw1'], 'fig')
    
    saveas(fig_nr1, [char(settings.plot_folders(2)) '/quat_body/q_raw1_' int2str(seq_nr)], 'fig')
    
    saveas(fig_nr2, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/q_raw2'], 'fig')
    
    saveas(fig_nr2, [char(settings.plot_folders(2)) '/quat_body/q_raw2_' int2str(seq_nr)], 'fig')
    
    end


end

