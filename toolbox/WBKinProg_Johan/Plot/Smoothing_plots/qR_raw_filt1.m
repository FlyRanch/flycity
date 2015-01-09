function qR_raw_filt1(settings,seq_nr,qR1_r,qR2_r,qR3_r,qR4_r,qR1_f1,qR2_f1,qR3_f1,qR4_f1,t,fig_nr,save_on_off)

    % Plots the unfiltered left wing quaternion and the left wing
    % quaternion corrected for the filtered body motion.
    
    figure(fig_nr)
    subplot(4,1,1); plot(t,qR1_r,t,qR1_f1)
    title('Right wing quaternion with unfiltered and filtered body attitude.')
    xlabel('t [s]')
    ylabel('W1')
    subplot(4,1,2); plot(t,qR2_r,t,qR2_f1)
    xlabel('t [s]')
    ylabel('W2')
    subplot(4,1,3); plot(t,qR3_r,t,qR3_f1)
    xlabel('t [s]')
    ylabel('W3')
    subplot(4,1,4); plot(t,qR4_r,t,qR4_f1)
    xlabel('t [s]')
    ylabel('S')
    
    
    % Save plots
    
    if save_on_off == 1
    
    saveas(fig_nr, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/qR_raw_filt1'], 'fig')
    
    saveas(fig_nr, [char(settings.plot_folders(2)) '/quat_right/qR_raw_filt1_' int2str(seq_nr)], 'fig')
    
    end


end