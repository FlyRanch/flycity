function qR_filt1_filt2(settings,seq_nr,qR1_f1,qR2_f1,qR3_f1,qR4_f1,qR1_f2,qR2_f2,qR3_f2,qR4_f2,t,fig_nr,save_on_off)

    %Plots the filt1 quaternions and the filt2 quaternions.
    
    figure(fig_nr)
    subplot(4,1,1); plot(t,qR1_f1,t,qR1_f2)
    title('Right wing quaternion with filt1 and filt2')
    xlabel('t [s]')
    ylabel('W1')
    subplot(4,1,2); plot(t,qR2_f1,t,qR2_f2)
    xlabel('t [s]')
    ylabel('W2')
    subplot(4,1,3); plot(t,qR3_f1,t,qR3_f2)
    xlabel('t [s]')
    ylabel('W3')
    subplot(4,1,4); plot(t,qR4_f1,t,qR4_f2)
    xlabel('t [s]')
    ylabel('S')
    
    % Save plots
    
    if save_on_off == 1
    
    saveas(fig_nr, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/qR_filt1_filt2'], 'fig')
    
    saveas(fig_nr, [char(settings.plot_folders(2)) '/quat_right/qR_filt1_filt2_' int2str(seq_nr)], 'fig')
    
    end

end
