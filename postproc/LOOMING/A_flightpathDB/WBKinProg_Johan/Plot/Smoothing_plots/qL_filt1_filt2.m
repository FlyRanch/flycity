function qL_filt1_filt2(settings,seq_nr,qL1_f1,qL2_f1,qL3_f1,qL4_f1,qL1_f2,qL2_f2,qL3_f2,qL4_f2,t,fig_nr,save_on_off)

    %Plots the filt1 quaternions and the filt2 quaternions.
    
    figure(fig_nr)
    subplot(4,1,1); plot(t,qL1_f1,t,qL1_f2)
    title('Left wing quaternion with filt1 and filt2')
    xlabel('t [s]')
    ylabel('W1')
    subplot(4,1,2); plot(t,qL2_f1,t,qL2_f2)
    xlabel('t [s]')
    ylabel('W2')
    subplot(4,1,3); plot(t,qL3_f1,t,qL3_f2)
    xlabel('t [s]')
    ylabel('W3')
    subplot(4,1,4); plot(t,qL4_f1,t,qL4_f2)
    xlabel('t [s]')
    ylabel('S')

    % Save plots
    
    if save_on_off == 1
    
    saveas(fig_nr, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/qL_filt1_filt2'], 'fig')
    
    saveas(fig_nr, [char(settings.plot_folders(2)) '/quat_left/qL_filt1_filt2_' int2str(seq_nr)], 'fig')
    
    end

end

