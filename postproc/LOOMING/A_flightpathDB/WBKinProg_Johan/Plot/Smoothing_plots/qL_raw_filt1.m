function qL_raw_filt1(settings,seq_nr,qL1_r,qL2_r,qL3_r,qL4_r,qL1_f1,qL2_f1,qL3_f1,qL4_f1,t,fig_nr1, fig_nr2, save_on_off)

    % Plots the unfiltered left wing quaternion and the left wing
    % quaternion corrected for the filtered body motion.
    
    figure(fig_nr1)
    subplot(4,1,1); plot(t,qL1_r,t,qL1_f1)
    title('Left wing quaternion with raw and filtered body attitude.')
    xlabel('t [s]')
    ylabel('W1')
    subplot(4,1,2); plot(t,qL2_r,t,qL2_f1)
    xlabel('t [s]')
    ylabel('W2')
    subplot(4,1,3); plot(t,qL3_r,t,qL3_f1)
    xlabel('t [s]')
    ylabel('W3')
    subplot(4,1,4); plot(t,qL4_r,t,qL4_f1)
    xlabel('t [s]')
    ylabel('S')

    
    r_r = real(atan2(2*(qL4_r.*qL1_r+qL2_r.*qL3_r),(1-2*(qL1_r.^2+qL2_r.^2))));
    p_r = -real(asin(2*(qL4_r.*qL2_r-qL3_r.*qL1_r)));
    y_r = real(atan2(2*(qL4_r.*qL3_r+qL1_r.*qL2_r),1-2*(qL2_r.^2+qL3_r.^2)));
    
    r_f = real(atan2(2*(qL4_f1.*qL1_f1+qL2_f1.*qL3_f1),(1-2*(qL1_f1.^2+qL2_f1.^2))));
    p_f = -real(asin(2*(qL4_f1.*qL2_f1-qL3_f1.*qL1_f1)));
    y_f = real(atan2(2*(qL4_f1.*qL3_f1+qL1_f1.*qL2_f1),1-2*(qL2_f1.^2+qL3_f1.^2)));
    
    figure(fig_nr2)
    subplot(3,1,1); plot(t,radtodeg(r_r),t,radtodeg(r_f))
    title('Unfiltered and filtered wing kinematics with body mounted stroke plane for the left wing')
    xlabel('t [s]')
    ylabel('roll [deg]')
    subplot(3,1,2); plot(t,radtodeg(p_r),t,radtodeg(p_f))
    xlabel('t [s]')
    ylabel('pitch [deg]')
    subplot(3,1,3); plot(t,radtodeg(y_r),t,radtodeg(y_f))
    xlabel('t [s]')
    ylabel('yaw [deg]')
    
    % Save plots
    
    if save_on_off == 1
    
    saveas(fig_nr1, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/qL_raw_filt1_1'], 'fig')
    
    saveas(fig_nr1, [char(settings.plot_folders(2)) '/quat_left/qL_raw_filt1_1_' int2str(seq_nr)], 'fig')
    
    saveas(fig_nr2, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/qL_raw_filt1_2'], 'fig')
    
    saveas(fig_nr2, [char(settings.plot_folders(2)) '/quat_left/qL_raw_filt1_2_' int2str(seq_nr)], 'fig')
    
    end


end

