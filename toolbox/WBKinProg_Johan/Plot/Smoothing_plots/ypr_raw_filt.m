function ypr_raw_filt(settings, seq_nr, qb1_r, qb2_r, qb3_r, qb4_r, qb1_f, qb2_f, qb3_f, qb4_f, t, fig_nr, save_on_off)
    
    %Plot roll, pitch and yaw in global reference frame of the flies body.
        
    % Calculate the yaw, pitch and roll
    
    r_r = atan2(-2*(qb4_r.*qb1_r+qb2_r.*qb3_r),-(1-2*(qb1_r.^2+qb2_r.^2)));
    p_r = -asin(2*(qb4_r.*qb2_r-qb3_r.*qb1_r));
    y_r = atan2(2*(qb4_r.*qb3_r+qb1_r.*qb2_r),1-2*(qb2_r.^2+qb3_r.^2));
    
    r_f = atan2(-2*(qb4_f.*qb1_f+qb2_f.*qb3_f),-(1-2*(qb1_f.^2+qb2_f.^2)));
    p_f = -asin(2*(qb4_f.*qb2_f-qb3_f.*qb1_f));
    y_f = atan2(2*(qb4_f.*qb3_f+qb1_f.*qb2_f),1-2*(qb2_f.^2+qb3_f.^2));
    
%     r_r = atan2(2*(qb4_r.*qb1_r+qb2_r.*qb3_r),(1-2*(qb1_r.^2+qb2_r.^2)));
%     p_r = asin(2*(qb4_r.*qb2_r-qb3_r.*qb1_r));
%     y_r = atan2(2*(qb4_r.*qb3_r+qb1_r.*qb2_r),1-2*(qb2_r.^2+qb3_r.^2));
%     
%     r_f = atan2(2*(qb4_f.*qb1_f+qb2_f.*qb3_f),(1-2*(qb1_f.^2+qb2_f.^2)));
%     p_f = asin(2*(qb4_f.*qb2_f-qb3_f.*qb1_f));
%     y_f = atan2(2*(qb4_f.*qb3_f+qb1_f.*qb2_f),1-2*(qb2_f.^2+qb3_f.^2));
    
    figure(fig_nr)
    subplot(3,1,1); plot(t,rad2deg(r_r),t,rad2deg(r_f))
    title('Filtered and unfiltered yaw, pitch and roll in global ref.')
    xlabel('t [s]')
    ylabel('roll [rad]')
    subplot(3,1,2); plot(t,rad2deg(p_r),t,rad2deg(p_f))
    xlabel('t [s]')
    ylabel('pitch [rad]')
    subplot(3,1,3); plot(t,rad2deg(y_r),t,rad2deg(y_f))
    xlabel('t [s]')
    ylabel('yaw [rad]')
    
    % Save plots
    
    if save_on_off == 1
    
    saveas(fig_nr, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/ypr_body'], 'fig')
    
    saveas(fig_nr, [char(settings.plot_folders(2)) '/ypr_body/ypr_body_' int2str(seq_nr)], 'fig')
    
    end

end

