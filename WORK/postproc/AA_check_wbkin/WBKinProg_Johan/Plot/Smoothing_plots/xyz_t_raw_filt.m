function xyz_t_raw_filt(settings, seq_nr, x_r, y_r, z_r, x_f, y_f, z_f, t, fig_nr1, fig_nr2, save_on_off)

    % Program that generates plots which compares the raw x,y,z data with
    % the filtered x,y,z data over time.
    
    figure(fig_nr1)
    subplot(3,1,1); plot(t,x_r,t,x_f)
    title('Filtered and unfiltered position over time')
    xlabel('t [s]')
    ylabel('x [mm]')
    subplot(3,1,2); plot(t,y_r,t,y_f)
    xlabel('t [s]')
    ylabel('y [mm]')
    subplot(3,1,3); plot(t,z_r,t,z_f)
    xlabel('t [s]')
    ylabel('z [mm]')


    figure(fig_nr2)
    subplot(3,1,1); plot(t,x_r-x_f,t,0)
    title('Deviation of unfiltered from filtered position')
    xlabel('t [s]')
    ylabel('x [mm]')
    subplot(3,1,2); plot(t,y_r-y_f,t,0)
    xlabel('t [s]')
    ylabel('y [mm]')
    subplot(3,1,3); plot(t,z_r-z_f,t,0)
    xlabel('t [s]')
    ylabel('z [mm]')
    
    
    % Save plots
    
    if save_on_off == 1
    
    saveas(fig_nr1, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/xyz1'], 'fig')
    
    saveas(fig_nr1, [char(settings.plot_folders(2)) '/xyz_body/xyz1_' int2str(seq_nr)], 'fig')
    
    saveas(fig_nr2, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/xyz2'], 'fig')
    
    saveas(fig_nr2, [char(settings.plot_folders(2)) '/xyz_body/xyz2_' int2str(seq_nr)], 'fig')
    
    end
    

end

