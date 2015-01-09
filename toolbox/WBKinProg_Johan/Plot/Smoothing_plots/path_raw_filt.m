function path_raw_filt(settings,seq_nr, x_r, y_r, z_r, x_f, y_f, z_f, u_f, v_f, w_f, fig_nr, save_on_off)

    % Program that generates a plot of the trajectory of the fly based on
    % the raw position data and the filtered position data.
    
    figure(fig_nr)
    plot3(x_r,y_r,z_r,x_f,y_f,z_f)
    title('Trajectory plot (unfiltered and filtered)')
    xlabel('x global [mm]')
    ylabel('y global [mm]')
    zlabel('z global [mm]')
    legend('Unfiltered','filtered')
    hold on
    quiver3(x_f(end),y_f(end),z_f(end),u_f(end),v_f(end),w_f(end),0.001)
    hold off
    axis equal
    grid on
    
    % Save plots
    
    if save_on_off == 1
    
    saveas(fig_nr, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/flight_path'], 'fig')
    
    saveas(fig_nr, [char(settings.plot_folders(2)) '/flight_path/flight_path_' int2str(seq_nr)], 'fig')
    
    end

end

