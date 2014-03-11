function [fitresult, gof] = plot_fourier_series(fitresult)


    % Plot fit with data.
    figure( 'Name', 'wb_fourier_series_8th_order' );
    plot( fitresult, 'r','linewidth',2);
    % Label axes
    xlabel( 't_loc' );
    ylabel( 'y_loc' );
    grid on

