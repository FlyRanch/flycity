function [fitresult, gof, coeffs] = create_fourier_series_4thN8th_order(t_loc, y_loc,n_order,plotting)

%CREATEFIT(T_LOC,Y_LOC)
%  Create a fit.
%
%  Data for 'wb_fourier_series_8th_order' fit:
%      X Input : t_loc
%      Y Output: y_loc
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 21-May-2013 10:43:32

if n_order == 4
    
    % Fit: 'wb_fourier_series_4th_order'.
    [xData, yData] = prepareCurveData( t_loc, y_loc );

    % Set up fittype and options.
    ft = fittype( 'a0 + a1*cos(2*pi()*x) + b1*sin(2*pi()*x)  + a2*cos(2*2*pi()*x) + b2*sin(2*2*pi()*x)  + a3*cos(3*2*pi()*x) + b3*sin(3*2*pi()*x)  + a4*cos(4*2*pi()*x) + b4*sin(4*2*pi()*x)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( ft );
    opts.Algorithm = 'Levenberg-Marquardt';
    opts.Display = 'Off';
    opts.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf];
    opts.StartPoint = [0.365816176838171 0.763504640848813 0.627896379614169 0.771980385554245 0.93285357027882 0.972740854003014 0.192028349427775 0.138874202829155 0.696266337082995];
    opts.Upper = [Inf Inf Inf Inf Inf Inf Inf Inf Inf];

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );

    
elseif  n_order == 8

    % Fit: 'wb_fourier_series_8th_order'.
    [xData, yData] = prepareCurveData( t_loc, y_loc );

    % Set up fittype and options.
    ft = fittype( 'a0 + a1*cos(2*pi()*x) + b1*sin(2*pi()*x)  + a2*cos(2*2*pi()*x) + b2*sin(2*2*pi()*x)  + a3*cos(3*2*pi()*x) + b3*sin(3*2*pi()*x)  + a4*cos(4*2*pi()*x) + b4*sin(4*2*pi()*x)  + a5*cos(5*2*pi()*x) + b5*sin(5*2*pi()*x)  + a6*cos(6*2*pi()*x) + b6*sin(6*2*pi()*x)  + a7*cos(7*2*pi()*x) + b7*sin(7*2*pi()*x)  + a8*cos(8*2*pi()*x) + b8*sin(8*2*pi()*x)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( ft );
    opts.Algorithm = 'Levenberg-Marquardt';
    opts.Display = 'Off';
    opts.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf];
    opts.StartPoint = [0.647617630172684 0.679016754093202 0.635786710514084 0.945174113109401 0.208934922426023 0.709281702710545 0.236230576993797 0.119396247797306 0.607303940685635 0.450137696965896 0.458725493648868 0.661944751905652 0.77028551480366 0.350218013441105 0.662009598359135 0.416158589969796 0.841929152691309];
    opts.Upper = [Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf];

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );

end

if plotting == 1;
    % Plot fit with data.
    figure( 'Name', 'wb_fourier_series_8th_order' );
    plot( fitresult, xData, yData );
    % Label axes
    xlabel( 't_loc' );
    ylabel( 'y_loc' );
    grid on
end

coeffs = coeffvalues(fitresult);