function [fitresult, gof] = freq_S2S3AmpRatio_parab_fit_fsteady(S2S3AmpRatioFunc, f_wb, plot_on)
%CREATEFIT(S2S3AMPRATIOFUNC,F_WB)
%  Create a fit.
%
%  Data for 'freq_parab_fit_fsteady' fit:
%      X Input : S2S3AmpRatioFunc
%      Y Output: f_wb
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 16-Apr-2015 17:34:40


%% Fit: 'freq_parab_fit_fsteady'.
[xData, yData] = prepareCurveData( S2S3AmpRatioFunc, f_wb );

% Set up fittype and options.
ft = fittype( 'a*x^2+b*x+188.7-a-b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf];
opts.StartPoint = [0.306349472016557 0.508508655381127];
opts.Upper = [Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
if plot_on == 1
%     figure( 'Name', 'freq_parab_fit_fsteady' );
    h = plot( fitresult, xData, yData );
    legend( h, 'f_wb vs. S2S3AmpRatioFunc', 'freq_parab_fit_fsteady', 'Location', 'NorthEast' );
    % Label axes
    xlabel( 'S2S3AmpRatioFunc' );
    ylabel( 'f_wb' );
    grid on
end

