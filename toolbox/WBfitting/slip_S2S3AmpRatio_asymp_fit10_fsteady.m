function [fitresult, gof] = slip_S2S3AmpRatio_asymp_fit10_fsteady(S2S3AmpRatioFunc, slip_wb, plot_on)
%CREATEFIT(S2S3AMPRATIOFUNC,slip_wb)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : S2S3AmpRatioFunc
%      Y Output: slip_wb
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 17-Apr-2015 10:16:38


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( S2S3AmpRatioFunc, slip_wb );

% Set up fittype and options.
ft = fittype( '-a+a/x^10', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = -Inf;
opts.StartPoint = 0.117417650855806;
opts.Upper = Inf;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
if plot_on == 1
%     figure( 'Name', 'untitled fit 1' );
    h = plot( fitresult, xData, yData );
    legend( h, 'slip_wb vs. S2S3AmpRatioFunc', 'untitled fit 1', 'Location', 'NorthEast' );
    % Label axes
    xlabel( 'S2S3AmpRatioFunc' );
    ylabel( 'slip_wb' );
    grid on
end

