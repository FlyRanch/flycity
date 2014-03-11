function [y] = calc_val_fourier_series_4thN8th_order(x,coeffs,plotting)

% set coeffs
n_order = length(coeffs);
for n = 1:n_order
    if n < ceil(n_order/2)+1
        
        eval(['a',int2str(n-1), ' = coeffs(n);']);
        
    else
        
        eval(['b',int2str(n-ceil(n_order/2)), ' = coeffs(n);']);
        
    end
end
    
% calc values
if n_order == 9 % 4TH ORDER
    
    y = a0 + a1*cos(2*pi()*x) + b1*sin(2*pi()*x)  + a2*cos(2*2*pi()*x) + b2*sin(2*2*pi()*x)  + a3*cos(3*2*pi()*x) + b3*sin(3*2*pi()*x)  + a4*cos(4*2*pi()*x) + b4*sin(4*2*pi()*x);

elseif  n_order == 17   % 8TH ORDER

    y = a0 + a1*cos(2*pi()*x) + b1*sin(2*pi()*x)  + a2*cos(2*2*pi()*x) + b2*sin(2*2*pi()*x)  + a3*cos(3*2*pi()*x) + b3*sin(3*2*pi()*x)  + a4*cos(4*2*pi()*x) + b4*sin(4*2*pi()*x)  + a5*cos(5*2*pi()*x) + b5*sin(5*2*pi()*x)  + a6*cos(6*2*pi()*x) + b6*sin(6*2*pi()*x)  + a7*cos(7*2*pi()*x) + b7*sin(7*2*pi()*x)  + a8*cos(8*2*pi()*x) + b8*sin(8*2*pi()*x);

else
    y = 'WRONG ORDER'
end

% plot data
if exist('plotting') == 1
    if plotting == 1;
        % Plot fit values.
        plot(x,y,'-r','linewidth',1)
%         plot(x,y);
    %     % Label axes
    %     xlabel( 't' );
    %     ylabel( 'y' );
    %     grid on
    end
end
