function [ a_periodic ] = periodic_bc_func( a_fit, down_up )

    n_pol = (size(a_fit,1)-2)/2;
    
    
    % create 11 subsequent
    
    a_fit2 = [ a_fit a_fit a_fit a_fit a_fit a_fit a_fit a_fit a_fit a_fit a_fit ...
       a_fit a_fit a_fit a_fit a_fit a_fit a_fit a_fit a_fit a_fit ];
    
    % Weights:
    
    W_1       = ones(1,21)*(1/21);
    W_2       = ones(1,21)*(1/21);
    

    % Use function average_fit to obtain the average for the left and
    % right wing combined:

    [a_periodic_1, a_periodic_2] = average_fit(   a_fit2((1:(n_pol+1)),:) , ...
                                                        a_fit2((n_pol+2):(2*(n_pol+1)),:), ...
                                                         n_pol,W_1,W_2,down_up);
                                                     
    a_periodic = [ a_periodic_1; a_periodic_2];


end

