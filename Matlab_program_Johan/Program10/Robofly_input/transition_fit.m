function [ a_trans ] = transition_fit( a_dev, n_pol, down_up )


    % Create the deviation function of a transition wingbeat:
    
    data_points = 101;
    
    [ ~, X ] = Wingbeat_Legendre_matrix( n_pol, down_up, data_points, 0, 1, 0 );
    
    y = X*a_dev;
    
    y_ramp = [0:(1/(data_points-1)):1]';
    
    Y = y.*y_ramp;
    
    down_up_loc = round(down_up*data_points);
    
    wb_loc = [1 down_up_loc; ...
              down_up_loc data_points];
          
    a_fit = Piecewise_polynomial_fit(Y,n_pol,wb_loc);
    
    a_trans = [a_fit(:,1); a_fit(:,2)];


end

