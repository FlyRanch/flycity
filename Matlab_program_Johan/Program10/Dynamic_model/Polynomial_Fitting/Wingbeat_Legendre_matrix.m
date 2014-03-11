function [ t, X ] = Wingbeat_Legendre_matrix( n_pol, down_up, data_points, t_start, t_end, n_deriv )


    % Function that combines part 1 and part 2 of a wingbeat to one output
    % data vector and the time vector with uniform Delta t over the
    % wingbeat.


    
    delta_t = (t_end-t_start)/(data_points-1);
    
    t = t_start:delta_t:t_end;
    
    m1 = round(down_up*data_points);
    
    m2 = data_points-m1+1;
    
    v1 = -1:(2/(m1-1)):1;
    
    v2 = -1:(2/(m2-1)):1;
    
    PN1 = Legendre_polynomial(n_pol,n_deriv,v1);
    
    PN2 = Legendre_polynomial(n_pol,n_deriv,v2);
    
    X = zeros(data_points,(n_pol+1)*2);
    
    X(1:(m1-1),1:(n_pol+1)) = PN1(:,1:(m1-1),n_deriv+1)';
    
    X(m1,:) = [0.5*PN1(:,m1,n_deriv+1)' 0.5*PN2(:,1,n_deriv+1)'];
           
    X((m1+1):(m1+m2-1),(n_pol+2):(2*(n_pol+1))) = PN2(:,2:m2,n_deriv+1)';
    
end

