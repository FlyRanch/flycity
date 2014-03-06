function [ a_dev_wb ] = a_man_wb_FM( n_pol, a_dev_man, down_up )


    % Construct a maneuvering wingbeat for a certain set of forces and
    % moments on the principal axes of the stroke plane reference frame:
    
    a_1 = [a_dev_man(1:(n_pol+1))];
    
    a_2 = [a_dev_man((n_pol+2):(2*(n_pol+1)))];
    
    Weight1 = 1;
    
    Weight2 = 1;
    
    [ a_fit_1, a_fit_2] = average_fit(a_1,a_2,n_pol,Weight1,Weight2,down_up);
    
    a_dev_wb = [a_fit_1; a_fit_2];
    
%     [ t, X ] = Wingbeat_Legendre_matrix( n_pol, down_up, 100, 0, 1, 0 )
% 
%     
%     figure()
%     plot(t,X*a_dev_wb)
   
   
end

