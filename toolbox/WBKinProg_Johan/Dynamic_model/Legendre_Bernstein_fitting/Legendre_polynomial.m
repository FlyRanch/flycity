function [ p_x ] = Legendre_polynomial( n_pol, n_deriv, x )


    % Function that returns the values of a n_pol order Legendre polynomial or its
    % n_deriv derivative for a vector x.
    
    m = length(x);
    
    PN = zeros(n_pol+1, m, n_deriv+1);
    
    
    
    for i = 1:(n_deriv+1)
        
        if i == 1
        
            for j = 1:(n_pol+1)

                if j == 1

                    PN(j,:,i) = ones(m,1);

                elseif j == 2

                    PN(j,:,i) = x;

                elseif j > 2

                    PN(j,:,i) = ((2*(j-2)+1)/(j-1))*x.*PN(j-1,:,i)-((j-2)/(j-1))*PN(j-2,:,i);

                end
            end
        
        else
            
            for j = 1:(n_pol+1)

                if j == 1

                    PN(j,:,i) = zeros(m,1);

                elseif j == 2

                    PN(j,:,i) = PN(j-1,:,i-1);
                    
                elseif j == 3
                    
                    PN(j,:,i) = 3*PN(j-1,:,i-1);
                    
                elseif j > 3
                    
                    PN(j,:,i) = (2*(j-2)+1)*PN(j-1,:,i-1)+PN(j-2,:,i);

                end
            end
        
        end
       
    end
    
    
     p_x = PN;

end

