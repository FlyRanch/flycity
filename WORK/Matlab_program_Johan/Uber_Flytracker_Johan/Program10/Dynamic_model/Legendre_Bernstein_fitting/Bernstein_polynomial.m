function [ p_x ] = Bernstein_polynomial( n_pol, x )

    % Create a Bernstein VanderMonde matrix:
    
    m = length(x);
    
    N = n_pol;
    
    PN = zeros(N+1, m);
    
               
    for j = 1:N+1
                    
        for k = 1:m

            PN(j,k) = bbt(j-1,N,x(k));

        end
                    
    end
                

    

     

    
     p_x = PN';


end


% Shape function:

function[b] = bbt(i,n,eta)

    if i >= 0
        if n>=0
            if(n-i) >= 0 
                
                b = factorial(n)/(factorial(i)*factorial(n-i)) * eta^i * (1-eta)^(n-i);
                
            else
                b = 0;
            end
        else
            b = 0;
        end
    else
        b = 0;
    end
            
end

