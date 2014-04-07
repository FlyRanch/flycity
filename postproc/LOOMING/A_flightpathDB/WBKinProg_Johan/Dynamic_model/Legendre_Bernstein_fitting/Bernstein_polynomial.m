function [ p_x ] = Bernstein_polynomial( n_pol, n_deriv, x )

    % Create a Bernstein VanderMonde matrix:
    
    m = length(x);
    
    N = n_pol;
    
    PN = zeros(N+1, m, n_deriv+1);
    
    if n_deriv <= 4
    
        for i = 1:n_deriv+1
            
            if i == 1
                
                for j = 1:N+1
                    
                    for k = 1:m
                    
                        PN(j,k,i) = bbt(j-1,N,x(k));
                    
                    end
                    
                end
                
            elseif i == 2
                
                PN_t = zeros(N+1,m);
                
                for j = 1:N+1
                    
                    for k = 1:m
                    
                        PN_t(j,k) = bbt_1(j-1,N-1,x(k));
                        
                    end
                    
                end
                                
                PN(:,:,i) = PN_t;
                
                                
            elseif i == 3
                
                PN_t = zeros(N+1,m);
                
                for j = 1:N+1

                    for k = 1:m
                    
                        PN_t(j,k) = bbt_2(j-1,N-2,x(k));
                        
                    end
                                        
                end
                
                PN(:,:,i) = PN_t;
                
            elseif i == 4
                
                PN_t = zeros(N+1,m);
                
                for j = 1:N+1
                    
                     for k = 1:m
                    
                        PN_t(j,k) = bbt_3(j-1,N-3,x(k));
                        
                    end                   
                    
                end
                
                PN(:,:,i) = PN_t;
                
            elseif i == 5
                
                PN_t = zeros(N+1,m);
                
                for j = 1:N+1
                    
                     for k = 1:m
                    
                        PN_t(j,k) = bbt_4(j-1,N-4,x(k));
                        
                    end                        
                    
                end
                
                PN(:,:,i) = PN_t;
                
            end
            
        end
    
    else
        
        'nr of derivatives is too high'
        
    end
     
     p_x = zeros(N+1, m, n_deriv+1);
    
     p_x(:,2:end-1,:) = PN(:,2:end-1,:);


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

function[b] = bbt_1(i,n,eta)


        b = (n+1)*(bbt(i-1,n,eta)-bbt(i,n,eta));
        
    
end

function[b] = bbt_2(i,n,eta)


        b = (n+1)*(bbt_1(i-1,n,eta)-bbt_1(i,n,eta));
        
    
end

function[b] = bbt_3(i,n,eta)


        b = (n+1)*(bbt_2(i-1,n,eta)-bbt_2(i,n,eta));
        
    
end

function[b] = bbt_4(i,n,eta)


        b = (n+1)*(bbt_3(i-1,n,eta)-bbt_3(i,n,eta));
        
    
end