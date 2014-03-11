function [ b ] = lin_fit( a_dev, FM )

    n_pol = (size(a_dev,1)-2)/2;
    nr_wb = size(a_dev,2);

    % Linear fit:
    
    b_1 = zeros(n_pol*2+2,1);
    b_2 = zeros(n_pol*2+2,1);
    b_3 = zeros(n_pol*2+2,1);
    
    if isempty(FM) == 0 && isempty(a_dev) == 0
    
    for i = 1:(n_pol*2+2)
        
        X_p = FM;
        Y_p = a_dev(i,:);

        if nr_wb > 5
            
            p_t1 = polyfit(X_p,Y_p,1);
            b_1(i) = p_t1(1);
            b_2(i) = p_t1(2);
            Y_fit  = polyval(p_t1,X_p);
            mu_x   = mean(X_p);
            sx     = std(X_p);
            sy     = std(Y_p);
            b_3(i) = sum((X_p-mu_x).*(Y_p-Y_fit))/((nr_wb-1)*sx*sy);
            
        end
    
    end
    
    end
    
    b = [b_1 b_2 b_3];

end

