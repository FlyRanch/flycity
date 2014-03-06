function [ a_dev ] = a_dev_fit( n_pol, dev_fit_const, order_dev, cor_coeff )


    % Construct deviation fit as a function of aplied non-dimensional force
    % or non-dimensional moment:
    
    a_dev = zeros((n_pol+1)*2,1);
    
    for i = 1:length(a_dev)
        
        if dev_fit_const(i,3) >= cor_coeff
        
            if i <= order_dev

                a_dev(i) = dev_fit_const(i,1);

            elseif i >= (n_pol+2) && i <= (n_pol+1+order_dev)

                a_dev(i) = dev_fit_const(i,1);

            end
        
        end
        
    end


end

