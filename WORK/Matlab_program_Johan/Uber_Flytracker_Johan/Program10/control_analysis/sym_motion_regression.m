function [ dev_fit_p, dev_fit_n ] = sym_motion_regression( n_a, non_dim_matrix)


    n = length(non_dim_matrix(:,1));

    dev_fit_p = zeros(n_a,3);
    dev_fit_n = zeros(n_a,3);

    for i = 1:n_a
        
        X_p = [];
        Y_p = [];
        
        X_n = []; 
        Y_n = [];
        
        for j = 1:n
            if abs((non_dim_matrix(j,2+i)+non_dim_matrix(j,2+n_a+i))/2) <=1
            if (non_dim_matrix(j,1)+non_dim_matrix(j,2)) >= 0
                X_p = [X_p; non_dim_matrix(j,1)+non_dim_matrix(j,2)];
                Y_p = [Y_p; (non_dim_matrix(j,2+i)+non_dim_matrix(j,2+n_a+i))/2];
            elseif (non_dim_matrix(j,1)+non_dim_matrix(j,2)) < 0
                X_n = [X_n; non_dim_matrix(j,1)+non_dim_matrix(j,2)];
                Y_n = [Y_n; (non_dim_matrix(j,2+i)+non_dim_matrix(j,2+n_a+i))/2];
            end
            end
        end
        
        b_1 = 0;
        b_2 = 0;
        b_3 = 0;
        
        b_4 = 0;
        b_5 = 0;
        b_6 = 0;
        

        
%         if length(X_p) > 5
%             [b_1,b_2,b_3,~,~]=lsqfitgm_adjust(X_p,Y_p,0,0);
%         end
%         
%         if length(X_n) > 5
%             [b_4,b_5,b_6,~,~]=lsqfitgm_adjust(X_n,Y_n,0,0);
%         end
        
        if length(X_p) > 5
            p_t1 = polyfit(X_p,Y_p,1);
            b_1 = p_t1(1);
            b_2 = p_t1(2);
            Yp1 = polyval(p_t1,X_p);
            Ymean1 = mean(Y_p);
            b_3 = sum((Yp1-Ymean1).^2)/sum((Y_p-Ymean1).^2);
        end
        
        if length(X_n) > 5
            p_t2 = polyfit(X_n,Y_n,1);
            b_4 = p_t2(1);
            b_5 = p_t2(2);
            Yp2 = polyval(p_t2,X_n);
            Ymean2 = mean(Y_n);
            b_6 = sum((Yp2-Ymean2).^2)/sum((Y_n-Ymean2).^2);
        end
        
        
        dev_fit_p(i,:) = [b_1 b_2 b_3 ];
        dev_fit_n(i,:) = [b_4 b_5 b_6 ];
        
    end

end

