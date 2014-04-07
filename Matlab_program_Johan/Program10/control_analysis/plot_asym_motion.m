function plot_asym_motion( n_a, non_dim_matrix, fit_p, fit_n, fit_pn )

   for i = 1:n_a
        
        X_p = [];
        Y_p = [];
        
        X_n = [];
        Y_n = [];
        
        X_pn = [];
        Y_pn = [];
        
        for j = 1:length(non_dim_matrix(:,1))
            
            if abs(non_dim_matrix(j,2+i))<=1 && abs(non_dim_matrix(j,2+n_a+i))<=1
        
                X_p = [X_p; non_dim_matrix(j,1)+non_dim_matrix(j,2)];
                Y_p = [Y_p; non_dim_matrix(j,2+i)];

                X_n = [X_n; non_dim_matrix(j,1)+non_dim_matrix(j,2)];
                Y_n = [Y_n; non_dim_matrix(j,2+n_a+i)];

                X_pn = [X_pn; non_dim_matrix(j,1)+non_dim_matrix(j,2)];
                Y_pn = [Y_pn; non_dim_matrix(j,2+i)-non_dim_matrix(j,2+n_a+i)];
        
            end
        
        end

        
        figure()
        hold on
        plot(X_p,Y_p,'o','Color','r')
        plot(X_n,Y_n,'o','Color','g')
        plot(X_pn,Y_pn,'o','Color','b')
        if fit_p(i,3) >= 0.025
            plot(sort(X_p),fit_p(i,1)*sort(X_p)+fit_p(i,2),'r')
        end
        if fit_n(i,3) >= 0.025
            plot(sort(X_n),fit_n(i,1)*sort(X_n)+fit_n(i,2),'g')
        end
        if fit_pn(i,3) >= 0.025
            plot(sort(X_pn),fit_pn(i,1)*sort(X_pn)+fit_pn(i,2),'b')
        end
        hold off
        title((char([ 'coefficient a-' int2str(i)])))
        xlabel('non-dimensional roll moment')
        ylabel('deviation-coefficient')
        legend('left','right','left+right')
    
    end

end
