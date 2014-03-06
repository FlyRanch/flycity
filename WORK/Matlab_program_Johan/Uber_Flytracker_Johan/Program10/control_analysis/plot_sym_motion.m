function plot_sym_motion( n_a, non_dim_matrix, fit_p, fit_n )


    for i = 1:n_a
        
        x_p = [];
        y_p = [];
        
        x_n = [];
        y_n = [];
        
        for j = 1:length(non_dim_matrix(:,1))
            if abs((non_dim_matrix(j,2+i)+non_dim_matrix(j,2+n_a+i))/2) <=1
            if (non_dim_matrix(j,1)+non_dim_matrix(j,2)) >= 0
                x_p = [x_p; non_dim_matrix(j,1)+non_dim_matrix(j,2)];
                y_p = [y_p; (non_dim_matrix(j,2+i)+non_dim_matrix(j,2+n_a+i))/2];
            elseif (non_dim_matrix(j,1)+non_dim_matrix(j,2)) < 0
                x_n = [x_n; non_dim_matrix(j,1)+non_dim_matrix(j,2)];
                y_n = [y_n; (non_dim_matrix(j,2+i)+non_dim_matrix(j,2+n_a+i))/2];
            end
            end
        end
        
        
        figure()
        hold on
        plot(x_p,y_p,'o','Color','r')
        plot(x_n,y_n,'o','Color','g')
        if fit_p(i,3) >= 0.05
        plot(sort(x_p),fit_p(i,1)*sort(x_p)+fit_n(i,2),'r')
        end
        if fit_n(i,3) >= 0.05
        plot(sort(x_n),fit_n(i,1)*sort(x_n)+fit_n(i,2),'g')
        end
        hold off
        title((char([ 'coefficient a-' int2str(i)])))
        
    
%         figure()
%         hold on
%         plot(x_p,y_p,'o','Color','r')
%         plot(x_n,y_n,'o','Color','g')
% %         if abs(fit_p(i,3)) >= 0
%             plot(sort(x_p),fit_p(i,1)*sort(x_p)+fit_p(i,2),'r')
% %         end
% %         if abs(fit_n(i,3)) >= 0
%             plot(sort(x_n),fit_n(i,1)*sort(x_n)+fit_n(i,2),'g')
% %         end
%         hold off
%         title((char([ 'coefficient a-' int2str(i)])))

    end


end

