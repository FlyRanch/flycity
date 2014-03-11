function plot_maneuvering_wingbeats( n_pol_theta, n_pol_eta, n_pol_phi, a_avg_theta, a_avg_eta, a_avg_phi, a_dev_theta, a_dev_eta, a_dev_phi, FM_range, f, down_up, nr_wb ,plot_color)

    [ t, X_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up, 200, 0, 1/f, 0 );
    [ ~, X_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, down_up, 200, 0, 1/f, 0 );
    [ ~, X_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, down_up, 200, 0, 1/f, 0 );
    
     
    a_theta = zeros(2*(n_pol_theta+1),nr_wb);
    a_eta = zeros(2*(n_pol_eta+1),nr_wb);
    a_phi = zeros(2*(n_pol_phi+1),nr_wb);
    


    FM = [0:(FM_range/(nr_wb-1)):FM_range];
            

    
    for i = 1:nr_wb
    
        a_theta(:,i) = [a_avg_theta + a_dev_theta.*FM(i)];
        a_eta(:,i) = [a_avg_eta + a_dev_eta.*FM(i)];
        a_phi(:,i) = [a_avg_phi + a_dev_phi.*FM(i)];
    
    end

    

    hold on
    for i = 1:nr_wb
        plot(t,X_theta*a_theta(:,i),'Color',plot_color)
    end
    hold off
    

    hold on
    for i = 1:nr_wb
        plot(t,X_eta*a_eta(:,i),'Color',plot_color)
    end
    hold off
    

    hold on
    for i = 1:nr_wb
        plot(t,X_phi*a_phi(:,i),'Color',plot_color)
    end
    hold off

%     figure()
%     hold on
%     subplot(3,1,1);
%     hold on
%     for i = 1:nr_wb
%         plot(t,X_theta*a_theta(:,i),'Color',[0.5 0.5 0.5])
%     end
%     hold off
%     subplot(3,1,2);
%     hold on
%     for i = 1:nr_wb
%         plot(t,X_eta*a_eta(:,i),'Color',[0.5 0.5 0.5])
%     end
%     hold off
%     subplot(3,1,3);
%     hold on
%     for i = 1:nr_wb
%         plot(t,X_phi*a_phi(:,i),'Color',[0.5 0.5 0.5])
%     end
%     hold off    
%     hold off

end

