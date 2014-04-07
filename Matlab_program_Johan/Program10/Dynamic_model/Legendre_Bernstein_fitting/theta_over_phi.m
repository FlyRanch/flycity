function theta_over_phi( settings, pathDB )


    % plot downstroke and upstroke wingbeats seperately over a standardized
    % phi.

    dt = pathDB.t(2)-pathDB.t(1);
    
    nr_of_seq = size(pathDB.x,2);
    
    n_pol_theta = 12;
    n_pol_eta = 14;
    n_pol_phi = 10;
    
    for i =1:nr_of_seq

        seq_nr = i;

        start = find(isnan(pathDB.x(:,seq_nr))==0, 1 );
        stop = find(isnan(pathDB.x(:,seq_nr))==0, 1, 'last' );
        
        t = pathDB.t(start:stop);

        % Analyze the wing kinematics on the spherical plane:

        down_time = pathDB.down_time_L(:,:,seq_nr);
        up_time = pathDB.up_time_L(:,:,seq_nr);

        nr_wb = find(isnan(down_time(:,1))==0, 1, 'last' );

        theta_L = pathDB.theta_L(start:stop,seq_nr);
        eta_L = pathDB.eta_L(start:stop,seq_nr);
        phi_L = pathDB.phi_L(start:stop,seq_nr);

        theta_R = pathDB.theta_R(start:stop,seq_nr);
        eta_R = pathDB.eta_R(start:stop,seq_nr);
        phi_R = pathDB.phi_R(start:stop,seq_nr);

        figure()
        hold on
        for j = 1:nr_wb
            
            down_time_end = find(isnan(down_time(j,:))==0, 1, 'last' );
            
            plot(phi_L(down_time(j,1:down_time_end)),theta_L(down_time(j,1:down_time_end)),'Color',[0.5 0.5 0.5])
            
        end
        hold off
        
        figure()
        hold on
        for j = 1:nr_wb
       
            up_time_end =  find(isnan(up_time(j,:))==0, 1, 'last' );
            
            plot(phi_L(up_time(j,1:up_time_end)),theta_L(up_time(j,1:up_time_end)),'Color',[0.5 0.5 0.5])
            
        end
        hold off
        
        
        
        [a_fit,a_avg,f_avg,down_up,trigger_wb,ratio_1_2,ratio_1_2_avg] = Standard_wingbeat( settings, pathDB, i, n_pol_theta, n_pol_eta, n_pol_phi );
        
        figure()
        hold on
        for j = 1:nr_wb
            
            up_time_end =  find(isnan(up_time(j,:))==0, 1, 'last' );
            
            [ t_wb, X_wb ] = Wingbeat_Legendre_matrix( n_pol_theta, ratio_1_2(j), 100, t(down_time(j,1)), t(up_time(j,up_time_end)+1), 0 );
            
            plot(t_wb,X_wb*[a_fit.theta_L1(:,j); a_fit.theta_L2(:,j)])
            
        end
        hold off
        
        figure()
        hold on
        for j = 1:nr_wb
            
            up_time_end =  find(isnan(up_time(j,:))==0, 1, 'last' );
            
            [ t_wb, X_wb_dot ] = Wingbeat_Legendre_matrix( n_pol_theta, ratio_1_2(j), 100, t(down_time(j,1)), t(up_time(j,up_time_end)+1), 1 );
            
            plot(t_wb,X_wb_dot*[a_fit.theta_L1(:,j); a_fit.theta_L2(:,j)])
            
        end
        hold off
        
        
        [theta_pks_pos, theta_pks_pos_id] = findpeaks(theta_L,'minpeakdistance',12);
        
        [theta_pks_neg, theta_pks_neg_id] = findpeaks(-theta_L,'minpeakdistance',12);        
        
        [eta_pks_pos, eta_pks_pos_id] = findpeaks(eta_L,'minpeakdistance',12);
        
        [eta_pks_neg, eta_pks_neg_id] = findpeaks(-eta_L,'minpeakdistance',12);
        
        
        
        % Find intersections of eta with eta = 90 deg.
        
        eta_90_id = 0;
        
        for j = 2:(stop-start+1)
            
            if eta_L(j)>= pi/2 && eta_L(j-1) <= pi/2
                
                if eta_90_id == 0
                    
                    if abs((pi/2)-eta_L(j)) >= abs((pi/2)-eta_L(j-1))
                        
                        eta_90_id = j;
                        
                    elseif abs((pi/2)-eta_L(j)) < abs((pi/2)-eta_L(j-1))
                        
                        eta_90_id = j-1;
                        
                    end
                    
                else
                    
                    if abs((pi/2)-eta_L(j)) >= abs((pi/2)-eta_L(j-1))
                        
                        eta_90_id = [eta_90_id; j];
                        
                    elseif abs((pi/2)-eta_L(j)) < abs((pi/2)-eta_L(j-1))
                        
                        eta_90_id = [eta_90_id; j-1];
                        
                    end
                    
                end
            
            end
            
            if eta_L(j) < pi/2 && eta_L(j-1) >= pi/2
                
                if eta_90_id == 0
                    
                    if abs((pi/2)-eta_L(j)) >= abs((pi/2)-eta_L(j-1))
                        
                        eta_90_id = j;
                        
                    elseif abs((pi/2)-eta_L(j)) < abs((pi/2)-eta_L(j-1))
                        
                        eta_90_id = j-1;
                        
                    end
                    
                else
                    
                    if abs((pi/2)-eta_L(j)) >= abs((pi/2)-eta_L(j-1))
                        
                        eta_90_id = [eta_90_id; j];
                        
                    elseif abs((pi/2)-eta_L(j)) < abs((pi/2)-eta_L(j-1))
                        
                        eta_90_id = [eta_90_id; j-1];
                        
                    end
                    
                end
            
            end
            
        end

        eta_90_id
        
        figure()
        plot(t,eta_L)
        hold on
        plot(t(eta_pks_pos_id), eta_L(eta_pks_pos_id),'o','Color','r')
        plot(t(eta_pks_neg_id), eta_L(eta_pks_neg_id),'o','Color','r')
        plot(t(eta_90_id),eta_L(eta_90_id),'o','Color','g')
        hold off
        
        figure()
        hold on
        for j = 1:nr_wb
       
            down_time_end = find(isnan(down_time(j,:))==0, 1, 'last' );
            
            up_time_end =  find(isnan(up_time(j,:))==0, 1, 'last' );
            
            plot(phi_L(down_time(j,1:down_time_end)),theta_L(down_time(j,1:down_time_end)),'Color',[0.5 0.5 0.5])
            
            plot(phi_L(up_time(j,1:up_time_end)),theta_L(up_time(j,1:up_time_end)),'Color',[0.5 0.5 0.5])
            
            
        end
        plot(phi_L(eta_90_id),theta_L(eta_90_id),'o','Color','r') 
        hold off 
        
%         figure()
%         hold on
%         for j = 1:nr_wb
%        
%             down_time_end = find(isnan(down_time(j,:))==0, 1, 'last' );
%             
%             up_time_end =  find(isnan(up_time(j,:))==0, 1, 'last' );
%             
%             plot(phi_L(down_time(j,1:down_time_end)),theta_L(down_time(j,1:down_time_end)),'Color',[0.5 0.5 0.5])
%             
%             plot(phi_L(up_time(j,1:up_time_end)),theta_L(up_time(j,1:up_time_end)),'Color',[0.5 0.5 0.5])
%             
%             
%         end
%         plot(phi_L(theta_pks_pos_id),theta_L(theta_pks_pos_id),'o','Color','r') 
%         plot(phi_L(theta_pks_neg_id),theta_L(thteta_pks_neg_id),'o','Color','b')
%         hold off      


%         figure()
%         hold on
%         for j = 1:nr_wb
%        
%             down_time_end = find(isnan(down_time(j,:))==0, 1, 'last' );
%             
%             up_time_end =  find(isnan(up_time(j,:))==0, 1, 'last' );
%             
%             plot(phi_L(down_time(j,1:down_time_end)),theta_L(down_time(j,1:down_time_end)),'Color',[0.5 0.5 0.5])
%             
%             plot(phi_L(up_time(j,1:up_time_end)),theta_L(up_time(j,1:up_time_end)),'Color',[0.5 0.5 0.5])
%             
%             
%         end
%         plot(phi_L(eta_pks_pos_id),theta_L(eta_pks_pos_id),'o','Color','r') 
%         plot(phi_L(eta_pks_neg_id),theta_L(eta_pks_neg_id),'o','Color','b')
%         hold off      

        
        
%         figure()
%         hold on
%         for j = 1:nr_wb
%        
%             down_time_end = find(isnan(down_time(j,:))==0, 1, 'last' );
%             
%             up_time_end =  find(isnan(up_time(j,:))==0, 1, 'last' );
%             
%         
%             [~, pks_pos_id_down] = findpeaks(theta_L(down_time(j,1:down_time_end)),'minpeakdistance',15);
% 
%             [~, pks_neg_id_down] = findpeaks(-theta_L(down_time(j,1:down_time_end)),'minpeakdistance',15) ;
% 
%             [~, pks_pos_id_up] = findpeaks(theta_L(up_time(j,1:up_time_end)),'minpeakdistance',15);
% 
%             [~, pks_neg_id_up] = findpeaks(-theta_L(up_time(j,1:up_time_end)),'minpeakdistance',15);
%             
%             plot(phi_L(down_time(j,1:down_time_end)),theta_L(down_time(j,1:down_time_end)),'Color',[0.5 0.5 0.5])
%             
%             plot(phi_L(up_time(j,1:up_time_end)),theta_L(up_time(j,1:up_time_end)),'Color',[0.5 0.5 0.5])
%             
%             plot(phi_L(down_time(j,1)-1+pks_pos_id_down),theta_L(down_time(j,1)-1+pks_pos_id_down),'o','Color','r') 
%             plot(phi_L(down_time(j,1)-1+pks_neg_id_down),theta_L(down_time(j,1)-1+pks_neg_id_down),'o','Color','b')
%             plot(phi_L(up_time(j,1)-1+pks_pos_id_up),theta_L(up_time(j,1)-1+pks_pos_id_up),'o','Color','r') 
%             plot(phi_L(up_time(j,1)-1+pks_neg_id_up),theta_L(up_time(j,1)-1+pks_neg_id_up),'o','Color','b')
%         end
%         hold off      
%         




%         figure()
%         hold on
%         for j = 1:nr_wb
%        
%             down_time_end = find(isnan(down_time(j,:))==0, 1, 'last' );
%             
%             up_time_end =  find(isnan(up_time(j,:))==0, 1, 'last' );
%             
%             plot(phi_L(down_time(j,1:down_time_end)),theta_L(down_time(j,1:down_time_end)),'Color',[0.5 0.5 0.5])
%             
%             plot(phi_L(up_time(j,1:up_time_end)),theta_L(up_time(j,1:up_time_end)),'Color',[0.5 0.5 0.5])
%             
%             [ ~, theta_max_id_down] = findpeaks(theta_L(down_time(j,1:down_time_end)));
%             
%             [ ~, theta_max_id_up] = findpeaks(theta_L(up_time(j,1:up_time_end)));
%             
%             plot(phi_L(down_time(j,theta_max_id_down)),theta_L(down_time(j,theta_max_id_down)),'o','Color','r')
%             
%             plot(phi_L(up_time(j,theta_max_id_up)),theta_L(up_time(j,theta_max_id_up)),'o','Color','b')
%             
%         end
%         hold off        
        

        pause

    end

end
