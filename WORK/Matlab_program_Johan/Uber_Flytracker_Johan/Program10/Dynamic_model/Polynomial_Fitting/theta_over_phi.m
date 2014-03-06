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

%         figure()
%         hold on
%         for j = 1:nr_wb
%             
%             down_time_end = find(isnan(down_time(j,:))==0, 1, 'last' );
%             
%             plot(phi_L(down_time(j,1:down_time_end)),theta_L(down_time(j,1:down_time_end)),'Color',[0.5 0.5 0.5])
%             
%         end
%         hold off
%         
%         figure()
%         hold on
%         for j = 1:nr_wb
%        
%             up_time_end =  find(isnan(up_time(j,:))==0, 1, 'last' );
%             
%             plot(phi_L(up_time(j,1:up_time_end)),theta_L(up_time(j,1:up_time_end)),'Color',[0.5 0.5 0.5])
%             
%         end
%         hold off

        
        
        [a_fit,a_avg,f_avg,down_up,trigger_wb,down_up_ratio] = Standard_wingbeat( settings, pathDB, i, n_pol_theta, n_pol_eta, n_pol_phi );
        
%         figure()
%         hold on
%         for j = 1:nr_wb
%             
%             up_time_end =  find(isnan(up_time(j,:))==0, 1, 'last' );
%             
%             [ t_wb, X_wb ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_ratio(j), 100, t(down_time(j,1)), t(up_time(j,up_time_end)+1), 0 );
%             
%             plot(t_wb,X_wb*[a_fit.theta_L1(:,j); a_fit.theta_L2(:,j)])
%             
%         end
%         hold off
%         
%         figure()
%         hold on
%         for j = 1:nr_wb
%             
%             up_time_end =  find(isnan(up_time(j,:))==0, 1, 'last' );
%             
%             [ t_wb, X_wb_dot ] = Wingbeat_Legendre_matrix( n_pol_theta,down_up_ratio(j), 100, t(down_time(j,1)), t(up_time(j,up_time_end)+1), 1 );
%             
%             plot(t_wb,X_wb_dot*[a_fit.theta_L1(:,j); a_fit.theta_L2(:,j)])
%             
%         end
%         hold off
%         
%         
%         [theta_pks_pos, theta_pks_pos_id] = findpeaks(theta_L,'minpeakdistance',1);
%         
%         [theta_pks_neg, theta_pks_neg_id] = findpeaks(-theta_L,'minpeakdistance',1);        
%         
%         [eta_pks_pos, eta_pks_pos_id] = findpeaks(eta_L,'minpeakdistance',26);
%         
%         [eta_pks_neg, eta_pks_neg_id] = findpeaks(-eta_L,'minpeakdistance',26);
%         
%         
%         
%         % Find intersections of eta with eta = 90 deg.
%         
%         eta_90_id = 0;
%         
%         for j = 2:(stop-start+1)
%             
%             if eta_L(j)>= pi/2 && eta_L(j-1) <= pi/2
%                 
%                 if eta_90_id == 0
%                     
%                     if abs((pi/2)-eta_L(j)) >= abs((pi/2)-eta_L(j-1))
%                         
%                         eta_90_id = j;
%                         
%                     elseif abs((pi/2)-eta_L(j)) < abs((pi/2)-eta_L(j-1))
%                         
%                         eta_90_id = j-1;
%                         
%                     end
%                     
%                 else
%                     
%                     if abs((pi/2)-eta_L(j)) >= abs((pi/2)-eta_L(j-1))
%                         
%                         eta_90_id = [eta_90_id; j];
%                         
%                     elseif abs((pi/2)-eta_L(j)) < abs((pi/2)-eta_L(j-1))
%                         
%                         eta_90_id = [eta_90_id; j-1];
%                         
%                     end
%                     
%                 end
%             
%             end
%             
%             if eta_L(j) < pi/2 && eta_L(j-1) >= pi/2
%                 
%                 if eta_90_id == 0
%                     
%                     if abs((pi/2)-eta_L(j)) >= abs((pi/2)-eta_L(j-1))
%                         
%                         eta_90_id = j;
%                         
%                     elseif abs((pi/2)-eta_L(j)) < abs((pi/2)-eta_L(j-1))
%                         
%                         eta_90_id = j-1;
%                         
%                     end
%                     
%                 else
%                     
%                     if abs((pi/2)-eta_L(j)) >= abs((pi/2)-eta_L(j-1))
%                         
%                         eta_90_id = [eta_90_id; j];
%                         
%                     elseif abs((pi/2)-eta_L(j)) < abs((pi/2)-eta_L(j-1))
%                         
%                         eta_90_id = [eta_90_id; j-1];
%                         
%                     end
%                     
%                 end
%             
%             end
%             
%         end

%         eta_90_id
%         
%         figure()
%         plot(t,eta_L)
%         hold on
%         plot(t(eta_pks_pos_id), eta_L(eta_pks_pos_id),'o','Color','r')
%         plot(t(eta_pks_neg_id), eta_L(eta_pks_neg_id),'o','Color','r')
%         plot(t(eta_90_id),eta_L(eta_90_id),'o','Color','g')
%         hold off
%         
%         figure()
%         plot(t,theta_L)
%         hold on
%         plot(t(theta_pks_pos_id), theta_L(theta_pks_pos_id),'o','Color','r')
%         plot(t(theta_pks_neg_id), theta_L(theta_pks_neg_id),'o','Color','g')
%         hold off
%         
%         figure()
%         plot(t,theta_L)
%         hold on
%         plot(t(eta_pks_pos_id), theta_L(eta_pks_pos_id),'o','Color','r')
%         plot(t(eta_pks_neg_id), theta_L(eta_pks_neg_id),'o','Color','g')
%         hold off
        
        theta_L_points = nan(nr_wb,5);
        
        for k = 1:nr_wb
            
            [ t_wb, X_wb ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_ratio(k), 100, 0, 1, 0 );
            
            theta_L_wb = X_wb*[a_fit.theta_L1(:,k); a_fit.theta_L2(:,k)];
            
            if k > 1
                
                [ t_wb_1, X_wb_1 ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_ratio(k-1), 100, 0, 1, 0 );
            
                theta_L_wb_1 = X_wb_1*[a_fit.theta_L1(:,k-1); a_fit.theta_L2(:,k-1)];
                
                [A_max_L, A_max_L_id] = max(theta_L_wb_1(80:100));
                
                [A_max_R, A_max_R_id] = max(theta_L_wb(1:15));
                
                if A_max_L > A_max_R
                    
                    theta_L_points(k-1,5) = 79+A_max_L_id;
                    
                    theta_L_points(k,1) = NaN;
                    
                elseif A_max_L < A_max_R
                    
                    theta_L_points(k,1) = A_max_R_id;
                    
                    theta_L_points(k-1,5) = NaN;
                    
                elseif A_max_L == A_max_R
                    
                    theta_L_points(k,1) = A_max_R_id;
                    
                    theta_L_points(k-1,5) = 79+A_max_L_id;
                    
                elseif (79+A_max_L_id) == 100
                    
                    theta_L_points(k,1) = 1;
                    
                    theta_L_points(k-1,5) = 100;
                    
                end
                
                [~, B_t_id] = min(theta_L_wb(10:40));
            
                [~, C_t_id] = max(theta_L_wb(30:70));

                [~, D_t_id] = min(theta_L_wb(60:90));

                theta_L_points(k,2:4) = [ 9+B_t_id 29+C_t_id 59+D_t_id ];
                
            else
            
            [~, A_t_id] = min(theta_L_wb(1:15));
            
            [~, B_t_id] = min(theta_L_wb(10:40));
            
            [~, C_t_id] = max(theta_L_wb(30:70));
            
            [~, D_t_id] = min(theta_L_wb(60:90));
            
            theta_L_points(k,:) = [ A_t_id 9+B_t_id 29+C_t_id 59+D_t_id NaN];
            
            end
        
        end
        
        theta_L_points
       

        
        figure()
        hold on
        for j = 1:nr_wb
            
            [ t_wb, X_wb ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_ratio(j), 100, 0, 1, 0 );
            
            theta_L_wb = X_wb*[a_fit.theta_L1(:,j); a_fit.theta_L2(:,j)];
            
            plot(t_wb,theta_L_wb,'Color',[0.5 0.5 0.5])
            
            if isnan(theta_L_points(j,1))==0
                plot(t_wb(theta_L_points(j,1)),theta_L_wb(theta_L_points(j,1)),'o','Color','r')
            end
            
            if isnan(theta_L_points(j,2))==0
                plot(t_wb(theta_L_points(j,2)),theta_L_wb(theta_L_points(j,2)),'o','Color','g')
            end
            
            if isnan(theta_L_points(j,3))==0
                plot(t_wb(theta_L_points(j,3)),theta_L_wb(theta_L_points(j,3)),'o','Color','b')
            end
            
            if isnan(theta_L_points(j,4))==0
                plot(t_wb(theta_L_points(j,4)),theta_L_wb(theta_L_points(j,4)),'o','Color','m')
            end
            
            if isnan(theta_L_points(j,5))==0
                plot(t_wb(theta_L_points(j,5)),theta_L_wb(theta_L_points(j,5)),'o','Color','r')
            end
            
        end
        hold off
        
        
        figure()
        hold on
        for j = 1:nr_wb
            
            [ t_wb, X_wb_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_ratio(j), 100, 0, 1, 0 );
            
            [ ~, X_wb_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_ratio(j), 100, 0, 1, 0 );
            
            theta_L_wb = X_wb_theta*[a_fit.theta_L1(:,j); a_fit.theta_L2(:,j)];
            
            phi_L_wb = X_wb_phi*[a_fit.phi_L1(:,j); a_fit.phi_L2(:,j)];
            
            plot(phi_L_wb,theta_L_wb,'Color',[0.5 0.5 0.5])
            
            if isnan(theta_L_points(j,1))==0
                plot(phi_L_wb(theta_L_points(j,1)),theta_L_wb(theta_L_points(j,1)),'o','Color','r')
            end
            
            if isnan(theta_L_points(j,2))==0
                plot(phi_L_wb(theta_L_points(j,2)),theta_L_wb(theta_L_points(j,2)),'o','Color','g')
            end
            
            if isnan(theta_L_points(j,3))==0
                plot(phi_L_wb(theta_L_points(j,3)),theta_L_wb(theta_L_points(j,3)),'o','Color','b')
            end
            
            if isnan(theta_L_points(j,4))==0
                plot(phi_L_wb(theta_L_points(j,4)),theta_L_wb(theta_L_points(j,4)),'o','Color','m')
            end
            
            if isnan(theta_L_points(j,5))==0
                plot(phi_L_wb(theta_L_points(j,5)),theta_L_wb(theta_L_points(j,5)),'o','Color','r')
            end
            
        end
        hold off

        
        figure()
        hold on
        for j = 1:nr_wb
            
            [ t_wb, X_wb_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_ratio(j), 100, 0, 1, 0 );
            
            [ ~, X_wb_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_ratio(j), 100, 0, 1, 0 );
            
            theta_L_wb = X_wb_theta*[a_fit.theta_L1(:,j); a_fit.theta_L2(:,j)];
            
            phi_L_wb = X_wb_phi*[a_fit.phi_L1(:,j); a_fit.phi_L2(:,j)];
            
            plot(phi_L_wb,theta_L_wb,'r')
            
        end
        plot(phi_L,theta_L,'g')
        hold off
        
        
        pause
        
%         figure()
%         plot(t,-theta_L.*phi_L)
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
%         end
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


    end

end

