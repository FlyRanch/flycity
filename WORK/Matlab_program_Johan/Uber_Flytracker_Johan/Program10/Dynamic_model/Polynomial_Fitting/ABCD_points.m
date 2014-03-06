function [ ABCD ] = ABCD_points( a_fit, a_avg, down_up_ratio, down_up_avg, n_pol_theta, n_pol_eta, n_pol_phi, nr_wb )
    

    

    % Define four characteristic points of the wingstroke according to the
    % strokeplane deviation local maxima and minima:

    
    
    ABCD = {};
    
    
    
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
            
            [~, A_t_id] = max(theta_L_wb(1:15));
            
            [~, B_t_id] = min(theta_L_wb(10:40));
            
            [~, C_t_id] = max(theta_L_wb(30:70));
            
            [~, D_t_id] = min(theta_L_wb(60:90));
            
            theta_L_points(k,:) = [ A_t_id 9+B_t_id 29+C_t_id 59+D_t_id NaN];
            
            end
        
        end
        
        
    theta_R_points = nan(nr_wb,5);
        
        for k = 1:nr_wb
            
            [ t_wb, X_wb ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_ratio(k), 100, 0, 1, 0 );
            
            theta_R_wb = X_wb*[a_fit.theta_R1(:,k); a_fit.theta_R2(:,k)];
            
            if k > 1
                
                [ t_wb_1, X_wb_1 ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_ratio(k-1), 100, 0, 1, 0 );
            
                theta_R_wb_1 = X_wb_1*[a_fit.theta_R1(:,k-1); a_fit.theta_R2(:,k-1)];
                
                [A_max_L, A_max_L_id] = max(theta_R_wb_1(80:100));
                
                [A_max_R, A_max_R_id] = max(theta_R_wb(1:15));
                
                if A_max_L > A_max_R
                    
                    theta_R_points(k-1,5) = 79+A_max_L_id;
                    
                    theta_R_points(k,1) = NaN;
                    
                elseif A_max_L < A_max_R
                    
                    theta_R_points(k,1) = A_max_R_id;
                    
                    theta_R_points(k-1,5) = NaN;
                    
                elseif A_max_L == A_max_R
                    
                    theta_R_points(k,1) = A_max_R_id;
                    
                    theta_R_points(k-1,5) = 79+A_max_L_id;
                    
                elseif (79+A_max_L_id) == 100
                    
                    theta_R_points(k,1) = 1;
                    
                    theta_R_points(k-1,5) = 100;
                    
                end
                
                [~, B_t_id] = min(theta_R_wb(10:40));
            
                [~, C_t_id] = max(theta_R_wb(30:70));

                [~, D_t_id] = min(theta_R_wb(60:90));

                theta_R_points(k,2:4) = [ 9+B_t_id 29+C_t_id 59+D_t_id ];
                
            else
            
            [~, A_t_id] = max(theta_R_wb(1:15));
            
            [~, B_t_id] = min(theta_R_wb(10:40));
            
            [~, C_t_id] = max(theta_R_wb(30:70));
            
            [~, D_t_id] = min(theta_R_wb(60:90));
            
            theta_R_points(k,:) = [ A_t_id 9+B_t_id 29+C_t_id 59+D_t_id NaN];
            
            end
        
        end
        
        
        [ ~, X_avg_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_avg, 100, 0, 1, 0 );
        
        [ ~, X_avg_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_avg, 100, 0, 1, 0 );
        
        [ ~, X_avg_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_avg, 100, 0, 1, 0 );
        
        theta_L_avg = X_avg_theta*[a_avg.theta_L1; a_avg.theta_L2];
        
        theta_R_avg = X_avg_theta*[a_avg.theta_R1; a_avg.theta_R2];
        
        eta_L_avg = X_avg_eta*[a_avg.eta_L1; a_avg.eta_L2];
        
        eta_R_avg = X_avg_eta*[a_avg.eta_R1; a_avg.eta_R2];

        phi_L_avg = X_avg_phi*[a_avg.phi_L1; a_avg.phi_L2];
        
        phi_R_avg = X_avg_phi*[a_avg.phi_R1; a_avg.phi_R2];
        
        
        
            
        [~, A_L_id] = max(theta_L_avg(1:15));
           
        [~, B_L_id] = min(theta_L_avg(10:40));
            
        [~, C_L_id] = max(theta_L_avg(30:70));
            
        [~, D_L_id] = min(theta_L_avg(60:90));
            
        avg_L_points = [ A_L_id 9+B_L_id 29+C_L_id 59+D_L_id NaN];
            
        [~, A_R_id] = max(theta_R_avg(1:15));
           
        [~, B_R_id] = min(theta_R_avg(10:40));
            
        [~, C_R_id] = max(theta_R_avg(30:70));
           
        [~, D_R_id] = min(theta_R_avg(60:90));
           
        avg_R_points = [ A_R_id 9+B_R_id 29+C_R_id 59+D_R_id NaN];

        
        dev_L_theta = zeros(nr_wb,4);
        
        dev_L_eta = zeros(nr_wb,4);
        
        dev_L_phi = zeros(nr_wb,4);
        
        dev_R_theta = zeros(nr_wb,4);
        
        dev_R_eta = zeros(nr_wb,4);
        
        dev_R_phi = zeros(nr_wb,4);
        
        
        
        
        
        for k = 1:nr_wb
            
            [ ~, X_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_ratio(k), 100, 0, 1, 0 );
            
            [ ~, X_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_ratio(k), 100, 0, 1, 0 );
            
            [ ~, X_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_ratio(k), 100, 0, 1, 0);
             
            theta_L_wb = X_theta*[a_fit.theta_L1(:,k); a_fit.theta_L2(:,k)];
            
            eta_L_wb = X_eta*[a_fit.eta_L1(:,k); a_fit.eta_L2(:,k)];
            
            phi_L_wb = X_phi*[a_fit.phi_L1(:,k); a_fit.phi_L2(:,k)];
            
            if k > 1
            
                if isnan(theta_L_points(k,1))==0
                    dev_L_theta(k,1) = theta_L_wb(theta_L_points(k,1))-theta_L_avg(avg_L_points(1));
                    dev_L_eta(k,1) = eta_L_wb(theta_L_points(k,1))-eta_L_avg(avg_L_points(1));
                    dev_L_phi(k,1) = phi_L_wb(theta_L_points(k,1))-phi_L_avg(avg_L_points(1));
                end
                
                if isnan(theta_L_points(k-1,5))==0
                    dev_L_theta(k,1) = theta_L_wb(theta_L_points(k-1,5))-theta_L_avg(avg_L_points(1));
                    dev_L_eta(k,1) = eta_L_wb(theta_L_points(k-1,5))-eta_L_avg(avg_L_points(1));
                    dev_L_phi(k,1) = phi_L_wb(theta_L_points(k-1,5))-phi_L_avg(avg_L_points(1));
                end
                
                if isnan(theta_L_points(k,2))==0
                    dev_L_theta(k,2) = theta_L_wb(theta_L_points(k,2))-theta_L_avg(avg_L_points(2));
                    dev_L_eta(k,2) = eta_L_wb(theta_L_points(k,2))-eta_L_avg(avg_L_points(2));
                    dev_L_phi(k,2) = phi_L_wb(theta_L_points(k,2))-phi_L_avg(avg_L_points(2));
                end
                
                if isnan(theta_L_points(k,3))==0
                    dev_L_theta(k,3) = theta_L_wb(theta_L_points(k,3))-theta_L_avg(avg_L_points(3));
                    dev_L_eta(k,3) = eta_L_wb(theta_L_points(k,3))-eta_L_avg(avg_L_points(3));
                    dev_L_phi(k,3) = phi_L_wb(theta_L_points(k,3))-phi_L_avg(avg_L_points(3));
                end
                
                if isnan(theta_L_points(k,4))==0
                    dev_L_theta(k,4) = theta_L_wb(theta_L_points(k,4))-theta_L_avg(avg_L_points(4));
                    dev_L_eta(k,4) = eta_L_wb(theta_L_points(k,4))-eta_L_avg(avg_L_points(4));
                    dev_L_phi(k,4) = phi_L_wb(theta_L_points(k,4))-phi_L_avg(avg_L_points(4));
                end
                
            else
                
                if isnan(theta_L_points(k,1))==0
                    dev_L_theta(k,1) = theta_L_wb(theta_L_points(k,1))-theta_L_avg(avg_L_points(1));
                    dev_L_eta(k,1) = eta_L_wb(theta_L_points(k,1))-eta_L_avg(avg_L_points(1));
                    dev_L_phi(k,1) = phi_L_wb(theta_L_points(k,1))-phi_L_avg(avg_L_points(1));
                end
                
                if isnan(theta_L_points(k,2))==0
                    dev_L_theta(k,2) = theta_L_wb(theta_L_points(k,2))-theta_L_avg(avg_L_points(2));
                    dev_L_eta(k,2) = eta_L_wb(theta_L_points(k,2))-eta_L_avg(avg_L_points(2));
                    dev_L_phi(k,2) = phi_L_wb(theta_L_points(k,2))-phi_L_avg(avg_L_points(2));
                end
                
                if isnan(theta_L_points(k,3))==0
                    dev_L_theta(k,3) = theta_L_wb(theta_L_points(k,3))-theta_L_avg(avg_L_points(3));
                    dev_L_eta(k,3) = eta_L_wb(theta_L_points(k,3))-eta_L_avg(avg_L_points(3));
                    dev_L_phi(k,3) = phi_L_wb(theta_L_points(k,3))-phi_L_avg(avg_L_points(3));
                end
                
                if isnan(theta_L_points(k,4))==0
                    dev_L_theta(k,4) = theta_L_wb(theta_L_points(k,4))-theta_L_avg(avg_L_points(4));
                    dev_L_eta(k,4) = eta_L_wb(theta_L_points(k,4))-eta_L_avg(avg_L_points(4));
                    dev_L_phi(k,4) = phi_L_wb(theta_L_points(k,4))-phi_L_avg(avg_L_points(4));
                end
                
            end
        end
        
        for k = 1:nr_wb
            
            [ ~, X_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_ratio(k), 100, 0, 1, 0 );
            
            [ ~, X_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_ratio(k), 100, 0, 1, 0 );
            
            theta_R_wb = X_theta*[a_fit.theta_R1(:,k); a_fit.theta_R2(:,k)];
            
            eta_R_wb = X_eta*[a_fit.eta_R1(:,k); a_fit.eta_R2(:,k)];
            
            phi_R_wb = X_phi*[a_fit.phi_R1(:,k); a_fit.phi_R2(:,k)];
                        
            if k > 1
            
                if isnan(theta_R_points(k,1))==0
                    dev_R_theta(k,1) = theta_R_wb(theta_R_points(k,1))-theta_R_avg(avg_R_points(1));
                    dev_R_eta(k,1) = eta_R_wb(theta_R_points(k,1))-eta_R_avg(avg_R_points(1));
                    dev_R_phi(k,1) = phi_R_wb(theta_R_points(k,1))-phi_R_avg(avg_R_points(1));
                end
                
                if isnan(theta_R_points(k-1,5))==0
                    dev_R_theta(k,1) = theta_R_wb(theta_R_points(k-1,5))-theta_R_avg(avg_R_points(1));
                    dev_R_eta(k,1) = eta_R_wb(theta_R_points(k-1,5))-eta_R_avg(avg_R_points(1));
                    dev_R_phi(k,1) = phi_R_wb(theta_R_points(k-1,5))-phi_R_avg(avg_R_points(1));
                end
                
                if isnan(theta_R_points(k,2))==0
                    dev_R_theta(k,2) = theta_R_wb(theta_R_points(k,2))-theta_R_avg(avg_R_points(2));
                    dev_R_eta(k,2) = eta_R_wb(theta_R_points(k,2))-eta_R_avg(avg_R_points(2));
                    dev_R_phi(k,2) = phi_R_wb(theta_R_points(k,2))-phi_R_avg(avg_R_points(2));
                end
                
                if isnan(theta_R_points(k,3))==0
                    dev_R_theta(k,3) = theta_R_wb(theta_R_points(k,3))-theta_R_avg(avg_R_points(3));
                    dev_R_eta(k,3) = eta_R_wb(theta_R_points(k,3))-eta_R_avg(avg_R_points(3));
                    dev_R_phi(k,3) = phi_R_wb(theta_R_points(k,3))-phi_R_avg(avg_R_points(3));
                end
                
                if isnan(theta_R_points(k,4))==0
                    dev_R_theta(k,4) = theta_R_wb(theta_R_points(k,4))-theta_R_avg(avg_R_points(4));
                    dev_R_eta(k,4) = eta_R_wb(theta_R_points(k,4))-eta_R_avg(avg_R_points(4));
                    dev_R_phi(k,4) = phi_R_wb(theta_R_points(k,4))-phi_R_avg(avg_R_points(4));
                end
                
            else
                
                if isnan(theta_R_points(k,1))==0
                    dev_R_theta(k,1) = theta_R_wb(theta_R_points(k,1))-theta_R_avg(avg_R_points(1));
                    dev_R_eta(k,1) = eta_R_wb(theta_R_points(k,1))-eta_R_avg(avg_R_points(1));
                    dev_R_phi(k,1) = phi_R_wb(theta_R_points(k,1))-phi_R_avg(avg_R_points(1));
                end
                
                if isnan(theta_R_points(k,2))==0
                    dev_R_theta(k,2) = theta_R_wb(theta_R_points(k,2))-theta_R_avg(avg_R_points(2));
                    dev_R_eta(k,2) = eta_R_wb(theta_R_points(k,2))-eta_R_avg(avg_R_points(2));
                    dev_R_phi(k,2) = phi_R_wb(theta_R_points(k,2))-phi_R_avg(avg_R_points(2));
                end
                
                if isnan(theta_R_points(k,3))==0
                    dev_R_theta(k,3) = theta_R_wb(theta_R_points(k,3))-theta_R_avg(avg_R_points(3));
                    dev_R_eta(k,3) = eta_R_wb(theta_R_points(k,3))-eta_R_avg(avg_R_points(3));
                    dev_R_phi(k,3) = phi_R_wb(theta_R_points(k,3))-phi_R_avg(avg_R_points(3));
                end
                
                if isnan(theta_R_points(k,4))==0
                    dev_R_theta(k,4) = theta_R_wb(theta_R_points(k,4))-theta_R_avg(avg_R_points(4));
                    dev_R_eta(k,4) = eta_R_wb(theta_R_points(k,4))-eta_R_avg(avg_R_points(4));
                    dev_R_phi(k,4) = phi_R_wb(theta_R_points(k,4))-phi_R_avg(avg_R_points(4));
                end
                
            end
        end
                
        
%         figure()
%         hold on
%         for j = 1:nr_wb
%             
%             [ t_wb, X_wb ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_ratio(j), 100, 0, 1, 0 );
%             
%             theta_L_wb = X_wb*[a_fit.theta_L1(:,j); a_fit.theta_L2(:,j)];
%             
%             plot(t_wb,theta_L_wb,'Color',[0.5 0.5 0.5])
%             
%             if isnan(theta_L_points(j,1))==0
%                 plot(t_wb(theta_L_points(j,1)),theta_L_wb(theta_L_points(j,1)),'o','Color','r')
%             end
%             
%             if isnan(theta_L_points(j,2))==0
%                 plot(t_wb(theta_L_points(j,2)),theta_L_wb(theta_L_points(j,2)),'o','Color','g')
%             end
%             
%             if isnan(theta_L_points(j,3))==0
%                 plot(t_wb(theta_L_points(j,3)),theta_L_wb(theta_L_points(j,3)),'o','Color','b')
%             end
%             
%             if isnan(theta_L_points(j,4))==0
%                 plot(t_wb(theta_L_points(j,4)),theta_L_wb(theta_L_points(j,4)),'o','Color','m')
%             end
%             
%             if isnan(theta_L_points(j,5))==0
%                 plot(t_wb(theta_L_points(j,5)),theta_L_wb(theta_L_points(j,5)),'o','Color','r')
%             end
%             
%         end
%         hold off
%         
%         
%         figure()
%         hold on
%         for j = 1:nr_wb
%             
%             [ t_wb, X_wb_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_ratio(j), 100, 0, 1, 0 );
%             
%             [ ~, X_wb_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_ratio(j), 100, 0, 1, 0 );
%             
%             theta_L_wb = X_wb_theta*[a_fit.theta_L1(:,j); a_fit.theta_L2(:,j)];
%             
%             phi_L_wb = X_wb_phi*[a_fit.phi_L1(:,j); a_fit.phi_L2(:,j)];
%             
%             plot(phi_L_wb,theta_L_wb,'Color',[0.5 0.5 0.5])
%             
%             if isnan(theta_L_points(j,1))==0
%                 plot(phi_L_wb(theta_L_points(j,1)),theta_L_wb(theta_L_points(j,1)),'o','Color','r')
%             end
%             
%             if isnan(theta_L_points(j,2))==0
%                 plot(phi_L_wb(theta_L_points(j,2)),theta_L_wb(theta_L_points(j,2)),'o','Color','g')
%             end
%             
%             if isnan(theta_L_points(j,3))==0
%                 plot(phi_L_wb(theta_L_points(j,3)),theta_L_wb(theta_L_points(j,3)),'o','Color','b')
%             end
%             
%             if isnan(theta_L_points(j,4))==0
%                 plot(phi_L_wb(theta_L_points(j,4)),theta_L_wb(theta_L_points(j,4)),'o','Color','m')
%             end
%             
%             if isnan(theta_L_points(j,5))==0
%                 plot(phi_L_wb(theta_L_points(j,5)),theta_L_wb(theta_L_points(j,5)),'o','Color','r')
%             end
%             
%         end
%         hold off
% 
%         
% 
%         
%         
%         figure()
%         hold on
%         for j = 1:nr_wb
%             
%             [ t_wb, X_wb ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_ratio(j), 100, 0, 1, 0 );
%             
%             theta_R_wb = X_wb*[a_fit.theta_R1(:,j); a_fit.theta_R2(:,j)];
%             
%             plot(t_wb,theta_R_wb,'Color',[0.5 0.5 0.5])
%             
%             if isnan(theta_R_points(j,1))==0
%                 plot(t_wb(theta_R_points(j,1)),theta_R_wb(theta_R_points(j,1)),'o','Color','r')
%             end
%             
%             if isnan(theta_R_points(j,2))==0
%                 plot(t_wb(theta_R_points(j,2)),theta_R_wb(theta_R_points(j,2)),'o','Color','g')
%             end
%             
%             if isnan(theta_R_points(j,3))==0
%                 plot(t_wb(theta_R_points(j,3)),theta_R_wb(theta_R_points(j,3)),'o','Color','b')
%             end
%             
%             if isnan(theta_R_points(j,4))==0
%                 plot(t_wb(theta_R_points(j,4)),theta_R_wb(theta_R_points(j,4)),'o','Color','m')
%             end
%             
%             if isnan(theta_R_points(j,5))==0
%                 plot(t_wb(theta_R_points(j,5)),theta_R_wb(theta_R_points(j,5)),'o','Color','r')
%             end
%             
%         end
%         hold off
%         
%         
%         figure()
%         hold on
%         for j = 1:nr_wb
%             
%             [ t_wb, X_wb_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_ratio(j), 100, 0, 1, 0 );
%             
%             [ ~, X_wb_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_ratio(j), 100, 0, 1, 0 );
%             
%             theta_R_wb = X_wb_theta*[a_fit.theta_R1(:,j); a_fit.theta_R2(:,j)];
%             
%             phi_R_wb = X_wb_phi*[a_fit.phi_R1(:,j); a_fit.phi_R2(:,j)];
%             
%             plot(phi_R_wb,theta_R_wb,'Color',[0.5 0.5 0.5])
%             
%             if isnan(theta_R_points(j,1))==0
%                 plot(phi_R_wb(theta_R_points(j,1)),theta_R_wb(theta_R_points(j,1)),'o','Color','r')
%             end
%             
%             if isnan(theta_R_points(j,2))==0
%                 plot(phi_R_wb(theta_R_points(j,2)),theta_R_wb(theta_R_points(j,2)),'o','Color','g')
%             end
%             
%             if isnan(theta_R_points(j,3))==0
%                 plot(phi_R_wb(theta_R_points(j,3)),theta_R_wb(theta_R_points(j,3)),'o','Color','b')
%             end
%             
%             if isnan(theta_R_points(j,4))==0
%                 plot(phi_R_wb(theta_R_points(j,4)),theta_R_wb(theta_R_points(j,4)),'o','Color','m')
%             end
%             
%             if isnan(theta_R_points(j,5))==0
%                 plot(phi_R_wb(theta_R_points(j,5)),theta_R_wb(theta_R_points(j,5)),'o','Color','r')
%             end
%             
%         end
%         hold off
%         
%         
%         figure()
%         hold on
%         for j = 1:nr_wb
%             plot(phi_L_avg(1)+dev_L_phi(j,1),theta_L_avg(1)+dev_L_theta(j,1),'o','Color','r')
%             plot(phi_L_avg(2)+dev_L_phi(j,2),theta_L_avg(2)+dev_L_theta(j,2),'o','Color','g')
%             plot(phi_L_avg(3)+dev_L_phi(j,3),theta_L_avg(3)+dev_L_theta(j,3),'o','Color','b')
%             plot(phi_L_avg(4)+dev_L_phi(j,4),theta_L_avg(4)+dev_L_theta(j,4),'o','Color','m')
%         end
%         [ ~, X_avg_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_avg, 100, 0, 1, 0 ); 
%         [ ~, X_avg_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_avg, 100, 0, 1, 0 );
%         theta_avg =  X_avg_theta*[a_avg.theta_L1; a_avg.theta_L2];
%         phi_avg = X_avg_phi*[a_avg.phi_L1; a_avg.phi_L2];
%         plot(phi_avg,theta_avg,'Color',[0.5 0.5 0.5])
%         hold off
%         
% 
%         figure()
%         hold on
%         for j = 1:nr_wb
%             plot(phi_R_avg(1)+dev_R_phi(j,1),theta_R_avg(1)+dev_R_theta(j,1),'o','Color','r')
%             plot(phi_R_avg(2)+dev_R_phi(j,2),theta_R_avg(2)+dev_R_theta(j,2),'o','Color','g')
%             plot(phi_R_avg(3)+dev_R_phi(j,3),theta_R_avg(3)+dev_R_theta(j,3),'o','Color','b')
%             plot(phi_R_avg(4)+dev_R_phi(j,4),theta_R_avg(4)+dev_R_theta(j,4),'o','Color','m')
%         end
%         [ ~, X_avg_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_avg, 100, 0, 1, 0 ); 
%         [ ~, X_avg_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_avg, 100, 0, 1, 0 );
%         theta_avg =  X_avg_theta*[a_avg.theta_R1; a_avg.theta_R2];
%         phi_avg = X_avg_phi*[a_avg.phi_R1; a_avg.phi_R2];
%         plot(phi_avg,theta_avg,'Color',[0.5 0.5 0.5])
%         hold off

        
        ABCD.theta_L_points = theta_L_points;
        ABCD.theta_R_points = theta_R_points;
        ABCD.avg_L_points = avg_L_points;
        ABCD.avg_R_points = avg_R_points;
        ABCD.avg_L_theta = theta_L_avg;
        ABCD.avg_R_theta = theta_R_avg;
        ABCD.avg_L_eta = eta_L_avg;
        ABCD.avg_R_eta = eta_R_avg;
        ABCD.avg_L_phi = phi_L_avg;
        ABCD.avg_R_phi = phi_R_avg;
        ABCD.dev_L_theta = dev_L_theta;
        ABCD.dev_R_theta = dev_R_theta;
        ABCD.dev_L_eta = dev_L_eta;
        ABCD.dev_R_eta = dev_R_eta;
        ABCD.dev_L_phi = dev_L_phi;
        ABCD.dev_R_phi = dev_R_phi;
        




end

