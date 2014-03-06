function [a_sym, a_dev] = Maneuver_wingbeat(a_fit,a_avg,n_pol_theta,n_pol_eta,n_pol_phi,trigger_wb,ratio_1_2,ratio_1_2_avg)


    % Determine the deviation between the average left wingbeat and the
    % actual left wingbeats and the deviation between the average right
    % wingbeat and the actual right wingbeats:
    
    nr_wb = length(a_fit.theta_L1(1,:));
    
    a_dev = {};
        
    a_dev_theta_L1 = zeros(n_pol_theta+1,nr_wb);
    
    a_dev_theta_L2 = zeros(n_pol_theta+1,nr_wb);
    
    a_dev_eta_L1 = zeros(n_pol_eta+1,nr_wb);
    
    a_dev_eta_L2 = zeros(n_pol_eta+1,nr_wb);
    
    a_dev_phi_L1 = zeros(n_pol_phi+1,nr_wb);
    
    a_dev_phi_L2 = zeros(n_pol_phi+1,nr_wb);
    
    
    a_dev_theta_R1 = zeros(n_pol_theta+1,nr_wb);
    
    a_dev_theta_R2 = zeros(n_pol_theta+1,nr_wb);
    
    a_dev_eta_R1 = zeros(n_pol_eta+1,nr_wb);
    
    a_dev_eta_R2 = zeros(n_pol_eta+1,nr_wb);
    
    a_dev_phi_R1 = zeros(n_pol_phi+1,nr_wb);
    
    a_dev_phi_R2 = zeros(n_pol_phi+1,nr_wb);
    
    
    
    for i = 1:nr_wb
        
        a_dev_theta_L1(:,i) = a_fit.theta_L1(:,i)-a_avg.theta_L1;

        a_dev_theta_L2(:,i) = a_fit.theta_L2(:,i)-a_avg.theta_L2;

        a_dev_eta_L1(:,i) = a_fit.eta_L1(:,i)-a_avg.eta_L1;

        a_dev_eta_L2(:,i) = a_fit.eta_L2(:,i)-a_avg.eta_L2;

        a_dev_phi_L1(:,i) = a_fit.phi_L1(:,i)-a_avg.phi_L1;
        
        a_dev_phi_L2(:,i) = a_fit.phi_L2(:,i)-a_avg.phi_L2;


        a_dev_theta_R1(:,i) = a_fit.theta_R1(:,i)-a_avg.theta_R1;

        a_dev_theta_R2(:,i) = a_fit.theta_R2(:,i)-a_avg.theta_R2;

        a_dev_eta_R1(:,i) = a_fit.eta_R1(:,i)-a_avg.eta_R1;

        a_dev_eta_R2(:,i) = a_fit.eta_R2(:,i)-a_avg.eta_R2;

        a_dev_phi_R1(:,i) = a_fit.phi_R1(:,i)-a_avg.phi_R1;

        a_dev_phi_R2(:,i) = a_fit.phi_R2(:,i)-a_avg.phi_R2;   
        
    end
    
    a_dev.theta_L1 = a_dev_theta_L1;
    
    a_dev.theta_L2 = a_dev_theta_L2;
    
    a_dev.eta_L1 = a_dev_eta_L1;
    
    a_dev.eta_L2 = a_dev_eta_L2;
    
    a_dev.phi_L1 = a_dev_phi_L1;
    
    a_dev.phi_L2 = a_dev_phi_L2;
    
     
    a_dev.theta_R1 = a_dev_theta_R1;
    
    a_dev.theta_R2 = a_dev_theta_R2;
    
    a_dev.eta_R1 = a_dev_eta_R1;
    
    a_dev.eta_R2 = a_dev_eta_R2;
    
    a_dev.phi_R1 = a_dev_phi_R1;
    
    a_dev.phi_R2 = a_dev_phi_R2;
    
    
    
    % Add up the deviation to the average symmetric wingbeat to create a
    % 'symmetric' wingbeat:
    
    a_sym = {};
    
    a_sym_theta_L1 = zeros(n_pol_theta+1,nr_wb);
    
    a_sym_theta_L2 = zeros(n_pol_theta+1,nr_wb);
    
    a_sym_eta_L1 = zeros(n_pol_eta+1,nr_wb);
    
    a_sym_eta_L2 = zeros(n_pol_eta+1,nr_wb);
    
    a_sym_phi_L1 = zeros(n_pol_phi+1,nr_wb);
    
    a_sym_phi_L2 = zeros(n_pol_phi+1,nr_wb);
    
    
    a_sym_theta_R1 = zeros(n_pol_theta+1,nr_wb);
    
    a_sym_theta_R2 = zeros(n_pol_theta+1,nr_wb);
    
    a_sym_eta_R1 = zeros(n_pol_eta+1,nr_wb);
    
    a_sym_eta_R2 = zeros(n_pol_eta+1,nr_wb);
    
    a_sym_phi_R1 = zeros(n_pol_phi+1,nr_wb);
    
    a_sym_phi_R2 = zeros(n_pol_phi+1,nr_wb);    
    
    
    for i = 1:nr_wb
        
        a_sym_theta_L1(:,i) = a_avg.theta_LR1+a_dev_theta_L1(:,i);

        a_sym_theta_L2(:,i) = a_avg.theta_LR2+a_dev_theta_L2(:,i);

        a_sym_eta_L1(:,i) = a_avg.eta_LR1+a_dev_eta_L1(:,i);

        a_sym_eta_L2(:,i) = a_avg.eta_LR2+a_dev_eta_L2(:,i);

        a_sym_phi_L1(:,i) = a_avg.phi_LR1+a_dev_phi_L1(:,i);
        
        a_sym_phi_L2(:,i) = a_avg.phi_LR2+a_dev_phi_L2(:,i);


        a_sym_theta_R1(:,i) = a_avg.theta_LR1+a_dev_theta_L1(:,i);

        a_sym_theta_R2(:,i) = a_avg.theta_LR2+a_dev_theta_L2(:,i);

        a_sym_eta_R1(:,i) = a_avg.eta_LR1+a_dev_eta_L1(:,i);

        a_sym_eta_R2(:,i) = a_avg.eta_LR2+a_dev_eta_L2(:,i);

        a_sym_phi_R1(:,i) = a_avg.phi_LR1+a_dev_phi_L1(:,i);

        a_sym_phi_R2(:,i) = a_avg.phi_LR2+a_dev_phi_L2(:,i);   
        
    end    
    
    a_sym.theta_L1 = a_sym_theta_L1;
    
    a_sym.theta_L2 = a_sym_theta_L2;
    
    a_sym.eta_L1 = a_sym_eta_L1;
    
    a_sym.eta_L2 = a_sym_eta_L2;
    
    a_sym.phi_L1 = a_sym_phi_L1;
    
    a_sym.phi_L2 = a_sym_phi_L2;
    
    
    a_sym.theta_R1 = a_sym_theta_R1;
    
    a_sym.theta_R2 = a_sym_theta_R2;
    
    a_sym.eta_R1 = a_sym_eta_R1;
    
    a_sym.eta_R2 = a_sym_eta_R2;
    
    a_sym.phi_R1 = a_sym_phi_R1;
    
    a_sym.phi_R2 = a_sym_phi_R2;
    
    
%     % Determine the local maxima and minima of the average wingbeat:
%     
%     loc_max_min = {};
%     
%     loc_max_min_theta = 0;
%     
%     loc_max_min_eta = 0;
%     
%     loc_max_min_phi = 0;
% 
% %         
    [ ~, X_dot_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, ratio_1_2_avg, 500, 0, 1, 1 );
            
    [ ~, X_dot_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, ratio_1_2_avg, 500, 0, 1, 1 );
            
    [ t_1, X_dot_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, ratio_1_2_avg, 500, 0, 1, 1 );
        
    theta_avg_dot = X_dot_theta*[a_avg.theta_LR1; a_avg.theta_LR2];
        
    eta_avg_dot = X_dot_eta*[a_avg.eta_LR1; a_avg.eta_LR2];
        
    phi_avg_dot = X_dot_phi*[a_avg.phi_LR1; a_avg.phi_LR2];
    
    [ ~, X_ddot_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, ratio_1_2_avg, 500, 0, 1, 2 );
            
    [ ~, X_ddot_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, ratio_1_2_avg, 500, 0, 1, 2 );
            
    [ ~, X_ddot_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, ratio_1_2_avg, 500, 0, 1, 2 );
    
    
%         
%         
%     for j = 2:498
%        
%         if theta_avg_dot(j) > 0 && theta_avg_dot(j+1) < 0
%             
%             loc_max_min_theta = [loc_max_min_theta; t_1(j)];
%             
%         elseif theta_avg_dot(j) < 0 && theta_avg_dot(j+1) > 0
%             
%             loc_max_min_theta = [loc_max_min_theta; t_1(j)];
%             
%         end
%             
%         if eta_avg_dot(j) > 0 && eta_avg_dot(j+1) < 0
%             
%             loc_max_min_eta = [loc_max_min_eta; t_1(j)];
%             
%         elseif eta_avg_dot(j) < 0 && eta_avg_dot(j+1) > 0
%         
%             loc_max_min_eta = [loc_max_min_eta; t_1(j)];
%             
%         end
%         
%         if phi_avg_dot(j) > 0 && phi_avg_dot(j+1) < 0
%             
%             loc_max_min_phi = [loc_max_min_phi; t_1(j)];
%             
%         elseif phi_avg_dot(j) < 0 && phi_avg_dot(j+1) > 0
%         
%             loc_max_min_phi = [loc_max_min_phi; t_1(j)];
%             
%         end
%         
%         
%     end
%     
%     loc_max_min_theta = [loc_max_min_theta; 1];
% 
%     loc_max_min_eta = [loc_max_min_eta; 1];
%     
%     loc_max_min_phi = [loc_max_min_phi; 1];
%     
%     
%     loc_max_min_theta
%     
%     loc_max_min_eta
%     
%     loc_max_min_phi
%     
%     
%     loc_max_min.theta = loc_max_min_theta;
%     
%     loc_max_min.eta = loc_max_min_eta;
%     
%     loc_max_min.phi = loc_max_min_phi;
% 
% 
    X_theta = zeros(500,2*(n_pol_theta+1),nr_wb);
    
    X_eta = zeros(500,2*(n_pol_eta+1),nr_wb);
    
    X_phi = zeros(500,2*(n_pol_phi+1),nr_wb);
    
    for i = 1:nr_wb
    
    [ ~, X_theta_t ] = Wingbeat_Legendre_matrix( n_pol_theta, ratio_1_2(i), 500, 0, 1, 0 );
            
    [ ~, X_eta_t ] = Wingbeat_Legendre_matrix( n_pol_eta, ratio_1_2(i), 500, 0, 1, 0 );
            
    [ t_2, X_phi_t ] = Wingbeat_Legendre_matrix( n_pol_phi, ratio_1_2(i), 500, 0, 1, 0 );
    
    
    X_theta(:,:,i) = X_theta_t;
    
    X_eta(:,:,i) = X_eta_t;
    
    X_phi(:,:,i) = X_phi_t;
    
    end
%     
% %     Mn_theta = Melisoids( loc_max_min.theta, t_2);
% %     
% %     Mn_eta = Melisoids( loc_max_min.eta, t_2);
% %     
% %     Mn_phi = Melisoids( loc_max_min.phi, t_2);
%     
%     Mn_theta = Bernstein_polynomial( 20, t_2 );
%     
%     Mn_eta = Bernstein_polynomial( 20, t_2 );
%     
%     Mn_phi = Bernstein_polynomial( 20, t_2 );
% 
%     
%     % Fit the shape functions to the symmetric maneuver data:
%     
%     X_CST_theta = zeros(500,length(Mn_theta(1,:)),nr_wb);
%     
%     X_CST_eta = zeros(500,length(Mn_eta(1,:)),nr_wb);
%     
%     X_CST_phi = zeros(500,length(Mn_phi(1,:)),nr_wb);
%     
%     for i = 1:nr_wb
%     
%     X_CST_theta(:,:,i) = ((X_theta(:,:,i)*[a_avg.theta_LR1; a_avg.theta_LR2])*ones(1,length(Mn_theta(1,:)))).*Mn_theta;
%     
%     X_CST_eta(:,:,i) = ((X_eta(:,:,i)*[a_avg.eta_LR1; a_avg.eta_LR2])*ones(1,length(Mn_eta(1,:)))).*Mn_eta;
%     
%     X_CST_phi(:,:,i) = ((X_phi(:,:,i)*[a_avg.phi_LR1; a_avg.phi_LR2])*ones(1,length(Mn_phi(1,:)))).*Mn_phi;
%     
%     end
% 
%     
%     b_theta_L = zeros(length(Mn_theta(1,:)),nr_wb);
%     
%     b_eta_L = zeros(length(Mn_eta(1,:)),nr_wb);
%     
%     b_phi_L = zeros(length(Mn_phi(1,:)),nr_wb);
%     
%     b_theta_R = zeros(length(Mn_theta(1,:)),nr_wb);
%     
%     b_eta_R = zeros(length(Mn_eta(1,:)),nr_wb);
%     
%     b_phi_R = zeros(length(Mn_phi(1,:)),nr_wb);
%     
%     
%     for i = 1:nr_wb
%     
%     Y_theta_sym_L = X_theta(:,:,i)*[a_sym.theta_L1(:,i); a_sym.theta_L2(:,i)];
%     
%     Y_eta_sym_L = X_eta(:,:,i)*[a_sym.eta_L1(:,i); a_sym.eta_L2(:,i)];
%     
%     Y_phi_sym_L = X_phi(:,:,i)*[a_sym.phi_L1(:,i); a_sym.phi_L2(:,i)];
% 
%     Y_theta_sym_R = X_theta(:,:,i)*[a_sym.theta_R1(:,i); a_sym.theta_R2(:,i)];
%     
%     Y_eta_sym_R = X_eta(:,:,i)*[a_sym.eta_R1(:,i); a_sym.eta_R2(:,i)];
%     
%     Y_phi_sym_R = X_phi(:,:,i)*[a_sym.phi_R1(:,i); a_sym.phi_R2(:,i)];
%     
% 
%     
%     b_theta_L(:,i) = (X_CST_theta(:,:,i)'*X_CST_theta(:,:,i))\(X_CST_theta(:,:,i)'*Y_theta_sym_L);
%     
%     b_eta_L(:,i) = (X_CST_eta(:,:,i)'*X_CST_eta(:,:,i))\(X_CST_eta(:,:,i)'*Y_eta_sym_L);
%     
%     b_phi_L(:,i) = (X_CST_phi(:,:,i)'*X_CST_phi(:,:,i))\(X_CST_phi(:,:,i)'*Y_phi_sym_L);
%     
%     b_theta_R(:,i) = (X_CST_theta(:,:,i)'*X_CST_theta(:,:,i))\(X_CST_theta(:,:,i)'*Y_theta_sym_R);
%     
%     b_eta_R(:,i) = (X_CST_eta(:,:,i)'*X_CST_eta(:,:,i))\(X_CST_eta(:,:,i)'*Y_eta_sym_R);
%     
%     b_phi_R(:,i) = (X_CST_phi(:,:,i)'*X_CST_phi(:,:,i))\(X_CST_phi(:,:,i)'*Y_phi_sym_R);
%     
%     end
    
%     figure()
%     surf(b_theta_L)
%     title('theta_L')
%     
%     figure()
%     surf(b_eta_L)
%     title('eta_L')
%     
%     figure()
%     surf(b_phi_L)
%     title('phi_L')
%     
%     figure()
%     surf(b_theta_R)
%     title('theta_R')
%     
%     figure()
%     surf(b_eta_R)
%     title('eta_R')
%     
%     figure()
%     surf(b_phi_R)
%     title('phi_R')
    
%     % Plot CST fits:
%     surf_dev_theta_L = zeros(nr_wb,500);
%     surf_dev_theta_R = zeros(nr_wb,500);
%     surf_dev_eta_L = zeros(nr_wb,500);
%     surf_dev_eta_R = zeros(nr_wb,500);
%     surf_dev_phi_L = zeros(nr_wb,500);
%     surf_dev_phi_R = zeros(nr_wb,500);
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%         
%         surf_dev_theta_L(k,:) = X_theta(:,:,k)*[a_dev.theta_L1(:,k); a_dev.theta_L2(:,k)];
%         
%         plot(t_2,radtodeg(X_theta(:,:,k)*[a_dev.theta_L1(:,k); a_dev.theta_L2(:,k)]),'Color',[0.5 0.5 0.5])
%         
%         
%     end
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%         
%         surf_dev_theta_R(k,:) = X_theta(:,:,k)*[a_dev.theta_R1(:,k); a_dev.theta_R2(:,k)];
%         
%         plot(t_2,radtodeg(X_theta(:,:,k)*[a_dev.theta_R1(:,k); a_dev.theta_R2(:,k)]),'Color',[0.5 0.5 0.5])
%         
%         
%     end
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%         
%         surf_dev_eta_L(k,:) = X_eta(:,:,k)*[a_dev.eta_L1(:,k); a_dev.eta_L2(:,k)];
%         
%         plot(t_2,radtodeg(X_eta(:,:,k)*[a_dev.eta_L1(:,k); a_dev.eta_L2(:,k)]),'Color',[0.5 0.5 0.5])
%         
%         
%     end
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%         
%         surf_dev_eta_R(k,:) = X_eta(:,:,k)*[a_dev.eta_R1(:,k); a_dev.eta_R2(:,k)];
%         
%         plot(t_2,radtodeg(X_eta(:,:,k)*[a_dev.eta_R1(:,k); a_dev.eta_R2(:,k)]),'Color',[0.5 0.5 0.5])
%         
%         
%     end
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%         
%         surf_dev_phi_L(k,:) = X_phi(:,:,k)*[a_dev.phi_L1(:,k); a_dev.phi_L2(:,k)];
%         
%         plot(t_2,radtodeg(X_phi(:,:,k)*[a_dev.phi_L1(:,k); a_dev.phi_L2(:,k)]),'Color',[0.5 0.5 0.5])
%         
%         
%     end
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%         
%         surf_dev_phi_R(k,:) = X_phi(:,:,k)*[a_dev.phi_R1(:,k); a_dev.phi_R2(:,k)];
%         
%         plot(t_2,radtodeg(X_phi(:,:,k)*[a_dev.phi_R1(:,k); a_dev.phi_R2(:,k)]),'Color',[0.5 0.5 0.5])
%         
%         
%     end
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%         surf_dev_theta_L(k,:) = X_theta(:,:,k)*[a_dev.theta_L1(:,k); a_dev.theta_L2(:,k)];
%         
%         plot(t_2,radtodeg(X_theta(:,:,k)*[a_dev.theta_L1(:,k)+a_dev.theta_R1(:,k); a_dev.theta_L2(:,k)+a_dev.theta_R2(:,k)]),'Color',[0.5 0.5 0.5])
%         
%         
%     end
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%         
%         plot(t_2,radtodeg(X_eta(:,:,k)*[a_dev.eta_L1(:,k)+a_dev.eta_R1(:,k); a_dev.eta_L2(:,k)+a_dev.eta_R2(:,k)]),'Color',[0.5 0.5 0.5])
%         
%         
%     end
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%         
%         plot(t_2,radtodeg(X_phi(:,:,k)*[a_dev.phi_L1(:,k)+a_dev.phi_R1(:,k); a_dev.phi_L2(:,k)+a_dev.phi_R2(:,k)]),'Color',[0.5 0.5 0.5])
%         
%         
%     end
%     hold off

%     surf_dev_theta_dot_L = zeros(nr_wb,500);
%     surf_dev_theta_dot_R = zeros(nr_wb,500);
%     surf_dev_eta_dot_L = zeros(nr_wb,500);
%     surf_dev_eta_dot_R = zeros(nr_wb,500);
%     surf_dev_phi_dot_L = zeros(nr_wb,500);
%     surf_dev_phi_dot_R = zeros(nr_wb,500);
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%         
%         surf_dev_theta_dot_L(k,:) = X_dot_theta*[a_dev.theta_L1(:,k); a_dev.theta_L2(:,k)];
%         
%         plot(t_2,radtodeg(X_dot_theta*[a_dev.theta_L1(:,k); a_dev.theta_L2(:,k)]),'Color',[0.5 0.5 0.5])
%         
%         
%     end
%     title('dev theta_dot_L')
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%         
%         surf_dev_theta_dot_R(k,:) = X_dot_theta*[a_dev.theta_R1(:,k); a_dev.theta_R2(:,k)];
%         
%         plot(t_2,radtodeg(X_dot_theta*[a_dev.theta_R1(:,k); a_dev.theta_R2(:,k)]),'Color',[0.5 0.5 0.5])
%         
%         
%     end
%     title('dev theta_dot_R')
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%         
%         surf_dev_eta_dot_L(k,:) = X_dot_eta*[a_dev.eta_L1(:,k); a_dev.eta_L2(:,k)];
%         
%         plot(t_2,radtodeg(X_dot_eta*[a_dev.eta_L1(:,k); a_dev.eta_L2(:,k)]),'Color',[0.5 0.5 0.5])
%         
%         
%     end
%     title('dev eta_dot_L')
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%         
%         surf_dev_eta_dot_R(k,:) = X_dot_eta*[a_dev.eta_R1(:,k); a_dev.eta_R2(:,k)];
%         
%         plot(t_2,radtodeg(X_dot_eta*[a_dev.eta_R1(:,k); a_dev.eta_R2(:,k)]),'Color',[0.5 0.5 0.5])
%         
%         
%     end
%     title('dev eta_dot_R')
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%         
%         surf_dev_phi_dot_L(k,:) = X_dot_phi*[a_dev.phi_L1(:,k); a_dev.phi_L2(:,k)];
%         
%         plot(t_2,radtodeg(X_dot_phi*[a_dev.phi_L1(:,k); a_dev.phi_L2(:,k)]),'Color',[0.5 0.5 0.5])
%         
%         
%     end
%     title('dev phi_dot_L')
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%         
%         surf_dev_phi_dot_R(k,:) = X_dot_phi*[a_dev.phi_R1(:,k); a_dev.phi_R2(:,k)];
%         
%         plot(t_2,radtodeg(X_dot_phi*[a_dev.phi_R1(:,k); a_dev.phi_R2(:,k)]),'Color',[0.5 0.5 0.5])
%         
%         
%     end
%     title('dev phi_dot_R')
%     hold off
%     
%     figure()
%     surf(surf_dev_theta_L)
%     
%     figure()
%     surf(surf_dev_theta_R)
%     
%     figure()
%     surf(surf_dev_eta_L)
%     
%     figure()
%     surf(surf_dev_eta_R)
%     
%     figure()
%     surf(surf_dev_phi_L)
%     
%     figure()
%     surf(surf_dev_phi_R)
%     
%     
%     figure()
%     surf(surf_dev_theta_dot_L)
%     
%     figure()
%     surf(surf_dev_theta_dot_R)
%     
%     figure()
%     surf(surf_dev_eta_dot_L)
%     
%     figure()
%     surf(surf_dev_eta_dot_R)
%     
%     figure()
%     surf(surf_dev_phi_dot_L)
%     
%     figure()
%     surf(surf_dev_phi_dot_R)
%     
%     
% %     figure()
% %     hold on
% %     for k = 1:nr_wb
% %         
% %         plot(t_2,radtodeg(X_ddot_theta*[a_dev.theta_L1(:,k); a_dev.theta_L2(:,k)]),'Color',[0.5 0.5 0.5])
% %         
% %         
% %     end
% %     title('dev theta_ddot_L')
% %     hold off
% %     
% %     figure()
% %     hold on
% %     for k = 1:nr_wb
% %         
% %         plot(t_2,radtodeg(X_ddot_theta*[a_dev.theta_R1(:,k); a_dev.theta_R2(:,k)]),'Color',[0.5 0.5 0.5])
% %         
% %         
% %     end
% %     title('dev theta_ddot_R')
% %     hold off
% %     
% %     figure()
% %     hold on
% %     for k = 1:nr_wb
% %         
% %         plot(t_2,radtodeg(X_ddot_eta*[a_dev.eta_L1(:,k); a_dev.eta_L2(:,k)]),'Color',[0.5 0.5 0.5])
% %         
% %         
% %     end
% %     title('dev eta_ddot_L')
% %     hold off
% %     
% %     figure()
% %     hold on
% %     for k = 1:nr_wb
% %         
% %         plot(t_2,radtodeg(X_ddot_eta*[a_dev.eta_R1(:,k); a_dev.eta_R2(:,k)]),'Color',[0.5 0.5 0.5])
% %         
% %         
% %     end
% %     title('dev eta_ddot_R')
% %     hold off
% %     
% %     figure()
% %     hold on
% %     for k = 1:nr_wb
% %         
% %         plot(t_2,radtodeg(X_ddot_phi*[a_dev.phi_L1(:,k); a_dev.phi_L2(:,k)]),'Color',[0.5 0.5 0.5])
% %         
% %         
% %     end
% %     title('dev phi_ddot_L')
% %     hold off
% %     
% %     figure()
% %     hold on
% %     for k = 1:nr_wb
% %         
% %         plot(t_2,radtodeg(X_ddot_phi*[a_dev.phi_R1(:,k); a_dev.phi_R2(:,k)]),'Color',[0.5 0.5 0.5])
% %         
% %         
% %     end
% %     title('dev phi_ddot_R')
% %     hold off
%     
%     
% %     figure()
% %     surf(surf_dev_theta_L)
% %     
% %     
% %     
% %     for k = 1:nr_wb
% %         
% %         figure(k)
% %         plot(t_2,X_theta(:,:,k)*[a_sym.theta_L1(:,k); a_sym.theta_L2(:,k)],'b')
% %         hold on
% %         %plot(t_2,X_CST_theta(:,:,k)*b_theta_L(:,k),'r')
% %         plot(t_2,X_CST_theta(:,:,k)*ones(length(b_theta_L(:,1)),1),'g')
% %         legend('symmetric maneuver','average wingbeat')
% %         hold off
% %         
% %     end
%         
%     
%     pause
%     
%     close all

end

