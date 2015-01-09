 function [a_Bernstein] = Piecewise_B_spline_fit(Y,n_pol_Bern,n_pol_Leg,wb_loc,a_Leg)


    nr_wb = length(wb_loc);
    
    Y_length = wb_loc(nr_wb,2)-wb_loc(1,1)+1;
    
    % Determine local polynomial fit:
    
    n_deriv = 3;
    
    beta_local = zeros(nr_wb*(n_pol_Bern+1),1);
    
    for i = 1:nr_wb
        
        m = wb_loc(i,2)-wb_loc(i,1)+1;
        
        xt1 = 0:(1/(m-1)):1;
        
        xt2 = -1:(2/(m-1)):1;

        PN_t = Bernstein_polynomial( n_pol_Bern, n_deriv, xt1 );
        
        PN_Leg = Legendre_polynomial( n_pol_Leg, n_deriv, xt2 );
        
%         PN_Leg2 = Legendre_polynomial( n_pol_Bern, n_deriv, xt2 );
        
        X_local = diag(PN_Leg(:,:,1)'*a_Leg(:,i))*PN_t(:,:,1)'
        
        Y_local = Y((wb_loc(i,1):wb_loc(i,2)));
        
        t = ((i-1)*(n_pol_Bern+1)+1):(i*(n_pol_Bern+1));
        
        beta_local(t) = (X_local'*X_local)\(X_local'*Y_local);

        beta_local(t)
        
        figure()
        plot(X_local*beta_local(t),'b')
        hold on
        plot(Y_local,'r')
        hold off
        
        
%         figure()
%         plot(PN_Leg(:,:,2)'*beta_local(t),'b')


        
        pause
        
    end
    

    
    
%     % Set up global matrix X
%     
%     X = zeros(Y_length,nr_wb*(n_pol+1));
%     
%     for j = 1:nr_wb
%         
%         m = wb_loc(j,2)-wb_loc(j,1)+1;
%         
%         xt = 0:(1/(m-1)):1;
%         
%         PN_t = Bernstein_polynomial( n_pol, n_deriv, xt );
%         
%         a = (wb_loc(j,1)-wb_loc(1,1)+1):(wb_loc(j,2)-wb_loc(1,1)+1);
%         
%         b = ((n_pol+1)*(j-1)+1):((n_pol+1)*j);
%         
%         X(a,b) = PN_t(:,:,1)';
%         
%         if j > 1 && j < nr_wb
%             
%             m_1 = wb_loc(j,2)-wb_loc(j,1);
%         
%             xt_1 = 0:(1/(m_1-1)):1;
%         
%             PN_t_1 = Bernstein_polynomial( n_pol, n_deriv, xt_1 );        
%             
%             X((wb_loc(j,1)-wb_loc(1,1)+1),b) = 0.5*X((wb_loc(j,1)-wb_loc(1,1)+1),b);
%             
%             X((wb_loc(j,2)-wb_loc(1,1)+1),b) = 0.5*X((wb_loc(j,2)-wb_loc(1,1)+1),b);
%             
%         elseif j == 1
%          
% 
%             X((wb_loc(j,2)-wb_loc(1,1)+1),b) = 0.5*X((wb_loc(j,2)-wb_loc(1,1)+1),b);
%         
%         end
%             
%         if j == nr_wb
%                   
%             
%             X((wb_loc(j,1)-wb_loc(1,1)+1),b) = 0.5*X((wb_loc(j,1)-wb_loc(1,1)+1),b);
% 
%         end       
%         
%     end
%     
% 
%     
%     % set up matrix R
%     
%     R = zeros((n_deriv+1)*(nr_wb+1),nr_wb*(n_pol+1));
%     
%     R_L = zeros(n_deriv+1,n_pol+1);
%     
%     R_R = zeros(n_deriv+1,n_pol+1);
%     
%     PN_temp = Bernstein_polynomial(n_pol,n_deriv,0:1);
%     
%     for q = 1:(n_deriv+1)
%         
%         R_L(q,:) = PN_temp(:,1,q);
%         
%         R_R(q,:) = PN_temp(:,2,q);
%         
%     end
%     
%     for k = 1:(nr_wb+1)
%         
%         a = ((k-1)*(n_deriv+1)+1):(k*(n_deriv+1));
%         
%         if k == 1
%             
%             b = ((k-1)*(n_pol+1)+1):(k*(n_pol+1));
%             
%             R(a,b) = R_L;
%             
%         elseif k == nr_wb+1
%             
%             c = ((k-2)*(n_pol+1)+1):((k-1)*(n_pol+1));
%             
%             R(a,c) = R_R;
%             
%         else
%             
%             b = ((k-1)*(n_pol+1)+1):(k*(n_pol+1));
%             
%             c = ((k-2)*(n_pol+1)+1):((k-1)*(n_pol+1));
%             
%             R(a,b) = R_L;
%             
%             R(a,c) = -R_R;
%             
%         end
%         
%         
%     end
%     
%     
%     % Find restriction vector V_r:
%     
%     V_r = zeros((n_deriv+1)*(nr_wb+1),1);
%     
%     V_r(1) = PN_temp(:,1,1)'*beta_local(1:(n_pol+1));
%     V_r(2) = PN_temp(:,1,2)'*beta_local(1:(n_pol+1));
%     V_r(3) = PN_temp(:,1,3)'*beta_local(1:(n_pol+1));
%     V_r(4) = PN_temp(:,1,4)'*beta_local(1:(n_pol+1));
%     V_r((n_deriv+1)*nr_wb+1) = PN_temp(:,3,1)'*beta_local(((nr_wb-1)*(n_pol+1)+1):end);
%     V_r((n_deriv+1)*nr_wb+2) = PN_temp(:,3,2)'*beta_local(((nr_wb-1)*(n_pol+1)+1):end);
%     V_r((n_deriv+1)*nr_wb+3) = PN_temp(:,3,3)'*beta_local(((nr_wb-1)*(n_pol+1)+1):end);
%     V_r((n_deriv+1)*nr_wb+4) = PN_temp(:,3,4)'*beta_local(((nr_wb-1)*(n_pol+1)+1):end);
% 
%     
%     % Calculate restricted polynomial fit:
% 
%     
%     XTX_inv = inv(X'*X);
% 
%     
%     beta_star = beta_local - XTX_inv*R'*((R*XTX_inv*R')\(R*beta_local-V_r));
% 
%    
%     figure()
%     plot(radtodeg(X*beta_star),'b')
%     hold on
%     plot(radtodeg(Y(wb_loc(1,1):wb_loc(nr_wb,2))),'r')
%     legend('fit','data')
%     hold off
%     
%     figure()
%     plot(radtodeg(X_dot*beta_star),'b')
%   
%     figure()
%     plot(radtodeg(X_ddot*beta_star),'b')
%     
%     beta_element = zeros(n_pol+1,nr_wb);
%     
%     for l = 1:nr_wb
%         
%         beta_element(:,l) = beta_star(((l-1)*(n_pol+1)+1):(l*(n_pol+1)));
%         
%     end
%     
%     a_fit = beta_element;
%     
%     figure()
%     surf(beta_element)
        
     
    
end