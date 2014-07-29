function [ a_fit_1, a_fit_2] = average_fit(a_1,a_2,n_pol,wb_loc,Weight1,Weight2)
% function [ a_fit_1, a_fit_2] = average_fit(Y,n_pol,wb_loc,Weight1,Weight2)


    % Use a polynomial fit to construct an average wingbeat:
    
    nr_wb = length(wb_loc)/2;
    
    n_deriv = 4;
    
    % Determine local polynomial fit---------------------------------------
    
    
%     % Set up local Y matrix:
%     
%     Y_local1 = Y(wb_loc(1,1):wb_loc(1,2));
%     
%     Y_local2 = Y(wb_loc(2,1):wb_loc(2,2));
%     
%     for i = 3:(nr_wb*2)
%         
%         if mod(i,2) == 1
%             
%             Y_local1 = [ Y_local1; Y(wb_loc(i,1):wb_loc(i,2)) ];
%             
%         elseif mod(i,2) == 0
%             
%             Y_local2 = [ Y_local2; Y(wb_loc(i,1):wb_loc(i,2)) ];
%             
%         end
%         
%     end    
%     
%     
%     
%     % Determine matrix X and W
%     
%     X1 = zeros(length(Y_local1),n_pol+1);
%     
%     W1 = zeros(length(Y_local1),length(Y_local1));
%     
%     X2 = zeros(length(Y_local2),n_pol+1);
%     
%     W2 = zeros(length(Y_local2),length(Y_local2));
%     
%     i_temp1 = 0;
%     
%     i_temp2 = 0;
%     
%     for i = 1:(2*nr_wb)
%         
%         m = wb_loc(i,2)-wb_loc(i,1)+1;
% 
%         xt = -1:(2/(m-1)):1;
% 
%         PN_t1 = Legendre_polynomial( n_pol, n_deriv, xt );
%         
%         if mod(i,2) == 1
%             
%             a = (i_temp1+1):(i_temp1+m);
%             
%             X1(a,:) = PN_t1(:,:,1)';
%             
%             W1(a,a) = diag(ones(m,1)*Weight1(ceil(i/2)));
%             
%             i_temp1 = i_temp1+m;
%             
%         elseif mod(i,2) == 0
%             
%             a = (i_temp2+1):(i_temp2+m);
%             
%             X2(a,:) = PN_t1(:,:,1)';
%             
%             W2(a,a) = diag(ones(m,1)*Weight2(ceil(i/2)));            
% 
%             i_temp2 = i_temp2+m;
%         end
%         
%     end
%     
%     clear i_temp1 i_temp2
%     
%     
%     beta_local1 = (X1'*W1*X1)\(X1'*W1*Y_local1);
%     
%     beta_local2 = (X2'*W2*X2)\(X2'*W2*Y_local2);


beta_local1 = a_1*Weight1';

beta_local2 = a_2*Weight2';


%     figure()
%     hold on
%     for k = 1:(nr_wb*2)
%         t = (wb_loc(k,1)-wb_loc(1,1)+1):(wb_loc(k,2)-wb_loc(1,1)+1);
%         plot(t,Y(wb_loc(k,1):wb_loc(k,2)),'Color',[0.5 0.5 0.5])
%         if mod(k,2) == 1
%             m = wb_loc(k,2)-wb_loc(k,1)+1;
%             xt = -1:(2/(m-1)):1;
%             PN_t = Legendre_polynomial( n_pol, n_deriv, xt );
%             plot(t,PN_t(:,:,1)'*beta_local1);            
%         elseif mod(k,2) == 0
%             m = wb_loc(k,2)-wb_loc(k,1)+1;
%             xt = -1:(2/(m-1)):1;
%             PN_t = Legendre_polynomial( n_pol, n_deriv, xt );
%             plot(t,PN_t(:,:,1)'*beta_local2);
%         end
%     end
%     hold off
%     
%     pause

    
    %----------------------------------------------------------------------
    
    
    % Create matrix R------------------------------------------------------
    
    
    PN_t2 = Legendre_polynomial( n_pol, n_deriv, -1:1 );
    
    PN_L = zeros(n_deriv+1,n_pol+1);
    
    PN_R = zeros(n_deriv+1,n_pol+1);
    
    for k = 1:((n_deriv+1))
        
        PN_L(k,:) = PN_t2(:,1,k)';
    
        PN_R(k,:) = PN_t2(:,3,k)';
        
    end
    
    
    R = zeros((n_deriv+1)*22,22*(n_pol+1));
    
   
    R_L = PN_L;
    
    R_R = PN_R;   
    
    
    for l = 1:22
        
        a = ((l-1)*(n_deriv+1)+1):(l*(n_deriv+1));
        
        if l == 1
            
            b = ((l-1)*(n_pol+1)+1):(l*(n_pol+1));
            
            R(a,b) = R_L;
            
        elseif l == nr_wb+1
            
            c = ((l-2)*(n_pol+1)+1):((l-1)*(n_pol+1));
            
            R(a,c) = R_R;
            
        else
            
            b = ((l-1)*(n_pol+1)+1):(l*(n_pol+1));
            
            c = ((l-2)*(n_pol+1)+1):((l-1)*(n_pol+1));
            
            R(a,b) = R_L;
            
            R(a,c) = -R_R;
            
        end
        
        
    end
    

    
    beta_local = [beta_local1; beta_local2; beta_local1; beta_local2; beta_local1; beta_local2; beta_local1; ...
                  beta_local2; beta_local1; beta_local2; beta_local1; beta_local2; beta_local1; beta_local2; ...
                  beta_local1; beta_local2; beta_local1; beta_local2; beta_local1; beta_local2; beta_local1; beta_local2];

    
    % Apply restricted least square fit on the boundaries:
    
    V_r = zeros((n_deriv+1)*22,1);
    
    V_r(1:(n_deriv+1)) = 0.5*(PN_L*beta_local1+PN_R*beta_local2);
    
    V_r((21*(n_deriv+1)+1):(22*(n_deriv+1))) = 0.5*(PN_L*beta_local2+PN_R*beta_local1);
    
    
    PN_t3 = Legendre_polynomial( n_pol, n_deriv, -1:0.1:1 );
    
    X3 = PN_t3(:,:,1)';
    
    zero_mat = zeros(21,n_pol+1);
    
    X4 = [X3 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat X3 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat X3 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat X3 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat X3 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat X3 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X3 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X3 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X3 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X3 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X3 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X3 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X3 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X3 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X3 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X3 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X3 zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X3 zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X3 zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X3 zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X3 zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X3  ];
  

%     % Set up global matrix X
%     
%     X4 = zeros(20*22+1,22*(n_pol+1));
%     
%     for j = 1:22
%         
%         m = 21;
%         
%         xt = -1:(2/(m-1)):1;
%         
%         PN_t = Legendre_polynomial( n_pol, n_deriv, xt );
%         
%         a = (1+20*(j-1)):(1+20*j);
%         
%         b = ((n_pol+1)*(j-1)+1):((n_pol+1)*j);
%         
%         X4(a,b) = PN_t(:,:,1)';
%         
%         if j > 1 && j < 22
%             
%             m_1 = wb_loc(j,2)-wb_loc(j,1);
%         
%             xt_1 = -1:(2/(m_1-1)):1;
%         
%             PN_t_1 = Legendre_polynomial( n_pol, n_deriv, xt_1 );          
%             
%             X4((1+20*(j-1)),b) = 0.5*X4((1+20*(j-1)),b);
%             
%             X4((1+20*j),b) = 0.5*X4((1+20*j),b);
%             
%         elseif j == 1
%          
% 
%             X4((1+20*j),b) = 0.5*X4((1+20*j),b);
%         
%         end
%             
%         if j == 22
%                   
%             
%             X4((1+20*(j-1)),b) = 0.5*X4((1+20*(j-1)),b);
% 
%         end       
%         
%     end
    
    
    XTX_inv = inv(X4'*X4);

   
    beta_star = beta_local - XTX_inv*R'*((R*XTX_inv*R')\(R*beta_local-V_r));
    

    a_fit_1 = beta_star(((n_pol+1)*10+1):((n_pol+1)*11));
    
    
    a_fit_2 = beta_star(((n_pol+1)*11+1):((n_pol+1)*12));
        
    
    
%     PN_t4 = Legendre_polynomial( n_pol, n_deriv, -1:0.01:1 );
%     figure()
%     hold on
%     for k = 1:(nr_wb*2)
%         t = (wb_loc(k,1)-wb_loc(1,1)+1):(wb_loc(k,2)-wb_loc(1,1)+1);
%         plot(t,Y(wb_loc(k,1):wb_loc(k,2)),'Color',[0.5 0.5 0.5])
%         if mod(k,2) == 1
%             m = wb_loc(k,2)-wb_loc(k,1)+1;
%             xt = -1:(2/(m-1)):1;
%             PN_t = Legendre_polynomial( n_pol, n_deriv, xt );
%             plot(t,PN_t(:,:,1)'*a_fit_1);            
%         elseif mod(k,2) == 0
%             m = wb_loc(k,2)-wb_loc(k,1)+1;
%             xt = -1:(2/(m-1)):1;
%             PN_t = Legendre_polynomial( n_pol, n_deriv, xt );
%             plot(t,PN_t(:,:,1)'*a_fit_2);
%         end
%     end
%     hold off
%     
%     figure()
%     hold on
%     for k = 1:(nr_wb*2)
%         if mod(k,2) == 1
%             m = wb_loc(k,2)-wb_loc(k,1)+1;
%             plot(-1:(1/(m-1)):0,Y(wb_loc(k,1):wb_loc(k,2)),'Color',[0.5 0.5 0.5]);            
%         elseif mod(k,2) == 0
%             m = wb_loc(k,2)-wb_loc(k,1)+1;
%             plot(0:(1/(m-1)):1,Y(wb_loc(k,1):wb_loc(k,2)),'Color',[0.5 0.5 0.5]);
%         end
%     end
%     plot(-1:0.005:0,PN_t4(:,:,1)'*a_fit_1,'r')
%     plot(0:0.005:1,PN_t4(:,:,1)'*a_fit_2,'r')
%     hold off
% %     
% %     figure()
% %     hold on
% %     plot(-1:0.005:0,PN_t4(:,:,2)'*a_fit_1)
% %     plot(0:0.005:1,PN_t4(:,:,2)'*a_fit_2)
% %     hold off
% %     
% %     figure()
% %     hold on
% %     plot(-1:0.005:0,PN_t4(:,:,3)'*a_fit_1)
% %     plot(0:0.005:1,PN_t4(:,:,3)'*a_fit_2)
% %     hold off
%     
%     
%     
%     
%     pause
%     
%     close all
end


