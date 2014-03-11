function [ a_fit_1, a_fit_2] = average_fit_periodic_updowncut(a_1, a_2,n_pol,Rds)

n_wb = 11;
M = 2000;
m1 = round(M*Rds);
m2 = M-m1+1;
    


    % Use a polynomial fit to construct an average wingbeat:
        
    n_deriv = 4;
    
    % Determine local polynomial fit---------------------------------------

    beta_local1 = a_1';
    
    beta_local2 = a_2';


    
    %----------------------------------------------------------------------
    
    
    % Create matrix R------------------------------------------------------
    
    
    PN_t2 = Legendre_polynomial( n_pol, n_deriv, -1:1 );
    
    PN_L = zeros(n_deriv+1,n_pol+1);
    
    PN_R = zeros(n_deriv+1,n_pol+1);
    
    for k = 1:((n_deriv+1))
        
        PN_L(k,:) = PN_t2(:,1,k)';
    
        PN_R(k,:) = PN_t2(:,3,k)';
        
    end
    
    
%     R = zeros((n_deriv+1)*2*n_wb,2*n_wb*(n_pol+1));
%     
%    
%     R_L = PN_L;
%     
%     R_R = PN_R;   
%     
%     
%     for l = 1:(n_wb*2)
%         
%         a = ((l-1)*(n_deriv+1)+1):(l*(n_deriv+1));
%         
%         if l == 1
%             
%             b = ((l-1)*(n_pol+1)+1):(l*(n_pol+1));
%             
%             R(a,b) = R_L;
%             
%         elseif l == n_wb+1
%             
%             c = ((l-2)*(n_pol+1)+1):((l-1)*(n_pol+1));
%             
%             R(a,c) = R_R;
%             
%         else
%             
%             b = ((l-1)*(n_pol+1)+1):(l*(n_pol+1));
%             
%             c = ((l-2)*(n_pol+1)+1):((l-1)*(n_pol+1));
%             
%             R(a,b) = R_L;
%             
%             R(a,c) = -R_R;
%             
%         end
%         
%         
%     end

    R = zeros((n_deriv+1)*(2*n_wb+1),(n_pol+1)*2*n_wb);
    
    R_L = PN_L;
    
    R_R = PN_R;  
    
    for k = 1:(2*n_wb)
        
        a = ((k-1)*(n_deriv+1)+1):(k*(n_deriv+1));
        
        if k == 1
            
            c = 1:(n_pol+1);
            
            R(a,c) = R_L;
            
        elseif k == (n_wb*2)
            
            c = ((n_pol+1)*(2*n_wb-2)+1):((n_pol+1)*(2*n_wb-1));
            
            d = ((n_pol+1)*(2*n_wb-1)+1):((n_pol+1)*2*n_wb);
            
            e = (2*n_wb*(n_deriv+1)+1):((2*n_wb+1)*(n_deriv+1));
            
            R(a,c) = -R_R;
            
            R(a,d) = R_L;
            
            R(e,d) = R_R;
            
        elseif k > 1 && k < (n_wb*2)
            
            c = ((k-2)*(n_pol+1)+1):((k-1)*(n_pol+1));
            
            d = ((k-1)*(n_pol+1)+1):(k*(n_pol+1));
            
            R(a,c) = -R_R;
            
            R(a,d) = R_L;
            
        end
        
    end
   
    
    beta_local = [beta_local1; beta_local2; beta_local1; beta_local2; beta_local1; beta_local2; beta_local1; ...
                  beta_local2; beta_local1; beta_local2; beta_local1; beta_local2; beta_local1; beta_local2; ...
                  beta_local1; beta_local2; beta_local1; beta_local2; beta_local1; beta_local2; beta_local1; beta_local2];

    
    % Apply restricted least square fit on the boundaries:
    
    V_r = zeros((n_deriv+1)*(2*n_wb+1),1);
    
    
    V_r(1:(n_deriv+1)) = 0.5*(PN_L*beta_local1+PN_R*beta_local2);
    
    V_r(((2*n_wb)*(n_deriv+1)+1):((2*n_wb+1)*(n_deriv+1))) = 0.5*(PN_L*beta_local2+PN_R*beta_local1);
    
    
    % Set up global matrix X
    
%     X4 = zeros(20*n_wb+1,n_wb*(n_pol+1));
    X4 = zeros(M*n_wb-n_wb,2*n_wb*(n_pol+1));
    
    clear a b c d
    
    for j = 1:n_wb
        
%         m = 21;
%         
%         xt = -1:(2/(m-1)):1;
%         
%         PN_t = Legendre_polynomial( n_pol, n_deriv, xt );
        
        x1 = -1:(2/(m1-1)):1;
        x2 = -1:(2/(m2-1)):1;
        
        X1 = Legendre_polynomial( n_pol, n_deriv, x1 );
        X2 = Legendre_polynomial( n_pol, n_deriv, x2 );
        
        a = (1+(M-1)*(j-1)):(m1+(M-1)*(j-1));
        
        b = (m1+(M-1)*(j-1)):(m1-1+m2+(M-1)*(j-1));
        
        c = (2*(n_pol+1)*(j-1)+1):(2*(n_pol+1)*(j-1)+n_pol+1);
        
        d = (2*(n_pol+1)*(j-1)+n_pol+2):(2*(n_pol+1)*(j-1)+2*(n_pol+1));
        
        X4(a,c) = X1(:,:,1)';
        
        X4(b,d) = X2(:,:,1)';
        
        if j > 1 && j < n_wb
                     
            
            X4(a(1),c) = 0.5*X4(a(1),c);
            
            X4(a(end),c) = 0.5*X4(a(end),c);
            
            X4(b(1),d) = 0.5*X4(b(1),d);
            
            X4(b(end),d) = 0.5*X4(b(end),d);
            
        elseif j == 1
         

            X4(a(end),c) = 0.5*X4(a(end),c);
            
            X4(b(1),d) = 0.5*X4(b(1),d);
            
            X4(b(end),d) = 0.5*X4(b(end),d);
        
        end
            
        if j == n_wb
                  
            
            X4(a(1),c) = 0.5*X4(a(1),c);
            
            X4(a(end),c) = 0.5*X4(a(end),c);
            
            X4(b(1),d) = 0.5*X4(b(1),d);

        end       
        
    end


    XTX_inv = inv(X4'*X4);

   
    beta_star = beta_local - XTX_inv*R'*((R*XTX_inv*R')\(R*beta_local-V_r));
    

    a_fit_1 = beta_star(((n_pol+1)*10+1):((n_pol+1)*11));
    
    
    a_fit_2 = beta_star(((n_pol+1)*11+1):((n_pol+1)*12));
    
 

end


