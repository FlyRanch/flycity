function [ a_fit_1, a_fit_2] = average_fit(a_1,a_2,n_pol,Weight1,Weight2,down_up)

    % Use a polynomial fit to construct an average wingbeat:
    

    
    n_deriv = 4;
    
    nr_avg_wb = 20;
    
    n = 100;
    
    m1 = round(down_up*n);
    
    m2 = n-m1+1;
    

    
    
    
    % Determine average local polynomial fit-------------------------------

%     beta_local1 = mean(a_1,2);
%     
%     beta_local2 = mean(a_2,2);
% 
% 



    beta_local1 = a_1*Weight1';

    beta_local2 = a_2*Weight2';

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
    
    
    R = zeros((n_deriv+1)*(2*nr_avg_wb+1),(2*nr_avg_wb)*(n_pol+1));
    
   
    R_L = PN_L;
    
    R_R = PN_R;   
    
    
    for k = 1:(2*nr_avg_wb+1)
        
        a = ((k-1)*(n_deriv+1)+1):(k*(n_deriv+1));
        
        if k == 1
            
            b = ((k-1)*(n_pol+1)+1):(k*(n_pol+1));
            
            R(a,b) = R_L;
            
        elseif k == (2*nr_avg_wb+1)
            
            c = ((k-2)*(n_pol+1)+1):((k-1)*(n_pol+1));
            
            R(a,c) = R_R;
            
        elseif k > 1 && k < (2*nr_avg_wb+1)
            
            b = ((k-1)*(n_pol+1)+1):(k*(n_pol+1));
            
            c = ((k-2)*(n_pol+1)+1):((k-1)*(n_pol+1));
            
            R(a,b) = R_L;
            
            R(a,c) = -R_R;
            
        end
        
        
    end    
    
    
    beta_local = zeros(nr_avg_wb*2*(n_pol+1),1);
    
    for k = 1:nr_avg_wb

        beta_local(((k-1)*(2*(n_pol+1))+1):(k*2*(n_pol+1))) = [beta_local1; beta_local2];
        
    end
    
    % Apply restricted least square fit on the boundaries:
    
    V_r = zeros((n_deriv+1)*(2*nr_avg_wb+1),1);
    
    V_r(1:(n_deriv+1)) = 0.5*(PN_L*beta_local1+PN_R*beta_local2);
    
    V_r(((n_deriv+1)*(2*nr_avg_wb)+1):((n_deriv+1)*(2*nr_avg_wb+1))) = 0.5*(PN_L*beta_local2+PN_R*beta_local1);
    
    
    % Set up global matrix X
    
    X4 = zeros(n*2*nr_avg_wb-nr_avg_wb+1,2*nr_avg_wb*(n_pol+1));
    
    for j = 1:(2*nr_avg_wb)
        
        if mod(j,2) == 1
                
            xt = -1:(2/(m1-1)):1;

            PN_t = Legendre_polynomial( n_pol, n_deriv, xt );
            
            a = ((j-1)*n-(j-1)+1):((j-1)*n-(j-1)+m1);
            
            b = ((j-1)*(n_pol+1)+1):(j*(n_pol+1));
            
        elseif mod(j,2) == 0
            
            xt = -1:(2/(m2-1)):1;

            PN_t = Legendre_polynomial( n_pol, n_deriv, xt );
            
            a = ((j-1)*n-(j-1)+m1):(j*n-(j-1));
            
            b = ((j-1)*(n_pol+1)+1):(j*(n_pol+1));
            
        end
        
        X4(a,b) = PN_t(:,:,1)';
        
        if j > 1 && j < (2*nr_avg_wb)                 
            
            X4(a(1),b) = 0.5*X4(a(1),b);
            
            X4(a(end),b) = 0.5*X4(a(end),b);
            
        elseif j == 1
         
            X4(a(end),b) = 0.5*X4(a(end),b);
        
        end
            
        if j == (2*nr_avg_wb)
                  
            X4(a(1),b) = 0.5*X4(a(1),b);

        end       
        
    end
    
    XTX_inv = inv(X4'*X4);

   
    beta_star = beta_local - XTX_inv*R'*((R*XTX_inv*R')\(R*beta_local-V_r));
    
    e = round(nr_avg_wb/2)*2;

    a_fit_1 = beta_star(((n_pol+1)*e+1):((n_pol+1)*(e+1)));
    
    
    a_fit_2 = beta_star(((n_pol+1)*(e+1)+1):((n_pol+1)*(e+2)));
% 
%     a_fit_1 = beta_local1;
%     
%     a_fit_2 = beta_local2;

    
    
end


