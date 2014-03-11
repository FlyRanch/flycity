function [ a_fit_1, a_fit_2] = average_fit(a_1, a_2,n_pol)


    % Use a polynomial fit to construct an average wingbeat:
        
    n_deriv = 4;
    
    % Determine local polynomial fit---------------------------------------

    beta_local1 = mean(a_1,1)';
    
    beta_local2 = mean(a_2,1)';


    
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
            
        elseif l == 22+1
            
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
    

    % Set up global matrix X
    
    X4 = zeros(20*22+1,22*(n_pol+1));
    
    for j = 1:22
        
        m = 21;
        
        xt = -1:(2/(m-1)):1;
        
        PN_t = Legendre_polynomial( n_pol, n_deriv, xt );
        
        a = (1+20*(j-1)):(1+20*j);
        
        b = ((n_pol+1)*(j-1)+1):((n_pol+1)*j);
        
        X4(a,b) = PN_t(:,:,1)';
        
        if j > 1 && j < 22
            
            m_1 = 21;
        
            xt_1 = -1:(2/(m_1-1)):1;
        
            PN_t_1 = Legendre_polynomial( n_pol, n_deriv, xt_1 );          
            
            X4((1+20*(j-1)),b) = 0.5*X4((1+20*(j-1)),b);
            
            X4((1+20*j),b) = 0.5*X4((1+20*j),b);
            
        elseif j == 1
         

            X4((1+20*j),b) = 0.5*X4((1+20*j),b);
        
        end
            
        if j == 22
                  
            
            X4((1+20*(j-1)),b) = 0.5*X4((1+20*(j-1)),b);

        end       
        
    end


    XTX_inv = inv(X4'*X4);

   
    beta_star = beta_local - XTX_inv*R'*((R*XTX_inv*R')\(R*beta_local-V_r));
    

    a_fit_1 = beta_star(((n_pol+1)*10+1):((n_pol+1)*11));
    
    
    a_fit_2 = beta_star(((n_pol+1)*11+1):((n_pol+1)*12));

end


