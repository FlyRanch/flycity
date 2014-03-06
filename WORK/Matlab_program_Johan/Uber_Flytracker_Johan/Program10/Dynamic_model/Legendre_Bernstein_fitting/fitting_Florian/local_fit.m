function [ a_fit_1, a_fit_2, ratio_1_2 ] = local_fit(Y,n_pol)

        n_deriv = 3;
        
        
        
        % Make local fit
        
        m = length(Y);
        
        m1 = ceil((m+1)/2);
        
        m2 = m-m1+1;
        
        xt1 = -1:(2/(m1-1)):1;
        
        xt2 = -1:(2/(m2-1)):1;
        
        Y1 = Y(1:m1);
        
        Y2 = Y(m1:end);
        
        PN_t1 = Legendre_polynomial( n_pol, n_deriv, xt1 );
        
        PN_t2 = Legendre_polynomial( n_pol, n_deriv, xt2 );
        
        X1 = PN_t1(:,:,1)';
        
        X2 = PN_t2(:,:,1)';

        a_fit_local_1 = (X1'*X1)\(X1'*Y1);
        
        a_fit_local_2 = (X2'*X2)\(X2'*Y2);
        
        ratio_1_2 = m1/(m2);

        a_fit_local_1_mirror = a_fit_local_1;
        
        a_fit_local_2_mirror = a_fit_local_2;
        
%         a_fit_local_1_mirror(2:end) = -a_fit_local_1(2:end);
%         
%         a_fit_local_2_mirror(2:end) = -a_fit_local_2(2:end);

        a_fit_local_1_mirror = a_fit_local_1;
        
        a_fit_local_2_mirror = a_fit_local_2;

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


    
    beta_local = [a_fit_local_1_mirror; a_fit_local_2_mirror; a_fit_local_1; a_fit_local_2; a_fit_local_1_mirror; a_fit_local_2_mirror; a_fit_local_1; ...
                  a_fit_local_2; a_fit_local_1_mirror; a_fit_local_2_mirror; a_fit_local_1; a_fit_local_2; a_fit_local_1_mirror; a_fit_local_2_mirror; ...
                  a_fit_local_1; a_fit_local_2; a_fit_local_1_mirror; a_fit_local_2_mirror; a_fit_local_1; a_fit_local_2; a_fit_local_1_mirror; a_fit_local_2_mirror];

    
    % Apply restricted least square fit on the boundaries:
    
    V_r = zeros((n_deriv+1)*22,1);
    
    V_r(1:(n_deriv+1)) = 0.5*(PN_L*a_fit_local_1_mirror+PN_R*a_fit_local_2_mirror);
    
    V_r((21*(n_deriv+1)+1):(22*(n_deriv+1))) = 0.5*(PN_L*a_fit_local_2_mirror+PN_R*a_fit_local_1_mirror);
    

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

