 function [a_fit] = Piecewise_polynomial_fit(Y,n_pol,wb_loc,down_up_ratio)


    nr_wb = size(wb_loc,1);
    
    Y_length = wb_loc(nr_wb,2)-wb_loc(1,1)+1;
    
    % Determine local polynomial fit:
    
    beta_local = zeros(nr_wb*(n_pol+1),1);
    
    n_deriv = 4;
    
       
    for i = 1:nr_wb
        
        m = wb_loc(i,2)-wb_loc(i,1)+1;
        
        xt = -1:(2/(m-1)):1;
        
        down_up = down_up_ratio(i);
        
        PN_t = Legendre_polynomial( n_pol, n_deriv, xt, down_up);
        
        X_local = PN_t(:,:,1)';
        
        Y_local = Y((wb_loc(i,1):wb_loc(i,2)));
       
        t = ((i-1)*(n_pol+1)+1):(i*(n_pol+1));

        beta_local(t) = (X_local'*X_local)\(X_local'*Y_local);
        
    end
    
    
    % Set up global matrix X
    
    X = zeros(Y_length,nr_wb*(n_pol+1));
    
    for j = 1:nr_wb
        
        m = wb_loc(j,2)-wb_loc(j,1)+1;
        
        xt = -1:(2/(m-1)):1;
        
        down_up = down_up_ratio(j);
        
        PN_t = Legendre_polynomial( n_pol, n_deriv, xt, down_up);
        
        a = (wb_loc(j,1)-wb_loc(1,1)+1):(wb_loc(j,2)-wb_loc(1,1)+1);
        
        b = ((n_pol+1)*(j-1)+1):((n_pol+1)*j);
        
        X(a,b) = PN_t(:,:,1)';
        
        if j > 1 && j < nr_wb
            
            m_1 = wb_loc(j,2)-wb_loc(j,1);
        
            xt_1 = -1:(2/(m_1-1)):1;        
            
            X((wb_loc(j,1)-wb_loc(1,1)+1),b) = 0.5*X((wb_loc(j,1)-wb_loc(1,1)+1),b);
            
            X((wb_loc(j,2)-wb_loc(1,1)+1),b) = 0.5*X((wb_loc(j,2)-wb_loc(1,1)+1),b);
            
        elseif j == 1
         

            X((wb_loc(j,2)-wb_loc(1,1)+1),b) = 0.5*X((wb_loc(j,2)-wb_loc(1,1)+1),b);
        
        end
            
        if j == nr_wb
                  
            
            X((wb_loc(j,1)-wb_loc(1,1)+1),b) = 0.5*X((wb_loc(j,1)-wb_loc(1,1)+1),b);

        end       
        
    end
    

    
    % set up matrix R
    
    R = zeros((n_deriv+1)*(nr_wb+1),nr_wb*(n_pol+1));
    
    PN_temp_start       = zeros(n_deriv+1,n_pol+1);
    PN_temp_end         = zeros(n_deriv+1,n_pol+1);
    
    for k = 1:(nr_wb+1)
        
                
        a = ((k-1)*(n_deriv+1)+1):(k*(n_deriv+1));
        
                
        if k == 1
            
            down_up_L = down_up_ratio(k);
        
            R_L = zeros(n_deriv+1,n_pol+1);

            PN_temp_L = Legendre_polynomial(n_pol,n_deriv,-1:1,down_up_L);

            for q = 1:(n_deriv+1)

                R_L(q,:) = PN_temp_L(:,1,q);

            end
            
            PN_temp_start = R_L;
            
            b = ((k-1)*(n_pol+1)+1):(k*(n_pol+1));
            
            R(a,b) = R_L;
            
        elseif k == nr_wb+1
                        
            down_up_R = down_up_ratio(k-1);

            R_R = zeros(n_deriv+1,n_pol+1);
            
            PN_temp_R = Legendre_polynomial(n_pol,n_deriv,-1:1,down_up_R);

            for q = 1:(n_deriv+1)

                R_R(q,:) = PN_temp_R(:,3,q);

            end
            
            PN_temp_end = R_R;
            
            c = ((k-2)*(n_pol+1)+1):((k-1)*(n_pol+1));
            
            R(a,c) = R_R;
            
        else
            
            down_up_L = down_up_ratio(k);
            
            down_up_R = down_up_ratio(k-1);
        
            R_L = zeros(n_deriv+1,n_pol+1);

            R_R = zeros(n_deriv+1,n_pol+1);

            PN_temp_L = Legendre_polynomial(n_pol,n_deriv,-1:1,down_up_L);
            
            PN_temp_R = Legendre_polynomial(n_pol,n_deriv,-1:1,down_up_R);

            for q = 1:(n_deriv+1)

                R_L(q,:) = PN_temp_L(:,1,q);

                R_R(q,:) = PN_temp_R(:,3,q);

            end
            
            b = ((k-1)*(n_pol+1)+1):(k*(n_pol+1));
            
            c = ((k-2)*(n_pol+1)+1):((k-1)*(n_pol+1));
            
            R(a,b) = R_L;
            
            R(a,c) = -R_R;
            
        end
        
        
    end
    
    
    % Find restriction vector V_r:
    
    V_r = zeros((n_deriv+1)*(nr_wb+1),1);
    
    V_r(1) = PN_temp_start(1,:)*beta_local(1:(n_pol+1));
    V_r(2) = PN_temp_start(2,:)*beta_local(1:(n_pol+1));
    V_r(3) = PN_temp_start(3,:)*beta_local(1:(n_pol+1));
    V_r(4) = PN_temp_start(4,:)*beta_local(1:(n_pol+1));
    V_r(5) = PN_temp_start(5,:)*beta_local(1:(n_pol+1));
    V_r((n_deriv+1)*nr_wb+1) = PN_temp_end(1,:)*beta_local(((nr_wb-1)*(n_pol+1)+1):end);
    V_r((n_deriv+1)*nr_wb+2) = PN_temp_end(2,:)*beta_local(((nr_wb-1)*(n_pol+1)+1):end);
    V_r((n_deriv+1)*nr_wb+3) = PN_temp_end(3,:)*beta_local(((nr_wb-1)*(n_pol+1)+1):end);
    V_r((n_deriv+1)*nr_wb+4) = PN_temp_end(4,:)*beta_local(((nr_wb-1)*(n_pol+1)+1):end);
    V_r((n_deriv+1)*nr_wb+5) = PN_temp_end(5,:)*beta_local(((nr_wb-1)*(n_pol+1)+1):end);

    
    % Calculate restricted polynomial fit:

    
    XTX_inv = inv(X'*X);
    
    beta_star = beta_local - XTX_inv*R'*((R*XTX_inv*R')\(R*beta_local-V_r));
        
    
    beta_element = zeros(n_pol+1,nr_wb);
    
    for l = 1:nr_wb
        
        beta_element(:,l) = beta_star(((l-1)*(n_pol+1)+1):(l*(n_pol+1)));
        
    end
    
    a_fit = beta_element;
        
%      % Compute goodness of fit:
%         
%      
%      var_Y = var((Y(wb_loc(1,1):wb_loc(nr_wb,2))-X*beta_star))
%      
%      goodness = sum((((Y(wb_loc(1,1):wb_loc(nr_wb,2))-X*beta_star)).^2)./sqrt(var_Y))
% 
% %     goodness = kstest2(Y(wb_loc(1,1):wb_loc(nr_wb,2)),X*beta_star)
% 
%      nr_wb*(n_pol+1)-Y_length
%      
%      pause
     
     
    
end
