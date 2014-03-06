function [ a_fit ] = average_fit(Y,n_pol,wb_loc,Weight)


    % Use a polynomial fit to construct an average wingbeat:
    
    nr_wb = length(wb_loc);
    
    n_deriv = 3;
    
    % Determine local polynomial fit---------------------------------------
    
    % Set up local Y matrix:
    
    Y_local = zeros(wb_loc(end,2)-wb_loc(1,1)+nr_wb,1);

    
    for i = 1:nr_wb
        
        a = (wb_loc(i,1)-wb_loc(1,1)+i):(wb_loc(i,2)-wb_loc(1,1)+i);
        
        b = wb_loc(i,1):wb_loc(i,2);
        
        Y_local(a) = Y(b);
        
    end
    
    
    % Determine matrix X and W
    
    X = zeros(length(Y_local),n_pol+1);
    
    W = zeros(length(Y_local),length(Y_local));
    
    for j = 1:nr_wb
        
        m = wb_loc(j,2)-wb_loc(j,1)+1;
        
        xt = -1:(2/(m-1)):1;
        
        PN_t1 = Legendre_polynomial( n_pol, n_deriv, xt );
        
        a = (wb_loc(j,1)-wb_loc(1,1)+j):(wb_loc(j,2)-wb_loc(1,1)+j);
        
        X(a,:) = PN_t1(:,:,1)'; 
        
        W(a,a) = diag(ones(m,1)*Weight(j));
        
    end
    
    
    beta_local =  (X'*W*X)\(X'*W*Y_local);
    
    %----------------------------------------------------------------------
    
    
    % Create matrix R------------------------------------------------------
    
    
    PN_t2 = Legendre_polynomial( n_pol, n_deriv, -1:1 );
    
    PN_L = zeros(4,n_pol+1);
    
    PN_R = zeros(4,n_pol+1);
    
    for k = 1:((n_deriv+1))
        
        PN_L(k,:) = PN_t2(:,1,k)';
    
        PN_R(k,:) = PN_t2(:,3,k)';
        
    end
    
    
    R = zeros((n_deriv+1)*21,21*(n_pol+1));
    
   
    R_L = PN_L;
    
    R_R = PN_R;   
    
    
    for l = 1:21
        
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

    
    beta_local2 = [beta_local; beta_local; beta_local; beta_local; beta_local; beta_local; beta_local; ...
                   beta_local; beta_local; beta_local; beta_local; beta_local; beta_local; beta_local; ...
                   beta_local; beta_local; beta_local; beta_local; beta_local; beta_local; beta_local];
      
    
    % Apply restricted least square fit on the boundaries:
    
    V_r = zeros((n_deriv+1)*21,1);
    
    V_r(1:(n_deriv+1)) = 0.5*(PN_L+PN_R)*beta_local;
    
    V_r((20*(n_deriv+1)+1):(21*(n_deriv+1))) = 0.5*(PN_L+PN_R)*beta_local;
    
    X2 = X(1:(wb_loc(1,2)-wb_loc(1,1)+1),1:(n_pol+1));
    
    zero_mat = zeros(wb_loc(1,2)-wb_loc(1,1)+1,n_pol+1);
    
    X3 = [X2 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat X2 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat X2 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat X2 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat X2 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat X2 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X2 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X2 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X2 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X2 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X2 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X2 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X2 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X2 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X2 zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X2 zero_mat zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X2 zero_mat zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X2 zero_mat zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X2 zero_mat zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X2 zero_mat; ...
          zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat zero_mat X2];
    
    XTX_inv = inv(X3'*X3);

   
    beta_star = beta_local2 - XTX_inv*R'*((R*XTX_inv*R')\(R*beta_local2-V_r));
    
    a_fit = beta_star(((n_pol+1)*10+1):((n_pol+1)*11));
    
    %a_fit = beta_local;
    
    
    
end

