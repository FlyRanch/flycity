function [a_fit] = Piecewise_polynomial_fit4(Y,n_pol,wb_loc)


    nr_wb = length(wb_loc);
    
    Y_length = wb_loc(nr_wb,2)-wb_loc(1,1)+1;
    
    n_deriv = 1;

    
    % Create vector Y
    
    Y_wb = Y(wb_loc(1,1):wb_loc(nr_wb,2));    
    
    
    % Create matrix X1 and R
    
    X1 = zeros(Y_length,(n_pol+1)*(nr_wb+1));
    
    R = zeros((nr_wb+1)*(n_deriv+1),(n_pol+1)*(nr_wb+1));
        
%     for i = 1:nr_wb
%         
%         a = (wb_loc(i,1)-wb_loc(1,1)+1):(wb_loc(i,2)-wb_loc(1,1)+1);
%         b = ((i-1)*(n_pol+1)+1):(i*(n_pol+1));
%         c = (i*(n_pol+1)+1):((i+1)*(n_pol+1));
%         
%         if i == 1
%             
%             m1 = wb_loc(i,2)-wb_loc(i,1)+1;
%            
%             t = [ -1:(1/(m1-1)):1 ];             
%             
%             PN = Legendre_polynomial( n_pol, n_deriv, t );
%             
%             X_i_R = 0.5*PN(:,m1:end,1)';
%             
%             X_ip1_L = 0.5*PN(:,1:m1,1)';
%             
%             X_i_R(1,:) = 1.5*X_i_R(1,:);
%             
%             X_i_R(end,:) = 0.5*X_i_R(end,:);
%             
%             X_ip1_L(1,:) = 0.5*X_ip1_L(1,:);
%             
%             X_ip1_L(end,:) = 0.5*X_ip1_L(end,:);
%                     
%             
%             
%             X1(a,b) = X_i_R;
%             
%             X1(a,c) = X_ip1_L;
%             
%         elseif i == nr_wb
%             
%             m1 = wb_loc(i,2)-wb_loc(i,1)+1;
%            
%             t = [ -1:(1/(m1-1)):1 ];             
%             
%             PN = Legendre_polynomial( n_pol, n_deriv, t );
%             
%             X_i_R = 0.5*PN(:,m1:end,1)';
%             
%             X_ip1_L = 0.5*PN(:,1:m1,1)';
%             
%             X_i_R(1,:) = 0.5*X_i_R(1,:);
%             
%             X_i_R(end,:) = 0.5*X_i_R(end,:);
%             
%             X_ip1_L(1,:) = 0.5*X_ip1_L(1,:);
%             
%             X_ip1_L(end,:) = 1.5*X_ip1_L(end,:);
%             
%             
%             
%             X1(a,b) = X1(a,b)+X_i_R;
%             
%             X1(a,c) = X_ip1_L;
%             
%         else
%             
%             m1 = wb_loc(i,2)-wb_loc(i,1)+1;
%            
%             t = [ -1:(1/(m1-1)):1 ];             
%             
%             PN = Legendre_polynomial( n_pol, n_deriv, t );
%             
%             X_i_R = 0.5*PN(:,m1:end,1)';
%             
%             X_ip1_L = 0.5*PN(:,1:m1,1)';
%             
%             X_i_R(1,:) = 0.5*X_i_R(1,:);
%             
%             X_i_R(end,:) = 0.5*X_i_R(end,:);
%             
%             X_ip1_L(1,:) = 0.5*X_ip1_L(1,:);
%             
%             X_ip1_L(end,:) = 0.5*X_ip1_L(end,:);
%             
% 
%                       
%             X1(a,b) = X1(a,b)+X_i_R;
%                        
%             X1(a,c) = X_ip1_L;
%             
%             
%         end
%         
%     end

    for i = 1:nr_wb
        
        a = (wb_loc(i,1)-wb_loc(1,1)+1):(wb_loc(i,2)-wb_loc(1,1)+1);
        b = ((i-1)*(n_pol+1)+1):(i*(n_pol+1));
        c = (i*(n_pol+1)+1):((i+1)*(n_pol+1));
        
        if i == 1
            
            m1 = wb_loc(i,2)-wb_loc(i,1)+1;
           
            t = [ -1:(1/(m1-1)):1 ];             
            
            PN = Legendre_polynomial( n_pol, n_deriv, t );
            
            X_i_R = PN(:,m1:end,1)';
            
            X_ip1_L = -PN(:,1:m1,1)';
            
                    
            
            
            X1(a,b) = X_i_R;
            
            X1(a,c) = X_ip1_L;
            
        elseif i == nr_wb
            
            m1 = wb_loc(i,2)-wb_loc(i,1)+1;
           
            t = [ -1:(1/(m1-1)):1 ];             
            
            PN = Legendre_polynomial( n_pol, n_deriv, t );
            
            X_i_R = PN(:,m1:end,1)';
            
            X_ip1_L = -PN(:,1:m1,1)';
            
            
            
            X1(a,b) = X1(a,b)+X_i_R;
            
            X1(a,c) = X_ip1_L;
            
        else
            
            m1 = wb_loc(i,2)-wb_loc(i,1)+1;
           
            t = [ -1:(1/(m1-1)):1 ];             
            
            PN = Legendre_polynomial( n_pol, n_deriv, t );
            
            X_i_R = PN(:,m1:end,1)';
            
            X_ip1_L = -PN(:,1:m1,1)';
            
            

                      
            X1(a,b) = X1(a,b)+X_i_R;
                       
            X1(a,c) = X_ip1_L;
            
            
        end
        
    end
    
    X1
    
    PN_t = Legendre_polynomial( n_pol, n_deriv, -1:1:1 );
    
    R_L = zeros(n_deriv+1,n_pol+1);
    
    R_M = zeros(n_deriv+1,n_pol+1);
    
    R_R = zeros(n_deriv+1,n_pol+1);
    
    for j = 1:(n_deriv+1)
        
        R_L(j,:) = 0.5*PN_t(:,3,j)';
        
        R_M(j,:) = 0.5*PN_t(:,2,j)';
        
        R_R(j,:) = 0.5*PN_t(:,1,j)';
        
    end
    
    for k = 1:(nr_wb+1)
        
        if k == 1
            
            b = ((k-1)*(n_deriv+1)+1):(k*(n_deriv+1));            
            c = (k*(n_deriv+1)+1):((k+1)*(n_deriv+1));
            d = ((k-1)*(n_pol+1)+1):(k*(n_pol+1));
            
            R(b,d) = 1.5*R_M;
            R(c,d) = 0.5*R_R;
            
        elseif k == nr_wb+1
            
            a = ((k-2)*(n_deriv+1)+1):((k-1)*(n_deriv+1));
            b = ((k-1)*(n_deriv+1)+1):(k*(n_deriv+1));
            d = ((k-1)*(n_pol+1)+1):(k*(n_pol+1));
            
            R(a,d) = 0.5*R_L;
            R(b,d) = 1.5*R_M;
            
        else
            
            a = ((k-2)*(n_deriv+1)+1):((k-1)*(n_deriv+1));            
            b = ((k-1)*(n_deriv+1)+1):(k*(n_deriv+1));            
            c = (k*(n_deriv+1)+1):((k+1)*(n_deriv+1));
            d = ((k-1)*(n_pol+1)+1):(k*(n_pol+1));
            
            R(a,d) = 0.5*R_L;
            R(b,d) = R_M;
            R(c,d) = 0.5*R_R;
            
        end

    end

    
    % Create matrix X2

    
    
    
    % Solve restricted least squares fit

    XTX_inv = inv(X1'*X1);
    
    beta_local = XTX_inv*(X1'*Y_wb);

    beta_1 = beta_local - XTX_inv*R'*((R*XTX_inv*R')\(R*beta_local));

    
    figure()
    plot(X1*beta_local,'b')
    hold on
    plot(Y_wb,'r')
    hold off
    
    
    
    beta_2 = zeros((n_pol+1)*nr_wb,1);
    
    
    % Solve 


end

