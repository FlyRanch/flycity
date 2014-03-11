function [p_fit,max_var_fit] = Piecewise_polynomial_fit2(Y,n_pol,wb_loc)


    % Piecewise polynomial fit to wingbeat data:
 
    nr_wb = length(wb_loc);
    
    o_pol = 15;
    
    o_deriv = 0;
    
    Y_m = Y(wb_loc(1,1):wb_loc(nr_wb,2));
    
    
    % Construct the X matrix:
    
    X = zeros(length(Y_m),(nr_wb+1)*(o_pol+1));
    
    PN_local = nan(100,o_pol+1,nr_wb+1);
    
%     for i = 1:nr_wb
%         
% %         Y_temp = Y(wb_loc(i,1):wb_loc(i,2));
% 
%         Y_temp = Y(wb_loc(i,1):wb_loc(i,2));
%         
%         n = (wb_loc(i,1)-wb_loc(1,1)+1):(wb_loc(i,2)-wb_loc(1,1)+1);
%         
%         t = -1:(2/(length(Y_temp)-1)):1;
%         
%         k = ((i-1)*(o_pol+1)+1):(i*(o_pol+1));
%         
%         PN_t = Legendre_polynomial(o_pol,o_deriv,t);
%         
%         PN_local(1:length(Y_temp),1:(o_pol+1),i) = PN_t';
%         
%         X(n,k) = PN_t(:,:,1)';
%         
%         if i == 1
%             
%             X(n,k) = PN_t(:,:,1)';
%             
%             X(n(end),k) = -X(n(end),k);
%             
%         elseif i == nr_wb
%             
%             X(n,k) = PN_t(:,:,1)';
%             
%             X(n(1),k) = -X(n(1),k);
%             
%         else
%             
%             X(n,k) = PN_t(:,:,1)';
%             
%             X(n(1),k) = -X(n(1),k);
%             
%             X(n(end),k) = -X(n(end),k);
%             
%         end
%         
%         
%     end


    for i = (1:nr_wb+1)
        
        
        if i == 1

            Y_temp1 = Y(wb_loc(i,1):wb_loc(i,2));
            
            m1 = length(Y_temp1);
            
            t1 = -1:(2/(m1-1)):1;
            
            t2 = 1:1:m1;
            
            PN_t = Legendre_polynomial(o_pol,o_deriv,t1);
            
            a_t = (wb_loc(i,1)-wb_loc(1,1)+1):(wb_loc(i,2)-wb_loc(1,1)+1);
            
            b_t = ((i-1)*(o_pol+1)+1):(i*(o_pol+1));
            
            X(a_t,b_t) = PN_t(:,t2,1)';
            
%             X((wb_loc(i,2)-wb_loc(1,1)+1),b_t) = 0.5*X((wb_loc(i,2)-wb_loc(1,1)+1),b_t);
            
            PN_local(1:m1,:,i) = PN_t(:,t2,1)';
            
        elseif i == nr_wb+1

            Y_temp1 = Y(wb_loc(i-1,1):wb_loc(i-1,2));
            
            m1 = length(Y_temp1);
            
            t1 = -1:(2/(m1-1)):1;
            
            t2 = 1:1:m1;
            
            PN_t = Legendre_polynomial(o_pol,o_deriv,t1);
            
            a_t = (wb_loc(i-1,1)-wb_loc(1,1)+1):(wb_loc(i-1,2)-wb_loc(1,1)+1);
            
            b_t = ((i-1)*(o_pol+1)+1):(i*(o_pol+1));
            
            X(a_t,b_t) = PN_t(:,t2,1)';
            
%             X((wb_loc(i-1,1)-wb_loc(1,1)+1),b_t) = 0.5*X((wb_loc(i-1,1)-wb_loc(1,1)+1),b_t);
            
            PN_local(1:m1,:,i) = PN_t(:,t2,1)';
            
        else
            
            Y_temp1 = Y(wb_loc(i-1,1):wb_loc(i,2));
            
            Y_temp2 = Y(wb_loc(i-1,1):wb_loc(i-1,2));
            
            Y_temp3 = Y(wb_loc(i,1):wb_loc(i,2));
            
            m1 = length(Y_temp1);
                       
            m2 = length(Y_temp2);
            
            m3 = length(Y_temp3);
            
            t1 = [ -1:(1/(m2-1)):0 (1/(m3-1)):(1/(m3-1)):1 ];
            
            PN_t = Legendre_polynomial(o_pol,o_deriv,t1);
            
            a_t = (wb_loc(i-1,1)-wb_loc(1,1)+1):(wb_loc(i,2)-wb_loc(1,1)+1);
            
            b_t = ((i-1)*(o_pol+1)+1):((i)*(o_pol+1));
            
            X(a_t,b_t) = PN_t(:,:,1)';
            
%             X((wb_loc(i-1,1)-wb_loc(1,1)+1),b_t) = 0.5*X((wb_loc(i-1,1)-wb_loc(1,1)+1),b_t);
%             
%             X((wb_loc(i,2)-wb_loc(1,1)+1),b_t) = 0.5*X((wb_loc(i,2)-wb_loc(1,1)+1),b_t);
            
            PN_local(1:m1,:,i) = PN_t(:,:,1)';
            
        end

    end
        
    % Least squares fit:
    
    a_fit = (X'*X)\(X'*Y_m);
    
    % Split in local coefficients:
    
    a_local = zeros(o_pol+1,nr_wb+1);
    
    for j = 1:(nr_wb+1)
        
        k = ((j-1)*(o_pol+1)+1):(j*(o_pol+1));
        
        a_local(:,j) = a_fit(k);
        
    end

    
    figure()
    surf(a_local)
    
    % Plot all the local curves:
    
    for l = 1:nr_wb
        
        Y_temp = Y(wb_loc(l,1):wb_loc(l,2));        
        
        PN_end = find(isnan(PN_local(:,1,l))==0, 1, 'last' );
        
        PN_temp = PN_local(1:PN_end,:,l);
        
        t = -1:(2/(length(Y_temp)-1)):1;
        
        figure()
        plot(t,Y_temp,'b')
        hold on
        plot(t,PN_temp*a_local(:,l),'r')
        hold off
        
    end
    

        
    
    
    
    



end

