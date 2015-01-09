function [ a_fit_1, a_fit_2, ratio_1_2 ] = local_fit(Y,n_pol)

        n_deriv = 5;
        
        
        
        % Make local fit
        
        m = length(Y);
        
        m1 = ceil(length(Y)/2);
        
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
        
        ratio_1_2 = m1/m2;

        
        beta_local = [a_fit_local_1; a_fit_local_2];
% %         
% %         a_fit_1 = a_fit_local_1;
% %         
% %         a_fit_2 = a_fit_local_2;      
        
        % Create matrix X
        
        PN_t3 = Legendre_polynomial( n_pol, n_deriv, -1:0.01:1 );
    
        X3 = PN_t3(:,:,1)';
        
%         X4 = [X3(1:200,:) zeros(200,n_pol+1); ...
%               0.5*X3(201,:) 0.5*X3(1,:); ...
%               zeros(200,n_pol+1) X3(2:201,:)];
        
        X4 = [X3 zeros(201,n_pol+1); ...
              zeros(201,n_pol+1) X3];


          
        % Create matrix R
        
        PN_t4 = Legendre_polynomial( n_pol, n_deriv, -1:1 );

        PN_L = zeros(n_deriv+1,n_pol+1);

        PN_R = zeros(n_deriv+1,n_pol+1);

        for k = 1:(n_deriv+1)

            PN_L(k,:) = PN_t4(:,1,k)';

            PN_R(k,:) = PN_t4(:,3,k)';

        end


        R = [-PN_R PN_L];
        
        
        XTX_inv = inv(X4'*X4);


        beta_star = beta_local - XTX_inv*R'*((R*XTX_inv*R')\(R*beta_local));
        
        a_fit_1 = beta_star(1:(n_pol+1));
        
        a_fit_2 = beta_star((n_pol+2):(2*(n_pol+1)));
%         
%         a_fit_1 = a_fit_local_1;
%         
%         a_fit_2 = a_fit_local_2;        
        
        
%         figure()
%         plot(1:1:m1,PN_t1(:,:,1)'*a_fit_local_1,'b')
%         hold on
%         plot(1:1:m1,PN_t1(:,:,1)'*a_fit_1,'r')
%         plot(m1:1:m,PN_t2(:,:,1)'*a_fit_local_2,'b')
%         plot(m1:1:m,PN_t2(:,:,1)'*a_fit_2,'r')
%         plot(1:m,Y,'g')
%         hold off
% 
%         figure()
%         plot(1:1:m1,PN_t1(:,:,2)'*a_fit_local_1,'b')
%         hold on
%         plot(1:1:m1,PN_t1(:,:,2)'*a_fit_1,'r')
%         plot(m1:1:m,PN_t2(:,:,2)'*a_fit_local_2,'b')
%         plot(m1:1:m,PN_t2(:,:,2)'*a_fit_2,'r')
%         hold off
%         
%         figure()
%         plot(1:1:m1,PN_t1(:,:,3)'*a_fit_local_1,'b')
%         hold on
%         plot(1:1:m1,PN_t1(:,:,3)'*a_fit_1,'r')
%         plot(m1:1:m,PN_t2(:,:,3)'*a_fit_local_2,'b')
%         plot(m1:1:m,PN_t2(:,:,3)'*a_fit_2,'r')
%         hold off
%         
%         pause
        
        

end

