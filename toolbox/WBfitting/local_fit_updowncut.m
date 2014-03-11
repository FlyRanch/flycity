function [ a_fit_local_1, a_fit_local_2, ratio_1_2 ] = local_fit_updowncut(Y,n_pol,Rds)

        n_deriv = 5;
        
        
        
        % Make local fit
        
        m = length(Y);
        
        m1 = round(Rds*m);
        
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
        
