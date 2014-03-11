function val = WBfit_val(a_fit,n,deriv)

    n_pol = length(a_fit.var1)-1;
    PN = Legendre_polynomial(n_pol,2,-1:0.01:1);

    nr_wb = length(a_fit.wb1_wb2(:,1));
    
%     figure()
%     hold on
%     for k = 1:nr_wb
%         m1 =a_fit.wb1_wb2(k,1);
%         m2 = a_fit.wb1_wb2(k,2);
%         t1 = -1:((m1/m2)/(n-1)):(-1+(m1/m2));
%         t2 = (-1+(m1/m2)):((2-(m1/m2))/(n-1)):1;
%         plot(t1,(PN(:,:,1)'*a_fit.var1(k,:)'),'Color',[0.5 0.5 0.5])
%         plot(t2,(PN(:,:,1)'*a_fit.var2(k,:)'),'Color',[0.5 0.5 0.5])
%     end
    m1 =a_fit.wb1_wb2(1);
    m2 = a_fit.wb1_wb2(2);
    t1 = -1:((m1/m2)/round(n/2-1)):(-1+(m1/m2));
    t2 = (-1+(m1/m2)):((2-(m1/m2))/round(n/2-1)):1;
%     plot(t1,(PN(:,:,1)'*a_fit.var1),'r')
%     plot(t2,(PN(:,:,1)'*a_fit.var2),'r')    
    t1_mod=(t1+1)/2;
    t2_mod=(t2+1)/2;
    
    val1 = PN(:,:,deriv)'*a_fit.var1';
    val2 = PN(:,:,deriv)'*a_fit.var2';
    
    val = [val1 val2];

    %     hold off

