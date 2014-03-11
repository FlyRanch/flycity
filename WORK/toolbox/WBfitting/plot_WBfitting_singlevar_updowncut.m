function plot_WBfitting_singlevar_updowncut(a_fit,color_down,color_up)

    n_pol = length(a_fit.var1)-1;

    nr_wb = length(a_fit.wb1_wb2(:,1));
    
%     figure()
%     hold on
%     for k = 1:nr_wb
%         m1 =a_fit.wb1_wb2(k,1);
%         m2 = a_fit.wb1_wb2(k,2);
%         t1 = -1:((m1/m2)/200):(-1+(m1/m2));
%         t2 = (-1+(m1/m2)):((2-(m1/m2))/200):1;
%         plot(t1,(PN(:,:,1)'*a_fit.var1(k,:)'),'Color',[0.5 0.5 0.5])
%         plot(t2,(PN(:,:,1)'*a_fit.var2(k,:)'),'Color',[0.5 0.5 0.5])
%     end
    m1 =a_fit.wb1_wb2(1);
    m2 = a_fit.wb1_wb2(2);
    Rds = m1/(m1+m2)
    
    t1 = 0:(Rds/(m1-1)):Rds;
    t2 = Rds:(1-Rds)/(m2-1):1;
    
    PN1 = Legendre_polynomial(n_pol,2,-1:(2/(m1-1)):1);
    PN2 = Legendre_polynomial(n_pol,2,-1:(2/(m2-1)):1);
    
    y1 = (PN1(:,:,1)'*a_fit.var1');
    y2 = (PN2(:,:,1)'*a_fit.var2');
    
    plot(t1,y1,'color',color_down,'linewidth',1)
    hold on
    plot(t2,y2,'color',color_up,'linewidth',1)    
%     hold off

% test BCs
    dt = t1(2)-t1(1);

    ds1 = (y1(2)-y1(1))/2*(m1-1)
    ds2 = (y1(end)-y1(end-1))/2*(m1-1)

    us1 = (y2(2)-y2(1))/2*(m2-1)
    us2 = (y2(end)-y2(end-1))/2*(m2-1)

    dds1 = ds1-us2
    dds2 = ds2-us1
    
    dy1 = (PN1(:,:,2)'*a_fit.var1');
    dy2 = (PN2(:,:,2)'*a_fit.var2');
    
    dy_ds_1 = dy1(1)
    dy_ds_2 = dy1(end)
    
    dy_us_1 = dy2(1)
    dy_us_2 = dy2(end)
    


