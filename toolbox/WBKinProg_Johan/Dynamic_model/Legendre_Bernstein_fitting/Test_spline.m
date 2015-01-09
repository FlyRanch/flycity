function Test_spline(Y, n_pol )


    % Test b-spline fitting:
    
    t = 1:length(Y)
    
    % Splines
    pp = splinefit(t,Y,16,5)

    pp.breaks
    
    pp.coefs
    
   
    
    % Plot
    figure()
    y1 = ppval(pp,t);
    plot(t,radtodeg(y1),'r')
    hold on 
    plot(t,radtodeg(Y),'b')
    hold off

    figure()
    y4 = ppval(pp,1:0.01:length(Y));
    plot(1:0.01:length(Y),radtodeg(y4),'r')

    qq = ppdiff(pp,1)
    
    rr = ppdiff(pp,2)
    
    % Plot
    figure()
    y2 = ppval(qq,1:0.01:length(Y));
    plot(1:0.01:length(Y),y2,'r')

    
    % Plot
    figure()
    y3 = ppval(rr,1:0.01:length(Y));
    plot(1:0.01:length(Y),y3,'r')
    
    
end

