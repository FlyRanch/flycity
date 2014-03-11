X = x(isnan(x)==0);
Y = y(isnan(y)==0);

% major axis
[m,b,r,sm,sb]=lsqfitma(X,Y);
MA = [m,b];
xp = [min(X) max(X)];
yp = polyval(MA,xp);

% 95% conf int
XY = [X Y];
XY = sortrows(XY);
p=0;
clear X_ci Y_high Y_low
for n = 1:dm:(length(X)-dn)
    p=p+1;
    
    X_sub = XY(n:n+dn-1,1);
    Y_sub = XY(n:n+dn-1,2);
    Y_MA = polyval(MA,X_sub);
    dY = Y_sub - Y_MA;
    dYci = 1.96*std(dY)/sqrt(dn);
    
    X_ci(p,1) = nanmean(X_sub);
    Y_high(p,1) = nanmean(Y_MA) + dYci;
    Y_low(p,1) = nanmean(Y_MA) - dYci;
end

plot(X,Y,'.','color',color_now)
hold on
ciplot(Y_low,Y_high,X_ci,.5*color_now)
alpha(.25)
plot(xp,yp,'color',.5*color_now,'linewidth',3)
