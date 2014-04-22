X = x(isnan(x)==0);
Y = y(isnan(x)==0);

% % major axis
% [m,b,r,sm,sb]=lsqfitma(X,Y);
% % [m,b,r,sm,sb]=lsqfitgm(X,Y);
% % [m,b,r,sm,sb]=lsqbisec(X,Y);
% MA = [m,b];

% 1st order polyfit
MA = polyfit(X,Y,1);

xp = [min(X) max(X)];
yp = polyval(MA,xp);

% 95% conf int
XY = [X Y];
XY = sortrows(XY);
[n,x_mid] = hist(X,N);
dx = x_mid(2)-x_mid(1);
x_max = min(X);

clear X_ci Y_high Y_low
for n = 1:N
    x_min = x_max;
    x_max = x_max + dx;
    
    X_sub = X(X>x_min & X<x_max);
    Y_sub = Y(X>x_min & X<x_max);
    
    if isempty(X_sub)==0
        Y_MA = polyval(MA,X_sub);
%         dY = Y_sub - Y_MA;
%         dYci = 1.96*circ_std_deg_nonan(dY)/sqrt(length(X_sub));
        dYci = 1.96*nanstd(Y_sub)/sqrt(length(Y_sub));
        if dYci == 0 || isnan(dYci)==1
%         if isnan(dYci)==1
            dYci = 1;
        end
        
        
        if n==1
            X_ci(n,1) = min(X);
            Y_ci = polyval(MA,X_ci(n));
        elseif n==N
            X_ci(n,1) = max(X);
            Y_ci = polyval(MA,X_ci(n));
        else
            X_ci(n,1) = nanmean(X_sub);
            Y_ci = nanmean(Y_MA);
        end
        Y_high(n,1) = Y_ci + dYci;
        Y_low(n,1) = Y_ci - dYci;
    else
        
        if n==1
            X_ci(n,1) = min(X);
        elseif n==N
            X_ci(n,1) = max(X);
        else
            X_ci(n,1) = nanmean([x_min x_max]);
        end
        Y_ci = polyval(MA,X_ci(n));
        Y_high(n,1) = Y_ci+1;
        Y_low(n,1) = Y_ci-1;
    end
        
end


plot(X,Y,'.','color',color_now)
hold on
ciplot(Y_low,Y_high,X_ci,.5*color_now)
alpha(.25)
plot(xp,yp,'color',.5*color_now,'linewidth',3)
