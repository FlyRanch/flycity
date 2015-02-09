% calc angle_post
angle_pre_calc = angle_pre;
angle_post_calc = angle_post;

% plot
angle_pre_plot =   angle_pre_calc(angle_pre_calc>angle_pre_min & angle_pre_calc<angle_pre_max);
angle_post_plot = angle_post_calc(angle_pre_calc>angle_pre_min & angle_pre_calc<angle_pre_max);

% average
angle_preNangle_post = [angle_pre_calc angle_post_calc];

% remove nans
angle_preNangle_post_nans = [angle_pre_calc.*angle_post_calc];
angle_preNangle_post(isnan(angle_preNangle_post_nans)==1,:) = [];

angle_preNangle_post = sortrows(angle_preNangle_post);

angle_pre_sort = angle_preNangle_post(:,1);
angle_post_sort = angle_preNangle_post(:,2);

% circ mean & 95% conf int
% N=dn;
X = angle_pre_calc;
Y = angle_post_calc;

XY = [X Y];
XY = sortrows(XY);
X = XY(:,1);
Y = XY(:,2);

N=round(Nbin/2);
dx = (mean_max-mean_min)/(N-1);

x_mean = mean_min - dx;
clear X_mean Y_ci Y_high Y_low
for n = 1:N
    x_mean = x_mean + dx;
    x_min = x_mean - .5*dt;
    x_max = x_mean + .5*dt;
    
    X_mean(n,1) = x_mean;
    
    X_sub = X(X>x_min & X<x_max);
    Y_sub = Y(X>x_min & X<x_max);
    
    if isempty(X_sub)==0
        
        [Y_ci(n,1) Y_high(n,1) Y_low(n,1)] = circ_mean_deg_nonan(Y_sub);
    else
        Y_ci(n,1) = nan;
        Y_high(n,1) = nan;
        Y_low(n,1) = nan;
    end
        
end

angle_pre_mean = X_mean;

angle_post_mean = Y_ci;
angle_post_high = Y_high;
angle_post_low = Y_low;

angle_pre_conf = angle_pre_mean;
% angle_pre_conf(isnan(angle_post_high)==1)=[];
% angle_post_low(isnan(angle_post_high)==1)=[];
% angle_post_high(isnan(angle_post_high)==1)=[];
angle_post_low(isnan(angle_post_high)==1)=-360;
angle_post_high(isnan(angle_post_high)==1)=360;

angle_post_mean_wrap = angle_post_mean;
angle_post_mean_wrap(angle_post_mean_wrap<-90) = angle_post_mean_wrap(angle_post_mean_wrap<-90)+360;

% hold off
% plot(angle_pre,angle_post,'ok','MarkerFaceColor',plotcolor,'markersize',5)
% plot(angle_pre_plot,angle_post_plot,'ok','MarkerFaceColor',plotcolor,'markersize',2)
hold on
ciplot(angle_post_low,angle_post_high,angle_pre_conf,plotcolor)
% alpha(.25)
plot(angle_pre_mean,angle_post_mean_wrap,'-','linewidth',1,'color',plotcolor)
% plot(angle_pre_mean,angle_post_mean,'-k','linewidth',1)
set(gca,'xlim',[angle_pre_min angle_pre_max])
set(gca,'XTick',[angle_pre_min:(angle_pre_max-angle_pre_min)/2:angle_pre_max])








