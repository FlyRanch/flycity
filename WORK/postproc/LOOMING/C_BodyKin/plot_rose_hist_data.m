subplot(3,3,1)
polar(0, maxHistogramValue,'-k')
hold on
X = rose(deg2rad(stim_angle_F_now+90),n_bins_rose);
for i=1:length(X)
    set(X(i),'color','k','linewidth',1)
end
x = get(X,'Xdata');
y = get(X,'Ydata');
g=patch(x,y,plot_color);
alpha(.5)
title('F relative to stimulus','fontsize',10) 

subplot(3,3,4)
polar(0, maxHistogramValue,'-k')
hold on
X = rose(deg2rad(Fsp_roll_now+90),n_bins_rose);
for i=1:length(X)
    set(X(i),'color','k','linewidth',1)
end
x = get(X,'Xdata');
y = get(X,'Ydata');
g=patch(x,y,plot_color);
alpha(.5)
title('F-SP roll','fontsize',10) 

subplot(3,3,7)
polar(0, maxHistogramValue,'-k')
hold on
X = rose(deg2rad(Fsp_pitch_now),n_bins_rose);
for i=1:length(X)
    set(X(i),'color','k','linewidth',1)
end
x = get(X,'Xdata');
y = get(X,'Ydata');
g=patch(x,y,plot_color);
alpha(.5)
title('F-SP pitch','fontsize',10) 

subplot(3,3,2)
polar(0, maxHistogramValue,'-k')
hold on
X = rose(deg2rad(roll_now+90),n_bins_rose);
for i=1:length(X)
    set(X(i),'color','k','linewidth',1)
end
x = get(X,'Xdata');
y = get(X,'Ydata');
g=patch(x,y,plot_color);
alpha(.5)
title('roll','fontsize',10) 

subplot(3,3,5)
polar(0, maxHistogramValue,'-k')
hold on
X = rose(deg2rad(pitch_now),n_bins_rose);
% for i=1:length(X)
%     set(X(i),'color','k','linewidth',1)
% end
x = get(X,'Xdata');
y = get(X,'Ydata');
g=patch(x,y,plot_color);
alpha(.5)
title('pitch','fontsize',10) 

subplot(3,3,8)
polar(0, maxHistogramValue,'-k')
hold on
X = rose(deg2rad(yaw_now+90),n_bins_rose);
for i=1:length(X)
    set(X(i),'color','k','linewidth',1)
end
x = get(X,'Xdata');
y = get(X,'Ydata');
g=patch(x,y,plot_color);
alpha(.5)
title('yaw','fontsize',10) 

subplot(3,3,3)
hold on
bin_min = 0;
bin_max = 2.5;
bins = [bin_min+dh/2:dh:bin_max-dh/2];

h1 = hist(F_now,bins);
bar(bins,h1,'FaceColor',plot_color,'EdgeColor','k')

set(gca,'xlim',[bin_min bin_max])
% set(gca,'XTick',[bin_min:dh:bin_max])
alpha(.5)
title('F','fontsize',10) 

subplot(3,3,6)
hold on
bin_min = -2.5;
bin_max = 2.5;
bins = [bin_min+dh/2:dh:bin_max-dh/2];

h1 = hist(Fn_hor_now,bins);
bar(bins,h1,'FaceColor',plot_color,'EdgeColor','k')

set(gca,'xlim',[bin_min bin_max])
% set(gca,'XTick',[bin_min:dh:bin_max])
alpha(.5)
title('Fn','fontsize',10) 

subplot(3,3,9)
hold on
bin_min = -2.5;
bin_max = 2.5;
bins = [bin_min+dh/2:dh:bin_max-dh/2];

h1 = hist(Ft_hor_now,bins);
bar(bins,h1,'FaceColor',plot_color,'EdgeColor','k')

set(gca,'xlim',[bin_min bin_max])
% set(gca,'XTick',[bin_min:dh:bin_max])
alpha(.5)
title('Ft','fontsize',10) 
