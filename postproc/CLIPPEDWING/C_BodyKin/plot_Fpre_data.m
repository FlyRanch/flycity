subplot(3,2,1)
% polar(0, maxHistogramValue,'-k')
% hold on
X = rose(deg2rad(stim_angle_F_pre),n_bins_rose);
for i=1:length(X)
    set(X(i),'color','k','linewidth',1)
end
x = get(X,'Xdata');
y = get(X,'Ydata');
g=patch(x,y,plot_color);
alpha(.5)
title('F relative to stimulus','fontsize',10) 
hold on

subplot(3,2,3)
polar(0, maxHistogramValue,'-k')
hold on
X = rose(deg2rad(Fsp_pitch_pre),n_bins_rose);
for i=1:length(X)
    set(X(i),'color','k','linewidth',1)
end
x = get(X,'Xdata');
y = get(X,'Ydata');
g=patch(x,y,plot_color);
alpha(.5)
title('F pitch relative to strokeplane','fontsize',10) 

subplot(3,2,5)
polar(0, maxHistogramValue,'-k')
hold on
X = rose(deg2rad(Fsp_roll_pre),n_bins_rose);
for i=1:length(X)
    set(X(i),'color','k','linewidth',1)
end
x = get(X,'Xdata');
y = get(X,'Ydata');
g=patch(x,y,plot_color);
alpha(.5)
title('F roll relative to strokeplane','fontsize',10) 

subplot(3,2,2)
% polar(0, maxHistogramValue,'-k')
% hold on
X = rose(deg2rad(slip_pre),n_bins_rose);
for i=1:length(X)
    set(X(i),'color','k','linewidth',1)
end
x = get(X,'Xdata');
y = get(X,'Ydata');
g=patch(x,y,plot_color);
alpha(.5)
title('slip','fontsize',10) 
hold on

subplot(3,2,4)
% polar(0, maxHistogramValue,'-k')
% hold on
X = rose(deg2rad(pitch_pre),n_bins_rose);
% for i=1:length(X)
%     set(X(i),'color','k','linewidth',1)
% end
x = get(X,'Xdata');
y = get(X,'Ydata');
g=patch(x,y,plot_color);
alpha(.5)
title('pitch','fontsize',10) 
hold on

subplot(3,2,6)
% polar(0, maxHistogramValue,'-k')
% hold on
X = rose(deg2rad(roll_pre),n_bins_rose);
% for i=1:length(X)
%     set(X(i),'color','k','linewidth',1)
% end
x = get(X,'Xdata');
y = get(X,'Ydata');
g=patch(x,y,plot_color);
alpha(.5)
title('roll','fontsize',10) 
hold on
