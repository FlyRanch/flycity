subplot(2,2,1)
polar(0, maxHistogramValue,'-k')
hold on
Z = rose(deg2rad(heading_post(teta_max==16 & OFF == 0)),36);
for i=1:length(Z)
    set(Z(i),'color','k','linewidth',1)
end
x = get(Z,'Xdata');
y = get(Z,'Ydata');
g=patch(x,y,'r');

Y = rose(deg2rad(heading_post(teta_max==16 & OFF == 1)),36);
for i=1:length(Y)
    set(Y(i),'color','k','linewidth',1)
end 
x = get(Y,'Xdata');
y = get(Y,'Ydata');
g=patch(x,y,'b');
title('max optical angle = 16deg','fontsize',18) 
alpha(.5)

subplot(2,2,2)
polar(0, maxHistogramValue,'-k')
hold on
Z = rose(deg2rad(heading_post(teta_max==32 & OFF == 0)),36);
for i=1:length(Z)
    set(Z(i),'color','k','linewidth',1)
end
x = get(Z,'Xdata');
y = get(Z,'Ydata');
g=patch(x,y,'r');

Y = rose(deg2rad(heading_post(teta_max==32 & OFF == 1)),36);
for i=1:length(Y)
    set(Y(i),'color','k','linewidth',1)
end 
x = get(Y,'Xdata');
y = get(Y,'Ydata');
g=patch(x,y,'b');
title('max optical angle = 32deg','fontsize',18) 
alpha(.5)

subplot(2,2,3)
polar(0, maxHistogramValue,'-k')
hold on
Z = rose(deg2rad(heading_post(teta_max==64 & OFF == 0)),36);
for i=1:length(Z)
    set(Z(i),'color','k','linewidth',1)
end
x = get(Z,'Xdata');
y = get(Z,'Ydata');
g=patch(x,y,'r');

Y = rose(deg2rad(heading_post(teta_max==64 & OFF == 1)),36);
for i=1:length(Y)
    set(Y(i),'color','k','linewidth',1)
end 
x = get(Y,'Xdata');
y = get(Y,'Ydata');
g=patch(x,y,'b');
title('max optical angle = 64deg','fontsize',18) 
alpha(.5)

subplot(2,2,4)
polar(0, maxHistogramValue,'-k')
hold on
Z = rose(deg2rad(heading_post(teta_max==165 & OFF == 0)),36);
for i=1:length(Z)
    set(Z(i),'color','k','linewidth',1)
end
x = get(Z,'Xdata');
y = get(Z,'Ydata');
g=patch(x,y,'r');

Y = rose(deg2rad(heading_post(teta_max==165 & OFF == 1)),36);
for i=1:length(Y)
    set(Y(i),'color','k','linewidth',1)
end 
x = get(Y,'Xdata');
y = get(Y,'Ydata');
g=patch(x,y,'b');
title('max optical angle = 165deg','fontsize',18) 
alpha(.5)
