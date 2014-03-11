
Ax_first = A_hor_first .* cosd(Adir_first);
Ay_first = A_hor_first .* sind(Adir_first);

Ax_post = A_hor_post .* cosd(Adir_post);
Ay_post = A_hor_post .* sind(Adir_post);

Ax_mean = A_hor_mean .* cosd(Adir_mean);
Ay_mean = A_hor_mean .* sind(Adir_mean);

Ax_max = A_hor_max .* cosd(Adir_Ahormax);
Ay_max = A_hor_max .* sind(Adir_Ahormax);

maxHistogramValue = (max(A_hor_max));
maxHistogramValue = 20;

subplot(2,3,1)
% polar(0, maxHistogramValue,'-r')
hold on
compass(-Ay_first,Ax_first,'r');

subplot(2,3,2)
polar(0, maxHistogramValue,'-r')
hold on
compass(-Ay_max,Ax_max,'r');

subplot(2,3,3)
polar(0, maxHistogramValue,'-r')
hold on
compass(-Ay_mean,Ax_mean,'r');

% figure
subplot(2,3,4)
polar(0, maxHistogramValue,'-r')
hold on
Y = rose(deg2rad(Adir_first+90),36);
% Y = rose(deg2rad(Adir_body_first),36);
for i=1:length(Y)
    set(Y(i),'color','k','linewidth',1)
end
x = get(Y,'Xdata');
y = get(Y,'Ydata');
g=patch(x,y,'r');
alpha(.5)

subplot(2,3,5)
polar(0, maxHistogramValue,'-r')
hold on
X = rose(deg2rad(Adir_Ahormax+90),36);
for i=1:length(X)
    set(X(i),'color','k','linewidth',1)
end
x = get(X,'Xdata');
y = get(X,'Ydata');
g=patch(x,y,'r');
alpha(.5)

subplot(2,3,6)
polar(0, maxHistogramValue,'-r')
hold on
X = rose(deg2rad(Adir_mean+90),36);
for i=1:length(X)
    set(X(i),'color','k','linewidth',1)
end
x = get(X,'Xdata');
y = get(X,'Ydata');
g=patch(x,y,'r');
alpha(.5)


