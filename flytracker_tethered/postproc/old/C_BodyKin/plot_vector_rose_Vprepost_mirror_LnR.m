figure
% maxHistogramValue = ceil(max(sqrt(u_post.^2 + v_post.^2)));
maxHistogramValue = floor(max(sqrt(u_post.^2 + v_post.^2)))-1;

subplot(1,2,1)
polar(0, maxHistogramValue,'-k')
% plot(0,0)
hold on
% plot(0,0,'r')
% legend('pre','post')
X = compass_zeroup(u_post,v_post,'b');
for i=1:length(X)
    set(X(i),'color','b','linewidth',1.5)
end
hold on

Z = compass(u_post_l,v_post_l,'r');
for i=1:length(Z)
    set(Z(i),'color','r','linewidth',1.5)
end

Y = compass(u_pre,v_pre,'k');
for i=1:length(Y)
    set(Y(i),'color','k','linewidth',1.5)
end 

% figure
subplot(1,2,2)

X = rose(deg2rad(Vdir_post+90),36);
for i=1:length(X)
    set(X(i),'color','k','linewidth',1.5)
end
x = get(X,'Xdata');
y = get(X,'Ydata');
g=patch(x,y,'b');
alpha(.5)
hold on

Z = rose(deg2rad(-Vdir_post+90),36);
for i=1:length(Z)
    set(Z(i),'color','k','linewidth',1.5)
end
x = get(Z,'Xdata');
y = get(Z,'Ydata');
g=patch(x,y,'r');
alpha(.5)

% Y = rose(deg2rad([Vdir_pre+90;-Vdir_pre+90]),36);
Y = rose(deg2rad([Vdir_pre+90]),36);
for i=1:length(Y)
    set(Y(i),'color','k','linewidth',1.5)
end 
x = get(Y,'Xdata');
y = get(Y,'Ydata');
g=patch(x,y,[.5 .5 .5]);
alpha(.5)

compass_zeroup(0, 0.01,'-k');
