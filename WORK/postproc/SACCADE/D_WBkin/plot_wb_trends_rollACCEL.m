% ROLL acceleration VS DEV L-R
figure
subplot(3,2,1)
title('Downstroke L-R')
hold on
subplot(3,2,2)
title('Upstroke L-R')
hold on

subplot(3,2,1)
x = roll_dot_dot_mean_wb(:);
y = rad2deg(ddev_max_udsPREV(:));
color_now = [1 0 0];
p = polyfit(x(isnan(x)==0),y(isnan(y)==0),1);
xp = [min(x) max(x)];
yp = polyval(p,xp);
plot(x,y,'.','color',color_now)
plot(xp,yp,'color',.5*color_now,'linewidth',3)
% xlabel('roll acceleration')
ylabel('L-R @ max dev')
axis([-300000 150000 -45 45])
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,2)
x = roll_dot_dot_mean_wb(:);
y = rad2deg(ddev_max_dus(:));
color_now = [1 0 0];
p = polyfit(x(isnan(x)==0),y(isnan(y)==0),1);
xp = [min(x) max(x)];
yp = polyval(p,xp);
plot(x,y,'.','color',color_now)
plot(xp,yp,'color',.5*color_now,'linewidth',3)
% xlabel('roll acceleration')
% ylabel('Ddev max')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,3)
x = roll_dot_dot_mean_wb(:);
y = rad2deg(ddev_min_ds(:));
color_now = [0 0 1];
p = polyfit(x(isnan(x)==0),y(isnan(y)==0),1);
xp = [min(x) max(x)];
yp = polyval(p,xp);
plot(x,y,'.','color',color_now)
hold on
plot(xp,yp,'color',.5*color_now,'linewidth',3)
% xlabel('roll acceleration')
ylabel('L-R @ min dev')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,4)
x = roll_dot_dot_mean_wb(:);
y = rad2deg(ddev_min_us(:));
color_now = [0 0 1];
p = polyfit(x(isnan(x)==0),y(isnan(y)==0),1);
xp = [min(x) max(x)];
yp = polyval(p,xp);
plot(x,y,'.','color',color_now)
hold on
plot(xp,yp,'color',.5*color_now,'linewidth',3)
% xlabel('roll acceleration')
% ylabel('Ddev min')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,5)
x = roll_dot_dot_mean_wb(:);
y = rad2deg(dAdev_ds(:));
color_now = [.5 .5 .5];
p = polyfit(x(isnan(x)==0),y(isnan(y)==0),1);
xp = [min(x) max(x)];
yp = polyval(p,xp);
plot(x,y,'.','color',color_now)
hold on
plot(xp,yp,'color',.5*color_now,'linewidth',3)
xlabel('roll acceleration')
ylabel('Adev diff')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,6)
x = roll_dot_dot_mean_wb(:);
y = rad2deg(dAdev_us(:));
color_now = [.5 .5 .5];
p = polyfit(x(isnan(x)==0),y(isnan(y)==0),1);
xp = [min(x) max(x)];
yp = polyval(p,xp);
plot(x,y,'.','color',color_now)
hold on
plot(xp,yp,'color',.5*color_now,'linewidth',3)
xlabel('roll acceleration')
% ylabel('dAdev')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

% ROLL acceleration VS PITCH L-R & STROKE L-R
figure
subplot(3,2,1)
title('wing stroke L-R')
hold on
subplot(3,2,2)
title('wing pitch L-R')
hold on

subplot(3,2,1)
x = roll_dot_dot_mean_wb(:);
y = rad2deg(dstroke_max_wb(:));
color_now = [1 0 0];
p = polyfit(x(isnan(x)==0),y(isnan(y)==0),1);
xp = [min(x) max(x)];
yp = polyval(p,xp);
plot(x,y,'.','color',color_now)
hold on
plot(xp,yp,'color',.5*color_now,'linewidth',3)
% xlabel('roll acceleration')
ylabel('L-R @ maximum')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,3)
x = roll_dot_dot_mean_wb(:);
y = rad2deg(dstroke_min_wb(:));
color_now = [0 0 1];
p = polyfit(x(isnan(x)==0),y(isnan(y)==0),1);
xp = [min(x) max(x)];
yp = polyval(p,xp);
plot(x,y,'.','color',color_now)
hold on
plot(xp,yp,'color',.5*color_now,'linewidth',3)
% xlabel('roll acceleration')
ylabel('L-R @ minimum')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,5)
x = roll_dot_dot_mean_wb(:);
y = rad2deg(dAstroke_wb(:));
color_now = [.5 .5 .5];
p = polyfit(x(isnan(x)==0),y(isnan(y)==0),1);
xp = [min(x) max(x)];
yp = polyval(p,xp);
plot(x,y,'.','color',color_now)
hold on
plot(xp,yp,'color',.5*color_now,'linewidth',3)
xlabel('roll acceleration')
ylabel('A diff')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,2)
x = roll_dot_dot_mean_wb(:);
y = rad2deg(dpitch_max_wb(:));
color_now = [1 0 0];
p = polyfit(x(isnan(x)==0),y(isnan(y)==0),1);
xp = [min(x) max(x)];
yp = polyval(p,xp);
plot(x,y,'.','color',color_now)
hold on
plot(xp,yp,'color',.5*color_now,'linewidth',3)
% xlabel('roll acceleration')
% ylabel('Dpitch max')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,4)
x = roll_dot_dot_mean_wb(:);
y = rad2deg(dpitch_min_wb(:));
color_now = [0 0 1];
p = polyfit(x(isnan(x)==0),y(isnan(y)==0),1);
xp = [min(x) max(x)];
yp = polyval(p,xp);
plot(x,y,'.','color',color_now)
hold on
plot(xp,yp,'color',.5*color_now,'linewidth',3)
% xlabel('roll acceleration')
% ylabel('Dpitch min')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,6)
x = roll_dot_dot_mean_wb(:);
y = rad2deg(dApitch_wb(:));
color_now = [.5 .5 .5];
p = polyfit(x(isnan(x)==0),y(isnan(y)==0),1);
xp = [min(x) max(x)];
yp = polyval(p,xp);
plot(x,y,'.','color',color_now)
hold on
plot(xp,yp,'color',.5*color_now,'linewidth',3)
xlabel('roll acceleration')
% ylabel('dApitch')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

