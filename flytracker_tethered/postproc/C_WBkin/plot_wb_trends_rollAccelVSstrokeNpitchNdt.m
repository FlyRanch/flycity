% roll acceleration VS PITCH L-R & STROKE L-R
figure
subplot(3,3,1)
title('wing stroke L-R')
hold on
subplot(3,3,2)
title('wing pitch L-R')
hold on
subplot(3,3,3)
title('wingstroke timing L-R')
hold on


subplot(3,3,1)
x = roll_dot_dot_mean_wb(:);
y = rad2deg(dstroke_max_wb(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd
% xlabel('roll acceleration')
ylabel('stroke max diff')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,3,4)
x = roll_dot_dot_mean_wb(:);
y = rad2deg(dstroke_min_wb(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd
% xlabel('roll acceleration')
ylabel('stroke min diff')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,3,7)
x = roll_dot_dot_mean_wb(:);
y = rad2deg(dAstroke_wb(:));
color_now = [.5 .5 .5];
plot_MAnCI_constwindow_circstd
xlabel('roll acceleration')
ylabel('dAstroke')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,3,2)
x = roll_dot_dot_mean_wb(:);
y = rad2deg(dpitch_max_wb(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd
% xlabel('roll acceleration')
ylabel('max pitch diff')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,3,5)
x = roll_dot_dot_mean_wb(:);
y = rad2deg(dpitch_min_wb(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd
% xlabel('roll acceleration')
ylabel('min pitch diff')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,3,8)
x = roll_dot_dot_mean_wb(:);
y = rad2deg(dApitch_wb(:));
color_now = [.5 .5 .5];
plot_MAnCI_constwindow_circstd
xlabel('roll acceleration')
ylabel('dApitch')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

% roll acceleration VS ddt & dRds L-R
x = roll_dot_dot_mean_wb(:);

subplot(3,3,3)
y = ddt_ds(:);
color_now = [1 0 0];
plot_MAnCI_constwindow
% xlabel('roll acceleration')
ylabel('down dt diff')
axis([-300000 150000 -.001 .001]) 
set(gca,'YTick',-.005:.001:.005,'fontsize',8)
grid on

subplot(3,3,6)
y = ddt_us(:);
color_now = [0 0 1];
plot_MAnCI_constwindow
% xlabel('roll acceleration')
ylabel('up dt diff')
axis([-300000 150000 -.001 .001]) 
set(gca,'YTick',-.005:.001:.005,'fontsize',8)
grid on

subplot(3,3,9)
Rds_L = dt_ds_L ./ dt_wb_L;
Rds_R = dt_ds_R ./ dt_wb_R;
dRds = Rds_L - Rds_R;
y = dRds(:);
color_now = [0.5 0.5 0.5];
plot_MAnCI_constwindow
xlabel('roll acceleration')
ylabel('downstroke ratio diff')
axis([-300000 150000 -.15 .15]) 
set(gca,'YTick',-.15:.05:.15,'fontsize',8)
grid on
