% roll acceleration VS PITCH L-R & STROKE L-R
figure
subplot(3,1,1)
title('L-R difference')
hold on

subplot(3,1,1)
x = roll_dot_dot_mean_wb(:);
y = ddt_ds(:);
color_now = [1 0 0];
plot_MAnCI_constwindow
% xlabel('roll acceleration')
ylabel('downstroke dt diff L-R')
axis([-300000 150000 -.001 .001]) 
set(gca,'YTick',-.005:.001:.005,'fontsize',8)
grid on

subplot(3,1,2)
x = roll_dot_dot_mean_wb(:);
y = ddt_us(:);
color_now = [0 0 1];
plot_MAnCI_constwindow
% xlabel('roll acceleration')
ylabel('upstroke dt')
axis([-300000 150000 -.001 .001]) 
set(gca,'YTick',-.005:.001:.005,'fontsize',8)
grid on

subplot(3,1,3)
Rds_L = dt_ds_L ./ dt_wb_L;
Rds_R = dt_ds_R ./ dt_wb_R;
dRds = Rds_L - Rds_R;
x = roll_dot_dot_mean_wb(:);
y = dRds(:);
color_now = [0 0 1];
plot_MAnCI_constwindow
% xlabel('roll acceleration')
ylabel('downstroke ratio')
axis([-300000 150000 -.15 .15]) 
set(gca,'YTick',-.15:.05:.15,'fontsize',8)
grid on

