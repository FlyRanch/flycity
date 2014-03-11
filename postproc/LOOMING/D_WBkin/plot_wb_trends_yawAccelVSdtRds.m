% yaw acceleration VS ddt & dRds (L-R)
x = yaw_dot_dot_mean_wb(:);

figure
subplot(3,1,1)
title('L-R difference')
hold on

subplot(3,1,1)
y = ddt_ds(:);
color_now = [1 0 0];
plot_MAnCI_constwindow
% xlabel('roll acceleration')
ylabel('down dt diff')
axis([-100000 200000 -.001 .001]) 
set(gca,'YTick',-.005:.001:.005,'fontsize',8)
grid on

subplot(3,1,2)
y = ddt_us(:);
color_now = [0 0 1];
plot_MAnCI_constwindow
% xlabel('roll acceleration')
ylabel('up dt diff')
axis([-100000 200000 -.001 .001]) 
set(gca,'YTick',-.005:.001:.005,'fontsize',8)
grid on

subplot(3,1,3)
Rds_L = dt_ds_L ./ dt_wb_L;
Rds_R = dt_ds_R ./ dt_wb_R;
dRds = Rds_L - Rds_R;
y = dRds(:);
color_now = [0 0 1];
plot_MAnCI_constwindow
% xlabel('roll acceleration')
ylabel('downstroke ratio diff')
axis([-100000 200000 -.15 .15]) 
set(gca,'YTick',-.15:.05:.15,'fontsize',8)
grid on

