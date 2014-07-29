% roll acceleration VS PITCH L-R & STROKE L-R
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
plot_MAnCI_constwindow_circstd
% xlabel('roll acceleration')
ylabel('stroke max diff')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,3)
x = roll_dot_dot_mean_wb(:);
y = rad2deg(dstroke_min_wb(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd
% xlabel('roll acceleration')
ylabel('stroke min diff')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,5)
x = roll_dot_dot_mean_wb(:);
y = rad2deg(dAstroke_wb(:));
color_now = [.5 .5 .5];
plot_MAnCI_constwindow_circstd
xlabel('roll acceleration')
ylabel('dAstroke')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,2)
x = roll_dot_dot_mean_wb(:);
y = rad2deg(dpitch_mid_ds(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd
% xlabel('roll acceleration')
ylabel('pitch diff @ mid down')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,4)
x = roll_dot_dot_mean_wb(:);
y = rad2deg(dpitch_mid_us(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd
% xlabel('roll acceleration')
ylabel('pitch diff @ mid up')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,6)
x = roll_dot_dot_mean_wb(:);
y = rad2deg(dApitch_mid(:));
color_now = [.5 .5 .5];
plot_MAnCI_constwindow_circstd
xlabel('roll acceleration')
ylabel('dAmidpitch')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

