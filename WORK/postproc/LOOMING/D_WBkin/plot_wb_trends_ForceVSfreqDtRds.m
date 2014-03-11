% pitch acceleration VS flapfreq dt & Rds
x = F_mean_wb(:);

figure
subplot(2,2,1)
% title('L-R difference')
hold on

subplot(2,2,1)
y = dt_ds_L(:);
color_now = [1 0 0];
plot_MAnCI_constwindow
y = dt_ds_R(:);
color_now = [0 0 1];
plot_MAnCI_constwindow
% xlabel('roll acceleration')
ylabel('downstroke dt')
axis([0.5 2.5 0 .005]) 
set(gca,'XTick',0:.5:3,'fontsize',8)
set(gca,'YTick',-.005:.0025:.005,'fontsize',8)
grid on

subplot(2,2,3)
y = dt_us_L(:);
color_now = [1 0 0];
plot_MAnCI_constwindow
y = dt_us_R(:);
color_now = [0 0 1];
plot_MAnCI_constwindow
% xlabel('roll acceleration')
ylabel('downstroke dt')
axis([0.5 2.5 0 .005]) 
set(gca,'XTick',0:.5:3,'fontsize',8)
set(gca,'YTick',-.005:.0025:.005,'fontsize',8)
grid on


subplot(2,2,2)
y = f_wb_L(:);
color_now = [1 0 0];
plot_MAnCI_constwindow
y = f_wb_R(:);
color_now = [0 0 1];
plot_MAnCI_constwindow
% xlabel('roll acceleration')
ylabel('flap freq')
axis([0.5 2.5 100 300]) 
set(gca,'XTick',0:.5:3,'fontsize',8)
set(gca,'YTick',0:100:500,'fontsize',8)
grid on


subplot(2,2,4)
Rds_L = dt_ds_L ./ dt_wb_L;
Rds_R = dt_ds_R ./ dt_wb_R;
y = Rds_L(:);
color_now = [1 0 0];
plot_MAnCI_constwindow
y = Rds_R(:);
color_now = [0 0 1];
plot_MAnCI_constwindow
% xlabel('roll acceleration')
ylabel('downstroke ratio')
axis([0.5 2.5 .25 .75]) 
set(gca,'XTick',0:.5:3,'fontsize',8)
set(gca,'YTick',0:.25:1,'fontsize',8)
grid on

