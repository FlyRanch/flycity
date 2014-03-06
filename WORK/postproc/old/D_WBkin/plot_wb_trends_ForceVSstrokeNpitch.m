% normalized force VS PITCH L&R & STROKE L-R
figure
subplot(3,2,1)
title('wing stroke')
hold on
subplot(3,2,2)
title('wing pitch')
hold on

subplot(3,2,1)
x = F_mean_wb(:);
y = rad2deg(stroke_max_wb_L(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd

x = F_mean_wb(:);
y = rad2deg(stroke_max_wb_R(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd

% xlabel('normalized force')
ylabel('stroke max')
axis([0.5 2.5 45 135]) 
set(gca,'YTick',-90:45:180,'fontsize',8)
grid on

subplot(3,2,3)
x = F_mean_wb(:);
y = rad2deg(stroke_min_wb_L(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd

x = F_mean_wb(:);
y = rad2deg(stroke_min_wb_R(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd

% xlabel('normalized force')
ylabel('stroke min')
axis([0.5 2.5 -105 -15]) 
set(gca,'YTick',-180:45:90,'fontsize',8)
grid on

subplot(3,2,5)
x = F_mean_wb(:);
y = rad2deg(Astroke_wb_L(:));
color_now = [.5 .5 .5];
plot_MAnCI_constwindow_circstd

x = F_mean_wb(:);
y = rad2deg(Astroke_wb_R(:));
color_now = [.25 .25 .25];
plot_MAnCI_constwindow_circstd

xlabel('normalized force')
ylabel('Astroke')
axis([0.5 2.5 90 180]) 
set(gca,'YTick',-90:45:180,'fontsize',8)
grid on

subplot(3,2,2)
x = F_mean_wb(:);
y = rad2deg(pitch_max_wb_L(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd

x = F_mean_wb(:);
y = rad2deg(pitch_max_wb_R(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd

% xlabel('normalized force')
ylabel('pitch max')
axis([0.5 2.5 120 210]) 
set(gca,'YTick',-180:45:180,'fontsize',8)
grid on

subplot(3,2,4)
x = F_mean_wb(:);
y = rad2deg(pitch_min_wb_L(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd

x = F_mean_wb(:);
y = rad2deg(pitch_min_wb_R(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd

% xlabel('normalized force')
ylabel('pitch min')
axis([0.5 2.5 -30 75]) 
set(gca,'YTick',-180:45:180,'fontsize',8)
grid on

subplot(3,2,6)
x = F_mean_wb(:);
y = rad2deg(Apitch_wb_L(:));
color_now = [.5 .5 .5];
plot_MAnCI_constwindow_circstd

x = F_mean_wb(:);
y = rad2deg(Apitch_wb_R(:));
color_now = [.25 .25 .25];
plot_MAnCI_constwindow_circstd

xlabel('normalized force')
ylabel('Apitch')
axis([0.5 2.5 90 180]) 
set(gca,'YTick',-180:45:180,'fontsize',8)
grid on

