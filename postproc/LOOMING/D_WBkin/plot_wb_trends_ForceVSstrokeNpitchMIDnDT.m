% normalized force VS PITCH L&R & STROKE L-R
figure
subplot(3,3,1)
title('wing stroke')
hold on
subplot(3,3,2)
title('wing pitch')
hold on
subplot(3,3,3)
title('wingstroke timing L-R')
hold on

% flight force
x = F_mean_wb(:);

subplot(3,3,1)
y = rad2deg(stroke_max_wb_L(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd
MA_L = MA;

y = rad2deg(stroke_max_wb_R(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd
MA = mean([MA_L;MA]);

% xlabel('normalized force')
ylabel('stroke max')
axis([0.5 2.5 45 135]) 
set(gca,'YTick',-90:45:180,'fontsize',8)
grid on
fit_force_strokemax = MA;

subplot(3,3,4)
y = rad2deg(stroke_min_wb_L(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd
MA_L = MA;

y = rad2deg(stroke_min_wb_R(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd
MA = mean([MA_L;MA]);

% xlabel('normalized force')
ylabel('stroke min')
axis([0.5 2.5 -105 -15]) 
set(gca,'YTick',-180:45:90,'fontsize',8)
grid on
fit_force_strokemin = MA;

subplot(3,3,7)
y = rad2deg(Astroke_wb_L(:));
color_now = [.5 .5 .5];
plot_MAnCI_constwindow_circstd
MA_L = MA;

y = rad2deg(Astroke_wb_R(:));
color_now = [.25 .25 .25];
plot_MAnCI_constwindow_circstd
MA = mean([MA_L;MA]);

xlabel('normalized force')
ylabel('Astroke')
axis([0.5 2.5 90 180]) 
set(gca,'YTick',-90:45:180,'fontsize',8)
grid on
fit_force_Astroke = MA;

subplot(3,3,2)
y = rad2deg(pitch_mid_ds_L(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd
MA_L = MA;

y = rad2deg(pitch_mid_ds_R(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd
MA = mean([MA_L;MA]);

% xlabel('normalized force')
ylabel('pitch @ mid ds')
axis([0.5 2.5 90 180]) 
set(gca,'YTick',-180:45:180,'fontsize',8)
grid on
fit_force_pitchmid_ds = MA;

subplot(3,3,5)
y = rad2deg(pitch_mid_us_L(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd
MA_L = MA;

y = rad2deg(pitch_mid_us_R(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd
MA = mean([MA_L;MA]);

% xlabel('normalized force')
ylabel('pitch @ mid us')
axis([0.5 2.5 0 90]) 
set(gca,'YTick',-180:45:180,'fontsize',8)
grid on
fit_force_pitchmid_us = MA;

subplot(3,3,8)
y = rad2deg(Apitch_mid_L(:));
color_now = [.5 .5 .5];
plot_MAnCI_constwindow_circstd
MA_L = MA;

y = rad2deg(Apitch_mid_R(:));
color_now = [.25 .25 .25];
plot_MAnCI_constwindow_circstd
MA = mean([MA_L;MA]);

xlabel('normalized force')
ylabel('Amidpitch')
axis([0.5 2.5 45 135]) 
set(gca,'YTick',-180:45:180,'fontsize',8)
grid on
fit_force_Apitchmid = MA;


% flapfreq dt & Rds
subplot(3,3,3)
y = dt_ds_L(:);
color_now = [0 0 1];
plot_MAnCI_constwindow
MA_L = MA;

y = dt_ds_R(:);
color_now = [1 0 0];
plot_MAnCI_constwindow
MA = mean([MA_L;MA]);

% xlabel('roll acceleration')
ylabel('downstroke dt')
axis([0.5 2.5 0 .005]) 
set(gca,'XTick',0:.5:3,'fontsize',8)
set(gca,'YTick',-.005:.0025:.005,'fontsize',8)
grid on
fit_force_dt_ds = MA;

subplot(3,3,6)
y = dt_us_L(:);
color_now = [0 0 1];
plot_MAnCI_constwindow
MA_L = MA;

y = dt_us_R(:);
color_now = [1 0 0];
plot_MAnCI_constwindow
MA = mean([MA_L;MA]);

% xlabel('roll acceleration')
ylabel('downstroke dt')
axis([0.5 2.5 0 .005]) 
set(gca,'XTick',0:.5:3,'fontsize',8)
set(gca,'YTick',-.005:.0025:.005,'fontsize',8)
grid on
fit_force_dt_us = MA;

subplot(3,3,9)
y = f_wb_L(:);
color_now = [.5 .5 .5];
plot_MAnCI_constwindow
MA_L = MA;

y = f_wb_R(:);
color_now = [0.25 0.25 0.25];
plot_MAnCI_constwindow
MA = mean([MA_L;MA]);

xlabel('normalized force')
ylabel('flap freq')
axis([0.5 2.5 100 300]) 
set(gca,'XTick',0:.5:3,'fontsize',8)
set(gca,'YTick',0:100:500,'fontsize',8)
grid on
fit_force_freq = MA;
