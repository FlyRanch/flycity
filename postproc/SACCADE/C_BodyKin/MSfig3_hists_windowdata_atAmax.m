% close all
figure
subplot(3,2,1)
% title('strokeplane angles','fontsize',10)
hold on
subplot(3,2,2)
% title('body euler angles','fontsize',10)
hold on
subplot(3,2,3)
hold on
subplot(3,2,4)
hold on
subplot(3,2,5)
hold on
subplot(3,2,6)
hold on



%% roll@Amax
hist_var = calc_value(roll_plot,n_Amax);
xmin=-45;
xmax=90;
dx=5;

subplot(3,2,1)
hist(hist_var,[xmin-dx/2:dx:xmax+dx/2])
% hist(hist_var,[xmin:dx:xmax])
h = findobj(gca,'Type','patch');
set(h,'FaceColor',grey_color)
xlim([xmin xmax])
set(gca,'XTick',xmin:xmax/4:xmax)

% hold on,plot(hist_var,[1:length(hist_var)]/3,'.')

%% pitch@Amax
hist_var = calc_value(pitch_plot,n_Amax);
xmin=-45;
xmax=90;
dx=5;

subplot(3,2,3)
hist(hist_var,[xmin-dx/2:dx:xmax+dx/2])
% hist(hist_var,[xmin:dx:xmax])
h = findobj(gca,'Type','patch');
set(h,'FaceColor',grey_color)
xlim([xmin xmax])
set(gca,'XTick',xmin:xmax/4:xmax)

% hold on,plot(hist_var,[1:length(hist_var)]/3,'.')

%% yaw@Amax
hist_var = calc_value(yaw_plot,n_Amax);
xmin=-45;
xmax=90;
dx=5;

subplot(3,2,5)
hist(hist_var,[xmin-dx/2:dx:xmax+dx/2])
% hist(hist_var,[xmin:dx:xmax])
h = findobj(gca,'Type','patch');
set(h,'FaceColor',grey_color)
xlim([xmin xmax])
set(gca,'XTick',xmin:xmax/4:xmax)

% hold on,plot(hist_var,[1:length(hist_var)]/3,'.')

%% angle accel @Amax
t_min = t_start;
t_max = t_stop;
t_window_subset = t_window(find(t_window>t_min & t_window<t_max));

xmin=-5;
xmax=5;
dx=.5;

%% rollAccel
hist_var = calc_value(roll_dot_dot_plot,n_Amax)/freq_steady^2;

subplot(3,2,2)
hist(hist_var,[xmin-dx/2:dx:xmax+dx/2])
% hist(hist_var,[xmin:dx:xmax])
h = findobj(gca,'Type','patch');
set(h,'FaceColor',grey_color)
xlim([xmin xmax])
set(gca,'XTick',xmin:xmax/4:xmax)

% hold on,plot(hist_var,[1:length(hist_var)]/3,'.')

%% pitchAccel
hist_var = calc_value(pitch_dot_dot_plot,n_Amax)/freq_steady^2;

subplot(3,2,4)
hist(hist_var,[xmin-dx/2:dx:xmax+dx/2])
% hist(hist_var,[xmin:dx:xmax])
h = findobj(gca,'Type','patch');
set(h,'FaceColor',grey_color)
xlim([xmin xmax])
set(gca,'XTick',xmin:xmax/4:xmax)

% hold on,plot(hist_var,[1:length(hist_var)]/3,'.')

%% yawAccel
hist_var = calc_value(yaw_dot_dot_plot,n_Amax)/freq_steady^2;

subplot(3,2,6)
hist(hist_var,[xmin-dx/2:dx:xmax+dx/2])
% hist(hist_var,[xmin:dx:xmax])
h = findobj(gca,'Type','patch');
set(h,'FaceColor',grey_color)
xlim([xmin xmax])
set(gca,'XTick',xmin:xmax/4:xmax)

% hold on,plot(hist_var,[1:length(hist_var)]/3,'.')

