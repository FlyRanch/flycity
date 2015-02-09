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


%% body angles @-20<t<40ms
t_min = t_start;
t_max = t_stop;
t_window_subset = t_window(find(t_window>t_min & t_window<t_max));

xmin=-45;
xmax=90;
dx=5;


%% droll@Amax
hist_var = droll_plot(find(t_window>t_min & t_window<t_max));

subplot(3,2,1)
hist(hist_var,[xmin-dx/2:dx:xmax+dx/2])
% hist(hist_var,[xmin:dx:xmax])
h = findobj(gca,'Type','patch');
set(h,'FaceColor',grey_color)
xlim([xmin xmax])
set(gca,'XTick',xmin:xmax/4:xmax)

% hold on,plot(hist_var,[1:length(hist_var)]/3,'.')

%% dpitch@Amax
hist_var = dpitch_plot(find(t_window>t_min & t_window<t_max));

subplot(3,2,3)
hist(hist_var,[xmin-dx/2:dx:xmax+dx/2])
% hist(hist_var,[xmin:dx:xmax])
h = findobj(gca,'Type','patch');
set(h,'FaceColor',grey_color)
xlim([xmin xmax])
set(gca,'XTick',xmin:xmax/4:xmax)

% hold on,plot(hist_var,[1:length(hist_var)]/3,'.')

%% yaw@Amax
hist_var = yaw_plot(find(t_window>t_min & t_window<t_max));

subplot(3,2,5)
hist(hist_var,[xmin-dx/2:dx:xmax+dx/2])
% hist(hist_var,[xmin:dx:xmax])
h = findobj(gca,'Type','patch');
set(h,'FaceColor',grey_color)
xlim([xmin xmax])
set(gca,'XTick',xmin:xmax/4:xmax)

% hold on,plot(hist_var,[1:length(hist_var)]/3,'.')

%% angle accel @-20<t<40ms
t_min = t_start;
t_max = t_stop;
t_window_subset = t_window(find(t_window>t_min & t_window<t_max));

xmin=-5;
xmax=5;
dx=.5;

%% rollAccel
hist_var = roll_dot_dot_plot(find(t_window>t_min & t_window<t_max))/freq_steady^2;

subplot(3,2,2)
hist(hist_var,[xmin-dx/2:dx:xmax+dx/2])
% hist(hist_var,[xmin:dx:xmax])
h = findobj(gca,'Type','patch');
set(h,'FaceColor',grey_color)
xlim([xmin xmax])
set(gca,'XTick',xmin:xmax/4:xmax)

% hold on,plot(hist_var,[1:length(hist_var)]/3,'.')

%% pitchAccel
hist_var = pitch_dot_dot_plot(find(t_window>t_min & t_window<t_max))/freq_steady^2;

subplot(3,2,4)
hist(hist_var,[xmin-dx/2:dx:xmax+dx/2])
% hist(hist_var,[xmin:dx:xmax])
h = findobj(gca,'Type','patch');
set(h,'FaceColor',grey_color)
xlim([xmin xmax])
set(gca,'XTick',xmin:xmax/4:xmax)

% hold on,plot(hist_var,[1:length(hist_var)]/3,'.')

%% yawAccel@Amax
hist_var = yaw_dot_dot_plot(find(t_window>t_min & t_window<t_max))/freq_steady^2;

subplot(3,2,6)
hist(hist_var,[xmin-dx/2:dx:xmax+dx/2])
% hist(hist_var,[xmin:dx:xmax])
h = findobj(gca,'Type','patch');
set(h,'FaceColor',grey_color)
xlim([xmin xmax])
set(gca,'XTick',xmin:xmax/4:xmax)

% hold on,plot(hist_var,[1:length(hist_var)]/3,'.')

