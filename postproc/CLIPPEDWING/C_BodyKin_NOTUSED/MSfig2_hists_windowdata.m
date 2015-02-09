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


%% escape heading
t_min = t_stop - dt;
t_max = t_stop + dt;
t_window_subset = t_window(find(t_window>t_min & t_window<t_max));

hist_var = stim_angle_vel_plot(find(t_window>t_min & t_window<t_max));
% hist_var = stim_angle_vel_post;

xmin=-180;
xmax=180;
dx=10;

subplot(3,2,1)
hist(hist_var,[xmin-dx/2:dx:xmax+dx/2])
% hist(hist_var,[xmin:dx:xmax])
h = findobj(gca,'Type','patch');
set(h,'FaceColor',grey_color)
xlim([xmin xmax])
set(gca,'XTick',xmin:xmax/2:xmax)

% hold on,plot(hist_var,[1:length(hist_var)]/3,'.')

%% escape dV
t_min = t_stop - dt;
t_max = t_stop + dt;
t_window_subset = t_window(find(t_window>t_min & t_window<t_max));

hist_var = dV(find(t_window>t_min & t_window<t_max));
% hist_var = dV_post;

xmin=-.25;
xmax=.5;
dx=.025;

subplot(3,2,3)
hist(hist_var,[xmin-dx/2:dx:xmax+dx/2])
% hist(hist_var,[xmin:dx:xmax])
h = findobj(gca,'Type','patch');
set(h,'FaceColor',grey_color)
xlim([xmin xmax])
set(gca,'XTick',xmin:xmax/2:xmax)

% hold on,plot(hist_var,[1:length(hist_var)]/3,'.')

%% escape Fhor
% t_min = t_stop - dt;
% t_max = t_stop + dt;
% t_window_subset = t_window(find(t_window>t_min & t_window<t_max));
% 
% hist_var = stim_angle_accel_plot(find(t_window>t_min & t_window<t_max));
% hist_var = stim_angle_accel_post;
% 
% xmin=-180;
% xmax=180;
% dx=10;
% 
% subplot(3,2,5)
% hist(hist_var,[xmin-dx/2:dx:xmax+dx/2])
% % hist(hist_var,[xmin:dx:xmax])
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor',grey_color)
% xlim([xmin xmax])
% set(gca,'XTick',xmin:xmax/2:xmax)
% 
% % hold on,plot(hist_var,[1:length(hist_var)]/3,'.')

%% Fhor at 0<t<40ms
t_min = 0;
t_max = t_stop;

t_window_subset = t_window(find(t_window>t_min & t_window<t_max));
hist_var = stim_angle_accel_plot(find(t_window>t_min & t_window<t_max));

xmin=-180;
xmax=180;
dx=10;

subplot(3,2,5)
hist(hist_var,[xmin-dx/2:dx:xmax+dx/2])
% hist(hist_var,[xmin:dx:xmax])
h = findobj(gca,'Type','patch');
set(h,'FaceColor',grey_color)
xlim([xmin xmax])
set(gca,'XTick',xmin:xmax/2:xmax)

% hold on,plot(hist_var,[1:length(hist_var)]/3,'.')


%% F/Mg max
hist_var = F_max;
xmin=-.5;
xmax=2;
dx=.1;

subplot(3,2,2)
hist(hist_var,[xmin-dx/2:dx:xmax+dx/2])
% hist(hist_var,[xmin:dx:xmax])
h = findobj(gca,'Type','patch');
set(h,'FaceColor',grey_color)
xlim([xmin xmax])
set(gca,'XTick',xmin:xmax/4:xmax)

% hold on,plot(hist_var,[1:length(hist_var)]/3,'.')




%% Froll at -20<t<40ms
t_min = t_start;
t_max = t_stop;

t_window_subset = t_window(find(t_window>t_min & t_window<t_max));
hist_var = Fsp_roll_plot(find(t_window>t_min & t_window<t_max));

xmin=-30;
xmax=30;
dx=2;

subplot(3,2,4)
hist(hist_var,[xmin-dx/2:dx:xmax+dx/2])
% hist(hist_var,[xmin:dx:xmax])
h = findobj(gca,'Type','patch');
set(h,'FaceColor',grey_color)
xlim([xmin xmax])
set(gca,'XTick',xmin:xmax/2:xmax)

% hold on,plot(hist_var,[1:length(hist_var)]/3,'.')


%% Fpitch at -20<t<40ms
t_min = t_start;
t_max = t_stop;

t_window_subset = t_window(find(t_window>t_min & t_window<t_max));
hist_var = Fsp_pitch_plot(find(t_window>t_min & t_window<t_max));

xmin=-30;
xmax=30;
dx=2;

subplot(3,2,6)
hist(hist_var,[xmin-dx/2:dx:xmax+dx/2])
% hist(hist_var,[xmin:dx:xmax])
h = findobj(gca,'Type','patch');
set(h,'FaceColor',grey_color)
xlim([xmin xmax])
set(gca,'XTick',xmin:xmax/2:xmax)

% hold on,plot(hist_var,[1:length(hist_var)]/3,'.')






