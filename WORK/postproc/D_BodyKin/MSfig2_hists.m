% close all
figure
subplot(3,3,1)
% title('strokeplane angles','fontsize',10)
hold on
subplot(3,3,2)
% title('body euler angles','fontsize',10)
hold on
subplot(3,3,3)
hold on
subplot(3,3,4)
hold on
subplot(3,3,5)
hold on
subplot(3,3,6)
hold on
subplot(3,3,7)
hold on
subplot(3,3,8)
hold on
subplot(3,3,9)
hold on


%% escape heading
hist_var = stim_angle_vel_post;
xmin=-180;
xmax=180;
dx=10;

subplot(3,3,1)
hist(hist_var,[xmin+dx/2:dx:xmax-dx/2])
% hist(hist_var,[xmin:dx:xmax])
h = findobj(gca,'Type','patch');
set(h,'FaceColor',grey_color)
xlim([xmin xmax])
set(gca,'XTick',xmin:xmax/2:xmax)

% hold on,plot(hist_var,[1:length(hist_var)]/3,'.')

%% escape dV
hist_var = dV_post;
xmin=-.25;
xmax=.5;
dx=.025;

subplot(3,3,4)
hist(hist_var,[xmin+dx/2:dx:xmax-dx/2])
% hist(hist_var,[xmin:dx:xmax])
h = findobj(gca,'Type','patch');
set(h,'FaceColor',grey_color)
xlim([xmin xmax])
set(gca,'XTick',xmin:xmax/2:xmax)

% hold on,plot(hist_var,[1:length(hist_var)]/3,'.')

%% escape Fhor
hist_var = stim_angle_accel_post;
xmin=-180;
xmax=180;
dx=10;

subplot(3,3,7)
hist(hist_var,[xmin+dx/2:dx:xmax-dx/2])
% hist(hist_var,[xmin:dx:xmax])
h = findobj(gca,'Type','patch');
set(h,'FaceColor',grey_color)
xlim([xmin xmax])
set(gca,'XTick',xmin:xmax/2:xmax)

% hold on,plot(hist_var,[1:length(hist_var)]/3,'.')


%% escape Fhor
hist_var = F_max;
xmin=-.5;
xmax=2;
dx=.075;

subplot(3,3,2)
hist(hist_var,[xmin+dx/2:dx:xmax-dx/2])
% hist(hist_var,[xmin:dx:xmax])
h = findobj(gca,'Type','patch');
set(h,'FaceColor',grey_color)
xlim([xmin xmax])
set(gca,'XTick',xmin:xmax/4:xmax)

% hold on,plot(hist_var,[1:length(hist_var)]/3,'.')









