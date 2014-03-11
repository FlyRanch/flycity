% plot wingbeat modifiactions for roll accel
colmap = abs(colormap(gray)-1);

%% pitch

t_hist = t_pitch(:);
t_bins = t_pitch(:);

biny_min = -180;
biny_max = 180;

y_min = -45;
y_max = 45;

% wing stroke
var = Dstroke_pitchaccel(:);
subplot(3,2,1)

plot_WBmean_heatmap
hold on
plot_WBmeanCI

% xlabel('time','fontsize',10) 
ylabel('wing stroke','fontsize',10) 
title('pitch accel mod','fontsize',10) 

% wing pitch
var = Dpitch_pitchaccel(:);
subplot(3,2,3)

plot_WBmean_heatmap
hold on
plot_WBmeanCI

% xlabel('time','fontsize',10) 
ylabel('wing pitch','fontsize',10) 
% title('pitch accel mod left')

% dev
var = Ddev_pitchaccel(:);
subplot(3,2,5)

plot_WBmean_heatmap
hold on
plot_WBmeanCI

xlabel('norm time','fontsize',10) 
ylabel('deviation','fontsize',10) 
% title('pitch accel left')




%% Force

t_hist = t_F(:);
t_bins = t_F(:);

biny_min = -180;
biny_max = 180;

y_min = -45;
y_max = 45;

% stroke
var = Dstroke_F(:);
subplot(3,2,2)

plot_WBmean_heatmap
hold on
plot_WBmeanCI

% xlabel('time','fontsize',10) 
% ylabel('wing stroke','fontsize',10) 
title('flight force mod','fontsize',10) 

% pitch
var = Dpitch_F(:);
subplot(3,2,4)

plot_WBmean_heatmap
hold on
plot_WBmeanCI

% xlabel('time','fontsize',10) 
% ylabel('wing pitch','fontsize',10) 
% title('F')

% dev
var = Ddev_F(:);
subplot(3,2,6)

plot_WBmean_heatmap
hold on
plot_WBmeanCI

xlabel('norm time','fontsize',10) 
% ylabel('deviation','fontsize',10) 
% title('F')

