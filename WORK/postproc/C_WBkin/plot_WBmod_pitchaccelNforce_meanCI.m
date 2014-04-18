% plot wingbeat modifiactions for pitch & F
%% pitch
t_bins = t_pitch_LR(:);
y_min = -45;
y_max = 45;

% wing stroke
var = Dstroke_pitchaccel(:);
subplot(3,2,1)
hold on
plot_WBmeanCI
% ylabel('wing stroke','fontsize',10) 
% title('pitch acceleration','fontsize',10) 

% wing pitch
var = Dpitch_pitchaccel(:);
subplot(3,2,3)
hold on
plot_WBmeanCI
% ylabel('wing pitch','fontsize',10) 

% dev
var = Ddev_pitchaccel(:);
subplot(3,2,5)
hold on
plot_WBmeanCI
% ylabel('wing deviation','fontsize',10) 




%% Force

t_bins = t_F_LR(:);
y_min = -45;
y_max = 45;

% wing stroke
var = Dstroke_F(:);
subplot(3,2,2)
hold on
plot_WBmeanCI
% title('flight force','fontsize',10) 

% wing pitch
var = Dpitch_F(:);
subplot(3,2,4)
hold on
plot_WBmeanCI

% dev
var = Ddev_F(:);
subplot(3,2,6)
hold on
plot_WBmeanCI

