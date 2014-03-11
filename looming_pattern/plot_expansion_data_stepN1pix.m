
clear
clc
close all

% camera frame time
fps_cam = 7500;
frame_end_cam = 5588;
trigger_frame_cam = 2795;
frames_cam = [1:frame_end_cam]';
t_cam = (frames_cam-trigger_frame_cam)/fps_cam;

% load('cmap_grey2darkred.mat')

%% stepwise slow
load('expansion_steps_Rlimited_exactpixval.mat')

clear t_new teta_new
dt = .0001
t = [exp_steps(1,4)-exp_steps(2,4);exp_steps(:,4)];
teta = [0;exp_steps(:,3)];

for i=2:length(t);
    j=2*i-1;
    t_new(j-1,1)=t(i);
    t_new(j,1)=t(i)+dt;
    teta_new(j:j+1,1)=teta(i);
end
% t_new(1,1)=t_new(1,1)-dt;
% t_new(j+1,1)=t(end)+2*dt;

% trigger delay (see labbook 20131203)
t_new = t_new + 0.0074;
t_new(1,1)=t_cam(1);
t_new(j+1,1)=t_cam(end);

teta_cam = interp1(t_new,teta_new,t_cam);

figure(1)
plot(t_new,teta_new,'.-r')
hold on
plot(t_cam,teta_cam,'.-b')


% figure(2)
% Y = polar(deg2rad(teta_new(1:14)/2),20*t_new(1:14)+10,'r-o');

%% stepwise fast
load('expansion_steps_DOUBLESPEED_looped_exactpixval.mat')

clear t_new teta_new
t = [exp_steps(1,4)-exp_steps(2,4);exp_steps(:,4)];
teta = [0;exp_steps(:,3)];

for i=2:length(t);
    j=2*i-1;
    t_new(j-1:j,1)=t(i);
    teta_new(j:j+1,1)=teta(i);
end
t_new(j+1,1)=t(end);

% trigger delay (see labbook 20131203)
t_new = t_new + 0.0074;

figure(1)
hold on
plot(t_new,teta_new,'g.-')

% figure(2)
% hold on
% Y = polar(deg2rad(teta_new(1:14)/2),20*t_new(1:14)+10,'g-o');

%% continuous fast 165deg
freq = 565;
teta_max = 165;
pixmax = 89;

pix = [[1:2:pixmax] pixmax*ones(1,300)]';
pix(1)=0;

dt = 1/freq;
pixperdeg = 24*8/360;
degperpix = 1/pixperdeg;
teta = pix*degperpix;

clear t
for i = 1:length(teta)
    t(i,1)=(i-1)*dt
end

clear t_new teta_new
for i=2:length(t);
    j=2*i-1;
    t_new(j-1:j,1)=t(i);
    teta_new(j:j+1,1)=teta(i);
end
t_new(j+1,1)=t(end);

% trigger delay (see labbook 20131203)
t_new = t_new + 0.0020;

figure(1)
hold on
plot(t_new,teta_new,'b.-')

% figure(2)
% hold on
% Y = polar(deg2rad(teta_new(1:14)/2),20*t_new(1:14)+10,'g-o');

%% continuous mediumspeed 165deg
freq = 340;
teta_max = 165;
pixmax = 89;

pix = [[1:2:pixmax] pixmax*ones(1,300)]';
pix(1)=0;

dt = 1/freq;
pixperdeg = 24*8/360;
degperpix = 1/pixperdeg;
teta = pix*degperpix;

clear t
for i = 1:length(teta)
    t(i,1)=(i-1)*dt
end

clear t_new teta_new
for i=2:length(t);
    j=2*i-1;
    t_new(j-1:j,1)=t(i);
    teta_new(j:j+1,1)=teta(i);
end
t_new(j+1,1)=t(end);

% trigger delay (see labbook 20131203)
t_new = t_new + 0.0041;

figure(1)
hold on
plot(t_new,teta_new,'k.-')

% figure(2)
% hold on
% Y = polar(deg2rad(teta_new(1:14)/2),20*t_new(1:14)+10,'g-o');



%% continuous fast 64deg
freq = 565;
teta_max = 64;
pixmax = 35;

pix = [[1:2:pixmax] pixmax*ones(1,300)]';
pix(1)=0;

dt = 1/freq;
pixperdeg = 24*8/360;
degperpix = 1/pixperdeg;
teta = pix*degperpix;

clear t
for i = 1:length(teta)
    t(i,1)=(i-1)*dt
end

clear t_new teta_new
for i=2:length(t);
    j=2*i-1;
    t_new(j-1:j,1)=t(i);
    teta_new(j:j+1,1)=teta(i);
end
t_new(j+1,1)=t(end);

% trigger delay (see labbook 20131203)
t_new = t_new + 0.0020;

figure(1)
hold on
plot(t_new,teta_new,'c.-')

% figure(2)
% hold on
% Y = polar(deg2rad(teta_new(1:14)/2),20*t_new(1:14)+10,'g-o');



figure(1)
xlabel('time')
ylabel('optical angle')
set(gca,'xlim',[-0.02 .2],'ylim',[0 180])
set(gca,'XTick',[0:.05:.2])
set(gca,'YTick',[0:45:180])
saveas(gca,'expansion_angle_vs_time.fig')
saveas(gca,'expansion_angle_vs_time.jpg')

