

% F VS DEV
figure
subplot(3,2,1)
title('Downstroke Deviation')
hold on
subplot(3,2,2)
title('Upstroke Deviation')
hold on

subplot(3,2,1)
plot(F_mean_wb(:),rad2deg(dev_max_udsPREV_L(:)),'.')
hold on
plot(F_mean_wb(:),rad2deg(dev_max_udsPREV_R(:)),'.r')
% xlabel('F')
ylabel('max dev')
axis([.5 2.5 -15 45])            
set(gca,'YTick',-90:15:90,'fontsize',8)

subplot(3,2,2)
plot(F_mean_wb(:),rad2deg(dev_max_dus_L(:)),'.')
hold on
plot(F_mean_wb(:),rad2deg(dev_max_dus_R(:)),'.r')
% xlabel('F')
% ylabel('max dev')
axis([.5 2.5 -15 45])            
set(gca,'YTick',-90:15:90,'fontsize',8)

subplot(3,2,3)
plot(F_mean_wb(:),rad2deg(dev_min_ds_L(:)),'.')
hold on
plot(F_mean_wb(:),rad2deg(dev_min_ds_R(:)),'.r')
% xlabel('F')
ylabel('min dev')
axis([.5 2.5 -45 15])            
set(gca,'YTick',-90:15:90,'fontsize',8)

subplot(3,2,4)
plot(F_mean_wb(:),rad2deg(dev_min_us_L(:)),'.')
hold on
plot(F_mean_wb(:),rad2deg(dev_min_us_R(:)),'.r')
% xlabel('F')
% ylabel('min dev')
axis([.5 2.5 -45 15])            
set(gca,'YTick',-90:15:90,'fontsize',8)

subplot(3,2,5)
plot(F_mean_wb(:),rad2deg(Adev_ds_L(:)),'.')
hold on
plot(F_mean_wb(:),rad2deg(Adev_ds_R(:)),'.r')
xlabel('F')
ylabel('Adev')
axis([.5 2.5 -0 45])            
set(gca,'YTick',-90:15:90,'fontsize',8)

subplot(3,2,6)
plot(F_mean_wb(:),rad2deg(Adev_us_L(:)),'.')
hold on
plot(F_mean_wb(:),rad2deg(Adev_us_R(:)),'.r')
xlabel('F')
% ylabel('Adev')
axis([.5 2.5 -0 45])            
set(gca,'YTick',-90:15:90,'fontsize',8)

% F VS wingPITCH&STROKE
figure
subplot(3,2,1)
title('Wing Pitch')
hold on
subplot(3,2,2)
title('Wing Stroke')
hold on

subplot(3,2,1)
plot(F_mean_wb(:),rad2deg(pitch_max_wb_L(:)),'.')
hold on
plot(F_mean_wb(:),rad2deg(pitch_max_wb_R(:)),'.r')
% xlabel('F')
ylabel('max wingpitch')
axis([.5 2.5 90 225])            
set(gca,'YTick',-90:45:225,'fontsize',8)

subplot(3,2,2)
plot(F_mean_wb(:),rad2deg(stroke_max_wb_L(:)),'.')
hold on
plot(F_mean_wb(:),rad2deg(stroke_max_wb_R(:)),'.r')
% xlabel('F')
ylabel('max wingstroke')
axis([.5 2.5 45 135])            
set(gca,'YTick',-90:45:225,'fontsize',8)

subplot(3,2,3)
plot(F_mean_wb(:),rad2deg(pitch_min_wb_L(:)),'.')
hold on
plot(F_mean_wb(:),rad2deg(pitch_min_wb_R(:)),'.r')
% xlabel('F')
ylabel('min wingpitch')
axis([.5 2.5 -45 90])            
set(gca,'YTick',-90:45:225,'fontsize',8)

subplot(3,2,4)
plot(F_mean_wb(:),rad2deg(stroke_min_wb_L(:)),'.')
hold on
plot(F_mean_wb(:),rad2deg(stroke_min_wb_R(:)),'.r')
% xlabel('F')
ylabel('min wingstroke')
axis([.5 2.5 -135 0])            
set(gca,'YTick',-180:45:225,'fontsize',8)

subplot(3,2,5)
plot(F_mean_wb(:),rad2deg(Apitch_wb_L(:)),'.')
hold on
plot(F_mean_wb(:),rad2deg(Apitch_wb_R(:)),'.r')
xlabel('F')
ylabel('Apitch')
axis([.5 2.5 45 225])            
set(gca,'YTick',-180:45:225,'fontsize',8)

subplot(3,2,6)
plot(F_mean_wb(:),rad2deg(Astroke_wb_L(:)),'.')
hold on
plot(F_mean_wb(:),rad2deg(Astroke_wb_R(:)),'.r')
xlabel('F')
ylabel('Astroke')
axis([.5 2.5 90 225])            
set(gca,'YTick',-180:45:225,'fontsize',8)

