

% F VS DEV
figure
subplot(3,2,1)
plot(F_mean_wb(:),rad2deg(dev_max_dus_L(:)),'.')
hold on
plot(F_mean_wb(:),rad2deg(dev_max_dus_R(:)),'.r')
% xlabel('F')
title('max dev down2up')
axis([.5 2.5 -30 60])            

subplot(3,2,2)
plot(F_mean_wb(:),rad2deg(dev_max_udsPREV_L(:)),'.')
hold on
plot(F_mean_wb(:),rad2deg(dev_max_udsPREV_R(:)),'.r')
% xlabel('F')
title('max dev up2down')
axis([.5 2.5 -30 60])            

subplot(3,2,3)
plot(F_mean_wb(:),rad2deg(dev_min_ds_L(:)),'.')
hold on
plot(F_mean_wb(:),rad2deg(dev_min_ds_R(:)),'.r')
xlabel('F')
title('min dev downstroke')
axis([.5 2.5 -45 15])            

subplot(3,2,4)
plot(F_mean_wb(:),rad2deg(dev_min_us_L(:)),'.')
hold on
plot(F_mean_wb(:),rad2deg(dev_min_us_R(:)),'.r')
xlabel('F')
title('min dev upstroke')
axis([.5 2.5 -45 15])            

subplot(3,2,5)
plot(F_mean_wb(:),rad2deg(Adev_ds_L(:)),'.')
hold on
plot(F_mean_wb(:),rad2deg(Adev_ds_R(:)),'.r')
xlabel('F')
title('Adev downstroke')
axis([.5 2.5 -0 50])            

subplot(3,2,6)
plot(F_mean_wb(:),rad2deg(Adev_us_L(:)),'.')
hold on
plot(F_mean_wb(:),rad2deg(Adev_us_R(:)),'.r')
xlabel('F')
title('Adev upstroke')
axis([.5 2.5 -0 50])            

% F VS wingPITCH&STROKE
figure
subplot(3,2,1)
plot(F_mean_wb(:),rad2deg(pitch_max_wb_L(:)),'.')
hold on
plot(F_mean_wb(:),rad2deg(pitch_max_wb_R(:)),'.r')
% xlabel('F')
title('max wingpitch')
axis([.5 2.5 120 210])            

subplot(3,2,2)
plot(F_mean_wb(:),rad2deg(stroke_max_wb_L(:)),'.')
hold on
plot(F_mean_wb(:),rad2deg(stroke_max_wb_R(:)),'.r')
% xlabel('F')
title('max wingstroke')
axis([.5 2.5 45 135])            

subplot(3,2,3)
plot(F_mean_wb(:),rad2deg(pitch_min_wb_L(:)),'.')
hold on
plot(F_mean_wb(:),rad2deg(pitch_min_wb_R(:)),'.r')
% xlabel('F')
title('min wingpitch')
axis([.5 2.5 -45 45])            

subplot(3,2,4)
plot(F_mean_wb(:),rad2deg(stroke_min_wb_L(:)),'.')
hold on
plot(F_mean_wb(:),rad2deg(stroke_min_wb_R(:)),'.r')
% xlabel('F')
title('min wingstroke')
axis([.5 2.5 -105 -15])            

subplot(3,2,5)
plot(F_mean_wb(:),rad2deg(Apitch_wb_L(:)),'.')
hold on
plot(F_mean_wb(:),rad2deg(Apitch_wb_R(:)),'.r')
xlabel('F')
title('Apitch')
axis([.5 2.5 90 180])            

subplot(3,2,6)
plot(F_mean_wb(:),rad2deg(Astroke_wb_L(:)),'.')
hold on
plot(F_mean_wb(:),rad2deg(Astroke_wb_R(:)),'.r')
xlabel('F')
title('Astroke')
axis([.5 2.5 90 180])            


