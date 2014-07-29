% YAWacceleration VS DEV L-R
figure
subplot(3,2,1)
plot(yaw_dot_dot_mean_wb(:),rad2deg(ddev_max_dus(:)),'.')
% xlabel('yaw acceleration')
title('Ddevmax down2up L-R')
axis([-100000 200000 -50 50])            

subplot(3,2,2)
plot(yaw_dot_dot_mean_wb(:),rad2deg(ddev_max_udsPREV(:)),'.')
% xlabel('yaw acceleration')
title('Ddevmax up2down L-R')
axis([-100000 200000 -50 50])            

subplot(3,2,3)
plot(yaw_dot_dot_mean_wb(:),rad2deg(ddev_min_ds(:)),'.r')
% xlabel('yaw acceleration')
title('Ddevmin down L-R')
axis([-100000 200000 -50 50])            

subplot(3,2,4)
plot(yaw_dot_dot_mean_wb(:),rad2deg(ddev_min_us(:)),'.r')
% xlabel('yaw acceleration')
title('Ddevmin up L-R')
axis([-100000 200000 -50 50])            

subplot(3,2,5)
plot(yaw_dot_dot_mean_wb(:),rad2deg(dAdev_ds(:)),'.k')
xlabel('yaw acceleration')
title('dAdev downstroke L-R')
axis([-100000 200000 -50 50])            

subplot(3,2,6)
plot(yaw_dot_dot_mean_wb(:),rad2deg(dAdev_us(:)),'.k')
xlabel('yaw acceleration')
title('dAdev upstroke L-R')
axis([-100000 200000 -50 50])            

% YAW acceleration VS PITCH L-R
figure
subplot(2,2,1)
plot(yaw_dot_dot_mean_wb(:),rad2deg(dpitch_max_ds(:)),'.')
% xlabel('yaw acceleration')
title('Dpitchmax downstroke L-R')
axis([-100000 200000 -50 50])            

subplot(2,2,2)
plot(yaw_dot_dot_mean_wb(:),rad2deg(dpitch_min_us(:)),'.')
% xlabel('yaw acceleration')
title('Dpitchmin upstroke L-R')
axis([-100000 200000 -50 50])            

subplot(2,2,3)
plot(yaw_dot_dot_mean_wb(:),rad2deg(dApitch_ds(:)),'.k')
xlabel('yaw acceleration')
title('dApitch downstroke L-R')
axis([-100000 200000 -50 50])            

subplot(2,2,4)
plot(yaw_dot_dot_mean_wb(:),rad2deg(dApitch_us(:)),'.k')
xlabel('yaw acceleration')
title('dApitch upstroke L-R')
axis([-100000 200000 -50 50])            

% YAW acceleration VS STROKE L-R
figure
subplot(2,2,1)
plot(yaw_dot_dot_mean_wb(:),rad2deg(dstroke_max_wb(:)),'.')
% xlabel('yaw acceleration')
title('Dstrokemax L-R')
axis([-100000 200000 -30 30])            

subplot(2,2,2)
plot(yaw_dot_dot_mean_wb(:),rad2deg(dstroke_min_wb(:)),'.')
% xlabel('yaw acceleration')
title('Dstrokemin L-R')
axis([-100000 200000 -30 30])            

subplot(2,2,3)
plot(yaw_dot_dot_mean_wb(:),rad2deg(dAstroke_ds(:)),'.k')
xlabel('yaw acceleration')
title('dAstroke downstroke L-R')
axis([-100000 200000 -30 30])            

subplot(2,2,4)
plot(yaw_dot_dot_mean_wb(:),rad2deg(dAstroke_us(:)),'.k')
xlabel('yaw acceleration')
title('dAstroke upstroke L-R')
axis([-100000 200000 -30 30])            
