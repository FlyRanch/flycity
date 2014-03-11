% ROLL acceleration VS DEV L-R
figure
subplot(3,2,1)
plot(roll_dot_dot_mean_wb(:),rad2deg(ddev_max_dus(:)),'.')
% xlabel('roll acceleration')
title('Ddevmax down2up L-R')
axis([-300000 100000 -50 50])            

subplot(3,2,2)
plot(roll_dot_dot_mean_wb(:),rad2deg(ddev_max_udsPREV(:)),'.')
% xlabel('roll acceleration')
title('Ddevmax up2down L-R')
axis([-300000 100000 -50 50])            

subplot(3,2,3)
plot(roll_dot_dot_mean_wb(:),rad2deg(ddev_min_ds(:)),'.r')
% xlabel('roll acceleration')
title('Ddevmin down L-R')
axis([-300000 100000 -50 50])            

subplot(3,2,4)
plot(roll_dot_dot_mean_wb(:),rad2deg(ddev_min_us(:)),'.r')
% xlabel('roll acceleration')
title('Ddevmin up L-R')
axis([-300000 100000 -50 50])            

subplot(3,2,5)
plot(roll_dot_dot_mean_wb(:),rad2deg(dAdev_ds(:)),'.k')
xlabel('roll acceleration')
title('dAdev downstroke L-R')
axis([-300000 100000      -50 50])            

subplot(3,2,6)
plot(roll_dot_dot_mean_wb(:),rad2deg(dAdev_us(:)),'.k')
xlabel('roll acceleration')
title('dAdev upstroke L-R')
axis([-300000 100000      -50 50])            
