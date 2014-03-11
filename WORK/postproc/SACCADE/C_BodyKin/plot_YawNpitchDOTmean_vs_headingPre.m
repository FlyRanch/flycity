subplot(2,1,1)

plotcolor = 'r'
angle_pre = heading_pre;
angle_post = pitch_dot_mean;
% mirror angles
for i = 1:length(angle_pre)
    if An_hor_max(i) < 0
        angle_pre(i) = -angle_pre(i);
    end
end
plot_Xangle_pre_post_mean_extendedsection

xlabel('initial heading','fontsize',10) 
ylabel('mean pitch rate','fontsize',10) 
set(gca,'xlim',[-225 45],'ylim',[-2000 4000])
set(gca,'XTick',[-180:90:180])
set(gca,'YTick',[-2000:2000:4000],'fontsize',8) 
grid on
% axis equal

subplot(2,1,2)

plotcolor = 'r'
angle_pre = heading_pre;
angle_post = yaw_dot_mean;
% mirror angles
% angle_post = -abs(angle_post);
for i = 1:length(angle_pre)
    if An_hor_max(i) < 0
        angle_pre(i) = -angle_pre(i);
        angle_post(i) = -angle_post(i);
    end
end
plot_Xangle_pre_post_mean_extendedsection

xlabel('initial heading','fontsize',10) 
ylabel('mean yaw rate','fontsize',10) 
set(gca,'xlim',[-225 45],'ylim',[-2000 4000])
set(gca,'XTick',[-180:90:180])
set(gca,'YTick',[-2000:2000:4000],'fontsize',8) 
grid on
% axis equal

