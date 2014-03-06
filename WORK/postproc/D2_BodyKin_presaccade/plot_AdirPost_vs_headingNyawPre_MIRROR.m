subplot(1,2,1)
plotcolor = 'r'
angle_pre = heading_pre;
angle_post = Adir_post;
% mirror angles
for i = 1:length(angle_pre)
    if An_hor_max(i) < 0
        angle_pre(i) = -angle_pre(i);
        angle_post(i) = -angle_post(i);
    end
%     if An_hor_max(i) < 0  && settings.expansion.HorPos(i) == 0
%         angle_pre(i) = -angle_pre(i);
%         angle_post(i) = -angle_post(i);
%     elseif An_hor_max(i) > 0  && settings.expansion.HorPos(i) == 180
%         angle_pre(i) = -angle_pre(i);
%         angle_post(i) = -angle_post(i);
%     end
end
% plot_angle_pre_post_csaps
% plot_angle_pre_post_circmean
plot_angle_pre_post_circmean_extendedsection

xlabel('initial heading','fontsize',10) 
ylabel('escape A direction','fontsize',10) 
% set(gca,'xlim',[0 180],'ylim',[-180 180])
set(gca,'xlim',[-225 45],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 
grid on

subplot(1,2,2)
plotcolor = 'b'
angle_pre = yaw_pre;
angle_post = Adir_post;
% mirror angles
for i = 1:length(angle_pre)
    if An_hor_max(i) < 0
        angle_pre(i) = -angle_pre(i);
        angle_post(i) = -angle_post(i);
    end
%     if An_hor_max(i) < 0  && settings.expansion.HorPos(i) == 0
%         angle_pre(i) = -angle_pre(i);
%         angle_post(i) = -angle_post(i);
%     elseif An_hor_max(i) > 0  && settings.expansion.HorPos(i) == 180
%         angle_pre(i) = -angle_pre(i);
%         angle_post(i) = -angle_post(i);
%     end
end
% plot_angle_pre_post_csaps
% plot_angle_pre_post_circmean
plot_angle_pre_post_circmean_extendedsection

xlabel('initial yaw','fontsize',10) 
ylabel('escape A direction','fontsize',10) 
% set(gca,'xlim',[0 180],'ylim',[-180 180])
set(gca,'xlim',[-225 45],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 
grid on