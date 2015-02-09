subplot(2,2,1)
plot(360,1,'ok','MarkerFaceColor','r','markersize',5)
hold on
plot(360,1,'ok','MarkerFaceColor','b','markersize',5)
legend('right turn','left turn','location','ne')

plotcolor = 'r'
angle_pre = heading_pre;
angle_post = heading_post;
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
plot_angle_pre_post_circmean

hold on
plotcolor = 'b'
angle_pre = heading_pre;
angle_post = heading_post;
% mirror angles
for i = 1:length(angle_pre)
    if An_hor_max(i) > 0
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
plot_angle_pre_post_circmean

% xlabel('initial heading','fontsize',10) 
ylabel('escape heading','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-180 180]) 
% set(gca,'xlim',[0 180],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180],'xticklabel',[])
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 
grid on

subplot(2,2,3)
plotcolor = 'r'
angle_pre = heading_pre;
angle_post = yaw_post;
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
plot_angle_pre_post_circmean

plotcolor = 'b'
angle_pre = heading_pre;
angle_post = yaw_post;
% mirror angles
for i = 1:length(angle_pre)
    if An_hor_max(i) > 0
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
plot_angle_pre_post_circmean

xlabel('initial heading','fontsize',10) 
ylabel('escape yaw','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-180 180]) 
% set(gca,'xlim',[0 180],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 
grid on

subplot(2,2,2)
plotcolor = 'r'
angle_pre = yaw_pre;
angle_post = heading_post;
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
plot_angle_pre_post_circmean

plotcolor = 'b'
angle_pre = yaw_pre;
angle_post = heading_post;
% mirror angles
for i = 1:length(angle_pre)
    if An_hor_max(i) > 0
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
plot_angle_pre_post_circmean

% xlabel('initial yaw','fontsize',10) 
% ylabel('escape heading','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-180 180]) 
% set(gca,'xlim',[0 180],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180],'xticklabel',[])
set(gca,'YTick',[-180 -90 0 90 180],'yticklabel',[],'fontsize',8) 
grid on

subplot(2,2,4)
plotcolor = 'r'
angle_pre = yaw_pre;
angle_post = yaw_post;
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
plot_angle_pre_post_circmean

plotcolor = 'b'
angle_pre = yaw_pre;
angle_post = yaw_post;
% mirror angles
for i = 1:length(angle_pre)
    if An_hor_max(i) > 0
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
plot_angle_pre_post_circmean

xlabel('initial yaw','fontsize',10) 
% ylabel('escape yaw','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-180 180]) 
% set(gca,'xlim',[0 180],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180],'yticklabel',[],'fontsize',8) 
grid on
