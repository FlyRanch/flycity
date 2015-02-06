subplot(3,1,1)
% plot(360,360,'ok','MarkerFaceColor','r','markersize',5)
% hold on
% plot(360,360,'ok','MarkerFaceColor','b','markersize',5)
% legend('@Amax','mean')
% 
% plotcolor = 'r'
% angle_pre = heading_pre;
% angle_post = roll_Amax;
% % mirror angles
% angle_post = abs(angle_post);
% for i = 1:length(angle_pre)
%     if An_hor_max(i) < 0
%         angle_pre(i) = -angle_pre(i);
%     end
% end
% plot_angle_pre_post_circmean_extendedsection

plotcolor = 'b'
angle_pre = heading_pre;
angle_post = roll_mean;
% mirror angles
angle_post = abs(angle_post);
for i = 1:length(angle_pre)
    if An_hor_max(i) < 0
        angle_pre(i) = -angle_pre(i);
    end
end
plot_angle_pre_post_circmean_extendedsection

xlabel('initial heading','fontsize',10) 
ylabel('roll','fontsize',10) 
set(gca,'xlim',[-225 45],'ylim',[-45 90])
set(gca,'XTick',[-180:90:180])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
grid on
% axis equal

subplot(3,1,2)

% plotcolor = 'r'
% angle_pre = heading_pre;
% angle_post = pitch_Amax - pitch_steady_mean;
% % mirror angles
% for i = 1:length(angle_pre)
%     if An_hor_max(i) < 0
%         angle_pre(i) = -angle_pre(i);
%     end
% end
% plot_angle_pre_post_circmean_extendedsection

plotcolor = 'b'
angle_pre = heading_pre;
% angle_post = pitch_mean - pitch_steady_mean;
angle_post = pitch_mean;
% mirror angles
for i = 1:length(angle_pre)
    if An_hor_max(i) < 0
        angle_pre(i) = -angle_pre(i);
    end
end
plot_angle_pre_post_circmean_extendedsection

xlabel('initial heading','fontsize',10) 
ylabel('Dpitch','fontsize',10) 
set(gca,'xlim',[-225 45],'ylim',[-45 90])
set(gca,'XTick',[-180:90:180])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
grid on
% axis equal

subplot(3,1,3)

% plotcolor = 'r'
% angle_pre = heading_pre;
% angle_post = yaw_Amax;
% % mirror angles
% angle_post = -abs(angle_post);
% for i = 1:length(angle_pre)
%     if An_hor_max(i) < 0
%         angle_pre(i) = -angle_pre(i);
%     end
% end
% plot_angle_pre_post_circmean_extendedsection

plotcolor = 'b'
angle_pre = heading_pre;
angle_post = yaw_mean;
% mirror angles
% angle_post = -abs(angle_post);
for i = 1:length(angle_pre)
    if An_hor_max(i) < 0
        angle_pre(i) = -angle_pre(i);
        angle_post(i) = -angle_post(i);
    end
end
plot_angle_pre_post_circmean_extendedsection

xlabel('initial heading','fontsize',10) 
ylabel('yaw','fontsize',10) 
set(gca,'xlim',[-225 45],'ylim',[-45 90])
set(gca,'XTick',[-180:90:180])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
grid on
% axis equal

