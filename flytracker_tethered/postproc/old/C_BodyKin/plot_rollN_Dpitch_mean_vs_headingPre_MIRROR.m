subplot(2,1,1)
plotcolor = 'r'
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
ylabel('mean roll','fontsize',10) 
set(gca,'xlim',[-225 45],'ylim',[0 90])
set(gca,'XTick',[-180:90:180])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
grid on
% axis equal

subplot(2,1,2)
plotcolor = 'b'
angle_pre = heading_pre;
angle_post = pitch_mean - pitch_steady_mean;
% mirror angles
for i = 1:length(angle_pre)
    if An_hor_max(i) < 0
        angle_pre(i) = -angle_pre(i);
    end
end

plot_angle_pre_post_circmean_extendedsection

xlabel('initial heading','fontsize',10) 
ylabel('meanDpitch','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-180 180]) 
% set(gca,'xlim',[0 180],'ylim',[-180 180])
set(gca,'xlim',[-225 45],'ylim',[-45 45])
set(gca,'XTick',[-180:90:180])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
grid on
% axis equal

