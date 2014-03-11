% roll rate timelines
figure
for i=1:80
    subplot(3,2,1)
%     hold off
    plot(t_wb_L(:,i)-t_pre(i),roll_dot_dot_mean_wb(:,i),'-k')
    hold on
    subplot(3,2,3)
%     hold off
    plot(t_wb_L(:,i)-t_pre(i),pitch_dot_dot_mean_wb(:,i),'-k')
    hold on
    subplot(3,2,5)
%     hold off
    plot(t_wb_L(:,i)-t_pre(i),yaw_dot_dot_mean_wb(:,i),'-k')
    hold on
    
    
    
    subplot(3,2,2)
%     hold off
    plot(t_wb_L(:,i)-t_pre(i),rad2deg(dev_max_dus_L(:,i)),'-b')
    hold on
    plot(t_wb_R(:,i)-t_pre(i),rad2deg(dev_max_dus_R(:,i)),'-r')
    subplot(3,2,4)
%     hold off
    plot(t_wb_L(:,i)-t_pre(i),rad2deg(ddev_max_dus(:,i)),'-k')
    hold on
    subplot(3,2,6)
%     hold off
    plot(t_wb_L(:,i)-t_pre(i),rad2deg(dAdev_us(:,i)),'-k')
    hold on
    
    for j = 1:150
        if isnan(ddev_max_dus(j,i)) == 0
            subplot(3,2,1)
            plot(t_wb_L(j,i)-t_pre(i),roll_dot_dot_mean_wb(j,i),'.','color',settings.cmap_k(IDX_wb_L(j,i),:),'markersize',20)
            subplot(3,2,3)
            plot(t_wb_L(j,i)-t_pre(i),pitch_dot_dot_mean_wb(j,i),'.','color',settings.cmap_k(IDX_wb_L(j,i),:),'markersize',20)
            subplot(3,2,5)
            plot(t_wb_L(j,i)-t_pre(i),yaw_dot_dot_mean_wb(j,i),'.','color',settings.cmap_k(IDX_wb_L(j,i),:),'markersize',20)

%             subplot(3,2,2)
%             plot(t_wb_L(j,i)-t_pre(i),rad2deg(dev_max_dus_L(j,i)),'.','color',settings.cmap_k(IDX_wb_L(j,i),:),'markersize',5)
%             plot(t_wb_R(j,i)-t_pre(i),rad2deg(dev_max_dus_R(j,i)),'.','color',settings.cmap_k(IDX_wb_R(j,i),:),'markersize',5)
            subplot(3,2,4)
            plot(t_wb_L(j,i)-t_pre(i),rad2deg(ddev_max_dus(j,i)),'.','color',settings.cmap_k(IDX_wb_L(j,i),:),'markersize',20)
            subplot(3,2,6)
            plot(t_wb_L(j,i)-t_pre(i),rad2deg(dAdev_us(j,i)),'.','color',settings.cmap_k(IDX_wb_L(j,i),:),'markersize',20)
        end
    end
%         pause
end

