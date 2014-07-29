%% plot flightpath timeline tshift headingstart

% close all
figure
subplot(3,4,1)
title('downstroke max')
hold on
subplot(3,4,2)
title('downstroke min')
hold on
subplot(3,4,3)
title('upstroke max')
hold on
subplot(3,4,4)
title('upstroke min')
hold on
subplot(3,4,5)
hold on
subplot(3,4,6)
hold on
subplot(3,4,7)
hold on
subplot(3,4,8)
hold on
subplot(3,4,9)
hold on
subplot(3,4,10)
hold on
subplot(3,4,11)
hold on
subplot(3,4,12)
hold on

for i=1:size(stim_angle_vel,2)
    size(stim_angle_vel,2) - i
    
    if isnan(color_var(i)) == 0 && color_var(i)~=0
        subplot(3,4,1)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(dev_max_udsPREV_L(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,4,5)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(dev_max_udsPREV_R(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,4,2)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(dev_min_ds_L(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,4,6)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(dev_min_ds_R(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        
        subplot(3,4,3)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(dev_max_dus_L(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,4,7)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(dev_max_dus_R(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,4,4)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(dev_min_us_L(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,4,8)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(dev_min_us_R(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        
        subplot(3,4,9)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(ddev_max_udsPREV(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,4,10)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(ddev_min_ds(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        
        subplot(3,4,11)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(ddev_max_dus(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,4,12)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(ddev_min_us(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
    end
end

subplot(3,4,1)
axis([t_start t_stop -45 45])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-45:45:90,'fontsize',8)
%     xlabel('time','fontsize',10)
ylabel('left deviation','fontsize',10) 
    grid on
subplot(3,4,2)
axis([t_start t_stop -45 45])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-45:45:90,'fontsize',8)
%     xlabel('time','fontsize',10)
% ylabel('left deviation','fontsize',10) 
    grid on
subplot(3,4,3)
axis([t_start t_stop -45 45])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-45:45:90,'fontsize',8)
%     xlabel('time','fontsize',10)
% ylabel('left deviation','fontsize',10) 
    grid on
subplot(3,4,4)
axis([t_start t_stop -45 45])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-45:45:90,'fontsize',8)
%     xlabel('time','fontsize',10)
% ylabel('left deviation','fontsize',10) 
    grid on
    
subplot(3,4,5)
axis([t_start t_stop -45 45])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-45:45:45,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('right deviation','fontsize',10)
    grid on
subplot(3,4,6)
axis([t_start t_stop -45 45])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-45:45:45,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('right deviation','fontsize',10)
    grid on
subplot(3,4,7)
axis([t_start t_stop -45 45])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-45:45:45,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('right deviation','fontsize',10)
    grid on
subplot(3,4,8)
axis([t_start t_stop -45 45])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-45:45:45,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('right deviation','fontsize',10)
    grid on
    
subplot(3,4,9)
axis([t_start t_stop -45 45])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-45:45:45,'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('dev diff L-R','fontsize',10)
    grid on
subplot(3,4,10)
axis([t_start t_stop -45 45])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-45:45:45,'fontsize',8)
    xlabel('time','fontsize',10)
%     ylabel('dev diff L-R','fontsize',10)
    grid on
subplot(3,4,11)
axis([t_start t_stop -45 45])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-45:45:45,'fontsize',8)
    xlabel('time','fontsize',10)
%     ylabel('dev diff L-R','fontsize',10)
    grid on
subplot(3,4,12)
axis([t_start t_stop -45 45])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-45:45:45,'fontsize',8)
    xlabel('time','fontsize',10)
%     ylabel('dev diff L-R','fontsize',10)
    grid on

