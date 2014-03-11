%% plot flightpath timeline tshift headingstart

% close all
figure
subplot(3,2,1)
title('strokeplane angles','fontsize',10)
hold on
subplot(3,2,2)
title('body euler angles','fontsize',10)
hold on
subplot(3,2,3)
hold on
subplot(3,2,4)
hold on
subplot(3,2,5)
hold on
subplot(3,2,6)
hold on

for i=1:size(stim_angle_vel,2)
    size(stim_angle_vel,2) - i
    
    if isnan(color_var(i)) == 0 && color_var(i)~=0
        subplot(3,2,1)
        plot(t-t_shift(i),roll_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,2,3)
        plot(t-t_shift(i),pitch_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,2,5)
        plot(t-t_shift(i),yaw_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        
        subplot(3,2,2)
        plot(t-t_shift(i),roll_global_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,2,4)
        plot(t-t_shift(i),pitch_global_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,2,6)
        plot(t-t_shift(i),slip_global_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
    end
end

subplot(3,2,1)
axis([t_start t_stop -135 136])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('roll','fontsize',10)
    grid on
subplot(3,2,3)
axis([t_start t_stop -45 135])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',[0;90;180],'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('pitch','fontsize',10)
    grid on
subplot(3,2,5)
axis([t_start t_stop -90 90])
    set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('yaw','fontsize',10)
    grid on

subplot(3,2,2)
axis([t_start t_stop -135 135])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('euler roll','fontsize',10)
    grid on
subplot(3,2,4)
axis([t_start t_stop -45 135])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',[0;90;180],'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('euler pitch','fontsize',10)
    grid on
subplot(3,2,6)
axis([t_start t_stop -90 90])
    set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
    xlabel('time','fontsize',10)
%     ylabel('euler slip','fontsize',10)
    grid on


