%% plot flightpath timeline tshift headingstart

% close all
figure
subplot(3,2,1)
hold on
subplot(3,2,3)
hold on
subplot(3,2,5)
hold on
subplot(3,2,2)
hold on
subplot(3,2,4)
hold on
subplot(3,2,6)
hold on

for i=1:size(stim_angle_vel,2)
    size(stim_angle_vel,2) - i
    
    if isnan(color_var(i)) == 0 && color_var(i)~=0
        subplot(3,2,1)
        plot(t-t_shift(i),roll_dot_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,2,3)
        plot(t-t_shift(i),pitch_dot_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,2,5)
        plot(t-t_shift(i),yaw_dot_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)

        subplot(3,2,2)
        plot(t-t_shift(i),roll_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,2,4)
        plot(t-t_shift(i),pitch_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,2,6)
        plot(t-t_shift(i),yaw_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
    end
end

subplot(3,2,1)
axis([t_start t_stop -4000 4000])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',[-4000:4000:4000],'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('roll rate','fontsize',10)
    grid on
subplot(3,2,3)
axis([t_start t_stop -1000 3000])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',[-1000;0;3000],'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('pitch rate','fontsize',10)
    grid on
subplot(3,2,5)
axis([t_start t_stop -2000 2000])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',[-2000;0;2000],'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('yaw rate','fontsize',10)
    grid on

subplot(3,2,2)
axis([t_start t_stop -90 90])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('roll','fontsize',10)
    grid on
subplot(3,2,4)
axis([t_start t_stop -45 135])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',[0;90;180],'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('pitch','fontsize',10)
    grid on
subplot(3,2,6)
axis([t_start t_stop -90 90])
    set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('yaw','fontsize',10)
    grid on


