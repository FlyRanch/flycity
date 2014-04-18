%% plot flightpath timeline tshift headingstart

% close all
figure
subplot(2,2,1)
hold on
subplot(2,2,2)
hold on
subplot(2,2,3)
hold on
subplot(2,2,4)
hold on

for j=1:length(seqs)
    i=seqs(j)
    size(stim_angle_vel,2) - i
    
    if isnan(color_var(i)) == 0 && color_var(i)~=0
        subplot(2,2,1)
        plot(t-t_shift(i),roll_dot_dot_norm_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(2,2,2)
        plot(t-t_shift(i),pitch_dot_dot_norm_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(2,2,3)
        plot(t-t_shift(i),yaw_dot_dot_norm_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(2,2,4)
        plot(t-t_shift(i),F_norm_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
    end
end

subplot(2,2,1)
axis([t_start t_stop -3 3])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',-3:3:3,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('MOD value roll accel','fontsize',10)
    grid on

subplot(2,2,2)
axis([t_start t_stop -3 3])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',-3:3:3,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('MOD value pitch accel','fontsize',10)
    grid on

subplot(2,2,3)
axis([t_start t_stop -3 3])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',-3:3:3,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('MOD value yaw accel','fontsize',10)
    grid on

subplot(2,2,4)
axis([t_start t_stop -1 2])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',-1:1:2,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('MOD value flight force','fontsize',10)
    grid on


