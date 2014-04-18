%% plot flightpath timeline tshift headingstart

% close all
figure
subplot(3,3,1)
% title('strokeplane angles','fontsize',10)
hold on
subplot(3,3,2)
% title('body euler angles','fontsize',10)
hold on
subplot(3,3,3)
hold on
subplot(3,3,4)
hold on
subplot(3,3,5)
hold on
subplot(3,3,6)
hold on
subplot(3,3,7)
hold on
subplot(3,3,8)
hold on
subplot(3,3,9)
hold on

for i=1:size(stim_angle_vel,2)
    size(stim_angle_vel,2) - i
    
    if isnan(color_var(i)) == 0 && color_var(i)~=0
        subplot(3,3,1)
        plot(t_wb_L(:,i)-t_shift(i),roll_dot_dot_mean_wb_norm(:,i),'-','color',grey_color,'linewidth',.25)
        subplot(3,3,4)
        plot(t_wb_L(:,i)-t_shift(i),pitch_dot_dot_mean_wb_norm(:,i),'-','color',grey_color,'linewidth',.25)
        subplot(3,3,7)
        plot(t_wb_L(:,i)-t_shift(i),yaw_dot_dot_mean_wb_norm(:,i),'-','color',grey_color,'linewidth',.25)
        
    end
end

for j=1:length(subset_seqs)
    i=subset_seqs(j);
    size(stim_angle_vel,2) - i;
    
    if isnan(color_var(i)) == 0 && color_var(i)~=0
        subplot(3,3,1)
        plot(t_wb_L(:,i)-t_shift(i),roll_dot_dot_mean_wb_norm(:,i),'-o','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines,'MarkerFaceColor',cmap_plot(color_var(i),:),'markersize',2)
        subplot(3,3,4)
        plot(t_wb_L(:,i)-t_shift(i),pitch_dot_dot_mean_wb_norm(:,i),'-o','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines,'MarkerFaceColor',cmap_plot(color_var(i),:),'markersize',2)
        subplot(3,3,7)
        plot(t_wb_L(:,i)-t_shift(i),yaw_dot_dot_mean_wb_norm(:,i),'-o','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines,'MarkerFaceColor',cmap_plot(color_var(i),:),'markersize',2)
        
%         subplot(3,3,1)
%         plot(t_wb_L(:,i)-t_shift(i),roll_dot_dot_mean_wb_norm(:,i),'-.','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
%         subplot(3,3,4)
%         plot(t_wb_L(:,i)-t_shift(i),pitch_dot_dot_mean_wb_norm(:,i),'-.','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
%         subplot(3,3,7)
%         plot(t_wb_L(:,i)-t_shift(i),yaw_dot_dot_mean_wb_norm(:,i),'-.','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
    end
end

subplot(3,3,1)
axis([t_start t_stop -7.5 7.5])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',-5:2.5:5,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('roll accel norm','fontsize',10)
%    grid on 
subplot(3,3,4)
axis([t_start t_stop -7.5 7.5])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',-5:2.5:5,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('pitch accel norm','fontsize',10)
%    grid on 
subplot(3,3,7)
axis([t_start t_stop -7.5 7.5])
    set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'YTick',-5:2.5:5,'fontsize',8)
    xlabel('time [s]','fontsize',10)
    ylabel('yaw accel norm','fontsize',10)
%    grid on 
