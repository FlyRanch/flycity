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

for j=1:length(subset_seqs)
    i=subset_seqs(j);
    size(stim_angle_vel,2) - i;
    
    if isnan(color_var(i)) == 0 && color_var(i)~=0
        subplot(3,3,1)
        plot(t-t_shift(i),roll_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,3,4)
        plot(t-t_shift(i),pitch_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,3,7)
        plot(t-t_shift(i),yaw_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
    end
end

subplot(3,3,1)
axis([t_start t_stop -45 90])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',-90:45:90,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('roll [deg]','fontsize',10)
%    grid on 
subplot(3,3,4)
axis([t_start t_stop -45 90])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',-90:45:90,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('pitch [deg]','fontsize',10)
%    grid on 
subplot(3,3,7)
axis([t_start t_stop -45 90])
    set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'YTick',-90:45:90,'fontsize',8)
    xlabel('time [s]','fontsize',10)
    ylabel('yaw [deg]','fontsize',10)
%    grid on 
