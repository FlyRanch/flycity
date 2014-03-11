%% plot flightpath timeline tshift headingstart

% close all
figure
subplot(3,4,1)
title('body roll')
hold on
subplot(3,4,2)
title('body yaw')
hold on
subplot(3,4,3)
title('wing pitch')
hold on
subplot(3,4,4)
title('wing stroke')
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
        plot(t_wb_L(:,i)-t_shift(i),roll_mean_wb(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,4,5)
        plot(t_wb_L(:,i)-t_shift(i),roll_dot_mean_wb(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,4,9)
        plot(t_wb_L(:,i)-t_shift(i),roll_dot_dot_mean_wb(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        
        subplot(3,4,2)
        plot(t_wb_L(:,i)-t_shift(i),yaw_mean_wb(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,4,6)
        plot(t_wb_L(:,i)-t_shift(i),yaw_dot_mean_wb(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,4,10)
        plot(t_wb_L(:,i)-t_shift(i),yaw_dot_dot_mean_wb(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        
        subplot(3,4,3)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(pitch_max_wb_L(:,i)),'-','color','b','linewidth',linewidth_timelines)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(pitch_max_wb_R(:,i)),'-','color','r','linewidth',linewidth_timelines)
        subplot(3,4,7)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(pitch_min_wb_L(:,i)),'-','color','b','linewidth',linewidth_timelines)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(pitch_min_wb_R(:,i)),'-','color','r','linewidth',linewidth_timelines)
        subplot(3,4,11)
%         plot(t_wb_L(:,i)-t_shift(i),rad2deg(dpitch_mean_wb(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(dApitch_wb(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        
        subplot(3,4,4)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(stroke_max_wb_L(:,i)),'-','color','b','linewidth',linewidth_timelines)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(stroke_max_wb_R(:,i)),'-','color','r','linewidth',linewidth_timelines)
        subplot(3,4,8)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(stroke_min_wb_L(:,i)),'-','color','b','linewidth',linewidth_timelines)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(stroke_min_wb_R(:,i)),'-','color','r','linewidth',linewidth_timelines)
        subplot(3,4,12)
%         plot(t_wb_L(:,i)-t_shift(i),rad2deg(dstroke_mean_wb(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(dAstroke_wb(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
    end
end

subplot(3,4,1)
axis([t_start t_stop -45 135])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('roll','fontsize',10)
    grid on
subplot(3,4,5)
axis([t_start t_stop -4000 4000])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',[-4000:4000:4000],'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('roll rate','fontsize',10)
    grid on
subplot(3,4,9)
axis([t_start t_stop -300000 200000])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'YTick',[-300000;0;200000],'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('roll acceleration','fontsize',10)
    grid on

subplot(3,4,2)
axis([t_start t_stop -45 135])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('yaw','fontsize',10)
    grid on
subplot(3,4,6)
axis([t_start t_stop -4000 4000])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',[-4000:4000:4000],'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('yaw rate','fontsize',10)
    grid on
subplot(3,4,10)
axis([t_start t_stop -300000 200000])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'YTick',[-300000;0;200000],'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('yaw acceleration','fontsize',10)
    grid on

subplot(3,4,3)
axis([t_start t_stop 90 225])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',0:45:180,'fontsize',8)
%     xlabel('time','fontsize',10)
ylabel('max pitch','fontsize',10) 
    grid on
subplot(3,4,7)
axis([t_start t_stop -45 90])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',0:45:180,'fontsize',8)
%     xlabel('time','fontsize',10)
ylabel('min pitch','fontsize',10) 
    grid on
subplot(3,4,11)
axis([t_start t_stop -90 90])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-90:45:90,'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('dApitch L-R','fontsize',10)
    grid on
    

subplot(3,4,4)
axis([t_start t_stop 0 135])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-90:45:180,'fontsize',8)
%     xlabel('time','fontsize',10)
ylabel('max stroke','fontsize',10) 
    grid on
subplot(3,4,8)
axis([t_start t_stop -135 0])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-180:45:180,'fontsize',8)
%     xlabel('time','fontsize',10)
ylabel('min stroke','fontsize',10) 
    grid on
subplot(3,4,12)
axis([t_start t_stop -45 45])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-90:45:90,'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('dAstroke L-R','fontsize',10)
    grid on
    
