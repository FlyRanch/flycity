%% plot flightpath timeline tshift headingstart

% close all
figure
subplot(3,4,1)
title('body pitch')
hold on
subplot(3,4,2)
title('body F')
hold on
subplot(3,4,3)
title('downstroke deviation')
hold on
subplot(3,4,4)
title('upstroke deviation')
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
        plot(t_wb_L(:,i)-t_shift(i),pitch_mean_wb(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,4,5)
        plot(t_wb_L(:,i)-t_shift(i),pitch_dot_mean_wb(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,4,9)
        plot(t_wb_L(:,i)-t_shift(i),pitch_dot_dot_mean_wb(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        
        subplot(3,4,2)
        plot(t_wb_L(:,i)-t_shift(i),F_mean_wb(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
%         subplot(3,4,6)
%         plot(t_wb_L(:,i)-t_shift(i),Fsp_roll_mean_wb(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
%         subplot(3,4,10)
%         plot(t_wb_L(:,i)-t_shift(i),Fsp_pitch_mean_wb(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        
        subplot(3,4,3)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(dev_max_udsPREV_L(:,i)),'-','color','b','linewidth',linewidth_timelines)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(dev_max_udsPREV_R(:,i)),'-','color','r','linewidth',linewidth_timelines)
        subplot(3,4,7)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(dev_min_ds_L(:,i)),'-','color','b','linewidth',linewidth_timelines)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(dev_min_ds_R(:,i)),'-','color','r','linewidth',linewidth_timelines)
        subplot(3,4,11)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(dAdev_ds(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        
        subplot(3,4,4)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(dev_max_dus_L(:,i)),'-','color','b','linewidth',linewidth_timelines)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(dev_max_dus_R(:,i)),'-','color','r','linewidth',linewidth_timelines)
        subplot(3,4,8)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(dev_min_us_L(:,i)),'-','color','b','linewidth',linewidth_timelines)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(dev_min_us_R(:,i)),'-','color','r','linewidth',linewidth_timelines)
        subplot(3,4,12)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(dAdev_us(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
    end
end

subplot(3,4,1)
axis([t_start t_stop -45 135])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('pitch','fontsize',10)
    grid on
subplot(3,4,5)
axis([t_start t_stop -4000 4000])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',[-4000:4000:4000],'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('pitch rate','fontsize',10)
    grid on
subplot(3,4,9)
axis([t_start t_stop -300000 200000])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'YTick',[-300000;0;200000],'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('pitch acceleration','fontsize',10)
    grid on

subplot(3,4,2)
axis([t_start t_stop -45 135])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('F','fontsize',10)
    grid on
subplot(3,4,6)
axis([t_start t_stop -4000 4000])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop,'XTickLabel',[]) 
    set(gca,'YTick',[-4000:4000:4000],'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('F roll','fontsize',10)
    grid on
subplot(3,4,10)
axis([t_start t_stop -300000 200000])
%     set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'XTick',t_start:-t_start:t_stop) 
    set(gca,'YTick',[-300000;0;200000],'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('F pitch','fontsize',10)
    grid on

subplot(3,4,3)
axis([t_start t_stop -45 45])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-45:45:90,'fontsize',8)
%     xlabel('time','fontsize',10)
ylabel('max deviation','fontsize',10) 
    grid on
subplot(3,4,7)
axis([t_start t_stop -45 45])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-45:45:45,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('min deviation','fontsize',10)
    grid on
subplot(3,4,11)
axis([t_start t_stop -45 45])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-45:45:45,'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('Adev diff L-R','fontsize',10)
    grid on

subplot(3,4,4)
axis([t_start t_stop -45 45])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-45:45:90,'fontsize',8)
%     xlabel('time','fontsize',10)
ylabel('max deviation','fontsize',10) 
    grid on
subplot(3,4,8)
axis([t_start t_stop -45 45])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-45:45:45,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('min deviation','fontsize',10)
    grid on
subplot(3,4,12)
axis([t_start t_stop -45 45])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-45:45:45,'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('Adev diff L-R','fontsize',10)
    grid on
    
