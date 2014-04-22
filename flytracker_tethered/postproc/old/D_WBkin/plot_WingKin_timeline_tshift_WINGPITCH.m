%% plot flightpath timeline tshift headingstart

% close all
figure
subplot(3,3,1)
title('wingpitch max')
hold on
subplot(3,3,2)
title('wingpitch min')
hold on
subplot(3,3,3)
title('A wingpitch')
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
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(pitch_max_wb_L(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,3,2)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(pitch_min_wb_L(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,3,3)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(Apitch_wb_L(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        
        subplot(3,3,4)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(pitch_max_wb_R(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,3,5)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(pitch_min_wb_R(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,3,6)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(Apitch_wb_R(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        
        subplot(3,3,7)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(dpitch_max_wb(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,3,8)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(dpitch_min_wb(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,3,9)
        plot(t_wb_L(:,i)-t_shift(i),rad2deg(dApitch_wb(:,i)),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
    end
end

subplot(3,3,1)
axis([t_start t_stop 90 225])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-45:45:270,'fontsize',8)
%     xlabel('time','fontsize',10)
ylabel('left wing','fontsize',10) 
    grid on
subplot(3,3,2)
axis([t_start t_stop -45 90])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-45:45:270,'fontsize',8)
%     xlabel('time','fontsize',10)
% ylabel('left wing','fontsize',10) 
    grid on
subplot(3,3,3)
axis([t_start t_stop 90 225])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-45:45:270,'fontsize',8)
%     xlabel('time','fontsize',10)
% ylabel('left wing','fontsize',10) 
    grid on
    
subplot(3,3,4)
axis([t_start t_stop 90 225])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-45:45:270,'fontsize',8)
%     xlabel('time','fontsize',10)
ylabel('right wing','fontsize',10) 
    grid on
subplot(3,3,5)
axis([t_start t_stop -45 90])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-45:45:270,'fontsize',8)
%     xlabel('time','fontsize',10)
% ylabel('right wing','fontsize',10) 
    grid on
subplot(3,3,6)
axis([t_start t_stop 90 225])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-45:45:270,'fontsize',8)
%     xlabel('time','fontsize',10)
% ylabel('right wing','fontsize',10) 
    grid on
    
    
subplot(3,3,7)
axis([t_start t_stop -90 90])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-90:45:90,'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('pitch diff L-R','fontsize',10)
    grid on
subplot(3,3,8)
axis([t_start t_stop -90 90])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-90:45:90,'fontsize',8)
    xlabel('time','fontsize',10)
%     ylabel('pitch diff L-R','fontsize',10)
    grid on
subplot(3,3,9)
axis([t_start t_stop -90 90])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-90:45:90,'fontsize',8)
    xlabel('time','fontsize',10)
%     ylabel('pitch diff L-R','fontsize',10)
    grid on
