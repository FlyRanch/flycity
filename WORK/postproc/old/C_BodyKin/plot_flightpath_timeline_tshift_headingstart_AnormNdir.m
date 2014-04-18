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
        plot(t-t_shift(i),stim_angle_vel_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,2,3)
        plot(t-t_shift(i),dV_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,2,5)
        plot(t-t_shift(i),stim_angle_accel_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        
        subplot(3,2,2)
        plot(t-t_shift(i),A_hor_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,2,4)
        plot(t-t_shift(i),An_hor_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,2,6)
        plot(t-t_shift(i),At_hor_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
    end
end

subplot(3,2,1)
axis([t_start t_stop -180 180])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-180:180:180,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('heading','fontsize',10)
    grid on
subplot(3,2,3)
axis([t_start t_stop -.25 .5])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-.25:.25:1,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('dV','fontsize',10)
    grid on
subplot(3,2,5)
axis([t_start t_stop -180 180])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-180:180:180,'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('A direction','fontsize',10)
    grid on

subplot(3,2,2)
axis([t_start t_stop 0 20])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',[0:10:20],'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('A','fontsize',10)
    grid on
subplot(3,2,4)
axis([t_start t_stop -20 20])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',[-20:20:20],'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('An','fontsize',10)
    grid on
subplot(3,2,6)
axis([t_start t_stop -20 20])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',[-20:20:20],'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('At','fontsize',10)
    grid on

