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
        plot(t-t_shift(i),stim_angle_yaw_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,2,5)
        plot(t-t_shift(i),slip_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,2,2)
        plot(t-t_shift(i),dV_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,2,4)
        plot(t-t_shift(i),An_hor_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,2,6)
        plot(t-t_shift(i),At_hor_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
    end
end

subplot(3,2,1)
axis([-.025 .05 -180 180])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-180:180:180,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('heading','fontsize',18)
    grid on

subplot(3,2,3)
axis([-.025 .05 -180 180])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-180:180:180,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('yaw','fontsize',18)
    grid on

subplot(3,2,5)
axis([-.025 .05 -90 90])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-90:90:90,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('slip','fontsize',18)
    grid on


subplot(3,2,2)
axis([-.025 .05 -.25 .5])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-.25:.25:1,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('dV','fontsize',18)
    grid on

subplot(3,2,4)
axis([-.025 .05 -5 20])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',[-5;0;20],'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('An','fontsize',18)
    grid on

subplot(3,2,6)
axis([-.025 .05 -15 20])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',[-15;0;20],'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('At','fontsize',18)
    grid on

