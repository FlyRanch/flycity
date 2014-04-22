%% plot flightpath timeline tshift clusters

% close all
figure
subplot(3,3,1)
hold on
subplot(3,3,3)
hold on
subplot(3,3,5)
hold on
subplot(3,3,2)
hold on
subplot(3,3,4)
hold on
subplot(3,3,6)
hold on
subplot(3,3,7)
hold on
subplot(3,3,8)
hold on
subplot(3,3,9)
hold on

for i=1:size(stim_angle_vel_plot,2)
    size(stim_angle_vel_plot,2) - i
    for j=2:size(stim_angle_vel_plot,1)
        if isnan(IDX_plot(j,i)) == 0 && isnan(IDX_plot(j-1,i)) == 0
            
            subplot(3,3,1)
            plot([t_plot(j-1)-t_shift(i) t_plot(j)-t_shift(i)],[stim_angle_vel_plot(j-1,i) stim_angle_vel_plot(j,i)],'-','color',cmap_k(IDX_plot(j,i),:),'linewidth',linewidth_timelines)
            subplot(3,3,4)
            plot([t_plot(j-1)-t_shift(i) t_plot(j)-t_shift(i)],[stim_angle_yaw_plot(j-1,i) stim_angle_yaw_plot(j,i)],'-','color',cmap_k(IDX_plot(j,i),:),'linewidth',linewidth_timelines)
            subplot(3,3,7)
            plot([t_plot(j-1)-t_shift(i) t_plot(j)-t_shift(i)],[stim_angle_accel_plot(j-1,i) stim_angle_accel_plot(j,i)],'-','color',cmap_k(IDX_plot(j,i),:),'linewidth',linewidth_timelines)

            subplot(3,3,2)
            plot([t_plot(j-1)-t_shift(i) t_plot(j)-t_shift(i)],[roll_plot(j-1,i) roll_plot(j,i)],'-','color',cmap_k(IDX_plot(j,i),:),'linewidth',linewidth_timelines)
            subplot(3,3,5)
            plot([t_plot(j-1)-t_shift(i) t_plot(j)-t_shift(i)],[slip_plot(j-1,i) slip_plot(j,i)],'-','color',cmap_k(IDX_plot(j,i),:),'linewidth',linewidth_timelines)
            subplot(3,3,8)
            plot([t_plot(j-1)-t_shift(i) t_plot(j)-t_shift(i)],[pitch_plot(j-1,i) pitch_plot(j,i)],'-','color',cmap_k(IDX_plot(j,i),:),'linewidth',linewidth_timelines)
            
            subplot(3,3,3)
            plot([t_plot(j-1)-t_shift(i) t_plot(j)-t_shift(i)],[dV_plot(j-1,i) dV_plot(j,i)],'-','color',cmap_k(IDX_plot(j,i),:),'linewidth',linewidth_timelines)
            subplot(3,3,6)
            plot([t_plot(j-1)-t_shift(i) t_plot(j)-t_shift(i)],[An_hor_plot(j-1,i) An_hor_plot(j,i)],'-','color',cmap_k(IDX_plot(j,i),:),'linewidth',linewidth_timelines)
            subplot(3,3,9)
            plot([t_plot(j-1)-t_shift(i) t_plot(j)-t_shift(i)],[At_hor_plot(j-1,i) At_hor_plot(j,i)],'-','color',cmap_k(IDX_plot(j,i),:),'linewidth',linewidth_timelines)
        end
    end
end

subplot(3,3,1)
axis([-.025 .05 -180 180])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-180:180:180,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('heading','fontsize',18)
    grid on
subplot(3,3,4)
axis([-.025 .05 -180 180])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-180:180:180,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('yaw','fontsize',18)
    grid on
subplot(3,3,7)
axis([-.025 .05 -180 180])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-180:180:180,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('A direction','fontsize',18)
    grid on

subplot(3,3,5)
axis([-.025 .05 -90 90])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-90:90:90,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('slip','fontsize',18)
    grid on
subplot(3,3,5)
axis([-.025 .05 -90 90])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-90:90:90,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('roll','fontsize',18)
    grid on
subplot(3,3,8)
axis([-.025 .05 -90 90])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-90:90:90,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('pitch','fontsize',18)
    grid on


subplot(3,3,3)
axis([-.025 .05 -.25 .5])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-.25:.25:1,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('dV','fontsize',18)
    grid on
subplot(3,3,6)
axis([-.025 .05 -5 20])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',[-5;0;20],'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('An','fontsize',18)
    grid on
subplot(3,3,9)
axis([-.025 .05 -15 20])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',[-15;0;20],'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('At','fontsize',18)
    grid on

