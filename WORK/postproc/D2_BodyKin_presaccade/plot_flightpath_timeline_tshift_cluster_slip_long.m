%% plot flightpath timeline tshift clusters

close all
subplot(6,1,1)
hold on
subplot(6,1,3)
hold on
subplot(6,1,5)
hold on
subplot(6,1,2)
hold on
subplot(6,1,4)
hold on
subplot(6,1,6)
hold on

for i=1:size(stim_angle_vel_plot,2)
    size(stim_angle_vel_plot,2) - i
    for j=2:size(stim_angle_vel_plot,1)
        if isnan(IDX_plot(j,i)) == 0 && isnan(IDX_plot(j-1,i)) == 0
            subplot(6,1,1)
            plot([t_plot(j-1)-t_shift(i) t_plot(j)-t_shift(i)],[stim_angle_vel_plot(j-1,i) stim_angle_vel_plot(j,i)],'-','color',cmap_k(IDX_plot(j,i),:))
            subplot(6,1,3)
            plot([t_plot(j-1)-t_shift(i) t_plot(j)-t_shift(i)],[stim_angle_yaw_plot(j-1,i) stim_angle_yaw_plot(j,i)],'-','color',cmap_k(IDX_plot(j,i),:))
            subplot(6,1,5)
            plot([t_plot(j-1)-t_shift(i) t_plot(j)-t_shift(i)],[slip_plot(j-1,i) slip_plot(j,i)],'-','color',cmap_k(IDX_plot(j,i),:))
            subplot(6,1,2)
            plot([t_plot(j-1)-t_shift(i) t_plot(j)-t_shift(i)],[V_plot(j-1,i) V_plot(j,i)],'-','color',cmap_k(IDX_plot(j,i),:))
            subplot(6,1,4)
            plot([t_plot(j-1)-t_shift(i) t_plot(j)-t_shift(i)],[An_hor_plot(j-1,i) An_hor_plot(j,i)],'-','color',cmap_k(IDX_plot(j,i),:))
            subplot(6,1,6)
            plot([t_plot(j-1)-t_shift(i) t_plot(j)-t_shift(i)],[At_hor_plot(j-1,i) At_hor_plot(j,i)],'-','color',cmap_k(IDX_plot(j,i),:))
        end
    end
end

subplot(6,1,1)
axis([-.3 .15 -180 180])
    set(gca,'XTick',-0.5:.15:.5) 
    set(gca,'YTick',-180:180:180,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('heading','fontsize',18)
    grid on

subplot(6,1,3)
axis([-.3 .15 -180 180])
    set(gca,'XTick',-0.5:.15:.5) 
    set(gca,'YTick',-180:180:180,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('yaw','fontsize',18)
    grid on

subplot(6,1,5)
axis([-.3 .15 -90 90])
    set(gca,'XTick',-0.5:.15:.5) 
    set(gca,'YTick',-90:90:90,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('slip','fontsize',18)
    grid on


subplot(6,1,2)
axis([-.3 .15 0 .8])
    set(gca,'XTick',-0.5:.15:.5) 
    set(gca,'YTick',0:.4:1,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('V','fontsize',18)
    grid on

subplot(6,1,4)
axis([-.3 .15 0 15])
    set(gca,'XTick',-0.5:.15:.5) 
    set(gca,'YTick',0:5:15,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('An','fontsize',18)
    grid on

subplot(6,1,6)
axis([-.3 .15 -5 15])
    set(gca,'XTick',-0.5:.15:.5) 
    set(gca,'YTick',[-5;0;15],'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('At','fontsize',18)
    grid on

