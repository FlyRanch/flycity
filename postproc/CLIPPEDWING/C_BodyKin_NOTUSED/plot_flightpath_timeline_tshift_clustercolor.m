%% plot flightpath timeline tshift clusters

close all
skip = 100;
for i=1:size(V,2)
    for j=skip+1:skip:size(V,1)
        if isnan(IDX(j,i)) == 0 && isnan(IDX(j-skip,i)) == 0
            subplot(3,2,1)
            hold on
            plot([t(j-skip)-t_shift(i) t(j)-t_shift(i)],[teta(j-skip,i) teta(j,i)],'-','color',cmap_k(IDX(j,i),:))
            subplot(3,2,3)
            hold on
            plot([t(j-skip)-t_shift(i) t(j)-t_shift(i)],[stim_angle_vel_plot(j-skip,i) stim_angle_vel_plot(j,i)],'-','color',cmap_k(IDX(j,i),:))
            subplot(3,2,5)
            hold on
            plot([t(j-skip)-t_shift(i) t(j)-t_shift(i)],[stim_angle_yaw_plot(j-skip,i) stim_angle_yaw_plot(j,i)],'-','color',cmap_k(IDX(j,i),:))
            subplot(3,2,2)
            hold on
            plot([t(j-skip)-t_shift(i) t(j)-t_shift(i)],[V(j-skip,i) V(j,i)],'-','color',cmap_k(IDX(j,i),:))
            subplot(3,2,4)
            hold on
            plot([t(j-skip)-t_shift(i) t(j)-t_shift(i)],[An_hor(j-skip,i) An_hor(j,i)],'-','color',cmap_k(IDX(j,i),:))
            subplot(3,2,6)
            hold on
            plot([t(j-skip)-t_shift(i) t(j)-t_shift(i)],[At_hor(j-skip,i) At_hor(j,i)],'-','color',cmap_k(IDX(j,i),:))
        end
    end
end



subplot(3,2,1)
axis([0-nanmean(t_shift) .15-nanmean(t_shift) 0 180])
    axis([0 .15 0 180])
    set(gca,'XTick',0:.05:.15) 
    set(gca,'YTick',0:180:180,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('optical angle','fontsize',18)
    grid on

subplot(3,2,3)
axis([0-nanmean(t_shift) .15-nanmean(t_shift) -180 180])
    set(gca,'XTick',0:.05:.15) 
    set(gca,'YTick',-180:180:180,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('heading','fontsize',18)
    grid on

subplot(3,2,5)
axis([0-nanmean(t_shift) .15-nanmean(t_shift) -180 180])
    set(gca,'XTick',0:.05:.15) 
    set(gca,'YTick',-180:180:180,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('yaw','fontsize',18)
    grid on


subplot(3,2,2)
axis([0-nanmean(t_shift) .15-nanmean(t_shift) 0 1])
    set(gca,'XTick',0:.05:.15) 
    set(gca,'YTick',0:.4:1,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('V','fontsize',18)
    grid on

subplot(3,2,4)
axis([0-nanmean(t_shift) .15-nanmean(t_shift) 0 15])
    set(gca,'XTick',0:.05:.15) 
    set(gca,'YTick',0:5:15,'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('An','fontsize',18)
    grid on

subplot(3,2,6)
axis([0-nanmean(t_shift) .15-nanmean(t_shift) -5 15])
    set(gca,'XTick',0:.05:.15) 
    set(gca,'YTick',[-5;0;15],'fontsize',12)
    xlabel('time','fontsize',18)
    ylabel('At','fontsize',18)
    grid on

