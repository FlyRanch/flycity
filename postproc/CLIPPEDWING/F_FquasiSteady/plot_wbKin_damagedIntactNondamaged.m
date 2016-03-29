    subplot(2,2,1)
    hold on
    plot(t_norm,stroke_damaged_all(:,1),'-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,stroke_intact_all(:,end),'-','color',[0 0 1],'linewidth',2)
    plot(t_norm,stroke_damaged_all(:,end),'-','color',[1 0 0],'linewidth',2)
    
    xlabel('wingbeat cycle')
    ylabel('stroke angle')
    axis([0 1 -90 90])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-90:90:90)    
    
    subplot(2,2,2)
    hold on
    plot(t_norm,dev_damaged_all(:,1),'-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,dev_intact_all(:,end),'-','color',[0 0 1],'linewidth',2)
    plot(t_norm,dev_damaged_all(:,end),'-','color',[1 0 0],'linewidth',2)
    
    xlabel('wingbeat cycle')
    ylabel('deviation angle')
    axis([0 1 -90 90])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-90:90:90)    
    
    subplot(2,2,3)
    hold on
    plot(t_norm,rot_damaged_all(:,1)-90,'-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,rot_intact_all(:,end)-90,'-','color',[0 0 1],'linewidth',2)
    plot(t_norm,rot_damaged_all(:,end)-90,'-','color',[1 0 0],'linewidth',2)
    
    xlabel('wingbeat cycle')
    ylabel('rotation angle')
    axis([0 1 -90 90])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-90:90:90)    
