    subplot(2,2,1)
    hold on
    plot([0 1],[-2*mean(Fz_intactwing_steady(:,end))    -2*mean(Fz_intactwing_steady(:,end))],   '-','color',[.5 .5 .5],'linewidth',1)

    plot([0 1],[-mean(Fz_intactwing_all(:,end))    -mean(Fz_intactwing_all(:,end))],   '--','color',[0 0 1],'linewidth',2)
    plot([0 1],[-mean(Fz_damagedwing_all(:,end))   -mean(Fz_damagedwing_all(:,end))],  '--','color',[1 0 0],'linewidth',2)
    plot([0 1],[-mean(Fz_damagedwing_all(:,end))-mean(Fz_intactwing_all(:,end)) -mean(Fz_damagedwing_all(:,end))-mean(Fz_intactwing_all(:,end))],  '--','color',[0 .5 0],'linewidth',2)

    plot(t_norm,-Fz_intactwing_all(:,end),    '-','color',[0 0 1],'linewidth',2)
    plot(t_norm,-Fz_damagedwing_all(:,end),   '-','color',[1 0 0],'linewidth',2)
    plot(t_norm,-Fz_damagedwing_all(:,end)-Fz_intactwing_all(:,end),   '-','color',[0 .5 0],'linewidth',2)


    xlabel('wingbeat cycle')
    ylabel('normalized vertical force Fz/mg')
    axis([0 1 0 2]) 
    set(gca,'xtick',0:1)
    set(gca,'ytick',-4:1)

    subplot(2,2,2)
    hold on
    plot([0 1],[0 0],   '-','color',[.5 .5 .5],'linewidth',1)

    plot([0 1],[mean(Mx_intactwing_all(:,end))    mean(Mx_intactwing_all(:,end))],   '--','color',[0 0 1],'linewidth',2)
    plot([0 1],[mean(Mx_damagedwing_all(:,end))   mean(Mx_damagedwing_all(:,end))],  '--','color',[1 0 0],'linewidth',2)
    plot([0 1],[mean(Mx_damagedwing_all(:,end))+mean(Mx_intactwing_all(:,end)) mean(Mx_damagedwing_all(:,end))+mean(Mx_intactwing_all(:,end))],  '--','color',[0 .5 0],'linewidth',2)
    
    plot(t_norm,Mx_damagedwing_all(:,end),   '-','color',[1 0 0],'linewidth',2)
    plot(t_norm,Mx_intactwing_all(:,end),    '-','color',[0 0 1],'linewidth',2)
    plot(t_norm,Mx_damagedwing_all(:,end)+Mx_intactwing_all(:,end),   '-','color',[0 .5 0],'linewidth',2)

    xlabel('wingbeat cycle')
    ylabel('normalized roll torque Mx/mgl')
    axis([0 1 -1 1]) 
    set(gca,'xtick',0:1)
    set(gca,'ytick',-4:1)

