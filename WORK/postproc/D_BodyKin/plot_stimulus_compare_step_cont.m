figure
plot(step_fast.t_cam,step_fast.teta_cam,'r','linewidth',2)
hold on
plot(cont_fast_165.t_cam,cont_fast_165.teta_cam,'b','linewidth',2)
plot(cont_fast_64_off.t_cam,cont_fast_64_off.teta_cam,'b','linewidth',2)



xlabel('time','fontsize',18)
ylabel('optical angle','fontsize',18)
set(gca,'xlim',[0 .15],'ylim',[0 180])
set(gca,'XTick',[0:.05:.15],'fontsize',12)
set(gca,'YTick',[0:90:180],'fontsize',12)


saveas(gca,['expansion_step_vs_cont.fig'])
saveas(gca,['expansion_step_vs_cont.png'])
plot2svg(['expansion_step_vs_cont.svg'])

