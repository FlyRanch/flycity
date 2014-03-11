n=40;

figure
subplot(3,5,1)
title('wing stroke')
hold on
subplot(3,5,2)
title('wing pitch')
hold on
subplot(3,5,3)
title('stroke deviation')
hold on
subplot(3,5,4)
title('angle of attack')
hold on
subplot(3,5,5)
title('wing speed')
hold on
subplot(3,5,6)
hold on
subplot(3,5,7)
hold on
subplot(3,5,8)
hold on
subplot(3,5,9)
hold on
subplot(3,5,10)
hold on
subplot(3,5,11)
hold on
subplot(3,5,12)
hold on
subplot(3,5,13)
hold on
subplot(3,5,14)
hold on
subplot(3,5,15)
hold on

            
            subplot(3,5,1)
            plot(t_steady_bins,stroke_wb_L_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)
            subplot(3,5,6)
            plot(t_steady_bins,stroke_wb_R_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)
            subplot(3,5,11)
            plot(t_steady_bins,Dstroke_wb_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)

            subplot(3,5,2)
            plot(t_steady_bins,pitch_wb_L_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)
            subplot(3,5,7)
            plot(t_steady_bins,pitch_wb_R_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)
            subplot(3,5,12)
            plot(t_steady_bins,Dpitch_wb_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)

            subplot(3,5,3)
            plot(t_steady_bins,dev_wb_L_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)
            subplot(3,5,8)
            plot(t_steady_bins,dev_wb_R_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)
            subplot(3,5,13)
            plot(t_steady_bins,Ddev_wb_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)

             subplot(3,5,4)
            plot(t_steady_bins,aoa_wb_L_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)
            subplot(3,5,9)
            plot(t_steady_bins,aoa_wb_R_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)
            subplot(3,5,14)
            plot(t_steady_bins,Daoa_wb_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)

            subplot(3,5,5)
            plot(t_steady_bins,U_wb_L_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)
            subplot(3,5,10)
            plot(t_steady_bins,U_wb_R_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)
            subplot(3,5,15)
            plot(t_steady_bins,DU_wb_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)

            
% plot_WBfunc_csaps_LRdiff
plot_WBfunc_circmeanCI


subplot(3,5,1)
axis([0 1 -90 90])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('Left wing','fontsize',10)
    grid on
subplot(3,5,6)
axis([0 1 -90 90])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('Right wing','fontsize',10)
    grid on
subplot(3,5,11)
axis([0 1 -45 45])
    set(gca,'XTick',0:.5:1) 
    set(gca,'YTick',-90:45:90,'fontsize',8)
    xlabel('normalized time','fontsize',10)
    ylabel('Left - Right','fontsize',10)
    grid on
            

subplot(3,5,2)
axis([0 1 0 180])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:180,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
    grid on
subplot(3,5,7)
axis([0 1 0 180])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:180,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
    grid on
subplot(3,5,12)
axis([0 1 -90 90])
    set(gca,'XTick',0:.5:1) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
    xlabel('normalized time','fontsize',10)
%     ylabel('Left - Right','fontsize',10)
    grid on
            
subplot(3,5,3)
axis([0 1 -30 30])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:30:180,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
    grid on
subplot(3,5,8)
axis([0 1 -30 30])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:30:180,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
    grid on
subplot(3,5,13)
axis([0 1 -45 45])
    set(gca,'XTick',0:.5:1) 
    set(gca,'YTick',-90:45:90,'fontsize',8)
    xlabel('normalized time','fontsize',10)
%     ylabel('Left - Right','fontsize',10)
    grid on
            
    

subplot(3,5,4)
axis([0 1 -90 90])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
    grid on
subplot(3,5,9)
axis([0 1 -90 90])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Right wing','fontsize',10)
    grid on
subplot(3,5,14)
axis([0 1 -90 90])
    set(gca,'XTick',0:.5:1) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
    xlabel('normalized time','fontsize',10)
%     ylabel('Left - Right','fontsize',10)
    grid on
            
subplot(3,5,5)
axis([0 1 0 6])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',0:3:6,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
    grid on
subplot(3,5,10)
axis([0 1 0 6])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',0:3:6,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
    grid on
subplot(3,5,15)
axis([0 1 -3 3])
    set(gca,'XTick',0:.5:1) 
    set(gca,'YTick',-3:3:3,'fontsize',8)
    xlabel('normalized time','fontsize',10)
%     ylabel('Left - Right','fontsize',10)
    grid on
            


