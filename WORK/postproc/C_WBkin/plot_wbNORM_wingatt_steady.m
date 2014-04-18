figure
subplot(3,1,1)
title('wing stroke')
hold on
subplot(3,1,2)
title('wing pitch')
hold on
subplot(3,1,3)
title('stroke deviation')
hold on

            
            subplot(3,1,1)
            plot(t_wb_steady_bins,stroke_wb_L_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)
            plot(t_wb_steady_bins,stroke_wb_R_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)

            subplot(3,1,2)
            plot(t_wb_steady_bins,pitch_wb_L_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)
            plot(t_wb_steady_bins,pitch_wb_R_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)

            subplot(3,1,3)
            plot(t_wb_steady_bins,dev_wb_L_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)
            plot(t_wb_steady_bins,dev_wb_R_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)

%% mean & std            
binx = t_wb_steady_bins(:,1);
t_bins = [t_wb_steady_bins t_wb_steady_bins];           


            subplot(3,1,1)
            var = [stroke_wb_L_steady_bins stroke_wb_R_steady_bins];
            plot_WBmeanstd
            
            subplot(3,1,2)
            var = [pitch_wb_L_steady_bins pitch_wb_R_steady_bins];
            plot_WBmeanstd

            subplot(3,1,3)
            var = [dev_wb_L_steady_bins dev_wb_R_steady_bins];
            plot_WBmeanstd

%% polynomials
subplot(3,1,1)
plot_WBfitting_singlevar_updowncut(stroke_steady_fit_binmean_periodic,'r','g');

subplot(3,1,2)
plot_WBfitting_singlevar_updowncut(pitch_steady_fit_binmean_periodic,'r','g');

subplot(3,1,3)
plot_WBfitting_singlevar_updowncut(dev_steady_fit_binmean_periodic,'r','g');



subplot(3,1,1)
axis([0 1 -90 90])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('Left wing','fontsize',10)
    grid on

subplot(3,1,2)
axis([0 1 0 180])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:180,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
    grid on
            
subplot(3,1,3)
axis([0 1 -30 30])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:30:180,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
    grid on