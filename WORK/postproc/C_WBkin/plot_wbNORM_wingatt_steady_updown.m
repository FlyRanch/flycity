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
%             plot(t_steady_bins,stroke_wb_L_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)
%             plot(t_steady_bins,stroke_wb_R_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)
            plot(t_ds_L_steady,stroke_ds_L_steady,'-','color',color_code_now,'linewidth',linewidth_timelines)
            plot(t_ds_R_steady,stroke_ds_R_steady,'-','color',color_code_now,'linewidth',linewidth_timelines)

            plot(t_us_L_steady,stroke_us_L_steady,'-','color',color_code_now,'linewidth',linewidth_timelines)
            plot(t_us_R_steady,stroke_us_R_steady,'-','color',color_code_now,'linewidth',linewidth_timelines)

            subplot(3,1,2)
%             plot(t_steady_bins,pitch_wb_L_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)
%             plot(t_steady_bins,pitch_wb_R_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)

            plot(t_ds_L_steady,pitch_ds_L_steady,'-','color',color_code_now,'linewidth',linewidth_timelines)
            plot(t_ds_R_steady,pitch_ds_R_steady,'-','color',color_code_now,'linewidth',linewidth_timelines)

            plot(t_us_L_steady,pitch_us_L_steady,'-','color',color_code_now,'linewidth',linewidth_timelines)
            plot(t_us_R_steady,pitch_us_R_steady,'-','color',color_code_now,'linewidth',linewidth_timelines)

            subplot(3,1,3)
%             plot(t_steady_bins,dev_wb_L_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)
%             plot(t_steady_bins,dev_wb_R_steady_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)

            plot(t_ds_L_steady,dev_ds_L_steady,'-','color',color_code_now,'linewidth',linewidth_timelines)
            plot(t_ds_R_steady,dev_ds_R_steady,'-','color',color_code_now,'linewidth',linewidth_timelines)

            plot(t_us_L_steady,dev_us_L_steady,'-','color',color_code_now,'linewidth',linewidth_timelines)
            plot(t_us_R_steady,dev_us_R_steady,'-','color',color_code_now,'linewidth',linewidth_timelines)

            
% plot_WBfunc_csaps_LRdiff
% plot_WBfunc_circmeanCI


            subplot(3,1,1)
            t_bins = [t_ds_steady_bins t_ds_steady_bins];           
            var = [stroke_ds_L_steady_bins stroke_ds_R_steady_bins];
            binx = t_bins(:,1);
            plot_WBmeanstd
            
            t_bins = [t_us_steady_bins t_us_steady_bins];           
            var = [stroke_us_L_steady_bins stroke_us_R_steady_bins];
            binx = t_bins(:,1);
            plot_WBmeanstd
            
            subplot(3,1,2)
            t_bins = [t_ds_steady_bins t_ds_steady_bins];           
            var = [pitch_ds_L_steady_bins pitch_ds_R_steady_bins];
            binx = t_bins(:,1);
            plot_WBmeanstd
            
            t_bins = [t_us_steady_bins t_us_steady_bins];           
            var = [pitch_us_L_steady_bins pitch_us_R_steady_bins];
            binx = t_bins(:,1);
            plot_WBmeanstd

            subplot(3,1,3)
            t_bins = [t_ds_steady_bins t_ds_steady_bins];           
            var = [dev_ds_L_steady_bins dev_ds_R_steady_bins];
            binx = t_bins(:,1);
            plot_WBmeanstd
            
            t_bins = [t_us_steady_bins t_us_steady_bins];           
            var = [dev_us_L_steady_bins dev_us_R_steady_bins];
            binx = t_bins(:,1);
            plot_WBmeanstd




% subplot(3,1,1)
% plot_WBfitting_singlevar(stroke_steady_fit_binmean_periodic,'r');
% 
% subplot(3,1,2)
% plot_WBfitting_singlevar(pitch_steady_fit_binmean_periodic,'r');
% 
% subplot(3,1,3)
% plot_WBfitting_singlevar(dev_steady_fit_binmean_periodic,'r');



subplot(3,1,1)
axis([-1 1 -90 90])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('Left wing','fontsize',10)
    grid on

subplot(3,1,2)
axis([-1 1 0 180])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:180,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
    grid on
            
subplot(3,1,3)
axis([-1 1 -30 30])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:30:180,'fontsize',8)
%     xlabel('time','fontsize',10)
%     ylabel('Left wing','fontsize',10)
    grid on
