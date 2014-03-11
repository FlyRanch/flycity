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
            plot(t_wb_Fenhance_bins,strokeMOD_wb_L_Fenhance_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)
            plot(t_wb_Fenhance_bins,strokeMOD_wb_R_Fenhance_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)

            subplot(3,1,2)
            plot(t_wb_Fenhance_bins,pitchMOD_wb_L_Fenhance_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)
            plot(t_wb_Fenhance_bins,pitchMOD_wb_R_Fenhance_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)

            subplot(3,1,3)
            plot(t_wb_Fenhance_bins,devMOD_wb_L_Fenhance_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)
            plot(t_wb_Fenhance_bins,devMOD_wb_R_Fenhance_bins,'-','color',color_code_now,'linewidth',linewidth_timelines)

            
%% plot mean & std OR CI
binx = t_wb_Fenhance_bins(:,1);
t_bins = [t_wb_Fenhance_bins t_wb_Fenhance_bins];           

subplot(3,1,1)
var = [strokeMOD_wb_L_Fenhance_bins strokeMOD_wb_R_Fenhance_bins];
plot_WBmeanstd
% plot_WBmeanCI

subplot(3,1,2)
var = [pitchMOD_wb_L_Fenhance_bins pitchMOD_wb_R_Fenhance_bins];
plot_WBmeanstd
% plot_WBmeanCI

subplot(3,1,3)
var = [devMOD_wb_L_Fenhance_bins devMOD_wb_R_Fenhance_bins];
plot_WBmeanstd
% plot_WBmeanCI

%% plot polinomials
subplot(3,1,1)
plot_WBfitting_singlevar_updowncut(strokeMOD_Fenhance_fit_binmean_periodic,'r','g');

subplot(3,1,2)
plot_WBfitting_singlevar_updowncut(pitchMOD_Fenhance_fit_binmean_periodic,'r','g');

subplot(3,1,3)
plot_WBfitting_singlevar_updowncut(devMOD_Fenhance_fit_binmean_periodic,'r','g');



subplot(3,1,1)
axis([0 1 -30 30])
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
    set(gca,'YTick',-90:90:90,'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('Left wing','fontsize',10)
    grid on

subplot(3,1,2)
axis([0 1 -30 30])
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
