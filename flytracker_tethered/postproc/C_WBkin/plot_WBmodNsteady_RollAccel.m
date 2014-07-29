figure
subplot(3,3,1)
title('wing stroke')
hold on
subplot(3,3,2)
title('wing pitch')
hold on
subplot(3,3,3)
title('stroke deviation')
hold on
subplot(3,3,4)
hold on
subplot(3,3,5)
hold on
subplot(3,3,6)
hold on
subplot(3,3,7)
hold on
subplot(3,3,8)
hold on
subplot(3,3,9)
hold on

color_band = [.5 .5 .5];

%% steady WB
% plot mean & std OR CI
binx = t_wb_steady_bins(:,1);
t_bins = [t_wb_steady_bins t_wb_steady_bins];           


            subplot(3,3,1)
            var = [stroke_wb_L_steady_bins stroke_wb_R_steady_bins];
            plot_WBmeanstd
%             plot_WBmeanCI
            
            subplot(3,3,2)
            var = [pitch_wb_L_steady_bins pitch_wb_R_steady_bins];
            plot_WBmeanstd
%             plot_WBmeanCI

            subplot(3,3,3)
            var = [dev_wb_L_steady_bins dev_wb_R_steady_bins];
            plot_WBmeanstd
%             plot_WBmeanCI

            subplot(3,3,4)
            var = [stroke_wb_L_steady_bins stroke_wb_R_steady_bins];
            plot_WBmeanstd
%             plot_WBmeanCI
            
            subplot(3,3,5)
            var = [pitch_wb_L_steady_bins pitch_wb_R_steady_bins];
            plot_WBmeanstd
%             plot_WBmeanCI

            subplot(3,3,6)
            var = [dev_wb_L_steady_bins dev_wb_R_steady_bins];
            plot_WBmeanstd
%             plot_WBmeanCI


            subplot(3,3,7)
            var = [stroke_wb_L_steady_bins stroke_wb_R_steady_bins];
            plot_WBmeanstd
%             plot_WBmeanCI
            
            subplot(3,3,8)
            var = [pitch_wb_L_steady_bins pitch_wb_R_steady_bins];
            plot_WBmeanstd
%             plot_WBmeanCI

            subplot(3,3,9)
            var = [dev_wb_L_steady_bins dev_wb_R_steady_bins];
            plot_WBmeanstd
%             plot_WBmeanCI


% fourier series
subplot(3,3,1)
stroke_steady = calc_val_fourier_series_4thN8th_order(binx,stroke_steady_fourier_coeffs_binmean,0);
plot(binx,stroke_steady,'--k','linewidth',2)
% plot(stroke_steady_fourier_fit_binmean);

subplot(3,3,2)
pitch_steady = calc_val_fourier_series_4thN8th_order(binx,pitch_steady_fourier_coeffs_binmean,0);
plot(binx,pitch_steady,'--k','linewidth',2)
% plot(pitch_steady_fourier_fit_binmean);

subplot(3,3,3)
dev_steady = calc_val_fourier_series_4thN8th_order(binx,dev_steady_fourier_coeffs_binmean,0);
plot(binx,dev_steady,'--k','linewidth',2)
% plot(dev_steady_fourier_fit_binmean);

subplot(3,3,4)
stroke_steady = calc_val_fourier_series_4thN8th_order(binx,stroke_steady_fourier_coeffs_binmean,0);
plot(binx,stroke_steady,'--k','linewidth',2)
% plot(stroke_steady_fourier_fit_binmean);

subplot(3,3,5)
pitch_steady = calc_val_fourier_series_4thN8th_order(binx,pitch_steady_fourier_coeffs_binmean,0);
plot(binx,pitch_steady,'--k','linewidth',2)
% plot(pitch_steady_fourier_fit_binmean);

subplot(3,3,6)
dev_steady = calc_val_fourier_series_4thN8th_order(binx,dev_steady_fourier_coeffs_binmean,0);
plot(binx,dev_steady,'--k','linewidth',2)
% plot(dev_steady_fourier_fit_binmean);

subplot(3,3,7)
stroke_steady = calc_val_fourier_series_4thN8th_order(binx,stroke_steady_fourier_coeffs_binmean,0);
plot(binx,stroke_steady,'--k','linewidth',2)
% plot(stroke_steady_fourier_fit_binmean);

subplot(3,3,8)
pitch_steady = calc_val_fourier_series_4thN8th_order(binx,pitch_steady_fourier_coeffs_binmean,0);
plot(binx,pitch_steady,'--k','linewidth',2)
% plot(pitch_steady_fourier_fit_binmean);

subplot(3,3,9)
dev_steady = calc_val_fourier_series_4thN8th_order(binx,dev_steady_fourier_coeffs_binmean,0);
plot(binx,dev_steady,'--k','linewidth',2)
% plot(dev_steady_fourier_fit_binmean);



%% MOD wb

%% fourier series

subplot(3,3,1)
strokeMOD_up = calc_val_fourier_series_4thN8th_order(binx,strokeMOD_up_RollAccel_fourier_coeffs_binmean,0);
plot(binx,stroke_steady+strokeMOD_up,'-r','linewidth',2)

legend off
xlabel([],'fontsize',10) 
ylabel('upward wing','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 

subplot(3,3,2)
pitchMOD_up = calc_val_fourier_series_4thN8th_order(binx,pitchMOD_up_RollAccel_fourier_coeffs_binmean,0);
plot(binx,pitch_steady+pitchMOD_up,'-r','linewidth',2)

legend off
xlabel([],'fontsize',10) 
ylabel([],'fontsize',10) 
biny_min = 0;
biny_max = 180;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 

subplot(3,3,3)
devMOD_up = calc_val_fourier_series_4thN8th_order(binx,devMOD_up_RollAccel_fourier_coeffs_binmean,0);
plot(binx,dev_steady+devMOD_up,'-r','linewidth',2)

legend off
xlabel([],'fontsize',10) 
ylabel([],'fontsize',10) 
biny_min = -30;
biny_max = 30;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 

subplot(3,3,4)
strokeMOD_down = calc_val_fourier_series_4thN8th_order(binx,strokeMOD_down_RollAccel_fourier_coeffs_binmean,0);
plot(binx,stroke_steady+strokeMOD_down,'-r','linewidth',2)

legend off
xlabel([],'fontsize',10) 
ylabel('downward wing','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 

subplot(3,3,5)
pitchMOD_down = calc_val_fourier_series_4thN8th_order(binx,pitchMOD_down_RollAccel_fourier_coeffs_binmean,0);
plot(binx,pitch_steady+pitchMOD_down,'-r','linewidth',2)

legend off
xlabel([],'fontsize',10) 
ylabel([],'fontsize',10) 
biny_min = 0;
biny_max = 180;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 

subplot(3,3,6)
devMOD_down = calc_val_fourier_series_4thN8th_order(binx,devMOD_down_RollAccel_fourier_coeffs_binmean,0);
plot(binx,dev_steady+devMOD_down,'-r','linewidth',2)

legend off
xlabel([],'fontsize',10) 
ylabel([],'fontsize',10) 
biny_min = -30;
biny_max = 30;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 

subplot(3,3,7)
DstrokeMOD = calc_val_fourier_series_4thN8th_order(binx,DstrokeMOD_RollAccel_fourier_coeffs_binmean,0);
plot(binx,stroke_steady+DstrokeMOD,'-r','linewidth',2)

legend off
xlabel('time','fontsize',10) 
ylabel('upward - downward','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(3,3,8)
DpitchMOD = calc_val_fourier_series_4thN8th_order(binx,DpitchMOD_RollAccel_fourier_coeffs_binmean,0);
plot(binx,pitch_steady+DpitchMOD,'-r','linewidth',2)

legend off
xlabel('time','fontsize',10) 
ylabel([],'fontsize',10) 
biny_min = 0;
biny_max = 180;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(3,3,9)
DdevMOD = calc_val_fourier_series_4thN8th_order(binx,DdevMOD_RollAccel_fourier_coeffs_binmean,0);
plot(binx,dev_steady+DdevMOD,'-r','linewidth',2)

legend off
xlabel('time','fontsize',10) 
ylabel([],'fontsize',10) 
biny_min = -30;
biny_max = 30;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 







