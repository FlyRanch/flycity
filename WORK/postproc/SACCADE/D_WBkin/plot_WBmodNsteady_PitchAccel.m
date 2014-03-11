figure
subplot(3,1,1)
title('Steady & Pitch Accel wingbeat','fontsize',12)
hold on
subplot(3,1,2)
% title('wing pitch')
hold on
subplot(3,1,3)
% title('stroke deviation')
hold on

color_band = [.5 .5 .5];

%% steady WB
% mean & std            
binx = t_wb_steady_bins(:,1);
t_bins = [t_wb_steady_bins t_wb_steady_bins];           


            subplot(3,1,1)
            var = [stroke_wb_L_steady_bins stroke_wb_R_steady_bins];
            plot_WBmeanstd
%             plot_WBmeanCI
            
            subplot(3,1,2)
            var = [pitch_wb_L_steady_bins pitch_wb_R_steady_bins];
            plot_WBmeanstd
%             plot_WBmeanCI

            subplot(3,1,3)
            var = [dev_wb_L_steady_bins dev_wb_R_steady_bins];
            plot_WBmeanstd
%             plot_WBmeanCI

% fourier series
subplot(3,1,1)
stroke_steady = calc_val_fourier_series_4thN8th_order(binx,stroke_steady_fourier_coeffs_binmean,0);
plot(binx,stroke_steady,'--k','linewidth',2)
% plot(stroke_steady_fourier_fit_binmean);
legend off
xlabel([],'fontsize',10) 
ylabel('wing stroke','fontsize',10) 

subplot(3,1,2)
pitch_steady = calc_val_fourier_series_4thN8th_order(binx,pitch_steady_fourier_coeffs_binmean,0);
plot(binx,pitch_steady,'--k','linewidth',2)
% plot(pitch_steady_fourier_fit_binmean);
legend off
xlabel([],'fontsize',10) 
ylabel('wing pitch','fontsize',10) 

subplot(3,1,3)
dev_steady = calc_val_fourier_series_4thN8th_order(binx,dev_steady_fourier_coeffs_binmean,0);
plot(binx,dev_steady,'--k','linewidth',2)
% plot(dev_steady_fourier_fit_binmean);
legend off
xlabel('normalized time','fontsize',10) 
ylabel('wing deviation','fontsize',10) 



%% MOD wb

subplot(3,1,1)
strokeMOD = calc_val_fourier_series_4thN8th_order(binx,strokeMOD_PitchAccel_fourier_coeffs_binmean,0);
plot(binx,stroke_steady+strokeMOD,'-r','linewidth',2)
plot(binx,stroke_steady-strokeMOD,'-b','linewidth',2)
% plot(strokeMOD_PitchAccel_fourier_fit_binmean);

legend off
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
xlabel([],'fontsize',10) 
ylabel('wing stroke','fontsize',10) 

subplot(3,1,2)
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,pitchMOD_PitchAccel_fourier_coeffs_binmean,0);
plot(binx,pitch_steady+pitchMOD,'-r','linewidth',2)
plot(binx,pitch_steady-pitchMOD,'-b','linewidth',2)
% plot(pitchMOD_PitchAccel_fourier_fit_binmean);

legend off
biny_min = 0;
biny_max = 180;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
xlabel([],'fontsize',10) 
ylabel('wing pitch','fontsize',10) 

subplot(3,1,3)
devMOD = calc_val_fourier_series_4thN8th_order(binx,devMOD_PitchAccel_fourier_coeffs_binmean,0);
plot(binx,dev_steady+devMOD,'-r','linewidth',2)
plot(binx,dev_steady-devMOD,'-b','linewidth',2)
% plot(devMOD_PitchAccel_fourier_fit_binmean);

legend off
biny_min = -30;
biny_max = 30;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 
xlabel('normalized time','fontsize',10) 
ylabel('wing deviation','fontsize',10) 





