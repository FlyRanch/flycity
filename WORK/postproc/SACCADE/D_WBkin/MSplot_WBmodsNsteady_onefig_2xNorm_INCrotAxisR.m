figure
subplot(1,3,1)
% title('wing stroke')
hold on
subplot(1,3,2)
% title('wing pitch')
hold on
subplot(1,3,3)
% title('stroke deviation')
hold on

color_band = [.5 .5 .5];

%% steady WB

subplot(1,3,1)
stroke_steady = calc_val_fourier_series_4thN8th_order(binx,stroke_steady_fourier_coeffs_binmean,0);
plot(binx,stroke_steady,'--k','linewidth',2)
% plot(stroke_steady_fourier_fit_binmean);

subplot(1,3,2)
pitch_steady = calc_val_fourier_series_4thN8th_order(binx,pitch_steady_fourier_coeffs_binmean,0);
plot(binx,pitch_steady,'--k','linewidth',2)
% plot(pitch_steady_fourier_fit_binmean);

subplot(1,3,3)
dev_steady = calc_val_fourier_series_4thN8th_order(binx,dev_steady_fourier_coeffs_binmean,0);
plot(binx,dev_steady,'--k','linewidth',2)
% plot(dev_steady_fourier_fit_binmean);

%% MOD wb

% Force
subplot(1,3,1)
strokeMOD = calc_val_fourier_series_4thN8th_order(binx,strokeMOD_Fenhance_fourier_coeffs_binmean,0);
plot(binx,stroke_steady+2*strokeMOD,'-','color',[1,.5,0],'linewidth',2)

subplot(1,3,2)
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,pitchMOD_Fenhance_fourier_coeffs_binmean,0);
plot(binx,pitch_steady+2*pitchMOD,'-','color',[1,.5,0],'linewidth',2)

subplot(1,3,3)
devMOD = calc_val_fourier_series_4thN8th_order(binx,devMOD_Fenhance_fourier_coeffs_binmean,0);
plot(binx,dev_steady+2*devMOD,'-','color',[1,.5,0],'linewidth',2)

% ROLL up-down
subplot(1,3,1)
DstrokeMOD = calc_val_fourier_series_4thN8th_order(binx,DstrokeMOD_RollAccel_fourier_coeffs_binmean,0);
plot(binx,stroke_steady+2*DstrokeMOD,'-g','linewidth',2)

subplot(1,3,2)
DpitchMOD = calc_val_fourier_series_4thN8th_order(binx,DpitchMOD_RollAccel_fourier_coeffs_binmean,0);
plot(binx,pitch_steady+2*DpitchMOD,'-g','linewidth',2)

subplot(1,3,3)
DdevMOD = calc_val_fourier_series_4thN8th_order(binx,DdevMOD_RollAccel_fourier_coeffs_binmean,0);
plot(binx,dev_steady+2*DdevMOD,'-g','linewidth',2)

% Pitch up&down
subplot(1,3,1)
strokeMOD = calc_val_fourier_series_4thN8th_order(binx,strokeMOD_PitchAccel_fourier_coeffs_binmean,0);
plot(binx,stroke_steady+2*strokeMOD,'-r','linewidth',2)
plot(binx,stroke_steady-2*strokeMOD,'-b','linewidth',2)

subplot(1,3,2)
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,pitchMOD_PitchAccel_fourier_coeffs_binmean,0);
plot(binx,pitch_steady+2*pitchMOD,'-r','linewidth',2)
plot(binx,pitch_steady-2*pitchMOD,'-b','linewidth',2)

subplot(1,3,3)
devMOD = calc_val_fourier_series_4thN8th_order(binx,devMOD_PitchAccel_fourier_coeffs_binmean,0);
plot(binx,dev_steady+2*devMOD,'-r','linewidth',2)
plot(binx,dev_steady-2*devMOD,'-b','linewidth',2)

% YAW fwd-rwd
subplot(1,3,1)
DstrokeMOD = calc_val_fourier_series_4thN8th_order(binx,DstrokeMOD_YawAccel_fourier_coeffs_binmean,0);
plot(binx,stroke_steady+2*DstrokeMOD,'-c','linewidth',2)

subplot(1,3,2)
DpitchMOD = calc_val_fourier_series_4thN8th_order(binx,DpitchMOD_YawAccel_fourier_coeffs_binmean,0);
plot(binx,pitch_steady+2*DpitchMOD,'-c','linewidth',2)

subplot(1,3,3)
DdevMOD = calc_val_fourier_series_4thN8th_order(binx,DdevMOD_YawAccel_fourier_coeffs_binmean,0);
plot(binx,dev_steady+2*DdevMOD,'-c','linewidth',2)

% DrotAxisR
subplot(1,3,1)
DstrokeMOD = calc_val_fourier_series_4thN8th_order(binx,DstrokeMOD_rotAxisRaccel_fourier_coeffs_binmean,0);
plot(binx,stroke_steady+2*DstrokeMOD,'-m','linewidth',2)

subplot(1,3,2)
DpitchMOD = calc_val_fourier_series_4thN8th_order(binx,DpitchMOD_rotAxisRaccel_fourier_coeffs_binmean,0);
plot(binx,pitch_steady+2*DpitchMOD,'-m','linewidth',2)

subplot(1,3,3)
DdevMOD = calc_val_fourier_series_4thN8th_order(binx,DdevMOD_rotAxisRaccel_fourier_coeffs_binmean,0);
plot(binx,dev_steady+2*DdevMOD,'-m','linewidth',2)

%%
subplot(1,3,1)
xlabel('t* (-)','fontsize',10) 
ylabel('stroke','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(1,3,2)
xlabel('t* (-)','fontsize',10) 
ylabel('pitch','fontsize',10) 
biny_min = 0;
biny_max = 180;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(1,3,3)
legend('steady','F/Mg','roll accel','pitch up','pitch down','yaw accel','rotAxisR accel')
xlabel('t* (-)','fontsize',10) 
ylabel('dev','fontsize',10) 
biny_min = -30;
biny_max = 30;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 



