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

imax = 200;

%% steady WB
stroke_steady = calc_val_fourier_series_4thN8th_order(binx,stroke_steady_fourier_coeffs_binmean,0);
pitch_steady = calc_val_fourier_series_4thN8th_order(binx,pitch_steady_fourier_coeffs_binmean,0);
dev_steady = calc_val_fourier_series_4thN8th_order(binx,dev_steady_fourier_coeffs_binmean,0);


%% Roll MOD wb
strokeMOD = calc_val_fourier_series_4thN8th_order(binx,DstrokeMOD_RollAccel_fourier_coeffs_binmean,0);
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,DpitchMOD_RollAccel_fourier_coeffs_binmean,0);
devMOD = calc_val_fourier_series_4thN8th_order(binx,DdevMOD_RollAccel_fourier_coeffs_binmean,0);

% color_map = autumn(imax);
color_map(1:imax,1) = 0;
% color_map(1:imax,1) = [.5:-.5/(imax-1):0];
% color_map(1:imax,1) = [.5:.5/(imax-1):1];

color_map(1:imax,2) = 0;
% color_map(1:imax,2) = [.5:-.5/(imax-1):0];

color_map(1:imax,3) = 0;
color_map(1:imax,3) = [0:1/(imax-1):1];
% color_map(1:imax,3) = [.5:-.5/(imax-1):0];
% color_map(1:imax,3) = [.5:.5/(imax-1):1];


% for i=1:imax
for i=1:imax-1
    k=(kmax*i)/imax;
    
    subplot(1,3,1)
    plot(binx,stroke_steady+k*strokeMOD,'-','color',color_map(i,:),'linewidth',.25)
    subplot(1,3,2)
    plot(binx,pitch_steady+k*pitchMOD,'-','color',color_map(i,:),'linewidth',.25)
    subplot(1,3,3)
    plot(binx,dev_steady+k*devMOD,'-','color',color_map(i,:),'linewidth',1)
end

% subplot(1,3,1)
% plot(binx,stroke_steady,'-','color','k','linewidth',1)
% subplot(1,3,2)
% plot(binx,pitch_steady,'-','color','k','linewidth',1)
% subplot(1,3,3)
% plot(binx,dev_steady,'-','color','k','linewidth',1)

subplot(1,3,1)
plot(binx,stroke_steady+kmax*strokeMOD,'-','color','k','linewidth',1)
subplot(1,3,2)
plot(binx,pitch_steady+kmax*pitchMOD,'-','color','k','linewidth',1)
subplot(1,3,3)
plot(binx,dev_steady+kmax*devMOD,'-','color','k','linewidth',1)

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
xlabel('t* (-)','fontsize',10) 
ylabel('dev','fontsize',10) 
biny_min = -30;
biny_max = 30;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 


