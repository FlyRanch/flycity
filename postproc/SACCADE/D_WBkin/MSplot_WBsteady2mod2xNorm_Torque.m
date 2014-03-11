figure
subplot(3,3,1)
% title('F/Mg')
hold on
subplot(3,3,2)
% title('Mroll')
hold on
subplot(3,3,3)
% title('Mpitch')
hold on
subplot(3,3,4)
% title('wing stroke')
hold on
subplot(3,3,5)
% title('wing pitch')
hold on
subplot(3,3,6)
% title('stroke deviation')
hold on
subplot(3,3,7)
% title('wing stroke')
hold on
subplot(3,3,8)
% title('wing pitch')
hold on
subplot(3,3,9)
% title('stroke deviation')
hold on

imax = 17;
kmax = 2;

%% steady WB
stroke_steady = calc_val_fourier_series_4thN8th_order(binx,stroke_steady_fourier_coeffs_binmean,0);
pitch_steady = calc_val_fourier_series_4thN8th_order(binx,pitch_steady_fourier_coeffs_binmean,0);
dev_steady = calc_val_fourier_series_4thN8th_order(binx,dev_steady_fourier_coeffs_binmean,0);


%% Force MOD wb
strokeMOD = calc_val_fourier_series_4thN8th_order(binx,strokeMOD_Fenhance_fourier_coeffs_binmean,0);
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,pitchMOD_Fenhance_fourier_coeffs_binmean,0);
devMOD = calc_val_fourier_series_4thN8th_order(binx,devMOD_Fenhance_fourier_coeffs_binmean,0);

% color_map = autumn(imax);
color_map(1:imax,1) = 0;
% color_map(1:imax,1) = [.5:-.5/(imax-1):0];
% color_map(1:imax,1) = [.5:.5/(imax-1):1];
% color_map(1:imax,1) = [0:.5/(imax-1):.5];

% color_map(1:imax,2) = 0;
% color_map(1:imax,2) = [.5:-.5/(imax-1):0];
% color_map(1:imax,2) = [0:.25/(imax-1):.25];
color_map(1:imax,2) = [0:.5/(imax-1):.5];

color_map(1:imax,3) = 0;
color_map(1:imax,3) = [0:1/(imax-1):1];
% color_map(1:imax,3) = [.5:-.5/(imax-1):0];
% color_map(1:imax,3) = [.5:.5/(imax-1):1];


for i=1:imax
    k=kmax*(i-1)/(imax-1)
    
    subplot(3,3,1)
    plot(binx,stroke_steady+k*strokeMOD,'-','color',color_map(i,:),'linewidth',1)
    subplot(3,3,4)
    plot(binx,pitch_steady+k*pitchMOD,'-','color',color_map(i,:),'linewidth',1)
    subplot(3,3,7)
    plot(binx,dev_steady+k*devMOD,'-','color',color_map(i,:),'linewidth',1)
end

% % subplot(3,3,1)
% % plot(binx,stroke_steady,'-','color','k','linewidth',1)
% % subplot(3,3,2)
% % plot(binx,pitch_steady,'-','color','k','linewidth',1)
% % subplot(3,3,3)
% % plot(binx,dev_steady,'-','color','k','linewidth',1)
% 
% subplot(3,3,1)
% plot(binx,stroke_steady+kmax*strokeMOD,'-','color','k','linewidth',.25)
% subplot(3,3,2)
% plot(binx,pitch_steady+kmax*pitchMOD,'-','color','k','linewidth',.25)
% subplot(3,3,3)
% plot(binx,dev_steady+kmax*devMOD,'-','color','k','linewidth',.25)

%% Roll MOD wb
strokeMOD = calc_val_fourier_series_4thN8th_order(binx,DstrokeMOD_RollTorque_fourier_coeffs_binmean,0);
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,DpitchMOD_RollTorque_fourier_coeffs_binmean,0);
devMOD = calc_val_fourier_series_4thN8th_order(binx,DdevMOD_RollTorque_fourier_coeffs_binmean,0);

% color_map = autumn(imax);
color_map(1:imax,1) = 0;
% color_map(1:imax,1) = [.5:-.5/(imax-1):0];
% color_map(1:imax,1) = [.5:.5/(imax-1):1];

% color_map(1:imax,2) = 0;
% color_map(1:imax,2) = [.5:-.5/(imax-1):0];
% color_map(1:imax,2) = [0:.25/(imax-1):.25];
color_map(1:imax,2) = [0:.5/(imax-1):.5];

color_map(1:imax,3) = 0;
color_map(1:imax,3) = [0:1/(imax-1):1];
% color_map(1:imax,3) = [.5:-.5/(imax-1):0];
% color_map(1:imax,3) = [.5:.5/(imax-1):1];


for i=1:imax
    k=kmax*(i-1)/(imax-1)
    
    subplot(3,3,2)
    plot(binx,stroke_steady+k*strokeMOD,'-','color',color_map(i,:),'linewidth',1)
    subplot(3,3,5)
    plot(binx,pitch_steady+k*pitchMOD,'-','color',color_map(i,:),'linewidth',1)
    subplot(3,3,8)
    plot(binx,dev_steady+k*devMOD,'-','color',color_map(i,:),'linewidth',1)
end
% 
% % subplot(3,3,1)
% % plot(binx,stroke_steady,'-','color','k','linewidth',1)
% % subplot(3,3,2)
% % plot(binx,pitch_steady,'-','color','k','linewidth',1)
% % subplot(3,3,3)
% % plot(binx,dev_steady,'-','color','k','linewidth',1)
% 
% subplot(3,3,4)
% plot(binx,stroke_steady+kmax*strokeMOD,'-','color','k','linewidth',.25)
% subplot(3,3,5)
% plot(binx,pitch_steady+kmax*pitchMOD,'-','color','k','linewidth',.25)
% subplot(3,3,6)
% plot(binx,dev_steady+kmax*devMOD,'-','color','k','linewidth',.25)

%% Pitch MOD wb
strokeMOD = calc_val_fourier_series_4thN8th_order(binx,strokeMOD_PitchTorque_fourier_coeffs_binmean,0);
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,pitchMOD_PitchTorque_fourier_coeffs_binmean,0);
devMOD = calc_val_fourier_series_4thN8th_order(binx,devMOD_PitchTorque_fourier_coeffs_binmean,0);

%% pitch up
% color_map = autumn(imax);
color_map(1:imax,1) = 0;
% color_map(1:imax,1) = [.5:-.5/(imax-1):0];
% color_map(1:imax,1) = [.5:.5/(imax-1):1];

% color_map(1:imax,2) = 0;
% color_map(1:imax,2) = [.5:-.5/(imax-1):0];
% color_map(1:imax,2) = [0:.25/(imax-1):.25];
color_map(1:imax,2) = [0:.5/(imax-1):.5];

color_map(1:imax,3) = 0;
color_map(1:imax,3) = [0:1/(imax-1):1];
% color_map(1:imax,3) = [.5:-.5/(imax-1):0];
% color_map(1:imax,3) = [.5:.5/(imax-1):1];


for i=1:imax
% for i=1:imax-1
    k=(kmax*(i-1))/(imax-1);
    
    subplot(3,3,3)
    plot(binx,stroke_steady+k*strokeMOD,'-','color',color_map(i,:),'linewidth',.25)
    subplot(3,3,6)
    plot(binx,pitch_steady+k*pitchMOD,'-','color',color_map(i,:),'linewidth',.25)
    subplot(3,3,9)
    plot(binx,dev_steady+k*devMOD,'-','color',color_map(i,:),'linewidth',.5)
end

%% pitch down
% color_map = autumn(imax);
% color_map(1:imax,1) = 0;
color_map(1:imax,1) = [0:1/(imax-1):1];
% color_map(1:imax,1) = [.5:-.5/(imax-1):0];
% color_map(1:imax,1) = [.5:.5/(imax-1):1];

% color_map(1:imax,2) = 0;
% color_map(1:imax,2) = [.5:-.5/(imax-1):0];
% color_map(1:imax,2) = [0:.5/(imax-1):.5];
color_map(1:imax,2) = [0:.1/(imax-1):.1];

% color_map(1:imax,3) = 0;
% color_map(1:imax,3) = [.5:-.5/(imax-1):0];
% color_map(1:imax,3) = [.5:.5/(imax-1):1];
color_map(1:imax,3) = [0:.1/(imax-1):.1];

for i=1:imax
% for i=1:imax-1
    k=(kmax*(i-1))/(imax-1);
    
    subplot(3,3,3)
    plot(binx,stroke_steady-k*strokeMOD,'-','color',color_map(i,:),'linewidth',.25)
    subplot(3,3,6)
    plot(binx,pitch_steady-k*pitchMOD,'-','color',color_map(i,:),'linewidth',.25)
    subplot(3,3,9)
    plot(binx,dev_steady-k*devMOD,'-','color',color_map(i,:),'linewidth',.5)
end

% % subplot(3,3,1)
% % plot(binx,stroke_steady,'-','color','k','linewidth',1)
% % subplot(3,3,2)
% % plot(binx,pitch_steady,'-','color','k','linewidth',1)
% % subplot(3,3,3)
% % plot(binx,dev_steady,'-','color','k','linewidth',1)
% 
% subplot(3,3,7)
% plot(binx,stroke_steady+kmax*strokeMOD,'-','color','k','linewidth',.25)
% subplot(3,3,8)
% plot(binx,pitch_steady+kmax*pitchMOD,'-','color','k','linewidth',.25)
% subplot(3,3,9)
% plot(binx,dev_steady+kmax*devMOD,'-','color','k','linewidth',.25)
% 
% subplot(3,3,7)
% plot(binx,stroke_steady-kmax*strokeMOD,'-','color','k','linewidth',.25)
% subplot(3,3,8)
% plot(binx,pitch_steady-kmax*pitchMOD,'-','color','k','linewidth',.25)
% subplot(3,3,9)
% plot(binx,dev_steady-kmax*devMOD,'-','color','k','linewidth',.25)

%%
subplot(3,3,1)
xlabel('t* (-)','fontsize',10) 
ylabel('stroke','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(3,3,4)
xlabel('t* (-)','fontsize',10) 
ylabel('pitch','fontsize',10) 
biny_min = 0;
biny_max = 180;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(3,3,7)
xlabel('t* (-)','fontsize',10) 
ylabel('dev','fontsize',10) 
biny_min = -30;
biny_max = 30;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 



%%
subplot(3,3,2)
xlabel('t* (-)','fontsize',10) 
ylabel('stroke','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(3,3,5)
xlabel('t* (-)','fontsize',10) 
ylabel('pitch','fontsize',10) 
biny_min = 0;
biny_max = 180;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(3,3,8)
xlabel('t* (-)','fontsize',10) 
ylabel('dev','fontsize',10) 
biny_min = -30;
biny_max = 30;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 


%%
subplot(3,3,3)
xlabel('t* (-)','fontsize',10) 
ylabel('stroke','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(3,3,6)
xlabel('t* (-)','fontsize',10) 
ylabel('pitch','fontsize',10) 
biny_min = 0;
biny_max = 180;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(3,3,9)
xlabel('t* (-)','fontsize',10) 
ylabel('dev','fontsize',10) 
biny_min = -30;
biny_max = 30;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 


