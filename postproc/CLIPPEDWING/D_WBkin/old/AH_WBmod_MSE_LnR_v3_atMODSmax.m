clc
clear
close all

name = 'WBdataset_all_steadyNmods_TorqueNorm.mat'
load(name)
load('norm_data_torque.mat')

% correlation variables
Dt_deg = 10
Dt = Dt_deg/360
Da = 10

Nt = 100 % Nx2
Na = 100 % Nx2

dt = Dt/Nt % norm time
dt_deg = 360*dt
da = Da/Na % deg

%% CONSTANT WING
% steady WB
stroke_steady = calc_val_fourier_series_4thN8th_order(binx,stroke_steady_fourier_coeffs_binmean,0);
pitch_steady  = calc_val_fourier_series_4thN8th_order(binx,pitch_steady_fourier_coeffs_binmean,0);
dev_steady    = calc_val_fourier_series_4thN8th_order(binx,dev_steady_fourier_coeffs_binmean,0);

% fwd wing 
stroke_Myaw_FWD = stroke_steady + max(YawMods)*calc_val_fourier_series_4thN8th_order(binx,strokeMOD_fwd_YawTorque_fourier_coeffs_binmean,0);
pitch_Myaw_FWD  = pitch_steady  + max(YawMods)*calc_val_fourier_series_4thN8th_order(binx,pitchMOD_fwd_YawTorque_fourier_coeffs_binmean,0);
dev_Myaw_FWD    = dev_steady    + max(YawMods)*calc_val_fourier_series_4thN8th_order(binx,devMOD_fwd_YawTorque_fourier_coeffs_binmean,0);

stroke_MaxisR_L = stroke_steady + max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx,strokeMOD_L_TorqueAxisR_fourier_coeffs_binmean,0);
pitch_MaxisR_L  = pitch_steady  + max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx,pitchMOD_L_TorqueAxisR_fourier_coeffs_binmean,0);
dev_MaxisR_L    = dev_steady    + max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx,devMOD_L_TorqueAxisR_fourier_coeffs_binmean,0);

% rwd wing 
stroke_Myaw_RWD0 = stroke_steady + max(YawMods)*calc_val_fourier_series_4thN8th_order(binx,strokeMOD_rwd_YawTorque_fourier_coeffs_binmean,0);
pitch_Myaw_RWD0  = pitch_steady  + max(YawMods)*calc_val_fourier_series_4thN8th_order(binx,pitchMOD_rwd_YawTorque_fourier_coeffs_binmean,0);
dev_Myaw_RWD0   = dev_steady    + max(YawMods)*calc_val_fourier_series_4thN8th_order(binx,devMOD_rwd_YawTorque_fourier_coeffs_binmean,0);

stroke_MaxisR_R0 = stroke_steady + max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx,strokeMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);
pitch_MaxisR_R0  = pitch_steady  + max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx,pitchMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);
dev_MaxisR_R0    = dev_steady    + max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx,devMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);

% dt loop
for i=1:2*Nt+1
    count = 2*Nt+1 -i
    dt_now = (i-Nt-1)*dt;
    dt_MSE(i,1) = dt_now;

    %% VARIABLE WING
    % steady WB
    stroke_steady_dt = calc_val_fourier_series_4thN8th_order(binx-dt_now,stroke_steady_fourier_coeffs_binmean,0);
    pitch_steady_dt  = calc_val_fourier_series_4thN8th_order(binx-dt_now,pitch_steady_fourier_coeffs_binmean,0);
    dev_steady_dt    = calc_val_fourier_series_4thN8th_order(binx-dt_now,dev_steady_fourier_coeffs_binmean,0);
    
    % rwd wing
    stroke_Myaw_RWD = stroke_steady_dt + max(YawMods)*calc_val_fourier_series_4thN8th_order(binx-dt_now,strokeMOD_rwd_YawTorque_fourier_coeffs_binmean,0);
    pitch_Myaw_RWD  = pitch_steady_dt  + max(YawMods)*calc_val_fourier_series_4thN8th_order(binx-dt_now,pitchMOD_rwd_YawTorque_fourier_coeffs_binmean,0);
    dev_Myaw_RWD    = dev_steady_dt    + max(YawMods)*calc_val_fourier_series_4thN8th_order(binx-dt_now,devMOD_rwd_YawTorque_fourier_coeffs_binmean,0);

    stroke_MaxisR_R = stroke_steady_dt + max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx-dt_now,strokeMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);
    pitch_MaxisR_R  = pitch_steady_dt  + max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx-dt_now,pitchMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);
    dev_MaxisR_R    = dev_steady_dt    + max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx-dt_now,devMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);

    % d angle loop
    for j=1:2*Na+1
        da_now = (j-Na-1)*da;
        da_MSE(j,1) = da_now;
        
        % WINGKIN DIFF
        Dstroke_Myaw = stroke_Myaw_FWD - stroke_Myaw_RWD - da_now;
        Dpitch_Myaw  = pitch_Myaw_FWD  - pitch_Myaw_RWD  - da_now;
        Ddev_Myaw    = dev_Myaw_FWD    - dev_Myaw_RWD    - da_now;
        
        Dstroke_MaxisR = stroke_MaxisR_L - stroke_MaxisR_R - da_now;
        Dpitch_MaxisR  = pitch_MaxisR_L  - pitch_MaxisR_R  - da_now;
        Ddev_MaxisR    = dev_MaxisR_L    - dev_MaxisR_R    - da_now;
        
        MSE_stroke_Myaw(i,j) = mse(Dstroke_Myaw);
        MSE_pitch_Myaw(i,j) = mse(Dpitch_Myaw);
        MSE_dev_Myaw(i,j) = mse(Ddev_Myaw);
        
        MSE_stroke_MaxisR(i,j) = mse(Dstroke_MaxisR);
        MSE_pitch_MaxisR(i,j) = mse(Dpitch_MaxisR);
        MSE_dev_MaxisR(i,j) = mse(Ddev_MaxisR);
    end
end

dt_deg_MSE = dt_MSE*360; % phase shift in deg;
[DA_MSE,DT_DEG_MSE] = meshgrid(da_MSE,dt_deg_MSE);
[DA_MSE,DT_MSE] = meshgrid(da_MSE,dt_MSE);

% min MSE
MSE_stroke_Myaw_MSEmin = min(MSE_stroke_Myaw(:));
dt_stroke_Myaw_MSEmin = DT_MSE(MSE_stroke_Myaw==MSE_stroke_Myaw_MSEmin);
phase_stroke_Myaw_MSEmin = DT_DEG_MSE(MSE_stroke_Myaw==MSE_stroke_Myaw_MSEmin);
dstroke_Myaw_MSEmin = DA_MSE(MSE_stroke_Myaw==MSE_stroke_Myaw_MSEmin);

MSE_pitch_Myaw_MSEmin = min(MSE_pitch_Myaw(:));
dt_pitch_Myaw_MSEmin = DT_MSE(MSE_pitch_Myaw==MSE_pitch_Myaw_MSEmin);
phase_pitch_Myaw_MSEmin = DT_DEG_MSE(MSE_pitch_Myaw==MSE_pitch_Myaw_MSEmin);
dpitch_Myaw_MSEmin = DA_MSE(MSE_pitch_Myaw==MSE_pitch_Myaw_MSEmin);

MSE_dev_Myaw_MSEmin = min(MSE_dev_Myaw(:));
dt_dev_Myaw_MSEmin = DT_MSE(MSE_dev_Myaw==MSE_dev_Myaw_MSEmin);
phase_dev_Myaw_MSEmin = DT_DEG_MSE(MSE_dev_Myaw==MSE_dev_Myaw_MSEmin);
ddev_Myaw_MSEmin = DA_MSE(MSE_dev_Myaw==MSE_dev_Myaw_MSEmin);

MSE_stroke_MaxisR_MSEmin = min(MSE_stroke_MaxisR(:));
dt_stroke_MaxisR_MSEmin = DT_MSE(MSE_stroke_MaxisR==MSE_stroke_MaxisR_MSEmin);
phase_stroke_MaxisR_MSEmin = DT_DEG_MSE(MSE_stroke_MaxisR==MSE_stroke_MaxisR_MSEmin);
dstroke_MaxisR_MSEmin = DA_MSE(MSE_stroke_MaxisR==MSE_stroke_MaxisR_MSEmin);

MSE_pitch_MaxisR_MSEmin = min(MSE_pitch_MaxisR(:));
dt_pitch_MaxisR_MSEmin = DT_MSE(MSE_pitch_MaxisR==MSE_pitch_MaxisR_MSEmin);
phase_pitch_MaxisR_MSEmin = DT_DEG_MSE(MSE_pitch_MaxisR==MSE_pitch_MaxisR_MSEmin);
dpitch_MaxisR_MSEmin = DA_MSE(MSE_pitch_MaxisR==MSE_pitch_MaxisR_MSEmin);

MSE_dev_MaxisR_MSEmin = min(MSE_dev_MaxisR(:));
dt_dev_MaxisR_MSEmin = DT_MSE(MSE_dev_MaxisR==MSE_dev_MaxisR_MSEmin);
phase_dev_MaxisR_MSEmin = DT_DEG_MSE(MSE_dev_MaxisR==MSE_dev_MaxisR_MSEmin);
ddev_MaxisR_MSEmin = DA_MSE(MSE_dev_MaxisR==MSE_dev_MaxisR_MSEmin);

%% WBs at MSEmin

% Yaw
stroke_steady_YawDt = calc_val_fourier_series_4thN8th_order(binx-dt_stroke_Myaw_MSEmin,stroke_steady_fourier_coeffs_binmean,0);
pitch_steady_YawDt  = calc_val_fourier_series_4thN8th_order(binx-dt_pitch_Myaw_MSEmin,pitch_steady_fourier_coeffs_binmean,0);
dev_steady_YawDt    = calc_val_fourier_series_4thN8th_order(binx-dt_dev_Myaw_MSEmin,dev_steady_fourier_coeffs_binmean,0);

stroke_Myaw_RWD_MSEmin = stroke_steady_YawDt + max(YawMods)*calc_val_fourier_series_4thN8th_order(binx-dt_stroke_Myaw_MSEmin,strokeMOD_rwd_YawTorque_fourier_coeffs_binmean,0);
pitch_Myaw_RWD_MSEmin  = pitch_steady_YawDt  + max(YawMods)*calc_val_fourier_series_4thN8th_order(binx-dt_pitch_Myaw_MSEmin,pitchMOD_rwd_YawTorque_fourier_coeffs_binmean,0);
dev_Myaw_RWD_MSEmin    = dev_steady_YawDt    + max(YawMods)*calc_val_fourier_series_4thN8th_order(binx-dt_dev_Myaw_MSEmin,devMOD_rwd_YawTorque_fourier_coeffs_binmean,0);

stroke_Myaw_RWD_MSEmin = stroke_Myaw_RWD_MSEmin + dstroke_Myaw_MSEmin;
pitch_Myaw_RWD_MSEmin  = pitch_Myaw_RWD_MSEmin  + dpitch_Myaw_MSEmin;
dev_Myaw_RWD_MSEmin    = dev_Myaw_RWD_MSEmin    + ddev_Myaw_MSEmin;

% axisR
stroke_steady_axisRdt = calc_val_fourier_series_4thN8th_order(binx-dt_stroke_MaxisR_MSEmin,stroke_steady_fourier_coeffs_binmean,0);
pitch_steady_axisRdt  = calc_val_fourier_series_4thN8th_order(binx-dt_pitch_MaxisR_MSEmin,pitch_steady_fourier_coeffs_binmean,0);
dev_steady_axisRdt    = calc_val_fourier_series_4thN8th_order(binx-dt_dev_MaxisR_MSEmin,dev_steady_fourier_coeffs_binmean,0);

stroke_MaxisR_R_MSEmin = stroke_steady_axisRdt + max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx-dt_stroke_MaxisR_MSEmin,strokeMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);
pitch_MaxisR_R_MSEmin  = pitch_steady_axisRdt  + max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx-dt_pitch_MaxisR_MSEmin,pitchMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);
dev_MaxisR_R_MSEmin    = dev_steady_axisRdt    + max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx-dt_dev_MaxisR_MSEmin,devMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);

stroke_MaxisR_R_MSEmin = stroke_MaxisR_R_MSEmin + dstroke_MaxisR_MSEmin;
pitch_MaxisR_R_MSEmin  = pitch_MaxisR_R_MSEmin  + dpitch_MaxisR_MSEmin;
dev_MaxisR_R_MSEmin    = dev_MaxisR_R_MSEmin    + ddev_MaxisR_MSEmin;

%% save data
save('WBkin_LnR_correlations_atMODmax.mat','binx',...
'dt_deg_MSE','dt_MSE','da_MSE','DA_MSE','DT_DEG_MSE','DT_MSE',...
'MSE_stroke_Myaw_MSEmin','MSE_stroke_Myaw','dt_stroke_Myaw_MSEmin','phase_stroke_Myaw_MSEmin','dstroke_Myaw_MSEmin',...
'MSE_pitch_Myaw_MSEmin','MSE_pitch_Myaw','dt_pitch_Myaw_MSEmin','phase_pitch_Myaw_MSEmin','dpitch_Myaw_MSEmin',...
'MSE_dev_Myaw_MSEmin','MSE_dev_Myaw','dt_dev_Myaw_MSEmin','phase_dev_Myaw_MSEmin','ddev_Myaw_MSEmin',...
'MSE_stroke_MaxisR_MSEmin','MSE_stroke_MaxisR','dt_stroke_MaxisR_MSEmin','phase_stroke_MaxisR_MSEmin','dstroke_MaxisR_MSEmin',...
'MSE_pitch_MaxisR_MSEmin','MSE_pitch_MaxisR','dt_pitch_MaxisR_MSEmin','phase_pitch_MaxisR_MSEmin','dpitch_MaxisR_MSEmin',...
'MSE_dev_MaxisR_MSEmin','MSE_dev_MaxisR','dt_dev_MaxisR_MSEmin','phase_dev_MaxisR_MSEmin','ddev_MaxisR_MSEmin',...
'stroke_Myaw_FWD','stroke_Myaw_RWD0','stroke_Myaw_RWD_MSEmin',...
'pitch_Myaw_FWD','pitch_Myaw_RWD0','pitch_Myaw_RWD_MSEmin',...
'dev_Myaw_FWD','dev_Myaw_RWD0','dev_Myaw_RWD_MSEmin',...
'stroke_MaxisR_L','stroke_MaxisR_R0','stroke_MaxisR_R_MSEmin',...
'pitch_MaxisR_L','pitch_MaxisR_R0','pitch_MaxisR_R_MSEmin',...
'dev_MaxisR_L','dev_MaxisR_R0','dev_MaxisR_R_MSEmin')

%% plots
mkdir('MSplot_WBkin_LnR_correlations')
cd('MSplot_WBkin_LnR_correlations')

%% WBs at MSEmin
figure
subplot(3,3,1)
hold on
subplot(3,3,2)
hold on
subplot(3,3,3)
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

% yaw @ MSEmin
subplot(3,3,1)
plot(binx,stroke_Myaw_FWD,'-','color',[1 .5 0],'linewidth',1)
plot(binx,stroke_Myaw_RWD_MSEmin,'-','color',[0 .5 1],'linewidth',1)

plot(binx,pitch_Myaw_FWD-90,'-','color',[1 .5 0],'linewidth',1)
plot(binx,pitch_Myaw_RWD_MSEmin-90,'-','color',[0 .5 1],'linewidth',1)

plot(binx,dev_Myaw_FWD,'-','color',[1 .5 0],'linewidth',1)
plot(binx,dev_Myaw_RWD_MSEmin,'-','color',[0 .5 1],'linewidth',1)

% axisR @ MSEmin
subplot(3,3,2)
plot(binx,stroke_MaxisR_L,'-','color',[1 .5 0],'linewidth',1)
plot(binx,stroke_MaxisR_R_MSEmin,'-','color',[0 .5 1],'linewidth',1)

plot(binx,pitch_MaxisR_L-90,'-','color',[1 .5 0],'linewidth',1)
plot(binx,pitch_MaxisR_R_MSEmin-90,'-','color',[0 .5 1],'linewidth',1)

plot(binx,dev_MaxisR_L,'-','color',[1 .5 0],'linewidth',1)
plot(binx,dev_MaxisR_R_MSEmin,'-','color',[0 .5 1],'linewidth',1)

% yaw original
subplot(3,3,4)
plot(binx,stroke_Myaw_FWD,'-','color',[1 .5 0],'linewidth',1)
plot(binx,stroke_Myaw_RWD0,'-','color',[0 .5 1],'linewidth',1)

plot(binx,pitch_Myaw_FWD-90,'-','color',[1 .5 0],'linewidth',1)
plot(binx,pitch_Myaw_RWD0-90,'-','color',[0 .5 1],'linewidth',1)

plot(binx,dev_Myaw_FWD,'-','color',[1 .5 0],'linewidth',1)
plot(binx,dev_Myaw_RWD0,'-','color',[0 .5 1],'linewidth',1)

% axisR original
subplot(3,3,5)
plot(binx,stroke_MaxisR_L,'-','color',[1 .5 0],'linewidth',1)
plot(binx,stroke_MaxisR_R0,'-','color',[0 .5 1],'linewidth',1)

plot(binx,pitch_MaxisR_L-90,'-','color',[1 .5 0],'linewidth',1)
plot(binx,pitch_MaxisR_R0-90,'-','color',[0 .5 1],'linewidth',1)

plot(binx,dev_MaxisR_L,'-','color',[1 .5 0],'linewidth',1)
plot(binx,dev_MaxisR_R0,'-','color',[0 .5 1],'linewidth',1)

% yaw original & @ MSEmin
subplot(3,3,7)
plot(binx,stroke_Myaw_FWD,'-','color',[1 .5 0],'linewidth',1)
plot(binx,stroke_Myaw_RWD_MSEmin,'-','color',[0 .5 1],'linewidth',1)
plot(binx,stroke_Myaw_RWD0,'-','color',[.5 .5 .5],'linewidth',1)

plot(binx,pitch_Myaw_FWD-90,'-','color',[1 .5 0],'linewidth',1)
plot(binx,pitch_Myaw_RWD_MSEmin-90,'-','color',[0 .5 1],'linewidth',1)
plot(binx,pitch_Myaw_RWD0-90,'-','color',[.5 .5 .5],'linewidth',1)

plot(binx,dev_Myaw_FWD,'-','color',[1 .5 0],'linewidth',1)
plot(binx,dev_Myaw_RWD_MSEmin,'-','color',[0 .5 1],'linewidth',1)
plot(binx,dev_Myaw_RWD0,'-','color',[.5 .5 .5],'linewidth',1)

% axisR original & @ MSEmin
subplot(3,3,8)
plot(binx,stroke_MaxisR_L,'-','color',[1 .5 0],'linewidth',1)
plot(binx,stroke_MaxisR_R_MSEmin,'-','color',[0 .5 1],'linewidth',1)
plot(binx,stroke_MaxisR_R0,'-','color',[.5 .5 .5],'linewidth',1)

plot(binx,pitch_MaxisR_L-90,'-','color',[1 .5 0],'linewidth',1)
plot(binx,pitch_MaxisR_R_MSEmin-90,'-','color',[0 .5 1],'linewidth',1)
plot(binx,pitch_MaxisR_R0-90,'-','color',[.5 .5 .5],'linewidth',1)

plot(binx,dev_MaxisR_L,'-','color',[1 .5 0],'linewidth',1)
plot(binx,dev_MaxisR_R_MSEmin,'-','color',[0 .5 1],'linewidth',1)
plot(binx,dev_MaxisR_R0,'-','color',[.5 .5 .5],'linewidth',1)

% legends
subplot(3,3,1)
% xlabel('t* (-)','fontsize',10) 
ylabel('yaw kin @ MSEmin','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(3,3,2)
% xlabel('t* (-)','fontsize',10) 
ylabel('Raxis up kin @ MSEmin','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(3,3,4)
% xlabel('t* (-)','fontsize',10) 
ylabel('yaw kin @ 3std','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(3,3,5)
% xlabel('t* (-)','fontsize',10) 
ylabel('Raxis up kin @ 3std','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(3,3,7)
xlabel('t* (-)','fontsize',10) 
ylabel('yaw kin @ 3std','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(3,3,8)
xlabel('t* (-)','fontsize',10) 
ylabel('Raxis up kin @ 3std','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

saveas(gca,['MSplot_WBkinMODmax_LnR_atMSEmin.fig'])
saveas(gca,['MSplot_WBkinMODmax_LnR_atMSEmin.png'])
plot2svg(['MSplot_WBkinMODmax_LnR_atMSEmin.svg'])

%% MSE plot
MSEmax = 100;

figure
colormap gray

subplot(3,2,1)
% surf(dt_deg_MSE, da_MSE, (MSE_stroke_Myaw'),'EdgeColor','none')
% shading interp
% view(2)
% hold on
% plot3([-10 10],[0 0],[10 10],':w')
% plot3([0 0],[-10 10],[10 10],':w')
contourf(dt_deg_MSE, da_MSE, (MSE_stroke_Myaw'),50,'LineColor','none');
set(gca,'clim',[0 MSEmax])
hold on
plot(phase_stroke_Myaw_MSEmin,dstroke_Myaw_MSEmin,'.w')
plot([-10 10],[0 0],':w')
plot([0 0],[-10 10],':w')
title('ln MSE yaw torque')
ylabel('stroke kin angle')
axis equal
axis tight

subplot(3,2,3)
% surf(dt_deg_MSE, da_MSE, (MSE_pitch_Myaw'),'EdgeColor','none')
% shading interp
% view(2)
% hold on
% plot3([-10 10],[0 0],[10 10],':w')
% plot3([0 0],[-10 10],[10 10],':w')
contourf(dt_deg_MSE, da_MSE, (MSE_pitch_Myaw'),50,'LineColor','none');
set(gca,'clim',[0 MSEmax])
hold on
plot(phase_pitch_Myaw_MSEmin,dpitch_Myaw_MSEmin,'.w')
plot([-10 10],[0 0],':w')
plot([0 0],[-10 10],':w')
% title('yaw rotation')
ylabel('rotation kin angle')
axis equal
axis tight

subplot(3,2,5)
% surf(dt_deg_MSE, da_MSE, (MSE_dev_Myaw'),'EdgeColor','none')
% shading interp
% view(2)
% hold on
% plot3([-10 10],[0 0],[10 10],':w')
% plot3([0 0],[-10 10],[10 10],':w')
contourf(dt_deg_MSE, da_MSE, (MSE_dev_Myaw'),50,'LineColor','none');
set(gca,'clim',[0 MSEmax])
hold on
plot(phase_dev_Myaw_MSEmin,ddev_Myaw_MSEmin,'.w')
plot([-10 10],[0 0],':w')
plot([0 0],[-10 10],':w')
% title('yaw dev')
ylabel('deviation kin angle')
xlabel('phase')
axis equal
axis tight

subplot(3,2,2)
% surf(dt_deg_MSE, da_MSE, (MSE_stroke_MaxisR'),'EdgeColor','none')
% shading interp
% view(2)
% hold on
% plot3([-10 10],[0 0],[10 10],':w')
% plot3([0 0],[-10 10],[10 10],':w')
contourf(dt_deg_MSE, da_MSE, (MSE_stroke_MaxisR'),50,'LineColor','none');
set(gca,'clim',[0 MSEmax])
hold on
plot(phase_stroke_MaxisR_MSEmin,dstroke_MaxisR_MSEmin,'.w')
plot([-10 10],[0 0],':w')
plot([0 0],[-10 10],':w')
title('ln MSE axisR torque')
% ylabel('kin angle')
axis equal
axis tight

subplot(3,2,4)
% surf(dt_deg_MSE, da_MSE, (MSE_pitch_MaxisR'),'EdgeColor','none')
% shading interp
% view(2)
% hold on
% plot3([-10 10],[0 0],[10 10],':w')
% plot3([0 0],[-10 10],[10 10],':w')
contourf(dt_deg_MSE, da_MSE, (MSE_pitch_MaxisR'),50,'LineColor','none');
set(gca,'clim',[0 MSEmax])
hold on
plot(phase_pitch_MaxisR_MSEmin,dpitch_MaxisR_MSEmin,'.w')
plot([-10 10],[0 0],':w')
plot([0 0],[-10 10],':w')
% title('axisR rotation')
% ylabel('kin angle')
axis equal
axis tight

subplot(3,2,6)
% surf(dt_deg_MSE, da_MSE, (MSE_dev_MaxisR'),'EdgeColor','none')
% shading interp
% view(2)
% hold on
% plot3([-10 10],[0 0],[10 10],':w')
% plot3([0 0],[-10 10],[10 10],':w')
contourf(dt_deg_MSE, da_MSE, (MSE_dev_MaxisR'),50,'LineColor','none');
set(gca,'clim',[0 MSEmax])
hold on
plot(phase_dev_MaxisR_MSEmin,ddev_MaxisR_MSEmin,'.w')
plot([-10 10],[0 0],':w')
plot([0 0],[-10 10],':w')
% title('axisR dev')
% ylabel('kin angle')
xlabel('phase')
axis equal
axis tight
colorbar

saveas(gca,['MSplot_MSE_WBkinMODmax_LnR_MSEmax',num2str(MSEmax),'.fig'])
saveas(gca,['MSplot_MSE_WBkinMODmax_LnR_MSEmax',num2str(MSEmax),'.png'])
plot2svg(['MSplot_MSE_WBkinMODmax_LnR_MSEmax',num2str(MSEmax),'.svg'])

%% log MSE plot
figure
colormap gray
logMSEmax = 2;

subplot(3,2,1)
% surf(dt_deg_MSE, da_MSE, log10(MSE_stroke_Myaw'),'EdgeColor','none')
% shading interp
% view(2)
% hold on
% plot3([-10 10],[0 0],[10 10],':w')
% plot3([0 0],[-10 10],[10 10],':w')
contourf(dt_deg_MSE, da_MSE, log10(MSE_stroke_Myaw'),50,'LineColor','none');
set(gca,'clim',[0 logMSEmax])
hold on
plot(phase_stroke_Myaw_MSEmin,dstroke_Myaw_MSEmin,'.w')
plot([-10 10],[0 0],':w')
plot([0 0],[-10 10],':w')
title('log MSE yaw torque')
ylabel('stroke kin angle')
axis equal
axis tight

subplot(3,2,3)
% surf(dt_deg_MSE, da_MSE, log10(MSE_pitch_Myaw'),'EdgeColor','none')
% shading interp
% view(2)
% hold on
% plot3([-10 10],[0 0],[10 10],':w')
% plot3([0 0],[-10 10],[10 10],':w')
contourf(dt_deg_MSE, da_MSE, log10(MSE_pitch_Myaw'),50,'LineColor','none');
set(gca,'clim',[0 logMSEmax])
hold on
plot(phase_pitch_Myaw_MSEmin,dpitch_Myaw_MSEmin,'.w')
plot([-10 10],[0 0],':w')
plot([0 0],[-10 10],':w')
% title('yaw rotation')
ylabel('rotation kin angle')
axis equal
axis tight

subplot(3,2,5)
% surf(dt_deg_MSE, da_MSE, log10(MSE_dev_Myaw'),'EdgeColor','none')
% shading interp
% view(2)
% hold on
% plot3([-10 10],[0 0],[10 10],':w')
% plot3([0 0],[-10 10],[10 10],':w')
contourf(dt_deg_MSE, da_MSE, log10(MSE_dev_Myaw'),50,'LineColor','none');
set(gca,'clim',[0 logMSEmax])
hold on
plot(phase_dev_Myaw_MSEmin,ddev_Myaw_MSEmin,'.w')
plot([-10 10],[0 0],':w')
plot([0 0],[-10 10],':w')
% title('yaw dev')
ylabel('deviation kin angle')
xlabel('phase')
axis equal
axis tight

subplot(3,2,2)
% surf(dt_deg_MSE, da_MSE, log10(MSE_stroke_MaxisR'),'EdgeColor','none')
% shading interp
% view(2)
% hold on
% plot3([-10 10],[0 0],[10 10],':w')
% plot3([0 0],[-10 10],[10 10],':w')
contourf(dt_deg_MSE, da_MSE, log10(MSE_stroke_MaxisR'),50,'LineColor','none');
set(gca,'clim',[0 logMSEmax])
hold on
plot(phase_stroke_MaxisR_MSEmin,dstroke_MaxisR_MSEmin,'.w')
plot([-10 10],[0 0],':w')
plot([0 0],[-10 10],':w')
title('log MSE axisR torque')
% ylabel('kin angle')
axis equal
axis tight

subplot(3,2,4)
% surf(dt_deg_MSE, da_MSE, log10(MSE_pitch_MaxisR'),'EdgeColor','none')
% shading interp
% view(2)
% hold on
% plot3([-10 10],[0 0],[10 10],':w')
% plot3([0 0],[-10 10],[10 10],':w')
contourf(dt_deg_MSE, da_MSE, log10(MSE_pitch_MaxisR'),50,'LineColor','none');
set(gca,'clim',[0 logMSEmax])
hold on
plot(phase_pitch_MaxisR_MSEmin,dpitch_MaxisR_MSEmin,'.w')
plot([-10 10],[0 0],':w')
plot([0 0],[-10 10],':w')
% title('axisR rotation')
% ylabel('kin angle')
axis equal
axis tight

subplot(3,2,6)
% surf(dt_deg_MSE, da_MSE, log10(MSE_dev_MaxisR'),'EdgeColor','none')
% shading interp
% view(2)
% hold on
% plot3([-10 10],[0 0],[10 10],':w')
% plot3([0 0],[-10 10],[10 10],':w')
contourf(dt_deg_MSE, da_MSE, log10(MSE_dev_MaxisR'),50,'LineColor','none');
set(gca,'clim',[0 logMSEmax])
hold on
plot(phase_dev_MaxisR_MSEmin,ddev_MaxisR_MSEmin,'.w')
plot([-10 10],[0 0],':w')
plot([0 0],[-10 10],':w')
% title('axisR dev')
% ylabel('kin angle')
xlabel('phase')
axis equal
axis tight
colorbar

saveas(gca,['MSplot_logMSE_WBkinMODmax_LnR_logMSEmax',num2str(logMSEmax),'.fig'])
saveas(gca,['MSplot_logMSE_WBkinMODmax_LnR_logMSEmax',num2str(logMSEmax),'.png'])
plot2svg(['MSplot_logMSE_WBkinMODmax_LnR_logMSEmax',num2str(logMSEmax),'.svg'])

cd ..


