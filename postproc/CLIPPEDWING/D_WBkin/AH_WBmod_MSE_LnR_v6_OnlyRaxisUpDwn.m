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

Nt = 25 % Nx2
Na = 25 % Nx2

dt = Dt/Nt % norm time
dt_deg = 360*dt
da = Da/Na % deg

%% plot steady wb with dt & da
da_steady = 5;
dt_steady = 5;

% steady WB
stroke_steady = calc_val_fourier_series_4thN8th_order(binx,stroke_steady_fourier_coeffs_binmean,0);
pitch_steady  = calc_val_fourier_series_4thN8th_order(binx,pitch_steady_fourier_coeffs_binmean,0);
dev_steady    = calc_val_fourier_series_4thN8th_order(binx,dev_steady_fourier_coeffs_binmean,0);

stroke_steady_dt    =              calc_val_fourier_series_4thN8th_order(binx-dt_steady/360,stroke_steady_fourier_coeffs_binmean,0);
pitch_steady_dadt   = da_steady +  calc_val_fourier_series_4thN8th_order(binx-dt_steady/360,pitch_steady_fourier_coeffs_binmean,0);
dev_steady_da       = da_steady +  calc_val_fourier_series_4thN8th_order(binx,dev_steady_fourier_coeffs_binmean,0);

%% plots
mkdir('MSplot_WBkin_LnR_correlations')
cd('MSplot_WBkin_LnR_correlations')

%% steady WBs & WBs with Dt & Da
% figure
% subplot(3,3,1)
% hold on
% 
% plot(360*binx,stroke_steady,'-k','linewidth',1)
% plot(360*binx,stroke_steady_dt,'-','color',[.5 .5 .5],'linewidth',1)
% 
% plot(360*binx,pitch_steady-90,'-k','linewidth',1)
% plot(360*binx,pitch_steady_dadt-90,'-','color',[.5 .5 .5],'linewidth',1)
% 
% plot(360*binx,dev_steady,'-k','linewidth',1)
% plot(360*binx,dev_steady_da,'-','color',[.5 .5 .5],'linewidth',1)
% 
% plot([180 185],[0 5],'-k')
% xlabel('t* (-)','fontsize',10) 
% ylabel('kin','fontsize',10) 
% biny_min = -90;
% biny_max = 90;
% axis equal
% axis([0 360 biny_min biny_max])
% set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
%     set(gca,'XTick',0:180:360) 
% 
% saveas(gca,['MSplot_steadyWB_Dt',num2str(da_steady),'deg_Da',num2str(da_steady),'deg.fig'])
% saveas(gca,['MSplot_steadyWB_Dt',num2str(da_steady),'deg_Da',num2str(da_steady),'deg.png'])
% plot2svg(['MSplot_steadyWB_Dt',num2str(da_steady),'deg_Da',num2str(da_steady),'deg.svg'])

cd ..

%% CONSTANT WING
% steady WB
stroke_steady = calc_val_fourier_series_4thN8th_order(binx,stroke_steady_fourier_coeffs_binmean,0);
pitch_steady  = calc_val_fourier_series_4thN8th_order(binx,pitch_steady_fourier_coeffs_binmean,0);
dev_steady    = calc_val_fourier_series_4thN8th_order(binx,dev_steady_fourier_coeffs_binmean,0);

% fwd wing 
% stroke_Myaw_FWD = stroke_steady + max(YawMods)*calc_val_fourier_series_4thN8th_order(binx,strokeMOD_fwd_YawTorque_fourier_coeffs_binmean,0);
% pitch_Myaw_FWD  = pitch_steady  + max(YawMods)*calc_val_fourier_series_4thN8th_order(binx,pitchMOD_fwd_YawTorque_fourier_coeffs_binmean,0);
% dev_Myaw_FWD    = dev_steady    + max(YawMods)*calc_val_fourier_series_4thN8th_order(binx,devMOD_fwd_YawTorque_fourier_coeffs_binmean,0);

stroke_MaxisR_L = stroke_steady + max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx,strokeMOD_L_TorqueAxisR_fourier_coeffs_binmean,0);
pitch_MaxisR_L  = pitch_steady  + max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx,pitchMOD_L_TorqueAxisR_fourier_coeffs_binmean,0);
dev_MaxisR_L    = dev_steady    + max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx,devMOD_L_TorqueAxisR_fourier_coeffs_binmean,0);

stroke_MaxisRdwn_L = stroke_steady - max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx,strokeMOD_L_TorqueAxisR_fourier_coeffs_binmean,0);
pitch_MaxisRdwn_L  = pitch_steady  - max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx,pitchMOD_L_TorqueAxisR_fourier_coeffs_binmean,0);
dev_MaxisRdwn_L    = dev_steady    - max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx,devMOD_L_TorqueAxisR_fourier_coeffs_binmean,0);

% rwd wing 
% stroke_Myaw_RWD0 = stroke_steady + max(YawMods)*calc_val_fourier_series_4thN8th_order(binx,strokeMOD_rwd_YawTorque_fourier_coeffs_binmean,0);
% pitch_Myaw_RWD0  = pitch_steady  + max(YawMods)*calc_val_fourier_series_4thN8th_order(binx,pitchMOD_rwd_YawTorque_fourier_coeffs_binmean,0);
% dev_Myaw_RWD0   = dev_steady    + max(YawMods)*calc_val_fourier_series_4thN8th_order(binx,devMOD_rwd_YawTorque_fourier_coeffs_binmean,0);

stroke_MaxisR_R0 = stroke_steady + max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx,strokeMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);
pitch_MaxisR_R0  = pitch_steady  + max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx,pitchMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);
dev_MaxisR_R0    = dev_steady    + max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx,devMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);

stroke_MaxisRdwn_R0 = stroke_steady - max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx,strokeMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);
pitch_MaxisRdwn_R0  = pitch_steady  - max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx,pitchMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);
dev_MaxisRdwn_R0    = dev_steady    - max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx,devMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);

%% dt loop

% if matlabpool('size')==0
%     matlabpool open 8
% end
% parfor i=1:2*Nt+1

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
%     stroke_Myaw_RWD = stroke_steady_dt + max(YawMods)*calc_val_fourier_series_4thN8th_order(binx-dt_now,strokeMOD_rwd_YawTorque_fourier_coeffs_binmean,0);
%     pitch_Myaw_RWD  = pitch_steady_dt  + max(YawMods)*calc_val_fourier_series_4thN8th_order(binx-dt_now,pitchMOD_rwd_YawTorque_fourier_coeffs_binmean,0);
%     dev_Myaw_RWD    = dev_steady_dt    + max(YawMods)*calc_val_fourier_series_4thN8th_order(binx-dt_now,devMOD_rwd_YawTorque_fourier_coeffs_binmean,0);

    stroke_MaxisR_R = stroke_steady_dt + max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx-dt_now,strokeMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);
    pitch_MaxisR_R  = pitch_steady_dt  + max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx-dt_now,pitchMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);
    dev_MaxisR_R    = dev_steady_dt    + max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx-dt_now,devMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);

    stroke_MaxisRdwn_R = stroke_steady_dt - max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx-dt_now,strokeMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);
    pitch_MaxisRdwn_R  = pitch_steady_dt  - max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx-dt_now,pitchMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);
    dev_MaxisRdwn_R    = dev_steady_dt    - max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx-dt_now,devMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);

    
    % d angle loop
    for j=1:2*Na+1
        da_now = (j-Na-1)*da;
        da_MSE(j,1) = da_now;
        
        % WINGKIN DIFF
%         Dstroke_Myaw = stroke_Myaw_FWD - stroke_Myaw_RWD - da_now;
%         Dpitch_Myaw  = pitch_Myaw_FWD  - pitch_Myaw_RWD  - da_now;
%         Ddev_Myaw    = dev_Myaw_FWD    - dev_Myaw_RWD    - da_now;
        
        Dstroke_MaxisR = stroke_MaxisR_L - stroke_MaxisR_R - da_now;
        Dpitch_MaxisR  = pitch_MaxisR_L  - pitch_MaxisR_R  - da_now;
        Ddev_MaxisR    = dev_MaxisR_L    - dev_MaxisR_R    - da_now;
        
        Dstroke_MaxisRdwn = stroke_MaxisRdwn_L - stroke_MaxisRdwn_R - da_now;
        Dpitch_MaxisRdwn  = pitch_MaxisRdwn_L  - pitch_MaxisRdwn_R  - da_now;
        Ddev_MaxisRdwn    = dev_MaxisRdwn_L    - dev_MaxisRdwn_R    - da_now;
        
%         MSE_stroke_Myaw(i,j) = mse(Dstroke_Myaw);
%         MSE_pitch_Myaw(i,j) = mse(Dpitch_Myaw);
%         MSE_dev_Myaw(i,j) = mse(Ddev_Myaw);
        
        MSE_stroke_MaxisR(i,j) = mse(Dstroke_MaxisR);
        MSE_pitch_MaxisR(i,j)  = mse(Dpitch_MaxisR);
        MSE_dev_MaxisR(i,j)    = mse(Ddev_MaxisR);
        
        MSE_stroke_MaxisRdwn(i,j) = mse(Dstroke_MaxisRdwn);
        MSE_pitch_MaxisRdwn(i,j)  = mse(Dpitch_MaxisRdwn);
        MSE_dev_MaxisRdwn(i,j)    = mse(Ddev_MaxisRdwn);
    end
end

dt_deg_MSE = dt_MSE*360; % phase shift in deg;
[DA_MSE,DT_DEG_MSE] = meshgrid(da_MSE,dt_deg_MSE);
[DA_MSE,DT_MSE] = meshgrid(da_MSE,dt_MSE);

%% min MSE yaw
% MSE_stroke_Myaw_MSEmin   = min(MSE_stroke_Myaw(:));
% dt_stroke_Myaw_MSEmin    = DT_MSE(MSE_stroke_Myaw==MSE_stroke_Myaw_MSEmin);
% phase_stroke_Myaw_MSEmin = DT_DEG_MSE(MSE_stroke_Myaw==MSE_stroke_Myaw_MSEmin);
% dstroke_Myaw_MSEmin      = DA_MSE(MSE_stroke_Myaw==MSE_stroke_Myaw_MSEmin);
% 
% MSE_pitch_Myaw_MSEmin   = min(MSE_pitch_Myaw(:));
% dt_pitch_Myaw_MSEmin    = DT_MSE(MSE_pitch_Myaw==MSE_pitch_Myaw_MSEmin);
% phase_pitch_Myaw_MSEmin = DT_DEG_MSE(MSE_pitch_Myaw==MSE_pitch_Myaw_MSEmin);
% dpitch_Myaw_MSEmin      = DA_MSE(MSE_pitch_Myaw==MSE_pitch_Myaw_MSEmin);
% 
% MSE_dev_Myaw_MSEmin   = min(MSE_dev_Myaw(:));
% dt_dev_Myaw_MSEmin    = DT_MSE(MSE_dev_Myaw==MSE_dev_Myaw_MSEmin);
% phase_dev_Myaw_MSEmin = DT_DEG_MSE(MSE_dev_Myaw==MSE_dev_Myaw_MSEmin);
% ddev_Myaw_MSEmin      = DA_MSE(MSE_dev_Myaw==MSE_dev_Myaw_MSEmin);

%% min MSE axisR
MSE_stroke_MaxisR_MSEmin   = min(MSE_stroke_MaxisR(:));
dt_stroke_MaxisR_MSEmin    = DT_MSE(MSE_stroke_MaxisR==MSE_stroke_MaxisR_MSEmin);
phase_stroke_MaxisR_MSEmin = DT_DEG_MSE(MSE_stroke_MaxisR==MSE_stroke_MaxisR_MSEmin);
dstroke_MaxisR_MSEmin      = DA_MSE(MSE_stroke_MaxisR==MSE_stroke_MaxisR_MSEmin);

MSE_pitch_MaxisR_MSEmin   = min(MSE_pitch_MaxisR(:));
dt_pitch_MaxisR_MSEmin    = DT_MSE(MSE_pitch_MaxisR==MSE_pitch_MaxisR_MSEmin);
phase_pitch_MaxisR_MSEmin = DT_DEG_MSE(MSE_pitch_MaxisR==MSE_pitch_MaxisR_MSEmin);
dpitch_MaxisR_MSEmin      = DA_MSE(MSE_pitch_MaxisR==MSE_pitch_MaxisR_MSEmin);

MSE_dev_MaxisR_MSEmin   = min(MSE_dev_MaxisR(:));
dt_dev_MaxisR_MSEmin    = DT_MSE(MSE_dev_MaxisR==MSE_dev_MaxisR_MSEmin);
phase_dev_MaxisR_MSEmin = DT_DEG_MSE(MSE_dev_MaxisR==MSE_dev_MaxisR_MSEmin);
ddev_MaxisR_MSEmin      = DA_MSE(MSE_dev_MaxisR==MSE_dev_MaxisR_MSEmin);

%% min MSE axisRdwn
MSE_stroke_MaxisRdwn_MSEmin   = min(MSE_stroke_MaxisRdwn(:));
dt_stroke_MaxisRdwn_MSEmin    = DT_MSE(MSE_stroke_MaxisRdwn==MSE_stroke_MaxisRdwn_MSEmin);
phase_stroke_MaxisRdwn_MSEmin = DT_DEG_MSE(MSE_stroke_MaxisRdwn==MSE_stroke_MaxisRdwn_MSEmin);
dstroke_MaxisRdwn_MSEmin      = DA_MSE(MSE_stroke_MaxisRdwn==MSE_stroke_MaxisRdwn_MSEmin);

MSE_pitch_MaxisRdwn_MSEmin   = min(MSE_pitch_MaxisRdwn(:));
dt_pitch_MaxisRdwn_MSEmin    = DT_MSE(MSE_pitch_MaxisRdwn==MSE_pitch_MaxisRdwn_MSEmin);
phase_pitch_MaxisRdwn_MSEmin = DT_DEG_MSE(MSE_pitch_MaxisRdwn==MSE_pitch_MaxisRdwn_MSEmin);
dpitch_MaxisRdwn_MSEmin      = DA_MSE(MSE_pitch_MaxisRdwn==MSE_pitch_MaxisRdwn_MSEmin);

MSE_dev_MaxisRdwn_MSEmin   = min(MSE_dev_MaxisRdwn(:));
dt_dev_MaxisRdwn_MSEmin    = DT_MSE(MSE_dev_MaxisRdwn==MSE_dev_MaxisRdwn_MSEmin);
phase_dev_MaxisRdwn_MSEmin = DT_DEG_MSE(MSE_dev_MaxisRdwn==MSE_dev_MaxisRdwn_MSEmin);
ddev_MaxisRdwn_MSEmin      = DA_MSE(MSE_dev_MaxisRdwn==MSE_dev_MaxisRdwn_MSEmin);

%% WBs at MSEmin

% % Yaw
% stroke_steady_YawDt = calc_val_fourier_series_4thN8th_order(binx-dt_stroke_Myaw_MSEmin,stroke_steady_fourier_coeffs_binmean,0);
% pitch_steady_YawDt  = calc_val_fourier_series_4thN8th_order(binx-dt_pitch_Myaw_MSEmin,pitch_steady_fourier_coeffs_binmean,0);
% dev_steady_YawDt    = calc_val_fourier_series_4thN8th_order(binx-dt_dev_Myaw_MSEmin,dev_steady_fourier_coeffs_binmean,0);
% 
% stroke_Myaw_RWD_MSEmin = stroke_steady_YawDt + max(YawMods)*calc_val_fourier_series_4thN8th_order(binx-dt_stroke_Myaw_MSEmin,strokeMOD_rwd_YawTorque_fourier_coeffs_binmean,0);
% pitch_Myaw_RWD_MSEmin  = pitch_steady_YawDt  + max(YawMods)*calc_val_fourier_series_4thN8th_order(binx-dt_pitch_Myaw_MSEmin,pitchMOD_rwd_YawTorque_fourier_coeffs_binmean,0);
% dev_Myaw_RWD_MSEmin    = dev_steady_YawDt    + max(YawMods)*calc_val_fourier_series_4thN8th_order(binx-dt_dev_Myaw_MSEmin,devMOD_rwd_YawTorque_fourier_coeffs_binmean,0);
% 
% stroke_Myaw_RWD_MSEmin = stroke_Myaw_RWD_MSEmin + dstroke_Myaw_MSEmin;
% pitch_Myaw_RWD_MSEmin  = pitch_Myaw_RWD_MSEmin  + dpitch_Myaw_MSEmin;
% dev_Myaw_RWD_MSEmin    = dev_Myaw_RWD_MSEmin    + ddev_Myaw_MSEmin;

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

% axisRdwn
stroke_steady_axisRdwndt = calc_val_fourier_series_4thN8th_order(binx-dt_stroke_MaxisRdwn_MSEmin,stroke_steady_fourier_coeffs_binmean,0);
pitch_steady_axisRdwndt  = calc_val_fourier_series_4thN8th_order(binx-dt_pitch_MaxisRdwn_MSEmin,pitch_steady_fourier_coeffs_binmean,0);
dev_steady_axisRdwndt    = calc_val_fourier_series_4thN8th_order(binx-dt_dev_MaxisRdwn_MSEmin,dev_steady_fourier_coeffs_binmean,0);

stroke_MaxisRdwn_R_MSEmin = stroke_steady_axisRdwndt - max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx-dt_stroke_MaxisRdwn_MSEmin,strokeMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);
pitch_MaxisRdwn_R_MSEmin  = pitch_steady_axisRdwndt  - max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx-dt_pitch_MaxisRdwn_MSEmin,pitchMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);
dev_MaxisRdwn_R_MSEmin    = dev_steady_axisRdwndt    - max(RaxisMods)*calc_val_fourier_series_4thN8th_order(binx-dt_dev_MaxisRdwn_MSEmin,devMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);

stroke_MaxisRdwn_R_MSEmin = stroke_MaxisRdwn_R_MSEmin + dstroke_MaxisRdwn_MSEmin;
pitch_MaxisRdwn_R_MSEmin  = pitch_MaxisRdwn_R_MSEmin  + dpitch_MaxisRdwn_MSEmin;
dev_MaxisRdwn_R_MSEmin    = dev_MaxisRdwn_R_MSEmin    + ddev_MaxisRdwn_MSEmin;

%% save data
% save('WBkin_LnR_correlations_atMODmax.mat','binx',...
% 'dt_deg_MSE','dt_MSE','da_MSE','DA_MSE','DT_DEG_MSE','DT_MSE',...
% 'MSE_stroke_Myaw_MSEmin','MSE_stroke_Myaw','dt_stroke_Myaw_MSEmin','phase_stroke_Myaw_MSEmin','dstroke_Myaw_MSEmin',...
% 'MSE_pitch_Myaw_MSEmin','MSE_pitch_Myaw','dt_pitch_Myaw_MSEmin','phase_pitch_Myaw_MSEmin','dpitch_Myaw_MSEmin',...
% 'MSE_dev_Myaw_MSEmin','MSE_dev_Myaw','dt_dev_Myaw_MSEmin','phase_dev_Myaw_MSEmin','ddev_Myaw_MSEmin',...
% 'MSE_stroke_MaxisR_MSEmin','MSE_stroke_MaxisR','dt_stroke_MaxisR_MSEmin','phase_stroke_MaxisR_MSEmin','dstroke_MaxisR_MSEmin',...
% 'MSE_pitch_MaxisR_MSEmin','MSE_pitch_MaxisR','dt_pitch_MaxisR_MSEmin','phase_pitch_MaxisR_MSEmin','dpitch_MaxisR_MSEmin',...
% 'MSE_dev_MaxisR_MSEmin','MSE_dev_MaxisR','dt_dev_MaxisR_MSEmin','phase_dev_MaxisR_MSEmin','ddev_MaxisR_MSEmin',...
% 'MSE_stroke_MaxisRdwn_MSEmin','MSE_stroke_MaxisRdwn','dt_stroke_MaxisRdwn_MSEmin','phase_stroke_MaxisRdwn_MSEmin','dstroke_MaxisRdwn_MSEmin',...
% 'MSE_pitch_MaxisRdwn_MSEmin','MSE_pitch_MaxisRdwn','dt_pitch_MaxisRdwn_MSEmin','phase_pitch_MaxisRdwn_MSEmin','dpitch_MaxisRdwn_MSEmin',...
% 'MSE_dev_MaxisRdwn_MSEmin','MSE_dev_MaxisRdwn','dt_dev_MaxisRdwn_MSEmin','phase_dev_MaxisRdwn_MSEmin','ddev_MaxisRdwn_MSEmin',...
% 'stroke_Myaw_FWD','stroke_Myaw_RWD0','stroke_Myaw_RWD_MSEmin',...
% 'pitch_Myaw_FWD','pitch_Myaw_RWD0','pitch_Myaw_RWD_MSEmin',...
% 'dev_Myaw_FWD','dev_Myaw_RWD0','dev_Myaw_RWD_MSEmin',...
% 'stroke_MaxisR_L','stroke_MaxisR_R0','stroke_MaxisR_R_MSEmin',...
% 'pitch_MaxisR_L','pitch_MaxisR_R0','pitch_MaxisR_R_MSEmin',...
% 'dev_MaxisR_L','dev_MaxisR_R0','dev_MaxisR_R_MSEmin',...
% 'stroke_MaxisRdwn_L','stroke_MaxisRdwn_R0','stroke_MaxisRdwn_R_MSEmin',...
% 'pitch_MaxisRdwn_L','pitch_MaxisRdwn_R0','pitch_MaxisRdwn_R_MSEmin',...
% 'dev_MaxisRdwn_L','dev_MaxisRdwn_R0','dev_MaxisRdwn_R_MSEmin')

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

% % yaw @ MSEmin
% subplot(3,3,1)
% plot(binx,stroke_Myaw_FWD,'-','color',[1 .5 0],'linewidth',1)
% plot(binx,stroke_Myaw_RWD_MSEmin,'-','color',[0 .5 1],'linewidth',1)
% 
% plot(binx,pitch_Myaw_FWD-90,'-','color',[1 .5 0],'linewidth',1)
% plot(binx,pitch_Myaw_RWD_MSEmin-90,'-','color',[0 .5 1],'linewidth',1)
% 
% plot(binx,dev_Myaw_FWD,'-','color',[1 .5 0],'linewidth',1)
% plot(binx,dev_Myaw_RWD_MSEmin,'-','color',[0 .5 1],'linewidth',1)

% axisR @ MSEmin
subplot(3,3,2)
plot(binx,stroke_MaxisR_L,'-','color',[1 .5 0],'linewidth',1)
plot(binx,stroke_MaxisR_R_MSEmin,'-','color',[0 .5 1],'linewidth',1)

plot(binx,pitch_MaxisR_L-90,'-','color',[1 .5 0],'linewidth',1)
plot(binx,pitch_MaxisR_R_MSEmin-90,'-','color',[0 .5 1],'linewidth',1)

plot(binx,dev_MaxisR_L,'-','color',[1 .5 0],'linewidth',1)
plot(binx,dev_MaxisR_R_MSEmin,'-','color',[0 .5 1],'linewidth',1)

% axisRdwn @ MSEmin
subplot(3,3,3)
plot(binx,stroke_MaxisRdwn_L,'-','color',[1 .5 0],'linewidth',1)
plot(binx,stroke_MaxisRdwn_R_MSEmin,'-','color',[0 .5 1],'linewidth',1)

plot(binx,pitch_MaxisRdwn_L-90,'-','color',[1 .5 0],'linewidth',1)
plot(binx,pitch_MaxisRdwn_R_MSEmin-90,'-','color',[0 .5 1],'linewidth',1)

plot(binx,dev_MaxisRdwn_L,'-','color',[1 .5 0],'linewidth',1)
plot(binx,dev_MaxisRdwn_R_MSEmin,'-','color',[0 .5 1],'linewidth',1)

% % yaw original
% subplot(3,3,4)
% plot(binx,stroke_Myaw_FWD,'-','color',[1 .5 0],'linewidth',1)
% plot(binx,stroke_Myaw_RWD0,'-','color',[0 .5 1],'linewidth',1)
% 
% plot(binx,pitch_Myaw_FWD-90,'-','color',[1 .5 0],'linewidth',1)
% plot(binx,pitch_Myaw_RWD0-90,'-','color',[0 .5 1],'linewidth',1)
% 
% plot(binx,dev_Myaw_FWD,'-','color',[1 .5 0],'linewidth',1)
% plot(binx,dev_Myaw_RWD0,'-','color',[0 .5 1],'linewidth',1)

% axisR original
subplot(3,3,5)
plot(binx,stroke_MaxisR_L,'-','color',[1 .5 0],'linewidth',1)
plot(binx,stroke_MaxisR_R0,'-','color',[0 .5 1],'linewidth',1)

plot(binx,pitch_MaxisR_L-90,'-','color',[1 .5 0],'linewidth',1)
plot(binx,pitch_MaxisR_R0-90,'-','color',[0 .5 1],'linewidth',1)

plot(binx,dev_MaxisR_L,'-','color',[1 .5 0],'linewidth',1)
plot(binx,dev_MaxisR_R0,'-','color',[0 .5 1],'linewidth',1)

% axisRdwn original
subplot(3,3,6)
plot(binx,stroke_MaxisRdwn_L,'-','color',[1 .5 0],'linewidth',1)
plot(binx,stroke_MaxisRdwn_R0,'-','color',[0 .5 1],'linewidth',1)

plot(binx,pitch_MaxisRdwn_L-90,'-','color',[1 .5 0],'linewidth',1)
plot(binx,pitch_MaxisRdwn_R0-90,'-','color',[0 .5 1],'linewidth',1)

plot(binx,dev_MaxisRdwn_L,'-','color',[1 .5 0],'linewidth',1)
plot(binx,dev_MaxisRdwn_R0,'-','color',[0 .5 1],'linewidth',1)

% % yaw original & @ MSEmin
% subplot(3,3,7)
% plot(binx,stroke_Myaw_FWD,'-','color',[1 .5 0],'linewidth',1)
% plot(binx,stroke_Myaw_RWD_MSEmin,'-','color',[0 .5 1],'linewidth',1)
% plot(binx,stroke_Myaw_RWD0,'-','color',[.5 .5 .5],'linewidth',1)
% 
% plot(binx,pitch_Myaw_FWD-90,'-','color',[1 .5 0],'linewidth',1)
% plot(binx,pitch_Myaw_RWD_MSEmin-90,'-','color',[0 .5 1],'linewidth',1)
% plot(binx,pitch_Myaw_RWD0-90,'-','color',[.5 .5 .5],'linewidth',1)
% 
% plot(binx,dev_Myaw_FWD,'-','color',[1 .5 0],'linewidth',1)
% plot(binx,dev_Myaw_RWD_MSEmin,'-','color',[0 .5 1],'linewidth',1)
% plot(binx,dev_Myaw_RWD0,'-','color',[.5 .5 .5],'linewidth',1)

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

% axisRdwn original & @ MSEmin
subplot(3,3,9)
plot(binx,stroke_MaxisRdwn_L,'-','color',[1 .5 0],'linewidth',1)
plot(binx,stroke_MaxisRdwn_R_MSEmin,'-','color',[0 .5 1],'linewidth',1)
plot(binx,stroke_MaxisRdwn_R0,'-','color',[.5 .5 .5],'linewidth',1)

plot(binx,pitch_MaxisRdwn_L-90,'-','color',[1 .5 0],'linewidth',1)
plot(binx,pitch_MaxisRdwn_R_MSEmin-90,'-','color',[0 .5 1],'linewidth',1)
plot(binx,pitch_MaxisRdwn_R0-90,'-','color',[.5 .5 .5],'linewidth',1)

plot(binx,dev_MaxisRdwn_L,'-','color',[1 .5 0],'linewidth',1)
plot(binx,dev_MaxisRdwn_R_MSEmin,'-','color',[0 .5 1],'linewidth',1)
plot(binx,dev_MaxisRdwn_R0,'-','color',[.5 .5 .5],'linewidth',1)

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

subplot(3,3,3)
% xlabel('t* (-)','fontsize',10) 
ylabel('Raxis down kin @ MSEmin','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(3,3,4)
% xlabel('t* (-)','fontsize',10) 
ylabel('yaw kin @ MODmax','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(3,3,5)
% xlabel('t* (-)','fontsize',10) 
ylabel('Raxis up kin @ MODmax','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(3,3,6)
% xlabel('t* (-)','fontsize',10) 
ylabel('Raxis down kin @ MODmax','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(3,3,7)
xlabel('t* (-)','fontsize',10) 
ylabel('yaw kin @ MSEmin','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(3,3,8)
xlabel('t* (-)','fontsize',10) 
ylabel('Raxis up kin @ MSEmin','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(3,3,9)
xlabel('t* (-)','fontsize',10) 
ylabel('Raxis down kin @ MSEmin','fontsize',10) 
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
steps = 10;

figure
colormap gray

% subplot(3,3,1)
% % surf(dt_deg_MSE, da_MSE, (MSE_stroke_Myaw'),'EdgeColor','none')
% % shading interp
% % view(2)
% % hold on
% % plot3([-10 10],[0 0],[10 10],':w')
% % plot3([0 0],[-10 10],[10 10],':w')
% contourf(dt_deg_MSE, da_MSE, (MSE_stroke_Myaw'),steps,'LineColor','none');
% set(gca,'clim',[0 MSEmax])
% hold on
% plot(phase_stroke_Myaw_MSEmin,dstroke_Myaw_MSEmin,'.w')
% plot([-10 10],[0 0],':w')
% plot([0 0],[-10 10],':w')
% title('ln MSE yaw torque')
% ylabel('stroke kin angle')
% axis equal
% axis tight
% 
% subplot(3,3,4)
% % surf(dt_deg_MSE, da_MSE, (MSE_pitch_Myaw'),'EdgeColor','none')
% % shading interp
% % view(2)
% % hold on
% % plot3([-10 10],[0 0],[10 10],':w')
% % plot3([0 0],[-10 10],[10 10],':w')
% contourf(dt_deg_MSE, da_MSE, (MSE_pitch_Myaw'),steps,'LineColor','none');
% set(gca,'clim',[0 MSEmax])
% hold on
% plot(phase_pitch_Myaw_MSEmin,dpitch_Myaw_MSEmin,'.w')
% plot([-10 10],[0 0],':w')
% plot([0 0],[-10 10],':w')
% % title('yaw rotation')
% ylabel('rotation kin angle')
% axis equal
% axis tight
% 
% subplot(3,3,7)
% % surf(dt_deg_MSE, da_MSE, (MSE_dev_Myaw'),'EdgeColor','none')
% % shading interp
% % view(2)
% % hold on
% % plot3([-10 10],[0 0],[10 10],':w')
% % plot3([0 0],[-10 10],[10 10],':w')
% contourf(dt_deg_MSE, da_MSE, (MSE_dev_Myaw'),steps,'LineColor','none');
% set(gca,'clim',[0 MSEmax])
% hold on
% plot(phase_dev_Myaw_MSEmin,ddev_Myaw_MSEmin,'.w')
% plot([-10 10],[0 0],':w')
% plot([0 0],[-10 10],':w')
% % title('yaw dev')
% ylabel('deviation kin angle')
% xlabel('phase')
% axis equal
% axis tight

subplot(3,3,2)
% surf(dt_deg_MSE, da_MSE, (MSE_stroke_MaxisR'),'EdgeColor','none')
% shading interp
% view(2)
% hold on
% plot3([-10 10],[0 0],[10 10],':w')
% plot3([0 0],[-10 10],[10 10],':w')
contourf(dt_deg_MSE, da_MSE, (MSE_stroke_MaxisR'),steps,'LineColor','none');
set(gca,'clim',[0 MSEmax])
hold on
plot(phase_stroke_MaxisR_MSEmin,dstroke_MaxisR_MSEmin,'.w')
plot([-10 10],[0 0],':w')
plot([0 0],[-10 10],':w')
title('ln MSE axisR torque')
% ylabel('kin angle')
axis equal
axis tight

subplot(3,3,5)
% surf(dt_deg_MSE, da_MSE, (MSE_pitch_MaxisR'),'EdgeColor','none')
% shading interp
% view(2)
% hold on
% plot3([-10 10],[0 0],[10 10],':w')
% plot3([0 0],[-10 10],[10 10],':w')
contourf(dt_deg_MSE, da_MSE, (MSE_pitch_MaxisR'),steps,'LineColor','none');
set(gca,'clim',[0 MSEmax])
hold on
plot(phase_pitch_MaxisR_MSEmin,dpitch_MaxisR_MSEmin,'.w')
plot([-10 10],[0 0],':w')
plot([0 0],[-10 10],':w')
% title('axisR rotation')
% ylabel('kin angle')
axis equal
axis tight

subplot(3,3,8)
% surf(dt_deg_MSE, da_MSE, (MSE_dev_MaxisR'),'EdgeColor','none')
% shading interp
% view(2)
% hold on
% plot3([-10 10],[0 0],[10 10],':w')
% plot3([0 0],[-10 10],[10 10],':w')
contourf(dt_deg_MSE, da_MSE, (MSE_dev_MaxisR'),steps,'LineColor','none');
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

subplot(3,3,3)
% surf(dt_deg_MSE, da_MSE, (MSE_stroke_MaxisRdwn'),'EdgeColor','none')
% shading interp
% view(2)
% hold on
% plot3([-10 10],[0 0],[10 10],':w')
% plot3([0 0],[-10 10],[10 10],':w')
contourf(dt_deg_MSE, da_MSE, (MSE_stroke_MaxisRdwn'),steps,'LineColor','none');
set(gca,'clim',[0 MSEmax])
hold on
plot(phase_stroke_MaxisRdwn_MSEmin,dstroke_MaxisRdwn_MSEmin,'.w')
plot([-10 10],[0 0],':w')
plot([0 0],[-10 10],':w')
title('ln MSE axisRdwn torque')
% ylabel('kin angle')
axis equal
axis tight

subplot(3,3,6)
% surf(dt_deg_MSE, da_MSE, (MSE_pitch_MaxisRdwn'),'EdgeColor','none')
% shading interp
% view(2)
% hold on
% plot3([-10 10],[0 0],[10 10],':w')
% plot3([0 0],[-10 10],[10 10],':w')
contourf(dt_deg_MSE, da_MSE, (MSE_pitch_MaxisRdwn'),steps,'LineColor','none');
set(gca,'clim',[0 MSEmax])
hold on
plot(phase_pitch_MaxisRdwn_MSEmin,dpitch_MaxisRdwn_MSEmin,'.w')
plot([-10 10],[0 0],':w')
plot([0 0],[-10 10],':w')
% title('axisRdwn rotation')
% ylabel('kin angle')
axis equal
axis tight

subplot(3,3,9)
% surf(dt_deg_MSE, da_MSE, (MSE_dev_MaxisRdwn'),'EdgeColor','none')
% shading interp
% view(2)
% hold on
% plot3([-10 10],[0 0],[10 10],':w')
% plot3([0 0],[-10 10],[10 10],':w')
contourf(dt_deg_MSE, da_MSE, (MSE_dev_MaxisRdwn'),steps,'LineColor','none');
set(gca,'clim',[0 MSEmax])
hold on
plot(phase_dev_MaxisRdwn_MSEmin,ddev_MaxisRdwn_MSEmin,'.w')
plot([-10 10],[0 0],':w')
plot([0 0],[-10 10],':w')
% title('axisRdwn dev')
% ylabel('kin angle')
xlabel('phase')
axis equal
axis tight

colorbar
saveas(gca,['MSplot_MSE_WBkinMODmax_LnR_MSEmax',num2str(MSEmax),'.fig'])
saveas(gca,['MSplot_MSE_WBkinMODmax_LnR_MSEmax',num2str(MSEmax),'.png'])
plot2svg(['MSplot_MSE_WBkinMODmax_LnR_MSEmax',num2str(MSEmax),'.svg'])

%% log MSE plot
logMSEmax = 2;
steps = 8;

figure
colormap gray

% subplot(3,3,1)
% % surf(dt_deg_MSE, da_MSE, log10(MSE_stroke_Myaw'),'EdgeColor','none')
% % shading interp
% % view(2)
% % hold on
% % plot3([-10 10],[0 0],[10 10],':w')
% % plot3([0 0],[-10 10],[10 10],':w')
% contourf(dt_deg_MSE, da_MSE, log10(MSE_stroke_Myaw'),steps,'LineColor','none');
% set(gca,'clim',[0 logMSEmax])
% hold on
% plot(phase_stroke_Myaw_MSEmin,dstroke_Myaw_MSEmin,'.w')
% plot([-10 10],[0 0],':w')
% plot([0 0],[-10 10],':w')
% title('ln MSE yaw torque')
% ylabel('stroke kin angle')
% axis equal
% axis tight
% 
% subplot(3,3,4)
% % surf(dt_deg_MSE, da_MSE, log10(MSE_pitch_Myaw'),'EdgeColor','none')
% % shading interp
% % view(2)
% % hold on
% % plot3([-10 10],[0 0],[10 10],':w')
% % plot3([0 0],[-10 10],[10 10],':w')
% contourf(dt_deg_MSE, da_MSE, log10(MSE_pitch_Myaw'),steps,'LineColor','none');
% set(gca,'clim',[0 logMSEmax])
% hold on
% plot(phase_pitch_Myaw_MSEmin,dpitch_Myaw_MSEmin,'.w')
% plot([-10 10],[0 0],':w')
% plot([0 0],[-10 10],':w')
% % title('yaw rotation')
% ylabel('rotation kin angle')
% axis equal
% axis tight
% 
% subplot(3,3,7)
% % surf(dt_deg_MSE, da_MSE, log10(MSE_dev_Myaw'),'EdgeColor','none')
% % shading interp
% % view(2)
% % hold on
% % plot3([-10 10],[0 0],[10 10],':w')
% % plot3([0 0],[-10 10],[10 10],':w')
% contourf(dt_deg_MSE, da_MSE, log10(MSE_dev_Myaw'),steps,'LineColor','none');
% set(gca,'clim',[0 logMSEmax])
% hold on
% plot(phase_dev_Myaw_MSEmin,ddev_Myaw_MSEmin,'.w')
% plot([-10 10],[0 0],':w')
% plot([0 0],[-10 10],':w')
% % title('yaw dev')
% ylabel('deviation kin angle')
% xlabel('phase')
% axis equal
% axis tight

subplot(3,3,2)
% surf(dt_deg_MSE, da_MSE, log10(MSE_stroke_MaxisR'),'EdgeColor','none')
% shading interp
% view(2)
% hold on
% plot3([-10 10],[0 0],[10 10],':w')
% plot3([0 0],[-10 10],[10 10],':w')
contourf(dt_deg_MSE, da_MSE, log10(MSE_stroke_MaxisR'),steps,'LineColor','none');
set(gca,'clim',[0 logMSEmax])
hold on
plot(phase_stroke_MaxisR_MSEmin,dstroke_MaxisR_MSEmin,'.w')
plot([-10 10],[0 0],':w')
plot([0 0],[-10 10],':w')
title('ln MSE axisR torque')
% ylabel('kin angle')
axis equal
axis tight

subplot(3,3,5)
% surf(dt_deg_MSE, da_MSE, log10(MSE_pitch_MaxisR'),'EdgeColor','none')
% shading interp
% view(2)
% hold on
% plot3([-10 10],[0 0],[10 10],':w')
% plot3([0 0],[-10 10],[10 10],':w')
contourf(dt_deg_MSE, da_MSE, log10(MSE_pitch_MaxisR'),steps,'LineColor','none');
set(gca,'clim',[0 logMSEmax])
hold on
plot(phase_pitch_MaxisR_MSEmin,dpitch_MaxisR_MSEmin,'.w')
plot([-10 10],[0 0],':w')
plot([0 0],[-10 10],':w')
% title('axisR rotation')
% ylabel('kin angle')
axis equal
axis tight

subplot(3,3,8)
% surf(dt_deg_MSE, da_MSE, log10(MSE_dev_MaxisR'),'EdgeColor','none')
% shading interp
% view(2)
% hold on
% plot3([-10 10],[0 0],[10 10],':w')
% plot3([0 0],[-10 10],[10 10],':w')
contourf(dt_deg_MSE, da_MSE, log10(MSE_dev_MaxisR'),steps,'LineColor','none');
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

subplot(3,3,3)
% surf(dt_deg_MSE, da_MSE, log10(MSE_stroke_MaxisRdwn'),'EdgeColor','none')
% shading interp
% view(2)
% hold on
% plot3([-10 10],[0 0],[10 10],':w')
% plot3([0 0],[-10 10],[10 10],':w')
contourf(dt_deg_MSE, da_MSE, log10(MSE_stroke_MaxisRdwn'),steps,'LineColor','none');
set(gca,'clim',[0 logMSEmax])
hold on
plot(phase_stroke_MaxisRdwn_MSEmin,dstroke_MaxisRdwn_MSEmin,'.w')
plot([-10 10],[0 0],':w')
plot([0 0],[-10 10],':w')
title('ln MSE axisRdwn torque')
% ylabel('kin angle')
axis equal
axis tight

subplot(3,3,6)
% surf(dt_deg_MSE, da_MSE, log10(MSE_pitch_MaxisRdwn'),'EdgeColor','none')
% shading interp
% view(2)
% hold on
% plot3([-10 10],[0 0],[10 10],':w')
% plot3([0 0],[-10 10],[10 10],':w')
contourf(dt_deg_MSE, da_MSE, log10(MSE_pitch_MaxisRdwn'),steps,'LineColor','none');
set(gca,'clim',[0 logMSEmax])
hold on
plot(phase_pitch_MaxisRdwn_MSEmin,dpitch_MaxisRdwn_MSEmin,'.w')
plot([-10 10],[0 0],':w')
plot([0 0],[-10 10],':w')
% title('axisRdwn rotation')
% ylabel('kin angle')
axis equal
axis tight

subplot(3,3,9)
% surf(dt_deg_MSE, da_MSE, log10(MSE_dev_MaxisRdwn'),'EdgeColor','none')
% shading interp
% view(2)
% hold on
% plot3([-10 10],[0 0],[10 10],':w')
% plot3([0 0],[-10 10],[10 10],':w')
contourf(dt_deg_MSE, da_MSE, log10(MSE_dev_MaxisRdwn'),steps,'LineColor','none');
set(gca,'clim',[0 logMSEmax])
hold on
plot(phase_dev_MaxisRdwn_MSEmin,ddev_MaxisRdwn_MSEmin,'.w')
plot([-10 10],[0 0],':w')
plot([0 0],[-10 10],':w')
% title('axisRdwn dev')
% ylabel('kin angle')
xlabel('phase')
axis equal
axis tight

colorbar
saveas(gca,['MSplot_logMSE_WBkinMODmax_LnR_logMSEmax',num2str(logMSEmax),'.fig'])
saveas(gca,['MSplot_logMSE_WBkinMODmax_LnR_logMSEmax',num2str(logMSEmax),'.png'])
plot2svg(['MSplot_logMSE_WBkinMODmax_LnR_logMSEmax',num2str(logMSEmax),'.svg'])

cd ..


