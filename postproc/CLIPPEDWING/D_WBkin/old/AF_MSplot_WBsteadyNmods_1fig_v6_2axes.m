clc
clear
close all

name = 'WBdataset_all_steadyNmods_TorqueNorm.mat'
load(name)

load('norm_data_torque.mat')

%% mod locations

% % low inertia (wingdisc)
% Fenhanses = [0 .3]'
% rollTorques = [0 .04]'
% pitchTorques = [0 .02]'
% yawTorques = [0 .07]'
% RaxisTorques = [0 .04]'
% LaxisTorques = [0 .015]'
% 
% FenhMods = Fenhanses/Fenhance_norm
% RollMods = rollTorques/Mroll_norm
% PitchMods = pitchTorques/Mpitch_norm
% YawMods = yawTorques/Myaw_norm
% RaxisMods = RaxisTorques/M_R_norm
% LaxisMods = LaxisTorques/M_L_norm
% 
% % save norm data
% save('norm_data_torque_DiskInertia.mat','f_wb_steady',...
%     'Fenhance_norm','Mroll_norm','Mpitch_norm','Myaw_norm','M_R_norm','M_L_norm',...
%     'Fenhanses','rollTorques','pitchTorques','yawTorques','RaxisTorques','LaxisTorques',...
%     'FenhMods','RollMods','PitchMods','YawMods','RaxisMods','LaxisMods')

% high inertia (escape data)
Fenhanses = [0 .26]'
rollTorques = [0 .034]'
pitchTorques = [0 .018]'
yawTorques = [0 .064]'
RaxisTorques = [0 .036]'
LaxisTorques = [0 .013]'
axis1Torques = [0 .035]'
axis2Torques = [0 .035]'
axis1normalTorques = [0 .016]'
axis2normalTorques = [0 .015]'

FenhMods = Fenhanses/Fenhance_norm
RollMods = rollTorques/Mroll_norm
PitchMods = pitchTorques/Mpitch_norm
YawMods = yawTorques/Myaw_norm
RaxisMods = RaxisTorques/M_R_norm
LaxisMods = LaxisTorques/M_L_norm
axis1Mods = axis1Torques/M_axis1_norm
axis2Mods = axis2Torques/M_axis2_norm
axis1normalMods = axis1normalTorques/M_axis1normal_norm
axis2normalMods = axis2normalTorques/M_axis2normal_norm

% save norm data
save('norm_data_torque.mat','f_wb_steady',...
    'Fenhance_norm','Mroll_norm','Mpitch_norm','Myaw_norm','M_R_norm','M_L_norm',...
    'M_axis1_norm','M_axis2_norm','M_axis1normal_norm','M_axis2normal_norm',...
    'Fenhanses','rollTorques','pitchTorques','yawTorques','RaxisTorques','LaxisTorques',...
    'axis1Torques','axis2Torques','axis1normalTorques','axis2normalTorques',...
    'FenhMods','RollMods','PitchMods','YawMods','RaxisMods','LaxisMods',...
    'axis1Mods','axis2Mods','axis1normalMods','axis2normalMods')

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

%% steady WB
stroke_steady = calc_val_fourier_series_4thN8th_order(binx,stroke_steady_fourier_coeffs_binmean,0);
pitch_steady = calc_val_fourier_series_4thN8th_order(binx,pitch_steady_fourier_coeffs_binmean,0);
dev_steady = calc_val_fourier_series_4thN8th_order(binx,dev_steady_fourier_coeffs_binmean,0);

%% Force MOD wb
mod_values = FenhMods;

strokeMOD = calc_val_fourier_series_4thN8th_order(binx,strokeMOD_Fenhance_fourier_coeffs_binmean,0);
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,pitchMOD_Fenhance_fourier_coeffs_binmean,0);
devMOD = calc_val_fourier_series_4thN8th_order(binx,devMOD_Fenhance_fourier_coeffs_binmean,0);

color_map = [0 0 0; 1 .5 0];
% color_map = [0 0 0; 0 1 1];
% color_map = [0 1 0; 1 .5 0; 1 0 0];

for i=1:length(mod_values)
    k=mod_values(i)
    
    subplot(3,3,3)
    plot(binx,stroke_steady+k*strokeMOD,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,4)
    plot(binx,pitch_steady+k*pitchMOD-90,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,7)
    plot(binx,dev_steady+k*devMOD,'-','color',color_map(i,:),'linewidth',1)
end

%% Roll MOD Dwb
mod_values = RollMods;

strokeMOD = calc_val_fourier_series_4thN8th_order(binx,DstrokeMOD_RollTorque_fourier_coeffs_binmean,0);
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,DpitchMOD_RollTorque_fourier_coeffs_binmean,0);
devMOD = calc_val_fourier_series_4thN8th_order(binx,DdevMOD_RollTorque_fourier_coeffs_binmean,0);

color_map = [0 0 0; 1 .5 0];
% color_map = [0 0 0; 0 1 1];
% color_map = [0 1 0; 1 .5 0; 1 0 0];

for i=1:length(mod_values)
    k=mod_values(i)
    
    subplot(3,3,4)
    plot(binx,stroke_steady+k*strokeMOD,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,5)
    plot(binx,pitch_steady+k*pitchMOD-90,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,8)
    plot(binx,dev_steady+k*devMOD,'-','color',color_map(i,:),'linewidth',1)
end

%% Roll MOD wb
mod_values = RollMods;

% up wing
strokeMOD = calc_val_fourier_series_4thN8th_order(binx,strokeMOD_up_RollTorque_fourier_coeffs_binmean,0);
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,pitchMOD_up_RollTorque_fourier_coeffs_binmean,0);
devMOD = calc_val_fourier_series_4thN8th_order(binx,devMOD_up_RollTorque_fourier_coeffs_binmean,0);

color_map = [0 0 0; 1 .5 0];
% color_map = [0 0 0; 0 1 1];
% color_map = [0 1 0; 1 .5 0; 1 0 0];

for i=1:length(mod_values)
    k=mod_values(i)
    
    subplot(3,3,5)
    plot(binx,stroke_steady+k*strokeMOD,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,5)
    plot(binx,pitch_steady+k*pitchMOD-90,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,8)
    plot(binx,dev_steady+k*devMOD,'-','color',color_map(i,:),'linewidth',1)
end

% down wing
strokeMOD = calc_val_fourier_series_4thN8th_order(binx,strokeMOD_down_RollTorque_fourier_coeffs_binmean,0);
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,pitchMOD_down_RollTorque_fourier_coeffs_binmean,0);
devMOD = calc_val_fourier_series_4thN8th_order(binx,devMOD_down_RollTorque_fourier_coeffs_binmean,0);

color_map = [0 0 0; 0 .5 1];
% color_map = [0 0 0; 0 1 1];
% color_map = [0 1 0; 1 .5 0; 1 0 0];

for i=1:length(mod_values)
    k=mod_values(i)
    
    subplot(3,3,5)
    plot(binx,stroke_steady+k*strokeMOD,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,5)
    plot(binx,pitch_steady+k*pitchMOD-90,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,8)
    plot(binx,dev_steady+k*devMOD,'-','color',color_map(i,:),'linewidth',1)
end

%% Yaw MOD Dwb
mod_values = YawMods;

strokeMOD = calc_val_fourier_series_4thN8th_order(binx,DstrokeMOD_YawTorque_fourier_coeffs_binmean,0);
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,DpitchMOD_YawTorque_fourier_coeffs_binmean,0);
devMOD = calc_val_fourier_series_4thN8th_order(binx,DdevMOD_YawTorque_fourier_coeffs_binmean,0);

color_map = [0 0 0; 1 .5 0];
% color_map = [0 0 0; 0 1 1];
% color_map = [0 1 0; 1 .5 0; 1 0 0];

for i=1:length(mod_values)
    k=mod_values(i)
    
    subplot(3,3,1)
    plot(binx,stroke_steady+k*strokeMOD,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,5)
    plot(binx,pitch_steady+k*pitchMOD-90,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,8)
    plot(binx,dev_steady+k*devMOD,'-','color',color_map(i,:),'linewidth',1)
end

%% Yaw MOD wb
mod_values = YawMods;

% fwd wing
strokeMOD = calc_val_fourier_series_4thN8th_order(binx,strokeMOD_fwd_YawTorque_fourier_coeffs_binmean,0);
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,pitchMOD_fwd_YawTorque_fourier_coeffs_binmean,0);
devMOD = calc_val_fourier_series_4thN8th_order(binx,devMOD_fwd_YawTorque_fourier_coeffs_binmean,0);

color_map = [0 0 0; 1 .5 0];
% color_map = [0 0 0; 0 1 1];
% color_map = [0 1 0; 1 .5 0; 1 0 0];

for i=1:length(mod_values)
    k=mod_values(i)
    
    subplot(3,3,2)
    plot(binx,stroke_steady+k*strokeMOD,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,5)
    plot(binx,pitch_steady+k*pitchMOD-90,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,8)
    plot(binx,dev_steady+k*devMOD,'-','color',color_map(i,:),'linewidth',1)
end

% rwd wing
strokeMOD = calc_val_fourier_series_4thN8th_order(binx,strokeMOD_rwd_YawTorque_fourier_coeffs_binmean,0);
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,pitchMOD_rwd_YawTorque_fourier_coeffs_binmean,0);
devMOD = calc_val_fourier_series_4thN8th_order(binx,devMOD_rwd_YawTorque_fourier_coeffs_binmean,0);

color_map = [0 0 0; 0 .5 1];
% color_map = [0 0 0; 0 1 1];
% color_map = [0 1 0; 1 .5 0; 1 0 0];

for i=1:length(mod_values)
    k=mod_values(i)
    
    subplot(3,3,2)
    plot(binx,stroke_steady+k*strokeMOD,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,5)
    plot(binx,pitch_steady+k*pitchMOD-90,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,8)
    plot(binx,dev_steady+k*devMOD,'-','color',color_map(i,:),'linewidth',1)
end

%% Pitch MOD wb
mod_values = PitchMods;

strokeMOD = calc_val_fourier_series_4thN8th_order(binx,strokeMOD_PitchTorque_fourier_coeffs_binmean,0);
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,pitchMOD_PitchTorque_fourier_coeffs_binmean,0);
devMOD = calc_val_fourier_series_4thN8th_order(binx,devMOD_PitchTorque_fourier_coeffs_binmean,0);

%% pitch up
color_map = [0 0 0; 1 .5 0];
% color_map = [0 0 0; 0 1 1];
% color_map = [0 1 0; 1 0 0];


for i=1:length(mod_values)
    k=mod_values(i)
    
    subplot(3,3,6)
    plot(binx,stroke_steady+k*strokeMOD,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,6)
    plot(binx,pitch_steady+k*pitchMOD-90,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,9)
    plot(binx,dev_steady+k*devMOD,'-','color',color_map(i,:),'linewidth',1)
end

%% pitch down
color_map = [0 0 0; 0 .5 1];
% color_map = [0 0 0; 0 1 1];
% color_map = [0 1 0; 0 0 1];

for i=1:length(mod_values)
    k=mod_values(i)
    
    subplot(3,3,6)
    plot(binx,stroke_steady-k*strokeMOD,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,6)
    plot(binx,pitch_steady-k*pitchMOD-90,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,9)
    plot(binx,dev_steady-k*devMOD,'-','color',color_map(i,:),'linewidth',1)
end

%% axisR MOD Dwb
mod_values = RaxisMods;

strokeMOD = calc_val_fourier_series_4thN8th_order(binx,DstrokeMOD_TorqueAxisR_fourier_coeffs_binmean,0);
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,DpitchMOD_TorqueAxisR_fourier_coeffs_binmean,0);
devMOD = calc_val_fourier_series_4thN8th_order(binx,DdevMOD_TorqueAxisR_fourier_coeffs_binmean,0);

%% axisR up
color_map = [0 0 0; 1 .5 0];
% color_map = [0 0 0; 0 1 1];
% color_map = [0 1 0; 1 0 0];


for i=1:length(mod_values)
    k=mod_values(i)
    
    subplot(3,3,7)
    plot(binx,stroke_steady+k*strokeMOD,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,6)
    plot(binx,pitch_steady+k*pitchMOD-90,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,9)
    plot(binx,dev_steady+k*devMOD,'-','color',color_map(i,:),'linewidth',1)
end

%% axisR down
color_map = [0 0 0; 0 .5 1];
% color_map = [0 0 0; 0 1 1];
% color_map = [0 1 0; 0 0 1];

for i=1:length(mod_values)
    k=mod_values(i)
    
    subplot(3,3,7)
    plot(binx,stroke_steady-k*strokeMOD,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,6)
    plot(binx,pitch_steady-k*pitchMOD-90,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,9)
    plot(binx,dev_steady-k*devMOD,'-','color',color_map(i,:),'linewidth',1)
end

%% Raxis MOD wb UP
mod_values = RaxisMods;

% fwd wing
strokeMOD = calc_val_fourier_series_4thN8th_order(binx,strokeMOD_L_TorqueAxisR_fourier_coeffs_binmean,0);
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,pitchMOD_L_TorqueAxisR_fourier_coeffs_binmean,0);
devMOD = calc_val_fourier_series_4thN8th_order(binx,devMOD_L_TorqueAxisR_fourier_coeffs_binmean,0);

color_map = [0 0 0; 1 .5 0];
% color_map = [0 0 0; 0 1 1];
% color_map = [0 1 0; 1 .5 0; 1 0 0];

for i=1:length(mod_values)
    k=mod_values(i)
    
    subplot(3,3,8)
    plot(binx,stroke_steady+k*strokeMOD,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,5)
    plot(binx,pitch_steady+k*pitchMOD-90,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,8)
    plot(binx,dev_steady+k*devMOD,'-','color',color_map(i,:),'linewidth',1)
end

% rwd wing
strokeMOD = calc_val_fourier_series_4thN8th_order(binx,strokeMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,pitchMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);
devMOD = calc_val_fourier_series_4thN8th_order(binx,devMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);

color_map = [0 0 0; 0 .5 1];
% color_map = [0 0 0; 0 1 1];
% color_map = [0 1 0; 1 .5 0; 1 0 0];

for i=1:length(mod_values)
    k=mod_values(i)
    
    subplot(3,3,8)
    plot(binx,stroke_steady+k*strokeMOD,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,5)
    plot(binx,pitch_steady+k*pitchMOD-90,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,8)
    plot(binx,dev_steady+k*devMOD,'-','color',color_map(i,:),'linewidth',1)
end

%% Raxis MOD wb DOWN
mod_values = RaxisMods;

% Left wing
strokeMOD = calc_val_fourier_series_4thN8th_order(binx,strokeMOD_L_TorqueAxisR_fourier_coeffs_binmean,0);
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,pitchMOD_L_TorqueAxisR_fourier_coeffs_binmean,0);
devMOD = calc_val_fourier_series_4thN8th_order(binx,devMOD_L_TorqueAxisR_fourier_coeffs_binmean,0);

color_map = [0 0 0; 1 .5 0];
% color_map = [0 0 0; 0 1 1];
% color_map = [0 1 0; 1 .5 0; 1 0 0];

for i=1:length(mod_values)
    k=mod_values(i)
    
    subplot(3,3,9)
    plot(binx,stroke_steady-k*strokeMOD,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,5)
    plot(binx,pitch_steady-k*pitchMOD-90,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,8)
    plot(binx,dev_steady-k*devMOD,'-','color',color_map(i,:),'linewidth',1)
end

% Right wing
strokeMOD = calc_val_fourier_series_4thN8th_order(binx,strokeMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,pitchMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);
devMOD = calc_val_fourier_series_4thN8th_order(binx,devMOD_R_TorqueAxisR_fourier_coeffs_binmean,0);

color_map = [0 0 0; 0 .5 1];
% color_map = [0 0 0; 0 1 1];
% color_map = [0 1 0; 1 .5 0; 1 0 0];

for i=1:length(mod_values)
    k=mod_values(i)
    
    subplot(3,3,9)
    plot(binx,stroke_steady-k*strokeMOD,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,5)
    plot(binx,pitch_steady-k*pitchMOD-90,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,8)
    plot(binx,dev_steady-k*devMOD,'-','color',color_map(i,:),'linewidth',1)
end


%%
subplot(3,3,1)
xlabel('t* (-)','fontsize',10) 
ylabel('yaw Dkin','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(3,3,2)
xlabel('t* (-)','fontsize',10) 
ylabel('yaw kin L&R','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(3,3,3)
xlabel('t* (-)','fontsize',10) 
ylabel('F kin','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(3,3,4)
xlabel('t* (-)','fontsize',10) 
ylabel('roll Dkin','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 
    
subplot(3,3,5)
xlabel('t* (-)','fontsize',10) 
ylabel('roll kin L&R','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 
    
subplot(3,3,6)
xlabel('t* (-)','fontsize',10) 
ylabel('pitch up&down kin','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 
    
subplot(3,3,7)
xlabel('t* (-)','fontsize',10) 
ylabel('Raxis up&down Dkin','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 
    
subplot(3,3,8)
xlabel('t* (-)','fontsize',10) 
ylabel('Raxis up kin','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 
    
subplot(3,3,9)
xlabel('t* (-)','fontsize',10) 
ylabel('Raxis down kin','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 
    
%% save plot
mkdir('MSfigs_WBkin_torque')
cd('MSfigs_WBkin_torque')
    saveas(gca,['MSplot_WBsteady2customMods_1fig.fig'])
    saveas(gca,['MSplot_WBsteady2customMods_1fig.png'])
    plot2svg(['MSplot_WBsteady2customMods_1fig.svg'])
cd ..

