clc
% clear
close all

name = 'WBdataset_all_steadyNmods_TorqueNorm.mat'
load(name)

% mod locations
Fenhanses = [0 .6]'
rollaccels = [0 3]'
pitchaccels = [0 2]'
yawaccels = [0 2]'

FenhMods = Fenhanses/Fenhance_norm
RollMods = rollaccels/rollaccel_norm
PitchMods = pitchaccels/pitchaccel_norm
YawMods = yawaccels/yawaccel_norm

% save norm data
save('norm_data.mat','rollaccel_norm','pitchaccel_norm','yawaccel_norm','Fenhance_norm',...
    'f_wb_steady','Fenhanses','rollaccels','pitchaccels','yawaccels',...
    'FenhMods','RollMods','PitchMods','YawMods')


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


%% steady WB
stroke_steady = calc_val_fourier_series_4thN8th_order(binx,stroke_steady_fourier_coeffs_binmean,0);
pitch_steady = calc_val_fourier_series_4thN8th_order(binx,pitch_steady_fourier_coeffs_binmean,0);
dev_steady = calc_val_fourier_series_4thN8th_order(binx,dev_steady_fourier_coeffs_binmean,0);


%% Force MOD wb
mod_values = FenhMods;

strokeMOD = calc_val_fourier_series_4thN8th_order(binx,strokeMOD_Fenhance_fourier_coeffs_binmean,0);
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,pitchMOD_Fenhance_fourier_coeffs_binmean,0);
devMOD = calc_val_fourier_series_4thN8th_order(binx,devMOD_Fenhance_fourier_coeffs_binmean,0);

color_map = [0 0 0; 0 1 1];
% color_map = [0 1 0; 1 .5 0; 1 0 0];

for i=1:length(mod_values)
    k=mod_values(i)
    
    subplot(3,3,1)
    plot(binx,stroke_steady+k*strokeMOD,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,4)
    plot(binx,pitch_steady+k*pitchMOD-90,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,7)
    plot(binx,dev_steady+k*devMOD,'-','color',color_map(i,:),'linewidth',1)
end

%% Roll MOD Dwb
mod_values = RollMods;

strokeMOD = calc_val_fourier_series_4thN8th_order(binx,DstrokeMOD_RollAccel_fourier_coeffs_binmean,0);
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,DpitchMOD_RollAccel_fourier_coeffs_binmean,0);
devMOD = calc_val_fourier_series_4thN8th_order(binx,DdevMOD_RollAccel_fourier_coeffs_binmean,0);

color_map = [0 0 0; 0 1 1];
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

%% Roll MOD wb
mod_values = RollMods;

% up wing
strokeMOD = calc_val_fourier_series_4thN8th_order(binx,strokeMOD_up_RollAccel_fourier_coeffs_binmean,0);
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,pitchMOD_up_RollAccel_fourier_coeffs_binmean,0);
devMOD = calc_val_fourier_series_4thN8th_order(binx,devMOD_up_RollAccel_fourier_coeffs_binmean,0);

color_map = [0 0 0; 0 1 1];
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
strokeMOD = calc_val_fourier_series_4thN8th_order(binx,strokeMOD_down_RollAccel_fourier_coeffs_binmean,0);
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,pitchMOD_down_RollAccel_fourier_coeffs_binmean,0);
devMOD = calc_val_fourier_series_4thN8th_order(binx,devMOD_down_RollAccel_fourier_coeffs_binmean,0);

color_map = [0 0 0; 1 0 0];
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

strokeMOD = calc_val_fourier_series_4thN8th_order(binx,DstrokeMOD_YawAccel_fourier_coeffs_binmean,0);
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,DpitchMOD_YawAccel_fourier_coeffs_binmean,0);
devMOD = calc_val_fourier_series_4thN8th_order(binx,DdevMOD_YawAccel_fourier_coeffs_binmean,0);

color_map = [0 0 0; 0 1 1];
% color_map = [0 1 0; 1 .5 0; 1 0 0];

for i=1:length(mod_values)
    k=mod_values(i)
    
    subplot(3,3,3)
    plot(binx,stroke_steady+k*strokeMOD,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,5)
    plot(binx,pitch_steady+k*pitchMOD-90,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,8)
    plot(binx,dev_steady+k*devMOD,'-','color',color_map(i,:),'linewidth',1)
end

%% Yaw MOD wb
mod_values = YawMods;

% fwd wing
strokeMOD = calc_val_fourier_series_4thN8th_order(binx,strokeMOD_fwd_YawAccel_fourier_coeffs_binmean,0);
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,pitchMOD_fwd_YawAccel_fourier_coeffs_binmean,0);
devMOD = calc_val_fourier_series_4thN8th_order(binx,devMOD_fwd_YawAccel_fourier_coeffs_binmean,0);

color_map = [0 0 0; 0 1 1];
% color_map = [0 1 0; 1 .5 0; 1 0 0];

for i=1:length(mod_values)
    k=mod_values(i)
    
    subplot(3,3,6)
    plot(binx,stroke_steady+k*strokeMOD,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,5)
    plot(binx,pitch_steady+k*pitchMOD-90,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,8)
    plot(binx,dev_steady+k*devMOD,'-','color',color_map(i,:),'linewidth',1)
end

% rwd wing
strokeMOD = calc_val_fourier_series_4thN8th_order(binx,strokeMOD_rwd_YawAccel_fourier_coeffs_binmean,0);
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,pitchMOD_rwd_YawAccel_fourier_coeffs_binmean,0);
devMOD = calc_val_fourier_series_4thN8th_order(binx,devMOD_rwd_YawAccel_fourier_coeffs_binmean,0);

color_map = [0 0 0; 1 0 0];
% color_map = [0 1 0; 1 .5 0; 1 0 0];

for i=1:length(mod_values)
    k=mod_values(i)
    
    subplot(3,3,6)
    plot(binx,stroke_steady+k*strokeMOD,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,5)
    plot(binx,pitch_steady+k*pitchMOD-90,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,8)
    plot(binx,dev_steady+k*devMOD,'-','color',color_map(i,:),'linewidth',1)
end

%% Pitch MOD wb
mod_values = PitchMods;

strokeMOD = calc_val_fourier_series_4thN8th_order(binx,strokeMOD_PitchAccel_fourier_coeffs_binmean,0);
pitchMOD = calc_val_fourier_series_4thN8th_order(binx,pitchMOD_PitchAccel_fourier_coeffs_binmean,0);
devMOD = calc_val_fourier_series_4thN8th_order(binx,devMOD_PitchAccel_fourier_coeffs_binmean,0);

%% pitch up
color_map = [0 0 0; 0 1 1];
% color_map = [0 1 0; 1 0 0];


for i=1:length(mod_values)
    k=mod_values(i)
    
    subplot(3,3,4)
    plot(binx,stroke_steady+k*strokeMOD,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,6)
    plot(binx,pitch_steady+k*pitchMOD-90,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,9)
    plot(binx,dev_steady+k*devMOD,'-','color',color_map(i,:),'linewidth',1)
end

%% pitch down
color_map = [0 0 0; 1 0 0];
% color_map = [0 1 0; 0 0 1];

for i=1:length(mod_values)
    k=mod_values(i)
    
    subplot(3,3,4)
    plot(binx,stroke_steady-k*strokeMOD,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,6)
    plot(binx,pitch_steady-k*pitchMOD-90,'-','color',color_map(i,:),'linewidth',1)
%     subplot(3,3,9)
    plot(binx,dev_steady-k*devMOD,'-','color',color_map(i,:),'linewidth',1)
end


%%
subplot(3,3,1)
xlabel('t* (-)','fontsize',10) 
ylabel('kin','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(3,3,2)
xlabel('t* (-)','fontsize',10) 
ylabel('kin','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(3,3,3)
xlabel('t* (-)','fontsize',10) 
ylabel('kin','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 

subplot(3,3,4)
xlabel('t* (-)','fontsize',10) 
ylabel('kin','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 
    
subplot(3,3,5)
xlabel('t* (-)','fontsize',10) 
ylabel('kin','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 
    
subplot(3,3,6)
xlabel('t* (-)','fontsize',10) 
ylabel('kin','fontsize',10) 
biny_min = -90;
biny_max = 90;
axis([0 1 biny_min biny_max])
set(gca,'YTick',[biny_min:(biny_max-biny_min)/2:biny_max],'fontsize',8) 
    set(gca,'XTick',0:.5:1) 
    
%% save plot
mkdir('MSfigs_WBkin')
cd('MSfigs_WBkin')
    saveas(gca,['MSplot_WBsteady2customMods_1fig.fig'])
    saveas(gca,['MSplot_WBsteady2customMods_1fig.png'])
    plot2svg(['MSplot_WBsteady2customMods_1fig.svg'])
cd ..

