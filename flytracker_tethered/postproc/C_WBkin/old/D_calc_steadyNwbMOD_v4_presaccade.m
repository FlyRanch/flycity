clc
clear
close all

addpath('/home/florian/Dropbox/WORK/toolbox')
addpath('/home/florian/Dropbox/WORK/toolbox/plot2svg')
addpath('/home/florian/Dropbox/WORK/toolbox/CircStat')
addpath('/home/florian/Dropbox/WORK/toolbox/flytracker')

loadname='WBdataset_all_3845WBs.mat'
loadname='WBdataset_all_1186WBs.mat'

steady_name = 'WBdataset_steady_1807WBs.mat'

calc_WB_STEADY = 0;
calc_WB_Fenhanced = 0;
calc_WB_PitchAccel = 0;
calc_WB_RollAccel = 0;
calc_WB_YawAccel = 0;

% calc_WB_STEADY = 1
calc_WB_Fenhanced = 1
calc_WB_PitchAccel = 1
calc_WB_RollAccel = 1
calc_WB_YawAccel = 1

mkdir('WBmod_figs')
load(loadname)
if calc_WB_STEADY == 0
if exist(steady_name) == 2
    load(steady_name)
end
end

%% constants
% damping effect

% inertia based on scaling 2^4
Iroll = 16*.1145e-12;
Ipitch = 16*.5060e-12;
Iyaw = 16*.4971e-12;

Croll = 29.21e-12;
Cpitch = 29.21e-12;
Cyaw = 29.21e-12;

k_damp_roll = Croll/Iroll;
k_damp_pitch = Cpitch/Ipitch;
k_damp_yaw = Cyaw/Iyaw;


%% settings
% limits (*std)
limit_steady = .5
limit_mod = .5
norm = 3

% fourier orders
    stroke_fourier_order = 4;
    pitch_fourier_order = 8;
    dev_fourier_order = 8;

    MOD_fourier_order = 8;

% number of polinomials
    n_pol_stroke = 12; % Order of used polynomials
    n_pol_pitch = 14; % Order of used polynomials
    n_pol_dev = 10; % Order of used polynomials 
    
    n_pol_MOD = 10;

plot_fits = 1;

n=200; % bins

linewidth_timelines = .5;
linewidth_meanWBs = 1;
skip = 50;

% heatmap resolution
nx = 1000;
ny = 100;

cmap_180 = jet(180);

% polyfit & 95% cof int settings
order = 3;
dn=20   % datapoints in bin
dm=20   % bin shift

color_code_now = [.5 .5 .5];
color_mid = [.25 .25 .25];
cmap = abs(colormap(gray)-1);
close all

%% calc & plot steady wb
% STEADY WB
if calc_WB_STEADY == 1
    
    calc_wbNORM_wingatt_steady
    
    cd('WBmod_figs')
%     plot_wbNORM_wingatt_steady
%     saveas(gca,['WBnorm_wingkin_',num2str(n_steady),'steadyWBs_limit',num2str(limit),'xSTD.fig'])
%     saveas(gca,['WBnorm_wingkin_',num2str(n_steady),'steadyWBs_limit',num2str(limit),'xSTD.png'])
%     plot2svg(['WBnorm_wingkin_',num2str(n_steady),'steadyWBs_limit',num2str(limit),'xSTD.svg'])
    
    plot_WBfunc_heatmap_steady_fourier
    saveas(gca,['WBnorm_wingkin_',num2str(n_steady),'steadyWBs_limit',num2str(limit_steady),'xSTD.fig'])
    saveas(gca,['WBnorm_wingkin_',num2str(n_steady),'steadyWBs_limit',num2str(limit_steady),'xSTD.png'])
    plot2svg(['WBnorm_wingkin_',num2str(n_steady),'steadyWBs_limit',num2str(limit_steady),'xSTD.svg'])
    
    plot_WBfreqNpitchV_hist_steady_LnR
    saveas(gca,['WBfreqNpitchV_wingkin_',num2str(n_steady),'steadyWBs_limit',num2str(limit_steady),'xSTD.fig'])
    saveas(gca,['WBfreqNpitchV_wingkin_',num2str(n_steady),'steadyWBs_limit',num2str(limit_steady),'xSTD.png'])
    plot2svg(['WBfreqNpitchV_wingkin_',num2str(n_steady),'steadyWBs_limit',num2str(limit_steady),'xSTD.svg'])
    cd ..
end

%% calc & plot MOD wb

% steady wb data
f_wb_steady = f_wb_steady_meanCIstd(:,1);

stroke_wb_steady = stroke_wb_steady_bins_meanCIstd(:,1);
stroke_ds_steady = stroke_ds_steady_bins_meanCIstd(:,1);
stroke_us_steady = stroke_us_steady_bins_meanCIstd(:,1);

pitch_wb_steady = pitch_wb_steady_bins_meanCIstd(:,1);
pitch_ds_steady = pitch_ds_steady_bins_meanCIstd(:,1);
pitch_us_steady = pitch_us_steady_bins_meanCIstd(:,1);

dev_wb_steady = dev_wb_steady_bins_meanCIstd(:,1);
dev_ds_steady = dev_ds_steady_bins_meanCIstd(:,1);
dev_us_steady = dev_us_steady_bins_meanCIstd(:,1);

% limits
rollaccel_limit_steady = limit_steady*nanstd(roll_dot_dot_mean_wb(:));
pitchaccel_limit_steady = limit_steady*nanstd(pitch_dot_dot_mean_wb(:));
yawaccel_limit_steady = limit_steady*nanstd(yaw_dot_dot_mean_wb(:));
Fenhance_limit_steady = limit_steady*nanstd(F_mean_wb(:));

rollaccel_limit_mod = limit_mod*nanstd(roll_dot_dot_mean_wb(:));
pitchaccel_limit_mod = limit_mod*nanstd(pitch_dot_dot_mean_wb(:));
yawaccel_limit_mod = limit_mod*nanstd(yaw_dot_dot_mean_wb(:));
Fenhance_limit_mod = limit_mod*nanstd(F_mean_wb(:));

% normalization value
rollaccel_norm = norm*nanstd(roll_dot_dot_mean_wb(:));
pitchaccel_norm = norm*nanstd(pitch_dot_dot_mean_wb(:));
yawaccel_norm = norm*nanstd(yaw_dot_dot_mean_wb(:));
Fenhance_norm = norm*nanstd(F_mean_wb(:));

% plot attitude ditributions inc thresholds
cd('WBmod_figs')

figure
h=hist(F_mean_wb(:),100);
hist(F_mean_wb(:),100);
hold on
plot(1+[Fenhance_limit_steady Fenhance_limit_steady],[0 max(h)],'-k','linewidth',2)
plot(1-[Fenhance_limit_steady Fenhance_limit_steady],[0 max(h)],'-k','linewidth',2)
plot(1+[Fenhance_limit_mod Fenhance_limit_mod],[0 max(h)],'-b','linewidth',2)
plot(1-[Fenhance_limit_mod Fenhance_limit_mod],[0 max(h)],'-b','linewidth',2)
plot(1+[Fenhance_norm Fenhance_norm],[0 max(h)],'-r','linewidth',2)
plot(1-[Fenhance_norm Fenhance_norm],[0 max(h)],'-r','linewidth',2)

saveas(gca,['WBmod_Fenhanced_disr_steadylimit',num2str(limit_steady),'_modlimit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
saveas(gca,['WBmod_Fenhanced_disr_steadylimit',num2str(limit_steady),'_modlimit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
plot2svg(['WBmod_Fenhanced_disr_steadylimit',num2str(limit_steady),'_modlimit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

figure
h=hist(roll_dot_dot_mean_wb(:),100);
hist(roll_dot_dot_mean_wb(:),100);
hold on
plot([rollaccel_limit_steady rollaccel_limit_steady],[0 max(h)],'-k','linewidth',2)
plot(-[rollaccel_limit_steady rollaccel_limit_steady],[0 max(h)],'-k','linewidth',2)
plot([rollaccel_limit_mod rollaccel_limit_mod],[0 max(h)],'-b','linewidth',2)
plot(-[rollaccel_limit_mod rollaccel_limit_mod],[0 max(h)],'-b','linewidth',2)
plot([rollaccel_norm rollaccel_norm],[0 max(h)],'-r','linewidth',2)
plot(-[rollaccel_norm rollaccel_norm],[0 max(h)],'-r','linewidth',2)

saveas(gca,['WBmod_RollAccel_disr_steadylimit',num2str(limit_steady),'_modlimit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
saveas(gca,['WBmod_RollAccel_disr_steadylimit',num2str(limit_steady),'_modlimit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
plot2svg(['WBmod_RollAccel_disr_steadylimit',num2str(limit_steady),'_modlimit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

figure
h=hist(pitch_dot_dot_mean_wb(:),100);
hist(pitch_dot_dot_mean_wb(:),100);
hold on
plot([pitchaccel_limit_steady pitchaccel_limit_steady],[0 max(h)],'-k','linewidth',2)
plot(-[pitchaccel_limit_steady pitchaccel_limit_steady],[0 max(h)],'-k','linewidth',2)
plot([pitchaccel_limit_mod pitchaccel_limit_mod],[0 max(h)],'-b','linewidth',2)
plot(-[pitchaccel_limit_mod pitchaccel_limit_mod],[0 max(h)],'-b','linewidth',2)
plot([pitchaccel_norm pitchaccel_norm],[0 max(h)],'-r','linewidth',2)
plot(-[pitchaccel_norm pitchaccel_norm],[0 max(h)],'-r','linewidth',2)

saveas(gca,['WBmod_PitchAccel_disr_steadylimit',num2str(limit_steady),'_modlimit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
saveas(gca,['WBmod_PitchAccel_disr_steadylimit',num2str(limit_steady),'_modlimit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
plot2svg(['WBmod_PitchAccel_disr_steadylimit',num2str(limit_steady),'_modlimit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

figure
h=hist(yaw_dot_dot_mean_wb(:),100);
hist(yaw_dot_dot_mean_wb(:),100);
hold on
plot([yawaccel_limit_steady yawaccel_limit_steady],[0 max(h)],'-k','linewidth',2)
plot(-[yawaccel_limit_steady yawaccel_limit_steady],[0 max(h)],'-k','linewidth',2)
plot([yawaccel_limit_mod yawaccel_limit_mod],[0 max(h)],'-b','linewidth',2)
plot(-[yawaccel_limit_mod yawaccel_limit_mod],[0 max(h)],'-b','linewidth',2)
plot([yawaccel_norm yawaccel_norm],[0 max(h)],'-r','linewidth',2)
plot(-[yawaccel_norm yawaccel_norm],[0 max(h)],'-r','linewidth',2)

saveas(gca,['WBmod_YawAccel_disr_steadylimit',num2str(limit_steady),'_modlimit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
saveas(gca,['WBmod_YawAccel_disr_steadylimit',num2str(limit_steady),'_modlimit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
plot2svg(['WBmod_YawAccel_disr_steadylimit',num2str(limit_steady),'_modlimit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

cd ..

%% Fenhanced WB
if calc_WB_Fenhanced == 1

    calc_wbNORM_wingatt_Fenhanced

    cd('WBmod_figs')
    
xmin = .5;
xmax = 1.5;

%     plot_wbMOD_wingatt_Fenhanced
%     saveas(gca,['WBmod_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
%     saveas(gca,['WBmod_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
%     plot2svg(['WBmod_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_heatmap_Fenhanced
    saveas(gca,['WBmod_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmodNsteady_Fenhanced
    saveas(gca,['WBmodNsteady_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmodNsteady_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmodNsteady_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_freq_hist_Fenhanced
    saveas(gca,['WBmod_freq_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_freq_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_freq_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_strokecoeff_Fenhanced
    saveas(gca,['WBmod_strokecoeff_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_strokecoeff_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_strokecoeff_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_pitchMIDcoeff_Fenhanced
    saveas(gca,['WBmod_pitchMIDcoeff_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_pitchMIDcoeff_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_pitchMIDcoeff_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devDScoeff_Fenhanced
    saveas(gca,['WBmod_devDScoeff_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devDScoeff_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devDScoeff_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devUScoeff_Fenhanced
    saveas(gca,['WBmod_devUScoeff_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devUScoeff_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devUScoeff_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    cd ..
end


%% Pitch Acceleration WB
if calc_WB_PitchAccel == 1

    calc_wbNORM_wingatt_PitchAccel

    cd('WBmod_figs')
xmin = -1e5;
xmax = -xmin;

%     plot_wbMOD_wingatt_PitchAccel
%     saveas(gca,['WBmod_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
%     saveas(gca,['WBmod_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
%     plot2svg(['WBmod_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_heatmap_PitchAccel
    saveas(gca,['WBmod_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmodNsteady_PitchAccel
    saveas(gca,['WBmodNsteady_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmodNsteady_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmodNsteady_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_freq_hist_PitchAccel
    saveas(gca,['WBmod_freq_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_freq_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_freq_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_strokecoeff_PitchAccel
    saveas(gca,['WBmod_strokecoeff_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_strokecoeff_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_strokecoeff_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_pitchMIDcoeff_PitchAccel
    saveas(gca,['WBmod_pitchMIDcoeff_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_pitchMIDcoeff_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_pitchMIDcoeff_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devDScoeff_PitchAccel
    saveas(gca,['WBmod_devDScoeff_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devDScoeff_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devDScoeff_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devUScoeff_PitchAccel
    saveas(gca,['WBmod_devUScoeff_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devUScoeff_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devUScoeff_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    cd ..
end



%% RollAccel WB
if calc_WB_RollAccel == 1

%     calc_wbNORM_wingatt_RollAccel
    calc_wbNORM_wingatt_RollAccel_updown

    cd('WBmod_figs')
xmin = 0;
xmax = 1e5;
    
%     plot_wbMOD_wingatt_RollAccel
%     saveas(gca,['WBmod_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
%     saveas(gca,['WBmod_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
%     plot2svg(['WBmod_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_heatmap_RollAccel
    saveas(gca,['WBmod_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmodNsteady_RollAccel
    saveas(gca,['WBmodNsteady_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmodNsteady_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmodNsteady_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_strokeMAXcoeff_RollAccel
    saveas(gca,['WBmod_strokeMAXcoeff_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_strokeMAXcoeff_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_strokeMAXcoeff_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_strokeMINcoeff_RollAccel
    saveas(gca,['WBmod_strokeMINcoeff_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_strokeMINcoeff_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_strokeMINcoeff_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_pitchMIDDScoeff_RollAccel
    saveas(gca,['WBmod_pitchMIDDScoeff_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_pitchMIDDScoeff_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_pitchMIDDScoeff_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_pitchMIDUScoeff_RollAccel
    saveas(gca,['WBmod_pitchMIDUScoeff_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_pitchMIDUScoeff_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_pitchMIDUScoeff_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devDSMAXcoeff_RollAccel
    saveas(gca,['WBmod_devDSMAXcoeff_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devDSMAXcoeff_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devDSMAXcoeff_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devDSMINcoeff_RollAccel
    saveas(gca,['WBmod_devDSMINcoeff_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devDSMINcoeff_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devDSMINcoeff_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devUSMAXcoeff_RollAccel
    saveas(gca,['WBmod_devUSMAXcoeff_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devUSMAXcoeff_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devUSMAXcoeff_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_devUSMINcoeff_RollAccel
    saveas(gca,['WBmod_devUSMINcoeff_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devUSMINcoeff_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devUSMINcoeff_RollAccel_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    cd ..
end



%% RollAccel WB
if calc_WB_YawAccel == 1

    calc_wbNORM_wingatt_YawAccel

    cd('WBmod_figs')
xmin = 0;
xmax = 1e5;
    
%     plot_wbMOD_wingatt_YawAccel
%     saveas(gca,['WBmod_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
%     saveas(gca,['WBmod_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
%     plot2svg(['WBmod_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_heatmap_YawAccel
    saveas(gca,['WBmod_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmodNsteady_YawAccel
    saveas(gca,['WBmodNsteady_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmodNsteady_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmodNsteady_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_strokeMAXcoeff_YawAccel
    saveas(gca,['WBmod_strokeMAXcoeff_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_strokeMAXcoeff_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_strokeMAXcoeff_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_strokeMINcoeff_YawAccel
    saveas(gca,['WBmod_strokeMINcoeff_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_strokeMINcoeff_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_strokeMINcoeff_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_pitchMIDDScoeff_YawAccel
    saveas(gca,['WBmod_pitchMIDDScoeff_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_pitchMIDDScoeff_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_pitchMIDDScoeff_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_pitchMIDUScoeff_YawAccel
    saveas(gca,['WBmod_pitchMIDUScoeff_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_pitchMIDUScoeff_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_pitchMIDUScoeff_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devDSMAXcoeff_YawAccel
    saveas(gca,['WBmod_devDSMAXcoeff_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devDSMAXcoeff_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devDSMAXcoeff_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devDSMINcoeff_YawAccel
    saveas(gca,['WBmod_devDSMINcoeff_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devDSMINcoeff_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devDSMINcoeff_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devUSMAXcoeff_YawAccel
    saveas(gca,['WBmod_devUSMAXcoeff_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devUSMAXcoeff_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devUSMAXcoeff_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_devUSMINcoeff_YawAccel
    saveas(gca,['WBmod_devUSMINcoeff_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devUSMINcoeff_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devUSMINcoeff_YawAccel_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    cd ..
end




