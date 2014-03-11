clc
clear
close all

addpath('/home/florian/Dropbox/WORK/toolbox')
addpath('/home/florian/Dropbox/WORK/toolbox/plot2svg')
addpath('/home/florian/Dropbox/WORK/toolbox/CircStat')
addpath('/home/florian/Dropbox/WORK/toolbox/flytracker')

loadname=dir('WBdataset_all_*')
loadname=loadname.name;

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

mkdir('WBmod_figs_damping')
load(loadname)
if calc_WB_STEADY == 0
    steady_name = dir('WBdataset_torquebased_steady_*')
if isempty(steady_name) == 0
    steady_name = steady_name.name;
    load(steady_name)
end
end

%% settings
% limits (*std)
limit_steady = .5
limit_mod = 1
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
    
%     calc_wbNORM_wingatt_steady
    calc_torquebased_wbNORM_wingatt_steady
    
    cd('WBmod_figs_damping')
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

% % limits
% rollaccel_limit_mod = limit_mod*nanstd(roll_dot_dot_mean_wb(:));
% pitchaccel_limit_mod = limit_mod*nanstd(pitch_dot_dot_mean_wb(:));
% yawaccel_limit_mod = limit_mod*nanstd(yaw_dot_dot_mean_wb(:));
% Fenhance_limit_mod = limit_mod*nanstd(F_mean_wb(:));
% 
% % normalization value
% rollaccel_norm = norm*nanstd(roll_dot_dot_mean_wb(:));
% pitchaccel_norm = norm*nanstd(pitch_dot_dot_mean_wb(:));
% yawaccel_norm = norm*nanstd(yaw_dot_dot_mean_wb(:));
% Fenhance_norm = norm*nanstd(F_mean_wb(:));

% TORQUE based constants
% limits
Mroll_limit_mod = limit_mod*nanstd(Mroll_mean_wb(:));
Mpitch_limit_mod = limit_mod*nanstd(Mpitch_mean_wb(:));
Myaw_limit_mod = limit_mod*nanstd(Myaw_mean_wb(:));
Fenhance_limit_mod = limit_mod*nanstd(F_mean_wb(:));

% normalization value
Mroll_norm = norm*nanstd(Mroll_mean_wb(:));
Mpitch_norm = norm*nanstd(Mpitch_mean_wb(:));
Myaw_norm = norm*nanstd(Myaw_mean_wb(:));
Fenhance_norm = norm*nanstd(F_mean_wb(:));

%% Fenhanced WB
if calc_WB_Fenhanced == 1

    calc_torquebased_wbNORM_wingatt_Fenhanced

    cd('WBmod_figs_damping')
    
xmin = .5;
xmax = 2.5;

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

    calc_torquebased_wbNORM_wingatt_PitchAccel

    cd('WBmod_figs_damping')
xmin = -.5;
xmax = -xmin;

%     plot_wbMOD_wingatt_PitchAccel
%     saveas(gca,['WBmod_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
%     saveas(gca,['WBmod_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
%     plot2svg(['WBmod_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_heatmap_Mpitch
    saveas(gca,['WBmod_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmodNsteady_Mpitch
    saveas(gca,['WBmodNsteady_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmodNsteady_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmodNsteady_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_freq_hist_Mpitch
    saveas(gca,['WBmod_freq_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_freq_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_freq_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_strokecoeff_Mpitch
    saveas(gca,['WBmod_strokecoeff_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_strokecoeff_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_strokecoeff_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_pitchMIDcoeff_Mpitch
    saveas(gca,['WBmod_pitchMIDcoeff_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_pitchMIDcoeff_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_pitchMIDcoeff_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devDScoeff_Mpitch
    saveas(gca,['WBmod_devDScoeff_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devDScoeff_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devDScoeff_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devUScoeff_Mpitch
    saveas(gca,['WBmod_devUScoeff_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devUScoeff_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devUScoeff_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    cd ..
end



%% RollAccel WB
if calc_WB_RollAccel == 1

%     calc_wbNORM_wingatt_RollAccel
%     calc_wbNORM_wingatt_RollAccel_updown
%     calc_wbNORM_wingatt_RollAccel_damping
    calc_torquebased_wbNORM_wingatt_RollAccel

    cd('WBmod_figs_damping')
xmin = 0;
xmax = 1;
    
    plot_WBmod_heatmap_Mroll
    saveas(gca,['WBmod_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmodNsteady_Mroll
    saveas(gca,['WBmodNsteady_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmodNsteady_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmodNsteady_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_strokeMAXcoeff_Mroll
    saveas(gca,['WBmod_strokeMAXcoeff_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_strokeMAXcoeff_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_strokeMAXcoeff_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_strokeMINcoeff_Mroll
    saveas(gca,['WBmod_strokeMINcoeff_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_strokeMINcoeff_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_strokeMINcoeff_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_pitchMIDDScoeff_Mroll
    saveas(gca,['WBmod_pitchMIDDScoeff_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_pitchMIDDScoeff_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_pitchMIDDScoeff_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_pitchMIDUScoeff_Mroll
    saveas(gca,['WBmod_pitchMIDUScoeff_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_pitchMIDUScoeff_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_pitchMIDUScoeff_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devDSMAXcoeff_Mroll
    saveas(gca,['WBmod_devDSMAXcoeff_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devDSMAXcoeff_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devDSMAXcoeff_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devDSMINcoeff_Mroll
    saveas(gca,['WBmod_devDSMINcoeff_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devDSMINcoeff_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devDSMINcoeff_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devUSMAXcoeff_Mroll
    saveas(gca,['WBmod_devUSMAXcoeff_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devUSMAXcoeff_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devUSMAXcoeff_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_devUSMINcoeff_Mroll
    saveas(gca,['WBmod_devUSMINcoeff_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devUSMINcoeff_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devUSMINcoeff_Mroll_',num2str(n_RollAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    cd ..
end



%% RollAccel WB
if calc_WB_YawAccel == 1

    calc_torquebased_wbNORM_wingatt_YawAccel

    cd('WBmod_figs_damping')
xmin = 0;
xmax = .25;
    
    plot_WBmod_heatmap_YawAccel
    saveas(gca,['WBmod_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmodNsteady_YawAccel
    saveas(gca,['WBmodNsteady_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmodNsteady_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmodNsteady_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_strokeMAXcoeff_Myaw
    saveas(gca,['WBmod_strokeMAXcoeff_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_strokeMAXcoeff_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_strokeMAXcoeff_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_strokeMINcoeff_Myaw
    saveas(gca,['WBmod_strokeMINcoeff_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_strokeMINcoeff_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_strokeMINcoeff_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_pitchMIDDScoeff_Myaw
    saveas(gca,['WBmod_pitchMIDDScoeff_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_pitchMIDDScoeff_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_pitchMIDDScoeff_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_pitchMIDUScoeff_Myaw
    saveas(gca,['WBmod_pitchMIDUScoeff_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_pitchMIDUScoeff_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_pitchMIDUScoeff_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devDSMAXcoeff_Myaw
    saveas(gca,['WBmod_devDSMAXcoeff_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devDSMAXcoeff_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devDSMAXcoeff_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devDSMINcoeff_Myaw
    saveas(gca,['WBmod_devDSMINcoeff_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devDSMINcoeff_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devDSMINcoeff_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devUSMAXcoeff_Myaw
    saveas(gca,['WBmod_devUSMAXcoeff_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devUSMAXcoeff_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devUSMAXcoeff_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_devUSMINcoeff_Myaw
    saveas(gca,['WBmod_devUSMINcoeff_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devUSMINcoeff_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devUSMINcoeff_Myaw_',num2str(n_YawAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    cd ..
end




