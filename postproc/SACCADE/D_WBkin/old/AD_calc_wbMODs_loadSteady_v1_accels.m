clc
clear
close all

addpath('/home/florian/Dropbox/WORK/toolbox')
addpath('/home/florian/Dropbox/WORK/toolbox/plot2svg')
addpath('/home/florian/Dropbox/WORK/toolbox/CircStat')
addpath('/home/florian/Dropbox/WORK/toolbox/flytracker')

loadname=dir('WBdataset_all_*')
loadname=loadname.name;
load(loadname)

steady_name=dir('WBdataset_steady_*')
steady_name=steady_name.name;
load(steady_name)

plot_on = 1;
% plot_on = 0;

save_on = 1;
% save_on = 0;

calc_WB_Fenhanced = 0;
calc_WB_PitchAccel_UP = 0;
calc_WB_PitchAccel_DOWN = 0;
calc_WB_PitchAccel = 0;
calc_WB_RollAccel = 0;
calc_WB_YawAccel = 0;

calc_WB_Fenhanced = 1
% calc_WB_PitchAccel_UP = 1;
% calc_WB_PitchAccel_DOWN = 1;
calc_WB_PitchAccel = 1
calc_WB_RollAccel = 1
calc_WB_YawAccel = 1

mkdir('WBmod_figs_accel')

%% settings
% limits (*std)
limit_mod = 1
norm = 3

% fourier orders
MOD_fourier_order = 8;

% number of polinomials
n_pol_MOD = 10;

% plot_fits = 0;
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

%% calc & plot MOD wb
% steady wb data
vel_wb_steady = V_steady_meanCIstd(:,1);
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
rollaccel_limit_mod = limit_mod*nanstd(roll_dot_dot_mean_wb(:));
pitchaccel_limit_mod = limit_mod*nanstd(pitch_dot_dot_mean_wb(:));
yawaccel_limit_mod = limit_mod*nanstd(yaw_dot_dot_mean_wb(:));
Fenhance_limit_mod = limit_mod*nanstd(F_mean_wb(:));

% normalization value
rollaccel_norm = norm*nanstd(roll_dot_dot_mean_wb(:));
pitchaccel_norm = norm*nanstd(pitch_dot_dot_mean_wb(:));
yawaccel_norm = norm*nanstd(yaw_dot_dot_mean_wb(:));
Fenhance_norm = norm*nanstd(F_mean_wb(:));

% % TORQUE based constants
% % limits
% Mroll_limit_mod = limit_mod*nanstd(Mroll_mean_wb(:));
% Mpitch_limit_mod = limit_mod*nanstd(Mpitch_mean_wb(:));
% Myaw_limit_mod = limit_mod*nanstd(Myaw_mean_wb(:));
% Fenhance_limit_mod = limit_mod*nanstd(F_mean_wb(:));
% 
% % normalization value
% Mroll_norm = norm*nanstd(Mroll_mean_wb(:));
% Mpitch_norm = norm*nanstd(Mpitch_mean_wb(:));
% Myaw_norm = norm*nanstd(Myaw_mean_wb(:));
% Fenhance_norm = norm*nanstd(F_mean_wb(:));

%% save norm_data
save('norm_data.mat','rollaccel_norm','pitchaccel_norm','yawaccel_norm','Fenhance_norm','f_wb_steady')

%% Fenhanced WB
if calc_WB_Fenhanced == 1

%     calc_torquebased_wbNORM_wingatt_Fenhanced
    calc_accelbased_wbNORM_wingatt_Fenhanced
    
    
    if plot_on ==1
    cd('WBmod_figs_accel')
    
xmin = .75;
xmax = 1.5;

%     plot_wbMOD_wingatt_Fenhanced
%     saveas(gca,['WBmod_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
%     saveas(gca,['WBmod_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
%     plot2svg(['WBmod_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

figure
h=hist(F_mean_wb,50);
hist(F_mean_wb,50)
hold on
plot([1+Fenhance_limit_mod 1+Fenhance_limit_mod],[0 max(h)],'-r')
title('F/g')
saveas(gca,['hist_Fnorm',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
saveas(gca,['hist_Fnorm',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
plot2svg(['hist_Fnorm',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_heatmap_Fenhanced
    saveas(gca,['WBmod_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmodNsteady_Fenhanced
    saveas(gca,['WBmodNsteady_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmodNsteady_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmodNsteady_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_vel_hist_Fenhanced
    saveas(gca,['WBmod_vel_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_vel_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_vel_Fenhanced_',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
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
end


%% Pitch Acceleration WB UP
if calc_WB_PitchAccel_UP == 1

%     calc_torquebased_wbNORM_wingatt_PitchAccel
%     calc_accelbased_wbNORM_wingatt_PitchAccel
    calc_accelbased_wbNORM_wingatt_PitchAccelUP

    if plot_on ==1

    cd('WBmod_figs_accel')
xmin = -ceil(max(abs(pitch_dot_dot_mean_wb)));
xmin = -1.5;
xmax = -xmin;


figure
h=hist(pitch_dot_dot_mean_wb,50);
hist(pitch_dot_dot_mean_wb,50)
hold on
plot([pitchaccel_limit_mod pitchaccel_limit_mod],[0 max(h)],'-r')
plot([-pitchaccel_limit_mod -pitchaccel_limit_mod],[0 max(h)],'-r')
title('pitch accel')
saveas(gca,['hist_pitchAccel',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
saveas(gca,['hist_pitchAccel',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
plot2svg(['hist_pitchAccel',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

%     plot_wbMOD_wingatt_PitchAccel
%     saveas(gca,['WBmod_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
%     saveas(gca,['WBmod_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
%     plot2svg(['WBmod_PitchAccel_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_heatmap_PitchAccel
    saveas(gca,['WBmod_PitchAccel_UP_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_PitchAccel_UP_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_PitchAccel_UP_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmodNsteady_PitchAccelUP
    saveas(gca,['WBmodNsteady_PitchAccel_UP_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmodNsteady_PitchAccel_UP_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmodNsteady_PitchAccel_UP_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_freq_hist_PitchAccel
    saveas(gca,['WBmod_freq_PitchAccel_UP_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_freq_PitchAccel_UP_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_freq_PitchAccel_UP_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_strokecoeff_PitchAccel
    saveas(gca,['WBmod_strokecoeff_PitchAccel_UP_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_strokecoeff_PitchAccel_UP_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_strokecoeff_PitchAccel_UP_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_pitchMIDcoeff_PitchAccel
    saveas(gca,['WBmod_pitchMIDcoeff_PitchAccel_UP_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_pitchMIDcoeff_PitchAccel_UP_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_pitchMIDcoeff_PitchAccel_UP_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devDScoeff_PitchAccel
    saveas(gca,['WBmod_devDScoeff_PitchAccel_UP_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devDScoeff_PitchAccel_UP_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devDScoeff_PitchAccel_UP_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devUScoeff_PitchAccel
    saveas(gca,['WBmod_devUScoeff_PitchAccel_UP_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devUScoeff_PitchAccel_UP_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devUScoeff_PitchAccel_UP_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    cd ..
    end
end

%% Pitch Acceleration WB DOWN
if calc_WB_PitchAccel_DOWN == 1

%     calc_torquebased_wbNORM_wingatt_PitchAccel
%     calc_accelbased_wbNORM_wingatt_PitchAccel
    calc_accelbased_wbNORM_wingatt_PitchAccelDOWN

    if plot_on ==1

    cd('WBmod_figs_accel')
    
xmin = -ceil(max(abs(pitch_dot_dot_mean_wb)));
xmin = -1.5;
xmax = -xmin;


figure
h=hist(pitch_dot_dot_mean_wb,50);
hist(pitch_dot_dot_mean_wb,50)
hold on
plot([pitchaccel_limit_mod pitchaccel_limit_mod],[0 max(h)],'-r')
plot([-pitchaccel_limit_mod -pitchaccel_limit_mod],[0 max(h)],'-r')
title('pitch accel')
saveas(gca,['hist_pitchAccel',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
saveas(gca,['hist_pitchAccel',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
plot2svg(['hist_pitchAccel',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

%     plot_wbMOD_wingatt_PitchAccel
%     saveas(gca,['WBmod_PitchAccel_DOWN_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
%     saveas(gca,['WBmod_PitchAccel_DOWN_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
%     plot2svg(['WBmod_PitchAccel_DOWN_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_heatmap_PitchAccel
    saveas(gca,['WBmod_PitchAccel_DOWN_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_PitchAccel_DOWN_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_PitchAccel_DOWN_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmodNsteady_PitchAccelDOWN
    saveas(gca,['WBmodNsteady_PitchAccel_DOWN_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmodNsteady_PitchAccel_DOWN_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmodNsteady_PitchAccel_DOWN_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_freq_hist_PitchAccel
    saveas(gca,['WBmod_freq_PitchAccel_DOWN_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_freq_PitchAccel_DOWN_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_freq_PitchAccel_DOWN_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_strokecoeff_PitchAccel
    saveas(gca,['WBmod_strokecoeff_PitchAccel_DOWN_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_strokecoeff_PitchAccel_DOWN_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_strokecoeff_PitchAccel_DOWN_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_pitchMIDcoeff_PitchAccel
    saveas(gca,['WBmod_pitchMIDcoeff_PitchAccel_DOWN_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_pitchMIDcoeff_PitchAccel_DOWN_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_pitchMIDcoeff_PitchAccel_DOWN_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devDScoeff_PitchAccel
    saveas(gca,['WBmod_devDScoeff_PitchAccel_DOWN_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devDScoeff_PitchAccel_DOWN_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devDScoeff_PitchAccel_DOWN_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devUScoeff_PitchAccel
    saveas(gca,['WBmod_devUScoeff_PitchAccel_DOWN_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devUScoeff_PitchAccel_DOWN_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devUScoeff_PitchAccel_DOWN_',num2str(n_PitchAccel),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    cd ..
    end
end

%% Pitch Acceleration WB
if calc_WB_PitchAccel == 1

%     calc_torquebased_wbNORM_wingatt_PitchAccel
    calc_accelbased_wbNORM_wingatt_PitchAccel

    if plot_on ==1

    cd('WBmod_figs_accel')
    
xmin = -ceil(max(abs(pitch_dot_dot_mean_wb)));
xmin = -1.5;
xmax = -xmin;


figure
h=hist(pitch_dot_dot_mean_wb,50);
hist(pitch_dot_dot_mean_wb,50)
hold on
plot([pitchaccel_limit_mod pitchaccel_limit_mod],[0 max(h)],'-r')
plot([-pitchaccel_limit_mod -pitchaccel_limit_mod],[0 max(h)],'-r')
title('pitch accel')
saveas(gca,['hist_pitchAccel',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
saveas(gca,['hist_pitchAccel',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
plot2svg(['hist_pitchAccel',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

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
end

%% RollAccel WB
if calc_WB_RollAccel == 1

%     calc_wbNORM_wingatt_RollAccel
%     calc_wbNORM_wingatt_RollAccel_updown
%     calc_wbNORM_wingatt_RollAccel_damping
%     calc_torquebased_wbNORM_wingatt_RollAccel
    calc_accelbased_wbNORM_wingatt_RollAccel

    if plot_on ==1
    cd('WBmod_figs_accel')

xmin = 0;
xmax = ceil(max(abs(roll_dot_dot_mean_wb)));

figure
h=hist(roll_dot_dot_mean_wb,50);
hist(roll_dot_dot_mean_wb,50)
hold on
plot([rollaccel_limit_mod rollaccel_limit_mod],[0 max(h)],'-r')
plot([-rollaccel_limit_mod -rollaccel_limit_mod],[0 max(h)],'-r')
title('roll accel')
saveas(gca,['hist_rollAccel',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
saveas(gca,['hist_rollAccel',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
plot2svg(['hist_rollAccel',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    
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
end



%% YawAccel WB
if calc_WB_YawAccel == 1

%     calc_torquebased_wbNORM_wingatt_YawAccel
    calc_accelbased_wbNORM_wingatt_YawAccel

    if plot_on ==1
    cd('WBmod_figs_accel')
xmin = 0;
xmax = ceil(max(abs(yaw_dot_dot_mean_wb)));

figure
h=hist(yaw_dot_dot_mean_wb,50);
hist(yaw_dot_dot_mean_wb,50)
hold on
plot([yawaccel_limit_mod yawaccel_limit_mod],[0 max(h)],'-r')
plot([-yawaccel_limit_mod -yawaccel_limit_mod],[0 max(h)],'-r')
title('yaw accel')
saveas(gca,['hist_yawAccel',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
saveas(gca,['hist_yawAccel',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
plot2svg(['hist_yawAccel',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    
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
end


%% save all
    save('WBdataset_all_steadyNmods_accelNorm.mat')

%% MS plots
mkdir('MSfigs_WBkin')
cd('MSfigs_WBkin')

MSplot_WBfunc_heatmap_steady_NOtick
    saveas(gca,['MSplot_WBfunc_heatmap_steady_NOtick.fig'])
    saveas(gca,['MSplot_WBfunc_heatmap_steady_NOtick.png'])
    plot2svg(['MSplot_WBfunc_heatmap_steady_NOtick.svg'])

MSplot_WBmodsNsteady_onefig_2xNorm
    saveas(gca,['MSplot_WBmodsNsteady_onefig.fig'])
    saveas(gca,['MSplot_WBmodsNsteady_onefig.png'])
    plot2svg(['MSplot_WBmodsNsteady_onefig.svg'])

MSplot_WBsteady2mod2xNorm_Fenhance
    saveas(gca,['MSplot_WBsteady2mod2xNorm_Fenhance.fig'])
    saveas(gca,['MSplot_WBsteady2mod2xNorm_Fenhance.png'])
    plot2svg(['MSplot_WBsteady2mod2xNorm_Fenhance.svg'])

MSplot_WBsteady2mod2xNorm_RollAccel
    saveas(gca,['MSplot_WBsteady2mod2xNorm_RollAccel.fig'])
    saveas(gca,['MSplot_WBsteady2mod2xNorm_RollAccel.png'])
    plot2svg(['MSplot_WBsteady2mod2xNorm_RollAccel.svg'])

MSplot_WBsteady2mod2xNorm_PitchAccel
    saveas(gca,['MSplot_WBsteady2mod2xNorm_PitchAccel.fig'])
    saveas(gca,['MSplot_WBsteady2mod2xNorm_PitchAccel.png'])
    plot2svg(['MSplot_WBsteady2mod2xNorm_PitchAccel.svg'])

MSplot_WBsteady2mod2xNorm
    saveas(gca,['MSplot_WBsteady2mod2xNorm.fig'])
    saveas(gca,['MSplot_WBsteady2mod2xNorm.png'])
    plot2svg(['MSplot_WBsteady2mod2xNorm.svg'])

MSplot_WBsteady2mod2xNorm_1fig
    saveas(gca,['MSplot_WBsteady2mod2xNorm_1fig.fig'])
    saveas(gca,['MSplot_WBsteady2mod2xNorm_1fig.png'])
    plot2svg(['MSplot_WBsteady2mod2xNorm_1fig.svg'])

cd ..


