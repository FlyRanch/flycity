clc
clear
close all

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
calc_WB_PitchTorque_UP = 0;
calc_WB_PitchTorque_DOWN = 0;
calc_WB_PitchTorque = 0;
calc_WB_RollTorque = 0;
calc_WB_YawTorque = 0;
calc_WB_TorqueAxisR = 0;

% calc_WB_Fenhanced = 1
% % calc_WB_PitchTorque_UP = 1;
% % calc_WB_PitchTorque_DOWN = 1;
% calc_WB_PitchTorque = 1
% calc_WB_RollTorque = 1
% calc_WB_YawTorque = 1
calc_WB_TorqueAxisR = 1

mkdir('WBmod_figs_torque')

%% settings
% limits (*std)
limit_mod = .5
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
Mroll_limit_mod = limit_mod*nanstd(Mroll_mean_wb(:));
Mpitch_limit_mod = limit_mod*nanstd(Mpitch_mean_wb(:));
Myaw_limit_mod = limit_mod*nanstd(Myaw_mean_wb(:));
Fenhance_limit_mod = limit_mod*nanstd(F_mean_wb(:));
M_R_limit_mod = limit_mod*nanstd(M_R_mean_wb(:));
M_L_limit_mod = limit_mod*nanstd(M_L_mean_wb(:));

% normalization value
Mroll_norm = norm*nanstd(Mroll_mean_wb(:));
Mpitch_norm = norm*nanstd(Mpitch_mean_wb(:));
Myaw_norm = norm*nanstd(Myaw_mean_wb(:));
Fenhance_norm = norm*nanstd(F_mean_wb(:));
M_R_norm = norm*nanstd(M_R_mean_wb(:));
M_L_norm = norm*nanstd(M_L_mean_wb(:));

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
save('norm_data_torque.mat','Mroll_norm','Mpitch_norm','Myaw_norm','Fenhance_norm','M_R_norm','M_L_norm','f_wb_steady')

%% Fenhanced WB
if calc_WB_Fenhanced == 1

%     calc_accelbased_wbNORM_wingatt_Fenhanced
    calc_torquebased_wbNORM_wingatt_Fenhanced
    
    
    if plot_on ==1
    cd('WBmod_figs_torque')
    
xmin = .75;
xmax = 1.5;

figure
h=hist(F_mean_wb,50);
hist(F_mean_wb,50)
hold on
plot([1+Fenhance_limit_mod 1+Fenhance_limit_mod],[0 max(h)],'-r')
title('F/g')
saveas(gca,['hist_Fnorm',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
saveas(gca,['hist_Fnorm',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
plot2svg(['hist_Fnorm',num2str(n_Fenhance),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

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


%% Pitch Torque WB UP
if calc_WB_PitchTorque_UP == 1

%     calc_torquebased_wbNORM_wingatt_PitchAccel
%     calc_accelbased_wbNORM_wingatt_PitchAccel
%     calc_accelbased_wbNORM_wingatt_PitchAccelUP
    calc_Torquebased_wbNORM_wingatt_PitchTorqueUP

    if plot_on ==1

    cd('WBmod_figs_torque')
xmin = -ceil(max(abs(Mpitch_mean_wb)));
xmin = -1.5;
xmax = -xmin;


figure
h=hist(Mpitch_mean_wb,50);
hist(Mpitch_mean_wb,50)
hold on
plot([Mpitch_limit_mod Mpitch_limit_mod],[0 max(h)],'-r')
plot([-Mpitch_limit_mod -Mpitch_limit_mod],[0 max(h)],'-r')
title('pitch Torque')
saveas(gca,['hist_pitchTorque',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
saveas(gca,['hist_pitchTorque',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
plot2svg(['hist_pitchTorque',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

%     plot_wbMOD_wingatt_PitchTorque
%     saveas(gca,['WBmod_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
%     saveas(gca,['WBmod_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
%     plot2svg(['WBmod_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_heatmap_PitchTorque
    saveas(gca,['WBmod_PitchTorque_UP_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_PitchTorque_UP_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_PitchTorque_UP_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmodNsteady_PitchTorqueUP
    saveas(gca,['WBmodNsteady_PitchTorque_UP_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmodNsteady_PitchTorque_UP_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmodNsteady_PitchTorque_UP_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_freq_hist_PitchTorque
    saveas(gca,['WBmod_freq_PitchTorque_UP_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_freq_PitchTorque_UP_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_freq_PitchTorque_UP_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_strokecoeff_PitchTorque
    saveas(gca,['WBmod_strokecoeff_PitchTorque_UP_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_strokecoeff_PitchTorque_UP_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_strokecoeff_PitchTorque_UP_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_pitchMIDcoeff_PitchTorque
    saveas(gca,['WBmod_pitchMIDcoeff_PitchTorque_UP_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_pitchMIDcoeff_PitchTorque_UP_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_pitchMIDcoeff_PitchTorque_UP_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devDScoeff_PitchTorque
    saveas(gca,['WBmod_devDScoeff_PitchTorque_UP_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devDScoeff_PitchTorque_UP_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devDScoeff_PitchTorque_UP_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devUScoeff_PitchTorque
    saveas(gca,['WBmod_devUScoeff_PitchTorque_UP_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devUScoeff_PitchTorque_UP_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devUScoeff_PitchTorque_UP_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    cd ..
    end
end

%% Pitch Torque WB DOWN
if calc_WB_PitchTorque_DOWN == 1

%     calc_torquebased_wbNORM_wingatt_PitchAccel
%     calc_accelbased_wbNORM_wingatt_PitchAccel
%     calc_accelbased_wbNORM_wingatt_PitchAccelDOWN
    calc_Torquebased_wbNORM_wingatt_PitchTorqueDOWN

    if plot_on ==1

    cd('WBmod_figs_torque')
    
xmin = -ceil(max(abs(Mpitch_mean_wb)));
xmin = -1.5;
xmax = -xmin;


figure
h=hist(Mpitch_mean_wb,50);
hist(Mpitch_mean_wb,50)
hold on
plot([Mpitch_limit_mod Mpitch_limit_mod],[0 max(h)],'-r')
plot([-Mpitch_limit_mod -Mpitch_limit_mod],[0 max(h)],'-r')
title('pitch Torque')
saveas(gca,['hist_pitchTorque',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
saveas(gca,['hist_pitchTorque',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
plot2svg(['hist_pitchTorque',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

%     plot_wbMOD_wingatt_PitchTorque
%     saveas(gca,['WBmod_PitchTorque_DOWN_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
%     saveas(gca,['WBmod_PitchTorque_DOWN_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
%     plot2svg(['WBmod_PitchTorque_DOWN_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_heatmap_PitchTorque
    saveas(gca,['WBmod_PitchTorque_DOWN_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_PitchTorque_DOWN_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_PitchTorque_DOWN_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmodNsteady_PitchTorqueDOWN
    saveas(gca,['WBmodNsteady_PitchTorque_DOWN_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmodNsteady_PitchTorque_DOWN_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmodNsteady_PitchTorque_DOWN_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_freq_hist_PitchTorque
    saveas(gca,['WBmod_freq_PitchTorque_DOWN_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_freq_PitchTorque_DOWN_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_freq_PitchTorque_DOWN_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_strokecoeff_PitchTorque
    saveas(gca,['WBmod_strokecoeff_PitchTorque_DOWN_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_strokecoeff_PitchTorque_DOWN_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_strokecoeff_PitchTorque_DOWN_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_pitchMIDcoeff_PitchTorque
    saveas(gca,['WBmod_pitchMIDcoeff_PitchTorque_DOWN_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_pitchMIDcoeff_PitchTorque_DOWN_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_pitchMIDcoeff_PitchTorque_DOWN_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devDScoeff_PitchTorque
    saveas(gca,['WBmod_devDScoeff_PitchTorque_DOWN_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devDScoeff_PitchTorque_DOWN_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devDScoeff_PitchTorque_DOWN_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devUScoeff_PitchTorque
    saveas(gca,['WBmod_devUScoeff_PitchTorque_DOWN_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devUScoeff_PitchTorque_DOWN_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devUScoeff_PitchTorque_DOWN_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    cd ..
    end
end

%% Pitch Torque WB
if calc_WB_PitchTorque == 1

%     calc_torquebased_wbNORM_wingatt_PitchAccel
%     calc_accelbased_wbNORM_wingatt_PitchAccel
    calc_Torquebased_wbNORM_wingatt_PitchTorque

    if plot_on ==1

    cd('WBmod_figs_torque')
    
xmin = -ceil(100*max(abs(Mpitch_mean_wb)))/100;
% xmin = -ceil(max(abs(Mpitch_mean_wb)));
% xmin = -.05;
xmax = -xmin;


figure
h=hist(Mpitch_mean_wb,50);
hist(Mpitch_mean_wb,50)
hold on
plot([Mpitch_limit_mod Mpitch_limit_mod],[0 max(h)],'-r')
plot([-Mpitch_limit_mod -Mpitch_limit_mod],[0 max(h)],'-r')
title('pitch Torque')
saveas(gca,['hist_pitchTorque',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
saveas(gca,['hist_pitchTorque',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
plot2svg(['hist_pitchTorque',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

%     plot_wbMOD_wingatt_PitchTorque
%     saveas(gca,['WBmod_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
%     saveas(gca,['WBmod_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
%     plot2svg(['WBmod_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_heatmap_PitchTorque
    saveas(gca,['WBmod_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmodNsteady_PitchTorque
    saveas(gca,['WBmodNsteady_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmodNsteady_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmodNsteady_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_freq_hist_PitchTorque
    saveas(gca,['WBmod_freq_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_freq_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_freq_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_strokecoeff_PitchTorque
    saveas(gca,['WBmod_strokecoeff_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_strokecoeff_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_strokecoeff_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_pitchMIDcoeff_PitchTorque
    saveas(gca,['WBmod_pitchMIDcoeff_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_pitchMIDcoeff_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_pitchMIDcoeff_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devDScoeff_PitchTorque
    saveas(gca,['WBmod_devDScoeff_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devDScoeff_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devDScoeff_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devUScoeff_PitchTorque
    saveas(gca,['WBmod_devUScoeff_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devUScoeff_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devUScoeff_PitchTorque_',num2str(n_PitchTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    cd ..
    end
end

%% RollTorque WB
if calc_WB_RollTorque == 1

%     calc_wbNORM_wingatt_RollAccel
%     calc_wbNORM_wingatt_RollAccel_updown
%     calc_wbNORM_wingatt_RollAccel_damping
%     calc_torquebased_wbNORM_wingatt_RollAccel
%     calc_accelbased_wbNORM_wingatt_RollAccel
    calc_Torquebased_wbNORM_wingatt_RollTorque

    if plot_on ==1
    cd('WBmod_figs_torque')

xmin = 0;
xmax = ceil(100*max(abs(Mroll_mean_wb)))/100;

figure
h=hist(Mroll_mean_wb,50);
hist(Mroll_mean_wb,50)
hold on
plot([Mroll_limit_mod Mroll_limit_mod],[0 max(h)],'-r')
plot([-Mroll_limit_mod -Mroll_limit_mod],[0 max(h)],'-r')
title('roll Torque')
saveas(gca,['hist_rollTorque',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
saveas(gca,['hist_rollTorque',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
plot2svg(['hist_rollTorque',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_heatmap_RollTorque
    saveas(gca,['WBmod_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmodNsteady_RollTorque
    saveas(gca,['WBmodNsteady_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmodNsteady_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmodNsteady_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_strokeMAXcoeff_RollTorque
    saveas(gca,['WBmod_strokeMAXcoeff_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_strokeMAXcoeff_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_strokeMAXcoeff_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_strokeMINcoeff_RollTorque
    saveas(gca,['WBmod_strokeMINcoeff_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_strokeMINcoeff_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_strokeMINcoeff_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_pitchMIDDScoeff_RollTorque
    saveas(gca,['WBmod_pitchMIDDScoeff_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_pitchMIDDScoeff_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_pitchMIDDScoeff_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_pitchMIDUScoeff_RollTorque
    saveas(gca,['WBmod_pitchMIDUScoeff_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_pitchMIDUScoeff_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_pitchMIDUScoeff_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devDSMAXcoeff_RollTorque
    saveas(gca,['WBmod_devDSMAXcoeff_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devDSMAXcoeff_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devDSMAXcoeff_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devDSMINcoeff_RollTorque
    saveas(gca,['WBmod_devDSMINcoeff_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devDSMINcoeff_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devDSMINcoeff_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devUSMAXcoeff_RollTorque
    saveas(gca,['WBmod_devUSMAXcoeff_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devUSMAXcoeff_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devUSMAXcoeff_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_devUSMINcoeff_RollTorque
    saveas(gca,['WBmod_devUSMINcoeff_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devUSMINcoeff_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devUSMINcoeff_RollTorque_',num2str(n_RollTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    cd ..
    end
end



%% YawTorque WB
if calc_WB_YawTorque == 1

%     calc_torquebased_wbNORM_wingatt_YawAccel
%     calc_accelbased_wbNORM_wingatt_YawAccel
    calc_Torquebased_wbNORM_wingatt_YawTorque_onlyPos

    if plot_on ==1
    cd('WBmod_figs_torque')
xmin = 0;
xmax = ceil(100*max(abs(Myaw_mean_wb)))/100;

figure
h=hist(Myaw_mean_wb,50);
hist(Myaw_mean_wb,50)
hold on
plot([Myaw_limit_mod Myaw_limit_mod],[0 max(h)],'-r')
% plot([-Myaw_limit_mod -Myaw_limit_mod],[0 max(h)],'-r')
title('yaw Torque')
saveas(gca,['hist_yawTorque',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
saveas(gca,['hist_yawTorque',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
plot2svg(['hist_yawTorque',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_heatmap_YawTorque
    saveas(gca,['WBmod_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmodNsteady_YawTorque
    saveas(gca,['WBmodNsteady_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmodNsteady_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmodNsteady_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_strokeMAXcoeff_YawTorque
    saveas(gca,['WBmod_strokeMAXcoeff_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_strokeMAXcoeff_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_strokeMAXcoeff_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_strokeMINcoeff_YawTorque
    saveas(gca,['WBmod_strokeMINcoeff_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_strokeMINcoeff_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_strokeMINcoeff_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_pitchMIDDScoeff_YawTorque
    saveas(gca,['WBmod_pitchMIDDScoeff_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_pitchMIDDScoeff_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_pitchMIDDScoeff_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_pitchMIDUScoeff_YawTorque
    saveas(gca,['WBmod_pitchMIDUScoeff_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_pitchMIDUScoeff_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_pitchMIDUScoeff_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devDSMAXcoeff_YawTorque
    saveas(gca,['WBmod_devDSMAXcoeff_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devDSMAXcoeff_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devDSMAXcoeff_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devDSMINcoeff_YawTorque
    saveas(gca,['WBmod_devDSMINcoeff_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devDSMINcoeff_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devDSMINcoeff_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devUSMAXcoeff_YawTorque
    saveas(gca,['WBmod_devUSMAXcoeff_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devUSMAXcoeff_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devUSMAXcoeff_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_devUSMINcoeff_YawTorque
    saveas(gca,['WBmod_devUSMINcoeff_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devUSMINcoeff_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devUSMINcoeff_YawTorque_',num2str(n_YawTorque),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    cd ..
    
    end
end

%% Torque around right rotation axis
if calc_WB_TorqueAxisR == 1

%     calc_Torquebased_wbNORM_wingatt_TorqueAxisR
    calc_Torquebased_wbNORM_wingatt_TorqueAxisR1POS

    if plot_on ==1

    cd('WBmod_figs_torque')
    
xmin = -ceil(100*max(abs(M_R_mean_wb)))/100;
% xmin = -1.5;
xmax = -xmin;


figure
h=hist(M_R_mean_wb,50);
hist(M_R_mean_wb,50)
hold on
plot([M_R_limit_mod M_R_limit_mod],[0 max(h)],'-r')
% plot([-M_R_limit_mod -M_R_limit_mod],[0 max(h)],'-r')
title('Torque axis R Torque')
saveas(gca,['hist_TorqueAxisR',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
saveas(gca,['hist_TorqueAxisR',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
plot2svg(['hist_TorqueAxisR',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_heatmap_TorqueAxisR
    saveas(gca,['WBmod_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmodNsteady_TorqueAxisR
    saveas(gca,['WBmodNsteady_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmodNsteady_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmodNsteady_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_strokeMAXcoeff_TorqueAxisR
    saveas(gca,['WBmod_strokeMAXcoeff_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_strokeMAXcoeff_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_strokeMAXcoeff_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_strokeMINcoeff_TorqueAxisR
    saveas(gca,['WBmod_strokeMINcoeff_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_strokeMINcoeff_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_strokeMINcoeff_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_pitchMIDDScoeff_TorqueAxisR
    saveas(gca,['WBmod_pitchMIDDScoeff_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_pitchMIDDScoeff_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_pitchMIDDScoeff_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_pitchMIDUScoeff_TorqueAxisR
    saveas(gca,['WBmod_pitchMIDUScoeff_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_pitchMIDUScoeff_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_pitchMIDUScoeff_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devDSMAXcoeff_TorqueAxisR
    saveas(gca,['WBmod_devDSMAXcoeff_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devDSMAXcoeff_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devDSMAXcoeff_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devDSMINcoeff_TorqueAxisR
    saveas(gca,['WBmod_devDSMINcoeff_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devDSMINcoeff_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devDSMINcoeff_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    plot_WBmod_devUSMAXcoeff_TorqueAxisR
    saveas(gca,['WBmod_devUSMAXcoeff_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devUSMAXcoeff_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devUSMAXcoeff_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])
    
    plot_WBmod_devUSMINcoeff_TorqueAxisR
    saveas(gca,['WBmod_devUSMINcoeff_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.fig'])
    saveas(gca,['WBmod_devUSMINcoeff_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.png'])
    plot2svg(['WBmod_devUSMINcoeff_TorqueAxisR_',num2str(n_TorqueAxisR),'WBs_limit',num2str(limit_mod),'_norm',num2str(norm),'xSTD.svg'])

    cd ..
    end
end

%% save all
    save('WBdataset_all_steadyNmods_TorqueNorm.mat')

%% MS plots
% % load('WBdataset_all_steadyNmods_TorqueNorm.mat')
% mkdir('MSfigs_WBkin_torque')
% cd('MSfigs_WBkin_torque')
% 
% MSplot_WBfunc_heatmap_steady_NOtick
%     saveas(gca,['MSplot_WBfunc_heatmap_steady_NOtick.fig'])
%     saveas(gca,['MSplot_WBfunc_heatmap_steady_NOtick.png'])
%     plot2svg(['MSplot_WBfunc_heatmap_steady_NOtick.svg'])
% 
% % MSplot_WBmodsNsteady_onefig_2xNorm_INCrotAxisR
% MSplot_WBmodsNsteady_onefig_2xNorm_TORQUE_incAxisR
%     saveas(gca,['MSplot_WBmodsNsteady_onefig.fig'])
%     saveas(gca,['MSplot_WBmodsNsteady_onefig.png'])
%     plot2svg(['MSplot_WBmodsNsteady_onefig.svg'])
% 
% MSplot_WBsteady2mod2xNorm_Fenhance
%     saveas(gca,['MSplot_WBsteady2mod2xNorm_Fenhance.fig'])
%     saveas(gca,['MSplot_WBsteady2mod2xNorm_Fenhance.png'])
%     plot2svg(['MSplot_WBsteady2mod2xNorm_Fenhance.svg'])
% 
% MSplot_WBsteady2mod2xNorm_RollTorque
%     saveas(gca,['MSplot_WBsteady2mod2xNorm_RollTorque.fig'])
%     saveas(gca,['MSplot_WBsteady2mod2xNorm_RollTorque.png'])
%     plot2svg(['MSplot_WBsteady2mod2xNorm_RollTorque.svg'])
% 
% MSplot_WBsteady2mod2xNorm_TorqueAxisR
%     saveas(gca,['MSplot_WBsteady2mod2xNorm_TorqueAxisR.fig'])
%     saveas(gca,['MSplot_WBsteady2mod2xNorm_TorqueAxisR.png'])
%     plot2svg(['MSplot_WBsteady2mod2xNorm_TorqueAxisR.svg'])
% 
% MSplot_WBsteady2mod2xNorm_Torque
%     saveas(gca,['MSplot_WBsteady2mod2xNorm.fig'])
%     saveas(gca,['MSplot_WBsteady2mod2xNorm.png'])
%     plot2svg(['MSplot_WBsteady2mod2xNorm.svg'])
% 
% MSplot_WBsteady2mod2xNorm_1fig_incTorqueAxisR
%     saveas(gca,['MSplot_WBsteady2mod2xNorm_1fig.fig'])
%     saveas(gca,['MSplot_WBsteady2mod2xNorm_1fig.png'])
%     plot2svg(['MSplot_WBsteady2mod2xNorm_1fig.svg'])
% 
% cd ..


