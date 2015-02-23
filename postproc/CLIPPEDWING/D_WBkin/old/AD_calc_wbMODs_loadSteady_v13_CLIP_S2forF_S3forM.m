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

% plot_fourier = 1;
plot_fourier = 0;

plot_fits = 1;
% plot_fits = 0;

calc_WB_dArea = 1;
% calc_WB_dArea = 0;

calc_WB_dS2 = 1;
% calc_WB_dS2 = 0;

%% settings
% cut-offs
% AreaRatio_cut = 1-nanstd(AreaRatio);
% S2Ratio_cut = 1-nanstd(S2Ratio);
AreaRatio_cut = .9;
S2Ratio_cut = .9;

AreaRatio_plot = .5;
S2Ratio_plot = .5;

% fourier orders
MOD_fourier_order = 8;

% number of polinomials
n_pol_MOD = 10;

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
% steady body kin data
vel_steady = V_steady_meanCIstd(:,1);
pitch_body_steady = pitch_global_steady_meanCIstd(:,1);

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
    
%% CALC CLIPPED WB MODS: 2nd MOMENT OF AREA DIFFERENCE for Forces
if calc_WB_dS2forF == 1
    
calc_wbNORM_wingatt_CLIPPED_dS2forF_cutoff

if plot_on == 1
mkdir('WBmod_figs_dS2')
cd('WBmod_figs_dS2')

xmin = floor(10*min(S2Ratio(:)))/10;
xmax = 1;

%% dS2 distribution
figure
h=hist(S2Ratio_dS2,10);
hist(S2Ratio_dS2,10)
hold on
title('S2Ratio of clipped and intact wings')
xlabel('S2Ratio')
ylabel('# of wingbeats')
xlim([xmin xmax])
saveas(gca,['hist_S2Ratio_',num2str(n_dS2),'WBs.fig'])
saveas(gca,['hist_S2Ratio_',num2str(n_dS2),'WBs.png'])
plot2svg(['hist_S2Ratio_',num2str(n_dS2),'WBs.svg'])

%% body kin MODs
    plot_BODYmod_vel_hist_dS2
    saveas(gca,['WBmod_vel_dS2_',num2str(n_dS2),'WBs.fig'])
    saveas(gca,['WBmod_vel_dS2_',num2str(n_dS2),'WBs.png'])
    plot2svg(['WBmod_vel_dS2_',num2str(n_dS2),'WBs.svg'])

    plot_BODYmod_slip_hist_dS2
    saveas(gca,['WBmod_slip_dS2_',num2str(n_dS2),'WBs.fig'])
    saveas(gca,['WBmod_slip_dS2_',num2str(n_dS2),'WBs.png'])
    plot2svg(['WBmod_slip_dS2_',num2str(n_dS2),'WBs.svg'])

    plot_BODYmod_pitch_hist_dS2
    saveas(gca,['WBmod_pitch_dS2_',num2str(n_dS2),'WBs.fig'])
    saveas(gca,['WBmod_pitch_dS2_',num2str(n_dS2),'WBs.png'])
    plot2svg(['WBmod_pitch_dS2_',num2str(n_dS2),'WBs.svg'])

    plot_BODYmod_roll_hist_dS2
    saveas(gca,['WBmod_roll_dS2_',num2str(n_dS2),'WBs.fig'])
    saveas(gca,['WBmod_roll_dS2_',num2str(n_dS2),'WBs.png'])
    plot2svg(['WBmod_roll_dS2_',num2str(n_dS2),'WBs.svg'])

    plot_BODYmod_Fsp_pitch_hist_dS2
    saveas(gca,['WBmod_Fsp_pitch_dS2_',num2str(n_dS2),'WBs.fig'])
    saveas(gca,['WBmod_Fsp_pitch_dS2_',num2str(n_dS2),'WBs.png'])
    plot2svg(['WBmod_Fsp_pitch_dS2_',num2str(n_dS2),'WBs.svg'])

    plot_BODYmod_Fsp_roll_hist_dS2
    saveas(gca,['WBmod_Fsp_roll_dS2_',num2str(n_dS2),'WBs.fig'])
    saveas(gca,['WBmod_Fsp_roll_dS2_',num2str(n_dS2),'WBs.png'])
    plot2svg(['WBmod_Fsp_roll_dS2_',num2str(n_dS2),'WBs.svg'])

%% WB kin MODs
    plot_WBmod_freq_hist_dS2
    saveas(gca,['WBmod_freq_dS2_',num2str(n_dS2),'WBs.fig'])
    saveas(gca,['WBmod_freq_dS2_',num2str(n_dS2),'WBs.png'])
    plot2svg(['WBmod_freq_dS2_',num2str(n_dS2),'WBs.svg'])

%     plot_WBmod_heatmap_dS2
    plot_WBmod_heatmap_dS2_plotval
    saveas(gca,['WBmod_dS2_',num2str(n_dS2),'WBs_2ndMomRatio',num2str(S2Ratio_plot),'.fig'])
    saveas(gca,['WBmod_dS2_',num2str(n_dS2),'WBs_2ndMomRatio',num2str(S2Ratio_plot),'.png'])
    plot2svg(['WBmod_dS2_',num2str(n_dS2),'WBs_2ndMomRatio',num2str(S2Ratio_plot),'.svg'])
    
%     plot_WBmodNsteady_dS2
    plot_WBmodNsteady_dS2_plotval
    saveas(gca,['WBmodNsteady_dS2_',num2str(n_dS2),'WBs_2ndMomRatio',num2str(S2Ratio_plot),'.fig'])
    saveas(gca,['WBmodNsteady_dS2_',num2str(n_dS2),'WBs_2ndMomRatio',num2str(S2Ratio_plot),'.png'])
    plot2svg(['WBmodNsteady_dS2_',num2str(n_dS2),'WBs_2ndMomRatio',num2str(S2Ratio_plot),'.svg'])
    
    plot_WBmod_strokeMAXcoeff_dS2
    saveas(gca,['WBmod_strokeMAXcoeff_dS2_',num2str(n_dS2),'WBs.fig'])
    saveas(gca,['WBmod_strokeMAXcoeff_dS2_',num2str(n_dS2),'WBs.png'])
    plot2svg(['WBmod_strokeMAXcoeff_dS2_',num2str(n_dS2),'WBs.svg'])

    plot_WBmod_strokeMINcoeff_dS2
    saveas(gca,['WBmod_strokeMINcoeff_dS2_',num2str(n_dS2),'WBs.fig'])
    saveas(gca,['WBmod_strokeMINcoeff_dS2_',num2str(n_dS2),'WBs.png'])
    plot2svg(['WBmod_strokeMINcoeff_dS2_',num2str(n_dS2),'WBs.svg'])

    plot_WBmod_pitchMIDDScoeff_dS2
    saveas(gca,['WBmod_pitchMIDDScoeff_dS2_',num2str(n_dS2),'WBs.fig'])
    saveas(gca,['WBmod_pitchMIDDScoeff_dS2_',num2str(n_dS2),'WBs.png'])
    plot2svg(['WBmod_pitchMIDDScoeff_dS2_',num2str(n_dS2),'WBs.svg'])

    plot_WBmod_pitchMIDUScoeff_dS2
    saveas(gca,['WBmod_pitchMIDUScoeff_dS2_',num2str(n_dS2),'WBs.fig'])
    saveas(gca,['WBmod_pitchMIDUScoeff_dS2_',num2str(n_dS2),'WBs.png'])
    plot2svg(['WBmod_pitchMIDUScoeff_dS2_',num2str(n_dS2),'WBs.svg'])

    plot_WBmod_devDSMAXcoeff_dS2
    saveas(gca,['WBmod_devDSMAXcoeff_dS2_',num2str(n_dS2),'WBs.fig'])
    saveas(gca,['WBmod_devDSMAXcoeff_dS2_',num2str(n_dS2),'WBs.png'])
    plot2svg(['WBmod_devDSMAXcoeff_dS2_',num2str(n_dS2),'WBs.svg'])

    plot_WBmod_devDSMINcoeff_dS2
    saveas(gca,['WBmod_devDSMINcoeff_dS2_',num2str(n_dS2),'WBs.fig'])
    saveas(gca,['WBmod_devDSMINcoeff_dS2_',num2str(n_dS2),'WBs.png'])
    plot2svg(['WBmod_devDSMINcoeff_dS2_',num2str(n_dS2),'WBs.svg'])

    plot_WBmod_devUSMAXcoeff_dS2
    saveas(gca,['WBmod_devUSMAXcoeff_dS2_',num2str(n_dS2),'WBs.fig'])
    saveas(gca,['WBmod_devUSMAXcoeff_dS2_',num2str(n_dS2),'WBs.png'])
    plot2svg(['WBmod_devUSMAXcoeff_dS2_',num2str(n_dS2),'WBs.svg'])
    
    plot_WBmod_devUSMINcoeff_dS2
    saveas(gca,['WBmod_devUSMINcoeff_dS2_',num2str(n_dS2),'WBs.fig'])
    saveas(gca,['WBmod_devUSMINcoeff_dS2_',num2str(n_dS2),'WBs.png'])
    plot2svg(['WBmod_devUSMINcoeff_dS2_',num2str(n_dS2),'WBs.svg'])
cd ..
end
end

%% CALC CLIPPED WB MODS: 3nd MOMENT OF AREA DIFFERENCE for Torques
if calc_WB_dS3forM == 1
    
calc_wbNORM_wingatt_CLIPPED_dS3forM_cutoff

if plot_on == 1
mkdir('WBmod_figs_dS3')
cd('WBmod_figs_dS3')

xmin = floor(10*min(S3Ratio(:)))/10;
xmax = 1;

%% dS3 distribution
figure
h=hist(S3Ratio_dS3,10);
hist(S3Ratio_dS3,10)
hold on
title('S3Ratio of clipped and intact wings')
xlabel('S3Ratio')
ylabel('# of wingbeats')
xlim([xmin xmax])
saveas(gca,['hist_S3Ratio_',num2str(n_dS3),'WBs.fig'])
saveas(gca,['hist_S3Ratio_',num2str(n_dS3),'WBs.png'])
plot2svg(['hist_S3Ratio_',num2str(n_dS3),'WBs.svg'])

%% body kin MODs
    plot_BODYmod_vel_hist_dS3
    saveas(gca,['WBmod_vel_dS3_',num2str(n_dS3),'WBs.fig'])
    saveas(gca,['WBmod_vel_dS3_',num2str(n_dS3),'WBs.png'])
    plot2svg(['WBmod_vel_dS3_',num2str(n_dS3),'WBs.svg'])

    plot_BODYmod_slip_hist_dS3
    saveas(gca,['WBmod_slip_dS3_',num2str(n_dS3),'WBs.fig'])
    saveas(gca,['WBmod_slip_dS3_',num2str(n_dS3),'WBs.png'])
    plot2svg(['WBmod_slip_dS3_',num2str(n_dS3),'WBs.svg'])

    plot_BODYmod_pitch_hist_dS3
    saveas(gca,['WBmod_pitch_dS3_',num2str(n_dS3),'WBs.fig'])
    saveas(gca,['WBmod_pitch_dS3_',num2str(n_dS3),'WBs.png'])
    plot2svg(['WBmod_pitch_dS3_',num2str(n_dS3),'WBs.svg'])

    plot_BODYmod_roll_hist_dS3
    saveas(gca,['WBmod_roll_dS3_',num2str(n_dS3),'WBs.fig'])
    saveas(gca,['WBmod_roll_dS3_',num2str(n_dS3),'WBs.png'])
    plot2svg(['WBmod_roll_dS3_',num2str(n_dS3),'WBs.svg'])

    plot_BODYmod_Fsp_pitch_hist_dS3
    saveas(gca,['WBmod_Fsp_pitch_dS3_',num2str(n_dS3),'WBs.fig'])
    saveas(gca,['WBmod_Fsp_pitch_dS3_',num2str(n_dS3),'WBs.png'])
    plot2svg(['WBmod_Fsp_pitch_dS3_',num2str(n_dS3),'WBs.svg'])

    plot_BODYmod_Fsp_roll_hist_dS3
    saveas(gca,['WBmod_Fsp_roll_dS3_',num2str(n_dS3),'WBs.fig'])
    saveas(gca,['WBmod_Fsp_roll_dS3_',num2str(n_dS3),'WBs.png'])
    plot2svg(['WBmod_Fsp_roll_dS3_',num2str(n_dS3),'WBs.svg'])

%% WB kin MODs
    plot_WBmod_freq_hist_dS3
    saveas(gca,['WBmod_freq_dS3_',num2str(n_dS3),'WBs.fig'])
    saveas(gca,['WBmod_freq_dS3_',num2str(n_dS3),'WBs.png'])
    plot2svg(['WBmod_freq_dS3_',num2str(n_dS3),'WBs.svg'])

%     plot_WBmod_heatmap_dS3
    plot_WBmod_heatmap_dS3_plotval
    saveas(gca,['WBmod_dS3_',num2str(n_dS3),'WBs_2ndMomRatio',num2str(S3Ratio_plot),'.fig'])
    saveas(gca,['WBmod_dS3_',num2str(n_dS3),'WBs_2ndMomRatio',num2str(S3Ratio_plot),'.png'])
    plot2svg(['WBmod_dS3_',num2str(n_dS3),'WBs_2ndMomRatio',num2str(S3Ratio_plot),'.svg'])
    
%     plot_WBmodNsteady_dS3
    plot_WBmodNsteady_dS3_plotval
    saveas(gca,['WBmodNsteady_dS3_',num2str(n_dS3),'WBs_2ndMomRatio',num2str(S3Ratio_plot),'.fig'])
    saveas(gca,['WBmodNsteady_dS3_',num2str(n_dS3),'WBs_2ndMomRatio',num2str(S3Ratio_plot),'.png'])
    plot2svg(['WBmodNsteady_dS3_',num2str(n_dS3),'WBs_2ndMomRatio',num2str(S3Ratio_plot),'.svg'])
    
    plot_WBmod_strokeMAXcoeff_dS3
    saveas(gca,['WBmod_strokeMAXcoeff_dS3_',num2str(n_dS3),'WBs.fig'])
    saveas(gca,['WBmod_strokeMAXcoeff_dS3_',num2str(n_dS3),'WBs.png'])
    plot2svg(['WBmod_strokeMAXcoeff_dS3_',num2str(n_dS3),'WBs.svg'])

    plot_WBmod_strokeMINcoeff_dS3
    saveas(gca,['WBmod_strokeMINcoeff_dS3_',num2str(n_dS3),'WBs.fig'])
    saveas(gca,['WBmod_strokeMINcoeff_dS3_',num2str(n_dS3),'WBs.png'])
    plot2svg(['WBmod_strokeMINcoeff_dS3_',num2str(n_dS3),'WBs.svg'])

    plot_WBmod_pitchMIDDScoeff_dS3
    saveas(gca,['WBmod_pitchMIDDScoeff_dS3_',num2str(n_dS3),'WBs.fig'])
    saveas(gca,['WBmod_pitchMIDDScoeff_dS3_',num2str(n_dS3),'WBs.png'])
    plot2svg(['WBmod_pitchMIDDScoeff_dS3_',num2str(n_dS3),'WBs.svg'])

    plot_WBmod_pitchMIDUScoeff_dS3
    saveas(gca,['WBmod_pitchMIDUScoeff_dS3_',num2str(n_dS3),'WBs.fig'])
    saveas(gca,['WBmod_pitchMIDUScoeff_dS3_',num2str(n_dS3),'WBs.png'])
    plot2svg(['WBmod_pitchMIDUScoeff_dS3_',num2str(n_dS3),'WBs.svg'])

    plot_WBmod_devDSMAXcoeff_dS3
    saveas(gca,['WBmod_devDSMAXcoeff_dS3_',num2str(n_dS3),'WBs.fig'])
    saveas(gca,['WBmod_devDSMAXcoeff_dS3_',num2str(n_dS3),'WBs.png'])
    plot2svg(['WBmod_devDSMAXcoeff_dS3_',num2str(n_dS3),'WBs.svg'])

    plot_WBmod_devDSMINcoeff_dS3
    saveas(gca,['WBmod_devDSMINcoeff_dS3_',num2str(n_dS3),'WBs.fig'])
    saveas(gca,['WBmod_devDSMINcoeff_dS3_',num2str(n_dS3),'WBs.png'])
    plot2svg(['WBmod_devDSMINcoeff_dS3_',num2str(n_dS3),'WBs.svg'])

    plot_WBmod_devUSMAXcoeff_dS3
    saveas(gca,['WBmod_devUSMAXcoeff_dS3_',num2str(n_dS3),'WBs.fig'])
    saveas(gca,['WBmod_devUSMAXcoeff_dS3_',num2str(n_dS3),'WBs.png'])
    plot2svg(['WBmod_devUSMAXcoeff_dS3_',num2str(n_dS3),'WBs.svg'])
    
    plot_WBmod_devUSMINcoeff_dS3
    saveas(gca,['WBmod_devUSMINcoeff_dS3_',num2str(n_dS3),'WBs.fig'])
    saveas(gca,['WBmod_devUSMINcoeff_dS3_',num2str(n_dS3),'WBs.png'])
    plot2svg(['WBmod_devUSMINcoeff_dS3_',num2str(n_dS3),'WBs.svg'])
cd ..
end
end

%% save all
if save_on == 1
    save('WBdataset_all_steadyNclipMods.mat')
end


