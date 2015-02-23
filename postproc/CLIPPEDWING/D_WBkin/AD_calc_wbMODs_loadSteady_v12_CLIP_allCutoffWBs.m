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

calc_WB_d2ndMom = 1;
% calc_WB_d2ndMom = 0;

%% settings
% cut-offs
% AreaRatio_cut = 1-nanstd(AreaRatio);
% SecondMomentRatio_cut = 1-nanstd(SecondMomentRatio);
AreaRatio_cut = .9;
SecondMomentRatio_cut = .9;

AreaRatio_plot = .5;
SecondMomentRatio_plot = .5;

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

%% CALC CLIPPED WB MODS: AREA DIFFERENCE
if calc_WB_dArea == 1
    
% calc_wbNORM_wingatt_CLIPPED_dArea
calc_wbNORM_wingatt_CLIPPED_dArea_cutoff
% calc_wbNORM_wingatt_CLIPPED_dArea_TEclip

if plot_on == 1
mkdir('WBmod_figs_dArea')
cd('WBmod_figs_dArea')

xmin = floor(10*min(AreaRatio(:)))/10;
xmax = 1;

%% dArea distribution
figure
h=hist(AreaRatio_dArea,10);
hist(AreaRatio_dArea,10)
hold on
title('AreaRatio of clipped and intact wings')
xlabel('AreaRatio')
ylabel('# of wingbeats')
xlim([xmin xmax])

saveas(gca,['hist_AreaRatio_',num2str(n_dArea),'WBs.fig'])
saveas(gca,['hist_AreaRatio_',num2str(n_dArea),'WBs.png'])
plot2svg(['hist_AreaRatio_',num2str(n_dArea),'WBs.svg'])

%% body kin MODs
    plot_BODYmod_vel_hist_dArea
    saveas(gca,['WBmod_vel_dArea_',num2str(n_dArea),'WBs.fig'])
    saveas(gca,['WBmod_vel_dArea_',num2str(n_dArea),'WBs.png'])
    plot2svg(['WBmod_vel_dArea_',num2str(n_dArea),'WBs.svg'])

    plot_BODYmod_slip_hist_dArea
    saveas(gca,['WBmod_slip_dArea_',num2str(n_dArea),'WBs.fig'])
    saveas(gca,['WBmod_slip_dArea_',num2str(n_dArea),'WBs.png'])
    plot2svg(['WBmod_slip_dArea_',num2str(n_dArea),'WBs.svg'])

    plot_BODYmod_pitch_hist_dArea
    saveas(gca,['WBmod_pitch_dArea_',num2str(n_dArea),'WBs.fig'])
    saveas(gca,['WBmod_pitch_dArea_',num2str(n_dArea),'WBs.png'])
    plot2svg(['WBmod_pitch_dArea_',num2str(n_dArea),'WBs.svg'])

    plot_BODYmod_roll_hist_dArea
    saveas(gca,['WBmod_roll_dArea_',num2str(n_dArea),'WBs.fig'])
    saveas(gca,['WBmod_roll_dArea_',num2str(n_dArea),'WBs.png'])
    plot2svg(['WBmod_roll_dArea_',num2str(n_dArea),'WBs.svg'])

    plot_BODYmod_Fsp_pitch_hist_dArea
    saveas(gca,['WBmod_Fsp_pitch_dArea_',num2str(n_dArea),'WBs.fig'])
    saveas(gca,['WBmod_Fsp_pitch_dArea_',num2str(n_dArea),'WBs.png'])
    plot2svg(['WBmod_Fsp_pitch_dArea_',num2str(n_dArea),'WBs.svg'])

    plot_BODYmod_Fsp_roll_hist_dArea
    saveas(gca,['WBmod_Fsp_roll_dArea_',num2str(n_dArea),'WBs.fig'])
    saveas(gca,['WBmod_Fsp_roll_dArea_',num2str(n_dArea),'WBs.png'])
    plot2svg(['WBmod_Fsp_roll_dArea_',num2str(n_dArea),'WBs.svg'])

%% WB kin MODs
    plot_WBmod_freq_hist_dArea
    saveas(gca,['WBmod_freq_dArea_',num2str(n_dArea),'WBs.fig'])
    saveas(gca,['WBmod_freq_dArea_',num2str(n_dArea),'WBs.png'])
    plot2svg(['WBmod_freq_dArea_',num2str(n_dArea),'WBs.svg'])

%     plot_WBmod_heatmap_dArea
    plot_WBmod_heatmap_dArea_plotval
    saveas(gca,['WBmod_dArea_',num2str(n_dArea),'WBs_AreaRatio',num2str(AreaRatio_plot),'.fig'])
    saveas(gca,['WBmod_dArea_',num2str(n_dArea),'WBs_AreaRatio',num2str(AreaRatio_plot),'.png'])
    plot2svg(['WBmod_dArea_',num2str(n_dArea),'WBs_AreaRatio',num2str(AreaRatio_plot),'.svg'])
    
%     plot_WBmodNsteady_dArea
    plot_WBmodNsteady_dArea_plotval
    saveas(gca,['WBmodNsteady_dArea_',num2str(n_dArea),'WBs_AreaRatio',num2str(AreaRatio_plot),'.fig'])
    saveas(gca,['WBmodNsteady_dArea_',num2str(n_dArea),'WBs_AreaRatio',num2str(AreaRatio_plot),'.png'])
    plot2svg(['WBmodNsteady_dArea_',num2str(n_dArea),'WBs_AreaRatio',num2str(AreaRatio_plot),'.svg'])
    
    plot_WBmod_strokeMAXcoeff_dArea
    saveas(gca,['WBmod_strokeMAXcoeff_dArea_',num2str(n_dArea),'WBs.fig'])
    saveas(gca,['WBmod_strokeMAXcoeff_dArea_',num2str(n_dArea),'WBs.png'])
    plot2svg(['WBmod_strokeMAXcoeff_dArea_',num2str(n_dArea),'WBs.svg'])

    plot_WBmod_strokeMINcoeff_dArea
    saveas(gca,['WBmod_strokeMINcoeff_dArea_',num2str(n_dArea),'WBs.fig'])
    saveas(gca,['WBmod_strokeMINcoeff_dArea_',num2str(n_dArea),'WBs.png'])
    plot2svg(['WBmod_strokeMINcoeff_dArea_',num2str(n_dArea),'WBs.svg'])

    plot_WBmod_pitchMIDDScoeff_dArea
    saveas(gca,['WBmod_pitchMIDDScoeff_dArea_',num2str(n_dArea),'WBs.fig'])
    saveas(gca,['WBmod_pitchMIDDScoeff_dArea_',num2str(n_dArea),'WBs.png'])
    plot2svg(['WBmod_pitchMIDDScoeff_dArea_',num2str(n_dArea),'WBs.svg'])

    plot_WBmod_pitchMIDUScoeff_dArea
    saveas(gca,['WBmod_pitchMIDUScoeff_dArea_',num2str(n_dArea),'WBs.fig'])
    saveas(gca,['WBmod_pitchMIDUScoeff_dArea_',num2str(n_dArea),'WBs.png'])
    plot2svg(['WBmod_pitchMIDUScoeff_dArea_',num2str(n_dArea),'WBs.svg'])

    plot_WBmod_devDSMAXcoeff_dArea
    saveas(gca,['WBmod_devDSMAXcoeff_dArea_',num2str(n_dArea),'WBs.fig'])
    saveas(gca,['WBmod_devDSMAXcoeff_dArea_',num2str(n_dArea),'WBs.png'])
    plot2svg(['WBmod_devDSMAXcoeff_dArea_',num2str(n_dArea),'WBs.svg'])

    plot_WBmod_devDSMINcoeff_dArea
    saveas(gca,['WBmod_devDSMINcoeff_dArea_',num2str(n_dArea),'WBs.fig'])
    saveas(gca,['WBmod_devDSMINcoeff_dArea_',num2str(n_dArea),'WBs.png'])
    plot2svg(['WBmod_devDSMINcoeff_dArea_',num2str(n_dArea),'WBs.svg'])

    plot_WBmod_devUSMAXcoeff_dArea
    saveas(gca,['WBmod_devUSMAXcoeff_dArea_',num2str(n_dArea),'WBs.fig'])
    saveas(gca,['WBmod_devUSMAXcoeff_dArea_',num2str(n_dArea),'WBs.png'])
    plot2svg(['WBmod_devUSMAXcoeff_dArea_',num2str(n_dArea),'WBs.svg'])
    
    plot_WBmod_devUSMINcoeff_dArea
    saveas(gca,['WBmod_devUSMINcoeff_dArea_',num2str(n_dArea),'WBs.fig'])
    saveas(gca,['WBmod_devUSMINcoeff_dArea_',num2str(n_dArea),'WBs.png'])
    plot2svg(['WBmod_devUSMINcoeff_dArea_',num2str(n_dArea),'WBs.svg'])
cd ..
end
end

%% CALC CLIPPED WB MODS: 2nd MOMENT OF AREA DIFFERENCE
if calc_WB_d2ndMom == 1
    
% calc_wbNORM_wingatt_CLIPPED_d2ndMom
calc_wbNORM_wingatt_CLIPPED_d2ndMom_cutoff
% calc_wbNORM_wingatt_CLIPPED_d2ndMom_TipClip

if plot_on == 1
mkdir('WBmod_figs_d2ndMom')
cd('WBmod_figs_d2ndMom')

xmin = floor(10*min(SecondMomentRatio(:)))/10;
xmax = 1;

%% d2ndMom distribution
figure
h=hist(SecondMomentRatio_d2ndMom,10);
hist(SecondMomentRatio_d2ndMom,10)
hold on
title('SecondMomentRatio of clipped and intact wings')
xlabel('SecondMomentRatio')
ylabel('# of wingbeats')
xlim([xmin xmax])
saveas(gca,['hist_SecondMomentRatio_',num2str(n_d2ndMom),'WBs.fig'])
saveas(gca,['hist_SecondMomentRatio_',num2str(n_d2ndMom),'WBs.png'])
plot2svg(['hist_SecondMomentRatio_',num2str(n_d2ndMom),'WBs.svg'])

%% body kin MODs
    plot_BODYmod_vel_hist_d2ndMom
    saveas(gca,['WBmod_vel_d2ndMom_',num2str(n_d2ndMom),'WBs.fig'])
    saveas(gca,['WBmod_vel_d2ndMom_',num2str(n_d2ndMom),'WBs.png'])
    plot2svg(['WBmod_vel_d2ndMom_',num2str(n_d2ndMom),'WBs.svg'])

    plot_BODYmod_slip_hist_d2ndMom
    saveas(gca,['WBmod_slip_d2ndMom_',num2str(n_d2ndMom),'WBs.fig'])
    saveas(gca,['WBmod_slip_d2ndMom_',num2str(n_d2ndMom),'WBs.png'])
    plot2svg(['WBmod_slip_d2ndMom_',num2str(n_d2ndMom),'WBs.svg'])

    plot_BODYmod_pitch_hist_d2ndMom
    saveas(gca,['WBmod_pitch_d2ndMom_',num2str(n_d2ndMom),'WBs.fig'])
    saveas(gca,['WBmod_pitch_d2ndMom_',num2str(n_d2ndMom),'WBs.png'])
    plot2svg(['WBmod_pitch_d2ndMom_',num2str(n_d2ndMom),'WBs.svg'])

    plot_BODYmod_roll_hist_d2ndMom
    saveas(gca,['WBmod_roll_d2ndMom_',num2str(n_d2ndMom),'WBs.fig'])
    saveas(gca,['WBmod_roll_d2ndMom_',num2str(n_d2ndMom),'WBs.png'])
    plot2svg(['WBmod_roll_d2ndMom_',num2str(n_d2ndMom),'WBs.svg'])

    plot_BODYmod_Fsp_pitch_hist_d2ndMom
    saveas(gca,['WBmod_Fsp_pitch_d2ndMom_',num2str(n_d2ndMom),'WBs.fig'])
    saveas(gca,['WBmod_Fsp_pitch_d2ndMom_',num2str(n_d2ndMom),'WBs.png'])
    plot2svg(['WBmod_Fsp_pitch_d2ndMom_',num2str(n_d2ndMom),'WBs.svg'])

    plot_BODYmod_Fsp_roll_hist_d2ndMom
    saveas(gca,['WBmod_Fsp_roll_d2ndMom_',num2str(n_d2ndMom),'WBs.fig'])
    saveas(gca,['WBmod_Fsp_roll_d2ndMom_',num2str(n_d2ndMom),'WBs.png'])
    plot2svg(['WBmod_Fsp_roll_d2ndMom_',num2str(n_d2ndMom),'WBs.svg'])

%% WB kin MODs
    plot_WBmod_freq_hist_d2ndMom
    saveas(gca,['WBmod_freq_d2ndMom_',num2str(n_d2ndMom),'WBs.fig'])
    saveas(gca,['WBmod_freq_d2ndMom_',num2str(n_d2ndMom),'WBs.png'])
    plot2svg(['WBmod_freq_d2ndMom_',num2str(n_d2ndMom),'WBs.svg'])

%     plot_WBmod_heatmap_d2ndMom
    plot_WBmod_heatmap_d2ndMom_plotval
    saveas(gca,['WBmod_d2ndMom_',num2str(n_d2ndMom),'WBs_2ndMomRatio',num2str(SecondMomentRatio_plot),'.fig'])
    saveas(gca,['WBmod_d2ndMom_',num2str(n_d2ndMom),'WBs_2ndMomRatio',num2str(SecondMomentRatio_plot),'.png'])
    plot2svg(['WBmod_d2ndMom_',num2str(n_d2ndMom),'WBs_2ndMomRatio',num2str(SecondMomentRatio_plot),'.svg'])
    
%     plot_WBmodNsteady_d2ndMom
    plot_WBmodNsteady_d2ndMom_plotval
    saveas(gca,['WBmodNsteady_d2ndMom_',num2str(n_d2ndMom),'WBs_2ndMomRatio',num2str(SecondMomentRatio_plot),'.fig'])
    saveas(gca,['WBmodNsteady_d2ndMom_',num2str(n_d2ndMom),'WBs_2ndMomRatio',num2str(SecondMomentRatio_plot),'.png'])
    plot2svg(['WBmodNsteady_d2ndMom_',num2str(n_d2ndMom),'WBs_2ndMomRatio',num2str(SecondMomentRatio_plot),'.svg'])
    
    plot_WBmod_strokeMAXcoeff_d2ndMom
    saveas(gca,['WBmod_strokeMAXcoeff_d2ndMom_',num2str(n_d2ndMom),'WBs.fig'])
    saveas(gca,['WBmod_strokeMAXcoeff_d2ndMom_',num2str(n_d2ndMom),'WBs.png'])
    plot2svg(['WBmod_strokeMAXcoeff_d2ndMom_',num2str(n_d2ndMom),'WBs.svg'])

    plot_WBmod_strokeMINcoeff_d2ndMom
    saveas(gca,['WBmod_strokeMINcoeff_d2ndMom_',num2str(n_d2ndMom),'WBs.fig'])
    saveas(gca,['WBmod_strokeMINcoeff_d2ndMom_',num2str(n_d2ndMom),'WBs.png'])
    plot2svg(['WBmod_strokeMINcoeff_d2ndMom_',num2str(n_d2ndMom),'WBs.svg'])

    plot_WBmod_pitchMIDDScoeff_d2ndMom
    saveas(gca,['WBmod_pitchMIDDScoeff_d2ndMom_',num2str(n_d2ndMom),'WBs.fig'])
    saveas(gca,['WBmod_pitchMIDDScoeff_d2ndMom_',num2str(n_d2ndMom),'WBs.png'])
    plot2svg(['WBmod_pitchMIDDScoeff_d2ndMom_',num2str(n_d2ndMom),'WBs.svg'])

    plot_WBmod_pitchMIDUScoeff_d2ndMom
    saveas(gca,['WBmod_pitchMIDUScoeff_d2ndMom_',num2str(n_d2ndMom),'WBs.fig'])
    saveas(gca,['WBmod_pitchMIDUScoeff_d2ndMom_',num2str(n_d2ndMom),'WBs.png'])
    plot2svg(['WBmod_pitchMIDUScoeff_d2ndMom_',num2str(n_d2ndMom),'WBs.svg'])

    plot_WBmod_devDSMAXcoeff_d2ndMom
    saveas(gca,['WBmod_devDSMAXcoeff_d2ndMom_',num2str(n_d2ndMom),'WBs.fig'])
    saveas(gca,['WBmod_devDSMAXcoeff_d2ndMom_',num2str(n_d2ndMom),'WBs.png'])
    plot2svg(['WBmod_devDSMAXcoeff_d2ndMom_',num2str(n_d2ndMom),'WBs.svg'])

    plot_WBmod_devDSMINcoeff_d2ndMom
    saveas(gca,['WBmod_devDSMINcoeff_d2ndMom_',num2str(n_d2ndMom),'WBs.fig'])
    saveas(gca,['WBmod_devDSMINcoeff_d2ndMom_',num2str(n_d2ndMom),'WBs.png'])
    plot2svg(['WBmod_devDSMINcoeff_d2ndMom_',num2str(n_d2ndMom),'WBs.svg'])

    plot_WBmod_devUSMAXcoeff_d2ndMom
    saveas(gca,['WBmod_devUSMAXcoeff_d2ndMom_',num2str(n_d2ndMom),'WBs.fig'])
    saveas(gca,['WBmod_devUSMAXcoeff_d2ndMom_',num2str(n_d2ndMom),'WBs.png'])
    plot2svg(['WBmod_devUSMAXcoeff_d2ndMom_',num2str(n_d2ndMom),'WBs.svg'])
    
    plot_WBmod_devUSMINcoeff_d2ndMom
    saveas(gca,['WBmod_devUSMINcoeff_d2ndMom_',num2str(n_d2ndMom),'WBs.fig'])
    saveas(gca,['WBmod_devUSMINcoeff_d2ndMom_',num2str(n_d2ndMom),'WBs.png'])
    plot2svg(['WBmod_devUSMINcoeff_d2ndMom_',num2str(n_d2ndMom),'WBs.svg'])
cd ..
end
end

%% save all
if save_on == 1
    save('WBdataset_all_steadyNclipMods.mat')
end


