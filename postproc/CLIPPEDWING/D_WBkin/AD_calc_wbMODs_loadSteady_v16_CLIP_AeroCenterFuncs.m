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

% calc_WB_AeroCenterFuncClipped = 1;
calc_WB_AeroCenterFuncClipped = 0;

% calc_WB_AeroCenterFuncIntact = 1;
calc_WB_AeroCenterFuncIntact = 0;

% calc_WB_AeroCenterFuncRatio = 1;
calc_WB_AeroCenterFuncRatio = 0;

calc_WB_AeroCenterFuncs = 1;
% calc_WB_AeroCenterFuncs = 0;

%% settings

% calc AeroCenterFuncs
AeroCenterFuncClipped = 1./(FirstMomentRatio+1);
AeroCenterFuncIntact = FirstMomentRatio./(FirstMomentRatio+1);
AeroCenterFuncRatio = 1./FirstMomentRatio;

% steady flight vals
FirstMomentRatio_steady = 1;
AeroCenterFuncClipped_steady = 1./(FirstMomentRatio_steady+1);
AeroCenterFuncIntact_steady = FirstMomentRatio_steady./(FirstMomentRatio_steady+1);
AeroCenterFuncRatio_steady = 1./FirstMomentRatio_steady;

% cut-offs
% FirstMomentRatio_cut = 0.9;
FirstMomentRatio_cut = 1-nanstd(FirstMomentRatio);
% FirstMomentRatio_cut = 1-2*nanstd(FirstMomentRatio);
AeroCenterFuncClipped_cut = 1./(FirstMomentRatio_cut+1);
AeroCenterFuncIntact_cut = FirstMomentRatio_cut./(FirstMomentRatio_cut+1);
AeroCenterFuncRatio_cut = 1./FirstMomentRatio_cut;

% plot vals
FirstMomentRatio_plot = 0.5;
AeroCenterFuncClipped_plot = 1./(FirstMomentRatio_plot+1);
AeroCenterFuncIntact_plot = FirstMomentRatio_plot./(FirstMomentRatio_plot+1);
AeroCenterFuncRatio_plot = 1./FirstMomentRatio_plot;

if plot_on == 1
    mkdir('WBmod_figs_AeroCenterFuncs')
    cd('WBmod_figs_AeroCenterFuncs')

    figure
    subplot(2,2,1)
    title('FirstMomentRatio')
    hold on
    hist(FirstMomentRatio)
    plot([FirstMomentRatio_steady FirstMomentRatio_steady],[0 200],'-ob')
    plot([FirstMomentRatio_cut FirstMomentRatio_cut],[0 200],'-og')
    plot([FirstMomentRatio_plot FirstMomentRatio_plot],[0 200],'-or')

    subplot(2,2,2)
    title('AeroCenterFunc Clipped wing')
    hold on
    hist(AeroCenterFuncClipped)
    plot([AeroCenterFuncClipped_steady AeroCenterFuncClipped_steady],[0 200],'-ob')
    plot([AeroCenterFuncClipped_cut AeroCenterFuncClipped_cut],[0 200],'-og')
    plot([AeroCenterFuncClipped_plot AeroCenterFuncClipped_plot],[0 200],'-or')

    subplot(2,2,3)
    title('AeroCenterFunc Intact wing')
    hold on
    hist(AeroCenterFuncIntact)
    plot([AeroCenterFuncIntact_steady AeroCenterFuncIntact_steady],[0 200],'-ob')
    plot([AeroCenterFuncIntact_cut AeroCenterFuncIntact_cut],[0 200],'-og')
    plot([AeroCenterFuncIntact_plot AeroCenterFuncIntact_plot],[0 200],'-or')

    subplot(2,2,4)
    title('AeroCenterFunc Ratio')
    hold on
    hist(AeroCenterFuncRatio)
    plot([AeroCenterFuncRatio_steady AeroCenterFuncRatio_steady],[0 200],'-ob')
    plot([AeroCenterFuncRatio_cut AeroCenterFuncRatio_cut],[0 200],'-og')
    plot([AeroCenterFuncRatio_plot AeroCenterFuncRatio_plot],[0 200],'-or')

    saveas(gca,['hist_AeroCenterFunc_allWBs.fig'])
    saveas(gca,['hist_AeroCenterFunc_allWBs.png'])
    plot2svg(['hist_AeroCenterFunc_allWBs.svg'])
    
    cd ..
end

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

%% steady body kin data
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

%% CALC CLIPPED WB MODS: Correlate WB kin with AeroCenter Functions
% filter wingbeats: FirstMomentRatio_cut = 1 - std(FirstMomentRatio)
% body kin: correlation with AeroCenterFuncRatio =1/r*
% clipped wing kin: correlation with AeroCenterFuncClipped = 1/(r*+1)
% intact wing kin: correlation with AeroCenterFuncIntact = r*/(r*+1)
% kin difference (intact wing - clipped wing): correlation with AeroCenterFuncRatio =1/r*

if calc_WB_AeroCenterFuncs == 1
    
%     calc_wbNORM_wingatt_CLIPPED_S1S2S3_AeroCenterFuncs_cutoff
    calc_wbNORM_wingatt_CLIPPED_S1S2S3_AeroCenterFuncs_NOcircMean

if plot_on == 1
mkdir('WBmod_figs_AeroCenterFuncs')
cd('WBmod_figs_AeroCenterFuncs')


%% AeroCenterFuncs distribution

figure
subplot(2,2,1)
title('FirstMomentRatio')
xlabel('AeroCenterFuncRatio')
ylabel('# of wingbeats')
xmin = floor(10*min(FirstMomentRatio_AeroCenterFuncs(:)))/10;
xmax = ceil(10*max(FirstMomentRatio_AeroCenterFuncs(:)))/10;
xlim([xmin xmax])
hold on
hist(FirstMomentRatio_AeroCenterFuncs)
plot([FirstMomentRatio_steady FirstMomentRatio_steady],[0 200],'-ob')
plot([FirstMomentRatio_cut FirstMomentRatio_cut],[0 200],'-og')
plot([FirstMomentRatio_plot FirstMomentRatio_plot],[0 200],'-or')

subplot(2,2,2)
title('AeroCenterFunc Clipped wing')
xlabel('AeroCenterFuncRatio')
ylabel('# of wingbeats')
xmin = floor(10*min(AeroCenterFuncClipped_AeroCenterFuncs(:)))/10;
xmax = ceil(10*max(AeroCenterFuncClipped_AeroCenterFuncs(:)))/10;
xlim([xmin xmax])
hold on
hist(AeroCenterFuncClipped_AeroCenterFuncs)
plot([AeroCenterFuncClipped_steady AeroCenterFuncClipped_steady],[0 200],'-ob')
plot([AeroCenterFuncClipped_cut AeroCenterFuncClipped_cut],[0 200],'-og')
plot([AeroCenterFuncClipped_plot AeroCenterFuncClipped_plot],[0 200],'-or')

subplot(2,2,3)
title('AeroCenterFunc Intact wing')
xlabel('AeroCenterFuncRatio')
ylabel('# of wingbeats')
xmin = floor(10*min(AeroCenterFuncIntact_AeroCenterFuncs(:)))/10;
xmax = ceil(10*max(AeroCenterFuncIntact_AeroCenterFuncs(:)))/10;
xlim([xmin xmax])
hold on
hist(AeroCenterFuncIntact_AeroCenterFuncs)
plot([AeroCenterFuncIntact_steady AeroCenterFuncIntact_steady],[0 200],'-ob')
plot([AeroCenterFuncIntact_cut AeroCenterFuncIntact_cut],[0 200],'-og')
plot([AeroCenterFuncIntact_plot AeroCenterFuncIntact_plot],[0 200],'-or')

subplot(2,2,4)
title('AeroCenterFunc Ratio')
xlabel('AeroCenterFuncRatio')
ylabel('# of wingbeats')
xmin = floor(10*min(AeroCenterFuncRatio_AeroCenterFuncs(:)))/10;
xmax = ceil(10*max(AeroCenterFuncRatio_AeroCenterFuncs(:)))/10;
xlim([xmin xmax])
hold on
hist(AeroCenterFuncRatio_AeroCenterFuncs)
plot([AeroCenterFuncRatio_steady AeroCenterFuncRatio_steady],[0 200],'-ob')
plot([AeroCenterFuncRatio_cut AeroCenterFuncRatio_cut],[0 200],'-og')
plot([AeroCenterFuncRatio_plot AeroCenterFuncRatio_plot],[0 200],'-or')

saveas(gca,['hist_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.fig'])
saveas(gca,['hist_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.png'])
plot2svg(['hist_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.svg'])

%% body kin MODs
    plot_BODYmod_vel_hist_AeroCenterFuncs
    saveas(gca,['WBmod_vel_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.fig'])
    saveas(gca,['WBmod_vel_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.png'])
    plot2svg(['WBmod_vel_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.svg'])

    plot_BODYmod_slip_hist_AeroCenterFuncs
    saveas(gca,['WBmod_slip_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.fig'])
    saveas(gca,['WBmod_slip_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.png'])
    plot2svg(['WBmod_slip_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.svg'])

    plot_BODYmod_pitch_hist_AeroCenterFuncs
    saveas(gca,['WBmod_pitch_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.fig'])
    saveas(gca,['WBmod_pitch_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.png'])
    plot2svg(['WBmod_pitch_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.svg'])

    plot_BODYmod_roll_hist_AeroCenterFuncs
    saveas(gca,['WBmod_roll_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.fig'])
    saveas(gca,['WBmod_roll_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.png'])
    plot2svg(['WBmod_roll_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.svg'])

    plot_BODYmod_Fsp_pitch_hist_AeroCenterFuncs
    saveas(gca,['WBmod_Fsp_pitch_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.fig'])
    saveas(gca,['WBmod_Fsp_pitch_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.png'])
    plot2svg(['WBmod_Fsp_pitch_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.svg'])

    plot_BODYmod_Fsp_roll_hist_AeroCenterFuncs
    saveas(gca,['WBmod_Fsp_roll_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.fig'])
    saveas(gca,['WBmod_Fsp_roll_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.png'])
    plot2svg(['WBmod_Fsp_roll_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.svg'])

%% WB kin MODs
    plot_WBmod_freq_hist_AeroCenterFuncs
    saveas(gca,['WBmod_freq_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.fig'])
    saveas(gca,['WBmod_freq_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.png'])
    plot2svg(['WBmod_freq_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.svg'])

%     plot_WBmod_heatmap_AeroCenterFuncs
    plot_WBmod_heatmap_AeroCenterFuncs_plotval
    saveas(gca,['WBmod_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs_1stMomRatio',num2str(FirstMomentRatio_plot),'.fig'])
    saveas(gca,['WBmod_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs_1stMomRatio',num2str(FirstMomentRatio_plot),'.png'])
    plot2svg(['WBmod_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs_1stMomRatio',num2str(FirstMomentRatio_plot),'.svg'])
    
%     plot_WBmodNsteady_AeroCenterFuncs
    plot_WBmodNsteady_AeroCenterFuncs_plotval
    saveas(gca,['WBmodNsteady_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs_1stMomRatio',num2str(FirstMomentRatio_plot),'.fig'])
    saveas(gca,['WBmodNsteady_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs_1stMomRatio',num2str(FirstMomentRatio_plot),'.png'])
    plot2svg(['WBmodNsteady_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs_1stMomRatio',num2str(FirstMomentRatio_plot),'.svg'])
    
    plot_WBmod_strokeMAXcoeff_AeroCenterFuncs
    saveas(gca,['WBmod_strokeMAXcoeff_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.fig'])
    saveas(gca,['WBmod_strokeMAXcoeff_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.png'])
    plot2svg(['WBmod_strokeMAXcoeff_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.svg'])

    plot_WBmod_strokeMINcoeff_AeroCenterFuncs
    saveas(gca,['WBmod_strokeMINcoeff_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.fig'])
    saveas(gca,['WBmod_strokeMINcoeff_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.png'])
    plot2svg(['WBmod_strokeMINcoeff_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.svg'])

    plot_WBmod_pitchMIDDScoeff_AeroCenterFuncs
    saveas(gca,['WBmod_pitchMIDDScoeff_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.fig'])
    saveas(gca,['WBmod_pitchMIDDScoeff_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.png'])
    plot2svg(['WBmod_pitchMIDDScoeff_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.svg'])

    plot_WBmod_pitchMIDUScoeff_AeroCenterFuncs
    saveas(gca,['WBmod_pitchMIDUScoeff_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.fig'])
    saveas(gca,['WBmod_pitchMIDUScoeff_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.png'])
    plot2svg(['WBmod_pitchMIDUScoeff_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.svg'])

    plot_WBmod_devDSMAXcoeff_AeroCenterFuncs
    saveas(gca,['WBmod_devDSMAXcoeff_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.fig'])
    saveas(gca,['WBmod_devDSMAXcoeff_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.png'])
    plot2svg(['WBmod_devDSMAXcoeff_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.svg'])

    plot_WBmod_devDSMINcoeff_AeroCenterFuncs
    saveas(gca,['WBmod_devDSMINcoeff_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.fig'])
    saveas(gca,['WBmod_devDSMINcoeff_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.png'])
    plot2svg(['WBmod_devDSMINcoeff_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.svg'])

    plot_WBmod_devUSMAXcoeff_AeroCenterFuncs
    saveas(gca,['WBmod_devUSMAXcoeff_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.fig'])
    saveas(gca,['WBmod_devUSMAXcoeff_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.png'])
    plot2svg(['WBmod_devUSMAXcoeff_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.svg'])
    
    plot_WBmod_devUSMINcoeff_AeroCenterFuncs
    saveas(gca,['WBmod_devUSMINcoeff_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.fig'])
    saveas(gca,['WBmod_devUSMINcoeff_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.png'])
    plot2svg(['WBmod_devUSMINcoeff_AeroCenterFuncs_',num2str(n_AeroCenterFuncs),'WBs.svg'])
cd ..
end
end


%% save all
if save_on == 1
    save('WBdataset_steadyNclipMods.mat')
end


