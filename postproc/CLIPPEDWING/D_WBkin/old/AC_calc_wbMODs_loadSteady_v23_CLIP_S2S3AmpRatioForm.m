clc
clear
close all

Eqname=dir('roboflyDB_CutAndIntactWing_EqSolved_AnFnM_vs_S2nS3_clippedFlyWBfreq*')
Eqname=Eqname.name;
load(Eqname)

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

calc_WB_S2S3AmpRatioFunc = 1;
% calc_WB_S2S3AmpRatioFunc = 1;

%% S2S3 Amplitude Ratio func settings
for i = 1: length(SecondMomentRatio)
%     counter = length(SecondMomentRatio)-i
    
    sol = subs(solAdAiRatio,S2,SecondMomentRatio(i));
    sol = subs(sol,S3,ThirdMomentRatio(i));
    S2S3AmpRatioFunc(i,1) = eval(sol);
end
    
    sol = subs(solAdAiRatio,S2,1);
    sol = subs(sol,S3,1);
    S2S3AmpRatioFunc_NONclipped = eval(sol);

% cut-offs
std_factor = 1.5;
S2S3AmpRatioFunc_cut = std_factor*nanstd(S2S3AmpRatioFunc);

%% plot values
% std_factor_plot = 3;
% S2S3AmpRatioFunc_plot = std_factor_plot*nanstd(S2S3AmpRatioFunc);
S2norm_plot = 0.5;
S3norm_plot = 0.5;

sol = subs(solAdAiRatio,S2,S2norm_plot);
sol = subs(sol,S3,S3norm_plot);
S2S3AmpRatioFunc_plot = eval(sol)-S2S3AmpRatioFunc_NONclipped;

figure
hist_all = S2S3AmpRatioFunc;
h=hist(hist_all,20);
hist(hist_all,20)
hold on
plot([S2S3AmpRatioFunc_NONclipped S2S3AmpRatioFunc_NONclipped],[0 max(h)])
plot([S2S3AmpRatioFunc_NONclipped+S2S3AmpRatioFunc_cut S2S3AmpRatioFunc_NONclipped+S2S3AmpRatioFunc_cut],[0 max(h)],'r')
plot([S2S3AmpRatioFunc_NONclipped+S2S3AmpRatioFunc_plot S2S3AmpRatioFunc_NONclipped+S2S3AmpRatioFunc_plot],[0 max(h)],'g')
legend('all wingbeats','steady','steady & threshold','plot value')
xlabel('S2S3AmpRatioFunc')
ylabel('# of wingbeats')

mkdir('clippedfly_steadyWBkin_param_figs')
cd('clippedfly_steadyWBkin_param_figs')

saveas(gca,['S2S3AmpRatioFunc_realfly_cutoff_n_plot_values.fig'])
saveas(gca,['S2S3AmpRatioFunc_realfly_cutoff_n_plot_values.png'])
plot2svg(['S2S3AmpRatioFunc_realfly_cutoff_n_plot_values.svg'])

cd ..

%% filter vars
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
% close all

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

%% CALC CLIPPED & INTACT WB MODS: S2S3AmpRatioFunc
if calc_WB_S2S3AmpRatioFunc == 1
    
    calc_wbNORM_wingatt_CLIPPED_S2S3AmpRatioFunc_NOcircmean

if plot_on == 1
mkdir('WBmod_figs_S2S3AmpRatioFunc')
cd('WBmod_figs_S2S3AmpRatioFunc')

% xmin = floor(100*min(S2S3AmpRatioFunc(:)))/100;
% xmax = ceil(100*max(S2S3AmpRatioFunc(:)))/100;

xmin = floor(100*min(S2S3AmpRatioFunc_NONclipped))/100;
xmax =  ceil(100*max(S2S3AmpRatioFunc_NONclipped+S2S3AmpRatioFunc_plot))/100;

%% body kin MODs
%     plot_BODYmod_vel_hist_S2S3AmpRatioFunc
%     saveas(gca,['WBmod_vel_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.fig'])
%     saveas(gca,['WBmod_vel_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.png'])
%     plot2svg(['WBmod_vel_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.svg'])
% 
%     plot_BODYmod_slip_hist_S2S3AmpRatioFunc
%     saveas(gca,['WBmod_slip_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.fig'])
%     saveas(gca,['WBmod_slip_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.png'])
%     plot2svg(['WBmod_slip_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.svg'])
% 
%     plot_BODYmod_pitch_hist_S2S3AmpRatioFunc
%     saveas(gca,['WBmod_pitch_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.fig'])
%     saveas(gca,['WBmod_pitch_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.png'])
%     plot2svg(['WBmod_pitch_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.svg'])
% 
%     plot_BODYmod_roll_hist_S2S3AmpRatioFunc
%     saveas(gca,['WBmod_roll_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.fig'])
%     saveas(gca,['WBmod_roll_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.png'])
%     plot2svg(['WBmod_roll_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.svg'])
% 
%     plot_BODYmod_Fsp_pitch_hist_S2S3AmpRatioFunc
%     saveas(gca,['WBmod_Fsp_pitch_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.fig'])
%     saveas(gca,['WBmod_Fsp_pitch_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.png'])
%     plot2svg(['WBmod_Fsp_pitch_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.svg'])
% 
%     plot_BODYmod_Fsp_roll_hist_S2S3AmpRatioFunc
%     saveas(gca,['WBmod_Fsp_roll_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.fig'])
%     saveas(gca,['WBmod_Fsp_roll_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.png'])
%     plot2svg(['WBmod_Fsp_roll_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.svg'])

%% WB kin MODs
    plot_WBmod_freq_hist_S2S3funcAmpRatioFunc
    saveas(gca,['WBmod_freq_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.fig'])
    saveas(gca,['WBmod_freq_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.png'])
    plot2svg(['WBmod_freq_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.svg'])

    plot_WBmod_heatmap_S2S3AmpRatioFunc
    saveas(gca,['WBmod_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs_S2S3AmpRatioFunc',num2str(S2S3AmpRatioFunc_plot),'.fig'])
    saveas(gca,['WBmod_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs_S2S3AmpRatioFunc',num2str(S2S3AmpRatioFunc_plot),'.png'])
    plot2svg(['WBmod_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs_S2S3AmpRatioFunc',num2str(S2S3AmpRatioFunc_plot),'.svg'])
    
    plot_WBmodNsteady_S2S3AmpRatioFunc
    saveas(gca,['WBmodNsteady_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs_S2S3AmpRatioFunc',num2str(S2S3AmpRatioFunc_plot),'.fig'])
    saveas(gca,['WBmodNsteady_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs_S2S3AmpRatioFunc',num2str(S2S3AmpRatioFunc_plot),'.png'])
    plot2svg(['WBmodNsteady_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs_S2S3AmpRatioFunc',num2str(S2S3AmpRatioFunc_plot),'.svg'])
    
    plot_WBmod_strokeMAXcoeff_S2S3AmpRatioFunc
    saveas(gca,['WBmod_strokeMAXcoeff_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.fig'])
    saveas(gca,['WBmod_strokeMAXcoeff_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.png'])
    plot2svg(['WBmod_strokeMAXcoeff_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.svg'])

    plot_WBmod_strokeMINcoeff_S2S3AmpRatioFunc
    saveas(gca,['WBmod_strokeMINcoeff_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.fig'])
    saveas(gca,['WBmod_strokeMINcoeff_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.png'])
    plot2svg(['WBmod_strokeMINcoeff_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.svg'])

    plot_WBmod_pitchMIDDScoeff_S2S3AmpRatioFunc
    saveas(gca,['WBmod_pitchMIDDScoeff_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.fig'])
    saveas(gca,['WBmod_pitchMIDDScoeff_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.png'])
    plot2svg(['WBmod_pitchMIDDScoeff_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.svg'])

    plot_WBmod_pitchMIDUScoeff_S2S3AmpRatioFunc
    saveas(gca,['WBmod_pitchMIDUScoeff_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.fig'])
    saveas(gca,['WBmod_pitchMIDUScoeff_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.png'])
    plot2svg(['WBmod_pitchMIDUScoeff_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.svg'])

    plot_WBmod_devDSMAXcoeff_S2S3AmpRatioFunc
    saveas(gca,['WBmod_devDSMAXcoeff_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.fig'])
    saveas(gca,['WBmod_devDSMAXcoeff_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.png'])
    plot2svg(['WBmod_devDSMAXcoeff_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.svg'])

    plot_WBmod_devDSMINcoeff_S2S3AmpRatioFunc
    saveas(gca,['WBmod_devDSMINcoeff_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.fig'])
    saveas(gca,['WBmod_devDSMINcoeff_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.png'])
    plot2svg(['WBmod_devDSMINcoeff_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.svg'])

    plot_WBmod_devUSMAXcoeff_S2S3AmpRatioFunc
    saveas(gca,['WBmod_devUSMAXcoeff_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.fig'])
    saveas(gca,['WBmod_devUSMAXcoeff_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.png'])
    plot2svg(['WBmod_devUSMAXcoeff_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.svg'])
    
    plot_WBmod_devUSMINcoeff_S2S3AmpRatioFunc
    saveas(gca,['WBmod_devUSMINcoeff_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.fig'])
    saveas(gca,['WBmod_devUSMINcoeff_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.png'])
    plot2svg(['WBmod_devUSMINcoeff_S2S3AmpRatioFunc_',num2str(n_S2S3AmpRatioFunc),'WBs.svg'])
cd ..
end
end

%% save all
if save_on == 1
    save('WBdataset_steadyNclipMods_S2S3AmpRatioFunc.mat')
end


