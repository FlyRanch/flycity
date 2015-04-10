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

calc_WB_S2S3funcIntact = 1;
% calc_WB_S2S3funcIntact = 0;

calc_WB_S2S3funcClipped = 1;
% calc_WB_S2S3funcClipped = 0;

%% S2S3 Amplitude func settings
% sol = subs(solAi,S2,SecondMomentRatio);
% sol = subs(sol,S3,ThirdMomentRatio);
% S2S3funcIntact = eval(sol);
% 
% sol = subs(solAd,S2,SecondMomentRatio);
% sol = subs(sol,S3,ThirdMomentRatio);
% S2S3funcClipped = eval(sol);
% 
for i = 1: length(SecondMomentRatio)
%     counter = length(SecondMomentRatio)-i
    
    sol = subs(solAi,S2,SecondMomentRatio(i));
    sol = subs(sol,S3,ThirdMomentRatio(i));
    S2S3funcIntact(i,1) = eval(sol);
    
    sol = subs(solAd,S2,SecondMomentRatio(i));
    sol = subs(sol,S3,ThirdMomentRatio(i));
    S2S3funcClipped(i,1) = eval(sol);
end
    
    sol = subs(solAi,S2,1);
    sol = subs(sol,S3,1);
    S2S3funcIntact_NONclipped = eval(sol);
    
    sol = subs(solAd,S2,1);
    sol = subs(sol,S3,1);
    S2S3funcClipped_NONclipped = eval(sol);

% cut-offs
std_factor = 1;
S2S3funcIntact_cut = std_factor*nanstd(S2S3funcIntact);
S2S3funcClipped_cut = std_factor*nanstd(S2S3funcClipped);

%% plot values
% std_factor_plot = 3;
% S2S3funcIntact_plot = std_factor_plot*nanstd(S2S3funcIntact);
% S2S3funcClipped_plot = std_factor_plot*nanstd(S2S3funcClipped);
S2norm_plot = 0.5;
S3norm_plot = 0.4;

sol = subs(solAi,S2,S2norm_plot);
sol = subs(sol,S3,S3norm_plot);
S2S3funcIntact_plot = S2S3funcIntact_NONclipped-eval(sol);

sol = subs(solAd,S2,S2norm_plot);
sol = subs(sol,S3,S3norm_plot);
S2S3funcClipped_plot = eval(sol)-S2S3funcClipped_NONclipped;

figure
plot(S2S3funcIntact,S2S3funcClipped,'o')
hold on
plot(S2S3funcIntact_NONclipped,S2S3funcClipped_NONclipped,'g*')
plot(S2S3funcIntact_NONclipped-S2S3funcIntact_plot,S2S3funcClipped_NONclipped+S2S3funcClipped_plot,'r*')
plot([S2S3funcIntact_NONclipped-S2S3funcIntact_cut S2S3funcIntact_NONclipped-S2S3funcIntact_cut],[S2S3funcClipped_NONclipped+S2S3funcClipped_plot S2S3funcClipped_NONclipped],'-g')
plot([S2S3funcIntact_NONclipped-S2S3funcIntact_plot S2S3funcIntact_NONclipped],[S2S3funcClipped_NONclipped+S2S3funcClipped_cut S2S3funcClipped_NONclipped+S2S3funcClipped_cut],'-g')
legend('all flies','steady non-clipped','plot value','cut-off boundaries')
xlabel('S2S3funcIntact')
ylabel('S2S3funcClipped')

mkdir('clippedfly_steadyWBkin_param_figs')
cd('clippedfly_steadyWBkin_param_figs')

saveas(gca,['S2S3AmpFuncs_realfly_cutoff_n_plot_values.fig'])
saveas(gca,['S2S3AmpFuncs_realfly_cutoff_n_plot_values.png'])
plot2svg(['S2S3AmpFuncs_realfly_cutoff_n_plot_values.svg'])

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

%% CALC CLIPPED WB MODS: S2S3funcIntact for intact wing
if calc_WB_S2S3funcIntact == 1
    
% % calc_wbNORM_wingatt_CLIPPED_S2S3funcIntact_cutoff
% calc_wbNORM_wingatt_CLIPPED_S2S3funcIntact_cutoff_NOcircmean
calc_wbNORM_wingatt_CLIPPED_S2S3AmpFuncIntact_NOcircmean
% calc_wbNORM_wingatt_CLIPPED_S2S3ForceFuncIntact_NOcircmean

if plot_on == 1
mkdir('WBmod_figs_S2S3AmpFuncIntact')
cd('WBmod_figs_S2S3AmpFuncIntact')

xmin = floor(100*min(S2S3funcIntact(:)))/100;
xmax = ceil(100*max(S2S3funcIntact(:)))/100;

%% d1stMom distribution
figure
hist_all = S2S3funcIntact;
hist_all(1:length(S2S3funcIntact_S2S3funcIntact),2) = S2S3funcIntact_S2S3funcIntact;
hist_all(hist_all==0)=nan;
h=hist(hist_all,10);
hist_all(1:max(h(:)),3) = S2S3funcIntact_NONclipped;
hist_all(hist_all==0)=nan;
hist(hist_all,10)
legend('all wingbeats','steady & threshold','steady non-clipped')
% h=hist(S2S3funcIntact_S2S3funcIntact,10);
% hist(S2S3funcIntact_S2S3funcIntact,10)
hold on
title('S2S3funcIntact of clipped and intact wings')
xlabel('S2S3funcIntact')
ylabel('# of wingbeats')
xlim([xmin xmax])
saveas(gca,['hist_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.fig'])
saveas(gca,['hist_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.png'])
plot2svg(['hist_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.svg'])

%% body kin MODs
%     plot_BODYmod_vel_hist_S2S3funcIntact
%     saveas(gca,['WBmod_vel_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.fig'])
%     saveas(gca,['WBmod_vel_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.png'])
%     plot2svg(['WBmod_vel_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.svg'])
% 
%     plot_BODYmod_slip_hist_S2S3funcIntact
%     saveas(gca,['WBmod_slip_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.fig'])
%     saveas(gca,['WBmod_slip_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.png'])
%     plot2svg(['WBmod_slip_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.svg'])
% 
%     plot_BODYmod_pitch_hist_S2S3funcIntact
%     saveas(gca,['WBmod_pitch_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.fig'])
%     saveas(gca,['WBmod_pitch_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.png'])
%     plot2svg(['WBmod_pitch_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.svg'])
% 
%     plot_BODYmod_roll_hist_S2S3funcIntact
%     saveas(gca,['WBmod_roll_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.fig'])
%     saveas(gca,['WBmod_roll_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.png'])
%     plot2svg(['WBmod_roll_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.svg'])
% 
%     plot_BODYmod_Fsp_pitch_hist_S2S3funcIntact
%     saveas(gca,['WBmod_Fsp_pitch_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.fig'])
%     saveas(gca,['WBmod_Fsp_pitch_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.png'])
%     plot2svg(['WBmod_Fsp_pitch_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.svg'])
% 
%     plot_BODYmod_Fsp_roll_hist_S2S3funcIntact
%     saveas(gca,['WBmod_Fsp_roll_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.fig'])
%     saveas(gca,['WBmod_Fsp_roll_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.png'])
%     plot2svg(['WBmod_Fsp_roll_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.svg'])

%% WB kin MODs
    plot_WBmod_freq_hist_S2S3funcIntact
    saveas(gca,['WBmod_freq_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.fig'])
    saveas(gca,['WBmod_freq_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.png'])
    plot2svg(['WBmod_freq_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.svg'])

%     plot_WBmod_heatmap_S2S3funcIntact
    plot_WBmod_heatmap_S2S3funcIntact_NEGplotval
    saveas(gca,['WBmod_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs_S2S3funcIntact',num2str(S2S3funcIntact_plot),'.fig'])
    saveas(gca,['WBmod_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs_S2S3funcIntact',num2str(S2S3funcIntact_plot),'.png'])
    plot2svg(['WBmod_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs_S2S3funcIntact',num2str(S2S3funcIntact_plot),'.svg'])
    
%     plot_WBmodNsteady_S2S3funcIntact
    plot_WBmodNsteady_S2S3funcIntact_NEGplotval
    saveas(gca,['WBmodNsteady_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs_S2S3funcIntact',num2str(S2S3funcIntact_plot),'.fig'])
    saveas(gca,['WBmodNsteady_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs_S2S3funcIntact',num2str(S2S3funcIntact_plot),'.png'])
    plot2svg(['WBmodNsteady_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs_S2S3funcIntact',num2str(S2S3funcIntact_plot),'.svg'])
    
    plot_WBmod_strokeMAXcoeff_S2S3funcIntact
    saveas(gca,['WBmod_strokeMAXcoeff_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.fig'])
    saveas(gca,['WBmod_strokeMAXcoeff_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.png'])
    plot2svg(['WBmod_strokeMAXcoeff_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.svg'])

    plot_WBmod_strokeMINcoeff_S2S3funcIntact
    saveas(gca,['WBmod_strokeMINcoeff_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.fig'])
    saveas(gca,['WBmod_strokeMINcoeff_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.png'])
    plot2svg(['WBmod_strokeMINcoeff_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.svg'])

    plot_WBmod_pitchMIDDScoeff_S2S3funcIntact
    saveas(gca,['WBmod_pitchMIDDScoeff_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.fig'])
    saveas(gca,['WBmod_pitchMIDDScoeff_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.png'])
    plot2svg(['WBmod_pitchMIDDScoeff_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.svg'])

    plot_WBmod_pitchMIDUScoeff_S2S3funcIntact
    saveas(gca,['WBmod_pitchMIDUScoeff_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.fig'])
    saveas(gca,['WBmod_pitchMIDUScoeff_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.png'])
    plot2svg(['WBmod_pitchMIDUScoeff_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.svg'])

    plot_WBmod_devDSMAXcoeff_S2S3funcIntact
    saveas(gca,['WBmod_devDSMAXcoeff_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.fig'])
    saveas(gca,['WBmod_devDSMAXcoeff_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.png'])
    plot2svg(['WBmod_devDSMAXcoeff_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.svg'])

    plot_WBmod_devDSMINcoeff_S2S3funcIntact
    saveas(gca,['WBmod_devDSMINcoeff_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.fig'])
    saveas(gca,['WBmod_devDSMINcoeff_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.png'])
    plot2svg(['WBmod_devDSMINcoeff_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.svg'])

    plot_WBmod_devUSMAXcoeff_S2S3funcIntact
    saveas(gca,['WBmod_devUSMAXcoeff_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.fig'])
    saveas(gca,['WBmod_devUSMAXcoeff_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.png'])
    plot2svg(['WBmod_devUSMAXcoeff_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.svg'])
    
    plot_WBmod_devUSMINcoeff_S2S3funcIntact
    saveas(gca,['WBmod_devUSMINcoeff_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.fig'])
    saveas(gca,['WBmod_devUSMINcoeff_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.png'])
    plot2svg(['WBmod_devUSMINcoeff_S2S3funcIntact_',num2str(n_S2S3funcIntact),'WBs.svg'])
cd ..
end
end

%% CALC CLIPPED WB MODS: S2S3funcClipped for CLIPPED wing
if calc_WB_S2S3funcClipped == 1
    
% % calc_wbNORM_wingatt_CLIPPED_S2S3funcClipped_cutoff
% calc_wbNORM_wingatt_CLIPPED_S2S3funcClipped_cutoff_NOcircmean
calc_wbNORM_wingatt_CLIPPED_S2S3AmpFuncClipped_NOcircmean
% calc_wbNORM_wingatt_CLIPPED_S2S3ForceFuncClipped_NOcircmean

if plot_on == 1
mkdir('WBmod_figs_S2S3AmpFuncClipped')
cd('WBmod_figs_S2S3AmpFuncClipped')

xmin = floor(100*min(S2S3funcClipped(:)))/100;
xmax = ceil(100*max(S2S3funcClipped(:)))/100;

%% d1stMom distribution
figure
hist_all = S2S3funcClipped;
hist_all(1:length(S2S3funcClipped_S2S3funcClipped),2) = S2S3funcClipped_S2S3funcClipped;
hist_all(hist_all==0)=nan;
h=hist(hist_all,10);
hist_all(1:max(h(:)),3) = S2S3funcClipped_NONclipped;
hist_all(hist_all==0)=nan;
hist(hist_all,10)
legend('all wingbeats','steady & threshold','steady non-clipped')
% h=hist(S2S3funcClipped_S2S3funcClipped,10);
% hist(S2S3funcClipped_S2S3funcClipped,10)
hold on
title('S2S3funcClipped of clipped and intact wings')
xlabel('S2S3funcClipped')
ylabel('# of wingbeats')
xlim([xmin xmax])
saveas(gca,['hist_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.fig'])
saveas(gca,['hist_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.png'])
plot2svg(['hist_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.svg'])

%% body kin MODs
%     plot_BODYmod_vel_hist_S2S3funcClipped
%     saveas(gca,['WBmod_vel_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.fig'])
%     saveas(gca,['WBmod_vel_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.png'])
%     plot2svg(['WBmod_vel_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.svg'])
% 
%     plot_BODYmod_slip_hist_S2S3funcClipped
%     saveas(gca,['WBmod_slip_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.fig'])
%     saveas(gca,['WBmod_slip_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.png'])
%     plot2svg(['WBmod_slip_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.svg'])
% 
%     plot_BODYmod_pitch_hist_S2S3funcClipped
%     saveas(gca,['WBmod_pitch_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.fig'])
%     saveas(gca,['WBmod_pitch_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.png'])
%     plot2svg(['WBmod_pitch_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.svg'])
% 
%     plot_BODYmod_roll_hist_S2S3funcClipped
%     saveas(gca,['WBmod_roll_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.fig'])
%     saveas(gca,['WBmod_roll_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.png'])
%     plot2svg(['WBmod_roll_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.svg'])
% 
%     plot_BODYmod_Fsp_pitch_hist_S2S3funcClipped
%     saveas(gca,['WBmod_Fsp_pitch_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.fig'])
%     saveas(gca,['WBmod_Fsp_pitch_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.png'])
%     plot2svg(['WBmod_Fsp_pitch_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.svg'])
% 
%     plot_BODYmod_Fsp_roll_hist_S2S3funcClipped
%     saveas(gca,['WBmod_Fsp_roll_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.fig'])
%     saveas(gca,['WBmod_Fsp_roll_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.png'])
%     plot2svg(['WBmod_Fsp_roll_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.svg'])

%% WB kin MODs
    plot_WBmod_freq_hist_S2S3funcClipped
    saveas(gca,['WBmod_freq_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.fig'])
    saveas(gca,['WBmod_freq_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.png'])
    plot2svg(['WBmod_freq_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.svg'])

%     plot_WBmod_heatmap_S2S3funcClipped
    plot_WBmod_heatmap_S2S3funcClipped_plotval
    saveas(gca,['WBmod_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs_S2S3funcClipped',num2str(S2S3funcClipped_plot),'.fig'])
    saveas(gca,['WBmod_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs_S2S3funcClipped',num2str(S2S3funcClipped_plot),'.png'])
    plot2svg(['WBmod_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs_S2S3funcClipped',num2str(S2S3funcClipped_plot),'.svg'])
    
%     plot_WBmodNsteady_S2S3funcClipped
    plot_WBmodNsteady_S2S3funcClipped_plotval
    saveas(gca,['WBmodNsteady_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs_S2S3funcClipped',num2str(S2S3funcClipped_plot),'.fig'])
    saveas(gca,['WBmodNsteady_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs_S2S3funcClipped',num2str(S2S3funcClipped_plot),'.png'])
    plot2svg(['WBmodNsteady_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs_S2S3funcClipped',num2str(S2S3funcClipped_plot),'.svg'])
    
    plot_WBmod_strokeMAXcoeff_S2S3funcClipped
    saveas(gca,['WBmod_strokeMAXcoeff_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.fig'])
    saveas(gca,['WBmod_strokeMAXcoeff_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.png'])
    plot2svg(['WBmod_strokeMAXcoeff_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.svg'])

    plot_WBmod_strokeMINcoeff_S2S3funcClipped
    saveas(gca,['WBmod_strokeMINcoeff_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.fig'])
    saveas(gca,['WBmod_strokeMINcoeff_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.png'])
    plot2svg(['WBmod_strokeMINcoeff_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.svg'])

    plot_WBmod_pitchMIDDScoeff_S2S3funcClipped
    saveas(gca,['WBmod_pitchMIDDScoeff_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.fig'])
    saveas(gca,['WBmod_pitchMIDDScoeff_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.png'])
    plot2svg(['WBmod_pitchMIDDScoeff_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.svg'])

    plot_WBmod_pitchMIDUScoeff_S2S3funcClipped
    saveas(gca,['WBmod_pitchMIDUScoeff_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.fig'])
    saveas(gca,['WBmod_pitchMIDUScoeff_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.png'])
    plot2svg(['WBmod_pitchMIDUScoeff_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.svg'])

    plot_WBmod_devDSMAXcoeff_S2S3funcClipped
    saveas(gca,['WBmod_devDSMAXcoeff_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.fig'])
    saveas(gca,['WBmod_devDSMAXcoeff_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.png'])
    plot2svg(['WBmod_devDSMAXcoeff_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.svg'])

    plot_WBmod_devDSMINcoeff_S2S3funcClipped
    saveas(gca,['WBmod_devDSMINcoeff_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.fig'])
    saveas(gca,['WBmod_devDSMINcoeff_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.png'])
    plot2svg(['WBmod_devDSMINcoeff_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.svg'])

    plot_WBmod_devUSMAXcoeff_S2S3funcClipped
    saveas(gca,['WBmod_devUSMAXcoeff_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.fig'])
    saveas(gca,['WBmod_devUSMAXcoeff_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.png'])
    plot2svg(['WBmod_devUSMAXcoeff_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.svg'])
    
    plot_WBmod_devUSMINcoeff_S2S3funcClipped
    saveas(gca,['WBmod_devUSMINcoeff_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.fig'])
    saveas(gca,['WBmod_devUSMINcoeff_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.png'])
    plot2svg(['WBmod_devUSMINcoeff_S2S3funcClipped_',num2str(n_S2S3funcClipped),'WBs.svg'])
cd ..
end
end

%% save all
if save_on == 1
    save('WBdataset_steadyNclipMods_S2S3AmpFuncs.mat')
end


