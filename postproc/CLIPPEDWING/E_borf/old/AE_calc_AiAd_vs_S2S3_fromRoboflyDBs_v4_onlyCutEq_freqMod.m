
clear
% clc
close all

%% steady wb freq
% const
A_min = .5;
A_max = 1.5;
S2_min = .5;
S2_max = 1;
S3_min = .4;
S3_max = 1;

F_min = -1;
F_max = 0;
M_min = -1;
M_max = 1;

Ftot_min = -1.5;
Ftot_max = -0.5;

Mtot_min = -1;
Mtot_max =  1;

% colormap: blue to white to red
cmap_surf=cbrewer('div','RdBu',100);
cmap_surf = flipud(cmap_surf);

% black at edges
cmap_surf(1,:) = [0 0 0];
cmap_surf(end,:) = [0 0 0];

plot_on = 0;
plot_on = 1;

%% load DB's
% fly data
var_file = dir('flyVar*');
load(var_file.name)
% steady Wb F&M&CoM data
load('steadyWB_FnMnCoM_data.mat')
% roboflyDB of Cut wing FnM vs Amp&S2&S3
load('roboflyDB_CutWing_FnM_vs_StrokeAmplitude_vs_S2nS3_INCcaliCF.mat')

%% fit param's
% cut wing vertical force vs Amp (parabolic fit) & S2 (linear fit)
Fz_Amp_S2_SurfFit_coeffs = coeffvalues(Fz_Amp_S2_SurfFit);

% cut wing roll torque vs Amp (parabolic fit) & S3 (linear fit)
Mx_Amp_S3_SurfFit_coeffs = coeffvalues(Mx_MinSteady_Amp_S3_SurfFit);

%% calc Amp's from vertical force and roll torque balance
calc_AiAd_FiFd_MiMd_AT_Fz0nMx0_onlyCutWingEq

%% plot Ai&Ad Fi&Fd Mi&Md @ weight support & zero roll torque
mkdir('figures_cutWing_robofly')
cd('figures_cutWing_robofly')

figure(1)
plot_AiAd_FiFd_MiMd_AT_weightSupport_n_RollTorqueEquilibrium
saveas(gcf,'AnFnM_vs_S2nS3_CutnIntactWing_steadyWBfreq.fig')
saveas(gcf,'AnFnM_vs_S2nS3_CutnIntactWing_steadyWBfreq.png')
% plot2svg('AnFnM_vs_S2nS3_CutnIntactWing_steadyWBfreq.svg')

figure(2)
plot_FtotalFiFd_MtotalMiMd_AT_weightSupport_n_RollEquilibrium
saveas(gcf,'FnMtotal_vs_S2nS3_steadyWBfreq.fig')
saveas(gcf,'FnMtotal_vs_S2nS3_steadyWBfreq.png')
% plot2svg('FnMtotal_vs_S2nS3_steadyWBfreq.svg')

cd ..

%% save data
save('roboflyDB_CutAndIntactWing_EqSolved_AnFnM_vs_S2nS3_steadyWBfreq.mat')







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% clipped fly wb freq

clear
clc
close all

% const
A_min = .5;
A_max = 1.5;
S2_min = .5;
S2_max = 1;
S3_min = .4;
S3_max = 1;

F_min = -1;
F_max = 0;
M_min = -1;
M_max = 1;

Ftot_min = -1.5;
Ftot_max = -0.5;

Mtot_min = -1;
Mtot_max =  1;

% colormap: blue to white to red
cmap_surf=cbrewer('div','RdBu',100);
cmap_surf = flipud(cmap_surf);

% black at edges
cmap_surf(1,:) = [0 0 0];
cmap_surf(end,:) = [0 0 0];

plot_on = 0;
plot_on = 1;

%% load DB's
% fly data
var_file = dir('flyVar*');
load(var_file.name)
% steady Wb F&M&CoM data
load('steadyWB_FnMnCoM_data.mat')
% roboflyDB of Non-cut wing (AmpReduce)
load('roboflyDB_NONcutWing_FnM_vs_ReducedAmpStrokeRatio_INCcaliCF_INCclippedWingFreq.mat')
% roboflyDB of Cut wing FnM vs Amp vs S2&S3
load('roboflyDB_CutWing_FnM_vs_StrokeAmplitude_vs_S2nS3_INCcaliCF.mat')


%% fit param's
% intact wing vertical force vs Amp (linear fit)
pFi1 = Fz_Amp_fit_freqMod(1);
pFi0 = Fz_Amp_fit_freqMod(2);

% intact wing roll torque vs Amp (linear fit)
pMi1 = MxMinSteady_Amp_fit_freqMod(1);
pMi0 = MxMinSteady_Amp_fit_freqMod(2);

% cut wing vertical force vs Amp (parabolic fit) & S2 (linear fit)
Fz_Amp_S2_SurfFit_coeffs = coeffvalues(Fz_Amp_S2_SurfFit_freqMod);

pFd00 = Fz_Amp_S2_SurfFit_coeffs(1);
pFd10 = Fz_Amp_S2_SurfFit_coeffs(2);
pFd01 = Fz_Amp_S2_SurfFit_coeffs(3);
pFd11 = Fz_Amp_S2_SurfFit_coeffs(4);
pFd02 = Fz_Amp_S2_SurfFit_coeffs(5);

% cut wing roll torque vs Amp (parabolic fit) & S3 (linear fit)
Mx_Amp_S3_SurfFit_coeffs = coeffvalues(Mx_MinSteady_Amp_S3_SurfFit_freqMod);

pMd00 = Mx_Amp_S3_SurfFit_coeffs(1);
pMd10 = Mx_Amp_S3_SurfFit_coeffs(2);
pMd01 = Mx_Amp_S3_SurfFit_coeffs(3);
pMd11 = Mx_Amp_S3_SurfFit_coeffs(4);
pMd02 = Mx_Amp_S3_SurfFit_coeffs(5);

% calc Amp's from vertical force and roll torque balance
calc_AiAd_FiFd_MiMd_AT_weightSupport_n_RollTorqueEquilibrium

% plot Ai&Ad Fi&Fd Mi&Md @ weight support & zero roll torque
mkdir('figures_cutWing_robofly')
cd('figures_cutWing_robofly')

figure(1)
plot_AiAd_FiFd_MiMd_AT_weightSupport_n_RollTorqueEquilibrium
saveas(gcf,'AnFnM_vs_S2nS3_CutnIntactWing_clippedFlyWBfreq.fig')
saveas(gcf,'AnFnM_vs_S2nS3_CutnIntactWing_clippedFlyWBfreq.png')
% plot2svg('AnFnM_vs_S2nS3_CutnIntactWing_clippedFlyWBfreq.svg')

figure(2)
plot_FtotalFiFd_MtotalMiMd_AT_weightSupport_n_RollEquilibrium
saveas(gcf,'FnMtotal_vs_S2nS3_clippedFlyWBfreq.fig')
saveas(gcf,'FnMtotal_vs_S2nS3_clippedFlyWBfreq.png')
% plot2svg('FnMtotal_vs_S2nS3_clippedFlyWBfreq.svg')

cd ..

%% save data
save('roboflyDB_CutAndIntactWing_EqSolved_AnFnM_vs_S2nS3_clippedFlyWBfreq.mat')










