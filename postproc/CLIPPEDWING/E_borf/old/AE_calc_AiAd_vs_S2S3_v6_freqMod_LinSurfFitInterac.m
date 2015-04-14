
%% steady wb freq

clear
clc
close all
warning off

% maximum possible wingbeat
load('WBdataset_steady_1603WBs.mat')
stroke_amp_steady = max(stroke_wb_steady_bins_meanCIstd(:,1))-min(stroke_wb_steady_bins_meanCIstd(:,1));
A_ratio_max = 180/stroke_amp_steady;

% const
A_min = 0.5;
A_max = 1.5;
AiAdRatio_min = 0.5;
AiAdRatio_max = 1.5;

S2_min = 0;
S2_max = 1;
S3_min = 0;
S3_max = 1;

F_min = -.75;
F_max = -.25;
M_min = -.5;
M_max = .5;

% colormap: blue to white to red
cmap_surf=cbrewer('div','RdBu',100);
cmap_surf = flipud(cmap_surf);

% black at Amp > 180deg
cmap_Aratio = cmap_surf;
n_Amax = round(99*(A_ratio_max-A_max+1)+1);
cmap_Aratio(n_Amax:end,:)=0;

% black at edges of cmap
cmap_edge = cmap_surf;
cmap_edge(1,:) = [0 0 0];
cmap_edge(end,:) = [0 0 0];

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
% roboflyDB of Cut wing FnM vs Amp&S2&S3 Linear Surface Fit INC interaction
load('roboflyDB_CutWing_FnM_vs_StrokeAmplitude_vs_S2nS3_INCcaliCF_LinSurfFitInterac.mat')

%% fit param's
% cut wing vertical force vs Amp (parabolic fit) & S2 (linear fit)
Fz_Amp_S2_SurfFit_coeffs = coeffvalues(Fz_Amp_S2_SurfFit);
Fz_Amp_fit_coeffs = Fz_Amp_fit;

% cut wing roll torque vs Amp (parabolic fit) & S3 (linear fit)
Mx_Amp_S3_SurfFit_coeffs = coeffvalues(Mx_MinSteady_Amp_S3_SurfFit);
MxMinSteady_Amp_fit_coeffs = MxMinSteady_Amp_fit;

%% calc Amp's from vertical force and roll torque balance
calc_AiAd_FiFd_MiMd_AT_Fz0nMx0

%% plot Ai&Ad Fi&Fd Mi&Md @ weight support & zero roll torque
figure(1)
plot_AiAd_AT_Fz0nMx0_LinSurfFitInterac_steadyWBfreq

figure(2)
plot_FtotalFiFd_MtotalMiMd_AT_AT_Fz0nMx0_LinSurfFitInterac

%% save data
save('roboflyDB_CutAndIntactWing_EqSolved_AnFnM_vs_S2nS3_steadyWBfreq.mat')







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% clipped fly wb freq

clear
% clc
% close all

% maximum possible wingbeat
load('WBdataset_steady_1603WBs.mat')
stroke_amp_steady = max(stroke_wb_steady_bins_meanCIstd(:,1))-min(stroke_wb_steady_bins_meanCIstd(:,1));
A_ratio_max = 180/stroke_amp_steady;

% const
A_min = 0.5;
A_max = 1.5;
AiAdRatio_min = 0.5;
AiAdRatio_max = 1.5;

S2_min = 0;
S2_max = 1;
S3_min = 0;
S3_max = 1;

F_min = -.75;
F_max = -.25;
M_min = -.5;
M_max = .5;

% colormap: blue to white to red
cmap_surf=cbrewer('div','RdBu',100);
cmap_surf = flipud(cmap_surf);

% black at Amp > 180deg
cmap_Aratio = cmap_surf;
n_Amax = round(99*(A_ratio_max-A_max+1)+1);
cmap_Aratio(n_Amax:end,:)=0;

% black at edges of cmap
cmap_edge = cmap_surf;
cmap_edge(1,:) = [0 0 0];
cmap_edge(end,:) = [0 0 0];

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
% roboflyDB of Cut wing FnM vs Amp&S2&S3 Linear Surface Fit INC interaction
load('roboflyDB_CutWing_FnM_vs_StrokeAmplitude_vs_S2nS3_INCcaliCF_LinSurfFitInterac.mat')

%% fit param's
% cut wing vertical force vs Amp (parabolic fit) & S2 (linear fit)
Fz_Amp_S2_SurfFit_coeffs = coeffvalues(Fz_Amp_S2_SurfFit_freqMod);
Fz_Amp_fit_coeffs = Fz_Amp_fit_freqMod;

% cut wing roll torque vs Amp (parabolic fit) & S3 (linear fit)
Mx_Amp_S3_SurfFit_coeffs = coeffvalues(Mx_MinSteady_Amp_S3_SurfFit_freqMod);
MxMinSteady_Amp_fit_coeffs = MxMinSteady_Amp_fit_freqMod;

%% calc Amp's from vertical force and roll torque balance
calc_AiAd_FiFd_MiMd_AT_Fz0nMx0

%% plot Ai&Ad Fi&Fd Mi&Md @ weight support & zero roll torque
figure(1)
plot_AiAd_AT_Fz0nMx0_LinSurfFitInterac_clippedFlyWBfreq

figure(3)
plot_FtotalFiFd_MtotalMiMd_AT_AT_Fz0nMx0_LinSurfFitInterac

%% save data
save('roboflyDB_CutAndIntactWing_EqSolved_AnFnM_vs_S2nS3_clippedFlyWBfreq.mat')

%% add real fly data
load('WBdataset_ClipNintact_wingbeat_kin.mat')
figure(1)
subplot(2,3,1)
hold on
for i = 1:length(Astroke_ratio_clip_mean)
    color_nr = round(98*(Astroke_ratio_clip_mean(i)-A_min)+1);
    if color_nr<1
        color_nr=1
    end
    plot3(S2_ratio_mean(i),S3_ratio_mean(i),2,'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',15)
end

subplot(2,3,2)
hold on
for i = 1:length(Astroke_ratio_intact_mean)
    color_nr = round(98*(Astroke_ratio_intact_mean(i)-A_min)+1);
    if color_nr<1
        color_nr=1
    end
    plot3(S2_ratio_mean(i),S3_ratio_mean(i),2,'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',10)
end

subplot(2,3,3)
hold on
for i = 1:length(Astroke_ratio_intact_clip_mean)
    color_nr = round(98*(Astroke_ratio_intact_clip_mean(i)-A_min)+1);
    if color_nr<1
        color_nr=1
    end
    plot3(S2_ratio_mean(i),S3_ratio_mean(i),2,'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',10)
end

subplot(2,3,4)
hold on
for i = 1:length(Astroke_ratio_clip_mean)
    color_nr = round(98*(Astroke_ratio_clip_mean(i)-A_min)+1);
    if color_nr<1
        color_nr=1
    end
    plot3(S2_ratio_mean(i),S3_ratio_mean(i),2,'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',8)
end

subplot(2,3,5)
hold on
for i = 1:length(Astroke_ratio_intact_mean)
    color_nr = round(98*(Astroke_ratio_intact_mean(i)-A_min)+1);
    if color_nr<1
        color_nr=1
    end
    plot3(S2_ratio_mean(i),S3_ratio_mean(i),2,'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',8)
end

subplot(2,3,6)
hold on
for i = 1:length(Astroke_ratio_intact_clip_mean)
    color_nr = round(98*(Astroke_ratio_intact_clip_mean(i)-A_min)+1);
    if color_nr<1
        color_nr=1
    end
    plot3(S2_ratio_mean(i),S3_ratio_mean(i),2,'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',8)
end

%% save plots
mkdir('figures_cutWing_robofly')
cd('figures_cutWing_robofly')

figure(1)
saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_clippedFlyWBfreq_n_steadyWBfreq_LinSurfFitInterac.fig')
saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_clippedFlyWBfreq_n_steadyWBfreq_LinSurfFitInterac.png')
saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_clippedFlyWBfreq_n_steadyWBfreq_LinSurfFitInterac.svg')
% plot2svg('AinAd_vs_S2nS3_CutnIntactWing_clippedFlyWBfreq_n_steadyWBfreq_LinSurfFitInterac.svg')

figure(2)
saveas(gcf,'FnMtotal_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.fig')
saveas(gcf,'FnMtotal_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.png')
saveas(gcf,'FnMtotal_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.svg')
% plot2svg('FnMtotal_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.svg')

figure(3)
saveas(gcf,'FnMtotal_vs_S2nS3_clippedFlyWBfreq_LinSurfFitInterac.fig')
saveas(gcf,'FnMtotal_vs_S2nS3_clippedFlyWBfreq_LinSurfFitInterac.png')
saveas(gcf,'FnMtotal_vs_S2nS3_clippedFlyWBfreq_LinSurfFitInterac.svg')
% plot2svg('FnMtotal_vs_S2nS3_clippedFlyWBfreq_LinSurfFitInterac.svg')

cd ..




