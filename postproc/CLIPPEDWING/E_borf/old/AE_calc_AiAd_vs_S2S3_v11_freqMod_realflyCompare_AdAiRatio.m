
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
AdAiRatio_min = 1;
AdAiRatio_max = 2;

S2_min = 0;
S2_max = 1;
S3_min = 0;
S3_max = 1;

F_min = -.75;
F_max = -.25;
M_min = -.5;
M_max = .5;

Ftot_min = -.75;
Ftot_max = -.25;
Mtot_min = -.5;
Mtot_max = .5;

Fd_min = -.75;
Fd_max = -.25;
Md_min = -.5;
Md_max = .5;

Fi_min = -.75;
Fi_max = -.25;
Mi_min = -.5;
Mi_max = .5;

% colormap: blue to white to red
cmap_surf=cbrewer('div','RdBu',100);
cmap_surf = flipud(cmap_surf);

% black at Amp > 180deg
cmap_Aratio = cmap_surf;
n_Amax = round(99/(A_max-A_min)*(A_ratio_max-A_min)+1);
if n_Amax <= length(cmap_Aratio)
    cmap_Aratio(n_Amax:end,:)=0;
end

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
calc_AdAi_FiFd_MiMd_AT_Fz0nMx0_AT_WBfreq
% calc_AdAi_FiFd_MiMd_AT_Fz0nMx0
% % calc_AiAd_FiFd_MiMd_AT_Fz0nMx0

%% plot Ai&Ad Fi&Fd Mi&Md @ weight support & zero roll torque
% plot_AiAd_AT_Fz0nMx0_LinSurfFitInterac_steadyWBfreq
plot_AdAi_AT_Fz0nMx0_LinSurfFitInterac_steadyWBfreq

figure(3)
plot_FtotalFiFd_MtotalMiMd_AT_Fz0nMx0_LinSurfFitInterac

figure(9)
plot_Fd_Md_AT_Fz0nMx0_LinSurfFitInterac

%% save data
save('roboflyDB_CutAndIntactWing_EqSolved_AnFnM_vs_S2nS3_steadyWBfreq.mat')

%% comparison with real fly data
load('WBdataset_ClipNintact_wingbeat_kin_S2S3ForceFunc.mat')

figure(2)
% subplot(1,2,2)
% title('stroke amplitude ratios @ non-clipped wingbeat freq')
% axis off

subplot(1,2,1)
hold on
for i = 1:length(Astroke_ratio_clip_mean)
    color_nr = round(99/(A_max-A_min)*(Astroke_ratio_clip_mean(i)-A_min)+1);
    if color_nr<1
        color_nr=1
    end
    plot(S2_ratio_mean(i),S3_ratio_mean(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',8)
end

axis equal
axis tight
axis([S2_min,S2_max,S3_min,S3_max])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
xlabel('S_2')
ylabel('S_3')

subplot(1,2,2)
hold on
for i = 1:length(Astroke_ratio_intact_mean)
    color_nr = round(99/(A_max-A_min)*(Astroke_ratio_intact_mean(i)-A_min)+1);
    if color_nr<1
        color_nr=1
    end
    plot(S2_ratio_mean(i),S3_ratio_mean(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',8)
end

axis equal
axis tight
axis([S2_min,S2_max,S3_min,S3_max])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
xlabel('S_2')
ylabel('S_3')

figure(8)
subplot(1,2,1)
hold on
for i = 1:length(Astroke_ratio_clip_intact_mean)
    color_nr = round(99/(AdAiRatio_max-AdAiRatio_min)*(Astroke_ratio_clip_intact_mean(i)-AdAiRatio_min)+1);
    if color_nr<1
        color_nr=1
    end
    plot(S2_ratio_mean(i),S3_ratio_mean(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',8)
end
axis equal
axis tight
axis([S2_min,S2_max,S3_min,S3_max])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
xlabel('S_2')
ylabel('S_3')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% clipped fly wb freq

clear
% clc
% close all

% maximum possible wingbeat
load('WBdataset_steady_1603WBs.mat')
stroke_amp_steady = max(stroke_wb_steady_bins_meanCIstd(:,1))-min(stroke_wb_steady_bins_meanCIstd(:,1));

A_ratio_max = 180/stroke_amp_steady;
% A_ratio0 = .72;

% const
AdAiRatio_min = 1;
AdAiRatio_max = 2;

S2_min = 0;
S2_max = 1;
S3_min = 0;
S3_max = 1;

Ftot_min = -2;
Ftot_max = 0;
Mtot_min = -.5;
Mtot_max = .5;

Fd_min = -.75;
Fd_max = -.25;
Md_min = -.5;
Md_max = .5;

Fi_min = -.75;
Fi_max = -.25;
Mi_min = -.5;
Mi_max = .5;

% colormap: blue to white to red
cmap_surf=cbrewer('div','RdBu',100);
cmap_surf = flipud(cmap_surf);

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
calc_AdAi_FiFd_MiMd_AT_Fz0nMx0_AT_WBfreq
% calc_AdAi_FiFd_MiMd_AT_Fz0nMx0
% % calc_AiAd_FiFd_MiMd_AT_Fz0nMx0

% colorbar centered at uncut equilibrium
A_min = 0;
A_max = 2*Aratio_intact;

A_max = A_ratio_max;
A_min = A_ratio_max - 2*(A_ratio_max-Aratio_intact);

% cmap_Aratio: centered at uncut equilibrium & black at Amp > 180deg
cmap_Aratio = cmap_surf;
n_Amax = round(99/(A_max-A_min)*(A_ratio_max-A_min)+1);
if n_Amax <= length(cmap_Aratio)
    cmap_Aratio(n_Amax:end,:)=0;
end

%% plot Ai&Ad Fi&Fd Mi&Md @ weight support & zero roll torque
% figure(4)
% plot_AiAd_AT_Fz0nMx0_LinSurfFitInterac_clippedFlyWBfreq
plot_AdAi_AT_Fz0nMx0_LinSurfFitInterac_clippedFlyWBfreq

figure(6)
plot_FtotalFiFd_MtotalMiMd_AT_AT_Fz0nMx0_LinSurfFitInterac

figure(10)
plot_Fd_Md_AT_Fz0nMx0_LinSurfFitInterac

%% save data
save('roboflyDB_CutAndIntactWing_EqSolved_AnFnM_vs_S2nS3_clippedFlyWBfreq.mat')

%% add real fly data
load('WBdataset_ClipNintact_wingbeat_kin_S2S3ForceFunc.mat')

figure(5)
% subplot(1,2,2)
% title('stroke amplitude ratios @ clipped fly wingbeat freq')
% axis off

subplot(1,2,1)
hold on
for i = 1:length(Astroke_ratio_clip_mean)
    color_nr = round(99/(A_max-A_min)*(Astroke_ratio_clip_mean(i)-A_min)+1);

    if color_nr<1
        color_nr=1
    end
    plot(S2_ratio_mean(i),S3_ratio_mean(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',8)
end
axis equal
axis tight
axis([S2_min,S2_max,S3_min,S3_max])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
xlabel('S_2')
ylabel('S_3')

subplot(1,2,2)
hold on
for i = 1:length(Astroke_ratio_intact_mean)
    color_nr = round(99/(A_max-A_min)*(Astroke_ratio_intact_mean(i)-A_min)+1);
    if color_nr<1
        color_nr=1
    end
    plot(S2_ratio_mean(i),S3_ratio_mean(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',8)
end
axis equal
axis tight
axis([S2_min,S2_max,S3_min,S3_max])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
xlabel('S_2')
ylabel('S_3')

figure(8)
subplot(1,2,2)
hold on
for i = 1:length(Astroke_ratio_clip_intact_mean)
    color_nr = round(99/(AdAiRatio_max-AdAiRatio_min)*(Astroke_ratio_clip_intact_mean(i)-AdAiRatio_min)+1);
    if color_nr<1
        color_nr=1
    end
    plot(S2_ratio_mean(i),S3_ratio_mean(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',8)
end
axis equal
axis tight
axis([S2_min,S2_max,S3_min,S3_max])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
xlabel('S_2')
ylabel('S_3')

%% save plots
mkdir('figures_cutWing_robofly')
cd('figures_cutWing_robofly')

figure(1)
saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_steadyWBfreq_LinSurfFitInterac.fig')
saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_steadyWBfreq_LinSurfFitInterac.png')
% % saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_steadyWBfreq_LinSurfFitInterac.svg')
% plot2svg('AinAd_vs_S2nS3_CutnIntactWing_steadyWBfreq_LinSurfFitInterac.svg')

figure(2)
saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_steadyWBfreq_realFlyData.fig')
saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_steadyWBfreq_realFlyData.png')
% % saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_steadyWBfreq_realFlyData.svg')
plot2svg('AinAd_vs_S2nS3_CutnIntactWing_steadyWBfreq_realFlyData.svg')

figure(3)
saveas(gcf,'FnMtotal_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.fig')
saveas(gcf,'FnMtotal_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.png')
% % saveas(gcf,'FnMtotal_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.svg')
% plot2svg('FnMtotal_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.svg')

figure(9)
saveas(gcf,'FnMd_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.fig')
saveas(gcf,'FnMd_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.png')
% % saveas(gcf,'FnMd_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.svg')
% plot2svg('FnMd_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.svg')





figure(4)
saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_cutWBfreq_LinSurfFitInterac.fig')
saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_cutWBfreq_LinSurfFitInterac.png')
% % saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_cutWBfreq_LinSurfFitInterac.svg')
% plot2svg('AinAd_vs_S2nS3_CutnIntactWing_cutWBfreq_LinSurfFitInterac.svg')

figure(5)
saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_cutWBfreq_realFlyData.fig')
saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_cutWBfreq_realFlyData.png')
% % saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_cutWBfreq_realFlyData.svg')
plot2svg('AinAd_vs_S2nS3_CutnIntactWing_cutWBfreq_realFlyData.svg')

figure(6)
saveas(gcf,'FnMtotal_vs_S2nS3_cutWBfreq_LinSurfFitInterac.fig')
saveas(gcf,'FnMtotal_vs_S2nS3_cutWBfreq_LinSurfFitInterac.png')
% % saveas(gcf,'FnMtotal_vs_S2nS3_cutWBfreq_LinSurfFitInterac.svg')
% plot2svg('FnMtotal_vs_S2nS3_cutWBfreq_LinSurfFitInterac.svg')

figure(10)
saveas(gcf,'FnMd_vs_S2nS3_cutWBfreq_LinSurfFitInterac.fig')
saveas(gcf,'FnMd_vs_S2nS3_cutWBfreq_LinSurfFitInterac.png')
% % saveas(gcf,'FnMd_vs_S2nS3_cutWBfreq_LinSurfFitInterac.svg')
% plot2svg('FnMd_vs_S2nS3_cutWBfreq_LinSurfFitInterac.svg')





figure(7)
saveas(gcf,'AdAiRatio_vs_S2nS3_CutnIntactWing_steadyNcutWBfreq_LinSurfFitInterac.fig')
saveas(gcf,'AdAiRatio_vs_S2nS3_CutnIntactWing_steadyNcutWBfreq_LinSurfFitInterac.png')
% % saveas(gcf,'AdAiRatio_vs_S2nS3_CutnIntactWing_steadyNcutWBfreq_LinSurfFitInterac.svg')
% plot2svg('AdAiRatio_vs_S2nS3_CutnIntactWing_steadyNcutWBfreq_LinSurfFitInterac.svg')

figure(8)
saveas(gcf,'AdAiRatio_vs_S2nS3_CutnIntactWing_steadyNcutWBfreq_realFlyData.fig')
saveas(gcf,'AdAiRatio_vs_S2nS3_CutnIntactWing_steadyNcutWBfreq_realFlyData.png')
% % saveas(gcf,'AdAiRatio_vs_S2nS3_CutnIntactWing_steadyNcutWBfreq_realFlyData.svg')
plot2svg('AdAiRatio_vs_S2nS3_CutnIntactWing_steadyNcutWBfreq_realFlyData.svg')

cd ..





