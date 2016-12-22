
%% steady wb freq

clear
clc
close all
warning off

% maximum possible wingbeat
load('WBdataset_steady_1603WBs.mat')
stroke_amp_steady = max(stroke_wb_steady_bins_meanCIstd(:,1))-min(stroke_wb_steady_bins_meanCIstd(:,1));
A_ratio_max = 178/stroke_amp_steady;

%% const
% A_min = 0.5;
% A_max = 1.5;
AdAiRatio_min = 1;
AdAiRatio_max = 2;

S2_min = 0;
S2_max = 1;
S3_min = 0;
S3_max = 1;

F_min = -1;
F_max = 0;
M_min = -.4;
M_max = .4;

Ftot_min = -2;
Ftot_max = 0;
Mtot_min = -.4;
Mtot_max = .4;

Fd_min = -1;
Fd_max = 0;
Md_min = -.4;
Md_max = .4;

Fi_min = -1;
Fi_max = 0;
Mi_min = -.4;
Mi_max = .4;

plot_on = 0;
plot_on = 1;

%% colormap for Aratio: blue to white to red
cmap_RdBu =cbrewer('div','RdBu',100);
cmap_BuRd = flipud(cmap_RdBu);

% black at edges of cmap
cmap_BuRd_edge = cmap_BuRd;
cmap_BuRd_edge(1,:) = [0 0 0];
cmap_BuRd_edge(end,:) = [0 0 0];

cmap_Aratio = cmap_BuRd;

%% colormap for A+: reds
% cmap_Reds =cbrewer('seq','Reds',100);
% cmap_AdAi = cmap_Reds;
% 
% cmap_YlOrRd =cbrewer('seq','YlOrRd',100);
% cmap_AdAi = cmap_YlOrRd;

cmap_hot =colormap(hot(100));
cmap_hot = flipud(cmap_hot);
cmap_AdAi = cmap_hot;

% %% colormap for A+: pink to white to green
% cmap_BrBG =cbrewer('div','BrBG',100);
% cmap_BGBr = flipud(cmap_BrBG);
% 
% % black at edges of cmap
% cmap_BGBr_edge = cmap_BGBr;
% cmap_BGBr_edge(1,:) = [0 0 0];
% cmap_BGBr_edge(end,:) = [0 0 0];
% 
% cmap_AdAi = cmap_BGBr;

%% colormap for Forces: purple to white to green
cmap_PRGn =cbrewer('div','PRGn',100);
cmap_GnPR = flipud(cmap_PRGn);

% black at edges of cmap
cmap_GnPR_edge = cmap_GnPR;
cmap_GnPR_edge(1,:) = [0 0 0];
cmap_GnPR_edge(end,:) = [0 0 0];

cmap_F = cmap_GnPR;
cmap_F_neg = cmap_PRGn;

%% colormap for Torques: brown to white to turquise
cmap_PiYG =cbrewer('div','PiYG',100);
cmap_YGPi = flipud(cmap_PiYG);

% black at edges of cmap
cmap_YGPi_edge = cmap_YGPi;
cmap_YGPi_edge(1,:) = [0 0 0];
cmap_YGPi_edge(end,:) = [0 0 0];

cmap_T = cmap_YGPi;

%% load DB's
% fly data
var_file = dir('flyVar*');
load(var_file.name)
% steady Wb F&M&CoM data
load('steadyWB_FnMnCoM_data.mat')

% % roboflyDB of Non-cut wing (AmpReduce)
% load('roboflyDB_NONcutWing_FnM_vs_ReducedAmpStrokeRatio_INCcaliCF_INCclippedWingFreq.mat')

% roboflyDB of Cut wing FnM vs Amp&S2&S3 Linear Surface Fit INC interaction
load('roboflyDB_CutWing_FnM_vs_StrokeAmplitude_vs_S2nS3_INCcaliCF_LinSurfFitInterac.mat')

%% fit param's
% cut wing vertical force vs Amp (parabolic fit) & S2 (linear fit)
Fz_Amp_S2_SurfFit_coeffs = coeffvalues(Fz_Amp_S2_SurfFit);
% Fz_Amp_fit_coeffs = Fz_Amp_fit;

% cut wing roll torque vs Amp (parabolic fit) & S3 (linear fit)
Mx_Amp_S3_SurfFit_coeffs = coeffvalues(Mx_MinSteady_Amp_S3_SurfFit);
% MxMinSteady_Amp_fit_coeffs = MxMinSteady_Amp_fit;

%% calc Amp's from vertical force and roll torque balance
calc_AdAi_FiFd_MiMd_AT_Fz0nMx0_AT_WBfreq_singleEq
% calc_AdAi_FiFd_MiMd_AT_Fz0nMx0_AT_WBfreq
% % calc_AdAi_FiFd_MiMd_AT_Fz0nMx0
% % % calc_AiAd_FiFd_MiMd_AT_Fz0nMx0

% colorbar centered at uncut equilibrium
A_max = A_ratio_max;
A_min = A_ratio_max - 2*(A_ratio_max-Aratio_intact);

% cmap_Aratio & cmap_AdAi: centered at uncut equilibrium & black at Amp > 178deg
n_Amax = round(99/(A_max-A_min)*(A_ratio_max-A_min)+1);
if n_Amax <= length(cmap_Aratio)
    cmap_Aratio(n_Amax:end,:)=0;
%     cmap_AdAi(n_Amax:end,:)=0;
end

%% plot Ai&Ad Fi&Fd Mi&Md @ weight support & zero roll torque
% plot_AiAd_AT_Fz0nMx0_LinSurfFitInterac_steadyWBfreq
fig_nr = 1
plot_AdAi_AT_Fz0nMx0_LinSurfFitInterac_steadyWBfreq

fig_nr = 2;
plot_FtotalFiFd_MtotalMiMd_AT_Fz0nMx0_LinSurfFitInterac

fig_nr = 3;
plot_Fd_Md_AT_Fz0nMx0_LinSurfFitInterac

figure(4)
plot_FdFlip_AT_Fz0nMx0_LinSurfFitInterac

%% save data
save('roboflyDB_CutAndIntactWing_EqSolved_AnFnM_vs_S2nS3_steadyWBfreq.mat')

%% comparison with real fly data
load('WBdataset_ClipNintact_wingbeat_kin_S2S3AmpFunc.mat')

figure(5)
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
    
    if clip_type_mean(i) > 1.5
        plot(S2_ratio_mean(i),S3_ratio_mean(i),'dk','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',8)
    else
        plot(S2_ratio_mean(i),S3_ratio_mean(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',8)
    end
        
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
    
    if clip_type_mean(i) > 1.5
        plot(S2_ratio_mean(i),S3_ratio_mean(i),'dk','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',8)
    else
        plot(S2_ratio_mean(i),S3_ratio_mean(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',8)
    end
end
axis equal
axis tight
axis([S2_min,S2_max,S3_min,S3_max])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
xlabel('S_2')
ylabel('S_3')

figure(6)
subplot(1,2,1)
hold on
for i = 1:length(Astroke_ratio_clip_intact_mean)
    color_nr = round(99/(AdAiRatio_max-AdAiRatio_min)*(Astroke_ratio_clip_intact_mean(i)-AdAiRatio_min)+1);
    if color_nr<1
        color_nr=1
    end
    
    if clip_type_mean(i) > 1.5
        plot(S2_ratio_mean(i),S3_ratio_mean(i),'dk','markerfacecolor',cmap_AdAi(color_nr,:),'markersize',8)
    else
        plot(S2_ratio_mean(i),S3_ratio_mean(i),'ok','markerfacecolor',cmap_AdAi(color_nr,:),'markersize',8)
    end
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

A_ratio_max = 178/stroke_amp_steady;
% A_ratio0 = .72;

%% const
% A_min = 0.5;
% A_max = 1.5;
AdAiRatio_min = 1;
AdAiRatio_max = 2;

S2_min = 0;
S2_max = 1;
S3_min = 0;
S3_max = 1;

F_min = -1;
F_max = 0;
M_min = -.4;
M_max = .4;

Ftot_min = -2;
Ftot_max = 0;
Mtot_min = -.4;
Mtot_max = .4;

Fd_min = -1;
Fd_max = 0;
Md_min = -.4;
Md_max = .4;

Fi_min = -1;
Fi_max = 0;
Mi_min = -.4;
Mi_max = .4;

plot_on = 0;
plot_on = 1;

%% colormap for Aratio: blue to white to red
cmap_RdBu =cbrewer('div','RdBu',100);
cmap_BuRd = flipud(cmap_RdBu);

% black at edges of cmap
cmap_BuRd_edge = cmap_BuRd;
cmap_BuRd_edge(1,:) = [0 0 0];
cmap_BuRd_edge(end,:) = [0 0 0];

cmap_Aratio = cmap_BuRd;

%% colormap for A+: reds
% cmap_Reds =cbrewer('seq','Reds',100);
% cmap_AdAi = cmap_Reds;
% 
% cmap_YlOrRd =cbrewer('seq','YlOrRd',100);
% cmap_AdAi = cmap_YlOrRd;

cmap_hot =colormap(hot(100));
cmap_hot = flipud(cmap_hot);
cmap_AdAi = cmap_hot;

% %% colormap for A+: pink to white to green
% cmap_BrBG =cbrewer('div','BrBG',100);
% cmap_BGBr = flipud(cmap_BrBG);
% 
% % black at edges of cmap
% cmap_BGBr_edge = cmap_BGBr;
% cmap_BGBr_edge(1,:) = [0 0 0];
% cmap_BGBr_edge(end,:) = [0 0 0];
% 
% cmap_AdAi = cmap_BGBr;

%% colormap for Forces: purple to white to green
cmap_PRGn =cbrewer('div','PRGn',100);
cmap_GnPR = flipud(cmap_PRGn);

% black at edges of cmap
cmap_GnPR_edge = cmap_GnPR;
cmap_GnPR_edge(1,:) = [0 0 0];
cmap_GnPR_edge(end,:) = [0 0 0];

cmap_F = cmap_GnPR;
cmap_F_neg = cmap_PRGn;

%% colormap for Torques: brown to white to turquise
cmap_PiYG =cbrewer('div','PiYG',100);
cmap_YGPi = flipud(cmap_PiYG);

% black at edges of cmap
cmap_YGPi_edge = cmap_YGPi;
cmap_YGPi_edge(1,:) = [0 0 0];
cmap_YGPi_edge(end,:) = [0 0 0];

cmap_T = cmap_YGPi;

%% load DB's
% fly data
var_file = dir('flyVar*');
load(var_file.name)
% steady Wb F&M&CoM data
load('steadyWB_FnMnCoM_data.mat')

% % roboflyDB of Non-cut wing (AmpReduce)
% load('roboflyDB_NONcutWing_FnM_vs_ReducedAmpStrokeRatio_INCcaliCF_INCclippedWingFreq.mat')

% roboflyDB of Cut wing FnM vs Amp&S2&S3 Linear Surface Fit INC interaction
load('roboflyDB_CutWing_FnM_vs_StrokeAmplitude_vs_S2nS3_INCcaliCF_LinSurfFitInterac.mat')

%% fit param's
% cut wing vertical force vs Amp (parabolic fit) & S2 (linear fit)
Fz_Amp_S2_SurfFit_coeffs = coeffvalues(Fz_Amp_S2_SurfFit_freqMod);
% Fz_Amp_fit_coeffs = Fz_Amp_fit_freqMod;

% cut wing roll torque vs Amp (parabolic fit) & S3 (linear fit)
Mx_Amp_S3_SurfFit_coeffs = coeffvalues(Mx_MinSteady_Amp_S3_SurfFit_freqMod);
% MxMinSteady_Amp_fit_coeffs = MxMinSteady_Amp_fit_freqMod;

%% calc Amp's from vertical force and roll torque balance
calc_AdAi_FiFd_MiMd_AT_Fz0nMx0_AT_WBfreq_singleEq
% calc_AdAi_FiFd_MiMd_AT_Fz0nMx0_AT_WBfreq
% % calc_AdAi_FiFd_MiMd_AT_Fz0nMx0
% % % calc_AiAd_FiFd_MiMd_AT_Fz0nMx0

% colorbar centered at uncut equilibrium
A_max = A_ratio_max;
A_min = A_ratio_max - 2*(A_ratio_max-Aratio_intact);

% cmap_Aratio & cmap_AdAi: centered at uncut equilibrium & black at Amp > 178deg
n_Amax = round(99/(A_max-A_min)*(A_ratio_max-A_min)+1);
if n_Amax <= length(cmap_Aratio)
    cmap_Aratio(n_Amax:end,:)=0;
%     cmap_AdAi(n_Amax:end,:)=0;
end

%% plot Ai&Ad Fi&Fd Mi&Md @ weight support & zero roll torque
% plot_AiAd_AT_Fz0nMx0_LinSurfFitInterac_clippedFlyWBfreq
fig_nr = 101
plot_AdAi_AT_Fz0nMx0_LinSurfFitInterac_clippedFlyWBfreq

fig_nr = 102;
plot_FtotalFiFd_MtotalMiMd_AT_Fz0nMx0_LinSurfFitInterac

fig_nr = 103;
plot_Fd_Md_AT_Fz0nMx0_LinSurfFitInterac

figure(104)
plot_FdFlip_AT_Fz0nMx0_LinSurfFitInterac

%% save data
save('roboflyDB_CutAndIntactWing_EqSolved_AnFnM_vs_S2nS3_clippedFlyWBfreq.mat')

%% add real fly data
load('WBdataset_ClipNintact_wingbeat_kin_S2S3AmpFunc.mat')

figure(105)
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
    
    if clip_type_mean(i) > 1.5
        plot(S2_ratio_mean(i),S3_ratio_mean(i),'dk','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',8)
    else
        plot(S2_ratio_mean(i),S3_ratio_mean(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',8)
    end
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
    
    if clip_type_mean(i) > 1.5
        plot(S2_ratio_mean(i),S3_ratio_mean(i),'dk','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',8)
    else
        plot(S2_ratio_mean(i),S3_ratio_mean(i),'ok','markerfacecolor',cmap_Aratio(color_nr,:),'markersize',8)
    end
end
axis equal
axis tight
axis([S2_min,S2_max,S3_min,S3_max])
set(gca,'xtick',S2_min:(S2_max-S2_min)/2:S2_max)
set(gca,'ytick',S3_min:(S3_max-S3_min)/2:S3_max)
xlabel('S_2')
ylabel('S_3')

figure(106)
subplot(1,2,2)
hold on
for i = 1:length(Astroke_ratio_clip_intact_mean)
    color_nr = round(99/(AdAiRatio_max-AdAiRatio_min)*(Astroke_ratio_clip_intact_mean(i)-AdAiRatio_min)+1);
    if color_nr<1
        color_nr=1
    end
    
    if clip_type_mean(i) > 1.5
        plot(S2_ratio_mean(i),S3_ratio_mean(i),'dk','markerfacecolor',cmap_AdAi(color_nr,:),'markersize',8)
    else
        plot(S2_ratio_mean(i),S3_ratio_mean(i),'ok','markerfacecolor',cmap_AdAi(color_nr,:),'markersize',8)
    end
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

%% steady WB freq
figure(1)
saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_steadyWBfreq_LinSurfFitInterac.fig')
saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_steadyWBfreq_LinSurfFitInterac.png')
% saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_steadyWBfreq_LinSurfFitInterac.svg')
plot2svg('AinAd_vs_S2nS3_CutnIntactWing_steadyWBfreq_LinSurfFitInterac.svg')

figure(11)
saveas(gcf,'AdAiRatio_vs_S2nS3_CutnIntactWing_steadyWBfreq_LinSurfFitInterac.fig')
saveas(gcf,'AdAiRatio_vs_S2nS3_CutnIntactWing_steadyWBfreq_LinSurfFitInterac.png')
% % saveas(gcf,'AdAiRatio_vs_S2nS3_CutnIntactWing_steadyWBfreq_LinSurfFitInterac.svg')
% plot2svg('AdAiRatio_vs_S2nS3_CutnIntactWing_steadyWBfreq_LinSurfFitInterac.svg')

figure(2)
saveas(gcf,'Fz_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.fig')
saveas(gcf,'Fz_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.png')
% saveas(gcf,'Fz_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.svg')
plot2svg('Fz_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.svg')

figure(12)
saveas(gcf,'Mx_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.fig')
saveas(gcf,'Mx_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.png')
% saveas(gcf,'Mx_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.svg')
plot2svg('Mx_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.svg')

figure(3)
saveas(gcf,'Fd_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.fig')
saveas(gcf,'Fd_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.png')
% saveas(gcf,'Fd_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.svg')
plot2svg('Fd_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.svg')

figure(13)
saveas(gcf,'Md_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.fig')
saveas(gcf,'Md_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.png')
% saveas(gcf,'Md_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.svg')
plot2svg('Md_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.svg')

figure(4)
saveas(gcf,'FdFlip_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.fig')
saveas(gcf,'FdFlip_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.png')
% saveas(gcf,'FdFlip_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.svg')
plot2svg('FdFlip_vs_S2nS3_steadyWBfreq_LinSurfFitInterac.svg')

figure(5)
saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_steadyWBfreq_realFlyData.fig')
saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_steadyWBfreq_realFlyData.png')
% saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_steadyWBfreq_realFlyData.svg')
plot2svg('AinAd_vs_S2nS3_CutnIntactWing_steadyWBfreq_realFlyData.svg')

figure(6)
saveas(gcf,'AdAiRatio_vs_S2nS3_CutnIntactWing_steadyNcutWBfreq_realFlyData.fig')
saveas(gcf,'AdAiRatio_vs_S2nS3_CutnIntactWing_steadyNcutWBfreq_realFlyData.png')
% % saveas(gcf,'AdAiRatio_vs_S2nS3_CutnIntactWing_steadyNcutWBfreq_realFlyData.svg')
plot2svg('AdAiRatio_vs_S2nS3_CutnIntactWing_steadyNcutWBfreq_realFlyData.svg')

%% cut WB freq
figure(101)
saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_cutWBfreq_LinSurfFitInterac.fig')
saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_cutWBfreq_LinSurfFitInterac.png')
% saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_cutWBfreq_LinSurfFitInterac.svg')
plot2svg('AinAd_vs_S2nS3_CutnIntactWing_cutWBfreq_LinSurfFitInterac.svg')

figure(111)
saveas(gcf,'AdAiRatio_vs_S2nS3_CutnIntactWing_cutWBfreq_LinSurfFitInterac.fig')
saveas(gcf,'AdAiRatio_vs_S2nS3_CutnIntactWing_cutWBfreq_LinSurfFitInterac.png')
% saveas(gcf,'AdAiRatio_vs_S2nS3_CutnIntactWing_cutWBfreq_LinSurfFitInterac.svg')
plot2svg('AdAiRatio_vs_S2nS3_CutnIntactWing_cutWBfreq_LinSurfFitInterac.svg')

figure(102)
saveas(gcf,'Fz_vs_S2nS3_cutWBfreq_LinSurfFitInterac.fig')
saveas(gcf,'Fz_vs_S2nS3_cutWBfreq_LinSurfFitInterac.png')
% saveas(gcf,'Fz_vs_S2nS3_cutWBfreq_LinSurfFitInterac.svg')
plot2svg('Fz_vs_S2nS3_cutWBfreq_LinSurfFitInterac.svg')

figure(112)
saveas(gcf,'Mx_vs_S2nS3_cutWBfreq_LinSurfFitInterac.fig')
saveas(gcf,'Mx_vs_S2nS3_cutWBfreq_LinSurfFitInterac.png')
% saveas(gcf,'Mx_vs_S2nS3_cutWBfreq_LinSurfFitInterac.svg')
% plot2svg('Mx_vs_S2nS3_cutWBfreq_LinSurfFitInterac.svg')

figure(103)
saveas(gcf,'Fd_vs_S2nS3_cutWBfreq_LinSurfFitInterac.fig')
saveas(gcf,'Fd_vs_S2nS3_cutWBfreq_LinSurfFitInterac.png')
% saveas(gcf,'Fd_vs_S2nS3_cutWBfreq_LinSurfFitInterac.svg')
plot2svg('Fd_vs_S2nS3_cutWBfreq_LinSurfFitInterac.svg')

figure(113)
saveas(gcf,'Md_vs_S2nS3_cutWBfreq_LinSurfFitInterac.fig')
saveas(gcf,'Md_vs_S2nS3_cutWBfreq_LinSurfFitInterac.png')
% saveas(gcf,'Md_vs_S2nS3_cutWBfreq_LinSurfFitInterac.svg')
plot2svg('Md_vs_S2nS3_cutWBfreq_LinSurfFitInterac.svg')

figure(104)
saveas(gcf,'FdFlip_vs_S2nS3_cutWBfreq_LinSurfFitInterac.fig')
saveas(gcf,'FdFlip_vs_S2nS3_cutWBfreq_LinSurfFitInterac.png')
% saveas(gcf,'FdFlip_vs_S2nS3_cutWBfreq_LinSurfFitInterac.svg')
plot2svg('FdFlip_vs_S2nS3_cutWBfreq_LinSurfFitInterac.svg')

figure(105)
saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_cutWBfreq_realFlyData.fig')
saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_cutWBfreq_realFlyData.png')
% saveas(gcf,'AinAd_vs_S2nS3_CutnIntactWing_cutWBfreq_realFlyData.svg')
plot2svg('AinAd_vs_S2nS3_CutnIntactWing_cutWBfreq_realFlyData.svg')

figure(106)
saveas(gcf,'AdAiRatio_vs_S2nS3_CutnIntactWing_steadyNcutWBfreq_realFlyData.fig')
saveas(gcf,'AdAiRatio_vs_S2nS3_CutnIntactWing_steadyNcutWBfreq_realFlyData.png')
% % saveas(gcf,'AdAiRatio_vs_S2nS3_CutnIntactWing_steadyNcutWBfreq_realFlyData.svg')
plot2svg('AdAiRatio_vs_S2nS3_CutnIntactWing_steadyNcutWBfreq_realFlyData.svg')

cd ..





