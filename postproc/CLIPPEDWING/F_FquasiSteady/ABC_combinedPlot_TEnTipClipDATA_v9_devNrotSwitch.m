clear;
clc;
close all
warning off

cmap_MODs_red = [1 1 1; 1 1 0; 1 .5 0; 1 0 0; .5 0 0; 0 0 0]; 
cmap_MODs_blue = [1 1 1; 0 1 1; 0 .5 1; 0 0 1; 0 0 .5; 0 0 0]; 

peakloc = 0;
rot_on = 1;

Aratio_max = 1.3;
Aratio_min = 1;

%% load data data
% freq & roll fit type
% fit_type = 1   % asymp10
fit_type = 2   % spline999;
% fit_type = 3   % linear;
% fit_type = 4   % power2;
    
% load('fit_type.mat')
% fit_type = fit_type+1

%% plot TIP cut wing data
    if fit_type == 1
        load(['allMODs_TipClip_indivFreqNrollFit_Asym10_roton',num2str(rot_on),'.mat'])
    elseif fit_type == 2
        load(['allMODs_TipClip_indivFreqNrollFit_spline999_roton',num2str(rot_on),'.mat'])
    elseif fit_type == 3
        load(['allMODs_TipClip_indivFreqNrollFit_linear_roton',num2str(rot_on),'.mat'])
    elseif fit_type == 4
        load(['allMODs_TipClip_indivFreqNrollFit_power2_roton',num2str(rot_on),'.mat'])
    end

    plot_on = 0;
    S2ratios_TEnTip = S2ratios;
    S3ratios_TEnTip = S3ratios;
    S2S3AmpRatioFuncs_TEnTip = S2S3AmpRatioFuncs;
    
    %% all&none
    Fx_damaged_mean_all_TEnTip = Fx_damaged_mean_all;
    Fy_damaged_mean_all_TEnTip = Fy_damaged_mean_all_rolled;
    Fz_damaged_mean_all_TEnTip = Fz_damaged_mean_all_rolled;
    
    Mx_damaged_mean_all_TEnTip = Mx_damaged_mean_all;
    My_damaged_mean_all_TEnTip = My_damaged_mean_all;
    Mz_damaged_mean_all_TEnTip = Mz_damaged_mean_all;
    
    Fx_damaged_mean_steady_TEnTip = Fx_damaged_mean_steady;
    Fy_damaged_mean_steady_TEnTip = Fy_damaged_mean_steady;
    Fz_damaged_mean_steady_TEnTip = Fz_damaged_mean_steady;
    
    Mx_damaged_mean_steady_TEnTip = Mx_damaged_mean_steady;
    My_damaged_mean_steady_TEnTip = My_damaged_mean_steady;
    Mz_damaged_mean_steady_TEnTip = Mz_damaged_mean_steady;
    
    %% MX: FREQ, ROLL, DEV, ROT, STROKE
    % freq
    Mx_damaged_mean_freqMOD_TEnTip = Mx_damaged_mean_freqMOD;
    My_damaged_mean_freqMOD_TEnTip = My_damaged_mean_freqMOD;
    Mz_damaged_mean_freqMOD_TEnTip = Mz_damaged_mean_freqMOD;
    
    %freq&rot
    Mx_damaged_mean_freqRotMOD_TEnTip = Mx_damaged_mean_freqRotMOD;
    My_damaged_mean_freqRotMOD_TEnTip = My_damaged_mean_freqRotMOD;
    Mz_damaged_mean_freqRotMOD_TEnTip = Mz_damaged_mean_freqRotMOD;
    
    %freq&dev
    Mx_damaged_mean_freqDevMOD_TEnTip = Mx_damaged_mean_freqDevMOD;
    My_damaged_mean_freqDevMOD_TEnTip = My_damaged_mean_freqDevMOD;
    Mz_damaged_mean_freqDevMOD_TEnTip = Mz_damaged_mean_freqDevMOD;
    
    %freq&rot&dev
    Mx_damaged_mean_freqRotDevMOD_TEnTip = Mx_damaged_mean_freqRotDevMOD;
    My_damaged_mean_freqRotDevMOD_TEnTip = My_damaged_mean_freqRotDevMOD;
    Mz_damaged_mean_freqRotDevMOD_TEnTip = Mz_damaged_mean_freqRotDevMOD;
    
    %% FZ: STROKE, ROT, DEV, ROLL, FREQ
    % stroke
    Fx_damaged_mean_StrokeMOD_TEnTip = Fx_damaged_mean_strokeMOD;
    Fy_damaged_mean_StrokeMOD_TEnTip = Fy_damaged_mean_strokeMOD;
    Fz_damaged_mean_StrokeMOD_TEnTip = Fz_damaged_mean_strokeMOD;
    
    Mx_damaged_mean_StrokeMOD_TEnTip = Mx_damaged_mean_strokeMOD;
    My_damaged_mean_StrokeMOD_TEnTip = My_damaged_mean_strokeMOD;
    Mz_damaged_mean_StrokeMOD_TEnTip = Mz_damaged_mean_strokeMOD;
    
    % stroke&dev
    Fx_damaged_mean_StrokeDevMOD_TEnTip = Fx_damaged_mean_strokeDevMOD;
    Fy_damaged_mean_StrokeDevMOD_TEnTip = Fy_damaged_mean_strokeDevMOD;
    Fz_damaged_mean_StrokeDevMOD_TEnTip = Fz_damaged_mean_strokeDevMOD;
    
    Mx_damaged_mean_StrokeDevMOD_TEnTip = Mx_damaged_mean_strokeDevMOD;
    My_damaged_mean_StrokeDevMOD_TEnTip = My_damaged_mean_strokeDevMOD;
    Mz_damaged_mean_StrokeDevMOD_TEnTip = Mz_damaged_mean_strokeDevMOD;
    
    % stroke&rot
    Fx_damaged_mean_StrokeRotMOD_TEnTip = Fx_damaged_mean_strokeRotMOD;
    Fy_damaged_mean_StrokeRotMOD_TEnTip = Fy_damaged_mean_strokeRotMOD;
    Fz_damaged_mean_StrokeRotMOD_TEnTip = Fz_damaged_mean_strokeRotMOD;
    
    Mx_damaged_mean_StrokeRotMOD_TEnTip = Mx_damaged_mean_strokeRotMOD;
    My_damaged_mean_StrokeRotMOD_TEnTip = My_damaged_mean_strokeRotMOD;
    Mz_damaged_mean_StrokeRotMOD_TEnTip = Mz_damaged_mean_strokeRotMOD;
    
    % stroke&dev&rot
    Fx_damaged_mean_StrokeDevRotMOD_TEnTip = Fx_damaged_mean_strokeDevRotMOD;
    Fy_damaged_mean_StrokeDevRotMOD_TEnTip = Fy_damaged_mean_strokeDevRotMOD;
    Fz_damaged_mean_StrokeDevRotMOD_TEnTip = Fz_damaged_mean_strokeDevRotMOD;
    
    Mx_damaged_mean_StrokeDevRotMOD_TEnTip = Mx_damaged_mean_strokeDevRotMOD;
    My_damaged_mean_StrokeDevRotMOD_TEnTip = My_damaged_mean_strokeDevRotMOD;
    Mz_damaged_mean_StrokeDevRotMOD_TEnTip = Mz_damaged_mean_strokeDevRotMOD;
    
    % stroke&dev&rot&roll
    Fx_damaged_mean_StrokeDevRotRollMOD_TEnTip = Fx_damaged_mean_strokeDevRotMOD;
    Fy_damaged_mean_StrokeDevRotRollMOD_TEnTip = Fy_damaged_mean_strokeDevRotMOD_rolled;
    Fz_damaged_mean_StrokeDevRotRollMOD_TEnTip = Fz_damaged_mean_strokeDevRotMOD_rolled;
    
    Mx_damaged_mean_StrokeDevRotRollMOD_TEnTip = Mx_damaged_mean_strokeDevRotMOD;
    My_damaged_mean_StrokeDevRotRollMOD_TEnTip = My_damaged_mean_strokeDevRotMOD;
    Mz_damaged_mean_StrokeDevRotRollMOD_TEnTip = Mz_damaged_mean_strokeDevRotMOD;
    
    %% PLOT ALL S2 & S3
    figure(1000)
    subplot(2,2,3)
    hold on
    for i = 1:length(S2S3AmpRatioFuncs)
        color_nr = round(99/(Aratio_max-Aratio_min)*(S2S3AmpRatioFuncs(i)-Aratio_min)+1);
        if color_nr<1
            color_nr=1
        elseif color_nr>size(cmap_AdAi,1)
            color_nr=size(cmap_AdAi,1)
        end

        plot(S2ratios(i),S3ratios(i),'ok','markerfacecolor',cmap_AdAi(color_nr,:),'markersize',10)
    end
    
    %% PLOT ALL MODs VS NO MODs
    figure(1)
    subplot(2,2,1)
    hold on
    plot(S2ratios,Fx_damaged_mean_steady,'ok','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S2ratios,-Fy_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','g')
    plot(S2ratios,-Fz_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','b')
    
    plot(S2ratios,Fx_damaged_mean_all,'ok','markersize',10,'markerfacecolor','r')
    plot(S2ratios,-Fy_damaged_mean_all_rolled,'ok','markersize',10,'markerfacecolor',[0 .5 0])
    plot(S2ratios,-Fz_damaged_mean_all_rolled,'ok','markersize',10,'markerfacecolor',[0 0 .5])
    
    subplot(2,2,2)
    hold on
    plot(S3ratios,Mx_damaged_mean_steady,'ok','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S3ratios,-My_damaged_mean_steady--My_damaged_mean_all(S3ratios==1),'ok','markersize',10,'markerfacecolor','g')
    plot(S3ratios,-Mz_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','b')
    
    plot(S3ratios,Mx_damaged_mean_all,'ok','markersize',10,'markerfacecolor','r')
    plot(S3ratios,-My_damaged_mean_all--My_damaged_mean_all(S3ratios==1),'ok','markersize',10,'markerfacecolor',[0 .5 0])
    plot(S3ratios,-Mz_damaged_mean_all,'ok','markersize',10,'markerfacecolor',[0 0 .5])
    
    %% compare MODs Mx COMBIMOD FREQ, ROLL, DEV, ROT, STROKE
    figure(4)
    subplot(2,2,2)
    hold on
    plot(S3ratios,Mx_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','w')
    plot(S3ratios,Mx_damaged_mean_freqMOD,'ok','markersize',10,'markerfacecolor','y')
    plot(S3ratios,Mx_damaged_mean_freqMOD,'ok','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S3ratios,Mx_damaged_mean_freqDevMOD,'ok','markersize',10,'markerfacecolor','r')
    plot(S3ratios,Mx_damaged_mean_freqRotDevMOD,'ok','markersize',10,'markerfacecolor',[.5 0 0])
    plot(S3ratios,Mx_damaged_mean_all,'ok','markersize',10,'markerfacecolor','k')
    
    %% compare MODs -Fz COMBIMOD STROKE, ROT, DEV, ROLL, FREQ
    figure(6)
    subplot(2,2,1)
    hold on
    plot(S2ratios,-Fz_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','w')
    plot(S2ratios,-Fz_damaged_mean_strokeMOD,'ok','markersize',10,'markerfacecolor','c')
    plot(S2ratios,-Fz_damaged_mean_strokeRotMOD,'ok','markersize',10,'markerfacecolor',[0 .5 1])
    plot(S2ratios,-Fz_damaged_mean_strokeDevRotMOD,'ok','markersize',10,'markerfacecolor','b')
    plot(S2ratios,-Fz_damaged_mean_strokeDevRotMOD_rolled,'ok','markersize',10,'markerfacecolor',[0 0 .5])
    plot(S2ratios,-Fz_damaged_mean_all_rolled,'ok','markersize',10,'markerfacecolor','k')
    
%% plot TE cut wing data
    if fit_type == 1
        load(['allMODs_TEclip_indivFreqNrollFit_Asym10_roton',num2str(rot_on),'.mat'])
    elseif fit_type == 2
        load(['allMODs_TEclip_indivFreqNrollFit_spline999_roton',num2str(rot_on),'.mat'])
    elseif fit_type == 3
        load(['allMODs_TEclip_indivFreqNrollFit_linear_roton',num2str(rot_on),'.mat'])
    elseif fit_type == 4
        load(['allMODs_TEclip_indivFreqNrollFit_power2_roton',num2str(rot_on),'.mat'])
    end

    S2ratios_TEnTip = [S2ratios_TEnTip;S2ratios];
    S3ratios_TEnTip = [S3ratios_TEnTip;S3ratios];
    S2S3AmpRatioFuncs_TEnTip = [S2S3AmpRatioFuncs_TEnTip;S2S3AmpRatioFuncs];

    %% ALL&NONE
    Fx_damaged_mean_all_TEnTip = [Fx_damaged_mean_all_TEnTip Fx_damaged_mean_all];
    Fy_damaged_mean_all_TEnTip = [Fy_damaged_mean_all_TEnTip Fy_damaged_mean_all_rolled];
    Fz_damaged_mean_all_TEnTip = [Fz_damaged_mean_all_TEnTip Fz_damaged_mean_all_rolled];
    
    Mx_damaged_mean_all_TEnTip = [Mx_damaged_mean_all_TEnTip Mx_damaged_mean_all];
    My_damaged_mean_all_TEnTip = [My_damaged_mean_all_TEnTip My_damaged_mean_all];
    Mz_damaged_mean_all_TEnTip = [Mz_damaged_mean_all_TEnTip Mz_damaged_mean_all];
    
    Fx_damaged_mean_steady_TEnTip = [Fx_damaged_mean_steady_TEnTip Fx_damaged_mean_steady];
    Fy_damaged_mean_steady_TEnTip = [Fy_damaged_mean_steady_TEnTip Fy_damaged_mean_steady];
    Fz_damaged_mean_steady_TEnTip = [Fz_damaged_mean_steady_TEnTip Fz_damaged_mean_steady];
    
    Mx_damaged_mean_steady_TEnTip = [Mx_damaged_mean_steady_TEnTip Mx_damaged_mean_steady];
    My_damaged_mean_steady_TEnTip = [My_damaged_mean_steady_TEnTip My_damaged_mean_steady];
    Mz_damaged_mean_steady_TEnTip = [Mz_damaged_mean_steady_TEnTip Mz_damaged_mean_steady];
    
    %% MX: FREQ, ROLL, DEV, ROT, STROKE
    % freq
    Mx_damaged_mean_freqMOD_TEnTip = [Mx_damaged_mean_freqMOD_TEnTip Mx_damaged_mean_freqMOD];
    My_damaged_mean_freqMOD_TEnTip = [My_damaged_mean_freqMOD_TEnTip My_damaged_mean_freqMOD];
    Mz_damaged_mean_freqMOD_TEnTip = [Mz_damaged_mean_freqMOD_TEnTip Mz_damaged_mean_freqMOD];
    
    % freqROT
    Mx_damaged_mean_freqRotMOD_TEnTip = [Mx_damaged_mean_freqRotMOD_TEnTip Mx_damaged_mean_freqRotMOD];
    My_damaged_mean_freqRotMOD_TEnTip = [My_damaged_mean_freqRotMOD_TEnTip My_damaged_mean_freqRotMOD];
    Mz_damaged_mean_freqRotMOD_TEnTip = [Mz_damaged_mean_freqRotMOD_TEnTip Mz_damaged_mean_freqRotMOD];
    
    % freqDEV
    Mx_damaged_mean_freqDevMOD_TEnTip = [Mx_damaged_mean_freqDevMOD_TEnTip Mx_damaged_mean_freqDevMOD];
    My_damaged_mean_freqDevMOD_TEnTip = [My_damaged_mean_freqDevMOD_TEnTip My_damaged_mean_freqDevMOD];
    Mz_damaged_mean_freqDevMOD_TEnTip = [Mz_damaged_mean_freqDevMOD_TEnTip Mz_damaged_mean_freqDevMOD];
    
    % freqROTDEV
    Mx_damaged_mean_freqRotDevMOD_TEnTip = [Mx_damaged_mean_freqRotDevMOD_TEnTip Mx_damaged_mean_freqRotDevMOD];
    My_damaged_mean_freqRotDevMOD_TEnTip = [My_damaged_mean_freqRotDevMOD_TEnTip My_damaged_mean_freqRotDevMOD];
    Mz_damaged_mean_freqRotDevMOD_TEnTip = [Mz_damaged_mean_freqRotDevMOD_TEnTip Mz_damaged_mean_freqRotDevMOD];
    
    %% STROKE, ROT, DEV, ROLL, FREQ
    % str
    Fx_damaged_mean_StrokeMOD_TEnTip = [Fx_damaged_mean_StrokeMOD_TEnTip Fx_damaged_mean_strokeMOD];
    Fy_damaged_mean_StrokeMOD_TEnTip = [Fy_damaged_mean_StrokeMOD_TEnTip Fy_damaged_mean_strokeMOD];
    Fz_damaged_mean_StrokeMOD_TEnTip = [Fz_damaged_mean_StrokeMOD_TEnTip Fz_damaged_mean_strokeMOD];
    
    Mx_damaged_mean_StrokeMOD_TEnTip = [Mx_damaged_mean_StrokeMOD_TEnTip Mx_damaged_mean_strokeMOD];
    My_damaged_mean_StrokeMOD_TEnTip = [My_damaged_mean_StrokeMOD_TEnTip My_damaged_mean_strokeMOD];
    Mz_damaged_mean_StrokeMOD_TEnTip = [Mz_damaged_mean_StrokeMOD_TEnTip Mz_damaged_mean_strokeMOD];

    % str&dev
    Fx_damaged_mean_StrokeDevMOD_TEnTip = [Fx_damaged_mean_StrokeDevMOD_TEnTip Fx_damaged_mean_strokeDevMOD];
    Fy_damaged_mean_StrokeDevMOD_TEnTip = [Fy_damaged_mean_StrokeDevMOD_TEnTip Fy_damaged_mean_strokeDevMOD];
    Fz_damaged_mean_StrokeDevMOD_TEnTip = [Fz_damaged_mean_StrokeDevMOD_TEnTip Fz_damaged_mean_strokeDevMOD];
    
    Mx_damaged_mean_StrokeDevMOD_TEnTip = [Mx_damaged_mean_StrokeDevMOD_TEnTip Mx_damaged_mean_strokeDevMOD];
    My_damaged_mean_StrokeDevMOD_TEnTip = [My_damaged_mean_StrokeDevMOD_TEnTip My_damaged_mean_strokeDevMOD];
    Mz_damaged_mean_StrokeDevMOD_TEnTip = [Mz_damaged_mean_StrokeDevMOD_TEnTip Mz_damaged_mean_strokeDevMOD];

    % str&rot
    Fx_damaged_mean_StrokeRotMOD_TEnTip = [Fx_damaged_mean_StrokeRotMOD_TEnTip Fx_damaged_mean_strokeRotMOD];
    Fy_damaged_mean_StrokeRotMOD_TEnTip = [Fy_damaged_mean_StrokeRotMOD_TEnTip Fy_damaged_mean_strokeRotMOD];
    Fz_damaged_mean_StrokeRotMOD_TEnTip = [Fz_damaged_mean_StrokeRotMOD_TEnTip Fz_damaged_mean_strokeRotMOD];
    
    Mx_damaged_mean_StrokeDevMOD_TEnTip = [Mx_damaged_mean_StrokeDevMOD_TEnTip Mx_damaged_mean_strokeDevMOD];
    My_damaged_mean_StrokeDevMOD_TEnTip = [My_damaged_mean_StrokeDevMOD_TEnTip My_damaged_mean_strokeDevMOD];
    Mz_damaged_mean_StrokeDevMOD_TEnTip = [Mz_damaged_mean_StrokeDevMOD_TEnTip Mz_damaged_mean_strokeDevMOD];

    % str&dev&rot
    Fx_damaged_mean_StrokeDevRotMOD_TEnTip = [Fx_damaged_mean_StrokeDevRotMOD_TEnTip Fx_damaged_mean_strokeDevRotMOD];
    Fy_damaged_mean_StrokeDevRotMOD_TEnTip = [Fy_damaged_mean_StrokeDevRotMOD_TEnTip Fy_damaged_mean_strokeDevRotMOD];
    Fz_damaged_mean_StrokeDevRotMOD_TEnTip = [Fz_damaged_mean_StrokeDevRotMOD_TEnTip Fz_damaged_mean_strokeDevRotMOD];
    
    Mx_damaged_mean_StrokeDevRotMOD_TEnTip = [Mx_damaged_mean_StrokeDevRotMOD_TEnTip Mx_damaged_mean_strokeDevRotMOD];
    My_damaged_mean_StrokeDevRotMOD_TEnTip = [My_damaged_mean_StrokeDevRotMOD_TEnTip My_damaged_mean_strokeDevRotMOD];
    Mz_damaged_mean_StrokeDevRotMOD_TEnTip = [Mz_damaged_mean_StrokeDevRotMOD_TEnTip Mz_damaged_mean_strokeDevRotMOD];

    % str&dev&rot&roll
    Fx_damaged_mean_StrokeDevRotRollMOD_TEnTip = [Fx_damaged_mean_StrokeDevRotRollMOD_TEnTip Fx_damaged_mean_strokeDevRotMOD];
    Fy_damaged_mean_StrokeDevRotRollMOD_TEnTip = [Fy_damaged_mean_StrokeDevRotRollMOD_TEnTip Fy_damaged_mean_strokeDevRotMOD_rolled];
    Fz_damaged_mean_StrokeDevRotRollMOD_TEnTip = [Fz_damaged_mean_StrokeDevRotRollMOD_TEnTip Fz_damaged_mean_strokeDevRotMOD_rolled];
    
    Mx_damaged_mean_StrokeDevRotRollMOD_TEnTip = [Mx_damaged_mean_StrokeDevRotRollMOD_TEnTip Mx_damaged_mean_strokeDevRotMOD];
    My_damaged_mean_StrokeDevRotRollMOD_TEnTip = [My_damaged_mean_StrokeDevRotRollMOD_TEnTip My_damaged_mean_strokeDevRotMOD];
    Mz_damaged_mean_StrokeDevRotRollMOD_TEnTip = [Mz_damaged_mean_StrokeDevRotRollMOD_TEnTip Mz_damaged_mean_strokeDevRotMOD];

    %% PLOT ALL S2 & S3
    figure(1000)
    subplot(2,2,3)
    hold on
    for i = 1:length(S2S3AmpRatioFuncs)
        color_nr = round(99/(Aratio_max-Aratio_min)*(S2S3AmpRatioFuncs(i)-Aratio_min)+1);
        if color_nr<1
            color_nr=1
        elseif color_nr>size(cmap_AdAi,1)
            color_nr=size(cmap_AdAi,1)
        end

        plot(S2ratios(i),S3ratios(i),'dk','markerfacecolor',cmap_AdAi(color_nr,:),'markersize',10)
    end
    
    plot(1,1,'sk','markerfacecolor','w','markersize',10)

%     legend('Tip clip','TE clip','location','NW')
    xlabel('S2 ratio')
    ylabel('S3 ratio')
    axis equal
    axis([0.5 1 0.5 1])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',0:.5:1)
    
    %% ALL MODs VS NO MODs
    figure(1)
    subplot(2,2,1)
    hold on
    plot(S2ratios,Fx_damaged_mean_steady,'dk','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S2ratios,-Fy_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','g')
    plot(S2ratios,-Fz_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','b')
    
    plot(S2ratios,Fx_damaged_mean_all,'dk','markersize',10,'markerfacecolor','r')
    plot(S2ratios,-Fy_damaged_mean_all_rolled,'dk','markersize',10,'markerfacecolor',[0 .5 0])
    plot(S2ratios,-Fz_damaged_mean_all_rolled,'dk','markersize',10,'markerfacecolor',[0 0 .5])
    
    subplot(2,2,2)
    hold on
    plot(S3ratios,Mx_damaged_mean_steady,'dk','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S3ratios,-My_damaged_mean_steady--My_damaged_mean_all(S3ratios==1),'dk','markersize',10,'markerfacecolor','g')
    plot(S3ratios,-Mz_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','b')
    
    plot(S3ratios,Mx_damaged_mean_all,'dk','markersize',10,'markerfacecolor','r')
    plot(S3ratios,-My_damaged_mean_all--My_damaged_mean_all(S3ratios==1),'dk','markersize',10,'markerfacecolor',[0 .5 0])
    plot(S3ratios,-Mz_damaged_mean_all,'dk','markersize',10,'markerfacecolor',[0 0 .5])
    
    subplot(2,2,1)
    xlabel('S2 ratio')
    ylabel('normalized forces F/mg')
    axis([0.5 1 0 1])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',0:1:1)
    
    subplot(2,2,2)
    xlabel('S3 ratio')
    ylabel('normalized torques T/mgl')
    axis([0.5 1 -.05 .15])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-.05:.05:.15)    
    
    subplot(2,2,3)
    hold on
    plot(0,0,'dk','markersize',10,'markerfacecolor',[1 .5 0])
    plot(0,0,'dk','markersize',10,'markerfacecolor','g')
    plot(0,0,'dk','markersize',10,'markerfacecolor','b')
    
    plot(0,0,'dk','markersize',10,'markerfacecolor','r')
    plot(0,0,'dk','markersize',10,'markerfacecolor',[0 .5 0])
    plot(0,0,'dk','markersize',10,'markerfacecolor',[0 0 .5])
    
    legend('x-Mod','y-Mod','z-Mod','x-steady','y-steady','z-steady','location','E')
    axis off
    
    %% compare MODs Mx COMBIMOD FREQ, ROLL, DEV, ROT, STROKE
    figure(4)
    subplot(2,2,2)
    hold on
    plot(S3ratios,Mx_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','w')
    plot(S3ratios,Mx_damaged_mean_freqMOD,'dk','markersize',10,'markerfacecolor','y')
    plot(S3ratios,Mx_damaged_mean_freqMOD,'dk','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S3ratios,Mx_damaged_mean_freqDevMOD,'dk','markersize',10,'markerfacecolor','r')
    plot(S3ratios,Mx_damaged_mean_freqRotDevMOD,'dk','markersize',10,'markerfacecolor',[.5 0 0])
    plot(S3ratios,Mx_damaged_mean_all,'dk','markersize',10,'markerfacecolor','k')

    xlabel('S3 ratio')
    ylabel('normalized roll torque Tx/mgl')
    axis([0.5 1 -.05 .15])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-.05:.05:.15)    
    
    subplot(2,2,3)
    hold on
    plot(0,0,'ok','markersize',10,'markerfacecolor','w')
    plot(0,0,'ok','markersize',10,'markerfacecolor','y')
    plot(0,0,'ok','markersize',10,'markerfacecolor',[1 .5 0])
    plot(0,0,'ok','markersize',10,'markerfacecolor','r')
    plot(0,0,'ok','markersize',10,'markerfacecolor',[.5 0 0])
    plot(0,0,'ok','markersize',10,'markerfacecolor','k')
    legend('normal','+freq','+roll','+dev','+rot','+stroke','location','E')
    colormap(cmap_MODs_red)
    colorbar
    axis off

    %% compare MODs -Fz & Mx COMBIMOD STROKE, ROT, DEV, ROLL, FREQ
    figure(6)
    subplot(2,2,1)
    hold on
    plot(S2ratios,-Fz_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','w')
    plot(S2ratios,-Fz_damaged_mean_strokeMOD,'dk','markersize',10,'markerfacecolor','c')
    plot(S2ratios,-Fz_damaged_mean_strokeRotMOD,'dk','markersize',10,'markerfacecolor',[0 .5 1])
    plot(S2ratios,-Fz_damaged_mean_strokeDevRotMOD,'dk','markersize',10,'markerfacecolor','b')
    plot(S2ratios,-Fz_damaged_mean_strokeDevRotMOD_rolled,'dk','markersize',10,'markerfacecolor',[0 0 .5])
    plot(S2ratios,-Fz_damaged_mean_all_rolled,'dk','markersize',10,'markerfacecolor','k')

    xlabel('S2 ratio')
    ylabel('normalized vertical force -Fz/mg')
    axis([0.5 1 .5 1.5])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-1.5:.5:1.5) 

    subplot(2,2,3)
    hold on
    plot(0,0,'ok','markersize',10,'markerfacecolor','w')
    plot(0,0,'ok','markersize',10,'markerfacecolor','c')
    plot(0,0,'ok','markersize',10,'markerfacecolor',[0 .5 1])
    plot(0,0,'ok','markersize',10,'markerfacecolor','b')
    plot(0,0,'ok','markersize',10,'markerfacecolor',[0 0 .5])
    plot(0,0,'ok','markersize',10,'markerfacecolor','k')
    legend('normal','+stroke','+rot','+dev','+roll','+freq','location','E')
    colormap(cmap_MODs_blue)
    colorbar
    axis off
    
    %% plot trendlines ALL&NONE
    figure(1)
    subplot(2,2,1)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,Fx_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color','r','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fy_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[0 .5 0],'linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[0 0 .5],'linewidth',2)

    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,Fx_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[1 .5 0],'linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fy_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color','g','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color','b','linewidth',2)

%     % adjusted
%     plot([min(S2ratios_TEnTip):.05:max(S2ratios_TEnTip) 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',2),[min(S2ratios_TEnTip):.05:max(S2ratios_TEnTip) 1]),'color',[0 0 .5],'linewidth',2)

    subplot(2,2,2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_all_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color','r','linewidth',2)
%     plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,-My_damaged_mean_all_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color',[0 .5 0],'linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,-Mz_damaged_mean_all_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color',[0 0 .5],'linewidth',2)

    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_steady_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color',[1 .5 0],'linewidth',2)
%     plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,-My_damaged_mean_steady_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color','g','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,-Mz_damaged_mean_steady_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color','b','linewidth',2)

    % adjusted
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,-My_damaged_mean_all_TEnTip'--My_damaged_mean_all(S3ratios==1),1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color',[0 .5 0],'linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,-My_damaged_mean_steady_TEnTip'--My_damaged_mean_steady(S3ratios==1),1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color','g','linewidth',2)

    %% plot trendlines MX: FREQ, ROLL, DEV, ROT, STROKE
    figure(4)
    subplot(2,2,2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_steady_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color',[.5 .5 .5],'linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color','y','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color',[1 .5 0],'linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqDevMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color','r','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqRotDevMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color',[.5 0 0],'linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_all_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color','k','linewidth',2)

    %% plot trendlines STROKE, DEV, ROT, ROLL, FREQ
    figure(6)
    subplot(2,2,1)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[.5 .5 .5],'linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_StrokeMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color','c','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_StrokeRotMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[0 .5 1],'linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_StrokeDevRotMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color','b','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_StrokeDevRotRollMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[0 0 .5],'linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color','k','linewidth',2)

%     % adjusted
%     plot([min(S2ratios_TEnTip):.05:max(S2ratios_TEnTip) 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',2),[min(S2ratios_TEnTip):.05:max(S2ratios_TEnTip) 1]),'color','k','linewidth',2)
% 

    %% percentage of roll equilibrium: FREQ, ROLL, DEV, ROT, STROKE
    Mx_steady = Mx_damaged_mean_steady_TEnTip(1);
    Mx_none = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_freq = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_freqDev = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqDevMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_freqRotDev = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqRotDevMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_freqRotDevStroke = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));

    freq2stroke.Mx_freq_perc_frds    = 100 / (-Mx_none) * (Mx_freq-Mx_none);
    freq2stroke.Mx_rot_perc_frds     = 100 / (-Mx_none) * (Mx_freqDev-Mx_freq);
    freq2stroke.Mx_dev_perc_frds     = 100 / (-Mx_none) * (Mx_freqRotDev-Mx_freqDev);
    freq2stroke.Mx_stroke_perc_frds  = 100 / (-Mx_none) * (Mx_freqRotDevStroke-Mx_freqRotDev);
    freq2stroke.Mx_all_perc_frds     = 100 / (-Mx_none) * (Mx_freqRotDevStroke-Mx_none);
    freq2stroke.Mx_sum_perc_frds     = freq2stroke.Mx_freq_perc_frds + freq2stroke.Mx_rot_perc_frds + freq2stroke.Mx_dev_perc_frds + freq2stroke.Mx_stroke_perc_frds;
    
    freq2stroke.Mx_normal = Mx_none;
    freq2stroke.Mx_WDR = Mx_freqRotDevStroke;
        
    My_steady = My_damaged_mean_steady_TEnTip(1);
    My_none = mean(polyval(polyfit(S3ratios_TEnTip,My_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)])) - My_steady;
    My_freq = mean(polyval(polyfit(S3ratios_TEnTip,My_damaged_mean_freqMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)])) - My_steady;
    My_freqDev = mean(polyval(polyfit(S3ratios_TEnTip,My_damaged_mean_freqDevMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)])) - My_steady;
    My_freqRotDev = mean(polyval(polyfit(S3ratios_TEnTip,My_damaged_mean_freqRotDevMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)])) - My_steady;
    My_freqRotDevStroke = mean(polyval(polyfit(S3ratios_TEnTip,My_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)])) - My_steady;

    freq2stroke.My_freq_perc_frds    = 100 / (-My_none) * (My_freq-My_none);
    freq2stroke.My_rot_perc_frds     = 100 / (-My_none) * (My_freqDev-My_freq);
    freq2stroke.My_dev_perc_frds     = 100 / (-My_none) * (My_freqRotDev-My_freqDev);
    freq2stroke.My_stroke_perc_frds  = 100 / (-My_none) * (My_freqRotDevStroke-My_freqRotDev);
    freq2stroke.My_all_perc_frds     = 100 / (-My_none) * (My_freqRotDevStroke-My_none);
    freq2stroke.My_sum_perc_frds     = freq2stroke.My_freq_perc_frds + freq2stroke.My_rot_perc_frds + freq2stroke.My_dev_perc_frds + freq2stroke.My_stroke_perc_frds;
    
    freq2stroke.My_normal = My_none;
    freq2stroke.My_WDR = My_freqRotDevStroke;
        
    Mz_steady = Mz_damaged_mean_steady_TEnTip(1);
    Mz_none = mean(polyval(polyfit(S3ratios_TEnTip,Mz_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mz_freq = mean(polyval(polyfit(S3ratios_TEnTip,Mz_damaged_mean_freqMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mz_freqDev = mean(polyval(polyfit(S3ratios_TEnTip,Mz_damaged_mean_freqDevMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mz_freqRotDev = mean(polyval(polyfit(S3ratios_TEnTip,Mz_damaged_mean_freqRotDevMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mz_freqRotDevStroke = mean(polyval(polyfit(S3ratios_TEnTip,Mz_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));

    freq2stroke.Mz_freq_perc_frds    = 100 / (-Mz_none) * (Mz_freq-Mz_none);
    freq2stroke.Mz_rot_perc_frds     = 100 / (-Mz_none) * (Mz_freqDev-Mz_freq);
    freq2stroke.Mz_dev_perc_frds     = 100 / (-Mz_none) * (Mz_freqRotDev-Mz_freqDev);
    freq2stroke.Mz_stroke_perc_frds  = 100 / (-Mz_none) * (Mz_freqRotDevStroke-Mz_freqRotDev);
    freq2stroke.Mz_all_perc_frds     = 100 / (-Mz_none) * (Mz_freqRotDevStroke-Mz_none);
    freq2stroke.Mz_sum_perc_frds     = freq2stroke.Mz_freq_perc_frds + freq2stroke.Mz_rot_perc_frds + freq2stroke.Mz_dev_perc_frds + freq2stroke.Mz_stroke_perc_frds;
    
    freq2stroke.Mz_normal = Mz_none;
    freq2stroke.Mz_WDR = Mz_freqRotDevStroke;
        
    %% percentage of weight support & sideways force: STROKE, DEV, ROT, ROLL, FREQ
    Fx_steady = -Fx_damaged_mean_steady_TEnTip(1);
    Fx_none = mean(polyval(polyfit(S2ratios_TEnTip,-Fx_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fx_Stroke = mean(polyval(polyfit(S2ratios_TEnTip,-Fx_damaged_mean_StrokeMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fx_StrokeRot = mean(polyval(polyfit(S2ratios_TEnTip,-Fx_damaged_mean_StrokeRotMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fx_StrokeDevRot = mean(polyval(polyfit(S2ratios_TEnTip,-Fx_damaged_mean_StrokeDevRotMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fx_StrokeDevRotRoll = mean(polyval(polyfit(S2ratios_TEnTip,-Fx_damaged_mean_StrokeDevRotRollMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fx_StrokeDevRotRollFreq = mean(polyval(polyfit(S2ratios_TEnTip,-Fx_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));

    stroke2freq.Fx_stroke_perc_sdrf  = 100* Fx_steady / (Fx_steady-Fx_none) * (Fx_Stroke-Fx_none);
    stroke2freq.Fx_dev_perc_sdrf     = 100* Fx_steady / (Fx_steady-Fx_none) * (Fx_StrokeRot-Fx_Stroke);
    stroke2freq.Fx_rot_perc_sdrf     = 100* Fx_steady / (Fx_steady-Fx_none) * (Fx_StrokeDevRot-Fx_StrokeRot);
    stroke2freq.Fx_roll_perc_sdrf    = 100* Fx_steady / (Fx_steady-Fx_none) * (Fx_StrokeDevRotRoll-Fx_StrokeDevRot);
    stroke2freq.Fx_freq_perc_sdrf    = 100* Fx_steady / (Fx_steady-Fx_none) * (Fx_StrokeDevRotRollFreq-Fx_StrokeDevRotRoll);
    stroke2freq.Fx_all_perc_sdrf     = 100* Fx_steady / (Fx_steady-Fx_none) * (Fx_StrokeDevRotRollFreq-Fx_none);
    stroke2freq.Fx_sum_perc_sdrf     = stroke2freq.Fx_stroke_perc_sdrf + stroke2freq.Fx_dev_perc_sdrf + stroke2freq.Fx_rot_perc_sdrf + stroke2freq.Fx_roll_perc_sdrf + stroke2freq.Fx_freq_perc_sdrf;
    
    stroke2freq.Fx_normal = Fx_none;
    stroke2freq.Fx_WDR = Fx_StrokeDevRotRollFreq;
    
    Fy_steady = -Fy_damaged_mean_steady_TEnTip(1);
    Fy_none = mean(polyval(polyfit(S2ratios_TEnTip,-Fy_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fy_Stroke = mean(polyval(polyfit(S2ratios_TEnTip,-Fy_damaged_mean_StrokeMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fy_StrokeRot = mean(polyval(polyfit(S2ratios_TEnTip,-Fy_damaged_mean_StrokeRotMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fy_StrokeDevRot = mean(polyval(polyfit(S2ratios_TEnTip,-Fy_damaged_mean_StrokeDevRotMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fy_StrokeDevRotRoll = mean(polyval(polyfit(S2ratios_TEnTip,-Fy_damaged_mean_StrokeDevRotRollMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fy_StrokeDevRotRollFreq = mean(polyval(polyfit(S2ratios_TEnTip,-Fy_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));

    stroke2freq.Fy_stroke_perc_sdrf  = 100 / (-Fy_none) * (Fy_Stroke-Fy_none);
    stroke2freq.Fy_dev_perc_sdrf     = 100 / (-Fy_none) * (Fy_StrokeRot-Fy_Stroke);
    stroke2freq.Fy_rot_perc_sdrf     = 100 / (-Fy_none) * (Fy_StrokeDevRot-Fy_StrokeRot);
    stroke2freq.Fy_roll_perc_sdrf    = 100 / (-Fy_none) * (Fy_StrokeDevRotRoll-Fy_StrokeDevRot);
    stroke2freq.Fy_freq_perc_sdrf    = 100 / (-Fy_none) * (Fy_StrokeDevRotRollFreq-Fy_StrokeDevRotRoll);
    stroke2freq.Fy_all_perc_sdrf     = 100 / (-Fy_none) * (Fy_StrokeDevRotRollFreq-Fy_none);
    stroke2freq.Fy_sum_perc_sdrf     = stroke2freq.Fy_stroke_perc_sdrf + stroke2freq.Fy_dev_perc_sdrf + stroke2freq.Fy_rot_perc_sdrf + stroke2freq.Fy_roll_perc_sdrf + stroke2freq.Fy_freq_perc_sdrf;
    
    stroke2freq.Fy_normal = Fy_none;
    stroke2freq.Fy_WDR = Fy_StrokeDevRotRollFreq;

    Fz_steady = -Fz_damaged_mean_steady_TEnTip(1);
    Fz_none = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_Stroke = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_StrokeMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_StrokeRot = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_StrokeRotMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_StrokeDevRot = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_StrokeDevRotMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_StrokeDevRotRoll = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_StrokeDevRotRollMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_StrokeDevRotRollFreq = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));

    stroke2freq.Fz_stroke_perc_sdrf  = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_Stroke-Fz_none);
    stroke2freq.Fz_dev_perc_sdrf     = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_StrokeRot-Fz_Stroke);
    stroke2freq.Fz_rot_perc_sdrf     = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_StrokeDevRot-Fz_StrokeRot);
    stroke2freq.Fz_roll_perc_sdrf    = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_StrokeDevRotRoll-Fz_StrokeDevRot);
    stroke2freq.Fz_freq_perc_sdrf    = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_StrokeDevRotRollFreq-Fz_StrokeDevRotRoll);
    stroke2freq.Fz_all_perc_sdrf     = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_StrokeDevRotRollFreq-Fz_none);
    stroke2freq.Fz_sum_perc_sdrf     = stroke2freq.Fz_stroke_perc_sdrf + stroke2freq.Fz_dev_perc_sdrf + stroke2freq.Fz_rot_perc_sdrf + stroke2freq.Fz_roll_perc_sdrf + stroke2freq.Fz_freq_perc_sdrf;
    
    stroke2freq.Fz_normal = Fz_none;
    stroke2freq.Fz_WDR = Fz_StrokeDevRotRollFreq;

    %% save plots
    %% save data
    if fit_type == 1
        save(['MODcontributionPerc2FznMx_indivFreqNrollFit_Asym10_roton',num2str(rot_on),'.mat'],'stroke2freq','freq2stroke')
        if plot_on == 1
            mkdir(['qsModel_FnM_indivFreqNrollFit_Asym10_roton',num2str(rot_on)])
            cd(['qsModel_FnM_indivFreqNrollFit_Asym10_roton',num2str(rot_on)])
        end
    elseif fit_type == 2
        save(['MODcontributionPerc2FznMx_indivFreqNrollFit_spline999_roton',num2str(rot_on),'.mat'],'stroke2freq','freq2stroke')
        if plot_on == 1
            mkdir(['qsModel_FnM_indivFreqNrollFit_spline999_roton',num2str(rot_on)])
            cd(['qsModel_FnM_indivFreqNrollFit_spline999_roton',num2str(rot_on)])
        end
    elseif fit_type == 3
        save(['MODcontributionPerc2FznMx_indivFreqNrollFit_linear_roton',num2str(rot_on),'.mat'],'stroke2freq','freq2stroke')
        if plot_on == 1
            mkdir(['qsModel_FnM_indivFreqNrollFit_linear_roton',num2str(rot_on)])
            cd(['qsModel_FnM_indivFreqNrollFit_linear_roton',num2str(rot_on)])
        end
    elseif fit_type == 4
        save(['MODcontributionPerc2FznMx_indivFreqNrollFit_power2_roton',num2str(rot_on),'.mat'],'stroke2freq','freq2stroke')
        if plot_on == 1
            mkdir(['qsModel_FnM_indivFreqNrollFit_power2_roton',num2str(rot_on)])
            cd(['qsModel_FnM_indivFreqNrollFit_power2_roton',num2str(rot_on)])
        end
    end
    
if plot_on == 1
    figure(1000)
        
    subplot(2,2,1)
    hold on
    plot(1,1,'sk','markerfacecolor','w','markersize',10)
    plot(1,1,'ok','markerfacecolor','w','markersize',10)
    plot(1,1,'dk','markerfacecolor','w','markersize',10)
    legend('undamaged','tip cut','TE cut')
    colormap(cmap_AdAi)
    caxis([Aratio_min Aratio_max])
    h = colorbar('location','northoutside'); 
    title(h,'A+')
    set(h,'xtick',Aratio_min:(Aratio_max-Aratio_min):Aratio_max)
    axis off

    saveas(gca,['S2vsS3vsAmpRatio.fig'])
    saveas(gca,['S2vsS3vsAmpRatio.png'])
    plot2svg(['S2vsS3vsAmpRatio.svg'])
    
    figure(1)
    saveas(gca,['FnM_WBmod_TEnTipClip_YnZflip.fig'])
    saveas(gca,['FnM_WBmod_TEnTipClip_YnZflip.png'])
    plot2svg(['FnM_WBmod_TEnTipClip_YnZflip.svg'])

    figure(4)
    saveas(gca,['Fz_Mx_WBmodCombies_freqDevRotStroke_TEnTipclip_YnZflip.fig'])
    saveas(gca,['Fz_Mx_WBmodCombies_freqDevRotStroke_TEnTipclip_YnZflip.png'])
    plot2svg(['Fz_Mx_WBmodCombies_freqDevRotStroke_TEnTipclip_YnZflip.svg'])

    figure(6)
    saveas(gca,['Fz_Mx_WBmodCombies_strokeRotDevRollFreq_TEnTipclip_YnZflip.fig'])
    saveas(gca,['Fz_Mx_WBmodCombies_strokeRotDevRollFreq_TEnTipclip_YnZflip.png'])
    plot2svg(['Fz_Mx_WBmodCombies_strokeRotDevRollFreq_TEnTipclip_YnZflip.svg'])

    cd ..
end    
    
    