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
% fit_type = 2   % spline999;
% fit_type = 3   % linear;
% fit_type = 4   % power2;
    
load('fit_type.mat')
fit_type = fit_type+1

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
    
    Mx_damaged_mean_all_TEnTip = Mx_damaged_mean_all_rolled;
    My_damaged_mean_all_TEnTip = My_damaged_mean_all;
    Mz_damaged_mean_all_TEnTip = Mz_damaged_mean_all;
    
    Fx_damaged_mean_steady_TEnTip = Fx_damaged_mean_steady;
    Fy_damaged_mean_steady_TEnTip = Fy_damaged_mean_steady;
    Fz_damaged_mean_steady_TEnTip = Fz_damaged_mean_steady;
    
    Mx_damaged_mean_steady_TEnTip = Mx_damaged_mean_steady;
    My_damaged_mean_steady_TEnTip = My_damaged_mean_steady;
    Mz_damaged_mean_steady_TEnTip = Mz_damaged_mean_steady;
    
    %% FREQ FIRST, ROLL LAST
    % freq
    Fx_damaged_mean_freqMOD_TEnTip = Fx_damaged_mean_freqMOD;
    Fy_damaged_mean_freqMOD_TEnTip = Fy_damaged_mean_freqMOD;
    Fz_damaged_mean_freqMOD_TEnTip = Fz_damaged_mean_freqMOD;
    
    Mx_damaged_mean_freqMOD_TEnTip = Mx_damaged_mean_freqMOD;
    My_damaged_mean_freqMOD_TEnTip = My_damaged_mean_freqMOD;
    Mz_damaged_mean_freqMOD_TEnTip = Mz_damaged_mean_freqMOD;
    
    %freq&rot
    Fx_damaged_mean_freqRotMOD_TEnTip = Fx_damaged_mean_freqRotMOD;
    Fy_damaged_mean_freqRotMOD_TEnTip = Fy_damaged_mean_freqRotMOD;
    Fz_damaged_mean_freqRotMOD_TEnTip = Fz_damaged_mean_freqRotMOD;
    
    Mx_damaged_mean_freqRotMOD_TEnTip = Mx_damaged_mean_freqRotMOD;
    My_damaged_mean_freqRotMOD_TEnTip = My_damaged_mean_freqRotMOD;
    Mz_damaged_mean_freqRotMOD_TEnTip = Mz_damaged_mean_freqRotMOD;
    
    %freq&rot&dev
    Fx_damaged_mean_freqRotDevMOD_TEnTip = Fx_damaged_mean_freqRotDevMOD;
    Fy_damaged_mean_freqRotDevMOD_TEnTip = Fy_damaged_mean_freqRotDevMOD;
    Fz_damaged_mean_freqRotDevMOD_TEnTip = Fz_damaged_mean_freqRotDevMOD;
    
    Mx_damaged_mean_freqRotDevMOD_TEnTip = Mx_damaged_mean_freqRotDevMOD;
    My_damaged_mean_freqRotDevMOD_TEnTip = My_damaged_mean_freqRotDevMOD;
    Mz_damaged_mean_freqRotDevMOD_TEnTip = Mz_damaged_mean_freqRotDevMOD;
    
    %freq&rot&dev&stroke
    Fx_damaged_mean_freqRotDevStrokeMOD_TEnTip = Fx_damaged_mean_all;
    Fy_damaged_mean_freqRotDevStrokeMOD_TEnTip = Fy_damaged_mean_all;
    Fz_damaged_mean_freqRotDevStrokeMOD_TEnTip = Fz_damaged_mean_all;
    
    Mx_damaged_mean_freqRotDevStrokeMOD_TEnTip = Mx_damaged_mean_all;
    My_damaged_mean_freqRotDevStrokeMOD_TEnTip = My_damaged_mean_all;
    Mz_damaged_mean_freqRotDevStrokeMOD_TEnTip = Mz_damaged_mean_all;
    
    %% ROLL FIRST, FREQ LAST
    % roll
    Fx_damaged_mean_rollMOD_TEnTip = Fx_damaged_mean_steady;
    Fy_damaged_mean_rollMOD_TEnTip = Fy_damaged_mean_steady_rolled;
    Fz_damaged_mean_rollMOD_TEnTip = Fz_damaged_mean_steady_rolled;
    
    Mx_damaged_mean_rollMOD_TEnTip = Mx_damaged_mean_steady_rolled;
    My_damaged_mean_rollMOD_TEnTip = My_damaged_mean_steady;
    Mz_damaged_mean_rollMOD_TEnTip = Mz_damaged_mean_steady;
    
    % roll&stroke
    Fx_damaged_mean_rollStrokeMOD_TEnTip = Fx_damaged_mean_strokeMOD;
    Fy_damaged_mean_rollStrokeMOD_TEnTip = Fy_damaged_mean_strokeMOD_rolled;
    Fz_damaged_mean_rollStrokeMOD_TEnTip = Fz_damaged_mean_strokeMOD_rolled;
    
    Mx_damaged_mean_rollStrokeMOD_TEnTip = Mx_damaged_mean_strokeMOD_rolled;
    My_damaged_mean_rollStrokeMOD_TEnTip = My_damaged_mean_strokeMOD;
    Mz_damaged_mean_rollStrokeMOD_TEnTip = Mz_damaged_mean_strokeMOD;
    
    % roll&stroke&dev
    Fx_damaged_mean_rollStrokeDevMOD_TEnTip = Fx_damaged_mean_strokeDevMOD;
    Fy_damaged_mean_rollStrokeDevMOD_TEnTip = Fy_damaged_mean_strokeDevMOD_rolled;
    Fz_damaged_mean_rollStrokeDevMOD_TEnTip = Fz_damaged_mean_strokeDevMOD_rolled;
    
    Mx_damaged_mean_rollStrokeDevMOD_TEnTip = Mx_damaged_mean_strokeDevMOD_rolled;
    My_damaged_mean_rollStrokeDevMOD_TEnTip = My_damaged_mean_strokeDevMOD;
    Mz_damaged_mean_rollStrokeDevMOD_TEnTip = Mz_damaged_mean_strokeDevMOD;
    
    % roll&stroke&dev&rot
    Fx_damaged_mean_rollStrokeDevRotMOD_TEnTip = Fx_damaged_mean_strokeDevRotMOD;
    Fy_damaged_mean_rollStrokeDevRotMOD_TEnTip = Fy_damaged_mean_strokeDevRotMOD_rolled;
    Fz_damaged_mean_rollStrokeDevRotMOD_TEnTip = Fz_damaged_mean_strokeDevRotMOD_rolled;
    
    Mx_damaged_mean_rollStrokeDevRotMOD_TEnTip = Mx_damaged_mean_strokeDevRotMOD_rolled;
    My_damaged_mean_rollStrokeDevRotMOD_TEnTip = My_damaged_mean_strokeDevRotMOD;
    Mz_damaged_mean_rollStrokeDevRotMOD_TEnTip = Mz_damaged_mean_strokeDevRotMOD;
    
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
    
    plot(S3ratios,Mx_damaged_mean_all_rolled,'ok','markersize',10,'markerfacecolor','r')
    plot(S3ratios,-My_damaged_mean_all--My_damaged_mean_all(S3ratios==1),'ok','markersize',10,'markerfacecolor',[0 .5 0])
    plot(S3ratios,-Mz_damaged_mean_all,'ok','markersize',10,'markerfacecolor',[0 0 .5])
    
    %% compare MODs -Fz & Mx COMBIMOD FREQ FIRST & ROLL LAST
    figure(4)
    subplot(2,2,1)
    hold on
    plot(S2ratios,-Fz_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','w')
    plot(S2ratios,-Fz_damaged_mean_freqMOD,'ok','markersize',10,'markerfacecolor','y')
    plot(S2ratios,-Fz_damaged_mean_freqRotMOD,'ok','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S2ratios,-Fz_damaged_mean_freqRotDevMOD,'ok','markersize',10,'markerfacecolor','r')
    plot(S2ratios,-Fz_damaged_mean_all,'ok','markersize',10,'markerfacecolor',[.5 0 0])
    plot(S2ratios,-Fz_damaged_mean_all_rolled,'ok','markersize',10,'markerfacecolor','k')
    
    subplot(2,2,2)
    hold on
    plot(S3ratios,Mx_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','w')
    plot(S3ratios,Mx_damaged_mean_freqMOD,'ok','markersize',10,'markerfacecolor','y')
    plot(S3ratios,Mx_damaged_mean_freqRotMOD,'ok','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S3ratios,Mx_damaged_mean_freqRotDevMOD,'ok','markersize',10,'markerfacecolor','r')
    plot(S3ratios,Mx_damaged_mean_all,'ok','markersize',10,'markerfacecolor',[.5 0 0])
    plot(S3ratios,Mx_damaged_mean_all_rolled,'ok','markersize',10,'markerfacecolor','k')
    
    %% compare MODs -Fz & Mx COMBIMOD ROLL FIRST & FREQ LAST
    figure(6)
    subplot(2,2,1)
    hold on
    plot(S2ratios,-Fz_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','w')
    plot(S2ratios,-Fz_damaged_mean_steady_rolled,'ok','markersize',10,'markerfacecolor','c')
    plot(S2ratios,-Fz_damaged_mean_strokeMOD_rolled,'ok','markersize',10,'markerfacecolor',[0 .5 1])
    plot(S2ratios,-Fz_damaged_mean_strokeDevMOD_rolled,'ok','markersize',10,'markerfacecolor','b')
    plot(S2ratios,-Fz_damaged_mean_strokeDevRotMOD_rolled,'ok','markersize',10,'markerfacecolor',[0 0 .5])
    plot(S2ratios,-Fz_damaged_mean_all_rolled,'ok','markersize',10,'markerfacecolor','k')
    
    subplot(2,2,2)
    hold on
    plot(S2ratios,Mx_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','w')
    plot(S2ratios,Mx_damaged_mean_steady_rolled,'ok','markersize',10,'markerfacecolor','c')
    plot(S2ratios,Mx_damaged_mean_strokeMOD_rolled,'ok','markersize',10,'markerfacecolor',[0 .5 1])
    plot(S2ratios,Mx_damaged_mean_strokeDevMOD_rolled,'ok','markersize',10,'markerfacecolor','b')
    plot(S2ratios,Mx_damaged_mean_strokeDevRotMOD_rolled,'ok','markersize',10,'markerfacecolor',[0 0 .5])
    plot(S2ratios,Mx_damaged_mean_all_rolled,'ok','markersize',10,'markerfacecolor','k')
    
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
    
    Mx_damaged_mean_all_TEnTip = [Mx_damaged_mean_all_TEnTip Mx_damaged_mean_all_rolled];
    My_damaged_mean_all_TEnTip = [My_damaged_mean_all_TEnTip My_damaged_mean_all];
    Mz_damaged_mean_all_TEnTip = [Mz_damaged_mean_all_TEnTip Mz_damaged_mean_all];
    
    Fx_damaged_mean_steady_TEnTip = [Fx_damaged_mean_steady_TEnTip Fx_damaged_mean_steady];
    Fy_damaged_mean_steady_TEnTip = [Fy_damaged_mean_steady_TEnTip Fy_damaged_mean_steady];
    Fz_damaged_mean_steady_TEnTip = [Fz_damaged_mean_steady_TEnTip Fz_damaged_mean_steady];
    
    Mx_damaged_mean_steady_TEnTip = [Mx_damaged_mean_steady_TEnTip Mx_damaged_mean_steady];
    My_damaged_mean_steady_TEnTip = [My_damaged_mean_steady_TEnTip My_damaged_mean_steady];
    Mz_damaged_mean_steady_TEnTip = [Mz_damaged_mean_steady_TEnTip Mz_damaged_mean_steady];
    
    %% FREQ FIRST, ROLL LAST
    % freq
    Fx_damaged_mean_freqMOD_TEnTip = [Fx_damaged_mean_freqMOD_TEnTip Fx_damaged_mean_freqMOD];
    Fy_damaged_mean_freqMOD_TEnTip = [Fy_damaged_mean_freqMOD_TEnTip Fy_damaged_mean_freqMOD];
    Fz_damaged_mean_freqMOD_TEnTip = [Fz_damaged_mean_freqMOD_TEnTip Fz_damaged_mean_freqMOD];
    
    Mx_damaged_mean_freqMOD_TEnTip = [Mx_damaged_mean_freqMOD_TEnTip Mx_damaged_mean_freqMOD];
    My_damaged_mean_freqMOD_TEnTip = [My_damaged_mean_freqMOD_TEnTip My_damaged_mean_freqMOD];
    Mz_damaged_mean_freqMOD_TEnTip = [Mz_damaged_mean_freqMOD_TEnTip Mz_damaged_mean_freqMOD];
    
    % freqROT
    Fx_damaged_mean_freqRotMOD_TEnTip = [Fx_damaged_mean_freqRotMOD_TEnTip Fx_damaged_mean_freqRotMOD];
    Fy_damaged_mean_freqRotMOD_TEnTip = [Fy_damaged_mean_freqRotMOD_TEnTip Fy_damaged_mean_freqRotMOD];
    Fz_damaged_mean_freqRotMOD_TEnTip = [Fz_damaged_mean_freqRotMOD_TEnTip Fz_damaged_mean_freqRotMOD];
    
    Mx_damaged_mean_freqRotMOD_TEnTip = [Mx_damaged_mean_freqRotMOD_TEnTip Mx_damaged_mean_freqRotMOD];
    My_damaged_mean_freqRotMOD_TEnTip = [My_damaged_mean_freqRotMOD_TEnTip My_damaged_mean_freqRotMOD];
    Mz_damaged_mean_freqRotMOD_TEnTip = [Mz_damaged_mean_freqRotMOD_TEnTip Mz_damaged_mean_freqRotMOD];
    
    % freqROTDEV
    Fx_damaged_mean_freqRotDevMOD_TEnTip = [Fx_damaged_mean_freqRotDevMOD_TEnTip Fx_damaged_mean_freqRotDevMOD];
    Fy_damaged_mean_freqRotDevMOD_TEnTip = [Fy_damaged_mean_freqRotDevMOD_TEnTip Fy_damaged_mean_freqRotDevMOD];
    Fz_damaged_mean_freqRotDevMOD_TEnTip = [Fz_damaged_mean_freqRotDevMOD_TEnTip Fz_damaged_mean_freqRotDevMOD];
    
    Mx_damaged_mean_freqRotDevMOD_TEnTip = [Mx_damaged_mean_freqRotDevMOD_TEnTip Mx_damaged_mean_freqRotDevMOD];
    My_damaged_mean_freqRotDevMOD_TEnTip = [My_damaged_mean_freqRotDevMOD_TEnTip My_damaged_mean_freqRotDevMOD];
    Mz_damaged_mean_freqRotDevMOD_TEnTip = [Mz_damaged_mean_freqRotDevMOD_TEnTip Mz_damaged_mean_freqRotDevMOD];
    
    % freqROTDEVSTR
    Fx_damaged_mean_freqRotDevStrokeMOD_TEnTip = [Fx_damaged_mean_freqRotDevStrokeMOD_TEnTip Fx_damaged_mean_all];
    Fy_damaged_mean_freqRotDevStrokeMOD_TEnTip = [Fy_damaged_mean_freqRotDevStrokeMOD_TEnTip Fy_damaged_mean_all];
    Fz_damaged_mean_freqRotDevStrokeMOD_TEnTip = [Fz_damaged_mean_freqRotDevStrokeMOD_TEnTip Fz_damaged_mean_all];
    
    Mx_damaged_mean_freqRotDevStrokeMOD_TEnTip = [Mx_damaged_mean_freqRotDevStrokeMOD_TEnTip Mx_damaged_mean_all];
    My_damaged_mean_freqRotDevStrokeMOD_TEnTip = [My_damaged_mean_freqRotDevStrokeMOD_TEnTip My_damaged_mean_all];
    Mz_damaged_mean_freqRotDevStrokeMOD_TEnTip = [Mz_damaged_mean_freqRotDevStrokeMOD_TEnTip Mz_damaged_mean_all];
    
    %% ROLL FIRST, FREQ LAST
    % roll
    Fx_damaged_mean_rollMOD_TEnTip = [Fx_damaged_mean_rollMOD_TEnTip Fx_damaged_mean_steady];
    Fy_damaged_mean_rollMOD_TEnTip = [Fy_damaged_mean_rollMOD_TEnTip Fy_damaged_mean_steady_rolled];
    Fz_damaged_mean_rollMOD_TEnTip = [Fz_damaged_mean_rollMOD_TEnTip Fz_damaged_mean_steady_rolled];
    
    Mx_damaged_mean_rollMOD_TEnTip = [Mx_damaged_mean_rollMOD_TEnTip Mx_damaged_mean_steady_rolled];
    My_damaged_mean_rollMOD_TEnTip = [My_damaged_mean_rollMOD_TEnTip My_damaged_mean_steady];
    Mz_damaged_mean_rollMOD_TEnTip = [Mz_damaged_mean_rollMOD_TEnTip Mz_damaged_mean_steady];
    
    % roll&str
    Fx_damaged_mean_rollStrokeMOD_TEnTip = [Fx_damaged_mean_rollStrokeMOD_TEnTip Fx_damaged_mean_strokeMOD];
    Fy_damaged_mean_rollStrokeMOD_TEnTip = [Fy_damaged_mean_rollStrokeMOD_TEnTip Fy_damaged_mean_strokeMOD_rolled];
    Fz_damaged_mean_rollStrokeMOD_TEnTip = [Fz_damaged_mean_rollStrokeMOD_TEnTip Fz_damaged_mean_strokeMOD_rolled];
    
    Mx_damaged_mean_rollStrokeMOD_TEnTip = [Mx_damaged_mean_rollStrokeMOD_TEnTip Mx_damaged_mean_strokeMOD_rolled];
    My_damaged_mean_rollStrokeMOD_TEnTip = [My_damaged_mean_rollStrokeMOD_TEnTip My_damaged_mean_strokeMOD];
    Mz_damaged_mean_rollStrokeMOD_TEnTip = [Mz_damaged_mean_rollStrokeMOD_TEnTip Mz_damaged_mean_strokeMOD];

    % roll&str&dev
    Fx_damaged_mean_rollStrokeDevMOD_TEnTip = [Fx_damaged_mean_rollStrokeDevMOD_TEnTip Fx_damaged_mean_strokeDevMOD];
    Fy_damaged_mean_rollStrokeDevMOD_TEnTip = [Fy_damaged_mean_rollStrokeDevMOD_TEnTip Fy_damaged_mean_strokeDevMOD_rolled];
    Fz_damaged_mean_rollStrokeDevMOD_TEnTip = [Fz_damaged_mean_rollStrokeDevMOD_TEnTip Fz_damaged_mean_strokeDevMOD_rolled];
    
    Mx_damaged_mean_rollStrokeDevMOD_TEnTip = [Mx_damaged_mean_rollStrokeDevMOD_TEnTip Mx_damaged_mean_strokeDevMOD_rolled];
    My_damaged_mean_rollStrokeDevMOD_TEnTip = [My_damaged_mean_rollStrokeDevMOD_TEnTip My_damaged_mean_strokeDevMOD];
    Mz_damaged_mean_rollStrokeDevMOD_TEnTip = [Mz_damaged_mean_rollStrokeDevMOD_TEnTip Mz_damaged_mean_strokeDevMOD];

    % roll&str&dev&rot
    Fx_damaged_mean_rollStrokeDevRotMOD_TEnTip = [Fx_damaged_mean_rollStrokeDevRotMOD_TEnTip Fx_damaged_mean_strokeDevRotMOD];
    Fy_damaged_mean_rollStrokeDevRotMOD_TEnTip = [Fy_damaged_mean_rollStrokeDevRotMOD_TEnTip Fy_damaged_mean_strokeDevRotMOD_rolled];
    Fz_damaged_mean_rollStrokeDevRotMOD_TEnTip = [Fz_damaged_mean_rollStrokeDevRotMOD_TEnTip Fz_damaged_mean_strokeDevRotMOD_rolled];
    
    Mx_damaged_mean_rollStrokeDevRotMOD_TEnTip = [Mx_damaged_mean_rollStrokeDevRotMOD_TEnTip Mx_damaged_mean_strokeDevRotMOD_rolled];
    My_damaged_mean_rollStrokeDevRotMOD_TEnTip = [My_damaged_mean_rollStrokeDevRotMOD_TEnTip My_damaged_mean_strokeDevRotMOD];
    Mz_damaged_mean_rollStrokeDevRotMOD_TEnTip = [Mz_damaged_mean_rollStrokeDevRotMOD_TEnTip Mz_damaged_mean_strokeDevRotMOD];

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
    
    plot(S3ratios,Mx_damaged_mean_all_rolled,'dk','markersize',10,'markerfacecolor','r')
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
    
    %% compare MODs -Fz & Mx COMBIMOD FREQ FIRST & ROLL LAST
    figure(4)
    subplot(2,2,1)
    hold on
    plot(S2ratios,-Fz_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','w')
    plot(S2ratios,-Fz_damaged_mean_freqMOD,'dk','markersize',10,'markerfacecolor','y')
    plot(S2ratios,-Fz_damaged_mean_freqRotMOD,'dk','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S2ratios,-Fz_damaged_mean_freqRotDevMOD,'dk','markersize',10,'markerfacecolor','r')
    plot(S2ratios,-Fz_damaged_mean_all,'dk','markersize',10,'markerfacecolor',[.5 0 0])
    plot(S2ratios,-Fz_damaged_mean_all_rolled,'dk','markersize',10,'markerfacecolor','k')
    
    subplot(2,2,2)
    hold on
    plot(S3ratios,Mx_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','w')
    plot(S3ratios,Mx_damaged_mean_freqMOD,'dk','markersize',10,'markerfacecolor','y')
    plot(S3ratios,Mx_damaged_mean_freqRotMOD,'dk','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S3ratios,Mx_damaged_mean_freqRotDevMOD,'dk','markersize',10,'markerfacecolor','r')
    plot(S3ratios,Mx_damaged_mean_all,'dk','markersize',10,'markerfacecolor',[.5 0 0])
    plot(S3ratios,Mx_damaged_mean_all_rolled,'dk','markersize',10,'markerfacecolor','k')

    subplot(2,2,1)
    xlabel('S2 ratio')
    ylabel('normalized vertical force -Fz/mg')
    axis([0.5 1 .5 1.5])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-1.5:.5:1.5) 

    subplot(2,2,2)
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
    legend('none','frequency','freq&rot','freq&rot&dev','freq&rot&dev&str','all','location','E')
    colormap(cmap_MODs_red)
    colorbar
    axis off

    %% compare MODs -Fz & Mx COMBIMOD ROLL FIRST & FREQ LAST
    figure(6)
    subplot(2,2,1)
    hold on
    plot(S2ratios,-Fz_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','w')
    plot(S2ratios,-Fz_damaged_mean_steady_rolled,'dk','markersize',10,'markerfacecolor','c')
    plot(S2ratios,-Fz_damaged_mean_strokeMOD_rolled,'dk','markersize',10,'markerfacecolor',[0 .5 1])
    plot(S2ratios,-Fz_damaged_mean_strokeDevMOD_rolled,'dk','markersize',10,'markerfacecolor','b')
    plot(S2ratios,-Fz_damaged_mean_strokeDevRotMOD_rolled,'dk','markersize',10,'markerfacecolor',[0 0 .5])
    plot(S2ratios,-Fz_damaged_mean_all_rolled,'dk','markersize',10,'markerfacecolor','k')
    
    subplot(2,2,2)
    hold on
    plot(S2ratios,Mx_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','w')
    plot(S2ratios,Mx_damaged_mean_steady_rolled,'dk','markersize',10,'markerfacecolor','c')
    plot(S2ratios,Mx_damaged_mean_strokeMOD_rolled,'dk','markersize',10,'markerfacecolor',[0 .5 1])
    plot(S2ratios,Mx_damaged_mean_strokeDevMOD_rolled,'dk','markersize',10,'markerfacecolor','b')
    plot(S2ratios,Mx_damaged_mean_strokeDevRotMOD_rolled,'dk','markersize',10,'markerfacecolor',[0 0 .5])
    plot(S2ratios,Mx_damaged_mean_all_rolled,'dk','markersize',10,'markerfacecolor','k')

    subplot(2,2,1)
    xlabel('S2 ratio')
    ylabel('normalized vertical force -Fz/mg')
    axis([0.5 1 .5 1.5])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-1.5:.5:1.5) 

    subplot(2,2,2)
    xlabel('S3 ratio')
    ylabel('normalized roll torque Tx/mgl')
    axis([0.5 1 -.05 .15])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-.05:.05:.15)    

    subplot(2,2,3)
    hold on
    plot(0,0,'ok','markersize',10,'markerfacecolor','w')
    plot(0,0,'ok','markersize',10,'markerfacecolor','c')
    plot(0,0,'ok','markersize',10,'markerfacecolor',[0 .5 1])
    plot(0,0,'ok','markersize',10,'markerfacecolor','b')
    plot(0,0,'ok','markersize',10,'markerfacecolor',[0 0 .5])
    plot(0,0,'ok','markersize',10,'markerfacecolor','k')
    legend('none','roll','roll&stroke','roll&str&dev','roll&str&dev&rot','all','location','E')
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

    % adjusted
    plot([min(S2ratios_TEnTip):.05:max(S2ratios_TEnTip) 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',2),[min(S2ratios_TEnTip):.05:max(S2ratios_TEnTip) 1]),'color',[0 0 .5],'linewidth',2)

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

    %% plot trendlines FREQ FIRST & ROLL LAST
    figure(4)
    subplot(2,2,1)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[.5 .5 .5],'linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color','y','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqRotMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[1 .5 0],'linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqRotDevMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color','r','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqRotDevStrokeMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[.5 0 0],'linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color','k','linewidth',2)

%     % adjusted
%     plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',2),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[.5 .5 .5],'linewidth',2)
%     plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqMOD_TEnTip',2),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color','b','linewidth',2)

    subplot(2,2,2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_steady_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color',[.5 .5 .5],'linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color','y','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqRotMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color',[1 .5 0],'linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqRotDevMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color','r','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqRotDevStrokeMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color',[.5 0 0],'linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_all_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color','k','linewidth',2)

    %% plot trendlines ROLL FIRST & FREQ LAST
    figure(6)
    subplot(2,2,1)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[.5 .5 .5],'linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_rollMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color','c','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_rollStrokeMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[0 .5 1],'linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_rollStrokeDevMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color','b','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_rollStrokeDevRotMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[0 0 .5],'linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color','k','linewidth',2)

    % adjusted
    plot([min(S2ratios_TEnTip):.05:max(S2ratios_TEnTip) 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',2),[min(S2ratios_TEnTip):.05:max(S2ratios_TEnTip) 1]),'color','k','linewidth',2)

    subplot(2,2,2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_steady_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color',[.5 .5 .5],'linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_rollMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color','c','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_rollStrokeMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color',[0 .5 1],'linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_rollStrokeDevMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color','b','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_rollStrokeDevRotMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color',[0 0 .5],'linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_all_TEnTip',1),[min(S3ratios_TEnTip) max(S2ratios_TEnTip)]),'color','k','linewidth',2)

    %% percentage of weight support and roll equilibrium FREQ FIRST & ROLL LAST
    Fz_steady = -Fz_damaged_mean_steady_TEnTip(1);
    Fz_none = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_freq = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_freqRot = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqRotMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_freqRotDev = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqRotDevMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_freqRotDevStroke = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqRotDevStrokeMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_freqRotDevStrokeRoll = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));

    freq2roll.Fz_freq_perc_frds    = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_freq-Fz_none)
    freq2roll.Fz_rot_perc_frds     = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_freqRot-Fz_freq)
    freq2roll.Fz_dev_perc_frds     = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_freqRotDev-Fz_freqRot)
    freq2roll.Fz_stroke_perc_frds  = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_freqRotDevStroke-Fz_freqRotDev)
    freq2roll.Fz_roll_perc_frds    = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_freqRotDevStrokeRoll-Fz_freqRotDevStroke)
    freq2roll.Fz_all_perc_frds     = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_freqRotDevStrokeRoll-Fz_none)
    freq2roll.Fz_sum_perc_frds     = freq2roll.Fz_freq_perc_frds + freq2roll.Fz_rot_perc_frds + freq2roll.Fz_dev_perc_frds + freq2roll.Fz_stroke_perc_frds + freq2roll.Fz_roll_perc_frds
    
    Mx_steady = Mx_damaged_mean_steady_TEnTip(1);
    Mx_none = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_freq = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_freqRot = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqRotMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_freqRotDev = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqRotDevMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_freqRotDevStroke = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqRotDevStrokeMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_freqRotDevStrokeRoll = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));

    freq2roll.Mx_freq_perc_frds    = 100 / (-Mx_none) * (Mx_freq-Mx_none)
    freq2roll.Mx_rot_perc_frds     = 100 / (-Mx_none) * (Mx_freqRot-Mx_freq)
    freq2roll.Mx_dev_perc_frds     = 100 / (-Mx_none) * (Mx_freqRotDev-Mx_freqRot)
    freq2roll.Mx_stroke_perc_frds  = 100 / (-Mx_none) * (Mx_freqRotDevStroke-Mx_freqRotDev)
    freq2roll.Mx_roll_perc_frds    = 100 / (-Mx_none) * (Mx_freqRotDevStrokeRoll-Mx_freqRotDevStroke)
    freq2roll.Mx_all_perc_frds     = 100 / (-Mx_none) * (Mx_freqRotDevStrokeRoll-Mx_none)
    freq2roll.Mx_sum_perc_frds     = freq2roll.Mx_freq_perc_frds + freq2roll.Mx_rot_perc_frds + freq2roll.Mx_dev_perc_frds + freq2roll.Mx_stroke_perc_frds + freq2roll.Mx_roll_perc_frds
        
    %% percentage of weight support and roll equilibrium ROLL FIRST & FREQ LAST
    Fz_steady = -Fz_damaged_mean_steady_TEnTip(1);
    Fz_none = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_roll = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_rollMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_rollStroke = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_rollStrokeMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_rollStrokeDev = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_rollStrokeDevMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_rollStrokeDevRot = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_rollStrokeDevRotMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_rollStrokeDevRotFreq = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));

    roll2freq.Fz_roll_perc_sdrf    = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_roll-Fz_none)
    roll2freq.Fz_stroke_perc_sdrf  = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_rollStroke-Fz_none)
    roll2freq.Fz_dev_perc_sdrf     = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_rollStrokeDev-Fz_rollStroke)
    roll2freq.Fz_rot_perc_sdrf     = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_rollStrokeDevRot-Fz_rollStrokeDev)
    roll2freq.Fz_freq_perc_sdrf    = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_rollStrokeDevRotFreq-Fz_rollStrokeDevRot)
    roll2freq.Fz_all_perc_sdrf     = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_rollStrokeDevRotFreq-Fz_none)
    roll2freq.Fz_sum_perc_sdrf     = roll2freq.Fz_freq_perc_sdrf + roll2freq.Fz_stroke_perc_sdrf + roll2freq.Fz_dev_perc_sdrf + roll2freq.Fz_rot_perc_sdrf + roll2freq.Fz_roll_perc_sdrf
    
    Mx_steady = Mx_damaged_mean_steady_TEnTip(1);
    Mx_none = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_roll = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_rollMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_rollStroke = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_rollStrokeMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_rollStrokeDev = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_rollStrokeDevMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_rollStrokeDevRot = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_rollStrokeDevRotMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_rollStrokeDevRotFreq = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));

    roll2freq.Mx_roll_perc_sdrf    = 100 / (-Mx_none) * (Mx_roll-Mx_none)
    roll2freq.Mx_stroke_perc_sdrf  = 100 / (-Mx_none) * (Mx_rollStroke-Mx_roll)
    roll2freq.Mx_dev_perc_sdrf     = 100 / (-Mx_none) * (Mx_rollStrokeDev-Mx_rollStroke)
    roll2freq.Mx_rot_perc_sdrf     = 100 / (-Mx_none) * (Mx_rollStrokeDevRot-Mx_rollStrokeDev)
    roll2freq.Mx_freq_perc_sdrf    = 100 / (-Mx_none) * (Mx_rollStrokeDevRotFreq-Mx_rollStrokeDevRot)
    roll2freq.Mx_all_perc_sdrf     = 100 / (-Mx_none) * (Mx_rollStrokeDevRotFreq-Mx_none)
    roll2freq.Mx_sum_perc_sdrf     = roll2freq.Mx_freq_perc_sdrf + roll2freq.Mx_stroke_perc_sdrf + roll2freq.Mx_dev_perc_sdrf + roll2freq.Mx_rot_perc_sdrf + roll2freq.Mx_roll_perc_sdrf
    
    %% save plots
    %% save data
    if fit_type == 1
        save(['MODcontributionPerc2FznMx_indivFreqNrollFit_Asym10_roton',num2str(rot_on),'.mat'],'roll2freq','freq2roll')
        if plot_on == 1
            mkdir(['qsModel_FnM_indivFreqNrollFit_Asym10_roton',num2str(rot_on)])
            cd(['qsModel_FnM_indivFreqNrollFit_Asym10_roton',num2str(rot_on)])
        end
    elseif fit_type == 2
        save(['MODcontributionPerc2FznMx_indivFreqNrollFit_spline999_roton',num2str(rot_on),'.mat'],'roll2freq','freq2roll')
        if plot_on == 1
            mkdir(['qsModel_FnM_indivFreqNrollFit_spline999_roton',num2str(rot_on)])
            cd(['qsModel_FnM_indivFreqNrollFit_spline999_roton',num2str(rot_on)])
        end
    elseif fit_type == 3
        save(['MODcontributionPerc2FznMx_indivFreqNrollFit_linear_roton',num2str(rot_on),'.mat'],'roll2freq','freq2roll')
        if plot_on == 1
            mkdir(['qsModel_FnM_indivFreqNrollFit_linear_roton',num2str(rot_on)])
            cd(['qsModel_FnM_indivFreqNrollFit_linear_roton',num2str(rot_on)])
        end
    elseif fit_type == 4
        save(['MODcontributionPerc2FznMx_indivFreqNrollFit_power2_roton',num2str(rot_on),'.mat'],'roll2freq','freq2roll')
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
    saveas(gca,['Fz_Mx_WBmodCombies_freqRotDevStroke_TEnTipclip_YnZflip.fig'])
    saveas(gca,['Fz_Mx_WBmodCombies_freqRotDevStroke_TEnTipclip_YnZflip.png'])
    plot2svg(['Fz_Mx_WBmodCombies_freqRotDevStroke_TEnTipclip_YnZflip.svg'])

    figure(6)
    saveas(gca,['Fz_Mx_WBmodCombies_strokeFirst_TEnTipclip_YnZflip.fig'])
    saveas(gca,['Fz_Mx_WBmodCombies_strokeFirst_TEnTipclip_YnZflip.png'])
    plot2svg(['Fz_Mx_WBmodCombies_strokeFirst_TEnTipclip_YnZflip.svg'])

    cd ..
end    
    
    