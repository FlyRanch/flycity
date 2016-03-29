clear;
clc;
close all
warning off

cmap_MODs = [1 1 1; 1 1 0; 1 .5 0; 1 0 0; .5 0 0]; 

freq_asymFitNr = 10;
peakloc = 0;

        %% plot TIP cut wing data
    load(['allMODs_TipClip_freqAsym',num2str(freq_asymFitNr),'.mat'])

    S2ratios_TEnTip = S2ratios;
    S3ratios_TEnTip = S3ratios;
    S2S3AmpRatioFuncs_TEnTip = S2S3AmpRatioFuncs;
    
    Fx_damaged_mean_all_TEnTip = Fx_damaged_mean_all;
    Fy_damaged_mean_all_TEnTip = Fy_damaged_mean_all;
    Fz_damaged_mean_all_TEnTip = Fz_damaged_mean_all;
    
    Mx_damaged_mean_all_TEnTip = Mx_damaged_mean_all;
    My_damaged_mean_all_TEnTip = My_damaged_mean_all;
    Mz_damaged_mean_all_TEnTip = Mz_damaged_mean_all;
    
    Fx_damaged_mean_steady_TEnTip = Fx_damaged_mean_steady;
    Fy_damaged_mean_steady_TEnTip = Fy_damaged_mean_steady;
    Fz_damaged_mean_steady_TEnTip = Fz_damaged_mean_steady;
    
    Mx_damaged_mean_steady_TEnTip = Mx_damaged_mean_steady;
    My_damaged_mean_steady_TEnTip = My_damaged_mean_steady;
    Mz_damaged_mean_steady_TEnTip = Mz_damaged_mean_steady;
    
    Fx_damaged_mean_freqMOD_TEnTip = Fx_damaged_mean_freqMOD;
    Fy_damaged_mean_freqMOD_TEnTip = Fy_damaged_mean_freqMOD;
    Fz_damaged_mean_freqMOD_TEnTip = Fz_damaged_mean_freqMOD;
    
    Mx_damaged_mean_freqMOD_TEnTip = Mx_damaged_mean_freqMOD;
    My_damaged_mean_freqMOD_TEnTip = My_damaged_mean_freqMOD;
    Mz_damaged_mean_freqMOD_TEnTip = Mz_damaged_mean_freqMOD;
    
    Fx_damaged_mean_strokeMOD_TEnTip = Fx_damaged_mean_strokeMOD;
    Fy_damaged_mean_strokeMOD_TEnTip = Fy_damaged_mean_strokeMOD;
    Fz_damaged_mean_strokeMOD_TEnTip = Fz_damaged_mean_strokeMOD;
    
    Mx_damaged_mean_strokeMOD_TEnTip = Mx_damaged_mean_strokeMOD;
    My_damaged_mean_strokeMOD_TEnTip = My_damaged_mean_strokeMOD;
    Mz_damaged_mean_strokeMOD_TEnTip = Mz_damaged_mean_strokeMOD;
    
    Fx_damaged_mean_devMOD_TEnTip = Fx_damaged_mean_devMOD;
    Fy_damaged_mean_devMOD_TEnTip = Fy_damaged_mean_devMOD;
    Fz_damaged_mean_devMOD_TEnTip = Fz_damaged_mean_devMOD;
    
    Mx_damaged_mean_devMOD_TEnTip = Mx_damaged_mean_devMOD;
    My_damaged_mean_devMOD_TEnTip = My_damaged_mean_devMOD;
    Mz_damaged_mean_devMOD_TEnTip = Mz_damaged_mean_devMOD;
    
    Fx_damaged_mean_rotMOD_TEnTip = Fx_damaged_mean_rotMOD;
    Fy_damaged_mean_rotMOD_TEnTip = Fy_damaged_mean_rotMOD;
    Fz_damaged_mean_rotMOD_TEnTip = Fz_damaged_mean_rotMOD;
    
    Mx_damaged_mean_rotMOD_TEnTip = Mx_damaged_mean_rotMOD;
    My_damaged_mean_rotMOD_TEnTip = My_damaged_mean_rotMOD;
    Mz_damaged_mean_rotMOD_TEnTip = Mz_damaged_mean_rotMOD;
    
    Fx_damaged_mean_freqStrokeMOD_TEnTip = Fx_damaged_mean_freqStrokeMOD;
    Fy_damaged_mean_freqStrokeMOD_TEnTip = Fy_damaged_mean_freqStrokeMOD;
    Fz_damaged_mean_freqStrokeMOD_TEnTip = Fz_damaged_mean_freqStrokeMOD;
    
    Mx_damaged_mean_freqStrokeMOD_TEnTip = Mx_damaged_mean_freqStrokeMOD;
    My_damaged_mean_freqStrokeMOD_TEnTip = My_damaged_mean_freqStrokeMOD;
    Mz_damaged_mean_freqStrokeMOD_TEnTip = Mz_damaged_mean_freqStrokeMOD;
    
    Fx_damaged_mean_freqStrokeDevMOD_TEnTip = Fx_damaged_mean_freqStrokeDevMOD;
    Fy_damaged_mean_freqStrokeDevMOD_TEnTip = Fy_damaged_mean_freqStrokeDevMOD;
    Fz_damaged_mean_freqStrokeDevMOD_TEnTip = Fz_damaged_mean_freqStrokeDevMOD;
    
    Mx_damaged_mean_freqStrokeDevMOD_TEnTip = Mx_damaged_mean_freqStrokeDevMOD;
    My_damaged_mean_freqStrokeDevMOD_TEnTip = My_damaged_mean_freqStrokeDevMOD;
    Mz_damaged_mean_freqStrokeDevMOD_TEnTip = Mz_damaged_mean_freqStrokeDevMOD;
    
    Fx_damaged_mean_strokeDevMOD_TEnTip = Fx_damaged_mean_strokeDevMOD;
    Fy_damaged_mean_strokeDevMOD_TEnTip = Fy_damaged_mean_strokeDevMOD;
    Fz_damaged_mean_strokeDevMOD_TEnTip = Fz_damaged_mean_strokeDevMOD;
    
    Mx_damaged_mean_strokeDevMOD_TEnTip = Mx_damaged_mean_strokeDevMOD;
    My_damaged_mean_strokeDevMOD_TEnTip = My_damaged_mean_strokeDevMOD;
    Mz_damaged_mean_strokeDevMOD_TEnTip = Mz_damaged_mean_strokeDevMOD;
    
    Fx_damaged_mean_strokeDevRotMOD_TEnTip = Fx_damaged_mean_strokeDevRotMOD;
    Fy_damaged_mean_strokeDevRotMOD_TEnTip = Fy_damaged_mean_strokeDevRotMOD;
    Fz_damaged_mean_strokeDevRotMOD_TEnTip = Fz_damaged_mean_strokeDevRotMOD;
    
    Mx_damaged_mean_strokeDevRotMOD_TEnTip = Mx_damaged_mean_strokeDevRotMOD;
    My_damaged_mean_strokeDevRotMOD_TEnTip = My_damaged_mean_strokeDevRotMOD;
    Mz_damaged_mean_strokeDevRotMOD_TEnTip = Mz_damaged_mean_strokeDevRotMOD;
    
    % ALL MODs VS NO MODs
    figure(1)
    subplot(2,2,1)
    hold on
    plot(S2ratios,Fx_damaged_mean_all,'ok','markersize',10,'markerfacecolor','b')
    plot(S2ratios,-Fy_damaged_mean_all,'ok','markersize',10,'markerfacecolor','r')
    plot(S2ratios,-Fz_damaged_mean_all,'ok','markersize',10,'markerfacecolor','g')
    
    plot(S2ratios,Fx_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','c')
    plot(S2ratios,-Fy_damaged_mean_steady,'ok','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S2ratios,-Fz_damaged_mean_steady,'ok','markersize',10,'markerfacecolor',[0 .5 0])
    
    subplot(2,2,2)
    hold on
    plot(S3ratios,Mx_damaged_mean_all,'ok','markersize',10,'markerfacecolor','b')
    plot(S3ratios,-My_damaged_mean_all--My_damaged_mean_all(S3ratios==1),'ok','markersize',10,'markerfacecolor','r')
    plot(S3ratios,-Mz_damaged_mean_all,'ok','markersize',10,'markerfacecolor','g')
    
    plot(S3ratios,Mx_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','c')
    plot(S3ratios,-My_damaged_mean_steady--My_damaged_mean_all(S3ratios==1),'ok','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S3ratios,-Mz_damaged_mean_steady,'ok','markersize',10,'markerfacecolor',[0 .5 0])
    
    % compare MODs -Fz & Mx
    figure(2)
    subplot(2,2,1)
    hold on
    plot(S2ratios,-Fz_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','k')
    plot(S2ratios,-Fz_damaged_mean_all,'ok','markersize',10,'markerfacecolor',[.5 .5 .5])
    plot(S2ratios,-Fz_damaged_mean_freqMOD,'ok','markersize',10,'markerfacecolor','b')
    plot(S2ratios,-Fz_damaged_mean_strokeMOD,'ok','markersize',10,'markerfacecolor','g')
    plot(S2ratios,-Fz_damaged_mean_devMOD,'ok','markersize',10,'markerfacecolor','r')
    plot(S2ratios,-Fz_damaged_mean_rotMOD,'ok','markersize',10,'markerfacecolor','c')
    
    subplot(2,2,2)
    hold on
    plot(S3ratios,Mx_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','k')
    plot(S3ratios,Mx_damaged_mean_all,'ok','markersize',10,'markerfacecolor',[.5 .5 .5])
    plot(S3ratios,Mx_damaged_mean_freqMOD,'ok','markersize',10,'markerfacecolor','b')
    plot(S3ratios,Mx_damaged_mean_strokeMOD,'ok','markersize',10,'markerfacecolor','g')
    plot(S3ratios,Mx_damaged_mean_devMOD,'ok','markersize',10,'markerfacecolor','r')
    plot(S3ratios,Mx_damaged_mean_rotMOD,'ok','markersize',10,'markerfacecolor','c')
    
    % compare MODs: -Fy & -Mz
    figure(3)
    subplot(2,2,1)
    hold on
    plot(S2ratios,-Fy_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','k')
    plot(S2ratios,-Fy_damaged_mean_all,'ok','markersize',10,'markerfacecolor',[.5 .5 .5])
    plot(S2ratios,-Fy_damaged_mean_freqMOD,'ok','markersize',10,'markerfacecolor','b')
    plot(S2ratios,-Fy_damaged_mean_strokeMOD,'ok','markersize',10,'markerfacecolor','g')
    plot(S2ratios,-Fy_damaged_mean_devMOD,'ok','markersize',10,'markerfacecolor','r')
    plot(S2ratios,-Fy_damaged_mean_rotMOD,'ok','markersize',10,'markerfacecolor','c')
    
    subplot(2,2,2)
    hold on
    plot(S3ratios,-Mz_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','k')
    plot(S3ratios,-Mz_damaged_mean_all,'ok','markersize',10,'markerfacecolor',[.5 .5 .5])
    plot(S3ratios,-Mz_damaged_mean_freqMOD,'ok','markersize',10,'markerfacecolor','b')
    plot(S3ratios,-Mz_damaged_mean_strokeMOD,'ok','markersize',10,'markerfacecolor','g')
    plot(S3ratios,-Mz_damaged_mean_devMOD,'ok','markersize',10,'markerfacecolor','r')
    plot(S3ratios,-Mz_damaged_mean_rotMOD,'ok','markersize',10,'markerfacecolor','c')
    
    % compare MODs -Fz & Mx COMBIMOD FREQ FIRST
    figure(4)
    subplot(2,2,1)
    hold on
    plot(S2ratios,-Fz_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','w')
    plot(S2ratios,-Fz_damaged_mean_freqMOD,'ok','markersize',10,'markerfacecolor','y')
    plot(S2ratios,-Fz_damaged_mean_freqStrokeMOD,'ok','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S2ratios,-Fz_damaged_mean_freqStrokeDevMOD,'ok','markersize',10,'markerfacecolor','r')
    plot(S2ratios,-Fz_damaged_mean_all,'ok','markersize',10,'markerfacecolor',[.5 0 0])
    
    subplot(2,2,2)
    hold on
    plot(S3ratios,Mx_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','w')
    plot(S3ratios,Mx_damaged_mean_freqMOD,'ok','markersize',10,'markerfacecolor','y')
    plot(S3ratios,Mx_damaged_mean_freqStrokeMOD,'ok','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S3ratios,Mx_damaged_mean_freqStrokeDevMOD,'ok','markersize',10,'markerfacecolor','r')
    plot(S3ratios,Mx_damaged_mean_all,'ok','markersize',10,'markerfacecolor',[.5 0 0])
    
    % compare MODs -Fz & Mx COMBIMOD STROKE FIRST
    figure(5)
    subplot(2,2,1)
    hold on
    plot(S2ratios,-Fz_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','w')
    plot(S2ratios,-Fz_damaged_mean_strokeMOD,'ok','markersize',10,'markerfacecolor','y')
    plot(S2ratios,-Fz_damaged_mean_strokeDevMOD,'ok','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S2ratios,-Fz_damaged_mean_strokeDevRotMOD,'ok','markersize',10,'markerfacecolor','r')
    plot(S2ratios,-Fz_damaged_mean_all,'ok','markersize',10,'markerfacecolor',[.5 0 0])
    
    subplot(2,2,2)
    hold on
    plot(S3ratios,Mx_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','w')
    plot(S3ratios,Mx_damaged_mean_strokeMOD,'ok','markersize',10,'markerfacecolor','y')
    plot(S3ratios,Mx_damaged_mean_strokeDevMOD,'ok','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S3ratios,Mx_damaged_mean_strokeDevRotMOD,'ok','markersize',10,'markerfacecolor','r')
    plot(S3ratios,Mx_damaged_mean_all,'ok','markersize',10,'markerfacecolor',[.5 0 0])
    
    %% plot TE cut wing data
    load(['allMODs_TEclip_freqAsym',num2str(freq_asymFitNr),'.mat'])

    S2ratios_TEnTip = [S2ratios_TEnTip;S2ratios];
    S3ratios_TEnTip = [S3ratios_TEnTip;S3ratios];
    S2S3AmpRatioFuncs_TEnTip = [S2S3AmpRatioFuncs_TEnTip;S2S3AmpRatioFuncs];

    Fx_damaged_mean_all_TEnTip = [Fx_damaged_mean_all_TEnTip Fx_damaged_mean_all];
    Fy_damaged_mean_all_TEnTip = [Fy_damaged_mean_all_TEnTip Fy_damaged_mean_all];
    Fz_damaged_mean_all_TEnTip = [Fz_damaged_mean_all_TEnTip Fz_damaged_mean_all];
    
    Mx_damaged_mean_all_TEnTip = [Mx_damaged_mean_all_TEnTip Mx_damaged_mean_all];
    My_damaged_mean_all_TEnTip = [My_damaged_mean_all_TEnTip My_damaged_mean_all];
    Mz_damaged_mean_all_TEnTip = [Mz_damaged_mean_all_TEnTip Mz_damaged_mean_all];
    
    Fx_damaged_mean_steady_TEnTip = [Fx_damaged_mean_steady_TEnTip Fx_damaged_mean_steady];
    Fy_damaged_mean_steady_TEnTip = [Fy_damaged_mean_steady_TEnTip Fy_damaged_mean_steady];
    Fz_damaged_mean_steady_TEnTip = [Fz_damaged_mean_steady_TEnTip Fz_damaged_mean_steady];
    
    Mx_damaged_mean_steady_TEnTip = [Mx_damaged_mean_steady_TEnTip Mx_damaged_mean_steady];
    My_damaged_mean_steady_TEnTip = [My_damaged_mean_steady_TEnTip My_damaged_mean_steady];
    Mz_damaged_mean_steady_TEnTip = [Mz_damaged_mean_steady_TEnTip Mz_damaged_mean_steady];
    
    Fx_damaged_mean_freqMOD_TEnTip = [Fx_damaged_mean_freqMOD_TEnTip Fx_damaged_mean_freqMOD];
    Fy_damaged_mean_freqMOD_TEnTip = [Fy_damaged_mean_freqMOD_TEnTip Fy_damaged_mean_freqMOD];
    Fz_damaged_mean_freqMOD_TEnTip = [Fz_damaged_mean_freqMOD_TEnTip Fz_damaged_mean_freqMOD];
    
    Mx_damaged_mean_freqMOD_TEnTip = [Mx_damaged_mean_freqMOD_TEnTip Mx_damaged_mean_freqMOD];
    My_damaged_mean_freqMOD_TEnTip = [My_damaged_mean_freqMOD_TEnTip My_damaged_mean_freqMOD];
    Mz_damaged_mean_freqMOD_TEnTip = [Mz_damaged_mean_freqMOD_TEnTip Mz_damaged_mean_freqMOD];
    
    Fx_damaged_mean_strokeMOD_TEnTip = [Fx_damaged_mean_strokeMOD_TEnTip Fx_damaged_mean_strokeMOD];
    Fy_damaged_mean_strokeMOD_TEnTip = [Fy_damaged_mean_strokeMOD_TEnTip Fy_damaged_mean_strokeMOD];
    Fz_damaged_mean_strokeMOD_TEnTip = [Fz_damaged_mean_strokeMOD_TEnTip Fz_damaged_mean_strokeMOD];
    
    Mx_damaged_mean_strokeMOD_TEnTip = [Mx_damaged_mean_strokeMOD_TEnTip Mx_damaged_mean_strokeMOD];
    My_damaged_mean_strokeMOD_TEnTip = [My_damaged_mean_strokeMOD_TEnTip My_damaged_mean_strokeMOD];
    Mz_damaged_mean_strokeMOD_TEnTip = [Mz_damaged_mean_strokeMOD_TEnTip Mz_damaged_mean_strokeMOD];
    
    Fx_damaged_mean_devMOD_TEnTip = [Fx_damaged_mean_devMOD_TEnTip Fx_damaged_mean_devMOD];
    Fy_damaged_mean_devMOD_TEnTip = [Fy_damaged_mean_devMOD_TEnTip Fy_damaged_mean_devMOD];
    Fz_damaged_mean_devMOD_TEnTip = [Fz_damaged_mean_devMOD_TEnTip Fz_damaged_mean_devMOD];
    
    Mx_damaged_mean_devMOD_TEnTip = [Mx_damaged_mean_devMOD_TEnTip Mx_damaged_mean_devMOD];
    My_damaged_mean_devMOD_TEnTip = [My_damaged_mean_devMOD_TEnTip My_damaged_mean_devMOD];
    Mz_damaged_mean_devMOD_TEnTip = [Mz_damaged_mean_devMOD_TEnTip Mz_damaged_mean_devMOD];
    
    Fx_damaged_mean_rotMOD_TEnTip = [Fx_damaged_mean_rotMOD_TEnTip Fx_damaged_mean_rotMOD];
    Fy_damaged_mean_rotMOD_TEnTip = [Fy_damaged_mean_rotMOD_TEnTip Fy_damaged_mean_rotMOD];
    Fz_damaged_mean_rotMOD_TEnTip = [Fz_damaged_mean_rotMOD_TEnTip Fz_damaged_mean_rotMOD];
    
    Mx_damaged_mean_rotMOD_TEnTip = [Mx_damaged_mean_rotMOD_TEnTip Mx_damaged_mean_rotMOD];
    My_damaged_mean_rotMOD_TEnTip = [My_damaged_mean_rotMOD_TEnTip My_damaged_mean_rotMOD];
    Mz_damaged_mean_rotMOD_TEnTip = [Mz_damaged_mean_rotMOD_TEnTip Mz_damaged_mean_rotMOD];
    
    Fx_damaged_mean_freqStrokeMOD_TEnTip = [Fx_damaged_mean_freqStrokeMOD_TEnTip Fx_damaged_mean_freqStrokeMOD];
    Fy_damaged_mean_freqStrokeMOD_TEnTip = [Fy_damaged_mean_freqStrokeMOD_TEnTip Fy_damaged_mean_freqStrokeMOD];
    Fz_damaged_mean_freqStrokeMOD_TEnTip = [Fz_damaged_mean_freqStrokeMOD_TEnTip Fz_damaged_mean_freqStrokeMOD];
    
    Mx_damaged_mean_freqStrokeMOD_TEnTip = [Mx_damaged_mean_freqStrokeMOD_TEnTip Mx_damaged_mean_freqStrokeMOD];
    My_damaged_mean_freqStrokeMOD_TEnTip = [My_damaged_mean_freqStrokeMOD_TEnTip My_damaged_mean_freqStrokeMOD];
    Mz_damaged_mean_freqStrokeMOD_TEnTip = [Mz_damaged_mean_freqStrokeMOD_TEnTip Mz_damaged_mean_freqStrokeMOD];
    
    Fx_damaged_mean_freqStrokeDevMOD_TEnTip = [Fx_damaged_mean_freqStrokeDevMOD_TEnTip Fx_damaged_mean_freqStrokeDevMOD];
    Fy_damaged_mean_freqStrokeDevMOD_TEnTip = [Fy_damaged_mean_freqStrokeDevMOD_TEnTip Fy_damaged_mean_freqStrokeDevMOD];
    Fz_damaged_mean_freqStrokeDevMOD_TEnTip = [Fz_damaged_mean_freqStrokeDevMOD_TEnTip Fz_damaged_mean_freqStrokeDevMOD];
    
    Mx_damaged_mean_freqStrokeDevMOD_TEnTip = [Mx_damaged_mean_freqStrokeDevMOD_TEnTip Mx_damaged_mean_freqStrokeDevMOD];
    My_damaged_mean_freqStrokeDevMOD_TEnTip = [My_damaged_mean_freqStrokeDevMOD_TEnTip My_damaged_mean_freqStrokeDevMOD];
    Mz_damaged_mean_freqStrokeDevMOD_TEnTip = [Mz_damaged_mean_freqStrokeDevMOD_TEnTip Mz_damaged_mean_freqStrokeDevMOD];
    
    Fx_damaged_mean_strokeDevMOD_TEnTip = [Fx_damaged_mean_strokeDevMOD_TEnTip Fx_damaged_mean_strokeDevMOD];
    Fy_damaged_mean_strokeDevMOD_TEnTip = [Fy_damaged_mean_strokeDevMOD_TEnTip Fy_damaged_mean_strokeDevMOD];
    Fz_damaged_mean_strokeDevMOD_TEnTip = [Fz_damaged_mean_strokeDevMOD_TEnTip Fz_damaged_mean_strokeDevMOD];
    
    Mx_damaged_mean_strokeDevMOD_TEnTip = [Mx_damaged_mean_strokeDevMOD_TEnTip Mx_damaged_mean_strokeDevMOD];
    My_damaged_mean_strokeDevMOD_TEnTip = [My_damaged_mean_strokeDevMOD_TEnTip My_damaged_mean_strokeDevMOD];
    Mz_damaged_mean_strokeDevMOD_TEnTip = [Mz_damaged_mean_strokeDevMOD_TEnTip Mz_damaged_mean_strokeDevMOD];
    
    Fx_damaged_mean_strokeDevRotMOD_TEnTip = [Fx_damaged_mean_strokeDevRotMOD_TEnTip Fx_damaged_mean_strokeDevRotMOD];
    Fy_damaged_mean_strokeDevRotMOD_TEnTip = [Fy_damaged_mean_strokeDevRotMOD_TEnTip Fy_damaged_mean_strokeDevRotMOD];
    Fz_damaged_mean_strokeDevRotMOD_TEnTip = [Fz_damaged_mean_strokeDevRotMOD_TEnTip Fz_damaged_mean_strokeDevRotMOD];
    
    Mx_damaged_mean_strokeDevRotMOD_TEnTip = [Mx_damaged_mean_strokeDevRotMOD_TEnTip Mx_damaged_mean_strokeDevRotMOD];
    My_damaged_mean_strokeDevRotMOD_TEnTip = [My_damaged_mean_strokeDevRotMOD_TEnTip My_damaged_mean_strokeDevRotMOD];
    Mz_damaged_mean_strokeDevRotMOD_TEnTip = [Mz_damaged_mean_strokeDevRotMOD_TEnTip Mz_damaged_mean_strokeDevRotMOD];

    % ALL MODs VS NO MODs
    figure(1)
    subplot(2,2,1)
    hold on
    plot(S2ratios,Fx_damaged_mean_all,'dk','markersize',10,'markerfacecolor','b')
    plot(S2ratios,-Fy_damaged_mean_all,'dk','markersize',10,'markerfacecolor','r')
    plot(S2ratios,-Fz_damaged_mean_all,'dk','markersize',10,'markerfacecolor','g')
    
    plot(S2ratios,Fx_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','c')
    plot(S2ratios,-Fy_damaged_mean_steady,'dk','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S2ratios,-Fz_damaged_mean_steady,'dk','markersize',10,'markerfacecolor',[0 .5 0])
    
    legend('x-Mod','y-Mod','z-Mod','x-steady','y-steady','z-steady','location','E')
    xlabel('S2 ratio')
    ylabel('normalized forces F/mg')
    axis([0.5 1 0 1])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',0:1:1)
    
    subplot(2,2,2)
    hold on
    plot(S3ratios,Mx_damaged_mean_all,'dk','markersize',10,'markerfacecolor','b')
    plot(S3ratios,-My_damaged_mean_all--My_damaged_mean_all(S3ratios==1),'dk','markersize',10,'markerfacecolor','r')
    plot(S3ratios,-Mz_damaged_mean_all,'dk','markersize',10,'markerfacecolor','g')
    
    plot(S3ratios,Mx_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','c')
    plot(S3ratios,-My_damaged_mean_steady--My_damaged_mean_all(S3ratios==1),'dk','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S3ratios,-Mz_damaged_mean_steady,'dk','markersize',10,'markerfacecolor',[0 .5 0])
    
    xlabel('S3 ratio')
    ylabel('normalized torques T/mgl')
    axis([0.5 1 -.05 .15])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-.05:.05:.15)    
    
    % compare MODs -Fz & Mx
    
    figure(2)
    subplot(2,2,1)
    hold on
    plot(S2ratios,-Fz_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','k')
    plot(S2ratios,-Fz_damaged_mean_all,'dk','markersize',10,'markerfacecolor',[.5 .5 .5])
    plot(S2ratios,-Fz_damaged_mean_freqMOD,'dk','markersize',10,'markerfacecolor','b')
    plot(S2ratios,-Fz_damaged_mean_strokeMOD,'dk','markersize',10,'markerfacecolor','g')
    plot(S2ratios,-Fz_damaged_mean_devMOD,'dk','markersize',10,'markerfacecolor','r')
    plot(S2ratios,-Fz_damaged_mean_rotMOD,'dk','markersize',10,'markerfacecolor','c')
    
    legend('steady','modified','frequency','stroke','deviation','rotation','location','SE')
    xlabel('S2 ratio')
    ylabel('normalized vertical force -Fz/mg')
    axis([0.5 1 .5 1.5])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-1.5:.5:1.5) 
    
    subplot(2,2,2)
    hold on
    plot(S3ratios,Mx_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','k')
    plot(S3ratios,Mx_damaged_mean_all,'dk','markersize',10,'markerfacecolor',[.5 .5 .5])
    plot(S3ratios,Mx_damaged_mean_freqMOD,'dk','markersize',10,'markerfacecolor','b')
    plot(S3ratios,Mx_damaged_mean_strokeMOD,'dk','markersize',10,'markerfacecolor','g')
    plot(S3ratios,Mx_damaged_mean_devMOD,'dk','markersize',10,'markerfacecolor','r')
    plot(S3ratios,Mx_damaged_mean_rotMOD,'dk','markersize',10,'markerfacecolor','c')
    
    xlabel('S3 ratio')
    ylabel('normalized roll torque Tx/mgl')
    axis([0.5 1 0 .2])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-.2:.2:.2)    

    % compare MODs: -Fy & -Mz
    
    figure(3)
    subplot(2,2,1)
    hold on
    plot(S2ratios,-Fy_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','k')
    plot(S2ratios,-Fy_damaged_mean_all,'dk','markersize',10,'markerfacecolor',[.5 .5 .5])
    plot(S2ratios,-Fy_damaged_mean_freqMOD,'dk','markersize',10,'markerfacecolor','b')
    plot(S2ratios,-Fy_damaged_mean_strokeMOD,'dk','markersize',10,'markerfacecolor','g')
    plot(S2ratios,-Fy_damaged_mean_devMOD,'dk','markersize',10,'markerfacecolor','r')
    plot(S2ratios,-Fy_damaged_mean_rotMOD,'dk','markersize',10,'markerfacecolor','c')
    
    legend('steady','modified','frequency','stroke','deviation','rotation','location','NE')
    xlabel('S2 ratio')
    ylabel('normalized sideways force -Fy/mg')
    axis([0.5 1 0 .5])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-1.5:.5:1.5) 
    
    subplot(2,2,2)
    hold on
    plot(S3ratios,-Mz_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','k')
    plot(S3ratios,-Mz_damaged_mean_all,'dk','markersize',10,'markerfacecolor',[.5 .5 .5])
    plot(S3ratios,-Mz_damaged_mean_freqMOD,'dk','markersize',10,'markerfacecolor','b')
    plot(S3ratios,-Mz_damaged_mean_strokeMOD,'dk','markersize',10,'markerfacecolor','g')
    plot(S3ratios,-Mz_damaged_mean_devMOD,'dk','markersize',10,'markerfacecolor','r')
    plot(S3ratios,-Mz_damaged_mean_rotMOD,'dk','markersize',10,'markerfacecolor','c')
    
    xlabel('S3 ratio')
    ylabel('normalized yaw torque Tz/mgl')
    axis([0.5 1 -.05 .15])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-.2:.05:.3)    
    
    % compare MODs -Fz & Mx COMBIMOD FREQ FIRST
    figure(4)
    subplot(2,2,1)
    hold on
    plot(S2ratios,-Fz_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','w')
    plot(S2ratios,-Fz_damaged_mean_freqMOD,'dk','markersize',10,'markerfacecolor','y')
    plot(S2ratios,-Fz_damaged_mean_freqStrokeMOD,'dk','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S2ratios,-Fz_damaged_mean_freqStrokeDevMOD,'dk','markersize',10,'markerfacecolor','r')
    plot(S2ratios,-Fz_damaged_mean_all,'dk','markersize',10,'markerfacecolor',[.5 0 0])
    
    xlabel('S2 ratio')
    ylabel('normalized vertical force -Fz/mg')
    axis([0.5 1 .5 1.5])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-1.5:.5:1.5) 

    subplot(2,2,2)
    hold on
    plot(S3ratios,Mx_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','w')
    plot(S3ratios,Mx_damaged_mean_freqMOD,'dk','markersize',10,'markerfacecolor','y')
    plot(S3ratios,Mx_damaged_mean_freqStrokeMOD,'dk','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S3ratios,Mx_damaged_mean_freqStrokeDevMOD,'dk','markersize',10,'markerfacecolor','r')
    plot(S3ratios,Mx_damaged_mean_all,'dk','markersize',10,'markerfacecolor',[.5 0 0])
    
    xlabel('S3 ratio')
    ylabel('normalized roll torque Tx/mgl')
    axis([0.5 1 0 .2])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-.2:.2:.2)    
    
    subplot(2,2,3)
    hold on
    plot(0,0,'ok','markersize',10,'markerfacecolor','w')
    plot(0,0,'ok','markersize',10,'markerfacecolor','y')
    plot(0,0,'ok','markersize',10,'markerfacecolor',[1 .5 0])
    plot(0,0,'ok','markersize',10,'markerfacecolor','r')
    plot(0,0,'ok','markersize',10,'markerfacecolor',[.5 0 0])
    legend('none','frequency','freq&stroke','freq&str&dev','all','location','NE')
    colormap(cmap_MODs)
    colorbar
    
    % compare MODs -Fz & Mx COMBIMOD STROKE FIRST
    figure(5)
    subplot(2,2,1)
    hold on
    plot(S2ratios,-Fz_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','w')
    plot(S2ratios,-Fz_damaged_mean_strokeMOD,'dk','markersize',10,'markerfacecolor','y')
    plot(S2ratios,-Fz_damaged_mean_strokeDevMOD,'dk','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S2ratios,-Fz_damaged_mean_strokeDevRotMOD,'dk','markersize',10,'markerfacecolor','r')
    plot(S2ratios,-Fz_damaged_mean_all,'dk','markersize',10,'markerfacecolor',[.5 0 0])
    
    xlabel('S2 ratio')
    ylabel('normalized vertical force -Fz/mg')
    axis([0.5 1 .5 1.5])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-1.5:.5:1.5) 

    subplot(2,2,2)
    hold on
    plot(S3ratios,Mx_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','w')
    plot(S3ratios,Mx_damaged_mean_strokeMOD,'dk','markersize',10,'markerfacecolor','y')
    plot(S3ratios,Mx_damaged_mean_strokeDevMOD,'dk','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S3ratios,Mx_damaged_mean_strokeDevRotMOD,'dk','markersize',10,'markerfacecolor','r')
    plot(S3ratios,Mx_damaged_mean_all,'dk','markersize',10,'markerfacecolor',[.5 0 0])
    
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
    legend('none','stroke','str&dev','str&dev&rot','all','location','NE')
    colormap(cmap_MODs)
    colorbar
    
    %% plot trendlines
    figure(1)
    subplot(2,2,1)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,Fx_damaged_mean_all_TEnTip',1),[.5 1]),'b','linewidth',2)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fy_damaged_mean_all_TEnTip',1),[.5 1]),'r','linewidth',2)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',1),[.5 1]),'g','linewidth',2)

    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,Fx_damaged_mean_steady_TEnTip',1),[.5 1]),'c','linewidth',2)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fy_damaged_mean_steady_TEnTip',1),[.5 1]),'color',[1 .5 0],'linewidth',2)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_steady_TEnTip',1),[.5 1]),'color',[0 .5 0],'linewidth',2)

    % adjusted
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',2),[.5 1]),'g','linewidth',2)

    subplot(2,2,2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_all_TEnTip',1),[.5 1]),'b','linewidth',2)
%     plot([.5 1],polyval(polyfit(S3ratios_TEnTip,-My_damaged_mean_all_TEnTip',1),[.5 1]),'r','linewidth',2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,-Mz_damaged_mean_all_TEnTip',1),[.5 1]),'g','linewidth',2)

    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_steady_TEnTip',1),[.5 1]),'c','linewidth',2)
%     plot([.5 1],polyval(polyfit(S3ratios_TEnTip,-My_damaged_mean_steady_TEnTip',1),[.5 1]),'color',[1 .5 0],'linewidth',2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,-Mz_damaged_mean_steady_TEnTip',1),[.5 1]),'color',[0 .5 0],'linewidth',2)

    % adjusted
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,-My_damaged_mean_all_TEnTip'--My_damaged_mean_all(S3ratios==1),1),[.5 1]),'r','linewidth',2)
%     plot([.5 1],polyval(polyfit(S3ratios_TEnTip,-Mz_damaged_mean_all_TEnTip',2),[.5 1]),'g','linewidth',2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,-My_damaged_mean_steady_TEnTip'--My_damaged_mean_steady(S3ratios==1),1),[.5 1]),'color',[1 .5 0],'linewidth',2)

    figure(2)
    subplot(2,2,1)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_steady_TEnTip',1),[.5 1]),'k','linewidth',2)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',1),[.5 1]),'color',[.5 .5 .5],'linewidth',2)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqMOD_TEnTip',1),[.5 1]),'b','linewidth',2)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_strokeMOD_TEnTip',1),[.5 1]),'g','linewidth',2)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_devMOD_TEnTip',1),[.5 1]),'r','linewidth',2)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_rotMOD_TEnTip',1),[.5 1]),'c','linewidth',2)

    % adjusted
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',2),[.5 1]),'color',[.5 .5 .5],'linewidth',2)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqMOD_TEnTip',2),[.5 1]),'b','linewidth',2)

    subplot(2,2,2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_steady_TEnTip',1),[.5 1]),'k','linewidth',2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_all_TEnTip',1),[.5 1]),'color',[.5 .5 .5],'linewidth',2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqMOD_TEnTip',1),[.5 1]),'b','linewidth',2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_strokeMOD_TEnTip',1),[.5 1]),'g','linewidth',2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_devMOD_TEnTip',1),[.5 1]),'r','linewidth',2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_rotMOD_TEnTip',1),[.5 1]),'c','linewidth',2)

    figure(3)
    subplot(2,2,1)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fy_damaged_mean_steady_TEnTip',1),[.5 1]),'k','linewidth',2)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fy_damaged_mean_all_TEnTip',1),[.5 1]),'color',[.5 .5 .5],'linewidth',2)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fy_damaged_mean_freqMOD_TEnTip',1),[.5 1]),'b','linewidth',2)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fy_damaged_mean_strokeMOD_TEnTip',1),[.5 1]),'g','linewidth',2)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fy_damaged_mean_devMOD_TEnTip',1),[.5 1]),'r','linewidth',2)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fy_damaged_mean_rotMOD_TEnTip',1),[.5 1]),'c','linewidth',2)

    subplot(2,2,2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,-Mz_damaged_mean_steady_TEnTip',1),[.5 1]),'k','linewidth',2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,-Mz_damaged_mean_all_TEnTip',1),[.5 1]),'color',[.5 .5 .5],'linewidth',2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,-Mz_damaged_mean_freqMOD_TEnTip',1),[.5 1]),'b','linewidth',2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,-Mz_damaged_mean_strokeMOD_TEnTip',1),[.5 1]),'g','linewidth',2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,-Mz_damaged_mean_devMOD_TEnTip',1),[.5 1]),'r','linewidth',2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,-Mz_damaged_mean_rotMOD_TEnTip',1),[.5 1]),'c','linewidth',2)

    figure(4)
    subplot(2,2,1)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_steady_TEnTip',1),[.5 1]),'k','linewidth',2)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqMOD_TEnTip',1),[.5 1]),'y','linewidth',2)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqStrokeMOD_TEnTip',1),[.5 1]),'color',[1 .5 0],'linewidth',2)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqStrokeDevMOD_TEnTip',1),[.5 1]),'r','linewidth',2)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',1),[.5 1]),'color',[.5 0 0],'linewidth',2)

%     % adjusted
%     plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',2),[.5 1]),'color',[.5 .5 .5],'linewidth',2)
%     plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqMOD_TEnTip',2),[.5 1]),'b','linewidth',2)

    subplot(2,2,2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_steady_TEnTip',1),[.5 1]),'k','linewidth',2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqMOD_TEnTip',1),[.5 1]),'y','linewidth',2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqStrokeMOD_TEnTip',1),[.5 1]),'color',[1 .5 0],'linewidth',2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqStrokeDevMOD_TEnTip',1),[.5 1]),'r','linewidth',2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_all_TEnTip',1),[.5 1]),'color',[.5 0 0],'linewidth',2)

    figure(5)
    subplot(2,2,1)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_steady_TEnTip',1),[.5 1]),'k','linewidth',2)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_strokeMOD_TEnTip',1),[.5 1]),'y','linewidth',2)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_strokeDevMOD_TEnTip',1),[.5 1]),'color',[1 .5 0],'linewidth',2)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_strokeDevRotMOD_TEnTip',1),[.5 1]),'r','linewidth',2)
    plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',1),[.5 1]),'color',[.5 0 0],'linewidth',2)

%     % adjusted
%     plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',2),[.5 1]),'color',[.5 .5 .5],'linewidth',2)
%     plot([.5 1],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqMOD_TEnTip',2),[.5 1]),'b','linewidth',2)

    subplot(2,2,2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_steady_TEnTip',1),[.5 1]),'k','linewidth',2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_strokeMOD_TEnTip',1),[.5 1]),'y','linewidth',2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_strokeDevMOD_TEnTip',1),[.5 1]),'color',[1 .5 0],'linewidth',2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_strokeDevRotMOD_TEnTip',1),[.5 1]),'r','linewidth',2)
    plot([.5 1],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_all_TEnTip',1),[.5 1]),'color',[.5 0 0],'linewidth',2)

    %% percentage of weight support and roll equilibrium FREQ FIRST
    Fz_steady = -Fz_damaged_mean_steady_TEnTip(1);
    Fz_none = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_steady_TEnTip',1),[.5 1]));
    Fz_freq = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqMOD_TEnTip',1),[.5 1]));
    Fz_freqStroke = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqStrokeMOD_TEnTip',1),[.5 1]));
    Fz_freqStrokeDev = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqStrokeDevMOD_TEnTip',1),[.5 1]));
    Fz_freqStrokeDevRot = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',1),[.5 1]));

    Fz_freq_perc    = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_freq-Fz_none)
    Fz_stroke_perc  = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_freqStroke-Fz_freq)
    Fz_dev_perc     = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_freqStrokeDev-Fz_freqStroke)
    Fz_rot_perc     = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_freqStrokeDevRot-Fz_freqStrokeDev)
    Fz_all_perc     = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_freqStrokeDevRot-Fz_none)
    Fz_sum_perc     = Fz_freq_perc + Fz_stroke_perc + Fz_dev_perc + Fz_rot_perc
    
    Mx_steady = Mx_damaged_mean_steady_TEnTip(1);
    Mx_none = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_steady_TEnTip',1),[.5 1]));
    Mx_freq = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqMOD_TEnTip',1),[.5 1]));
    Mx_freqStroke = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqStrokeMOD_TEnTip',1),[.5 1]));
    Mx_freqStrokeDev = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqStrokeDevMOD_TEnTip',1),[.5 1]));
    Mx_freqStrokeDevRot = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_all_TEnTip',1),[.5 1]));

    Mx_freq_perc    = 100 / (-Mx_none) * (Mx_freq-Mx_none)
    Mx_stroke_perc  = 100 / (-Mx_none) * (Mx_freqStroke-Mx_freq)
    Mx_dev_perc     = 100 / (-Mx_none) * (Mx_freqStrokeDev-Mx_freqStroke)
    Mx_rot_perc     = 100 / (-Mx_none) * (Mx_freqStrokeDevRot-Mx_freqStrokeDev)
    Mx_all_perc     = 100 / (-Mx_none) * (Mx_freqStrokeDevRot-Mx_none)
    Mx_sum_perc     = Mx_freq_perc + Mx_stroke_perc + Mx_dev_perc + Mx_rot_perc
    
    %% percentage of weight support and roll equilibrium STROKE FIRST
    Fz_steady = -Fz_damaged_mean_steady_TEnTip(1);
    Fz_none = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_steady_TEnTip',1),[.5 1]));
    Fz_stroke = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_strokeMOD_TEnTip',1),[.5 1]));
    Fz_strokeDev = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_strokeDevMOD_TEnTip',1),[.5 1]));
    Fz_strokeDevRot = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_strokeDevRotMOD_TEnTip',1),[.5 1]));
    Fz_strokeDevRotFreq = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',1),[.5 1]));

    Fz_stroke_perc    = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_stroke-Fz_none)
    Fz_dev_perc     = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_strokeDev-Fz_stroke)
    Fz_rot_perc     = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_strokeDevRot-Fz_strokeDev)
    Fz_freq_perc     = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_strokeDevRotFreq-Fz_strokeDevRot)
    Fz_all_perc     = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_strokeDevRotFreq-Fz_none)
    Fz_sum_perc     = Fz_freq_perc + Fz_stroke_perc + Fz_dev_perc + Fz_rot_perc
    
    Mx_steady = Mx_damaged_mean_steady_TEnTip(1);
    Mx_none = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_steady_TEnTip',1),[.5 1]));
    Mx_stroke = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_strokeMOD_TEnTip',1),[.5 1]));
    Mx_strokeDev = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_strokeDevMOD_TEnTip',1),[.5 1]));
    Mx_strokeDevRot = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_strokeDevRotMOD_TEnTip',1),[.5 1]));
    Mx_strokeDevRotFreq = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_all_TEnTip',1),[.5 1]));

    Mx_stroke_perc    = 100 / (-Mx_none) * (Mx_stroke-Mx_none)
    Mx_dev_perc     = 100 / (-Mx_none) * (Mx_strokeDev-Mx_stroke)
    Mx_rot_perc     = 100 / (-Mx_none) * (Mx_strokeDevRot-Mx_strokeDev)
    Mx_freq_perc     = 100 / (-Mx_none) * (Mx_strokeDevRotFreq-Mx_strokeDevRot)
    Mx_all_perc     = 100 / (-Mx_none) * (Mx_strokeDevRotFreq-Mx_none)
    Mx_sum_perc     = Mx_freq_perc + Mx_stroke_perc + Mx_dev_perc + Mx_rot_perc
    
    %% save plots
    
    mkdir('qsModel_FnM_TEnTipCut')
    cd('qsModel_FnM_TEnTipCut')
    
    figure(1)
    saveas(gca,['FnM_WBmod_TEnTipClip_YnZflip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['FnM_WBmod_TEnTipClip_YnZflip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['FnM_WBmod_TEnTipClip_YnZflip_asympFit',num2str(freq_asymFitNr),'.svg'])

    figure(2)
    saveas(gca,['Fz_Mx_WBmodComponents_TEnTipClip_YnZflip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['Fz_Mx_WBmodComponents_TEnTipClip_YnZflip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['Fz_Mx_WBmodComponents_TEnTipClip_YnZflip_asympFit',num2str(freq_asymFitNr),'.svg'])

    figure(3)
    saveas(gca,['Fy_Mz_WBmodComponents_TEnTipclip_YnZflip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['Fy_Mz_WBmodComponents_TEnTipclip_YnZflip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['Fy_Mz_WBmodComponents_TEnTipclip_YnZflip_asympFit',num2str(freq_asymFitNr),'.svg'])

    figure(4)
    saveas(gca,['Fz_Mx_WBmodCombies_freqFirst_TEnTipclip_YnZflip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['Fz_Mx_WBmodCombies_freqFirst_TEnTipclip_YnZflip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['Fz_Mx_WBmodCombies_freqFirst_TEnTipclip_YnZflip_asympFit',num2str(freq_asymFitNr),'.svg'])

    figure(5)
    subplot(2,2,3)
    colormap(cmap_MODs)
    colorbar
    saveas(gca,['Fz_Mx_WBmodCombies_strokeFirst_TEnTipclip_YnZflip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['Fz_Mx_WBmodCombies_strokeFirst_TEnTipclip_YnZflip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['Fz_Mx_WBmodCombies_strokeFirst_TEnTipclip_YnZflip_asympFit',num2str(freq_asymFitNr),'.svg'])

    cd ..
    
    
    