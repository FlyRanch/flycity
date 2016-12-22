clear;
clc;
close all
warning off

cmap_MODs_red = [1 1 1; 1 1 0; 1 .5 0; 1 0 0; .5 0 0]; 
cmap_MODs_blue = [1 1 1; 0 1 1; 0 .5 1; 0 0 1; 0 0 .5]; 

freq_asymFitNr = 10;
peakloc = 0;

Aratio_max = 1.3;
Aratio_min = 1;


        %% plot TIP cut wing data
    load(['allMODs_TipClip_freqAsym',num2str(freq_asymFitNr),'_roton1.mat'])

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
    
    Fx_damaged_mean_freqRotMOD_TEnTip = Fx_damaged_mean_freqRotMOD;
    Fy_damaged_mean_freqRotMOD_TEnTip = Fy_damaged_mean_freqRotMOD;
    Fz_damaged_mean_freqRotMOD_TEnTip = Fz_damaged_mean_freqRotMOD;
    
    Mx_damaged_mean_freqRotMOD_TEnTip = Mx_damaged_mean_freqRotMOD;
    My_damaged_mean_freqRotMOD_TEnTip = My_damaged_mean_freqRotMOD;
    Mz_damaged_mean_freqRotMOD_TEnTip = Mz_damaged_mean_freqRotMOD;
    
    Fx_damaged_mean_freqRotDevMOD_TEnTip = Fx_damaged_mean_freqRotDevMOD;
    Fy_damaged_mean_freqRotDevMOD_TEnTip = Fy_damaged_mean_freqRotDevMOD;
    Fz_damaged_mean_freqRotDevMOD_TEnTip = Fz_damaged_mean_freqRotDevMOD;
    
    Mx_damaged_mean_freqRotDevMOD_TEnTip = Mx_damaged_mean_freqRotDevMOD;
    My_damaged_mean_freqRotDevMOD_TEnTip = My_damaged_mean_freqRotDevMOD;
    Mz_damaged_mean_freqRotDevMOD_TEnTip = Mz_damaged_mean_freqRotDevMOD;
    
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
    
    % ALL S2 & S3
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
    
    % compare MODs -Fz & Mx COMBIMOD FREQ FIRST & STROKE LAST
    figure(4)
    subplot(2,2,1)
    hold on
    plot(S2ratios,-Fz_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','w')
    plot(S2ratios,-Fz_damaged_mean_freqMOD,'ok','markersize',10,'markerfacecolor','y')
    plot(S2ratios,-Fz_damaged_mean_freqRotMOD,'ok','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S2ratios,-Fz_damaged_mean_freqRotDevMOD,'ok','markersize',10,'markerfacecolor','r')
    plot(S2ratios,-Fz_damaged_mean_all,'ok','markersize',10,'markerfacecolor',[.5 0 0])
    
    subplot(2,2,2)
    hold on
    plot(S3ratios,Mx_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','w')
    plot(S3ratios,Mx_damaged_mean_freqMOD,'ok','markersize',10,'markerfacecolor','y')
    plot(S3ratios,Mx_damaged_mean_freqRotMOD,'ok','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S3ratios,Mx_damaged_mean_freqRotDevMOD,'ok','markersize',10,'markerfacecolor','r')
    plot(S3ratios,Mx_damaged_mean_all,'ok','markersize',10,'markerfacecolor',[.5 0 0])
    
    % compare MODs -Fz & Mx COMBIMOD FREQ FIRST & ROT LAST
    figure(5)
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
    figure(6)
    subplot(2,2,1)
    hold on
    plot(S2ratios,-Fz_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','w')
    plot(S2ratios,-Fz_damaged_mean_strokeMOD,'ok','markersize',10,'markerfacecolor','c')
    plot(S2ratios,-Fz_damaged_mean_strokeDevMOD,'ok','markersize',10,'markerfacecolor',[0 .5 1])
    plot(S2ratios,-Fz_damaged_mean_strokeDevRotMOD,'ok','markersize',10,'markerfacecolor','b')
    plot(S2ratios,-Fz_damaged_mean_all,'ok','markersize',10,'markerfacecolor',[0 0 .5])
    
    subplot(2,2,2)
    hold on
    plot(S3ratios,Mx_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','w')
    plot(S3ratios,Mx_damaged_mean_strokeMOD,'ok','markersize',10,'markerfacecolor','c')
    plot(S3ratios,Mx_damaged_mean_strokeDevMOD,'ok','markersize',10,'markerfacecolor',[0 .5 1])
    plot(S3ratios,Mx_damaged_mean_strokeDevRotMOD,'ok','markersize',10,'markerfacecolor','b')
    plot(S3ratios,Mx_damaged_mean_all,'ok','markersize',10,'markerfacecolor',[0 0 .5])
    
    %% plot TE cut wing data
    load(['allMODs_TEclip_freqAsym',num2str(freq_asymFitNr),'_roton1.mat'])

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
    
    Fx_damaged_mean_freqRotMOD_TEnTip = [Fx_damaged_mean_freqRotMOD_TEnTip Fx_damaged_mean_freqRotMOD];
    Fy_damaged_mean_freqRotMOD_TEnTip = [Fy_damaged_mean_freqRotMOD_TEnTip Fy_damaged_mean_freqRotMOD];
    Fz_damaged_mean_freqRotMOD_TEnTip = [Fz_damaged_mean_freqRotMOD_TEnTip Fz_damaged_mean_freqRotMOD];
    
    Mx_damaged_mean_freqRotMOD_TEnTip = [Mx_damaged_mean_freqRotMOD_TEnTip Mx_damaged_mean_freqRotMOD];
    My_damaged_mean_freqRotMOD_TEnTip = [My_damaged_mean_freqRotMOD_TEnTip My_damaged_mean_freqRotMOD];
    Mz_damaged_mean_freqRotMOD_TEnTip = [Mz_damaged_mean_freqRotMOD_TEnTip Mz_damaged_mean_freqRotMOD];
    
    Fx_damaged_mean_freqRotDevMOD_TEnTip = [Fx_damaged_mean_freqRotDevMOD_TEnTip Fx_damaged_mean_freqRotDevMOD];
    Fy_damaged_mean_freqRotDevMOD_TEnTip = [Fy_damaged_mean_freqRotDevMOD_TEnTip Fy_damaged_mean_freqRotDevMOD];
    Fz_damaged_mean_freqRotDevMOD_TEnTip = [Fz_damaged_mean_freqRotDevMOD_TEnTip Fz_damaged_mean_freqRotDevMOD];
    
    Mx_damaged_mean_freqRotDevMOD_TEnTip = [Mx_damaged_mean_freqRotDevMOD_TEnTip Mx_damaged_mean_freqRotDevMOD];
    My_damaged_mean_freqRotDevMOD_TEnTip = [My_damaged_mean_freqRotDevMOD_TEnTip My_damaged_mean_freqRotDevMOD];
    Mz_damaged_mean_freqRotDevMOD_TEnTip = [Mz_damaged_mean_freqRotDevMOD_TEnTip Mz_damaged_mean_freqRotDevMOD];
    
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

    % ALL S2 & S3
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
    
    % compare MODs -Fz & Mx COMBIMOD FREQ FIRST & STROKE LAST
    figure(4)
    subplot(2,2,1)
    hold on
    plot(S2ratios,-Fz_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','w')
    plot(S2ratios,-Fz_damaged_mean_freqMOD,'dk','markersize',10,'markerfacecolor','y')
    plot(S2ratios,-Fz_damaged_mean_freqRotMOD,'dk','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S2ratios,-Fz_damaged_mean_freqRotDevMOD,'dk','markersize',10,'markerfacecolor','r')
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
    plot(S3ratios,Mx_damaged_mean_freqRotMOD,'dk','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S3ratios,Mx_damaged_mean_freqRotDevMOD,'dk','markersize',10,'markerfacecolor','r')
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
    legend('none','frequency','freq&rot','freq&rot&dev','all','location','NE')
    colormap(cmap_MODs_red)
    colorbar
    
    % compare MODs -Fz & Mx COMBIMOD FREQ FIRST & ROT LAST
    figure(5)
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
    colormap(cmap_MODs_red)
    colorbar
    
    % compare MODs -Fz & Mx COMBIMOD STROKE FIRST
    figure(6)
    subplot(2,2,1)
    hold on
    plot(S2ratios,-Fz_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','w')
    plot(S2ratios,-Fz_damaged_mean_strokeMOD,'dk','markersize',10,'markerfacecolor','c')
    plot(S2ratios,-Fz_damaged_mean_strokeDevMOD,'dk','markersize',10,'markerfacecolor',[0 .5 1])
    plot(S2ratios,-Fz_damaged_mean_strokeDevRotMOD,'dk','markersize',10,'markerfacecolor','b')
    plot(S2ratios,-Fz_damaged_mean_all,'dk','markersize',10,'markerfacecolor',[0 0 .5])
    
    xlabel('S2 ratio')
    ylabel('normalized vertical force -Fz/mg')
    axis([0.5 1 .5 1.5])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-1.5:.5:1.5) 

    subplot(2,2,2)
    hold on
    plot(S3ratios,Mx_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','w')
    plot(S3ratios,Mx_damaged_mean_strokeMOD,'dk','markersize',10,'markerfacecolor','c')
    plot(S3ratios,Mx_damaged_mean_strokeDevMOD,'dk','markersize',10,'markerfacecolor',[0 .5 1])
    plot(S3ratios,Mx_damaged_mean_strokeDevRotMOD,'dk','markersize',10,'markerfacecolor','b')
    plot(S3ratios,Mx_damaged_mean_all,'dk','markersize',10,'markerfacecolor',[0 0 .5])
    
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
    legend('none','stroke','str&dev','str&dev&rot','all','location','NE')
    colormap(cmap_MODs_blue)
    colorbar
    
    %% plot trendlines
    figure(1)
    subplot(2,2,1)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,Fx_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'b','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fy_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'r','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'g','linewidth',2)

    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,Fx_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'c','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fy_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[1 .5 0],'linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[0 .5 0],'linewidth',2)

    % adjusted
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',2),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'g','linewidth',2)

    subplot(2,2,2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_all_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'b','linewidth',2)
%     plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,-My_damaged_mean_all_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'r','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,-Mz_damaged_mean_all_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'g','linewidth',2)

    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_steady_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'c','linewidth',2)
%     plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,-My_damaged_mean_steady_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color',[1 .5 0],'linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,-Mz_damaged_mean_steady_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color',[0 .5 0],'linewidth',2)

    % adjusted
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,-My_damaged_mean_all_TEnTip'--My_damaged_mean_all(S3ratios==1),1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'r','linewidth',2)
%     plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,-Mz_damaged_mean_all_TEnTip',2),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'g','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,-My_damaged_mean_steady_TEnTip'--My_damaged_mean_steady(S3ratios==1),1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color',[1 .5 0],'linewidth',2)

    figure(2)
    subplot(2,2,1)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'k','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[.5 .5 .5],'linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'b','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_strokeMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'g','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_devMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'r','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_rotMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'c','linewidth',2)

    % adjusted
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',2),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[.5 .5 .5],'linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqMOD_TEnTip',2),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'b','linewidth',2)

    subplot(2,2,2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_steady_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'k','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_all_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color',[.5 .5 .5],'linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'b','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_strokeMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'g','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_devMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'r','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_rotMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'c','linewidth',2)

    figure(3)
    subplot(2,2,1)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fy_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'k','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fy_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[.5 .5 .5],'linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fy_damaged_mean_freqMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'b','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fy_damaged_mean_strokeMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'g','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fy_damaged_mean_devMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'r','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fy_damaged_mean_rotMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'c','linewidth',2)

    subplot(2,2,2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,-Mz_damaged_mean_steady_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'k','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,-Mz_damaged_mean_all_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color',[.5 .5 .5],'linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,-Mz_damaged_mean_freqMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'b','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,-Mz_damaged_mean_strokeMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'g','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,-Mz_damaged_mean_devMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'r','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,-Mz_damaged_mean_rotMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'c','linewidth',2)

    figure(4)
    subplot(2,2,1)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'k','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'y','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqRotMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[1 .5 0],'linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqRotDevMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'r','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[.5 0 0],'linewidth',2)

%     % adjusted
%     plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',2),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[.5 .5 .5],'linewidth',2)
%     plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqMOD_TEnTip',2),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'b','linewidth',2)

    subplot(2,2,2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_steady_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'k','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'y','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqRotMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color',[1 .5 0],'linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqRotDevMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'r','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_all_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color',[.5 0 0],'linewidth',2)

    figure(5)
    subplot(2,2,1)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'k','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'y','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqStrokeMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[1 .5 0],'linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqStrokeDevMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'r','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[.5 0 0],'linewidth',2)

%     % adjusted
%     plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',2),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[.5 .5 .5],'linewidth',2)
%     plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqMOD_TEnTip',2),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'b','linewidth',2)

    subplot(2,2,2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_steady_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'k','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'y','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqStrokeMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color',[1 .5 0],'linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqStrokeDevMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'r','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_all_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color',[.5 0 0],'linewidth',2)

    figure(6)
    subplot(2,2,1)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'k','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_strokeMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'c','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_strokeDevMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[0 .5 1],'linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_strokeDevRotMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'b','linewidth',2)
    plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[0 0 .5],'linewidth',2)

%     % adjusted
%     plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',2),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[.5 .5 .5],'linewidth',2)
%     plot([min(S2ratios_TEnTip) max(S2ratios_TEnTip)],polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqMOD_TEnTip',2),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]),'b','linewidth',2)

    subplot(2,2,2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_steady_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'k','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_strokeMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'c','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_strokeDevMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'color',[0 .5 1],'linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_strokeDevRotMOD_TEnTip',1),[min(S3ratios_TEnTip) max(S3ratios_TEnTip)]),'b','linewidth',2)
    plot([min(S3ratios_TEnTip) max(S3ratios_TEnTip)],polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_all_TEnTip',1),[min(S3ratios_TEnTip) max(S2ratios_TEnTip)]),'color',[0 0 .5],'linewidth',2)

    %% percentage of weight support and roll equilibrium FREQ FIRST & STROKE LAST
    Fz_steady = -Fz_damaged_mean_steady_TEnTip(1);
    Fz_none = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_freq = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_freqRot = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqRotMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_freqRotDev = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqRotDevMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_freqRotDevStroke = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));

    Fz_freq_perc_frds    = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_freq-Fz_none)
    Fz_rot_perc_frds     = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_freqRot-Fz_freq)
    Fz_dev_perc_frds     = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_freqRotDev-Fz_freqRot)
    Fz_stroke_perc_frds  = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_freqRotDevStroke-Fz_freqRotDev)
    Fz_all_perc_frds     = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_freqRotDevStroke-Fz_none)
    Fz_sum_perc_frds     = Fz_freq_perc_frds + Fz_rot_perc_frds + Fz_dev_perc_frds + Fz_rot_perc_frds
    
    Mx_steady = Mx_damaged_mean_steady_TEnTip(1);
    Mx_none = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_freq = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_freqRot = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqRotMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_freqRotDev = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqRotDevMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_freqRotDevStroke = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));

    Mx_freq_perc_frds    = 100 / (-Mx_none) * (Mx_freq-Mx_none)
    Mx_rot_perc_frds     = 100 / (-Mx_none) * (Mx_freqRot-Mx_freq)
    Mx_dev_perc_frds     = 100 / (-Mx_none) * (Mx_freqRotDev-Mx_freqRot)
    Mx_stroke_perc_frds  = 100 / (-Mx_none) * (Mx_freqRotDevStroke-Mx_freqRotDev)
    Mx_all_perc_frds     = 100 / (-Mx_none) * (Mx_freqRotDevStroke-Mx_none)
    Mx_sum_perc_frds     = Mx_freq_perc_frds + Mx_rot_perc_frds + Mx_dev_perc_frds + Mx_stroke_perc_frds
    
    %% percentage of weight support and roll equilibrium FREQ FIRST & ROT LAST
    Fz_steady = -Fz_damaged_mean_steady_TEnTip(1);
    Fz_none = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_freq = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_freqStroke = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqStrokeMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_freqStrokeDev = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_freqStrokeDevMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_freqStrokeDevRot = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));

    Fz_freq_perc_fsdr    = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_freq-Fz_none)
    Fz_stroke_perc_fsdr  = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_freqStroke-Fz_freq)
    Fz_dev_perc_fsdr     = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_freqStrokeDev-Fz_freqStroke)
    Fz_rot_perc_fsdr     = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_freqStrokeDevRot-Fz_freqStrokeDev)
    Fz_all_perc_fsdr     = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_freqStrokeDevRot-Fz_none)
    Fz_sum_perc_fsdr     = Fz_freq_perc_fsdr + Fz_stroke_perc_fsdr + Fz_dev_perc_fsdr + Fz_rot_perc_fsdr
    
    Mx_steady = Mx_damaged_mean_steady_TEnTip(1);
    Mx_none = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_freq = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_freqStroke = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqStrokeMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_freqStrokeDev = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_freqStrokeDevMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_freqStrokeDevRot = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));

    Mx_freq_perc_fsdr    = 100 / (-Mx_none) * (Mx_freq-Mx_none)
    Mx_stroke_perc_fsdr  = 100 / (-Mx_none) * (Mx_freqStroke-Mx_freq)
    Mx_dev_perc_fsdr     = 100 / (-Mx_none) * (Mx_freqStrokeDev-Mx_freqStroke)
    Mx_rot_perc_fsdr     = 100 / (-Mx_none) * (Mx_freqStrokeDevRot-Mx_freqStrokeDev)
    Mx_all_perc_fsdr     = 100 / (-Mx_none) * (Mx_freqStrokeDevRot-Mx_none)
    Mx_sum_perc_fsdr     = Mx_freq_perc_fsdr + Mx_stroke_perc_fsdr + Mx_dev_perc_fsdr + Mx_rot_perc_fsdr
    
    %% percentage of weight support and roll equilibrium STROKE FIRST
    Fz_steady = -Fz_damaged_mean_steady_TEnTip(1);
    Fz_none = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_stroke = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_strokeMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_strokeDev = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_strokeDevMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_strokeDevRot = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_strokeDevRotMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Fz_strokeDevRotFreq = mean(polyval(polyfit(S2ratios_TEnTip,-Fz_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));

    Fz_stroke_perc_sdrf  = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_stroke-Fz_none)
    Fz_dev_perc_sdrf     = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_strokeDev-Fz_stroke)
    Fz_rot_perc_sdrf     = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_strokeDevRot-Fz_strokeDev)
    Fz_freq_perc_sdrf    = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_strokeDevRotFreq-Fz_strokeDevRot)
    Fz_all_perc_sdrf     = 100* Fz_steady / (Fz_steady-Fz_none) * (Fz_strokeDevRotFreq-Fz_none)
    Fz_sum_perc_sdrf     = Fz_freq_perc_sdrf + Fz_stroke_perc_sdrf + Fz_dev_perc_sdrf + Fz_rot_perc_sdrf
    
    Mx_steady = Mx_damaged_mean_steady_TEnTip(1);
    Mx_none = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_steady_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_stroke = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_strokeMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_strokeDev = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_strokeDevMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_strokeDevRot = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_strokeDevRotMOD_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));
    Mx_strokeDevRotFreq = mean(polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_all_TEnTip',1),[min(S2ratios_TEnTip) max(S2ratios_TEnTip)]));

    Mx_stroke_perc_sdrf  = 100 / (-Mx_none) * (Mx_stroke-Mx_none)
    Mx_dev_perc_sdrf     = 100 / (-Mx_none) * (Mx_strokeDev-Mx_stroke)
    Mx_rot_perc_sdrf     = 100 / (-Mx_none) * (Mx_strokeDevRot-Mx_strokeDev)
    Mx_freq_perc_sdrf    = 100 / (-Mx_none) * (Mx_strokeDevRotFreq-Mx_strokeDevRot)
    Mx_all_perc_sdrf     = 100 / (-Mx_none) * (Mx_strokeDevRotFreq-Mx_none)
    Mx_sum_perc_sdrf     = Mx_freq_perc_sdrf + Mx_stroke_perc_sdrf + Mx_dev_perc_sdrf + Mx_rot_perc_sdrf
    
    %% save plots
    
    mkdir('qsModel_FnM_TEnTipCut')
    cd('qsModel_FnM_TEnTipCut')
    
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
    saveas(gca,['Fz_Mx_WBmodCombies_freqRotDevStroke_TEnTipclip_YnZflip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['Fz_Mx_WBmodCombies_freqRotDevStroke_TEnTipclip_YnZflip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['Fz_Mx_WBmodCombies_freqRotDevStroke_TEnTipclip_YnZflip_asympFit',num2str(freq_asymFitNr),'.svg'])

    figure(5)
    saveas(gca,['Fz_Mx_WBmodCombies_freqStrokeDevRot_TEnTipclip_YnZflip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['Fz_Mx_WBmodCombies_freqStrokeDevRot_TEnTipclip_YnZflip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['Fz_Mx_WBmodCombies_freqStrokeDevRot_TEnTipclip_YnZflip_asympFit',num2str(freq_asymFitNr),'.svg'])

    figure(6)
    saveas(gca,['Fz_Mx_WBmodCombies_strokeFirst_TEnTipclip_YnZflip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['Fz_Mx_WBmodCombies_strokeFirst_TEnTipclip_YnZflip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['Fz_Mx_WBmodCombies_strokeFirst_TEnTipclip_YnZflip_asympFit',num2str(freq_asymFitNr),'.svg'])

    cd ..
    
    
    