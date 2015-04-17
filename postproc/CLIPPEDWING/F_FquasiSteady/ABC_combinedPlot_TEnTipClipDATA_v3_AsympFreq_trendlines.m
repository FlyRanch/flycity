clear;
clc;
close all
warning off


freq_asymFitNr = 10;
    

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
    
    % ALL MODs VS NO MODs
    figure(1)
    subplot(1,2,1)
    hold on
    plot(S2ratios,Fx_damaged_mean_all,'ok','markersize',10,'markerfacecolor','b')
    plot(S2ratios,Fy_damaged_mean_all,'ok','markersize',10,'markerfacecolor','r')
    plot(S2ratios,Fz_damaged_mean_all,'ok','markersize',10,'markerfacecolor','g')
    
    plot(S2ratios,Fx_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','c')
    plot(S2ratios,Fy_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','m')
    plot(S2ratios,Fz_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','y')
    
    subplot(1,2,2)
    hold on
    plot(S2ratios,Mx_damaged_mean_all,'ok','markersize',10,'markerfacecolor','b')
    plot(S2ratios,My_damaged_mean_all-My_damaged_mean_all(S2ratios==1),'ok','markersize',10,'markerfacecolor','r')
    plot(S2ratios,Mz_damaged_mean_all,'ok','markersize',10,'markerfacecolor','g')
    
    plot(S2ratios,Mx_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','c')
    plot(S2ratios,My_damaged_mean_steady-My_damaged_mean_all(S2ratios==1),'ok','markersize',10,'markerfacecolor','m')
    plot(S2ratios,Mz_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','y')
    
    % compare MODs Fz & Mx
    figure(2)
    subplot(1,2,1)
    hold on
    plot(S2ratios,Fz_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','w')
    plot(S2ratios,Fz_damaged_mean_all,'ok','markersize',10,'markerfacecolor',[.5 .5 .5])
    plot(S2ratios,Fz_damaged_mean_freqMOD,'ok','markersize',10,'markerfacecolor','b')
    plot(S2ratios,Fz_damaged_mean_strokeMOD,'ok','markersize',10,'markerfacecolor','g')
    plot(S2ratios,Fz_damaged_mean_devMOD,'ok','markersize',10,'markerfacecolor','r')
    plot(S2ratios,Fz_damaged_mean_rotMOD,'ok','markersize',10,'markerfacecolor','c')
    
    subplot(1,2,2)
    hold on
    plot(S2ratios,Mx_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','w')
    plot(S2ratios,Mx_damaged_mean_all,'ok','markersize',10,'markerfacecolor',[.5 .5 .5])
    plot(S2ratios,Mx_damaged_mean_freqMOD,'ok','markersize',10,'markerfacecolor','b')
    plot(S2ratios,Mx_damaged_mean_strokeMOD,'ok','markersize',10,'markerfacecolor','g')
    plot(S2ratios,Mx_damaged_mean_devMOD,'ok','markersize',10,'markerfacecolor','r')
    plot(S2ratios,Mx_damaged_mean_rotMOD,'ok','markersize',10,'markerfacecolor','c')
    
    % compare MODs: Fy & Mz
    figure(3)
    subplot(1,2,1)
    hold on
    plot(S2ratios,Fy_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','w')
    plot(S2ratios,Fy_damaged_mean_all,'ok','markersize',10,'markerfacecolor',[.5 .5 .5])
    plot(S2ratios,Fy_damaged_mean_freqMOD,'ok','markersize',10,'markerfacecolor','b')
    plot(S2ratios,Fy_damaged_mean_strokeMOD,'ok','markersize',10,'markerfacecolor','g')
    plot(S2ratios,Fy_damaged_mean_devMOD,'ok','markersize',10,'markerfacecolor','r')
    plot(S2ratios,Fy_damaged_mean_rotMOD,'ok','markersize',10,'markerfacecolor','c')
    
    subplot(1,2,2)
    hold on
    plot(S2ratios,Mz_damaged_mean_steady,'ok','markersize',10,'markerfacecolor','w')
    plot(S2ratios,Mz_damaged_mean_all,'ok','markersize',10,'markerfacecolor',[.5 .5 .5])
    plot(S2ratios,Mz_damaged_mean_freqMOD,'ok','markersize',10,'markerfacecolor','b')
    plot(S2ratios,Mz_damaged_mean_strokeMOD,'ok','markersize',10,'markerfacecolor','g')
    plot(S2ratios,Mz_damaged_mean_devMOD,'ok','markersize',10,'markerfacecolor','r')
    plot(S2ratios,Mz_damaged_mean_rotMOD,'ok','markersize',10,'markerfacecolor','c')
    
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
    

    % ALL MODs VS NO MODs
    figure(1)
    subplot(1,2,1)
    hold on
    plot(S2ratios,Fx_damaged_mean_all,'dk','markersize',10,'markerfacecolor','b')
    plot(S2ratios,Fy_damaged_mean_all,'dk','markersize',10,'markerfacecolor','r')
    plot(S2ratios,Fz_damaged_mean_all,'dk','markersize',10,'markerfacecolor','g')
    
    plot(S2ratios,Fx_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','c')
    plot(S2ratios,Fy_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','m')
    plot(S2ratios,Fz_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','y')
    
    legend('x-Mod','y-Mod','z-Mod','x-steady','y-steady','z-steady','location','E')
    xlabel('S2 ratio')
    ylabel('normalized forces F/mg')
    axis([0.5 1 -1.25 .25])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-1.5:.25:1)
    
    subplot(1,2,2)
    hold on
    plot(S2ratios,Mx_damaged_mean_all,'dk','markersize',10,'markerfacecolor','b')
    plot(S2ratios,My_damaged_mean_all-My_damaged_mean_all(S2ratios==1),'dk','markersize',10,'markerfacecolor','r')
    plot(S2ratios,Mz_damaged_mean_all,'dk','markersize',10,'markerfacecolor','g')
    
    plot(S2ratios,Mx_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','c')
    plot(S2ratios,My_damaged_mean_steady-My_damaged_mean_all(S2ratios==1),'dk','markersize',10,'markerfacecolor','m')
    plot(S2ratios,Mz_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','y')
    
    xlabel('S3 ratio')
    ylabel('normalized torques T/mgl')
    axis([0.5 1 -.1 .15])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-.15:.05:.2)    
    
    % compare MODs Fz & Mx
    
    figure(2)
    subplot(1,2,1)
    hold on
    plot(S2ratios,Fz_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','w')
    plot(S2ratios,Fz_damaged_mean_all,'dk','markersize',10,'markerfacecolor',[.5 .5 .5])
    plot(S2ratios,Fz_damaged_mean_freqMOD,'dk','markersize',10,'markerfacecolor','b')
    plot(S2ratios,Fz_damaged_mean_strokeMOD,'dk','markersize',10,'markerfacecolor','g')
    plot(S2ratios,Fz_damaged_mean_devMOD,'dk','markersize',10,'markerfacecolor','r')
    plot(S2ratios,Fz_damaged_mean_rotMOD,'dk','markersize',10,'markerfacecolor','c')
    
    legend('steady','modified','frequency','stroke','deviation','rotation','location','NE')
    xlabel('S2 ratio')
    ylabel('normalized vertical force Fz/mg')
    axis([0.5 1 -1 -.5])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-1.5:.25:1)
    
    subplot(1,2,2)
    hold on
    plot(S2ratios,Mx_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','w')
    plot(S2ratios,Mx_damaged_mean_all,'dk','markersize',10,'markerfacecolor',[.5 .5 .5])
    plot(S2ratios,Mx_damaged_mean_freqMOD,'dk','markersize',10,'markerfacecolor','b')
    plot(S2ratios,Mx_damaged_mean_strokeMOD,'dk','markersize',10,'markerfacecolor','g')
    plot(S2ratios,Mx_damaged_mean_devMOD,'dk','markersize',10,'markerfacecolor','r')
    plot(S2ratios,Mx_damaged_mean_rotMOD,'dk','markersize',10,'markerfacecolor','c')
    
    xlabel('S3 ratio')
    ylabel('normalized roll torque Tx/mgl')
    axis([0.5 1 0 .25])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-.25:.25:.3)    

    % compare MODs: Fy & Mz
    
    figure(3)
    subplot(1,2,1)
    hold on
    plot(S2ratios,Fy_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','w')
    plot(S2ratios,Fy_damaged_mean_all,'dk','markersize',10,'markerfacecolor',[.5 .5 .5])
    plot(S2ratios,Fy_damaged_mean_freqMOD,'dk','markersize',10,'markerfacecolor','b')
    plot(S2ratios,Fy_damaged_mean_strokeMOD,'dk','markersize',10,'markerfacecolor','g')
    plot(S2ratios,Fy_damaged_mean_devMOD,'dk','markersize',10,'markerfacecolor','r')
    plot(S2ratios,Fy_damaged_mean_rotMOD,'dk','markersize',10,'markerfacecolor','c')
    
    legend('steady','modified','frequency','stroke','deviation','rotation','location','SE')
    xlabel('S2 ratio')
    ylabel('normalized sideways force Fy/mg')
    axis([0.5 1 -.5 0])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-1.5:.25:1)
    
    subplot(1,2,2)
    hold on
    plot(S2ratios,Mz_damaged_mean_steady,'dk','markersize',10,'markerfacecolor','w')
    plot(S2ratios,Mz_damaged_mean_all,'dk','markersize',10,'markerfacecolor',[.5 .5 .5])
    plot(S2ratios,Mz_damaged_mean_freqMOD,'dk','markersize',10,'markerfacecolor','b')
    plot(S2ratios,Mz_damaged_mean_strokeMOD,'dk','markersize',10,'markerfacecolor','g')
    plot(S2ratios,Mz_damaged_mean_devMOD,'dk','markersize',10,'markerfacecolor','r')
    plot(S2ratios,Mz_damaged_mean_rotMOD,'dk','markersize',10,'markerfacecolor','c')
    
    xlabel('S3 ratio')
    ylabel('normalized yaw torque Tz/mgl')
    axis([0.5 1 -.15 .1])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-.2:.05:.3)    
    
    
    %% plot trendlines
    figure(1)
    subplot(1,2,1)
%     plot(sort(S2ratios_TEnTip),polyval(polyfit(S2ratios_TEnTip,Fx_damaged_mean_all_TEnTip',1),sort(S2ratios_TEnTip)),'b','linewidth',2)
    plot(sort(S2ratios_TEnTip),polyval(polyfit(S2ratios_TEnTip,Fy_damaged_mean_all_TEnTip',1),sort(S2ratios_TEnTip)),'r','linewidth',2)
    plot(sort(S2ratios_TEnTip),polyval(polyfit(S2ratios_TEnTip,Fz_damaged_mean_all_TEnTip',1),sort(S2ratios_TEnTip)),'g','linewidth',2)

    plot(sort(S2ratios_TEnTip),polyval(polyfit(S2ratios_TEnTip,Fx_damaged_mean_steady_TEnTip',1),sort(S2ratios_TEnTip)),'c','linewidth',2)
    plot(sort(S2ratios_TEnTip),polyval(polyfit(S2ratios_TEnTip,Fy_damaged_mean_steady_TEnTip',1),sort(S2ratios_TEnTip)),'m','linewidth',2)
    plot(sort(S2ratios_TEnTip),polyval(polyfit(S2ratios_TEnTip,Fz_damaged_mean_steady_TEnTip',1),sort(S2ratios_TEnTip)),'y','linewidth',2)

    % adjusted
    plot(sort(S2ratios_TEnTip),polyval(polyfit(S2ratios_TEnTip,Fz_damaged_mean_all_TEnTip',2),sort(S2ratios_TEnTip)),'g','linewidth',2)

    subplot(1,2,2)
    plot(sort(S3ratios_TEnTip),polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_all_TEnTip',1),sort(S3ratios_TEnTip)),'b','linewidth',2)
%     plot(sort(S3ratios_TEnTip),polyval(polyfit(S3ratios_TEnTip,My_damaged_mean_all_TEnTip',1),sort(S3ratios_TEnTip)),'r','linewidth',2)
%     plot(sort(S3ratios_TEnTip),polyval(polyfit(S3ratios_TEnTip,Mz_damaged_mean_all_TEnTip',1),sort(S3ratios_TEnTip)),'g','linewidth',2)

    plot(sort(S3ratios_TEnTip),polyval(polyfit(S3ratios_TEnTip,Mx_damaged_mean_steady_TEnTip',1),sort(S3ratios_TEnTip)),'c','linewidth',2)
%     plot(sort(S3ratios_TEnTip),polyval(polyfit(S3ratios_TEnTip,My_damaged_mean_steady_TEnTip',1),sort(S3ratios_TEnTip)),'m','linewidth',2)
    plot(sort(S3ratios_TEnTip),polyval(polyfit(S3ratios_TEnTip,Mz_damaged_mean_steady_TEnTip',1),sort(S3ratios_TEnTip)),'y','linewidth',2)

    % adjusted
    plot(sort(S3ratios_TEnTip),polyval(polyfit(S3ratios_TEnTip,My_damaged_mean_all_TEnTip'-My_damaged_mean_all(S2ratios==1),1),sort(S3ratios_TEnTip)),'r','linewidth',2)
    plot(sort(S3ratios_TEnTip),polyval(polyfit(S3ratios_TEnTip,Mz_damaged_mean_all_TEnTip',2),sort(S3ratios_TEnTip)),'g','linewidth',2)
    plot(sort(S3ratios_TEnTip),polyval(polyfit(S3ratios_TEnTip,My_damaged_mean_steady_TEnTip'-My_damaged_mean_steady(S2ratios==1),1),sort(S3ratios_TEnTip)),'m','linewidth',2)

    %% save plots
    
    mkdir('qsModel_FnM_TEnTipCut')
    cd('qsModel_FnM_TEnTipCut')
    
    figure(1)
    saveas(gca,['FnM_WBmod_TEnTipClip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['FnM_WBmod_TEnTipClip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['FnM_WBmod_TEnTipClip_asympFit',num2str(freq_asymFitNr),'.svg'])

    figure(2)
    saveas(gca,['Fz_Mx_WBmodComponents_TEnTipClip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['Fz_Mx_WBmodComponents_TEnTipClip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['Fz_Mx_WBmodComponents_TEnTipClip_asympFit',num2str(freq_asymFitNr),'.svg'])

    figure(3)
    saveas(gca,['Fy_Mz_WBmodComponents_TEnTipclip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['Fy_Mz_WBmodComponents_TEnTipclip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['Fy_Mz_WBmodComponents_TEnTipclip_asympFit',num2str(freq_asymFitNr),'.svg'])
    


    cd ..
    
    
    