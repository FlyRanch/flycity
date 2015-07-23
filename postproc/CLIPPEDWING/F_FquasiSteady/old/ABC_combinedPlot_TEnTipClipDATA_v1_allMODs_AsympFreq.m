clear;
clc;
close all
warning off


freq_asymFitNr = 10;
    

        %% plot TIP cut wing data
    load(['allMODs_TipClip_freqAsym',num2str(freq_asymFitNr),'.mat'])

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
    
    
    