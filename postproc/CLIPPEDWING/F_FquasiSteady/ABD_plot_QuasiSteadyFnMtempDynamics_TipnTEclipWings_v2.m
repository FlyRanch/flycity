clear;
clc;
close all
warning off

%% Tip clip    
load(['allMODs_TipClip_freqAsym10.mat'])
    
    figure(1)
    subplot(2,1,1)
    hold on
    
    Fz_mean_steady_nondamaged = mean(Fz_intactwing_steady(:,1)+Fz_intactwing_steady(:,end));
    Fz_mean_steady = mean(Fz_intactwing_steady(:,1)+Fz_damagedwing_steady(:,end));
    Fz_mean_all = mean(Fz_intactwing_all(:,1)+Fz_damagedwing_all(:,end));
    plot([0 1],[Fz_mean_steady_nondamaged Fz_mean_steady_nondamaged],'--k','linewidth',2)
    plot([0 1],[Fz_mean_steady Fz_mean_steady],'--','color',[.5 .5 .5],'linewidth',2)
    plot([0 1],[Fz_mean_all Fz_mean_all],'--r','linewidth',2)
    
    plot(t_norm,Fz_intactwing_steady(:,1),'-k','linewidth',2)
    plot(t_norm,Fz_damagedwing_steady(:,end),'-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,Fz_intactwing_all(:,end),'-','color',[1 .5 0],'linewidth',2)
    plot(t_norm,Fz_damagedwing_all(:,end),'-r','linewidth',2)

%     legend('non-damaged wing steady WB','damaged wing steady WB','non-damaged wing modified WB','damaged wing modified WB','','location','SW')
    xlabel('wingbeat cycle')
    ylabel('normalized vertical force Fz/mg')
    axis([0 1 -2 1])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-4:1)
    
    subplot(2,1,2)
    hold on
    
    Mx_mean_steady = mean(Mx_intactwing_steady(:,1)+Mx_damagedwing_steady(:,end));
    Mx_mean_all = mean(Mx_intactwing_all(:,1)+Mx_damagedwing_all(:,end));
    plot([0 1],[0 0],'--k','linewidth',2)
    plot([0 1],[Mx_mean_steady Mx_mean_steady],'--','color',[.5 .5 .5],'linewidth',2)
    plot([0 1],[Mx_mean_all Mx_mean_all],'--r','linewidth',2)

    plot(t_norm,Mx_intactwing_steady(:,1),'-k','linewidth',2)
    plot(t_norm,Mx_damagedwing_steady(:,end),'-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,Mx_intactwing_all(:,end),'-','color',[1 .5 0],'linewidth',2)
    plot(t_norm,Mx_damagedwing_all(:,end),'-r','linewidth',2)
    
    xlabel('wingbeat cycle')
    ylabel('normalized roll torque Tx/mgl')
    axis([0 1 -1 1])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-1:1)    

    
%% TE clip    
load(['allMODs_TEclip_freqAsym10.mat'])
    
%     figure(2)
    subplot(2,1,1)
    hold on
    
    Fz_mean_steady_nondamaged = mean(Fz_intactwing_steady(:,1)+Fz_intactwing_steady(:,end));
    Fz_mean_steady = mean(Fz_intactwing_steady(:,1)+Fz_damagedwing_steady(:,end));
    Fz_mean_all = mean(Fz_intactwing_all(:,1)+Fz_damagedwing_all(:,end));
%     plot([0 1],[Fz_mean_steady_nondamaged Fz_mean_steady_nondamaged],'--k','linewidth',2)
%     plot([0 1],[Fz_mean_steady Fz_mean_steady],'--','color',[.5 .5 .5],'linewidth',2)
    plot([0 1],[Fz_mean_all Fz_mean_all],'--b','linewidth',2)
    
%     plot(t_norm,Fz_intactwing_steady(:,1),'-k','linewidth',2)
%     plot(t_norm,Fz_damagedwing_steady(:,end),'-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,Fz_intactwing_all(:,end),'-c','linewidth',2)
    plot(t_norm,Fz_damagedwing_all(:,end),'-b','linewidth',2)

%     legend('non-damaged wing steady WB','damaged wing steady WB','non-damaged wing modified WB','damaged wing modified WB','','location','SW')
    xlabel('wingbeat cycle')
    ylabel('normalized vertical force Fz/mg')
    axis([0 1 -2 1])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-4:1)
    
    subplot(2,1,2)
    hold on
    
    Mx_mean_steady = mean(Mx_intactwing_steady(:,1)+Mx_damagedwing_steady(:,end));
    Mx_mean_all = mean(Mx_intactwing_all(:,1)+Mx_damagedwing_all(:,end));
%     plot([0 1],[0 0],'--k','linewidth',2)
%     plot([0 1],[Mx_mean_steady Mx_mean_steady],'--','color',[.5 .5 .5],'linewidth',2)
    plot([0 1],[Mx_mean_all Mx_mean_all],'--b','linewidth',2)

%     plot(t_norm,Mx_intactwing_steady(:,1),'-k','linewidth',2)
%     plot(t_norm,Mx_damagedwing_steady(:,end),'-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,Mx_intactwing_all(:,end),'-c','linewidth',2)
    plot(t_norm,Mx_damagedwing_all(:,end),'-b','linewidth',2)
    
    xlabel('wingbeat cycle')
    ylabel('normalized roll torque Tx/mgl')
    axis([0 1 -1 1])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-1:1)    

%% save figs
mkdir('qsModel_FnM_TempDynamics_TipCut')
cd('qsModel_FnM_TempDynamics_TipCut')
% 
% figure(1)
%     saveas(gca,['FzMx_TempDynamics_TipClip_asympFit',num2str(freq_asymFitNr),'.fig'])
%     saveas(gca,['FzMx_TempDynamics_TipClip_asympFit',num2str(freq_asymFitNr),'.png'])
%     plot2svg(['FzMx_TempDynamics_TipClip_asympFit',num2str(freq_asymFitNr),'.svg'])
% 
% figure(2)
    saveas(gca,['FzMx_TempDynamics_TipTEclip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['FzMx_TempDynamics_TipTEclip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['FzMx_TempDynamics_TipTEclip_asympFit',num2str(freq_asymFitNr),'.svg'])

    cd ..
    