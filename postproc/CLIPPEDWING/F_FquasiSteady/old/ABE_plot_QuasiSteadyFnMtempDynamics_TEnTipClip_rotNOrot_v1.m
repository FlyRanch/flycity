clear;
clc;
close all
warning off

%% TEclip ROT on
clear
figure
freq_asymFitNr = 10;
load(['allMODs_TEclip_freqAsym',num2str(freq_asymFitNr),'_roton1.mat'])

    subplot(2,1,1)
    hold on
    
    Fz_mean_steady_nondamaged = -mean(Fz_intactwing_steady(:,1)+Fz_intactwing_steady(:,end));
    Fz_mean_all = -mean(Fz_intactwing_all(:,end)+Fz_damagedwing_all(:,end));
%     plot([0 1],[Fz_mean_steady_nondamaged Fz_mean_steady_nondamaged],'--','color',[.5 .5 .5],'linewidth',2)
    plot([0 1],[Fz_mean_all Fz_mean_all],'--k','linewidth',2)
    
    plot(t_norm,-Fz_intactwing_steady(:,end), '-k','linewidth',2)
    plot(t_norm,-Fz_intactwing_all(:,end), '-','color',[0 0 .5],'linewidth',2)
    plot(t_norm,-Fz_damagedwing_all(:,end),'-','color',[.5 0 0],'linewidth',2)

%     legend('non-damaged wing steady WB','damaged wing steady WB','non-damaged wing modified WB','damaged wing modified WB','','location','SW')
    xlabel('wingbeat cycle')
    ylabel('normalized vertical force Fz/mg')
    axis([0 1 -1 2])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-4:1)
    
    subplot(2,1,2)
    hold on
    
    Mx_mean_steady_nondamaged = mean(Mx_intactwing_steady(:,1)-Mx_intactwing_steady(:,end));
    Mx_mean_all = mean(Mx_intactwing_all(:,end)+Mx_damagedwing_all(:,end));
%     plot([0 1],[Mx_mean_steady_nondamaged Mx_mean_steady_nondamaged],'--','color',[.5 .5 .5],'linewidth',2)
    plot([0 1],[Mx_mean_all Mx_mean_all],'--k','linewidth',2)
    
    plot(t_norm,Mx_intactwing_steady(:,end), '-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,-Mx_intactwing_steady(:,end), '-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,Mx_intactwing_all(:,end), '-','color',[0 0 .5],'linewidth',2)
    plot(t_norm,Mx_damagedwing_all(:,end),'-','color',[.5 0 0],'linewidth',2)
    
    xlabel('wingbeat cycle')
    ylabel('normalized roll torque Tx/mgl')
    axis([0 1 -1 1])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-1:1)    

%% TEclip ROT off
clear
freq_asymFitNr = 10;
load(['allMODs_TEclip_freqAsym',num2str(freq_asymFitNr),'_roton0.mat'])

    subplot(2,1,1)
    hold on
    
    Fz_mean_steady_nondamaged = -mean(Fz_intactwing_steady(:,1)+Fz_intactwing_steady(:,end));
    Fz_mean_all = -mean(Fz_intactwing_all(:,end)+Fz_damagedwing_all(:,end));
%     plot([0 1],[Fz_mean_steady_nondamaged Fz_mean_steady_nondamaged],'--','color',[.5 .5 .5],'linewidth',2)
    plot([0 1],[Fz_mean_all Fz_mean_all],'--','color',[.75 .75 .75],'linewidth',2)
    
    plot(t_norm,-Fz_intactwing_steady(:,end), '-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,-Fz_intactwing_all(:,end), '-','color',[0 .5 1],'linewidth',2)
    plot(t_norm,-Fz_damagedwing_all(:,end),'-','color',[1 .5 0],'linewidth',2)

%     legend('non-damaged wing steady WB','damaged wing steady WB','non-damaged wing modified WB','damaged wing modified WB','','location','SW')
    xlabel('wingbeat cycle')
    ylabel('normalized vertical force Fz/mg')
    axis([0 1 -1 2])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-4:1)
    
    subplot(2,1,2)
    hold on
    
    Mx_mean_steady_nondamaged = mean(Mx_intactwing_steady(:,1)-Mx_intactwing_steady(:,end));
    Mx_mean_all = mean(Mx_intactwing_all(:,end)+Mx_damagedwing_all(:,end));
%     plot([0 1],[Mx_mean_steady_nondamaged Mx_mean_steady_nondamaged],'--','color',[.5 .5 .5],'linewidth',2)
    plot([0 1],[Mx_mean_all Mx_mean_all],'--','color',[.75 .75 .75],'linewidth',2)
    
%     plot(t_norm,Mx_intactwing_steady(:,end), '-','color',[.5 .5 .5],'linewidth',2)
%     plot(t_norm,-Mx_intactwing_steady(:,end), '-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,Mx_intactwing_all(:,end), '-','color',[0 .5 1],'linewidth',2)
    plot(t_norm,Mx_damagedwing_all(:,end),'-','color',[1 .5 0],'linewidth',2)
    
    xlabel('wingbeat cycle')
    ylabel('normalized roll torque Tx/mgl')
    axis([0 1 -1 1])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-1:1)    

%% save figs
mkdir('qsModel_FnM_TempDynamics_TEnTipCut_rotNOrot')
cd('qsModel_FnM_TempDynamics_TEnTipCut_rotNOrot')

saveas(gca,['FzFlip_Mx_TempDynamics_TEclip_rotNOrot_asympFit',num2str(freq_asymFitNr),'.fig'])
saveas(gca,['FzFlip_Mx_TempDynamics_TEclip_rotNOrot_asympFit',num2str(freq_asymFitNr),'.png'])
plot2svg(  ['FzFlip_Mx_TempDynamics_TEclip_rotNOrot_asympFit',num2str(freq_asymFitNr),'.svg'])

cd ..

%% TipClip ROT on
figure
clear
freq_asymFitNr = 10;
load(['allMODs_TipClip_freqAsym',num2str(freq_asymFitNr),'_roton1.mat'])

    subplot(2,1,1)
    hold on
    
    Fz_mean_steady_nondamaged = -mean(Fz_intactwing_steady(:,1)+Fz_intactwing_steady(:,end));
    Fz_mean_all = -mean(Fz_intactwing_all(:,end)+Fz_damagedwing_all(:,end));
%     plot([0 1],[Fz_mean_steady_nondamaged Fz_mean_steady_nondamaged],'--','color',[.5 .5 .5],'linewidth',2)
    plot([0 1],[Fz_mean_all Fz_mean_all],'--k','linewidth',2)
    
    plot(t_norm,-Fz_intactwing_steady(:,end), '-k','linewidth',2)
    plot(t_norm,-Fz_intactwing_all(:,end), '-','color',[0 0 .5],'linewidth',2)
    plot(t_norm,-Fz_damagedwing_all(:,end),'-','color',[.5 0 0],'linewidth',2)

%     legend('non-damaged wing steady WB','damaged wing steady WB','non-damaged wing modified WB','damaged wing modified WB','','location','SW')
    xlabel('wingbeat cycle')
    ylabel('normalized vertical force Fz/mg')
    axis([0 1 -1 2])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-4:1)
    
    subplot(2,1,2)
    hold on
    
    Mx_mean_steady_nondamaged = mean(Mx_intactwing_steady(:,1)-Mx_intactwing_steady(:,end));
    Mx_mean_all = mean(Mx_intactwing_all(:,end)+Mx_damagedwing_all(:,end));
%     plot([0 1],[Mx_mean_steady_nondamaged Mx_mean_steady_nondamaged],'--','color',[.5 .5 .5],'linewidth',2)
    plot([0 1],[Mx_mean_all Mx_mean_all],'--k','linewidth',2)
    
    plot(t_norm,Mx_intactwing_steady(:,end), '-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,-Mx_intactwing_steady(:,end), '-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,Mx_intactwing_all(:,end), '-','color',[0 0 .5],'linewidth',2)
    plot(t_norm,Mx_damagedwing_all(:,end),'-','color',[.5 0 0],'linewidth',2)
    
    xlabel('wingbeat cycle')
    ylabel('normalized roll torque Tx/mgl')
    axis([0 1 -1 1])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-1:1)    

%% TipClip ROT off
clear
freq_asymFitNr = 10;
load(['allMODs_TipClip_freqAsym',num2str(freq_asymFitNr),'_roton0.mat'])

    subplot(2,1,1)
    hold on
    
    Fz_mean_steady_nondamaged = -mean(Fz_intactwing_steady(:,1)+Fz_intactwing_steady(:,end));
    Fz_mean_all = -mean(Fz_intactwing_all(:,end)+Fz_damagedwing_all(:,end));
%     plot([0 1],[Fz_mean_steady_nondamaged Fz_mean_steady_nondamaged],'--','color',[.5 .5 .5],'linewidth',2)
    plot([0 1],[Fz_mean_all Fz_mean_all],'--','color',[.75 .75 .75],'linewidth',2)
    
    plot(t_norm,-Fz_intactwing_steady(:,end), '-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,-Fz_intactwing_all(:,end), '-','color',[0 .5 1],'linewidth',2)
    plot(t_norm,-Fz_damagedwing_all(:,end),'-','color',[1 .5 0],'linewidth',2)

%     legend('non-damaged wing steady WB','damaged wing steady WB','non-damaged wing modified WB','damaged wing modified WB','','location','SW')
    xlabel('wingbeat cycle')
    ylabel('normalized vertical force Fz/mg')
    axis([0 1 -1 2])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-4:1)
    
    subplot(2,1,2)
    hold on
    
    Mx_mean_steady_nondamaged = mean(Mx_intactwing_steady(:,1)-Mx_intactwing_steady(:,end));
    Mx_mean_all = mean(Mx_intactwing_all(:,end)+Mx_damagedwing_all(:,end));
%     plot([0 1],[Mx_mean_steady_nondamaged Mx_mean_steady_nondamaged],'--','color',[.5 .5 .5],'linewidth',2)
    plot([0 1],[Mx_mean_all Mx_mean_all],'--','color',[.75 .75 .75],'linewidth',2)
    
%     plot(t_norm,Mx_intactwing_steady(:,end), '-','color',[.5 .5 .5],'linewidth',2)
%     plot(t_norm,-Mx_intactwing_steady(:,end), '-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,Mx_intactwing_all(:,end), '-','color',[0 .5 1],'linewidth',2)
    plot(t_norm,Mx_damagedwing_all(:,end),'-','color',[1 .5 0],'linewidth',2)
    
    xlabel('wingbeat cycle')
    ylabel('normalized roll torque Tx/mgl')
    axis([0 1 -1 1])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-1:1)    

%% save figs
mkdir('qsModel_FnM_TempDynamics_TEnTipCut_rotNOrot')
cd('qsModel_FnM_TempDynamics_TEnTipCut_rotNOrot')

saveas(gca,['FzFlip_Mx_TempDynamics_TipClip_rotNOrot_asympFit',num2str(freq_asymFitNr),'.fig'])
saveas(gca,['FzFlip_Mx_TempDynamics_TipClip_rotNOrot_asympFit',num2str(freq_asymFitNr),'.png'])
plot2svg(  ['FzFlip_Mx_TempDynamics_TipClip_rotNOrot_asympFit',num2str(freq_asymFitNr),'.svg'])

cd ..
    