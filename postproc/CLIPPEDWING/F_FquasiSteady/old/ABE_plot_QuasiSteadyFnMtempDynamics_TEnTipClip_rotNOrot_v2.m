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
    
    plot([0 1],[-mean(Fz_intactwing_steady(:,end)) -mean(Fz_intactwing_steady(:,end))],'--','color',[0 0 0],'linewidth',2)
    plot([0 1],[-mean(Fz_intactwing_all(:,end))    -mean(Fz_intactwing_all(:,end))],   '--','color',[0 0 1],'linewidth',2)
    plot([0 1],[-mean(Fz_damagedwing_all(:,end))   -mean(Fz_damagedwing_all(:,end))],  '--','color',[1 0 0],'linewidth',2)
    
    plot(t_norm,-Fz_intactwing_steady(:,end), '-','color',[0 0 0],'linewidth',2)
    plot(t_norm,-Fz_intactwing_all(:,end),    '-','color',[0 0 1],'linewidth',2)
    plot(t_norm,-Fz_damagedwing_all(:,end),   '-','color',[1 0 0],'linewidth',2)

%     legend('non-damaged wing steady WB','damaged wing steady WB','non-damaged wing modified WB','damaged wing modified WB','','location','SW')
    xlabel('wingbeat cycle')
    ylabel('normalized vertical force Fz/mg')
    axis([0 1 -1 2])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-4:1)
    
    subplot(2,1,2)
    hold on
    
    plot([0 1],[mean(Mx_intactwing_steady(:,end)) mean(Mx_intactwing_steady(:,end))],'--','color',[0 0 0],'linewidth',2)
    plot([0 1],[-mean(Mx_intactwing_steady(:,end)) -mean(Mx_intactwing_steady(:,end))],'--','color',[0 0 0],'linewidth',2)
    plot([0 1],[mean(Mx_intactwing_all(:,end))    mean(Mx_intactwing_all(:,end))],   '--','color',[0 0 1],'linewidth',2)
    plot([0 1],[mean(Mx_damagedwing_all(:,end))   mean(Mx_damagedwing_all(:,end))],  '--','color',[1 0 0],'linewidth',2)

    plot(t_norm,Mx_intactwing_steady(:,end), '-','color',[0 0 0],'linewidth',2)
    plot(t_norm,-Mx_intactwing_steady(:,end), '-','color',[0 0 0],'linewidth',2)
    plot(t_norm,Mx_intactwing_all(:,end),    '-','color',[0 0 1],'linewidth',2)
    plot(t_norm,Mx_damagedwing_all(:,end),   '-','color',[1 0 0],'linewidth',2)
    
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
    
    plot([0 1],[-mean(Fz_intactwing_steady(:,end)) -mean(Fz_intactwing_steady(:,end))],'--','color',[.5 .5 .5],'linewidth',2)
    plot([0 1],[-mean(Fz_intactwing_all(:,end))    -mean(Fz_intactwing_all(:,end))],   '--','color',[0 1 1],'linewidth',2)
    plot([0 1],[-mean(Fz_damagedwing_all(:,end))   -mean(Fz_damagedwing_all(:,end))],  '--','color',[1 .5 0],'linewidth',2)
    
    plot(t_norm,-Fz_intactwing_steady(:,end), '-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,-Fz_intactwing_all(:,end),    '-','color',[0 1 1],'linewidth',2)
    plot(t_norm,-Fz_damagedwing_all(:,end),   '-','color',[1 .5 0],'linewidth',2)

%     legend('non-damaged wing steady WB','damaged wing steady WB','non-damaged wing modified WB','damaged wing modified WB','','location','SW')
    xlabel('wingbeat cycle')
    ylabel('normalized vertical force Fz/mg')
    axis([0 1 -1 2])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-4:1)
    
    subplot(2,1,2)
    hold on
    
    plot([0 1],[mean(Mx_intactwing_steady(:,end)) mean(Mx_intactwing_steady(:,end))],'--','color',[.5 .5 .5],'linewidth',2)
    plot([0 1],[-mean(Mx_intactwing_steady(:,end)) -mean(Mx_intactwing_steady(:,end))],'--','color',[.5 .5 .5],'linewidth',2)
    plot([0 1],[mean(Mx_intactwing_all(:,end))    mean(Mx_intactwing_all(:,end))],   '--','color',[0 1 1],'linewidth',2)
    plot([0 1],[mean(Mx_damagedwing_all(:,end))   mean(Mx_damagedwing_all(:,end))],  '--','color',[1 .5 0],'linewidth',2)

    plot(t_norm,Mx_intactwing_steady(:,end), '-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,-Mx_intactwing_steady(:,end), '-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,Mx_intactwing_all(:,end),    '-','color',[0 1 1],'linewidth',2)
    plot(t_norm,Mx_damagedwing_all(:,end),   '-','color',[1 .5 0],'linewidth',2)

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
    
    plot([0 1],[-mean(Fz_intactwing_steady(:,end)) -mean(Fz_intactwing_steady(:,end))],'--','color',[0 0 0],'linewidth',2)
    plot([0 1],[-mean(Fz_intactwing_all(:,end))    -mean(Fz_intactwing_all(:,end))],   '--','color',[0 0 1],'linewidth',2)
    plot([0 1],[-mean(Fz_damagedwing_all(:,end))   -mean(Fz_damagedwing_all(:,end))],  '--','color',[1 0 0],'linewidth',2)
    
    plot(t_norm,-Fz_intactwing_steady(:,end), '-','color',[0 0 0],'linewidth',2)
    plot(t_norm,-Fz_intactwing_all(:,end),    '-','color',[0 0 1],'linewidth',2)
    plot(t_norm,-Fz_damagedwing_all(:,end),   '-','color',[1 0 0],'linewidth',2)

%     legend('non-damaged wing steady WB','damaged wing steady WB','non-damaged wing modified WB','damaged wing modified WB','','location','SW')
    xlabel('wingbeat cycle')
    ylabel('normalized vertical force Fz/mg')
    axis([0 1 -1 2])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-4:1)
    
    subplot(2,1,2)
    hold on
    
    plot([0 1],[mean(Mx_intactwing_steady(:,end)) mean(Mx_intactwing_steady(:,end))],'--','color',[0 0 0],'linewidth',2)
    plot([0 1],[-mean(Mx_intactwing_steady(:,end)) -mean(Mx_intactwing_steady(:,end))],'--','color',[0 0 0],'linewidth',2)
    plot([0 1],[mean(Mx_intactwing_all(:,end))    mean(Mx_intactwing_all(:,end))],   '--','color',[0 0 1],'linewidth',2)
    plot([0 1],[mean(Mx_damagedwing_all(:,end))   mean(Mx_damagedwing_all(:,end))],  '--','color',[1 0 0],'linewidth',2)

    plot(t_norm,Mx_intactwing_steady(:,end), '-','color',[0 0 0],'linewidth',2)
    plot(t_norm,-Mx_intactwing_steady(:,end), '-','color',[0 0 0],'linewidth',2)
    plot(t_norm,Mx_intactwing_all(:,end),    '-','color',[0 0 1],'linewidth',2)
    plot(t_norm,Mx_damagedwing_all(:,end),   '-','color',[1 0 0],'linewidth',2)
    
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
    
    plot([0 1],[-mean(Fz_intactwing_steady(:,end)) -mean(Fz_intactwing_steady(:,end))],'--','color',[.5 .5 .5],'linewidth',2)
    plot([0 1],[-mean(Fz_intactwing_all(:,end))    -mean(Fz_intactwing_all(:,end))],   '--','color',[0 1 1],'linewidth',2)
    plot([0 1],[-mean(Fz_damagedwing_all(:,end))   -mean(Fz_damagedwing_all(:,end))],  '--','color',[1 .5 0],'linewidth',2)
    
    plot(t_norm,-Fz_intactwing_steady(:,end), '-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,-Fz_intactwing_all(:,end),    '-','color',[0 1 1],'linewidth',2)
    plot(t_norm,-Fz_damagedwing_all(:,end),   '-','color',[1 .5 0],'linewidth',2)

%     legend('non-damaged wing steady WB','damaged wing steady WB','non-damaged wing modified WB','damaged wing modified WB','','location','SW')
    xlabel('wingbeat cycle')
    ylabel('normalized vertical force Fz/mg')
    axis([0 1 -1 2])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-4:1)
    
    subplot(2,1,2)
    hold on
    
    plot([0 1],[mean(Mx_intactwing_steady(:,end)) mean(Mx_intactwing_steady(:,end))],'--','color',[.5 .5 .5],'linewidth',2)
    plot([0 1],[-mean(Mx_intactwing_steady(:,end)) -mean(Mx_intactwing_steady(:,end))],'--','color',[.5 .5 .5],'linewidth',2)
    plot([0 1],[mean(Mx_intactwing_all(:,end))    mean(Mx_intactwing_all(:,end))],   '--','color',[0 1 1],'linewidth',2)
    plot([0 1],[mean(Mx_damagedwing_all(:,end))   mean(Mx_damagedwing_all(:,end))],  '--','color',[1 .5 0],'linewidth',2)

    plot(t_norm,Mx_intactwing_steady(:,end), '-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,-Mx_intactwing_steady(:,end), '-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,Mx_intactwing_all(:,end),    '-','color',[0 1 1],'linewidth',2)
    plot(t_norm,Mx_damagedwing_all(:,end),   '-','color',[1 .5 0],'linewidth',2)

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
    