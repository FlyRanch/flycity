clear;
clc;
close all
warning off


freq_asymFitNr = 10;
peakloc = 0;

%% Tip clip    
% load(['allMODs_TipClip_freqAsym',num2str(freq_asymFitNr),'_peakloc',num2str(peakloc),'.mat'])
load(['allMODs_TipClip_freqAsym',num2str(freq_asymFitNr),'_roton1.mat'])
    
    figure(1)
    subplot(2,1,1)
    hold on
    
    Fz_mean_steady_nondamaged = -mean(Fz_intactwing_steady(:,1)+Fz_intactwing_steady(:,end));
    Fz_mean_steady = -mean(Fz_intactwing_steady(:,1)+Fz_damagedwing_steady(:,end));
    Fz_mean_all = -mean(Fz_intactwing_all(:,1)+Fz_damagedwing_all(:,end));
    plot([0 1],[Fz_mean_steady_nondamaged Fz_mean_steady_nondamaged],'--','color',[0 0 .5],'linewidth',2)
    plot([0 1],[Fz_mean_steady Fz_mean_steady],'--','color',[0 .5 1],'linewidth',2)
    plot([0 1],[Fz_mean_all Fz_mean_all],'--r','linewidth',2)
    plot([0 1],[0 0],'-k','linewidth',1)
    
    plot(t_norm,-Fz_intactwing_steady(:,1),'-','color',[0 0 .5],'linewidth',2)
    plot(t_norm,-Fz_damagedwing_steady(:,end),'-','color',[0 .5 1],'linewidth',2)
    plot(t_norm,-Fz_intactwing_all(:,end),'-','color',[1 .5 0],'linewidth',2)
    plot(t_norm,-Fz_damagedwing_all(:,end),'-r','linewidth',2)

%     legend('non-damaged wing steady WB','damaged wing steady WB','non-damaged wing modified WB','damaged wing modified WB','','location','SW')
    xlabel('wingbeat cycle')
    ylabel('normalized vertical force Fz/mg')
    axis([0 1 -1 2])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-4:1)
    
    subplot(2,1,2)
    hold on
    
    Mx_mean_steady = mean(Mx_intactwing_steady(:,1)+Mx_damagedwing_steady(:,end));
    Mx_mean_all = mean(Mx_intactwing_all(:,1)+Mx_damagedwing_all(:,end));
    plot([0 1],[0 0],'--','color',[0 0 .5],'linewidth',2)
    plot([0 1],[Mx_mean_steady Mx_mean_steady],'--','color',[0 .5 1],'linewidth',2)
    plot([0 1],[Mx_mean_all Mx_mean_all],'--r','linewidth',2)

    plot(t_norm,Mx_intactwing_steady(:,1),'-','color',[0 0 .5],'linewidth',2)
    plot(t_norm,Mx_damagedwing_steady(:,end),'-','color',[0 .5 1],'linewidth',2)
    plot(t_norm,Mx_intactwing_all(:,end),'-','color',[1 .5 0],'linewidth',2)
    plot(t_norm,Mx_damagedwing_all(:,end),'-r','linewidth',2)
    
    xlabel('wingbeat cycle')
    ylabel('normalized roll torque Tx/mgl')
    axis([0 1 -1 1])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-1:1)    
    
    figure(3)
    subplot(2,1,1)
    hold on
    
    Fz_mean_steady_nondamaged = -mean(Fz_intactwing_steady(:,1)+Fz_intactwing_steady(:,end));
    Fz_mean_steady = -mean(Fz_intactwing_steady(:,1)+Fz_damagedwing_steady(:,end));
    Fz_mean_all = -mean(Fz_intactwing_all(:,1)+Fz_damagedwing_all(:,end));
    plot([0 1],[Fz_mean_steady_nondamaged Fz_mean_steady_nondamaged],'--k','linewidth',2)
    plot([0 1],[Fz_mean_steady Fz_mean_steady],'--','color',[.5 .5 .5],'linewidth',2)
    plot([0 1],[Fz_mean_all Fz_mean_all],'--r','linewidth',2)
    
    plot(t_norm,-Fz_intactwing_steady(:,1),'-k','linewidth',2)
    plot(t_norm,-Fz_damagedwing_steady(:,end),'-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,-Fz_intactwing_all(:,end),'-','color',[1 .5 0],'linewidth',2)
    plot(t_norm,-Fz_damagedwing_all(:,end),'-r','linewidth',2)

%     legend('non-damaged wing steady WB','damaged wing steady WB','non-damaged wing modified WB','damaged wing modified WB','','location','SW')
    xlabel('wingbeat cycle')
    ylabel('normalized vertical force Fz/mg')
    axis([0 1 -1 2])
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
    
    figure(4)
    subplot(3,1,1)
    hold on
    plot(t_norm,stroke_damaged_all(:,1),'-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,stroke_intact_all(:,end),'-','color',[1 .5 0],'linewidth',2)
    plot(t_norm,stroke_damaged_all(:,end),'-r','linewidth',2)
    
    xlabel('wingbeat cycle')
    ylabel('stroke angle')
    axis([0 1 -90 90])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-90:90:90)    
    
    subplot(3,1,2)
    hold on
    plot(t_norm,dev_damaged_all(:,1),'-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,dev_intact_all(:,end),'-','color',[1 .5 0],'linewidth',2)
    plot(t_norm,dev_damaged_all(:,end),'-r','linewidth',2)
    
    xlabel('wingbeat cycle')
    ylabel('deviation angle')
    axis([0 1 -90 90])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-90:90:90)    
    
    subplot(3,1,3)
    hold on
    plot(t_norm,rot_damaged_all(:,1)-90,'-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,rot_intact_all(:,end)-90,'-','color',[1 .5 0],'linewidth',2)
    plot(t_norm,rot_damaged_all(:,end)-90,'-r','linewidth',2)
    
    xlabel('wingbeat cycle')
    ylabel('rotation angle')
    axis([0 1 -90 90])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-90:90:90)    
    
    figure(6)
    subplot(3,1,1)
    hold on
    plot(t_norm,stroke_damaged_all(:,1),'-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,stroke_intact_all(:,end),'-','color',[1 .5 0],'linewidth',2)
    plot(t_norm,stroke_damaged_all(:,end),'-r','linewidth',2)
    
    xlabel('wingbeat cycle')
    ylabel('stroke angle')
    axis([0 1 -90 90])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-90:90:90)    
    
    subplot(3,1,2)
    hold on
    plot(t_norm,dev_damaged_all(:,1),'-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,dev_intact_all(:,end),'-','color',[1 .5 0],'linewidth',2)
    plot(t_norm,dev_damaged_all(:,end),'-r','linewidth',2)
    
    xlabel('wingbeat cycle')
    ylabel('deviation angle')
    axis([0 1 -90 90])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-90:90:90)    
    
    subplot(3,1,3)
    hold on
    plot(t_norm,rot_damaged_all(:,1)-90,'-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,rot_intact_all(:,end)-90,'-','color',[1 .5 0],'linewidth',2)
    plot(t_norm,rot_damaged_all(:,end)-90,'-r','linewidth',2)
    
    xlabel('wingbeat cycle')
    ylabel('rotation angle')
    axis([0 1 -90 90])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-90:90:90)    
    
    %% F&T wingbeat MODs
    figure(7)
    subplot(2,1,1)
    hold on
    
    Fz_steady = (Fz_intactwing_steady(:,end)+Fz_damagedwing_steady(:,end));
    Fz_all = (Fz_intactwing_all(:,end)+Fz_damagedwing_all(:,end));
    Fz_freqMOD = (Fz_intactwing_freqMOD(:,end)+Fz_damagedwing_freqMOD(:,end));
    Fz_strokeMOD = (Fz_intactwing_strokeMOD(:,end)+Fz_damagedwing_strokeMOD(:,end));
    Fz_devMOD = (Fz_intactwing_devMOD(:,end)+Fz_damagedwing_devMOD(:,end));
    Fz_rotMOD = (Fz_intactwing_rotMOD(:,end)+Fz_damagedwing_rotMOD(:,end));

    plot(t_norm,-Fz_all,'-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,-Fz_steady,'-k','linewidth',2)
    plot(t_norm,-Fz_freqMOD,'-b','linewidth',2)
    plot(t_norm,-Fz_strokeMOD,'-g','linewidth',2)
    plot(t_norm,-Fz_devMOD,'-r','linewidth',2)
    plot(t_norm,-Fz_rotMOD,'-c','linewidth',2)

    xlabel('wingbeat cycle')
    ylabel('normalized vertical force Fz/mg')
    axis([0 1 -1 3])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-4:4)
    
    subplot(2,1,2)
    hold on
    
    Mx_steady = (Mx_intactwing_steady(:,end)+Mx_damagedwing_steady(:,end));
    Mx_all = (Mx_intactwing_all(:,end)+Mx_damagedwing_all(:,end));
    Mx_freqMOD = (Mx_intactwing_freqMOD(:,end)+Mx_damagedwing_freqMOD(:,end));
    Mx_strokeMOD = (Mx_intactwing_strokeMOD(:,end)+Mx_damagedwing_strokeMOD(:,end));
    Mx_devMOD = (Mx_intactwing_devMOD(:,end)+Mx_damagedwing_devMOD(:,end));
    Mx_rotMOD = (Mx_intactwing_rotMOD(:,end)+Mx_damagedwing_rotMOD(:,end));

    plot(t_norm,Mx_all,'-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,Mx_steady,'-k','linewidth',2)
    plot(t_norm,Mx_freqMOD,'-b','linewidth',2)
    plot(t_norm,Mx_strokeMOD,'-g','linewidth',2)
    plot(t_norm,Mx_devMOD,'-r','linewidth',2)
    plot(t_norm,Mx_rotMOD,'-c','linewidth',2)
    
    
    legend('all','none','freq','stroke','dev','rot','location','NE')
    xlabel('wingbeat cycle')
    ylabel('normalized roll torque Tx/mgl')
    axis([0 1 -1 1])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-1:1)    
    
%% TE clip    
load(['allMODs_TEClip_freqAsym',num2str(freq_asymFitNr),'_roton1.mat'])

    figure(2)
    subplot(2,1,1)
    hold on
    
    Fz_mean_steady_nondamaged = -mean(Fz_intactwing_steady(:,1)+Fz_intactwing_steady(:,end));
    Fz_mean_steady = -mean(Fz_intactwing_steady(:,1)+Fz_damagedwing_steady(:,end));
    Fz_mean_all = -mean(Fz_intactwing_all(:,1)+Fz_damagedwing_all(:,end));
    plot([0 1],[Fz_mean_steady_nondamaged Fz_mean_steady_nondamaged],'--k','linewidth',2)
    plot([0 1],[Fz_mean_steady Fz_mean_steady],'--','color',[.5 .5 .5],'linewidth',2)
    plot([0 1],[Fz_mean_all Fz_mean_all],'--r','linewidth',2)
    
    plot(t_norm,-Fz_intactwing_steady(:,1),'-k','linewidth',2)
    plot(t_norm,-Fz_damagedwing_steady(:,end),'-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,-Fz_intactwing_all(:,end),'-','color',[1 .5 0],'linewidth',2)
    plot(t_norm,-Fz_damagedwing_all(:,end),'-r','linewidth',2)

%     legend('non-damaged wing steady WB','damaged wing steady WB','non-damaged wing modified WB','damaged wing modified WB','','location','SW')
    xlabel('wingbeat cycle')
    ylabel('normalized vertical force Fz/mg')
    axis([0 1 -1 2])
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
    
    figure(3)
    subplot(2,1,1)
    hold on
    
    Fz_mean_steady_nondamaged = -mean(Fz_intactwing_steady(:,1)+Fz_intactwing_steady(:,end));
    Fz_mean_steady = -mean(Fz_intactwing_steady(:,1)+Fz_damagedwing_steady(:,end));
    Fz_mean_all = -mean(Fz_intactwing_all(:,1)+Fz_damagedwing_all(:,end));
%     plot([0 1],[Fz_mean_steady_nondamaged Fz_mean_steady_nondamaged],'--k','linewidth',2)
%     plot([0 1],[Fz_mean_steady Fz_mean_steady],'--','color',[.5 .5 .5],'linewidth',2)
    plot([0 1],[Fz_mean_all Fz_mean_all],'--b','linewidth',2)
    
%     plot(t_norm,-Fz_intactwing_steady(:,1),'-k','linewidth',2)
%     plot(t_norm,-Fz_damagedwing_steady(:,end),'-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,-Fz_intactwing_all(:,end),'-c','linewidth',2)
    plot(t_norm,-Fz_damagedwing_all(:,end),'-b','linewidth',2)

%     legend('non-damaged wing steady WB','damaged wing steady WB','non-damaged wing modified WB','damaged wing modified WB','','location','SW')
    xlabel('wingbeat cycle')
    ylabel('normalized vertical force Fz/mg')
    axis([0 1 -1 2])
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
    

    figure(5)
    subplot(3,1,1)
    hold on
    plot(t_norm,stroke_damaged_all(:,1),'-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,stroke_intact_all(:,end),'-','color',[1 .5 0],'linewidth',2)
    plot(t_norm,stroke_damaged_all(:,end),'-r','linewidth',2)
    
    xlabel('wingbeat cycle')
    ylabel('stroke angle')
    axis([0 1 -90 90])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-90:90:90)    
    
    subplot(3,1,2)
    hold on
    plot(t_norm,dev_damaged_all(:,1),'-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,dev_intact_all(:,end),'-','color',[1 .5 0],'linewidth',2)
    plot(t_norm,dev_damaged_all(:,end),'-r','linewidth',2)
    
    xlabel('wingbeat cycle')
    ylabel('deviation angle')
    axis([0 1 -90 90])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-90:90:90)    
    
    subplot(3,1,3)
    hold on
    plot(t_norm,rot_damaged_all(:,1)-90,'-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,rot_intact_all(:,end)-90,'-','color',[1 .5 0],'linewidth',2)
    plot(t_norm,rot_damaged_all(:,end)-90,'-r','linewidth',2)
    
    xlabel('wingbeat cycle')
    ylabel('rotation angle')
    axis([0 1 -90 90])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-90:90:90)    

    figure(6)
    subplot(3,1,1)
    hold on
%     plot(t_norm,stroke_damaged_all(:,1),'-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,stroke_intact_all(:,end),'-c','linewidth',2)
    plot(t_norm,stroke_damaged_all(:,end),'-b','linewidth',2)
    
    xlabel('wingbeat cycle')
    ylabel('stroke angle')
    axis([0 1 -90 90])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-90:90:90)    
    
    subplot(3,1,2)
    hold on
%     plot(t_norm,dev_damaged_all(:,1),'-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,dev_intact_all(:,end),'-c','linewidth',2)
    plot(t_norm,dev_damaged_all(:,end),'-b','linewidth',2)
    
    xlabel('wingbeat cycle')
    ylabel('deviation angle')
    axis([0 1 -90 90])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-90:90:90)    
    
    subplot(3,1,3)
    hold on
    plot(t_norm,rot_damaged_all(:,1)-90,'-','color',[.5 .5 .5],'linewidth',2)
    plot(t_norm,rot_intact_all(:,end)-90,'-c','linewidth',2)
    plot(t_norm,rot_damaged_all(:,end)-90,'-b','linewidth',2)
    
    xlabel('wingbeat cycle')
    ylabel('rotation angle')
    axis([0 1 -90 90])
    set(gca,'xtick',0:1)
    set(gca,'ytick',-90:90:90)    

%% save figs
mkdir('qsModel_FnM_TempDynamics_TipTEcut')
cd('qsModel_FnM_TempDynamics_TipTEcut')

figure(1)
    saveas(gca,['FzFlip_Mx_TempDynamics_TipClip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['FzFlip_Mx_TempDynamics_TipClip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['FzFlip_Mx_TempDynamics_TipClip_asympFit',num2str(freq_asymFitNr),'.svg'])

figure(2)
    saveas(gca,['FzFlip_Mx_TempDynamics_TEclip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['FzFlip_Mx_TempDynamics_TEclip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['FzFlip_Mx_TempDynamics_TEclip_asympFit',num2str(freq_asymFitNr),'.svg'])
    
figure(3)
    saveas(gca,['FzFlip_Mx_TempDynamics_TipTEclip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['FzFlip_Mx_TempDynamics_TipTEclip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['FzFlip_Mx_TempDynamics_TipTEclip_asympFit',num2str(freq_asymFitNr),'.svg'])

figure(4)
    saveas(gca,['WBkinAngles_TempDynamics_TipClip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['WBkinAngles_TempDynamics_TipClip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['WBkinAngles_TempDynamics_TipClip_asympFit',num2str(freq_asymFitNr),'.svg'])

figure(5)
    saveas(gca,['WBkinAngles_TempDynamics_TEclip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['WBkinAngles_TempDynamics_TEclip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['WBkinAngles_TempDynamics_TEclip_asympFit',num2str(freq_asymFitNr),'.svg'])

figure(6)
    saveas(gca,['WBkinAngles_TempDynamics_TipTEclip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['WBkinAngles_TempDynamics_TipTEclip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['WBkinAngles_TempDynamics_TipTEclip_asympFit',num2str(freq_asymFitNr),'.svg'])

figure(7)
    saveas(gca,['FzFlip_Mx_wbMODs_TempDynamics_TipClip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['FzFlip_Mx_wbMODs_TempDynamics_TipClip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['FzFlip_Mx_wbMODs_TempDynamics_TipClip_asympFit',num2str(freq_asymFitNr),'.svg'])
    
    cd ..
    