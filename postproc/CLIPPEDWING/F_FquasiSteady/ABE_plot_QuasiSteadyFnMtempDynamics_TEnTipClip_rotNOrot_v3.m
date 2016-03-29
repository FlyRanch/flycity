clear;
clc;
close all
warning off

%% TEclip ROT on
clear
freq_asymFitNr = 10;
plot_nr = 1;
load(['allMODs_TEclip_freqAsym',num2str(freq_asymFitNr),'_roton1.mat'])

spanloc = .7;
Fscale = 4;
Fref = 2;
Vscale = 2;
Vref = 4;
dref = 10;
c_intact = 10;
c_damaged = 5;
Nsteps = 23;

% Fz rot & NO rot
figure(plot_nr)
plot_qsFtempDyn_rotON

% wb kin tempdyn
figure(10+plot_nr)
plot_wbKin_damagedIntactNondamaged

% wb kin & forces lollypops
figure(100+plot_nr)
plot_wingprofNforce_lollypops_oneplot

% wb kin & forces lollypops separate
Fscale = 8;
Fref = 2;
c_intact = 15;
c_damaged = 7.5;
Nsteps = 10;
plot_wingprofNforce_lollypops

% Fz & Mx
figure(1000+plot_nr)
plot_qsFnMtempDyn_rotON

%% TEclip ROT off
clear
freq_asymFitNr = 10;
plot_nr = 1;
load(['allMODs_TEclip_freqAsym',num2str(freq_asymFitNr),'_roton0.mat'])

% Fz rot & NO rot
figure(plot_nr)
plot_qsFtempDyn_rotOFF


%% TipClip ROT on
clear
freq_asymFitNr = 10;
plot_nr = 2;
load(['allMODs_TipClip_freqAsym',num2str(freq_asymFitNr),'_roton1.mat'])

spanloc = .7;
Fscale = 4;
Fref = 2;
Vscale = 2;
Vref = 4;
dref = 10;
c_intact = 10;
c_damaged = 10;
Nsteps = 23;

% Fz rot & NO rot
figure(plot_nr)
plot_qsFtempDyn_rotON

% wb kin tempdyn
figure(10+plot_nr)
plot_wbKin_damagedIntactNondamaged

% wb kin & forces lollypops
figure(100+plot_nr)
plot_wingprofNforce_lollypops_oneplot

% wb kin & forces lollypops separate
Fscale = 8;
Fref = 2;
c_intact = 15;
c_damaged = 15;
Nsteps = 10;
plot_wingprofNforce_lollypops

% Fz & Mx
figure(1000+plot_nr)
plot_qsFnMtempDyn_rotON

%% TEclip ROT off
clear
freq_asymFitNr = 10;
plot_nr = 2;
load(['allMODs_TEclip_freqAsym',num2str(freq_asymFitNr),'_roton0.mat'])

% Fz rot & NO rot
figure(plot_nr)
plot_qsFtempDyn_rotOFF

%% save figs
mkdir('qsModel_FnM_TempDynamics')
cd('qsModel_FnM_TempDynamics')

figure(1)
saveas(gca,['FzFlip_TempDynamics_TEclip_rotNOrot_asympFit',num2str(freq_asymFitNr),'.fig'])
saveas(gca,['FzFlip_TempDynamics_TEclip_rotNOrot_asympFit',num2str(freq_asymFitNr),'.png'])
plot2svg(  ['FzFlip_TempDynamics_TEclip_rotNOrot_asympFit',num2str(freq_asymFitNr),'.svg'])

figure(2)
saveas(gca,['FzFlip_TempDynamics_TipClip_rotNOrot_asympFit',num2str(freq_asymFitNr),'.fig'])
saveas(gca,['FzFlip_TempDynamics_TipClip_rotNOrot_asympFit',num2str(freq_asymFitNr),'.png'])
plot2svg(  ['FzFlip_TempDynamics_TipClip_rotNOrot_asympFit',num2str(freq_asymFitNr),'.svg'])

figure(11)
saveas(gca,['WBkin_TempDynamics_TEclip_asympFit',num2str(freq_asymFitNr),'.fig'])
saveas(gca,['WBkin_TempDynamics_TEclip_asympFit',num2str(freq_asymFitNr),'.png'])
plot2svg(  ['WBkin_TempDynamics_TEclip_asympFit',num2str(freq_asymFitNr),'.svg'])

figure(12)
saveas(gca,['WBkin_TempDynamics_TipClip_asympFit',num2str(freq_asymFitNr),'.fig'])
saveas(gca,['WBkin_TempDynamics_TipClip_asympFit',num2str(freq_asymFitNr),'.png'])
plot2svg(  ['WBkin_TempDynamics_TipClip_asympFit',num2str(freq_asymFitNr),'.svg'])

figure(101)
saveas(gca,['WBkinNforce_TempDynamics_TEclip_asympFit',num2str(freq_asymFitNr),'.fig'])
saveas(gca,['WBkinNforce_TempDynamics_TEclip_asympFit',num2str(freq_asymFitNr),'.png'])
plot2svg(  ['WBkinNforce_TempDynamics_TEclip_asympFit',num2str(freq_asymFitNr),'.svg'])

figure(102)
saveas(gca,['WBkinNforce_TempDynamics_TipClip_asympFit',num2str(freq_asymFitNr),'.fig'])
saveas(gca,['WBkinNforce_TempDynamics_TipClip_asympFit',num2str(freq_asymFitNr),'.png'])
plot2svg(  ['WBkinNforce_TempDynamics_TipClip_asympFit',num2str(freq_asymFitNr),'.svg'])

figure(103)
saveas(gca,['WBkinNforce_TempDynamics_TEclip_damaged_asympFit',num2str(freq_asymFitNr),'.fig'])
saveas(gca,['WBkinNforce_TempDynamics_TEclip_damaged_asympFit',num2str(freq_asymFitNr),'.png'])
plot2svg(  ['WBkinNforce_TempDynamics_TEclip_damaged_asympFit',num2str(freq_asymFitNr),'.svg'])

figure(104)
saveas(gca,['WBkinNforce_TempDynamics_TipClip_damaged_asympFit',num2str(freq_asymFitNr),'.fig'])
saveas(gca,['WBkinNforce_TempDynamics_TipClip_damaged_asympFit',num2str(freq_asymFitNr),'.png'])
plot2svg(  ['WBkinNforce_TempDynamics_TipClip_damaged_asympFit',num2str(freq_asymFitNr),'.svg'])

figure(105)
saveas(gca,['WBkinNforce_TempDynamics_TEclip_intact_asympFit',num2str(freq_asymFitNr),'.fig'])
saveas(gca,['WBkinNforce_TempDynamics_TEclip_intact_asympFit',num2str(freq_asymFitNr),'.png'])
plot2svg(  ['WBkinNforce_TempDynamics_TEclip_intact_asympFit',num2str(freq_asymFitNr),'.svg'])

figure(106)
saveas(gca,['WBkinNforce_TempDynamics_TipClip_intact_asympFit',num2str(freq_asymFitNr),'.fig'])
saveas(gca,['WBkinNforce_TempDynamics_TipClip_intact_asympFit',num2str(freq_asymFitNr),'.png'])
plot2svg(  ['WBkinNforce_TempDynamics_TipClip_intact_asympFit',num2str(freq_asymFitNr),'.svg'])

figure(1001)
saveas(gca,['FzMx_TempDynamics_TEclip_asympFit',num2str(freq_asymFitNr),'.fig'])
saveas(gca,['FzMx_TempDynamics_TEclip_asympFit',num2str(freq_asymFitNr),'.png'])
plot2svg(  ['FzMx_TempDynamics_TEclip_asympFit',num2str(freq_asymFitNr),'.svg'])

figure(1002)
saveas(gca,['FzMx_TempDynamics_TipClip_asympFit',num2str(freq_asymFitNr),'.fig'])
saveas(gca,['FzMx_TempDynamics_TipClip_asympFit',num2str(freq_asymFitNr),'.png'])
plot2svg(  ['FzMx_TempDynamics_TipClip_asympFit',num2str(freq_asymFitNr),'.svg'])

cd ..
    