clear;
clc;
close all
warning off

addpath('/home/florian/Dropbox/WORK/flytracker/flytracker');
addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/mex/');
addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/core/');
addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/results/');

loadname=dir('WBdataset_temporal_dynamics_*')
loadname = loadname.name;
load(loadname)

load('FnMqs_timeseries.mat')

mkdir('MSfigs_timeseries_KinNqsTorque')
cd('MSfigs_timeseries_KinNqsTorque')

%% freq
figure
subplot(3,1,1)
hold on
ciplot(f_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)-1.96*f_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),f_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)+1.96*f_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),t_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0))
plot(t_wb_seq_mean_all,f_wb_seq_mean_all,'-k.')

% xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('freq [1/sec]')
ylim([180 200])
set(gca,'XTick',[-.05:.05:.05],'XTickLabel',[])
set(gca,'YTick',[0:10:200])

%% Force
subplot(3,1,2)
hold on

% body dyn torque
ciplot(F_mean_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)-1.96*F_mean_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),F_mean_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)+1.96*F_mean_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),t_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0))
plot(t_wb_seq_mean_all,F_mean_wb_seq_mean_all,'-k.')

% QSmodel torque
F_trans_norm_timeseries_WBmean = sqrt(Fx_trans_norm_timeseries_WBmean.^2 + Fy_trans_norm_timeseries_WBmean.^2 + Fz_trans_norm_timeseries_WBmean.^2);
F_transNrot_norm_timeseries_WBmean = sqrt(Fx_transNrot_norm_timeseries_WBmean.^2 + Fy_transNrot_norm_timeseries_WBmean.^2 + Fz_transNrot_norm_timeseries_WBmean.^2);

plot(t_wb_seq_mean_all,F_trans_norm_timeseries_WBmean,'-m.')
plot(t_wb_seq_mean_all,F_transNrot_norm_timeseries_WBmean,'-g.')

% xlabel('time [sec]')
ylabel('F/mg')
xlim([-.05 .06]) 
ylim([.75 1.25])
set(gca,'XTick',[-.05:.05:.05])
set(gca,'YTick',[0 :.25:2])

saveas(gca,['MSfig_timeseries_freqNforce.fig'])
saveas(gca,['MSfig_timeseries_freqNforce.png'])
plot2svg(['MSfig_timeseries_freqNforce.svg'])

%% roll torque
figure
subplot(3,1,1)
hold on

% body dyn torque
ciplot(Mroll_mean_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)-1.96*Mroll_mean_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),Mroll_mean_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)+1.96*Mroll_mean_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),t_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0))
plot(t_wb_seq_mean_all,Mroll_mean_wb_seq_mean_all,'-k.')

% QSmodel torque
plot(t_wb_seq_mean_all,Mx_trans_norm_timeseries_WBmean,'-m.')
plot(t_wb_seq_mean_all,Mx_transNrot_norm_timeseries_WBmean,'-g.')

% xlabel('time [sec]')
ylabel('Troll/M(g+a)')
xlim([-.05 .06]) 
ylim([-.03 .03])
set(gca,'XTick',[-.05:.05:.05],'XTickLabel',[])
set(gca,'YTick',[-.03 :.03:1])

%% pitch torque
subplot(3,1,2)
hold on

% body dyn torque
ciplot(Mpitch_mean_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)-1.96*Mpitch_mean_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),Mpitch_mean_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)+1.96*Mpitch_mean_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),t_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0))
plot(t_wb_seq_mean_all,Mpitch_mean_wb_seq_mean_all,'-k.')

% QSmodel torque
plot(t_wb_seq_mean_all,My_trans_norm_timeseries_WBmean,'-m.')
plot(t_wb_seq_mean_all,My_transNrot_norm_timeseries_WBmean,'-g.')

% xlabel('time [sec]')
ylabel('Tpitch/M(g+a)')
xlim([-.05 .06]) 
ylim([-.03 .03])
set(gca,'XTick',[-.05:.05:.05],'XTickLabel',[])
set(gca,'YTick',[-.03 :.03:1])

%% yaw torque
subplot(3,1,3)
hold on

% body dyn torque
ciplot(Myaw_mean_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)-1.96*Myaw_mean_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),Myaw_mean_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)+1.96*Myaw_mean_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),t_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0))
plot(t_wb_seq_mean_all,Myaw_mean_wb_seq_mean_all,'-k.')

% QSmodel torque
plot(t_wb_seq_mean_all,Mz_trans_norm_timeseries_WBmean,'-m.')
plot(t_wb_seq_mean_all,Mz_transNrot_norm_timeseries_WBmean,'-g.')

xlabel('time [sec]')
ylabel('Tyaw/M(g+a)')
xlim([-.05 .06]) 
ylim([0 .06])
set(gca,'XTick',[-.05:.05:.05])
set(gca,'YTick',[-.03 :.03:1])

saveas(gca,['MSfig_timeseries_TrollNpitchNyaw.fig'])
saveas(gca,['MSfig_timeseries_TrollNpitchNyaw.png'])
plot2svg(['MSfig_timeseries_TrollNpitchNyaw.svg'])

%% Raxis torque
figure
subplot(3,1,1)
hold on

% body dyn torque
ciplot(M_R_mean_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)-1.96*M_R_mean_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),M_R_mean_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)+1.96*M_R_mean_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),t_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0))
plot(t_wb_seq_mean_all,M_R_mean_wb_seq_mean_all,'-k.')

% QSmodel torque
plot(t_wb_seq_mean_all,M_R_trans_norm_timeseries_WBmean,'-m.')
plot(t_wb_seq_mean_all,M_R_transNrot_norm_timeseries_WBmean,'-g.')

% xlabel('time [sec]')
ylabel('T_R/M(g+a)')
xlim([-.05 .06]) 
ylim([-.03 .03])
set(gca,'XTick',[-.05:.05:.05],'XTickLabel',[])
set(gca,'YTick',[-.03 :.03:1])

%% Laxis torque
subplot(3,1,2)
hold on

% body dyn torque
ciplot(M_L_mean_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)-1.96*M_L_mean_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),M_L_mean_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)+1.96*M_L_mean_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),t_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0))
plot(t_wb_seq_mean_all,M_L_mean_wb_seq_mean_all,'-k.')

% QSmodel torque
plot(t_wb_seq_mean_all,M_L_trans_norm_timeseries_WBmean,'-m.')
plot(t_wb_seq_mean_all,M_L_transNrot_norm_timeseries_WBmean,'-g.')

% xlabel('time [sec]')
ylabel('T_L/M(g+a)')
xlim([-.05 .06]) 
ylim([-.03 .03])
set(gca,'XTick',[-.05:.05:.05],'XTickLabel',[])
set(gca,'YTick',[-.03 :.03:1])

%% yaw torque
subplot(3,1,3)
hold on

% body dyn torque
ciplot(Myaw_mean_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)-1.96*Myaw_mean_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),Myaw_mean_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)+1.96*Myaw_mean_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),t_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0))
plot(t_wb_seq_mean_all,Myaw_mean_wb_seq_mean_all,'-k.')

% QSmodel torque
plot(t_wb_seq_mean_all,Mz_trans_norm_timeseries_WBmean,'-m.')
plot(t_wb_seq_mean_all,Mz_transNrot_norm_timeseries_WBmean,'-g.')

xlabel('time [sec]')
ylabel('Tyaw/M(g+a)')
xlim([-.05 .06]) 
ylim([0 .06])
set(gca,'XTick',[-.05:.05:.05])
set(gca,'YTick',[-.03 :.03:1])

saveas(gca,['MSfig_timeseries_Torque_RnLaxis_Tyaw.fig'])
saveas(gca,['MSfig_timeseries_Torque_RnLaxis_Tyaw.png'])
plot2svg(['MSfig_timeseries_Torque_RnLaxis_Tyaw.svg'])

%% wingbeat kin
figure

% stroke
subplot(3,1,1)
hold on

plot(t_wb_seq_bins_mean_all,stroke_wb_R_seq_bins_mean_all,'-b','linewidth',.5)
plot(t_wb_seq_bins_mean_all,stroke_wb_L_seq_bins_mean_all,'-r','linewidth',.5)

% xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('stroke angle [deg]')
ylim([-90 90])
set(gca,'XTick',[-.05:.05:.05],'XTickLabel',[])
set(gca,'YTick',[-90:45:90])

% deviation
subplot(3,1,2)
hold on

plot(t_wb_seq_bins_mean_all,dev_wb_R_seq_bins_mean_all,'-b','linewidth',.5)
plot(t_wb_seq_bins_mean_all,dev_wb_L_seq_bins_mean_all,'-r','linewidth',.5)

% xlabel('time [sec]')
xlim([-.05 .06])
ylabel('deviation angle [deg]')
ylim([-30 30])
set(gca,'XTick',[-.05:.05:.05],'XTickLabel',[])
set(gca,'YTick',[-90:30:90])

% rotation
subplot(3,1,3)
hold on
plot(t_wb_seq_bins_mean_all,pitch_wb_R_seq_bins_mean_all,'-b','linewidth',.5)
plot(t_wb_seq_bins_mean_all,pitch_wb_L_seq_bins_mean_all,'-r','linewidth',.5)

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('rotation angle [deg]')
ylim([0 180])
set(gca,'XTick',[-.05:.05:.05])
set(gca,'YTick',[-90:45:180])

saveas(gca,['MSfig_timeseries_WBangles.fig'])
saveas(gca,['MSfig_timeseries_WBangles.png'])
plot2svg(['MSfig_timeseries_WBangles.svg'])

%% WB parameters

%% stroke
figure

% L+R
subplot(3,1,1)
hold on
plot(t_wb_seq_mean_all,Astroke_wb_R_seq_bins,'-b.')
plot(t_wb_seq_mean_all,Astroke_wb_L_seq_bins,'-r.')
% plot(t_wb_seq_mean_all,Astroke_wb_seq_bins,'-g.'])

% xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('stroke angle amplitude [deg]')
ylim([125 135])
set(gca,'XTick',[-.05:.05:.05],'XTickLabel',[])
set(gca,'YTick',[-90:5:190])

% L-R
subplot(3,1,2)
hold on
plot(t_wb_seq_mean_all,-dAstroke_wb_seq_bins,'-g.','color',[.5 .5 .5])
plot(t_wb_seq_mean_all,Dstroke_wb_seq_bins_MAX,'-r.','color',[1 .5 0])
plot(t_wb_seq_mean_all,Dstroke_wb_seq_bins_MIN,'-b.','color',[0 .5 1])

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('stroke angle L-R [deg]')
ylim([-5 5])
set(gca,'XTick',[-.05:.05:.05])
set(gca,'YTick',[-90:5:90])

saveas(gca,['MSfig_timeseries_strokeParams.fig'])
saveas(gca,['MSfig_timeseries_strokeParams.png'])
plot2svg(['MSfig_timeseries_strokeParams.svg'])

%% pitch
figure

% Downstroke
subplot(3,1,1)
hold on
plot(t_wb_seq_mean_all,pitch_wb_R_seq_bins_MAXmidDS,'-b.')
plot(t_wb_seq_mean_all,pitch_wb_L_seq_bins_MAXmidDS,'-r.')
% plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MAXmidDS,'-g.')

% xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('rotation angle downstroke [deg]')
ylim([50 60])
set(gca,'XTick',[-.05:.05:.05],'XTickLabel',[])
set(gca,'YTick',[-90:5:90])

% Upstroke
subplot(3,1,2)
hold on
plot(t_wb_seq_mean_all,pitch_wb_R_seq_bins_MINmidUS,'-b.')
plot(t_wb_seq_mean_all,pitch_wb_L_seq_bins_MINmidUS,'-r.')
% plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MINmidUS,'-g.')

% xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('rotation angle upstroke [deg]')
ylim([-50 -40])
set(gca,'XTick',[-.05:.05:.05],'XTickLabel',[])
set(gca,'YTick',[-90:5:90])

% L-R
subplot(3,1,3)
hold on
plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MAXmidDS,'-r.','color',[1 .5 0])
plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MINmidUS,'-b.','color',[0 .5 1])
% plot(t_wb_seq_mean_all,-dApitch_wb_seq_bins,'-g.','color',[.5 .5 .5])

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('rotation angle L-R [deg]')
ylim([-5 5])
set(gca,'XTick',[-.05:.05:.05])
set(gca,'YTick',[-90:5:90])

saveas(gca,['MSfig_timeseries_rotationParams.fig'])
saveas(gca,['MSfig_timeseries_rotationParams.png'])
plot2svg(['MSfig_timeseries_rotationParams.svg'])

%% dev
figure

% Downstroke
subplot(3,1,1)
hold on
plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_DS,'-b.')
plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_DS,'-r.')

% xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('deviation amplitude downstroke [deg]')
ylim([10 20])
set(gca,'XTick',[-.05:.05:.05],'XTickLabel',[])
set(gca,'YTick',[-90:5:90])

% Upstroke
subplot(3,1,2)
hold on

plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_US,'-b.')
plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_US,'-r.')

% xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('deviation amplitude upstroke [deg]')
ylim([10 20])
set(gca,'XTick',[-.05:.05:.05],'XTickLabel',[])
set(gca,'YTick',[-90:5:90])

% L-R
subplot(3,1,3)
hold on
plot(t_wb_seq_mean_all,dAdev_wb_seq_bins_DS,'-b.','color',[0 .5 1])
plot(t_wb_seq_mean_all,dAdev_wb_seq_bins_US,'-r.','color',[1 .5 0])
% plot(t_wb_seq_mean_all,-dAdev_wb_seq_bins,'-g.','color',[.5 .5 .5])

xlabel('time [sec]')
xlim([-.05 .06]) 
ylabel('deviation amplitude L-R [deg]')
ylim([-5 5])
set(gca,'XTick',[-.05:.05:.05])
set(gca,'YTick',[-90:5:90])

saveas(gca,['MSfig_timeseries_deviationParams.fig'])
saveas(gca,['MSfig_timeseries_deviationParams.png'])
plot2svg(['MSfig_timeseries_deviationParams.svg'])



cd ..


















