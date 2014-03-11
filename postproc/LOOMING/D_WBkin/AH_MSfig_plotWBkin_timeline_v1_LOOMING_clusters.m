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
% 
% load('FnMqs_timeseries.mat')

mkdir('MSfigs_WBnBodyKin_TempDynamics')
cd('MSfigs_WBnBodyKin_TempDynamics')

%% freq
figure
subplot(3,1,1)
hold on
ciplot(f_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)-1.96*f_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),...
       f_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)+1.96*f_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),...
       t_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0))
plot(t_wb_seq_mean_all,f_wb_seq_mean_all,'-k.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('freq [1/sec]')
ylim([180 220])
set(gca,'XTick',[-.04:.02:.04])
set(gca,'YTick',[0:20:220])

subplot(3,1,2)
hold on
ciplot(f_wb_seq_mean_all_c1(isnan(f_wb_seq_mean_all_c1)==0)-1.96*f_wb_seq_ste_all_c1(isnan(f_wb_seq_mean_all_c1)==0),...
       f_wb_seq_mean_all_c1(isnan(f_wb_seq_mean_all_c1)==0)+1.96*f_wb_seq_ste_all_c1(isnan(f_wb_seq_mean_all_c1)==0),...
       t_wb_seq_mean_all(isnan(f_wb_seq_mean_all_c1)==0))
plot(t_wb_seq_mean_all,f_wb_seq_mean_all_c1,'-k.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('freq c1 [1/sec]')
ylim([180 220])
set(gca,'XTick',[-.04:.02:.04])
set(gca,'YTick',[0:20:220])

subplot(3,1,3)
hold on
ciplot(f_wb_seq_mean_all_c2(isnan(f_wb_seq_mean_all_c2)==0)-1.96*f_wb_seq_ste_all_c2(isnan(f_wb_seq_mean_all_c2)==0),...
       f_wb_seq_mean_all_c2(isnan(f_wb_seq_mean_all_c2)==0)+1.96*f_wb_seq_ste_all_c2(isnan(f_wb_seq_mean_all_c2)==0),...
       t_wb_seq_mean_all(isnan(f_wb_seq_mean_all_c2)==0))
plot(t_wb_seq_mean_all,f_wb_seq_mean_all_c2,'-k.')

xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('freq c2 [1/sec]')
ylim([180 220])
set(gca,'XTick',[-.04:.02:.04])
set(gca,'YTick',[0:20:220])

saveas(gca,['MSfig_timeseries_freq.fig'])
saveas(gca,['MSfig_timeseries_freq.png'])
plot2svg(['MSfig_timeseries_freq.svg'])

figure
subplot(3,1,1)
hold on
ciplot(f_wb_seq_mean_all_c3(isnan(f_wb_seq_mean_all_c3)==0)-1.96*f_wb_seq_ste_all_c3(isnan(f_wb_seq_mean_all_c3)==0),...
       f_wb_seq_mean_all_c3(isnan(f_wb_seq_mean_all_c3)==0)+1.96*f_wb_seq_ste_all_c3(isnan(f_wb_seq_mean_all_c3)==0),...
       t_wb_seq_mean_all(isnan(f_wb_seq_mean_all_c3)==0))
plot(t_wb_seq_mean_all,f_wb_seq_mean_all_c3,'-k.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('freq c3 [1/sec]')
ylim([180 220])
set(gca,'XTick',[-.04:.02:.04])
set(gca,'YTick',[0:20:220])

subplot(3,1,2)
hold on
ciplot(f_wb_seq_mean_all_c4(isnan(f_wb_seq_mean_all_c4)==0)-1.96*f_wb_seq_ste_all_c4(isnan(f_wb_seq_mean_all_c4)==0),...
       f_wb_seq_mean_all_c4(isnan(f_wb_seq_mean_all_c4)==0)+1.96*f_wb_seq_ste_all_c4(isnan(f_wb_seq_mean_all_c4)==0),...
       t_wb_seq_mean_all(isnan(f_wb_seq_mean_all_c4)==0))
plot(t_wb_seq_mean_all,f_wb_seq_mean_all_c4,'-k.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('freq c4 [1/sec]')
ylim([180 220])
set(gca,'XTick',[-.04:.02:.04])
set(gca,'YTick',[0:20:220])

saveas(gca,['MSfig_timeseries_freq2.fig'])
saveas(gca,['MSfig_timeseries_freq2.png'])
plot2svg(['MSfig_timeseries_freq2.svg'])

%% wingbeat kin
%% ALL
figure

% stroke
subplot(3,1,1)
hold on
% 
% ciplot(stroke_wb_L_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0)-1.96*stroke_wb_L_seq_bins_ste_all(isnan(t_wb_seq_bins_mean_all)==0),...
%     stroke_wb_L_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0)+1.96*stroke_wb_L_seq_bins_ste_all(isnan(t_wb_seq_bins_mean_all)==0),...
%     t_wb_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0),[0 .5 1])
% ciplot(stroke_wb_R_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0)-1.96*stroke_wb_R_seq_bins_ste_all(isnan(t_wb_seq_bins_mean_all)==0),...
%     stroke_wb_R_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0)+1.96*stroke_wb_R_seq_bins_ste_all(isnan(t_wb_seq_bins_mean_all)==0),...
%     t_wb_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0),[1 0 0])

plot(t_wb_seq_bins_mean_all,stroke_wb_L_seq_bins_mean_all,'-b','linewidth',.5)
plot(t_wb_seq_bins_mean_all,stroke_wb_R_seq_bins_mean_all,'-r','linewidth',.5)

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('stroke angle [deg]')
ylim([-90 90])
set(gca,'YTick',[-90:45:90])

% deviation
subplot(3,1,2)
hold on

% ciplot(dev_wb_L_seq_bins_mean_all-1.96*dev_wb_L_seq_bins_ste_all,...
%     dev_wb_L_seq_bins_mean_all+1.96*dev_wb_L_seq_bins_ste_all,t_wb_seq_bins_mean_all,[0 .5 1])
% ciplot(dev_wb_R_seq_bins_mean_all-1.96*dev_wb_R_seq_bins_ste_all,...
%     dev_wb_R_seq_bins_mean_all+1.96*dev_wb_R_seq_bins_ste_all,t_wb_seq_bins_mean_all,[1 0 0])
% 
plot(t_wb_seq_bins_mean_all,dev_wb_L_seq_bins_mean_all,'-b','linewidth',.5)
plot(t_wb_seq_bins_mean_all,dev_wb_R_seq_bins_mean_all,'-r','linewidth',.5)

% xlabel('time [sec]')
xlim([-.04 .04])
ylabel('deviation angle [deg]')
ylim([-30 15])
set(gca,'YTick',[-90:30:90])

% rotation
subplot(3,1,3)
hold on
% plot(t_wb_seq_bins_mean_all,pitch_wb_L_seq_bins_all,'b','linewidth',.5)
% plot(t_wb_seq_bins_mean_all,pitch_wb_R_seq_bins_all,'r','linewidth',.5)

plot(t_wb_seq_bins_mean_all,pitch_wb_L_seq_bins_mean_all,'-b','linewidth',.5)
plot(t_wb_seq_bins_mean_all,pitch_wb_R_seq_bins_mean_all,'-r','linewidth',.5)
% 
% ciplot(pitch_wb_L_seq_bins_mean_all-1.96*pitch_wb_L_seq_bins_ste_all,...
%     pitch_wb_L_seq_bins_mean_all+1.96*pitch_wb_L_seq_bins_ste_all,t_wb_seq_bins_mean_all,[0 .5 1])
% ciplot(pitch_wb_R_seq_bins_mean_all-1.96*pitch_wb_R_seq_bins_ste_all,...
%     pitch_wb_R_seq_bins_mean_all+1.96*pitch_wb_R_seq_bins_ste_all,t_wb_seq_bins_mean_all,[1 0 0])

xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('rotation angle [deg]')
ylim([0 180])
set(gca,'YTick',[-90:45:180])

saveas(gca,['MSfig_timeseries_WBangles.fig'])
saveas(gca,['MSfig_timeseries_WBangles.png'])
plot2svg(['MSfig_timeseries_WBangles.svg'])

%% c1
figure

% stroke
subplot(3,1,1)
hold on
% 
% ciplot(stroke_wb_L_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0)-1.96*stroke_wb_L_seq_bins_ste_all(isnan(t_wb_seq_bins_mean_all)==0),...
%     stroke_wb_L_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0)+1.96*stroke_wb_L_seq_bins_ste_all(isnan(t_wb_seq_bins_mean_all)==0),...
%     t_wb_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0),[0 .5 1])
% ciplot(stroke_wb_R_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0)-1.96*stroke_wb_R_seq_bins_ste_all(isnan(t_wb_seq_bins_mean_all)==0),...
%     stroke_wb_R_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0)+1.96*stroke_wb_R_seq_bins_ste_all(isnan(t_wb_seq_bins_mean_all)==0),...
%     t_wb_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0),[1 0 0])

plot(t_wb_seq_bins_mean_all,stroke_wb_L_seq_bins_mean_all_c1,'-b','linewidth',.5)
plot(t_wb_seq_bins_mean_all,stroke_wb_R_seq_bins_mean_all_c1,'-r','linewidth',.5)

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('stroke angle c1 [deg]')
ylim([-90 90])
set(gca,'YTick',[-90:45:90])

% deviation
subplot(3,1,2)
hold on

% ciplot(dev_wb_L_seq_bins_mean_all-1.96*dev_wb_L_seq_bins_ste_all,...
%     dev_wb_L_seq_bins_mean_all+1.96*dev_wb_L_seq_bins_ste_all,t_wb_seq_bins_mean_all,[0 .5 1])
% ciplot(dev_wb_R_seq_bins_mean_all-1.96*dev_wb_R_seq_bins_ste_all,...
%     dev_wb_R_seq_bins_mean_all+1.96*dev_wb_R_seq_bins_ste_all,t_wb_seq_bins_mean_all,[1 0 0])
% 
plot(t_wb_seq_bins_mean_all,dev_wb_L_seq_bins_mean_all_c1,'-b','linewidth',.5)
plot(t_wb_seq_bins_mean_all,dev_wb_R_seq_bins_mean_all_c1,'-r','linewidth',.5)

% xlabel('time [sec]')
xlim([-.04 .04])
ylabel('deviation angle c1 [deg]')
ylim([-30 15])
set(gca,'YTick',[-90:30:90])

% rotation
subplot(3,1,3)
hold on
% plot(t_wb_seq_bins_mean_all,pitch_wb_L_seq_bins_all,'b','linewidth',.5)
% plot(t_wb_seq_bins_mean_all,pitch_wb_R_seq_bins_all,'r','linewidth',.5)

plot(t_wb_seq_bins_mean_all,pitch_wb_L_seq_bins_mean_all_c1,'-b','linewidth',.5)
plot(t_wb_seq_bins_mean_all,pitch_wb_R_seq_bins_mean_all_c1,'-r','linewidth',.5)
% 
% ciplot(pitch_wb_L_seq_bins_mean_all-1.96*pitch_wb_L_seq_bins_ste_all,...
%     pitch_wb_L_seq_bins_mean_all+1.96*pitch_wb_L_seq_bins_ste_all,t_wb_seq_bins_mean_all,[0 .5 1])
% ciplot(pitch_wb_R_seq_bins_mean_all-1.96*pitch_wb_R_seq_bins_ste_all,...
%     pitch_wb_R_seq_bins_mean_all+1.96*pitch_wb_R_seq_bins_ste_all,t_wb_seq_bins_mean_all,[1 0 0])

xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('rotation angle c1 [deg]')
ylim([0 180])
set(gca,'YTick',[-90:45:180])

saveas(gca,['MSfig_timeseries_WBangles_c1.fig'])
saveas(gca,['MSfig_timeseries_WBangles_c1.png'])
plot2svg(['MSfig_timeseries_WBangles_c1.svg'])

%% c2
figure

% stroke
subplot(3,1,1)
hold on
% 
% ciplot(stroke_wb_L_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0)-1.96*stroke_wb_L_seq_bins_ste_all(isnan(t_wb_seq_bins_mean_all)==0),...
%     stroke_wb_L_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0)+1.96*stroke_wb_L_seq_bins_ste_all(isnan(t_wb_seq_bins_mean_all)==0),...
%     t_wb_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0),[0 .5 1])
% ciplot(stroke_wb_R_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0)-1.96*stroke_wb_R_seq_bins_ste_all(isnan(t_wb_seq_bins_mean_all)==0),...
%     stroke_wb_R_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0)+1.96*stroke_wb_R_seq_bins_ste_all(isnan(t_wb_seq_bins_mean_all)==0),...
%     t_wb_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0),[1 0 0])

plot(t_wb_seq_bins_mean_all,stroke_wb_L_seq_bins_mean_all_c2,'-b','linewidth',.5)
plot(t_wb_seq_bins_mean_all,stroke_wb_R_seq_bins_mean_all_c2,'-r','linewidth',.5)

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('stroke angle c2 [deg]')
ylim([-90 90])
set(gca,'YTick',[-90:45:90])

% deviation
subplot(3,1,2)
hold on

% ciplot(dev_wb_L_seq_bins_mean_all-1.96*dev_wb_L_seq_bins_ste_all,...
%     dev_wb_L_seq_bins_mean_all+1.96*dev_wb_L_seq_bins_ste_all,t_wb_seq_bins_mean_all,[0 .5 1])
% ciplot(dev_wb_R_seq_bins_mean_all-1.96*dev_wb_R_seq_bins_ste_all,...
%     dev_wb_R_seq_bins_mean_all+1.96*dev_wb_R_seq_bins_ste_all,t_wb_seq_bins_mean_all,[1 0 0])
% 
plot(t_wb_seq_bins_mean_all,dev_wb_L_seq_bins_mean_all_c2,'-b','linewidth',.5)
plot(t_wb_seq_bins_mean_all,dev_wb_R_seq_bins_mean_all_c2,'-r','linewidth',.5)

% xlabel('time [sec]')
xlim([-.04 .04])
ylabel('deviation angle c2 [deg]')
ylim([-30 15])
set(gca,'YTick',[-90:30:90])

% rotation
subplot(3,1,3)
hold on
% plot(t_wb_seq_bins_mean_all,pitch_wb_L_seq_bins_all,'b','linewidth',.5)
% plot(t_wb_seq_bins_mean_all,pitch_wb_R_seq_bins_all,'r','linewidth',.5)

plot(t_wb_seq_bins_mean_all,pitch_wb_L_seq_bins_mean_all_c2,'-b','linewidth',.5)
plot(t_wb_seq_bins_mean_all,pitch_wb_R_seq_bins_mean_all_c2,'-r','linewidth',.5)
% 
% ciplot(pitch_wb_L_seq_bins_mean_all-1.96*pitch_wb_L_seq_bins_ste_all,...
%     pitch_wb_L_seq_bins_mean_all+1.96*pitch_wb_L_seq_bins_ste_all,t_wb_seq_bins_mean_all,[0 .5 1])
% ciplot(pitch_wb_R_seq_bins_mean_all-1.96*pitch_wb_R_seq_bins_ste_all,...
%     pitch_wb_R_seq_bins_mean_all+1.96*pitch_wb_R_seq_bins_ste_all,t_wb_seq_bins_mean_all,[1 0 0])

xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('rotation angle c2 [deg]')
ylim([0 180])
set(gca,'YTick',[-90:45:180])

saveas(gca,['MSfig_timeseries_WBangles_c2.fig'])
saveas(gca,['MSfig_timeseries_WBangles_c2.png'])
plot2svg(['MSfig_timeseries_WBangles_c2.svg'])

%% c3
figure

% stroke
subplot(3,1,1)
hold on
% 
% ciplot(stroke_wb_L_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0)-1.96*stroke_wb_L_seq_bins_ste_all(isnan(t_wb_seq_bins_mean_all)==0),...
%     stroke_wb_L_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0)+1.96*stroke_wb_L_seq_bins_ste_all(isnan(t_wb_seq_bins_mean_all)==0),...
%     t_wb_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0),[0 .5 1])
% ciplot(stroke_wb_R_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0)-1.96*stroke_wb_R_seq_bins_ste_all(isnan(t_wb_seq_bins_mean_all)==0),...
%     stroke_wb_R_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0)+1.96*stroke_wb_R_seq_bins_ste_all(isnan(t_wb_seq_bins_mean_all)==0),...
%     t_wb_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0),[1 0 0])

plot(t_wb_seq_bins_mean_all,stroke_wb_L_seq_bins_mean_all_c3,'-b','linewidth',.5)
plot(t_wb_seq_bins_mean_all,stroke_wb_R_seq_bins_mean_all_c3,'-r','linewidth',.5)

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('stroke angle c3 [deg]')
ylim([-90 90])
set(gca,'YTick',[-90:45:90])

% deviation
subplot(3,1,2)
hold on

% ciplot(dev_wb_L_seq_bins_mean_all-1.96*dev_wb_L_seq_bins_ste_all,...
%     dev_wb_L_seq_bins_mean_all+1.96*dev_wb_L_seq_bins_ste_all,t_wb_seq_bins_mean_all,[0 .5 1])
% ciplot(dev_wb_R_seq_bins_mean_all-1.96*dev_wb_R_seq_bins_ste_all,...
%     dev_wb_R_seq_bins_mean_all+1.96*dev_wb_R_seq_bins_ste_all,t_wb_seq_bins_mean_all,[1 0 0])
% 
plot(t_wb_seq_bins_mean_all,dev_wb_L_seq_bins_mean_all_c3,'-b','linewidth',.5)
plot(t_wb_seq_bins_mean_all,dev_wb_R_seq_bins_mean_all_c3,'-r','linewidth',.5)

% xlabel('time [sec]')
xlim([-.04 .04])
ylabel('deviation angle c3 [deg]')
ylim([-30 15])
set(gca,'YTick',[-90:30:90])

% rotation
subplot(3,1,3)
hold on
% plot(t_wb_seq_bins_mean_all,pitch_wb_L_seq_bins_all,'b','linewidth',.5)
% plot(t_wb_seq_bins_mean_all,pitch_wb_R_seq_bins_all,'r','linewidth',.5)

plot(t_wb_seq_bins_mean_all,pitch_wb_L_seq_bins_mean_all_c3,'-b','linewidth',.5)
plot(t_wb_seq_bins_mean_all,pitch_wb_R_seq_bins_mean_all_c3,'-r','linewidth',.5)
% 
% ciplot(pitch_wb_L_seq_bins_mean_all-1.96*pitch_wb_L_seq_bins_ste_all,...
%     pitch_wb_L_seq_bins_mean_all+1.96*pitch_wb_L_seq_bins_ste_all,t_wb_seq_bins_mean_all,[0 .5 1])
% ciplot(pitch_wb_R_seq_bins_mean_all-1.96*pitch_wb_R_seq_bins_ste_all,...
%     pitch_wb_R_seq_bins_mean_all+1.96*pitch_wb_R_seq_bins_ste_all,t_wb_seq_bins_mean_all,[1 0 0])

xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('rotation angle c3 [deg]')
ylim([0 180])
set(gca,'YTick',[-90:45:180])

saveas(gca,['MSfig_timeseries_WBangles_c3.fig'])
saveas(gca,['MSfig_timeseries_WBangles_c3.png'])
plot2svg(['MSfig_timeseries_WBangles_c3.svg'])

%% c4
figure

% stroke
subplot(3,1,1)
hold on
% 
% ciplot(stroke_wb_L_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0)-1.96*stroke_wb_L_seq_bins_ste_all(isnan(t_wb_seq_bins_mean_all)==0),...
%     stroke_wb_L_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0)+1.96*stroke_wb_L_seq_bins_ste_all(isnan(t_wb_seq_bins_mean_all)==0),...
%     t_wb_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0),[0 .5 1])
% ciplot(stroke_wb_R_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0)-1.96*stroke_wb_R_seq_bins_ste_all(isnan(t_wb_seq_bins_mean_all)==0),...
%     stroke_wb_R_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0)+1.96*stroke_wb_R_seq_bins_ste_all(isnan(t_wb_seq_bins_mean_all)==0),...
%     t_wb_seq_bins_mean_all(isnan(t_wb_seq_bins_mean_all)==0),[1 0 0])

plot(t_wb_seq_bins_mean_all,stroke_wb_L_seq_bins_mean_all_c4,'-b','linewidth',.5)
plot(t_wb_seq_bins_mean_all,stroke_wb_R_seq_bins_mean_all_c4,'-r','linewidth',.5)

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('stroke angle c4 [deg]')
ylim([-90 90])
set(gca,'YTick',[-90:45:90])

% deviation
subplot(3,1,2)
hold on

% ciplot(dev_wb_L_seq_bins_mean_all-1.96*dev_wb_L_seq_bins_ste_all,...
%     dev_wb_L_seq_bins_mean_all+1.96*dev_wb_L_seq_bins_ste_all,t_wb_seq_bins_mean_all,[0 .5 1])
% ciplot(dev_wb_R_seq_bins_mean_all-1.96*dev_wb_R_seq_bins_ste_all,...
%     dev_wb_R_seq_bins_mean_all+1.96*dev_wb_R_seq_bins_ste_all,t_wb_seq_bins_mean_all,[1 0 0])
% 
plot(t_wb_seq_bins_mean_all,dev_wb_L_seq_bins_mean_all_c4,'-b','linewidth',.5)
plot(t_wb_seq_bins_mean_all,dev_wb_R_seq_bins_mean_all_c4,'-r','linewidth',.5)

% xlabel('time [sec]')
xlim([-.04 .04])
ylabel('deviation angle c4 [deg]')
ylim([-30 15])
set(gca,'YTick',[-90:30:90])

% rotation
subplot(3,1,3)
hold on
% plot(t_wb_seq_bins_mean_all,pitch_wb_L_seq_bins_all,'b','linewidth',.5)
% plot(t_wb_seq_bins_mean_all,pitch_wb_R_seq_bins_all,'r','linewidth',.5)

plot(t_wb_seq_bins_mean_all,pitch_wb_L_seq_bins_mean_all_c4,'-b','linewidth',.5)
plot(t_wb_seq_bins_mean_all,pitch_wb_R_seq_bins_mean_all_c4,'-r','linewidth',.5)
% 
% ciplot(pitch_wb_L_seq_bins_mean_all-1.96*pitch_wb_L_seq_bins_ste_all,...
%     pitch_wb_L_seq_bins_mean_all+1.96*pitch_wb_L_seq_bins_ste_all,t_wb_seq_bins_mean_all,[0 .5 1])
% ciplot(pitch_wb_R_seq_bins_mean_all-1.96*pitch_wb_R_seq_bins_ste_all,...
%     pitch_wb_R_seq_bins_mean_all+1.96*pitch_wb_R_seq_bins_ste_all,t_wb_seq_bins_mean_all,[1 0 0])

xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('rotation angle c4 [deg]')
ylim([0 180])
set(gca,'YTick',[-90:45:180])

saveas(gca,['MSfig_timeseries_WBangles_c4.fig'])
saveas(gca,['MSfig_timeseries_WBangles_c4.png'])
plot2svg(['MSfig_timeseries_WBangles_c4.svg'])


%% WB parameters cluster ALL

%% stroke
figure

% L+R
subplot(3,1,1)
hold on
% plot(t_wb_seq_mean_all,Astroke_wb_seq_bins,'-g.')
plot(t_wb_seq_mean_all,Astroke_wb_R_seq_bins,'-r.')
plot(t_wb_seq_mean_all,Astroke_wb_L_seq_bins,'-b.')
% plot(t_wb_seq_mean_all,Astroke_wb_seq_bins,'-g.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('stroke angle amplitude [deg]')
ylim([120 165])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:15:180])

% L+R
subplot(3,1,2)
hold on
% plot(t_wb_seq_mean_all,stroke_wb_seq_bins_MAX,'-k.','color',[1 .5 0])
% plot(t_wb_seq_mean_all,-stroke_wb_seq_bins_MIN,'-k.','color',[0 .5 1])
plot(t_wb_seq_mean_all,mean([stroke_wb_L_seq_bins_MIN,stroke_wb_L_seq_bins_MAX]')','-r.')
plot(t_wb_seq_mean_all,mean([stroke_wb_R_seq_bins_MIN,stroke_wb_R_seq_bins_MAX]')','-b.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('mean stroke angle [deg]')
ylim([0 15])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:7.5:180])

% L-R
subplot(3,1,3)
hold on
plot(t_wb_seq_mean_all,Dstroke_wb_seq_bins_MAX,'-k.','color',[1 .5 0])
plot(t_wb_seq_mean_all,Dstroke_wb_seq_bins_MIN,'-k.','color',[0 .5 1])
% plot(t_wb_seq_mean_all,-dAstroke_wb_seq_bins,'-g.')

xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('stroke angle L-R [deg]')
ylim([-7.5 7.5])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:7.5:90])

saveas(gca,['MSfig_timeseries_strokeParams.fig'])
saveas(gca,['MSfig_timeseries_strokeParams.png'])
plot2svg(['MSfig_timeseries_strokeParams.svg'])

%% pitch
figure

% Downstroke
subplot(3,1,1)
hold on
plot(t_wb_seq_mean_all,pitch_wb_R_seq_bins_MAXmidDS,'-r.')
plot(t_wb_seq_mean_all,pitch_wb_L_seq_bins_MAXmidDS,'-b.')
% plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MAXmidDS,'-g.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('rotation angle downstroke [deg]')
ylim([45 60])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:5:90])

% Upstroke
subplot(3,1,2)
hold on
plot(t_wb_seq_mean_all,pitch_wb_R_seq_bins_MINmidUS,'-r.')
plot(t_wb_seq_mean_all,pitch_wb_L_seq_bins_MINmidUS,'-b.')
% plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MINmidUS,'-g.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('rotation angle upstroke [deg]')
ylim([-55 -40])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:5:90])

% L-R
subplot(3,1,3)
hold on
plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MAXmidDS,'-k.','color',[1 .5 0])
plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MINmidUS,'-k.','color',[0 .5 1])
% plot(t_wb_seq_mean_all,-dApitch_wb_seq_bins,'-g.')

xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('rotation angle L-R [deg]')
ylim([-7.5 7.5])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:5:90])

saveas(gca,['MSfig_timeseries_rotationParams.fig'])
saveas(gca,['MSfig_timeseries_rotationParams.png'])
plot2svg(['MSfig_timeseries_rotationParams.svg'])

%% dev
figure

% Downstroke
subplot(3,1,1)
hold on
% plot(t_wb_seq_mean_all,dev_wb_R_seq_bins_MAXds,'-c.')
% plot(t_wb_seq_mean_all,dev_wb_L_seq_bins_MAXds,'-m.')
% % plot(t_wb_seq_mean_all,Ddev_wb_seq_bins_MAXds,'-g.')

plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_DS,'-r.')
plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_DS,'-b.')
% 
% plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_US,'-c.')
% plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_US,'-m.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('deviation amplitude downstroke [deg]')
ylim([15 30])
set(gca,'YTick',[-90:5:90])

% Upstroke
subplot(3,1,2)
hold on
% plot(t_wb_seq_mean_all,dev_wb_R_seq_bins_MAXds,'-c.')
% plot(t_wb_seq_mean_all,dev_wb_L_seq_bins_MAXds,'-m.')
% % plot(t_wb_seq_mean_all,Ddev_wb_seq_bins_MAXds,'-g.')
% 
% plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_DS,'-k.','color',[1 .5 0])
% plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_DS,'-k.','color',[0 .5 1])

plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_US,'-r.')
plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_US,'-b.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('deviation amplitude upstroke [deg]')
ylim([5 20])
set(gca,'YTick',[-90:5:90])

% L-R
subplot(3,1,3)
hold on
plot(t_wb_seq_mean_all,dAdev_wb_seq_bins_DS,'-k.','color',[1 .5 0])
plot(t_wb_seq_mean_all,dAdev_wb_seq_bins_US,'-k.','color',[0 .5 1])
% plot(t_wb_seq_mean_all,-dAdev_wb_seq_bins,'-g.')

xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('deviation amplitude L-R [deg]')
ylim([-7.5 7.5])
set(gca,'YTick',[-90:5:90])

saveas(gca,['MSfig_timeseries_deviationParams.fig'])
saveas(gca,['MSfig_timeseries_deviationParams.png'])
plot2svg(['MSfig_timeseries_deviationParams.svg'])

%% WB parameters cluster c1

%% stroke
figure

% L+R
subplot(3,1,1)
hold on
% plot(t_wb_seq_mean_all,Astroke_wb_seq_bins_c1,'-g.')
plot(t_wb_seq_mean_all,Astroke_wb_R_seq_bins_c1,'-r.')
plot(t_wb_seq_mean_all,Astroke_wb_L_seq_bins_c1,'-b.')
% plot(t_wb_seq_mean_all,Astroke_wb_seq_bins_c1,'-g.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('stroke angle amplitude [deg]')
ylim([120 165])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:15:180])

% L+R
subplot(3,1,2)
hold on
% plot(t_wb_seq_mean_all,stroke_wb_seq_bins_MAX_c1,'-k.','color',[1 .5 0])
% plot(t_wb_seq_mean_all,-stroke_wb_seq_bins_MIN_c1,'-k.','color',[0 .5 1])
plot(t_wb_seq_mean_all,mean([stroke_wb_L_seq_bins_MIN_c1,stroke_wb_L_seq_bins_MAX_c1]')','-r.')
plot(t_wb_seq_mean_all,mean([stroke_wb_R_seq_bins_MIN_c1,stroke_wb_R_seq_bins_MAX_c1]')','-b.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('mean stroke angle [deg]')
ylim([0 15])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:7.5:180])

% L-R
subplot(3,1,3)
hold on
plot(t_wb_seq_mean_all,Dstroke_wb_seq_bins_MAX_c1,'-k.','color',[1 .5 0])
plot(t_wb_seq_mean_all,Dstroke_wb_seq_bins_MIN_c1,'-k.','color',[0 .5 1])
% plot(t_wb_seq_mean_all,-dAstroke_wb_seq_bins_c1,'-g.')

xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('stroke angle L-R [deg]')
ylim([-7.5 7.5])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:7.5:90])

saveas(gca,['MSfig_timeseries_strokeParams_c1.fig'])
saveas(gca,['MSfig_timeseries_strokeParams_c1.png'])
plot2svg(['MSfig_timeseries_strokeParams_c1.svg'])

%% pitch
figure

% Downstroke
subplot(3,1,1)
hold on
plot(t_wb_seq_mean_all,pitch_wb_R_seq_bins_MAXmidDS_c1,'-r.')
plot(t_wb_seq_mean_all,pitch_wb_L_seq_bins_MAXmidDS_c1,'-b.')
% plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MAXmidDS_c1,'-g.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('rotation angle downstroke c1 [deg]')
ylim([45 60])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:5:90])

% Upstroke
subplot(3,1,2)
hold on
plot(t_wb_seq_mean_all,pitch_wb_R_seq_bins_MINmidUS_c1,'-r.')
plot(t_wb_seq_mean_all,pitch_wb_L_seq_bins_MINmidUS_c1,'-b.')
% plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MINmidUS_c1,'-g.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('rotation angle upstroke c1 [deg]')
ylim([-55 -40])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:5:90])

% L-R
subplot(3,1,3)
hold on
plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MAXmidDS_c1,'-k.','color',[1 .5 0])
plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MINmidUS_c1,'-k.','color',[0 .5 1])
% plot(t_wb_seq_mean_all,-dApitch_wb_seq_bins_c1,'-g.')

xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('rotation angle L-R c1 [deg]')
ylim([-7.5 7.5])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:5:90])

saveas(gca,['MSfig_timeseries_rotationParams_c1.fig'])
saveas(gca,['MSfig_timeseries_rotationParams_c1.png'])
plot2svg(['MSfig_timeseries_rotationParams_c1.svg'])

%% dev
figure

% Downstroke
subplot(3,1,1)
hold on
% plot(t_wb_seq_mean_all,dev_wb_R_seq_bins_MAXds_c1,'-c.')
% plot(t_wb_seq_mean_all,dev_wb_L_seq_bins_MAXds_c1,'-m.')
% % plot(t_wb_seq_mean_all,Ddev_wb_seq_bins_MAXds_c1,'-g.')

plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_DS_c1,'-r.')
plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_DS_c1,'-b.')
% 
% plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_US_c1,'-c.')
% plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_US_c1,'-m.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('deviation amplitude downstroke c1 [deg]')
ylim([15 30])
set(gca,'YTick',[-90:5:90])

% Upstroke
subplot(3,1,2)
hold on
% plot(t_wb_seq_mean_all,dev_wb_R_seq_bins_MAXds_c1,'-c.')
% plot(t_wb_seq_mean_all,dev_wb_L_seq_bins_MAXds_c1,'-m.')
% % plot(t_wb_seq_mean_all,Ddev_wb_seq_bins_MAXds_c1,'-g.')
% 
% plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_DS_c1,'-k.','color',[1 .5 0])
% plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_DS_c1,'-k.','color',[0 .5 1])

plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_US_c1,'-r.')
plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_US_c1,'-b.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('deviation amplitude upstroke c1 [deg]')
ylim([5 20])
set(gca,'YTick',[-90:5:90])

% L-R
subplot(3,1,3)
hold on
plot(t_wb_seq_mean_all,dAdev_wb_seq_bins_DS_c1,'-k.','color',[1 .5 0])
plot(t_wb_seq_mean_all,dAdev_wb_seq_bins_US_c1,'-k.','color',[0 .5 1])
% plot(t_wb_seq_mean_all,-dAdev_wb_seq_bins_c1,'-g.')

xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('deviation amplitude L-R c1 [deg]')
ylim([-7.5 7.5])
set(gca,'YTick',[-90:5:90])

saveas(gca,['MSfig_timeseries_deviationParams_c1.fig'])
saveas(gca,['MSfig_timeseries_deviationParams_c1.png'])
plot2svg(['MSfig_timeseries_deviationParams_c1.svg'])

%% WB parameters cluster c2

%% stroke
figure

% L+R
subplot(3,1,1)
hold on
% plot(t_wb_seq_mean_all,Astroke_wb_seq_bins_c2,'-g.')
plot(t_wb_seq_mean_all,Astroke_wb_R_seq_bins_c2,'-r.')
plot(t_wb_seq_mean_all,Astroke_wb_L_seq_bins_c2,'-b.')
% plot(t_wb_seq_mean_all,Astroke_wb_seq_bins_c2,'-g.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('stroke angle amplitude [deg]')
ylim([120 165])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:15:180])

% L+R
subplot(3,1,2)
hold on
% plot(t_wb_seq_mean_all,stroke_wb_seq_bins_MAX_c2,'-k.','color',[1 .5 0])
% plot(t_wb_seq_mean_all,-stroke_wb_seq_bins_MIN_c2,'-k.','color',[0 .5 1])
plot(t_wb_seq_mean_all,mean([stroke_wb_L_seq_bins_MIN_c2,stroke_wb_L_seq_bins_MAX_c2]')','-r.')
plot(t_wb_seq_mean_all,mean([stroke_wb_R_seq_bins_MIN_c2,stroke_wb_R_seq_bins_MAX_c2]')','-b.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('mean stroke angle [deg]')
ylim([0 15])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:7.5:180])

% L-R
subplot(3,1,3)
hold on
plot(t_wb_seq_mean_all,Dstroke_wb_seq_bins_MAX_c2,'-k.','color',[1 .5 0])
plot(t_wb_seq_mean_all,Dstroke_wb_seq_bins_MIN_c2,'-k.','color',[0 .5 1])
% plot(t_wb_seq_mean_all,-dAstroke_wb_seq_bins_c2,'-g.')

xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('stroke angle L-R [deg]')
ylim([-7.5 7.5])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:7.5:90])

saveas(gca,['MSfig_timeseries_strokeParams_c2.fig'])
saveas(gca,['MSfig_timeseries_strokeParams_c2.png'])
plot2svg(['MSfig_timeseries_strokeParams_c2.svg'])

%% pitch
figure

% Downstroke
subplot(3,1,1)
hold on
plot(t_wb_seq_mean_all,pitch_wb_R_seq_bins_MAXmidDS_c2,'-r.')
plot(t_wb_seq_mean_all,pitch_wb_L_seq_bins_MAXmidDS_c2,'-b.')
% plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MAXmidDS_c2,'-g.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('rotation angle downstroke c2 [deg]')
ylim([45 60])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:5:90])

% Upstroke
subplot(3,1,2)
hold on
plot(t_wb_seq_mean_all,pitch_wb_R_seq_bins_MINmidUS_c2,'-r.')
plot(t_wb_seq_mean_all,pitch_wb_L_seq_bins_MINmidUS_c2,'-b.')
% plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MINmidUS_c2,'-g.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('rotation angle upstroke c2 [deg]')
ylim([-55 -40])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:5:90])

% L-R
subplot(3,1,3)
hold on
plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MAXmidDS_c2,'-k.','color',[1 .5 0])
plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MINmidUS_c2,'-k.','color',[0 .5 1])
% plot(t_wb_seq_mean_all,-dApitch_wb_seq_bins_c2,'-g.')

xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('rotation angle L-R c2 [deg]')
ylim([-7.5 7.5])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:5:90])

saveas(gca,['MSfig_timeseries_rotationParams_c2.fig'])
saveas(gca,['MSfig_timeseries_rotationParams_c2.png'])
plot2svg(['MSfig_timeseries_rotationParams_c2.svg'])

%% dev
figure

% Downstroke
subplot(3,1,1)
hold on
% plot(t_wb_seq_mean_all,dev_wb_R_seq_bins_MAXds_c2,'-c.')
% plot(t_wb_seq_mean_all,dev_wb_L_seq_bins_MAXds_c2,'-m.')
% % plot(t_wb_seq_mean_all,Ddev_wb_seq_bins_MAXds_c2,'-g.')

plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_DS_c2,'-r.')
plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_DS_c2,'-b.')
% 
% plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_US_c2,'-c.')
% plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_US_c2,'-m.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('deviation amplitude downstroke c2 [deg]')
ylim([15 30])
set(gca,'YTick',[-90:5:90])

% Upstroke
subplot(3,1,2)
hold on
% plot(t_wb_seq_mean_all,dev_wb_R_seq_bins_MAXds_c2,'-c.')
% plot(t_wb_seq_mean_all,dev_wb_L_seq_bins_MAXds_c2,'-m.')
% % plot(t_wb_seq_mean_all,Ddev_wb_seq_bins_MAXds_c2,'-g.')
% 
% plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_DS_c2,'-k.','color',[1 .5 0])
% plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_DS_c2,'-k.','color',[0 .5 1])

plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_US_c2,'-r.')
plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_US_c2,'-b.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('deviation amplitude upstroke c2 [deg]')
ylim([5 20])
set(gca,'YTick',[-90:5:90])

% L-R
subplot(3,1,3)
hold on
plot(t_wb_seq_mean_all,dAdev_wb_seq_bins_DS_c2,'-k.','color',[1 .5 0])
plot(t_wb_seq_mean_all,dAdev_wb_seq_bins_US_c2,'-k.','color',[0 .5 1])
% plot(t_wb_seq_mean_all,-dAdev_wb_seq_bins_c2,'-g.')

xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('deviation amplitude L-R c2 [deg]')
ylim([-7.5 7.5])
set(gca,'YTick',[-90:5:90])

saveas(gca,['MSfig_timeseries_deviationParams_c2.fig'])
saveas(gca,['MSfig_timeseries_deviationParams_c2.png'])
plot2svg(['MSfig_timeseries_deviationParams_c2.svg'])

%% WB parameters cluster c3

%% stroke
figure

% L+R
subplot(3,1,1)
hold on
% plot(t_wb_seq_mean_all,Astroke_wb_seq_bins_c3,'-g.')
plot(t_wb_seq_mean_all,Astroke_wb_R_seq_bins_c3,'-r.')
plot(t_wb_seq_mean_all,Astroke_wb_L_seq_bins_c3,'-b.')
% plot(t_wb_seq_mean_all,Astroke_wb_seq_bins_c3,'-g.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('stroke angle amplitude [deg]')
ylim([120 165])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:15:180])

% L+R
subplot(3,1,2)
hold on
% plot(t_wb_seq_mean_all,stroke_wb_seq_bins_MAX_c3,'-k.','color',[1 .5 0])
% plot(t_wb_seq_mean_all,-stroke_wb_seq_bins_MIN_c3,'-k.','color',[0 .5 1])
plot(t_wb_seq_mean_all,mean([stroke_wb_L_seq_bins_MIN_c3,stroke_wb_L_seq_bins_MAX_c3]')','-r.')
plot(t_wb_seq_mean_all,mean([stroke_wb_R_seq_bins_MIN_c3,stroke_wb_R_seq_bins_MAX_c3]')','-b.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('mean stroke angle [deg]')
ylim([0 15])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:7.5:180])

% L-R
subplot(3,1,3)
hold on
plot(t_wb_seq_mean_all,Dstroke_wb_seq_bins_MAX_c3,'-k.','color',[1 .5 0])
plot(t_wb_seq_mean_all,Dstroke_wb_seq_bins_MIN_c3,'-k.','color',[0 .5 1])
% plot(t_wb_seq_mean_all,-dAstroke_wb_seq_bins_c3,'-g.')

xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('stroke angle L-R [deg]')
ylim([-7.5 7.5])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:7.5:90])

saveas(gca,['MSfig_timeseries_strokeParams_c3.fig'])
saveas(gca,['MSfig_timeseries_strokeParams_c3.png'])
plot2svg(['MSfig_timeseries_strokeParams_c3.svg'])

%% pitch
figure

% Downstroke
subplot(3,1,1)
hold on
plot(t_wb_seq_mean_all,pitch_wb_R_seq_bins_MAXmidDS_c3,'-r.')
plot(t_wb_seq_mean_all,pitch_wb_L_seq_bins_MAXmidDS_c3,'-b.')
% plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MAXmidDS_c3,'-g.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('rotation angle downstroke c3 [deg]')
ylim([45 60])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:5:90])

% Upstroke
subplot(3,1,2)
hold on
plot(t_wb_seq_mean_all,pitch_wb_R_seq_bins_MINmidUS_c3,'-r.')
plot(t_wb_seq_mean_all,pitch_wb_L_seq_bins_MINmidUS_c3,'-b.')
% plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MINmidUS_c3,'-g.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('rotation angle upstroke c3 [deg]')
ylim([-55 -40])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:5:90])

% L-R
subplot(3,1,3)
hold on
plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MAXmidDS_c3,'-k.','color',[1 .5 0])
plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MINmidUS_c3,'-k.','color',[0 .5 1])
% plot(t_wb_seq_mean_all,-dApitch_wb_seq_bins_c3,'-g.')

xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('rotation angle L-R c3 [deg]')
ylim([-7.5 7.5])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:5:90])

saveas(gca,['MSfig_timeseries_rotationParams_c3.fig'])
saveas(gca,['MSfig_timeseries_rotationParams_c3.png'])
plot2svg(['MSfig_timeseries_rotationParams_c3.svg'])

%% dev
figure

% Downstroke
subplot(3,1,1)
hold on
% plot(t_wb_seq_mean_all,dev_wb_R_seq_bins_MAXds_c3,'-c.')
% plot(t_wb_seq_mean_all,dev_wb_L_seq_bins_MAXds_c3,'-m.')
% % plot(t_wb_seq_mean_all,Ddev_wb_seq_bins_MAXds_c3,'-g.')

plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_DS_c3,'-r.')
plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_DS_c3,'-b.')
% 
% plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_US_c3,'-c.')
% plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_US_c3,'-m.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('deviation amplitude downstroke c3 [deg]')
ylim([15 30])
set(gca,'YTick',[-90:5:90])

% Upstroke
subplot(3,1,2)
hold on
% plot(t_wb_seq_mean_all,dev_wb_R_seq_bins_MAXds_c3,'-c.')
% plot(t_wb_seq_mean_all,dev_wb_L_seq_bins_MAXds_c3,'-m.')
% % plot(t_wb_seq_mean_all,Ddev_wb_seq_bins_MAXds_c3,'-g.')
% 
% plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_DS_c3,'-k.','color',[1 .5 0])
% plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_DS_c3,'-k.','color',[0 .5 1])

plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_US_c3,'-r.')
plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_US_c3,'-b.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('deviation amplitude upstroke c3 [deg]')
ylim([5 20])
set(gca,'YTick',[-90:5:90])

% L-R
subplot(3,1,3)
hold on
plot(t_wb_seq_mean_all,dAdev_wb_seq_bins_DS_c3,'-k.','color',[1 .5 0])
plot(t_wb_seq_mean_all,dAdev_wb_seq_bins_US_c3,'-k.','color',[0 .5 1])
% plot(t_wb_seq_mean_all,-dAdev_wb_seq_bins_c3,'-g.')

xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('deviation amplitude L-R c3 [deg]')
ylim([-7.5 7.5])
set(gca,'YTick',[-90:5:90])

saveas(gca,['MSfig_timeseries_deviationParams_c3.fig'])
saveas(gca,['MSfig_timeseries_deviationParams_c3.png'])
plot2svg(['MSfig_timeseries_deviationParams_c3.svg'])

%% WB parameters cluster c4

%% stroke
figure

% L+R
subplot(3,1,1)
hold on
% plot(t_wb_seq_mean_all,Astroke_wb_seq_bins_c4,'-g.')
plot(t_wb_seq_mean_all,Astroke_wb_R_seq_bins_c4,'-r.')
plot(t_wb_seq_mean_all,Astroke_wb_L_seq_bins_c4,'-b.')
% plot(t_wb_seq_mean_all,Astroke_wb_seq_bins_c4,'-g.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('stroke angle amplitude [deg]')
ylim([120 165])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:15:180])

% L+R
subplot(3,1,2)
hold on
% plot(t_wb_seq_mean_all,stroke_wb_seq_bins_MAX_c4,'-k.','color',[1 .5 0])
% plot(t_wb_seq_mean_all,-stroke_wb_seq_bins_MIN_c4,'-k.','color',[0 .5 1])
plot(t_wb_seq_mean_all,mean([stroke_wb_L_seq_bins_MIN_c4,stroke_wb_L_seq_bins_MAX_c4]')','-r.')
plot(t_wb_seq_mean_all,mean([stroke_wb_R_seq_bins_MIN_c4,stroke_wb_R_seq_bins_MAX_c4]')','-b.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('mean stroke angle [deg]')
ylim([0 15])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:7.5:180])

% L-R
subplot(3,1,3)
hold on
plot(t_wb_seq_mean_all,Dstroke_wb_seq_bins_MAX_c4,'-k.','color',[1 .5 0])
plot(t_wb_seq_mean_all,Dstroke_wb_seq_bins_MIN_c4,'-k.','color',[0 .5 1])
% plot(t_wb_seq_mean_all,-dAstroke_wb_seq_bins_c4,'-g.')

xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('stroke angle L-R [deg]')
ylim([-7.5 7.5])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:7.5:90])

saveas(gca,['MSfig_timeseries_strokeParams_c4.fig'])
saveas(gca,['MSfig_timeseries_strokeParams_c4.png'])
plot2svg(['MSfig_timeseries_strokeParams_c4.svg'])

%% pitch
figure

% Downstroke
subplot(3,1,1)
hold on
plot(t_wb_seq_mean_all,pitch_wb_R_seq_bins_MAXmidDS_c4,'-r.')
plot(t_wb_seq_mean_all,pitch_wb_L_seq_bins_MAXmidDS_c4,'-b.')
% plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MAXmidDS_c4,'-g.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('rotation angle downstroke c4 [deg]')
ylim([45 60])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:5:90])

% Upstroke
subplot(3,1,2)
hold on
plot(t_wb_seq_mean_all,pitch_wb_R_seq_bins_MINmidUS_c4,'-r.')
plot(t_wb_seq_mean_all,pitch_wb_L_seq_bins_MINmidUS_c4,'-b.')
% plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MINmidUS_c4,'-g.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('rotation angle upstroke c4 [deg]')
ylim([-55 -40])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:5:90])

% L-R
subplot(3,1,3)
hold on
plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MAXmidDS_c4,'-k.','color',[1 .5 0])
plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MINmidUS_c4,'-k.','color',[0 .5 1])
% plot(t_wb_seq_mean_all,-dApitch_wb_seq_bins_c4,'-g.')

xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('rotation angle L-R c4 [deg]')
ylim([-7.5 7.5])
set(gca,'XTick',[-.04:.04:.04])
set(gca,'YTick',[-90:5:90])

saveas(gca,['MSfig_timeseries_rotationParams_c4.fig'])
saveas(gca,['MSfig_timeseries_rotationParams_c4.png'])
plot2svg(['MSfig_timeseries_rotationParams_c4.svg'])

%% dev
figure

% Downstroke
subplot(3,1,1)
hold on
% plot(t_wb_seq_mean_all,dev_wb_R_seq_bins_MAXds_c4,'-c.')
% plot(t_wb_seq_mean_all,dev_wb_L_seq_bins_MAXds_c4,'-m.')
% % plot(t_wb_seq_mean_all,Ddev_wb_seq_bins_MAXds_c4,'-g.')

plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_DS_c4,'-r.')
plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_DS_c4,'-b.')
% 
% plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_US_c4,'-c.')
% plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_US_c4,'-m.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('deviation amplitude downstroke c4 [deg]')
ylim([15 30])
set(gca,'YTick',[-90:5:90])

% Upstroke
subplot(3,1,2)
hold on
% plot(t_wb_seq_mean_all,dev_wb_R_seq_bins_MAXds_c4,'-c.')
% plot(t_wb_seq_mean_all,dev_wb_L_seq_bins_MAXds_c4,'-m.')
% % plot(t_wb_seq_mean_all,Ddev_wb_seq_bins_MAXds_c4,'-g.')
% 
% plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_DS_c4,'-k.','color',[1 .5 0])
% plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_DS_c4,'-k.','color',[0 .5 1])

plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_US_c4,'-r.')
plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_US_c4,'-b.')

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('deviation amplitude upstroke c4 [deg]')
ylim([5 20])
set(gca,'YTick',[-90:5:90])

% L-R
subplot(3,1,3)
hold on
plot(t_wb_seq_mean_all,dAdev_wb_seq_bins_DS_c4,'-k.','color',[1 .5 0])
plot(t_wb_seq_mean_all,dAdev_wb_seq_bins_US_c4,'-k.','color',[0 .5 1])
% plot(t_wb_seq_mean_all,-dAdev_wb_seq_bins_c4,'-g.')

xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('deviation amplitude L-R c4 [deg]')
ylim([-7.5 7.5])
set(gca,'YTick',[-90:5:90])

saveas(gca,['MSfig_timeseries_deviationParams_c4.fig'])
saveas(gca,['MSfig_timeseries_deviationParams_c4.png'])
plot2svg(['MSfig_timeseries_deviationParams_c4.svg'])



%% ACCELS
% figure
% 
% %% Fy accel
% subplot(3,1,1)
% hold on
% 
% % body dyn accel
% ciplot(roll_dot_dot_mean_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)-1.96*roll_dot_dot_mean_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),roll_dot_dot_mean_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)+1.96*roll_dot_dot_mean_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),t_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0))
% plot(t_wb_seq_mean_all,roll_dot_dot_mean_wb_seq_mean_all,'-k.','color',[0 .5 1])
% 
% % xlabel('time [sec]')
% ylabel('Aroll')
% xlim([-.04 .04]) 
% ylim([-3 3])
% set(gca,'XTick',[-.04:.02:.04])
% set(gca,'YTick',[-3:3:3])

%% Angular ACCELS
figure

%% roll accel
subplot(3,1,1)
hold on

% body dyn accel
ciplot(roll_dot_dot_mean_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)-1.96*roll_dot_dot_mean_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),roll_dot_dot_mean_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)+1.96*roll_dot_dot_mean_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),t_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0))
plot(t_wb_seq_mean_all,roll_dot_dot_mean_wb_seq_mean_all,'-k.')

% xlabel('time [sec]')
ylabel('Aroll')
xlim([-.04 .04]) 
ylim([-3 3])
set(gca,'XTick',[-.04:.02:.04])
set(gca,'YTick',[-3:3:3])

%% pitch accel
subplot(3,1,2)
hold on

% body dyn accel
ciplot(pitch_dot_dot_mean_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)-1.96*pitch_dot_dot_mean_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),pitch_dot_dot_mean_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)+1.96*pitch_dot_dot_mean_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),t_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0))
plot(t_wb_seq_mean_all,pitch_dot_dot_mean_wb_seq_mean_all,'-k.')

% xlabel('time [sec]')
ylabel('Apitch')
xlim([-.04 .04]) 
ylim([-3 3])
set(gca,'XTick',[-.04:.02:.04])
set(gca,'YTick',[-3:3:3])

%% yaw accel
subplot(3,1,3)
hold on

% body dyn accel
ciplot(yaw_dot_dot_mean_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)-1.96*yaw_dot_dot_mean_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),yaw_dot_dot_mean_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0)+1.96*yaw_dot_dot_mean_wb_seq_ste_all(isnan(t_wb_seq_mean_all)==0),t_wb_seq_mean_all(isnan(t_wb_seq_mean_all)==0))
plot(t_wb_seq_mean_all,yaw_dot_dot_mean_wb_seq_mean_all,'-k.')

xlabel('time [sec]')
ylabel('Ayaw')
xlim([-.04 .04]) 
ylim([-3 3])
set(gca,'XTick',[-.04:.02:.04])
set(gca,'YTick',[-3:3:3])

saveas(gca,['MSfig_timeseries_accels.fig'])
saveas(gca,['MSfig_timeseries_accels.png'])
plot2svg(['MSfig_timeseries_accels.svg'])


%% Angular ACCELS c1
figure

%% roll accel
subplot(3,1,1)
hold on

% body dyn accel
ciplot(roll_dot_dot_mean_wb_seq_mean_all_c1(isnan(roll_dot_dot_mean_wb_seq_mean_all_c1)==0)-1.96*roll_dot_dot_mean_wb_seq_ste_all_c1(isnan(roll_dot_dot_mean_wb_seq_mean_all_c1)==0),...
    roll_dot_dot_mean_wb_seq_mean_all_c1(isnan(roll_dot_dot_mean_wb_seq_mean_all_c1)==0)+1.96*roll_dot_dot_mean_wb_seq_ste_all_c1(isnan(roll_dot_dot_mean_wb_seq_mean_all_c1)==0),...
    t_wb_seq_mean_all(isnan(roll_dot_dot_mean_wb_seq_mean_all_c1)==0))
plot(t_wb_seq_mean_all,roll_dot_dot_mean_wb_seq_mean_all_c1,'-k.')

% xlabel('time [sec]')
ylabel('Aroll')
xlim([-.04 .04]) 
ylim([-3 3])
set(gca,'XTick',[-.04:.02:.04])
set(gca,'YTick',[-3:3:3])

%% pitch accel
subplot(3,1,2)
hold on

% body dyn accel
ciplot(pitch_dot_dot_mean_wb_seq_mean_all_c1(isnan(roll_dot_dot_mean_wb_seq_mean_all_c1)==0)-1.96*pitch_dot_dot_mean_wb_seq_ste_all_c1(isnan(roll_dot_dot_mean_wb_seq_mean_all_c1)==0),...
    pitch_dot_dot_mean_wb_seq_mean_all_c1(isnan(roll_dot_dot_mean_wb_seq_mean_all_c1)==0)+1.96*pitch_dot_dot_mean_wb_seq_ste_all_c1(isnan(roll_dot_dot_mean_wb_seq_mean_all_c1)==0),...
    t_wb_seq_mean_all(isnan(roll_dot_dot_mean_wb_seq_mean_all_c1)==0))
plot(t_wb_seq_mean_all,pitch_dot_dot_mean_wb_seq_mean_all_c1,'-k.')

% xlabel('time [sec]')
ylabel('Apitch')
xlim([-.04 .04]) 
ylim([-3 3])
set(gca,'XTick',[-.04:.02:.04])
set(gca,'YTick',[-3:3:3])

%% yaw accel
subplot(3,1,3)
hold on

% body dyn accel
ciplot(yaw_dot_dot_mean_wb_seq_mean_all_c1(isnan(roll_dot_dot_mean_wb_seq_mean_all_c1)==0)-1.96*yaw_dot_dot_mean_wb_seq_ste_all_c1(isnan(roll_dot_dot_mean_wb_seq_mean_all_c1)==0),...
    yaw_dot_dot_mean_wb_seq_mean_all_c1(isnan(roll_dot_dot_mean_wb_seq_mean_all_c1)==0)+1.96*yaw_dot_dot_mean_wb_seq_ste_all_c1(isnan(roll_dot_dot_mean_wb_seq_mean_all_c1)==0),...
    t_wb_seq_mean_all(isnan(roll_dot_dot_mean_wb_seq_mean_all_c1)==0))
plot(t_wb_seq_mean_all,yaw_dot_dot_mean_wb_seq_mean_all_c1,'-k.')

xlabel('time [sec]')
ylabel('Ayaw')
xlim([-.04 .04]) 
ylim([-3 3])
set(gca,'XTick',[-.04:.02:.04])
set(gca,'YTick',[-3:3:3])

saveas(gca,['MSfig_timeseries_accels_c1.fig'])
saveas(gca,['MSfig_timeseries_accels_c1.png'])
plot2svg(['MSfig_timeseries_accels_c1.svg'])


%% Angular ACCELS c2
figure

%% roll accel
subplot(3,1,1)
hold on

% body dyn accel
ciplot(roll_dot_dot_mean_wb_seq_mean_all_c2(isnan(roll_dot_dot_mean_wb_seq_mean_all_c2)==0)-1.96*roll_dot_dot_mean_wb_seq_ste_all_c2(isnan(roll_dot_dot_mean_wb_seq_mean_all_c2)==0),...
    roll_dot_dot_mean_wb_seq_mean_all_c2(isnan(roll_dot_dot_mean_wb_seq_mean_all_c2)==0)+1.96*roll_dot_dot_mean_wb_seq_ste_all_c2(isnan(roll_dot_dot_mean_wb_seq_mean_all_c2)==0),...
    t_wb_seq_mean_all(isnan(roll_dot_dot_mean_wb_seq_mean_all_c2)==0))
plot(t_wb_seq_mean_all,roll_dot_dot_mean_wb_seq_mean_all_c2,'-k.')

% xlabel('time [sec]')
ylabel('Aroll')
xlim([-.04 .04]) 
ylim([-3 3])
set(gca,'XTick',[-.04:.02:.04])
set(gca,'YTick',[-3:3:3])

%% pitch accel
subplot(3,1,2)
hold on

% body dyn accel
ciplot(pitch_dot_dot_mean_wb_seq_mean_all_c2(isnan(roll_dot_dot_mean_wb_seq_mean_all_c2)==0)-1.96*pitch_dot_dot_mean_wb_seq_ste_all_c2(isnan(roll_dot_dot_mean_wb_seq_mean_all_c2)==0),...
    pitch_dot_dot_mean_wb_seq_mean_all_c2(isnan(roll_dot_dot_mean_wb_seq_mean_all_c2)==0)+1.96*pitch_dot_dot_mean_wb_seq_ste_all_c2(isnan(roll_dot_dot_mean_wb_seq_mean_all_c2)==0),...
    t_wb_seq_mean_all(isnan(roll_dot_dot_mean_wb_seq_mean_all_c2)==0))
plot(t_wb_seq_mean_all,pitch_dot_dot_mean_wb_seq_mean_all_c2,'-k.')

% xlabel('time [sec]')
ylabel('Apitch')
xlim([-.04 .04]) 
ylim([-3 3])
set(gca,'XTick',[-.04:.02:.04])
set(gca,'YTick',[-3:3:3])

%% yaw accel
subplot(3,1,3)
hold on

% body dyn accel
ciplot(yaw_dot_dot_mean_wb_seq_mean_all_c2(isnan(roll_dot_dot_mean_wb_seq_mean_all_c2)==0)-1.96*yaw_dot_dot_mean_wb_seq_ste_all_c2(isnan(roll_dot_dot_mean_wb_seq_mean_all_c2)==0),...
    yaw_dot_dot_mean_wb_seq_mean_all_c2(isnan(roll_dot_dot_mean_wb_seq_mean_all_c2)==0)+1.96*yaw_dot_dot_mean_wb_seq_ste_all_c2(isnan(roll_dot_dot_mean_wb_seq_mean_all_c2)==0),...
    t_wb_seq_mean_all(isnan(roll_dot_dot_mean_wb_seq_mean_all_c2)==0))
plot(t_wb_seq_mean_all,yaw_dot_dot_mean_wb_seq_mean_all_c2,'-k.')

xlabel('time [sec]')
ylabel('Ayaw')
xlim([-.04 .04]) 
ylim([-3 3])
set(gca,'XTick',[-.04:.02:.04])
set(gca,'YTick',[-3:3:3])

saveas(gca,['MSfig_timeseries_accels_c2.fig'])
saveas(gca,['MSfig_timeseries_accels_c2.png'])
plot2svg(['MSfig_timeseries_accels_c2.svg'])


%% Angular ACCELS c3
figure

%% roll accel
subplot(3,1,1)
hold on

% body dyn accel
ciplot(roll_dot_dot_mean_wb_seq_mean_all_c3(isnan(roll_dot_dot_mean_wb_seq_mean_all_c3)==0)-1.96*roll_dot_dot_mean_wb_seq_ste_all_c3(isnan(roll_dot_dot_mean_wb_seq_mean_all_c3)==0),...
    roll_dot_dot_mean_wb_seq_mean_all_c3(isnan(roll_dot_dot_mean_wb_seq_mean_all_c3)==0)+1.96*roll_dot_dot_mean_wb_seq_ste_all_c3(isnan(roll_dot_dot_mean_wb_seq_mean_all_c3)==0),...
    t_wb_seq_mean_all(isnan(roll_dot_dot_mean_wb_seq_mean_all_c3)==0))
plot(t_wb_seq_mean_all,roll_dot_dot_mean_wb_seq_mean_all_c3,'-k.')

% xlabel('time [sec]')
ylabel('Aroll')
xlim([-.04 .04]) 
ylim([-3 3])
set(gca,'XTick',[-.04:.02:.04])
set(gca,'YTick',[-3:3:3])

%% pitch accel
subplot(3,1,2)
hold on

% body dyn accel
ciplot(pitch_dot_dot_mean_wb_seq_mean_all_c3(isnan(roll_dot_dot_mean_wb_seq_mean_all_c3)==0)-1.96*pitch_dot_dot_mean_wb_seq_ste_all_c3(isnan(roll_dot_dot_mean_wb_seq_mean_all_c3)==0),...
    pitch_dot_dot_mean_wb_seq_mean_all_c3(isnan(roll_dot_dot_mean_wb_seq_mean_all_c3)==0)+1.96*pitch_dot_dot_mean_wb_seq_ste_all_c3(isnan(roll_dot_dot_mean_wb_seq_mean_all_c3)==0),...
    t_wb_seq_mean_all(isnan(roll_dot_dot_mean_wb_seq_mean_all_c3)==0))
plot(t_wb_seq_mean_all,pitch_dot_dot_mean_wb_seq_mean_all_c3,'-k.')

% xlabel('time [sec]')
ylabel('Apitch')
xlim([-.04 .04]) 
ylim([-3 3])
set(gca,'XTick',[-.04:.02:.04])
set(gca,'YTick',[-3:3:3])

%% yaw accel
subplot(3,1,3)
hold on

% body dyn accel
ciplot(yaw_dot_dot_mean_wb_seq_mean_all_c3(isnan(roll_dot_dot_mean_wb_seq_mean_all_c3)==0)-1.96*yaw_dot_dot_mean_wb_seq_ste_all_c3(isnan(roll_dot_dot_mean_wb_seq_mean_all_c3)==0),...
    yaw_dot_dot_mean_wb_seq_mean_all_c3(isnan(roll_dot_dot_mean_wb_seq_mean_all_c3)==0)+1.96*yaw_dot_dot_mean_wb_seq_ste_all_c3(isnan(roll_dot_dot_mean_wb_seq_mean_all_c3)==0),...
    t_wb_seq_mean_all(isnan(roll_dot_dot_mean_wb_seq_mean_all_c3)==0))
plot(t_wb_seq_mean_all,yaw_dot_dot_mean_wb_seq_mean_all_c3,'-k.')

xlabel('time [sec]')
ylabel('Ayaw')
xlim([-.04 .04]) 
ylim([-3 3])
set(gca,'XTick',[-.04:.02:.04])
set(gca,'YTick',[-3:3:3])

saveas(gca,['MSfig_timeseries_accels_c3.fig'])
saveas(gca,['MSfig_timeseries_accels_c3.png'])
plot2svg(['MSfig_timeseries_accels_c3.svg'])


%% Angular ACCELS c4
figure

%% roll accel
subplot(3,1,1)
hold on

% body dyn accel
ciplot(roll_dot_dot_mean_wb_seq_mean_all_c4(isnan(roll_dot_dot_mean_wb_seq_mean_all_c4)==0)-1.96*roll_dot_dot_mean_wb_seq_ste_all_c4(isnan(roll_dot_dot_mean_wb_seq_mean_all_c4)==0),...
    roll_dot_dot_mean_wb_seq_mean_all_c4(isnan(roll_dot_dot_mean_wb_seq_mean_all_c4)==0)+1.96*roll_dot_dot_mean_wb_seq_ste_all_c4(isnan(roll_dot_dot_mean_wb_seq_mean_all_c4)==0),...
    t_wb_seq_mean_all(isnan(roll_dot_dot_mean_wb_seq_mean_all_c4)==0))
plot(t_wb_seq_mean_all,roll_dot_dot_mean_wb_seq_mean_all_c4,'-k.')

% xlabel('time [sec]')
ylabel('Aroll')
xlim([-.04 .04]) 
ylim([-3 3])
set(gca,'XTick',[-.04:.02:.04])
set(gca,'YTick',[-3:3:3])

%% pitch accel
subplot(3,1,2)
hold on

% body dyn accel
ciplot(pitch_dot_dot_mean_wb_seq_mean_all_c4(isnan(roll_dot_dot_mean_wb_seq_mean_all_c4)==0)-1.96*pitch_dot_dot_mean_wb_seq_ste_all_c4(isnan(roll_dot_dot_mean_wb_seq_mean_all_c4)==0),...
    pitch_dot_dot_mean_wb_seq_mean_all_c4(isnan(roll_dot_dot_mean_wb_seq_mean_all_c4)==0)+1.96*pitch_dot_dot_mean_wb_seq_ste_all_c4(isnan(roll_dot_dot_mean_wb_seq_mean_all_c4)==0),...
    t_wb_seq_mean_all(isnan(roll_dot_dot_mean_wb_seq_mean_all_c4)==0))
plot(t_wb_seq_mean_all,pitch_dot_dot_mean_wb_seq_mean_all_c4,'-k.')

% xlabel('time [sec]')
ylabel('Apitch')
xlim([-.04 .04]) 
ylim([-3 3])
set(gca,'XTick',[-.04:.02:.04])
set(gca,'YTick',[-3:3:3])

%% yaw accel
subplot(3,1,3)
hold on

% body dyn accel
ciplot(yaw_dot_dot_mean_wb_seq_mean_all_c4(isnan(roll_dot_dot_mean_wb_seq_mean_all_c4)==0)-1.96*yaw_dot_dot_mean_wb_seq_ste_all_c4(isnan(roll_dot_dot_mean_wb_seq_mean_all_c4)==0),...
    yaw_dot_dot_mean_wb_seq_mean_all_c4(isnan(roll_dot_dot_mean_wb_seq_mean_all_c4)==0)+1.96*yaw_dot_dot_mean_wb_seq_ste_all_c4(isnan(roll_dot_dot_mean_wb_seq_mean_all_c4)==0),...
    t_wb_seq_mean_all(isnan(roll_dot_dot_mean_wb_seq_mean_all_c4)==0))
plot(t_wb_seq_mean_all,yaw_dot_dot_mean_wb_seq_mean_all_c4,'-k.')

xlabel('time [sec]')
ylabel('Ayaw')
xlim([-.04 .04]) 
ylim([-3 3])
set(gca,'XTick',[-.04:.02:.04])
set(gca,'YTick',[-3:3:3])

saveas(gca,['MSfig_timeseries_accels_c4.fig'])
saveas(gca,['MSfig_timeseries_accels_c4.png'])
plot2svg(['MSfig_timeseries_accels_c4.svg'])



cd ..


















