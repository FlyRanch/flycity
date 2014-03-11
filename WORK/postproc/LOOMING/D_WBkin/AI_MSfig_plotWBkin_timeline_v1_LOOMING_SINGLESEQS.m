% kinematics for seqs in fig 1d
% clear;
% clc;
% close all
% warning off
% 
% addpath('/home/florian/Dropbox/WORK/flytracker/flytracker');
% addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/mex/');
% addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/core/');
% addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/results/');
% 
% loadname=dir('WBdataset_temporal_dynamics_*')
% loadname = loadname.name;
% load(loadname)
% % 
% % load('FnMqs_timeseries.mat')
% 
mkdir('MSfigs_WBnBodyKin_TempDynamics')
cd('MSfigs_WBnBodyKin_TempDynamics')

for i = [11 3 9 5 1] % resp 1205#3 1120.4 1205.1 1129.3 1120.5 1120.1

close all

freq = f_wb_seq_all(i,:)';

stroke_L = stroke_wb_L_seq_bins_all(:,i);
stroke_R = stroke_wb_R_seq_bins_all(:,i);
dev_L = dev_wb_L_seq_bins_all(:,i);
dev_R = dev_wb_R_seq_bins_all(:,i);
pitch_L = pitch_wb_L_seq_bins_all(:,i);
pitch_R = pitch_wb_R_seq_bins_all(:,i);

stroke_L_pre(:,:) = stroke_wb_L_seq_bins_pre(:,i,:);
stroke_R_pre(:,:) = stroke_wb_R_seq_bins_pre(:,i,:);
stroke_L_post(:,:) = stroke_wb_L_seq_bins_post(:,i,:);
stroke_R_post(:,:) = stroke_wb_R_seq_bins_post(:,i,:);

stroke_L_bins = [stroke_L_pre stroke_L_post];
stroke_R_bins = [stroke_R_pre stroke_R_post];

dev_L_pre(:,:) = dev_wb_L_seq_bins_pre(:,i,:);
dev_R_pre(:,:) = dev_wb_R_seq_bins_pre(:,i,:);
dev_L_post(:,:) = dev_wb_L_seq_bins_post(:,i,:);
dev_R_post(:,:) = dev_wb_R_seq_bins_post(:,i,:);

dev_L_bins = [dev_L_pre dev_L_post];
dev_R_bins = [dev_R_pre dev_R_post];

pitch_L_pre(:,:) = pitch_wb_L_seq_bins_pre(:,i,:);
pitch_R_pre(:,:) = pitch_wb_R_seq_bins_pre(:,i,:);
pitch_L_post(:,:) = pitch_wb_L_seq_bins_post(:,i,:);
pitch_R_post(:,:) = pitch_wb_R_seq_bins_post(:,i,:);

pitch_L_bins = [pitch_L_pre pitch_L_post];
pitch_R_bins = [pitch_R_pre pitch_R_post];

%% wb params

%% stroke
% L
stroke_L_bins_MAX1 = max(stroke_L_bins(1:100,:))';
stroke_L_bins_MAX2 = max(stroke_L_bins(100:end,:))';
stroke_L_bins_MAX1 = [stroke_L_bins_MAX1;nan];
stroke_L_bins_MAX2 = [nan;stroke_L_bins_MAX2];
stroke_L_bins_MAX = [stroke_L_bins_MAX1 stroke_L_bins_MAX2];
stroke_L_bins_MAX = max(stroke_L_bins_MAX(1:end-1,:)')';

stroke_L_bins_MIN = min(stroke_L_bins)';
Astroke_L_bins = stroke_L_bins_MAX - stroke_L_bins_MIN;

% R
stroke_R_bins_MAX1 = max(stroke_R_bins(1:100,:))';
stroke_R_bins_MAX2 = max(stroke_R_bins(100:end,:))';
stroke_R_bins_MAX1 = [stroke_R_bins_MAX1;nan];
stroke_R_bins_MAX2 = [nan;stroke_R_bins_MAX2];
stroke_R_bins_MAX = [stroke_R_bins_MAX1 stroke_R_bins_MAX2];
stroke_R_bins_MAX = max(stroke_R_bins_MAX(1:end-1,:)')';

stroke_R_bins_MIN = min(stroke_R_bins)';
Astroke_R_bins = stroke_R_bins_MAX - stroke_R_bins_MIN;

% R+L
Astroke_bins = nanmean([Astroke_R_bins Astroke_L_bins]')';
stroke_bins_MAX = nanmean([stroke_R_bins_MAX stroke_L_bins_MAX]')';
stroke_bins_MIN = nanmean([stroke_R_bins_MIN stroke_L_bins_MIN]')';

% R-L
dAstroke_bins = Astroke_R_bins - Astroke_L_bins;
Dstroke_bins_MAX = stroke_R_bins_MAX - stroke_L_bins_MAX;
Dstroke_bins_MIN = stroke_R_bins_MIN - stroke_L_bins_MIN;

%% wing rotation
% L
pitch_L_bins_MAXmidDS = max(pitch_L_bins(40:100,:))'-90;
pitch_L_bins_MINmidUS = min(pitch_L_bins(150:end,:))'-90;
Apitch_L_bins = pitch_L_bins_MAXmidDS - pitch_L_bins_MINmidUS;

% R
pitch_R_bins_MAXmidDS = max(pitch_R_bins(40:100,:))'-90;
pitch_R_bins_MINmidUS = min(pitch_R_bins(150:end,:))'-90;
Apitch_R_bins = pitch_R_bins_MAXmidDS - pitch_R_bins_MINmidUS;

% R+L
Apitch_bins = nanmean([Apitch_R_bins Apitch_L_bins]')';
pitch_bins_MAXmidDS = nanmean([pitch_R_bins_MAXmidDS pitch_L_bins_MAXmidDS]')';
pitch_bins_MINmidUS = nanmean([pitch_R_bins_MINmidUS pitch_L_bins_MINmidUS]')';

% R-L
dApitch_bins = Apitch_R_bins - Apitch_L_bins;
Dpitch_bins_MAXmidDS = pitch_R_bins_MAXmidDS - pitch_L_bins_MAXmidDS;
Dpitch_bins_MINmidUS = pitch_R_bins_MINmidUS - pitch_L_bins_MINmidUS;

%% wing deviation
% L
dev_L_bins_MAXds = max(dev_L_bins(1:50,:))'-90;
dev_L_bins_MINds = min(dev_L_bins(1:100,:))'-90;
dev_L_bins_MAXus = max(dev_L_bins(50:150,:))'-90;
dev_L_bins_MINus = min(dev_L_bins(100:end,:))'-90;

Adev_L_bins_DS = dev_L_bins_MAXds - dev_L_bins_MINds;
Adev_L_bins_US = dev_L_bins_MAXus - dev_L_bins_MINus;

% R
dev_R_bins_MAXds = max(dev_R_bins(1:50,:))'-90;
dev_R_bins_MINds = min(dev_R_bins(1:100,:))'-90;
dev_R_bins_MAXus = max(dev_R_bins(50:150,:))'-90;
dev_R_bins_MINus = min(dev_R_bins(100:end,:))'-90;

Adev_R_bins_DS = dev_R_bins_MAXds - dev_R_bins_MINds;
Adev_R_bins_US = dev_R_bins_MAXus - dev_R_bins_MINus;

% R+L
Adev_bins_DS = nanmean([Adev_R_bins_DS Adev_L_bins_DS]')';
Adev_bins_US = nanmean([Adev_R_bins_US Adev_L_bins_US]')';

dev_bins_MAXds = nanmean([dev_R_bins_MAXds dev_L_bins_MAXds]')';
dev_bins_MINds = nanmean([dev_R_bins_MINds dev_L_bins_MINds]')';
dev_bins_MAXus = nanmean([dev_R_bins_MAXus dev_L_bins_MAXus]')';
dev_bins_MINus = nanmean([dev_R_bins_MINus dev_L_bins_MINus]')';

% R-L
dAdev_bins_DS = Adev_R_bins_DS - Adev_L_bins_DS;
dAdev_bins_US = Adev_R_bins_US - Adev_L_bins_US;

Ddev_bins_MAXds = dev_R_bins_MAXds - dev_L_bins_MAXds;
Ddev_bins_MINds = dev_R_bins_MINds - dev_L_bins_MINds;
Ddev_bins_MAXus = dev_R_bins_MAXus - dev_L_bins_MAXus;
Ddev_bins_MINus = dev_R_bins_MINus - dev_L_bins_MINus;

%% plot kin

% freq
figure
subplot(3,1,1)
hold on
plot(t_wb_seq_mean_all,freq,'-k.','linewidth',1)

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('freq [1/sec]')
ylim([180 220])
set(gca,'XTick',[-.04:.02:.04])
set(gca,'YTick',[0:20:220])

saveas(gca,['WBvsTimeSeqs_freq_seq',num2str(i),'.fig'])
saveas(gca,['WBvsTimeSeqs_freq_seq',num2str(i),'.png'])
plot2svg(['WBvsTimeSeqs_freq_seq',num2str(i),'.svg'])

% angles
figure
subplot(3,1,1)
hold on
plot(t_wb_seq_bins_mean_all,stroke_L,'linewidth',.25,'color',[0 .5 1])
plot(t_wb_seq_bins_mean_all,stroke_R,'linewidth',.25,'color',[1 0 0])

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('stroke angle [deg]')
ylim([-90 90])
set(gca,'YTick',[-90:45:90])

subplot(3,1,2)
hold on
plot(t_wb_seq_bins_mean_all,dev_L,'linewidth',.25,'color',[0 .5 1])
plot(t_wb_seq_bins_mean_all,dev_R,'linewidth',.25,'color',[1 0 0])

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('deviation angle [deg]')
ylim([-30 30])
set(gca,'YTick',[-90:30:90])

subplot(3,1,3)
hold on
plot(t_wb_seq_bins_mean_all,pitch_L,'linewidth',.25,'color',[0 .5 1])
plot(t_wb_seq_bins_mean_all,pitch_R,'linewidth',.25,'color',[1 0 0])

xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('rotation angle [deg]')
ylim([0 180])
set(gca,'YTick',[-90:45:180])

saveas(gca,['WBvsTimeSeqs_WBangles_seq',num2str(i),'.fig'])
saveas(gca,['WBvsTimeSeqs_WBangles_seq',num2str(i),'.png'])
plot2svg(['WBvsTimeSeqs_WBangles_seq',num2str(i),'.svg'])


%% WB parameters

%% stroke
figure

% L+R
subplot(3,1,1)
hold on
plot(t_wb_seq_mean_all,Astroke_R_bins,'-k.','color',[0 .5 1])
plot(t_wb_seq_mean_all,Astroke_L_bins,'-k.','color',[1 0 0])

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('stroke angle amplitude [deg]')
ylim([120 180])
set(gca,'YTick',[-90:30:190])

% L-R
subplot(3,1,2)
hold on
plot(t_wb_seq_mean_all,Dstroke_bins_MAX,'-k.')
plot(t_wb_seq_mean_all,Dstroke_bins_MIN,'-k.','color',[.5 .5 .5])
% plot(t_wb_seq_mean_all,-dAstroke_bins,'-g.')

xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('stroke angle L-R [deg]')
ylim([-7.5 7.5])
set(gca,'YTick',[-90:7.5:90])

saveas(gca,['MSfig_timeseries_strokeParams_seq',num2str(i),'.fig'])
saveas(gca,['MSfig_timeseries_strokeParams_seq',num2str(i),'.png'])
plot2svg(['MSfig_timeseries_strokeParams_seq',num2str(i),'.svg'])

%% pitch
figure

% Downstroke
subplot(3,1,1)
hold on
plot(t_wb_seq_mean_all,pitch_R_bins_MAXmidDS,'-c.','color',[0 .5 1])
plot(t_wb_seq_mean_all,pitch_L_bins_MAXmidDS,'-m.','color',[1 0 0])

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('rotation angle downstroke [deg]')
ylim([45 75])
set(gca,'YTick',[-90:15:90])

% Upstroke
subplot(3,1,2)
hold on
plot(t_wb_seq_mean_all,pitch_R_bins_MINmidUS,'-c.','color',[0 .5 1])
plot(t_wb_seq_mean_all,pitch_L_bins_MINmidUS,'-m.','color',[1 0 0])

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('rotation angle upstroke [deg]')
ylim([-60 -30])
set(gca,'YTick',[-90:15:90])

% L-R
subplot(3,1,3)
hold on
plot(t_wb_seq_mean_all,Dpitch_bins_MAXmidDS,'-k.')
plot(t_wb_seq_mean_all,Dpitch_bins_MINmidUS,'-k.','color',[.5 .5 .5])
% plot(t_wb_seq_mean_all,-dApitch_bins,'-g.')

xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('rotation angle L-R [deg]')
ylim([-7.5 7.5])
set(gca,'YTick',[-90:7.5:90])

saveas(gca,['MSfig_timeseries_rotParams_seq',num2str(i),'.fig'])
saveas(gca,['MSfig_timeseries_rotParams_seq',num2str(i),'.png'])
plot2svg(['MSfig_timeseries_rotParams_seq',num2str(i),'.svg'])

%% dev
figure

% Downstroke
subplot(3,1,1)
hold on
plot(t_wb_seq_mean_all,Adev_R_bins_DS,'-c.','color',[0 .5 1])
plot(t_wb_seq_mean_all,Adev_L_bins_DS,'-m.','color',[1 0 0])

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('deviation amplitude downstroke [deg]')
ylim([-15 45])
set(gca,'YTick',[-90:15:90])

% Upstroke
subplot(3,1,2)
hold on
plot(t_wb_seq_mean_all,Adev_R_bins_US,'-c.','color',[0 .5 1])
plot(t_wb_seq_mean_all,Adev_L_bins_US,'-m.','color',[1 0 0])

% xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('deviation amplitude upstroke [deg]')
ylim([-15 45])
set(gca,'YTick',[-90:15:90])

% L-R
subplot(3,1,3)
hold on
plot(t_wb_seq_mean_all,dAdev_bins_DS,'-k.')
plot(t_wb_seq_mean_all,dAdev_bins_US,'-k.','color',[.5 .5 .5])
% plot(t_wb_seq_mean_all,-dAdev_bins,'-g.')

xlabel('time [sec]')
xlim([-.04 .04]) 
ylabel('deviation amplitude L-R [deg]')
ylim([-30 30])
set(gca,'YTick',[-90:15:90])

saveas(gca,['MSfig_timeseries_devParams_seq',num2str(i),'.fig'])
saveas(gca,['MSfig_timeseries_devParams_seq',num2str(i),'.png'])
plot2svg(['MSfig_timeseries_devParams_seq',num2str(i),'.svg'])

end

cd ..