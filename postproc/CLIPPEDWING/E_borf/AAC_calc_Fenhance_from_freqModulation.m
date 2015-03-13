clear
clc
close all

dir_now = cd;

%% wingbeat frequency of clipped wing flies from CLIPPED_GOOD wingbeat DB
dir_flydata = 'E:\Florian\Dropbox\Projects\flytracker\CLIPPED\WORK\CLIPPED_GOOD';
flydata = 'WBdataset_all_1252WBs_S2nS3.mat';

cd(dir_flydata)
load(flydata)
cd(dir_now)

steady_wbs = find(steady_nr_mean_wb == 1);

f_wb = mean([f_wb_L f_wb_R]')';
freq_clip = f_wb(steady_wbs);
freq_clip_mean = mean(freq_clip);
freq_clip_std = std(freq_clip);
freq_clip_n = length(freq_clip);

freq_clip_mean = 195

S2ratio_clip = SecondMomentRatio(steady_wbs);
S3ratio_clip = ThirdMomentRatio(steady_wbs);

%% freq MOD - Fz enhancement fit data from Muijres et al 2014
freq_steady = 188.7;
MODfreqFnorm = 41.5;

Fnorm_clip_steady = (freq_clip_mean - freq_steady)/MODfreqFnorm + 1;

%% plot results
mkdir('figures_cutWing_robofly')
cd('figures_cutWing_robofly')

figure
subplot(2,1,1)
plot(S2ratio_clip,freq_clip,'ok','markerfacecolor','b')
hold on
plot([.4 1.2],[freq_steady freq_steady],'-k')
xlabel('S2 ratio')
ylabel('freq')
title('wingbeat frequency of clipped flies')
axis([0.4 1.2 150 250])

subplot(2,1,2)
plot(S3ratio_clip,freq_clip,'ok','markerfacecolor','r')
hold on
plot([.4 1.2],[freq_steady freq_steady],'-k')
xlabel('S3 ratio')
ylabel('freq')
title('wingbeat frequency of clipped flies')
axis([0.4 1.2 150 250])

saveas(gcf,['FnMvsAmpStrokeRatio_MyAtCoM_robofly_NONcutWing_INCcali_REDUCE_LinFit.fig'])
saveas(gcf,['FnMvsAmpStrokeRatio_MyAtCoM_robofly_NONcutWing_INCcali_REDUCE_LinFit.png'])
plot2svg(['FnMvsAmpStrokeRatio_MyAtCoM_robofly_NONcutWing_INCcali_REDUCE_LinFit.svg'])

cd ..

%% save data
save('Fenhance_by_freqMod_4SteadyWB.mat',...
    'freq_clip','freq_clip_mean','freq_clip_std','freq_clip_n',...
    'S2ratio_clip','S3ratio_clip',...
    'freq_steady','MODfreqFnorm','Fnorm_clip_steady')
