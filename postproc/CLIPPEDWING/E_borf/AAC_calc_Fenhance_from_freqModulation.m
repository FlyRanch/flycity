clear
clc
close all

dir_now = cd;

%% wingbeat frequency of clipped wing flies from CLIPPED_GOOD wingbeat DB
steady_db = 'WBdataset_steady_1603WBs.mat';
clippedfly_db = 'WBdataset_all_1252WBs_S2nS3.mat';

cd ..
cd('CLIPPED_GOOD')
dir_flydata = cd;
load(clippedfly_db)
cd(dir_now)

steady_wbs = find(steady_nr_mean_wb == 1);

f_wb = mean([f_wb_L f_wb_R]')';
freq_clip = f_wb(steady_wbs);
freq_clip_mean = mean(freq_clip);
freq_clip_std = std(freq_clip);
freq_clip_n = length(freq_clip);
freq_clip_ste = freq_clip_std/sqrt(freq_clip_n);

% freq_clip_mean = 195

S2ratio_clip = SecondMomentRatio(steady_wbs);
S3ratio_clip = ThirdMomentRatio(steady_wbs);

%% freq MOD - Fz enhancement fit data from Muijres et al 2014
load(steady_db)
freq_steady_mean = f_wb_steady_meanCIstd(1);
freq_steady_ste = f_wb_steady_meanCIstd(2);

MODfreqFnorm = 41.5;
Fnorm_clip_steady = (freq_clip_mean - freq_steady_mean)/MODfreqFnorm + 1;

%% plot results
mkdir('figures_cutWing_robofly')
cd('figures_cutWing_robofly')

figure
subplot(2,1,1)
plot(S2ratio_clip,freq_clip,'ok','markerfacecolor','b')
hold on

plot([.4 1.2],[freq_clip_mean freq_clip_mean],'-b')
ciplot([freq_clip_mean-freq_clip_ste freq_clip_mean-freq_clip_ste],[freq_clip_mean+freq_clip_ste freq_clip_mean+freq_clip_ste],[.4 1.2],'-b')

plot([.4 1.2],[freq_steady_mean freq_steady_mean],'-k')
ciplot([freq_steady_mean-freq_steady_ste freq_steady_mean-freq_steady_ste],[freq_steady_mean+freq_steady_ste freq_steady_mean+freq_steady_ste],[.4 1.2],'-k')

% alpha(.5)
xlabel('S2 ratio')
ylabel('freq')
title('wingbeat frequency of clipped flies')
axis([0.4 1.2 150 250])

subplot(2,1,2)
plot(S3ratio_clip,freq_clip,'ok','markerfacecolor','r')
hold on

plot([.4 1.2],[freq_clip_mean freq_clip_mean],'-r')
ciplot([freq_clip_mean-freq_clip_ste freq_clip_mean-freq_clip_ste],[freq_clip_mean+freq_clip_ste freq_clip_mean+freq_clip_ste],[.4 1.2],'-r')

plot([.4 1.2],[freq_steady_mean freq_steady_mean],'-k')
ciplot([freq_steady_mean-freq_steady_ste freq_steady_mean-freq_steady_ste],[freq_steady_mean+freq_steady_ste freq_steady_mean+freq_steady_ste],[.4 1.2],'-k')

% alpha(.5)
xlabel('S3 ratio')
ylabel('freq')
title('wingbeat frequency of clipped flies')
axis([0.4 1.2 150 250])

saveas(gcf,['wingbeat_frequency_of_clipped_flies_nSteadyMeanSte_vs_S2S3.fig'])
saveas(gcf,['wingbeat_frequency_of_clipped_flies_nSteadyMeanSte_vs_S2S3.png'])
% saveas(gcf,['wingbeat_frequency_of_clipped_flies_nSteadyMeanSte_vs_S2S3.svg'])
plot2svg(['wingbeat_frequency_of_clipped_flies_nSteadyMeanSte_vs_S2S3.svg'])

cd ..

%% save data
save('Fenhance_by_freqMod_4SteadyWB.mat',...
    'freq_clip','freq_clip_mean','freq_clip_std','freq_clip_n','freq_clip_ste',...
    'S2ratio_clip','S3ratio_clip',...
    'freq_steady_mean','freq_steady_ste','MODfreqFnorm','Fnorm_clip_steady')
