clc
clear
close all

addpath('/home/florian/Dropbox/WORK/toolbox')
addpath('/home/florian/Dropbox/WORK/toolbox/plot2svg')
addpath('/home/florian/Dropbox/WORK/toolbox/CircStat')
addpath('/home/florian/Dropbox/WORK/toolbox/flytracker')

loadname='WBdataset_all_3845WBs.mat'
steady_name = 'WBdataset_steady_1807WBs.mat'

%% constants

%% settings
% limits (*std)
limit_steady = .5

% fourier orders
    stroke_fourier_order = 4;
    pitch_fourier_order = 8;
    dev_fourier_order = 8;
    
n=200; % bins

linewidth_timelines = .5;
linewidth_meanWBs = 1;
skip = 50;

% heatmap resolution
nx = 1000;
ny = 100;

dn=20   % datapoints in bin
dm=20   % bin shift

color_code_now = [.5 .5 .5];
color_mid = [.25 .25 .25];
cmap = abs(colormap(gray)-1);
close all

%% calc & plot steady wb
    calc_wbNORM_wingatt_steady
    
    plot_WBfunc_heatmap_steady_fourier
    saveas(gca,['WBnorm_wingkin_',num2str(n_steady),'steadyWBs_limit',num2str(limit_steady),'xSTD_heatmap.fig'])
    saveas(gca,['WBnorm_wingkin_',num2str(n_steady),'steadyWBs_limit',num2str(limit_steady),'xSTD_heatmap.png'])
    plot2svg(['WBnorm_wingkin_',num2str(n_steady),'steadyWBs_limit',num2str(limit_steady),'xSTD_heatmap.svg'])
    
    plot_WBfreqNpitchV_hist_steady_LnR
    saveas(gca,['WBfreqNpitchV_wingkin_',num2str(n_steady),'steadyWBs_limit',num2str(limit_steady),'xSTD.fig'])
    saveas(gca,['WBfreqNpitchV_wingkin_',num2str(n_steady),'steadyWBs_limit',num2str(limit_steady),'xSTD.png'])
    plot2svg(['WBfreqNpitchV_wingkin_',num2str(n_steady),'steadyWBs_limit',num2str(limit_steady),'xSTD.svg'])
    cd ..
