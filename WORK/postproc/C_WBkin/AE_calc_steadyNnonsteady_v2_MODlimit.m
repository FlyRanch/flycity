clc
clear
close all

addpath('/home/florian/Dropbox/WORK/toolbox')
addpath('/home/florian/Dropbox/WORK/toolbox/plot2svg')
addpath('/home/florian/Dropbox/WORK/toolbox/CircStat')
addpath('/home/florian/Dropbox/WORK/toolbox/flytracker')

loadname=dir('WBdataset_all_*')
loadname=loadname.name;

plot_on = 1;
plot_on = 0;

save_on = 1;
% save_on = 0;

mkdir('MSfigs_WBkin')
load(loadname)

%% settings
% limits (*std)
limit_steady = .5
limit_mod = 1
norm = 3

% fourier orders
    stroke_fourier_order = 4;
    pitch_fourier_order = 8;
    dev_fourier_order = 8;

    MOD_fourier_order = 8;

% number of polinomials
    n_pol_stroke = 12; % Order of used polynomials
    n_pol_pitch = 14; % Order of used polynomials
    n_pol_dev = 10; % Order of used polynomials 
    
    n_pol_MOD = 10;

plot_fits = 1;

n=200; % bins

linewidth_timelines = .5;
linewidth_meanWBs = 1;
skip = 50;

% heatmap resolution
nx = 1000;
ny = 100;

cmap_180 = jet(180);

% polyfit & 95% cof int settings
order = 3;
dn=20   % datapoints in bin
dm=20   % bin shift

color_code_now = [.5 .5 .5];
color_mid = [.25 .25 .25];
cmap = abs(colormap(gray)-1);
close all

%% calc & plot steadyNnonsteady wb
    
%     calc_wbNORM_wingatt_steadyNnonsteady
    calc_wbNORM_wingatt_steadyNmod
    
if plot_on ==1
    cd('MSfigs_WBkin')
    
    plot_WBfunc_heatmap_steadyNmod
%     plot_WBfunc_heatmap_steady_NOfourier
%     plot_WBfunc_heatmap_NONsteady_NOfourier
%     plot_WBfunc_heatmap_ALL
%     axis off
% 
%     saveas(gca,['MSplot_WBsteadyNmod_heatmap.fig'])
%     saveas(gca,['MSplot_WBsteadyNmod_heatmap.png'])
%     plot2svg(['MSplot_WBsteadyNmod_heatmap.svg'])    
    cd ..
end

