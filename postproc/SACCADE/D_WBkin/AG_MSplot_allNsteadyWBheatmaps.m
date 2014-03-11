clc
clear
close all

name = 'WBdataset_all_steadyNmods_accelNorm.mat'
load(name)


plot_WBfunc_heatmap_steady_NOfourier
plot_WBfunc_heatmap_ALL
axis off

    saveas(gca,['MSplot_WBsteadyNall_heatmap.fig'])
    saveas(gca,['MSplot_WBsteadyNall_heatmap.png'])
    plot2svg(['MSplot_WBsteadyNall_heatmap.svg'])
