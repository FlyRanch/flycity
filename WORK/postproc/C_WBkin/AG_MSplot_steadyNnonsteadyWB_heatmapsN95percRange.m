clc
clear
close all

name = 'WBdataset_steadyNnonsteady.mat'
load(name)

mkdir('MSfigs_WBkin')
cd('MSfigs_WBkin')

% plot & save heatmap
plot_WBfunc_heatmap_steadyNmod

% plot & save 95% range
figure

subplot(3,3,1)
t_bin = t_wb_NONsteady_bins(:,1);
var = [stroke_wb_L_NONsteady_bins stroke_wb_R_NONsteady_bins];
plotcolor = 'c';
plot_WBbin_timelines_95percRange
% plot_WBbin_timelines_quartiles

t_bin = t_wb_steady_bins(:,1);
var = [stroke_wb_L_steady_bins stroke_wb_R_steady_bins];
plotcolor = 'r';
plot_WBbin_timelines_95percRange
% plot_WBbin_timelines_quartiles
% alpha(.25)
ylim([-90 90])

subplot(3,3,2)
t_bin = t_wb_NONsteady_bins(:,1);
var = [dev_wb_L_NONsteady_bins dev_wb_R_NONsteady_bins];
plotcolor = 'c';
plot_WBbin_timelines_95percRange
% plot_WBbin_timelines_quartiles

t_bin = t_wb_steady_bins(:,1);
var = [dev_wb_L_steady_bins dev_wb_R_steady_bins];
plotcolor = 'r';
plot_WBbin_timelines_95percRange
% plot_WBbin_timelines_quartiles
% alpha(.25)
ylim([-30 30])

subplot(3,3,3)
t_bin = t_wb_NONsteady_bins(:,1);
var = [pitch_wb_L_NONsteady_bins pitch_wb_R_NONsteady_bins];
plotcolor = 'c';
plot_WBbin_timelines_95percRange
% plot_WBbin_timelines_quartiles

t_bin = t_wb_steady_bins(:,1);
var = [pitch_wb_L_steady_bins pitch_wb_R_steady_bins];
plotcolor = 'r';
plot_WBbin_timelines_95percRange
% plot_WBbin_timelines_quartiles
% alpha(.25)
ylim([0 180])

saveas(gca,['MSplot_WBsteadyNmod_95percRange.fig'])
saveas(gca,['MSplot_WBsteadyNmod_95percRange.png'])
plot2svg(['MSplot_WBsteadyNmod_95percRange.svg'])    

% plot & save quartiles range
figure

subplot(3,3,1)
t_bin = t_wb_NONsteady_bins(:,1);
var = [stroke_wb_L_NONsteady_bins stroke_wb_R_NONsteady_bins];
plotcolor = 'c';
% plot_WBbin_timelines_95percRange
plot_WBbin_timelines_quartiles

t_bin = t_wb_steady_bins(:,1);
var = [stroke_wb_L_steady_bins stroke_wb_R_steady_bins];
plotcolor = 'r';
% plot_WBbin_timelines_95percRange
plot_WBbin_timelines_quartiles
% alpha(.25)
ylim([-90 90])

subplot(3,3,2)
t_bin = t_wb_NONsteady_bins(:,1);
var = [dev_wb_L_NONsteady_bins dev_wb_R_NONsteady_bins];
plotcolor = 'c';
% plot_WBbin_timelines_95percRange
plot_WBbin_timelines_quartiles

t_bin = t_wb_steady_bins(:,1);
var = [dev_wb_L_steady_bins dev_wb_R_steady_bins];
plotcolor = 'r';
% plot_WBbin_timelines_95percRange
plot_WBbin_timelines_quartiles
% alpha(.25)
ylim([-30 30])

subplot(3,3,3)
t_bin = t_wb_NONsteady_bins(:,1);
var = [pitch_wb_L_NONsteady_bins pitch_wb_R_NONsteady_bins];
plotcolor = 'c';
% plot_WBbin_timelines_95percRange
plot_WBbin_timelines_quartiles

t_bin = t_wb_steady_bins(:,1);
var = [pitch_wb_L_steady_bins pitch_wb_R_steady_bins];
plotcolor = 'r';
% plot_WBbin_timelines_95percRange
plot_WBbin_timelines_quartiles
% alpha(.25)
ylim([0 180])

saveas(gca,['MSplot_WBsteadyNmod_quartiles.fig'])
saveas(gca,['MSplot_WBsteadyNmod_quartiles.png'])
plot2svg(['MSplot_WBsteadyNmod_quartiles.svg'])    

cd ..