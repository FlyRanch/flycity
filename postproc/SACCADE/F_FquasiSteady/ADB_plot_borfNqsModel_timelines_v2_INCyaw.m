clc
clear
close all

plot_timelines_QSmodel_transONLY
% plot_timelines_borf_NOyaw
plot_timelines_borf_INCyaw_Mnorm


mkdir('MSfigs')
cd('MSfigs')

figure(1)
saveas(gca, 'MSfig_Fenhance_timelines_BORFnQSmodel_transONLY.fig')
saveas(gca, 'MSfig_Fenhance_timelines_BORFnQSmodel_transONLY.png')
plot2svg(['MSfig_Fenhance_timelines_BORFnQSmodel_transONLY.svg'])

figure(2)
saveas(gca, 'MSfig_Mroll_timelines_BORFnQSmodel_transONLY.fig')
saveas(gca, 'MSfig_Mroll_timelines_BORFnQSmodel_transONLY.png')
plot2svg(['MSfig_Mroll_timelines_BORFnQSmodel_transONLY.svg'])

figure(3)
saveas(gca, 'MSfig_Myaw_timelines_BORFnQSmodel_transONLY.fig')
saveas(gca, 'MSfig_Myaw_timelines_BORFnQSmodel_transONLY.png')
plot2svg(['MSfig_Myaw_timelines_BORFnQSmodel_transONLY.svg'])

figure(4)
saveas(gca, 'MSfig_Mpitch_timelines_BORFnQSmodel_transONLY.fig')
saveas(gca, 'MSfig_Mpitch_timelines_BORFnQSmodel_transONLY.png')
plot2svg(['MSfig_Mpitch_timelines_BORFnQSmodel_transONLY.svg'])

cd ..
