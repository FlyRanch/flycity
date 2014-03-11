clc
clear
close all

% plot_timelines_QSmodel_transONLY_rollREVERSE
plot_timelines_QSmodel_transNrot_rollREVERSE
% plot_timelines_borf_NOyaw
plot_timelines_borf_INCyaw_Mnorm
% plot_timelines_borf_INCyaw_Mnorm_butter


mkdir('MSfigs')
cd('MSfigs')

% % trans only
% figure(1)
% saveas(gca, 'MSfig_Fenhance_timelines_BORFnQSmodel_transONLY.fig')
% saveas(gca, 'MSfig_Fenhance_timelines_BORFnQSmodel_transONLY.png')
% plot2svg(['MSfig_Fenhance_timelines_BORFnQSmodel_transONLY.svg'])
% 
% figure(2)
% saveas(gca, 'MSfig_Mroll_timelines_BORFnQSmodel_transONLY.fig')
% saveas(gca, 'MSfig_Mroll_timelines_BORFnQSmodel_transONLY.png')
% plot2svg(['MSfig_Mroll_timelines_BORFnQSmodel_transONLY.svg'])
% 
% figure(3)
% saveas(gca, 'MSfig_Myaw_timelines_BORFnQSmodel_transONLY.fig')
% saveas(gca, 'MSfig_Myaw_timelines_BORFnQSmodel_transONLY.png')
% plot2svg(['MSfig_Myaw_timelines_BORFnQSmodel_transONLY.svg'])
% 
% figure(4)
% saveas(gca, 'MSfig_Mpitch_timelines_BORFnQSmodel_transONLY.fig')
% saveas(gca, 'MSfig_Mpitch_timelines_BORFnQSmodel_transONLY.png')
% plot2svg(['MSfig_Mpitch_timelines_BORFnQSmodel_transONLY.svg'])

% trans&rot
figure(1)
saveas(gca, 'MSfig_Fenhance_timelines_BORFnQSmodel_transNrot_NObutter.fig')
saveas(gca, 'MSfig_Fenhance_timelines_BORFnQSmodel_transNrot_NObutter.png')
plot2svg(['MSfig_Fenhance_timelines_BORFnQSmodel_transNrot_NObutter.svg'])

figure(2)
saveas(gca, 'MSfig_Mroll_timelines_BORFnQSmodel_transNrot_NObutter.fig')
saveas(gca, 'MSfig_Mroll_timelines_BORFnQSmodel_transNrot_NObutter.png')
plot2svg(['MSfig_Mroll_timelines_BORFnQSmodel_transNrot_NObutter.svg'])

figure(3)
saveas(gca, 'MSfig_Myaw_timelines_BORFnQSmodel_transNrot_NObutter.fig')
saveas(gca, 'MSfig_Myaw_timelines_BORFnQSmodel_transNrot_NObutter.png')
plot2svg(['MSfig_Myaw_timelines_BORFnQSmodel_transNrot_NObutter.svg'])

figure(4)
saveas(gca, 'MSfig_Mpitch_timelines_BORFnQSmodel_transNrot_NObutter.fig')
saveas(gca, 'MSfig_Mpitch_timelines_BORFnQSmodel_transNrot_NObutter.png')
plot2svg(['MSfig_Mpitch_timelines_BORFnQSmodel_transNrot_NObutter.svg'])

cd ..
