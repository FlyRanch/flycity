clc
clear
close all

% plot_qsFnMmean_transONLY_NOsumNtrend
plot_qsFnMmean_transNrot_NOsumNtrend
plot_mean_Borf_data_NOsumNtrend_INCyaw_Mnorm

mkdir('MSfigs')
cd('MSfigs')
% saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch_borfQSmodel.fig')
% saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch_borfQSmodel.png')
% plot2svg(['MSfig_Fenhance_Mroll_Myaw_Mpitch_borfQSmodel.svg'])

saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch_borfQSmodel_INCrot.fig')
saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch_borfQSmodel_INCrot.png')
plot2svg(['MSfig_Fenhance_Mroll_Myaw_Mpitch_borfQSmodel_INCrot.svg'])
cd ..
