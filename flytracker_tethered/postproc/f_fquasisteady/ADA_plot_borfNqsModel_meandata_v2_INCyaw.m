clc
clear
close all

% plot_mean_Borf_data_NOsumNtrend_NOyaw
plot_qsFnMmean_transONLY_NOsumNtrend
plot_mean_Borf_data_NOsumNtrend_INCyaw_Mnorm

mkdir('MSfigs')
cd('MSfigs')
saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch_borfQSmodel.fig')
saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch_borfQSmodel.png')
plot2svg(['MSfig_Fenhance_Mroll_Myaw_Mpitch_borfQSmodel.svg'])
cd ..
