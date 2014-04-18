% all MODs
figure
hold on
plot(mod_value_allNOfreq,Mx_mean_allNOfreq,'.b')
plot(mod_value_allNOfreq,My_mean_allNOfreq,'.g')
plot(mod_value_allNOfreq,Mz_mean_allNOfreq,'.r')
legend('Mroll','Mpitch','Myaw','location','nw')
grid on
xlabel('MOD value')
ylabel('M')
saveas(gca, [save_name,'_MODvalue_allNOfreqMOD_LnRsym.png'])

% % freq MODs
% figure
% hold on
% plot(freq_freq,Mx_mean_freq,'.b')
% plot(freq_freq,My_mean_freq,'.g')
% plot(freq_freq,Mz_mean_freq,'.r')
% legend('Mroll','Mpitch','Myaw','location','nw') 
% grid on
% xlabel('freq')
% ylabel('M')
% saveas(gca, [save_name,'_freq_freqMOD_LnRsym.png'])
% 
% figure
% hold on
% plot(mod_value_freq,Mx_mean_freq,'.b')
% plot(mod_value_freq,My_mean_freq,'.g')
% plot(mod_value_freq,Mz_mean_freq,'.r')
% legend('Mroll','Mpitch','Myaw','location','nw') 
% grid on
% xlabel('MOD value')
% ylabel('M')
% saveas(gca, [save_name,'_MODvalue_freqMOD_LnRsym.png'])

% % stroke MODs
% figure
% hold on
% plot(Astroke_stroke,Mx_mean_stroke,'.b')
% plot(Astroke_stroke,My_mean_stroke,'.g')
% plot(Astroke_stroke,Mz_mean_stroke,'.r')
% legend('Mroll','Mpitch','Myaw','location','nw') 
% grid on
% xlabel('Astroke')
% ylabel('M')
% saveas(gca, [save_name,'_Astroke_strokeMOD_LnRsym.png'])

figure
hold on
plot(mod_value_stroke,Mx_mean_stroke,'.b')
plot(mod_value_stroke,My_mean_stroke,'.g')
plot(mod_value_stroke,Mz_mean_stroke,'.r')
legend('Mroll','Mpitch','Myaw','location','nw') 
grid on
xlabel('MOD value')
ylabel('M')
saveas(gca, [save_name,'_MODvalue_strokeMOD_LnRsym.png'])

% % pitch MODs
% figure
% hold on
% plot(Apitch_pitch,Mx_mean_pitch,'.b')
% plot(Apitch_pitch,My_mean_pitch,'.g')
% plot(Apitch_pitch,Mz_mean_pitch,'.r')
% legend('Mroll','Mpitch','Myaw','location','nw') 
% grid on
% xlabel('Apitch')
% ylabel('M')
% saveas(gca, [save_name,'_Apitch_pitchMOD_LnRsym.png'])

figure
hold on
plot(mod_value_pitch,Mx_mean_pitch,'.b')
plot(mod_value_pitch,My_mean_pitch,'.g')
plot(mod_value_pitch,Mz_mean_pitch,'.r')
legend('Mroll','Mpitch','Myaw','location','nw') 
grid on
xlabel('MOD value')
ylabel('M')
saveas(gca, [save_name,'_MODvalue_pitchMOD_LnRsym.png'])

% % dev MODs
% figure
% hold on
% plot(Adev_dev,Mx_mean_dev,'.b')
% plot(Adev_dev,My_mean_dev,'.g')
% plot(Adev_dev,Mz_mean_dev,'.r')
% legend('Mroll','Mpitch','Myaw','location','nw') 
% grid on
% xlabel('Adev')
% ylabel('M')
% saveas(gca, [save_name,'_Adev_devMOD_LnRsym.png'])

figure
hold on
plot(mod_value_dev,Mx_mean_dev,'.b')
plot(mod_value_dev,My_mean_dev,'.g')
plot(mod_value_dev,Mz_mean_dev,'.r')
legend('Mroll','Mpitch','Myaw','location','nw') 
grid on
xlabel('MOD value')
ylabel('M')
saveas(gca, [save_name,'_MODvalue_devMOD_LnRsym.png'])

% compare MODs
figure
hold on
plot(mod_value_allNOfreq,Mx_mean_allNOfreq,'.k')
% plot(mod_value_freq,Mx_mean_freq,'.r')
plot(mod_value_stroke,Mx_mean_stroke,'.g')
plot(mod_value_pitch,Mx_mean_pitch,'.b')
plot(mod_value_dev,Mx_mean_dev,'.c')
plot(mod_value_allNOfreq(1:MODmin),Mx_mean_sum,'ok')

legend('all','stroke','pitch','dev','sum','location','nw')
grid on
xlabel('MOD value')
ylabel('Mroll')
saveas(gca, [save_name,'_MODvalue_compare_LnRsym.png'])
