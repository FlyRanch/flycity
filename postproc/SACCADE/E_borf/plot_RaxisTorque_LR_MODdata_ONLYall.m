% all MODs
figure
hold on

% plot(mod_value_allNOfreq,Mx_mean_allNOfreq,'-b.')
% plot(mod_value_allNOfreq,My_mean_allNOfreq,'-g.')
% plot(mod_value_allNOfreq,Mz_mean_allNOfreq,'-r.')
% plot(mod_value_allNOfreq,M_R_mean_allNOfreq,'-c.')
% plot(mod_value_allNOfreq,M_L_mean_allNOfreq,'-m.')

plot(mod_value_allNOfreq,Mx_mean_allNOfreq-Mx_steady,'-b.')
plot(mod_value_allNOfreq,My_mean_allNOfreq-My_steady,'-g.')
plot(mod_value_allNOfreq,Mz_mean_allNOfreq-Mz_steady,'-r.')
plot(mod_value_allNOfreq,M_R_mean_allNOfreq-M_R_steady,'-c.')
plot(mod_value_allNOfreq,M_L_mean_allNOfreq-M_L_steady,'-m.')

legend('Mroll','Mpitch','Myaw','M_R','M_L','location','nw')
grid on
xlabel('MOD value')
ylabel('M')
saveas(gca, [save_name,'_MODvalue_allNOfreqMOD.png'])

% figure
% hold on
% plot(mod_value_freq,Mx_mean_freq,'-b.')
% plot(mod_value_freq,My_mean_freq,'-g.')
% plot(mod_value_freq,Mz_mean_freq,'-r.')
% legend('Mroll','Mpitch','Myaw','location','nw') 
% grid on
% xlabel('MOD value')
% ylabel('M')
% saveas(gca, [save_name,'_MODvalue_freqMOD.png'])

% figure
% hold on
% 
% % plot(mod_value_stroke,Mx_mean_stroke,'-b.')
% % plot(mod_value_stroke,My_mean_stroke,'-g.')
% % plot(mod_value_stroke,Mz_mean_stroke,'-r.')
% % plot(mod_value_stroke,M_R_mean_stroke,'-c.')
% % plot(mod_value_stroke,M_L_mean_stroke,'-m.')
% 
% plot(mod_value_stroke,Mx_mean_stroke-Mx_steady,'-b.')
% plot(mod_value_stroke,My_mean_stroke-My_steady,'-g.')
% plot(mod_value_stroke,Mz_mean_stroke-Mz_steady,'-r.')
% plot(mod_value_stroke,M_R_mean_stroke-M_R_steady,'-c.')
% plot(mod_value_stroke,M_L_mean_stroke-M_L_steady,'-m.')
% 
% legend('Mroll','Mpitch','Myaw','M_R','M_L','location','nw')
% grid on
% xlabel('MOD value')
% ylabel('M')
% saveas(gca, [save_name,'_MODvalue_strokeMOD.png'])
% 
% figure
% hold on
% 
% % plot(mod_value_pitch,Mx_mean_pitch,'-b.')
% % plot(mod_value_pitch,My_mean_pitch,'-g.')
% % plot(mod_value_pitch,Mz_mean_pitch,'-r.')
% % plot(mod_value_pitch,M_R_mean_pitch,'-c.')
% % plot(mod_value_pitch,M_L_mean_pitch,'-m.')
% 
% plot(mod_value_pitch,Mx_mean_pitch-Mx_steady,'-b.')
% plot(mod_value_pitch,My_mean_pitch-My_steady,'-g.')
% plot(mod_value_pitch,Mz_mean_pitch-Mz_steady,'-r.')
% plot(mod_value_pitch,M_R_mean_pitch-M_R_steady,'-c.')
% plot(mod_value_pitch,M_L_mean_pitch-M_L_steady,'-m.')
% 
% legend('Mroll','Mpitch','Myaw','M_R','M_L','location','nw')
% grid on
% xlabel('MOD value')
% ylabel('M')
% saveas(gca, [save_name,'_MODvalue_pitchMOD.png'])
% 
% figure
% hold on
% 
% % plot(mod_value_dev,Mx_mean_dev,'-b.')
% % plot(mod_value_dev,My_mean_dev,'-g.')
% % plot(mod_value_dev,Mz_mean_dev,'-r.')
% % plot(mod_value_dev,M_R_mean_dev,'-c.')
% % plot(mod_value_dev,M_L_mean_dev,'-m.')
% 
% plot(mod_value_dev,Mx_mean_dev-Mx_steady,'-b.')
% plot(mod_value_dev,My_mean_dev-My_steady,'-g.')
% plot(mod_value_dev,Mz_mean_dev-Mz_steady,'-r.')
% plot(mod_value_dev,M_R_mean_dev-M_R_steady,'-c.')
% plot(mod_value_dev,M_L_mean_dev-M_L_steady,'-m.')
% 
% legend('Mroll','Mpitch','Myaw','M_R','M_L','location','nw')
% grid on
% xlabel('MOD value')
% ylabel('M')
% saveas(gca, [save_name,'_MODvalue_devMOD.png'])
% 
% % compare MODs
% figure
% hold on
% plot(mod_value_allNOfreq,M_R_mean_allNOfreq-M_R_steady,'-k.')
% % plot(mod_value_freq,M_R_mean_freq-M_R_steady,'-r.')
% plot(mod_value_stroke,M_R_mean_stroke-M_R_steady,'-g.')
% plot(mod_value_pitch,M_R_mean_pitch-M_R_steady,'-b.')
% plot(mod_value_dev,M_R_mean_dev-M_R_steady,'-c.')
% plot(mod_value_allNOfreq(1:MODmin),M_R_mean_sum-M_R_mean_sum_steady,'-ok')
% 
% legend('all','stroke','pitch','dev','sum','location','nw')
% grid on
% xlabel('MOD value')
% ylabel('M_R')
% saveas(gca, [save_name,'_MODvalue_compare.png'])

