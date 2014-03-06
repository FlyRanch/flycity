% all MODs
figure
hold on
plot(mod_value_all,F_mean_all/Mg_fly,'ok')
plot(mod_value_all,Fx_mean_all/Mg_fly,'.b')
plot(mod_value_all,Fy_mean_all/Mg_fly,'.g')
plot(mod_value_all,-Fz_mean_all/Mg_fly,'.r')
legend('Ftotal','Thrust','Fside','Lift','location','nw')
grid on
xlabel('MOD value')
ylabel('F/Mg')
saveas(gca, [save_name,'_MODvalue_allMOD.png'])

% % freq MODs
% figure
% hold on
% plot(freq_freq,F_mean_freq/Mg_fly,'ok')
% plot(freq_freq,Fx_mean_freq/Mg_fly,'.b')
% plot(freq_freq,Fy_mean_freq/Mg_fly,'.g')
% plot(freq_freq,-Fz_mean_freq/Mg_fly,'.r')
% legend('Ftotal','Thrust','Fside','Lift','location','nw')
% grid on
% xlabel('freq')
% ylabel('F/Mg')
% saveas(gca, [save_name,'_freq_freqMOD.png'])

figure
hold on
plot(mod_value_freq,F_mean_freq/Mg_fly,'ok')
plot(mod_value_freq,Fx_mean_freq/Mg_fly,'.b')
plot(mod_value_freq,Fy_mean_freq/Mg_fly,'.g')
plot(mod_value_freq,-Fz_mean_freq/Mg_fly,'.r')
legend('Ftotal','Thrust','Fside','Lift','location','nw')
grid on
xlabel('MOD value')
ylabel('F/Mg')
saveas(gca, [save_name,'_MODvalue_freqMOD.png'])

figure
hold on
plot(mod_value_vel,F_mean_vel/Mg_fly,'ok')
plot(mod_value_vel,Fx_mean_vel/Mg_fly,'.b')
plot(mod_value_vel,Fy_mean_vel/Mg_fly,'.g')
plot(mod_value_vel,-Fz_mean_vel/Mg_fly,'.r')
legend('Ftotal','Thrust','Fside','Lift','location','nw')
grid on
xlabel('MOD value')
ylabel('F/Mg')
saveas(gca, [save_name,'_MODvalue_velMOD.png'])

% % stroke MODs
% figure
% hold on
% plot(Astroke_stroke,F_mean_stroke/Mg_fly,'ok')
% plot(Astroke_stroke,Fx_mean_stroke/Mg_fly,'.b')
% plot(Astroke_stroke,Fy_mean_stroke/Mg_fly,'.g')
% plot(Astroke_stroke,-Fz_mean_stroke/Mg_fly,'.r')
% legend('Ftotal','Thrust','Fside','Lift','location','nw')
% grid on
% xlabel('Astroke')
% ylabel('F/Mg')
% saveas(gca, [save_name,'_Astroke_strokeMOD.png'])

figure
hold on
plot(mod_value_stroke,F_mean_stroke/Mg_fly,'ok')
plot(mod_value_stroke,Fx_mean_stroke/Mg_fly,'.b')
plot(mod_value_stroke,Fy_mean_stroke/Mg_fly,'.g')
plot(mod_value_stroke,-Fz_mean_stroke/Mg_fly,'.r')
legend('Ftotal','Thrust','Fside','Lift','location','nw')
grid on
xlabel('MOD value')
ylabel('F/Mg')
saveas(gca, [save_name,'_MODvalue_strokeMOD.png'])

% % pitch MODs
% figure
% hold on
% plot(Apitch_pitch,F_mean_pitch/Mg_fly,'ok')
% plot(Apitch_pitch,Fx_mean_pitch/Mg_fly,'.b')
% plot(Apitch_pitch,Fy_mean_pitch/Mg_fly,'.g')
% plot(Apitch_pitch,-Fz_mean_pitch/Mg_fly,'.r')
% legend('Ftotal','Thrust','Fside','Lift','location','nw')
% grid on
% xlabel('Apitch')
% ylabel('F/Mg')
% saveas(gca, [save_name,'_Apitch_pitchMOD.png'])
% 
figure
hold on
plot(mod_value_pitch,F_mean_pitch/Mg_fly,'ok')
plot(mod_value_pitch,Fx_mean_pitch/Mg_fly,'.b')
plot(mod_value_pitch,Fy_mean_pitch/Mg_fly,'.g')
plot(mod_value_pitch,-Fz_mean_pitch/Mg_fly,'.r')
legend('Ftotal','Thrust','Fside','Lift','location','nw')
grid on
xlabel('MOD value')
ylabel('F/Mg')
saveas(gca, [save_name,'_MODvalue_pitchMOD.png'])

% % dev MODs
% figure
% hold on
% plot(Adev_dev,F_mean_dev/Mg_fly,'ok')
% plot(Adev_dev,Fx_mean_dev/Mg_fly,'.b')
% plot(Adev_dev,Fy_mean_dev/Mg_fly,'.g')
% plot(Adev_dev,-Fz_mean_dev/Mg_fly,'.r')
% legend('Ftotal','Thrust','Fside','Lift','location','nw')
% grid on
% xlabel('Adev')
% ylabel('F/Mg')
% saveas(gca, [save_name,'_Adev_devMOD.png'])

figure
hold on
plot(mod_value_dev,F_mean_dev/Mg_fly,'ok')
plot(mod_value_dev,Fx_mean_dev/Mg_fly,'.b')
plot(mod_value_dev,Fy_mean_dev/Mg_fly,'.g')
plot(mod_value_dev,-Fz_mean_dev/Mg_fly,'.r')
legend('Ftotal','Thrust','Fside','Lift','location','nw')
grid on
xlabel('MOD value')
ylabel('F/Mg')
saveas(gca, [save_name,'_MODvalue_devMOD.png'])

% compare MODs
figure
hold on
plot(mod_value_all,F_mean_all/Mg_fly,'.k')
plot(mod_value_freq,F_mean_freq/Mg_fly,'.r')
plot(mod_value_vel,F_mean_vel/Mg_fly,'.m')
plot(mod_value_stroke,F_mean_stroke/Mg_fly,'.g')
plot(mod_value_pitch,F_mean_pitch/Mg_fly,'.b')
plot(mod_value_dev,F_mean_dev/Mg_fly,'.c')
plot(mod_value_all(1:MODmin),F_mean_sum/Mg_fly,'ok')
plot(mod_value_all(1:MODmin),F_mean_sum_nopitch/Mg_fly,'ob')

legend('all','freq','vel','stroke','pitch','dev','sum','no pitch','location','nw')
grid on
xlabel('MOD value')
ylabel('F/Mg')
saveas(gca, [save_name,'_MODvalue_compare.png'])
% 
% figure
% hold on
% plot(mod_value_all,F_mean_all/F_steady,'.k')
% plot(mod_value_freq,F_mean_freq/F_steady,'.r')
% plot(mod_value_stroke,F_mean_stroke/F_steady,'.g')
% plot(mod_value_pitch,F_mean_pitch/F_steady,'.b')
% plot(mod_value_dev,F_mean_dev/F_steady,'.c')
% plot(mod_value_all(1:MODmin),F_mean_sum/F_steady,'ok')
% plot(mod_value_all(1:MODmin),F_mean_sum_nopitch/F_steady,'ob')
% 
% legend('all','freq','stroke','pitch','dev','sum','no pitch','location','nw')
% grid on
% xlabel('MOD value')
% ylabel('F/Fsteady')
% saveas(gca, 'Fnorm_MODvalue_compare.png'])