%% compare MODs NORM minus Steady

% load('norm_data.mat')
figure
subplot(3,3,1)
hold on
plot(Fenhance_norm*mod_value_all+1,F_mean_all/Mg_fly,'-k','linewidth',1)
plot(Fenhance_norm*mod_value_freq+1,F_mean_freq/Mg_fly,'-k','linewidth',1)
plot(Fenhance_norm*mod_value_stroke+1,F_mean_stroke/Mg_fly,'-k','linewidth',1)
plot(Fenhance_norm*mod_value_pitch+1,F_mean_pitch/Mg_fly,'-k','linewidth',1)
plot(Fenhance_norm*mod_value_dev+1,F_mean_dev/Mg_fly,'-k','linewidth',1)
% plot(Fenhance_norm*mod_value_all(1:MODmin),F_mean_sum/Mg_fly,'-k','linewidth',1)

% all
for i = 2:2:length(mod_value_all)
    mod_now = mod_value_all(i);
    if mod_now>=0
        plot(Fenhance_norm*mod_now+1,F_mean_all(i)/Mg_fly,'o',...
            'markerfacecolor',colormap_black2blue(round(10*mod_now)+1,:),...
            'markeredgecolor','k','markersize',4)
    end
end

% freq
for i = 2:2:length(mod_value_freq)
    mod_now = mod_value_freq(i);
    if mod_now>=0
        plot(Fenhance_norm*mod_now+1,F_mean_freq(i)/Mg_fly,'v',...
            'markerfacecolor',colormap_black2blue(round(10*mod_now)+1,:),...
            'markeredgecolor','k','markersize',4)
    end
end

% stroke
for i = 1:2:length(mod_value_stroke)
    mod_now = mod_value_stroke(i);
    if mod_now>=0
        plot(Fenhance_norm*mod_now+1,F_mean_stroke(i)/Mg_fly,'s',...
            'markerfacecolor',colormap_black2blue(round(10*mod_now)+1,:),...
            'markeredgecolor','k','markersize',4)
    end
end

% pitch
for i = 2:2:length(mod_value_pitch)
    mod_now = mod_value_pitch(i);
    if mod_now>=0
        plot(Fenhance_norm*mod_now+1,F_mean_pitch(i)/Mg_fly,'^',...
            'markerfacecolor',colormap_black2blue(round(10*mod_now)+1,:),...
            'markeredgecolor','k','markersize',4)
    end
end

% dev
for i = 1:2:length(mod_value_dev)
    mod_now = mod_value_dev(i);
    if mod_now>=0
        plot(Fenhance_norm*mod_now+1,F_mean_dev(i)/Mg_fly,'d',...
            'markerfacecolor',colormap_black2blue(round(10*mod_now)+1,:),...
            'markeredgecolor','k','markersize',4)
    end
end

xlabel('F/Mg in (mNmm)','fontsize',10) 
ylabel('F/Mg out (mNmm)','fontsize',10)

x_min = 1;
x_max = 2;
y_min = .6;
y_max = 1.3;
dy = .1;

axis([x_min x_max y_min y_max])
set(gca,'XTick',[x_min:(x_max-x_min)/4:x_max],'fontsize',8) 
set(gca,'YTick',[y_min:dy:y_max],'fontsize',8) 
saveas(gca, 'Fenhance_MSfig.fig')
saveas(gca, 'Fenhance_MSfig.png')
plot2svg(['Fenhance_MSfig.svg'])







%% all MODs
% figure
% hold on
% plot(mod_value_all,F_mean_all/Mg_fly,'ok')
% plot(mod_value_all,Fx_mean_all/Mg_fly,'.b')
% plot(mod_value_all,Fy_mean_all/Mg_fly,'.g')
% plot(mod_value_all,-Fz_mean_all/Mg_fly,'.r')
% legend('Ftotal','Thrust','Fside','Lift','location','nw')
% grid on
% xlabel('MOD value')
% ylabel('F/Mg')
% saveas(gca, 'FMg_MODvalue_allMOD.png')
% 
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
% saveas(gca, 'FMg_freq_freqMOD.png')
% 
% figure
% hold on
% plot(mod_value_freq,F_mean_freq/Mg_fly,'ok')
% plot(mod_value_freq,Fx_mean_freq/Mg_fly,'.b')
% plot(mod_value_freq,Fy_mean_freq/Mg_fly,'.g')
% plot(mod_value_freq,-Fz_mean_freq/Mg_fly,'.r')
% legend('Ftotal','Thrust','Fside','Lift','location','nw')
% grid on
% xlabel('MOD value')
% ylabel('F/Mg')
% saveas(gca, 'FMg_MODvalue_freqMOD.png')
% 
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
% saveas(gca, 'FMg_Astroke_strokeMOD.png')
% 
% figure
% hold on
% plot(mod_value_stroke,F_mean_stroke/Mg_fly,'ok')
% plot(mod_value_stroke,Fx_mean_stroke/Mg_fly,'.b')
% plot(mod_value_stroke,Fy_mean_stroke/Mg_fly,'.g')
% plot(mod_value_stroke,-Fz_mean_stroke/Mg_fly,'.r')
% legend('Ftotal','Thrust','Fside','Lift','location','nw')
% grid on
% xlabel('MOD value')
% ylabel('F/Mg')
% saveas(gca, 'FMg_MODvalue_strokeMOD.png')
% 
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
% saveas(gca, 'FMg_Apitch_pitchMOD.png')
% 
% figure
% hold on
% plot(mod_value_pitch,F_mean_pitch/Mg_fly,'ok')
% plot(mod_value_pitch,Fx_mean_pitch/Mg_fly,'.b')
% plot(mod_value_pitch,Fy_mean_pitch/Mg_fly,'.g')
% plot(mod_value_pitch,-Fz_mean_pitch/Mg_fly,'.r')
% legend('Ftotal','Thrust','Fside','Lift','location','nw')
% grid on
% xlabel('MOD value')
% ylabel('F/Mg')
% saveas(gca, 'FMg_MODvalue_pitchMOD.png')
% 
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
% saveas(gca, 'FMg_Adev_devMOD.png')
% 
% figure
% hold on
% plot(mod_value_dev,F_mean_dev/Mg_fly,'ok')
% plot(mod_value_dev,Fx_mean_dev/Mg_fly,'.b')
% plot(mod_value_dev,Fy_mean_dev/Mg_fly,'.g')
% plot(mod_value_dev,-Fz_mean_dev/Mg_fly,'.r')
% legend('Ftotal','Thrust','Fside','Lift','location','nw')
% grid on
% xlabel('MOD value')
% ylabel('F/Mg')
% saveas(gca, 'FMg_MODvalue_devMOD.png')
% 
% % compare MODs
% figure
% hold on
% plot(mod_value_all,F_mean_all/Mg_fly,'.k')
% plot(mod_value_freq,F_mean_freq/Mg_fly,'.r')
% plot(mod_value_stroke,F_mean_stroke/Mg_fly,'.g')
% plot(mod_value_pitch,F_mean_pitch/Mg_fly,'.b')
% plot(mod_value_dev,F_mean_dev/Mg_fly,'.c')
% plot(mod_value_all(1:MODmin),F_mean_sum/Mg_fly,'ok')
% plot(mod_value_all(1:MODmin),F_mean_sum_nopitch/Mg_fly,'ob')
% 
% legend('all','freq','stroke','pitch','dev','sum','no pitch','location','nw')
% grid on
% xlabel('MOD value')
% ylabel('F/Mg')
% saveas(gca, 'FMg_MODvalue_compare.png')
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
% saveas(gca, 'Fnorm_MODvalue_compare.png')
% 
% 
