


%% compare MODs NORM minus Steady

load('norm_data.mat')
figure
subplot(3,3,1)
subplot(3,3,4)
hold on
plot(Mroll_norm*mod_value_all,Nm2mNmm*(Mx_mean_all-Mx_steady),'-k','linewidth',1)
plot(Mroll_norm*mod_value_freq,Nm2mNmm*(Mx_mean_freq-Mx_steady),'-k','linewidth',1)
plot(Mroll_norm*mod_value_stroke,Nm2mNmm*(Mx_mean_stroke-Mx_steady),'-k','linewidth',1)
plot(Mroll_norm*mod_value_pitch,Nm2mNmm*(Mx_mean_pitch-Mx_steady),'-k','linewidth',1)
plot(Mroll_norm*mod_value_dev,Nm2mNmm*(Mx_mean_dev-Mx_steady),'-k','linewidth',1)
% plot(Mroll_norm*mod_value_all(1:MODmin),Nm2mNmm*(Mx_mean_sum-Mx_steady),'-k','linewidth',1)

% all
for i = 2:2:length(mod_value_all)
    mod_now = mod_value_all(i);
    if mod_now>=0
        plot(Mroll_norm*mod_now,Nm2mNmm*(Mx_mean_all(i)-Mx_steady),'o',...
            'markerfacecolor',colormap_black2blue(round(10*mod_now)+1,:),...
            'markeredgecolor','k','markersize',4)
    end
end

% freq
for i = 2:2:length(mod_value_freq)
    mod_now = mod_value_freq(i);
    if mod_now>=0
        plot(Mroll_norm*mod_now,Nm2mNmm*(Mx_mean_freq(i)-Mx_steady),'v',...
            'markerfacecolor',colormap_black2blue(round(10*mod_now)+1,:),...
            'markeredgecolor','k','markersize',4)
    end
end

% stroke
for i = 2:2:length(mod_value_stroke)
    mod_now = mod_value_stroke(i);
    if mod_now>=0
        plot(Mroll_norm*mod_now,Nm2mNmm*(Mx_mean_stroke(i)-Mx_steady),'s',...
            'markerfacecolor',colormap_black2blue(round(10*mod_now)+1,:),...
            'markeredgecolor','k','markersize',4)
    end
end

% pitch
for i = 1:2:length(mod_value_pitch)
    mod_now = mod_value_pitch(i);
    if mod_now>=0
        plot(Mroll_norm*mod_now,Nm2mNmm*(Mx_mean_pitch(i)-Mx_steady),'^',...
            'markerfacecolor',colormap_black2blue(round(10*mod_now)+1,:),...
            'markeredgecolor','k','markersize',4)
    end
end

% dev
for i = 1:2:length(mod_value_dev)
    mod_now = mod_value_dev(i);
    if mod_now>=0
        plot(Mroll_norm*mod_now,Nm2mNmm*(Mx_mean_dev(i)-Mx_steady),'d',...
            'markerfacecolor',colormap_black2blue(round(10*mod_now)+1,:),...
            'markeredgecolor','k','markersize',4)
    end
end

xlabel('Mroll in (mNmm)','fontsize',10) 
ylabel('Mroll out (mNmm)','fontsize',10)

x_min = -.5e-3;
x_max = 4e-3;
y_min = -.5e-3;
y_max = 4e-3;
dy = 1e-3;

axis equal
axis([x_min x_max y_min y_max])
set(gca,'XTick',[0:dy:x_max],'fontsize',8) 
set(gca,'YTick',[0:dy:y_max],'fontsize',8) 

saveas(gca, 'Mroll_MSfig.fig')
saveas(gca, 'Mroll_MSfig.png')
plot2svg(['Mroll_MSfig.svg'])



%% all MODs
% figure
% hold on
% plot(mod_value_all,Nm2mNmm*Mx_mean_all,'.b')
% plot(mod_value_all,Nm2mNmm*My_mean_all,'.g')
% plot(mod_value_all,Nm2mNmm*Mz_mean_all,'.r')
% legend('Mroll','Mpitch','Myaw','location','nw')
% grid on
% xlabel('MOD value')
% ylabel('M (mNmm)')
% saveas(gca, 'Mrollaccel_MODvalue_allMOD.png')
% 
% % freq MODs
% figure
% hold on
% plot(freq_freq,Nm2mNmm*Mx_mean_freq,'.b')
% plot(freq_freq,Nm2mNmm*My_mean_freq,'.g')
% plot(freq_freq,Nm2mNmm*Mz_mean_freq,'.r')
% legend('Mroll','Mpitch','Myaw','location','nw') 
% grid on
% xlabel('freq')
% ylabel('M (mNmm)')
% saveas(gca, 'Mrollaccel_freq_freqMOD.png')
% 
% figure
% hold on
% plot(mod_value_freq,Nm2mNmm*Mx_mean_freq,'.b')
% plot(mod_value_freq,Nm2mNmm*My_mean_freq,'.g')
% plot(mod_value_freq,Nm2mNmm*Mz_mean_freq,'.r')
% legend('Mroll','Mpitch','Myaw','location','nw') 
% grid on
% xlabel('MOD value')
% ylabel('M (mNmm)')
% saveas(gca, 'Mrollaccel_MODvalue_freqMOD.png')
% 
% % stroke MODs
% figure
% hold on
% plot(Astroke_stroke,Nm2mNmm*Mx_mean_stroke,'.b')
% plot(Astroke_stroke,Nm2mNmm*My_mean_stroke,'.g')
% plot(Astroke_stroke,Nm2mNmm*Mz_mean_stroke,'.r')
% legend('Mroll','Mpitch','Myaw','location','nw') 
% grid on
% xlabel('Astroke')
% ylabel('M (mNmm)')
% saveas(gca, 'Mrollaccel_Astroke_strokeMOD.png')
% 
% figure
% hold on
% plot(mod_value_stroke,Nm2mNmm*Mx_mean_stroke,'.b')
% plot(mod_value_stroke,Nm2mNmm*My_mean_stroke,'.g')
% plot(mod_value_stroke,Nm2mNmm*Mz_mean_stroke,'.r')
% legend('Mroll','Mpitch','Myaw','location','nw') 
% grid on
% xlabel('MOD value')
% ylabel('M (mNmm)')
% saveas(gca, 'Mrollaccel_MODvalue_strokeMOD.png')
% 
% % pitch MODs
% figure
% hold on
% plot(Apitch_pitch,Nm2mNmm*Mx_mean_pitch,'.b')
% plot(Apitch_pitch,Nm2mNmm*My_mean_pitch,'.g')
% plot(Apitch_pitch,Nm2mNmm*Mz_mean_pitch,'.r')
% legend('Mroll','Mpitch','Myaw','location','nw') 
% grid on
% xlabel('Apitch')
% ylabel('M (mNmm)')
% saveas(gca, 'Mrollaccel_Apitch_pitchMOD.png')
% 
% figure
% hold on
% plot(mod_value_pitch,Nm2mNmm*Mx_mean_pitch,'.b')
% plot(mod_value_pitch,Nm2mNmm*My_mean_pitch,'.g')
% plot(mod_value_pitch,Nm2mNmm*Mz_mean_pitch,'.r')
% legend('Mroll','Mpitch','Myaw','location','nw') 
% grid on
% xlabel('MOD value')
% ylabel('M (mNmm)')
% saveas(gca, 'Mrollaccel_MODvalue_pitchMOD.png')
% 
% % dev MODs
% figure
% hold on
% plot(Adev_dev,Nm2mNmm*Mx_mean_dev,'.b')
% plot(Adev_dev,Nm2mNmm*My_mean_dev,'.g')
% plot(Adev_dev,Nm2mNmm*Mz_mean_dev,'.r')
% legend('Mroll','Mpitch','Myaw','location','nw') 
% grid on
% xlabel('Adev')
% ylabel('M (mNmm)')
% saveas(gca, 'Mrollaccel_Adev_devMOD.png')
% 
% figure
% hold on
% plot(mod_value_dev,Nm2mNmm*Mx_mean_dev,'.b')
% plot(mod_value_dev,Nm2mNmm*My_mean_dev,'.g')
% plot(mod_value_dev,Nm2mNmm*Mz_mean_dev,'.r')
% legend('Mroll','Mpitch','Myaw','location','nw') 
% grid on
% xlabel('MOD value')
% ylabel('M (mNmm)')
% saveas(gca, 'Mrollaccel_MODvalue_devMOD.png')
% 
% % compare MODs
% figure
% hold on
% plot(mod_value_all,Nm2mNmm*Mx_mean_all,'.k')
% plot(mod_value_freq,Nm2mNmm*Mx_mean_freq,'.r')
% plot(mod_value_stroke,Nm2mNmm*Mx_mean_stroke,'.g')
% plot(mod_value_pitch,Nm2mNmm*Mx_mean_pitch,'.b')
% plot(mod_value_dev,Nm2mNmm*Mx_mean_dev,'.c')
% plot(mod_value_all(1:MODmin),Nm2mNmm*Mx_mean_sum,'ok')
% 
% legend('all','freq','stroke','pitch','dev','sum','location','nw')
% grid on
% xlabel('MOD value')
% ylabel('Mroll (mNmm)')
% saveas(gca, 'Mrollaccel_MODvalue_compare.png')
% 
% % compare MODs NORM
% figure
% hold on
% plot(Mroll_norm*mod_value_all,Nm2mNmm*Mx_mean_all,'.-k')
% plot(Mroll_norm*mod_value_freq,Nm2mNmm*Mx_mean_freq,'.-r')
% plot(Mroll_norm*mod_value_stroke,Nm2mNmm*Mx_mean_stroke,'.-g')
% plot(Mroll_norm*mod_value_pitch,Nm2mNmm*Mx_mean_pitch,'.-b')
% plot(Mroll_norm*mod_value_dev,Nm2mNmm*Mx_mean_dev,'.-c')
% % plot(Mroll_norm*mod_value_all(1:MODmin),Nm2mNmm*Mx_mean_sum,'ok')
% 
% legend('all','freq','stroke','pitch','dev','sum','location','nw')
% grid on
% xlabel('MOD value')
% ylabel('Mroll (mNmm)')
% saveas(gca, 'Mrollaccel_MODvalue_compare.png')
% 
