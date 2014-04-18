% make borf plots mean F&M

clear
clc
close all

%% load MOD&norm data
load('MOD_norm_data.mat')

%% Fenhance
load(FDB_mean)

x_min = .8;
x_max = 1 + max(Fenhanses);
dx = max(Fenhanses)/2;
y_min = x_min;
y_max = x_max;
dy = dx;

dmod = mod_value_all(2)-mod_value_all(1);
subset_all = find(Fenhance_norm*mod_value_all>=x_min-1-dmod & Fenhance_norm*mod_value_all<=x_max-1+dmod );
subset_freq = find(Fenhance_norm*mod_value_freq>=x_min-1-dmod & Fenhance_norm*mod_value_freq<=x_max-1+dmod );
subset_stroke = find(Fenhance_norm*mod_value_stroke>=x_min-1-dmod & Fenhance_norm*mod_value_stroke<=x_max-1+dmod );
subset_pitch = find(Fenhance_norm*mod_value_pitch>=x_min-1-dmod & Fenhance_norm*mod_value_pitch<=x_max-1+dmod );
subset_dev = find(Fenhance_norm*mod_value_dev>=x_min-1-dmod & Fenhance_norm*mod_value_dev<=x_max-1+dmod );

figure
subplot(3,3,1)
hold on

plot(Fenhance_norm*mod_value_all(subset_all)+1,F_mean_all(subset_all)/Mg_fly,'-r','linewidth',1)
plot(Fenhance_norm*mod_value_freq(subset_freq)+1,F_mean_freq(subset_freq)/Mg_fly,'-g','linewidth',1)
plot(Fenhance_norm*mod_value_stroke(subset_stroke)+1,F_mean_stroke(subset_stroke)/Mg_fly,'-b','linewidth',1)
plot(Fenhance_norm*mod_value_pitch(subset_pitch)+1,F_mean_pitch(subset_pitch)/Mg_fly,'-c','linewidth',1)
plot(Fenhance_norm*mod_value_dev(subset_dev)+1,F_mean_dev(subset_dev)/Mg_fly,'-m','linewidth',1)
plot(Fenhance_norm*mod_value_all((subset_all))+1,F_mean_sum(subset_all)/Mg_fly,'--','color',[.5 .5 .5],'linewidth',1)

xlabel('a+g/g','fontsize',10) 
ylabel('F/Mg','fontsize',10)


% trendline
plot([x_min x_max],[x_min x_max],'--k')

axis equal
axis([x_min x_max y_min y_max])
set(gca,'XTick',[1:dx:x_max],'fontsize',8) 
set(gca,'YTick',[1:dy:y_max],'fontsize',8) 


%% Mroll
load(rollDB_mean)

x_min = -1;
x_max = max(rollaccels);
dx = max(rollaccels)/2;
y_min = -.25;
y_max = 1;
dy = .5;

dmod = mod_value_all(2)-mod_value_all(1);
subset_allNOfreq = find(rollaccel_norm*mod_value_allNOfreq>=x_min-dmod & rollaccel_norm*mod_value_allNOfreq<=x_max+dmod );
% subset_freq = find(rollaccel_norm*mod_value_freq>=x_min-dmod & rollaccel_norm*mod_value_freq<=x_max+dmod );
subset_stroke = find(rollaccel_norm*mod_value_stroke>=x_min-dmod & rollaccel_norm*mod_value_stroke<=x_max+dmod );
subset_pitch = find(rollaccel_norm*mod_value_pitch>=x_min-dmod & rollaccel_norm*mod_value_pitch<=x_max+dmod );
subset_dev = find(rollaccel_norm*mod_value_dev>=x_min-dmod & rollaccel_norm*mod_value_dev<=x_max+dmod );

% subplot(3,3,4)
% hold on
% plot(rollaccel_norm*mod_value_allNOfreq(subset_allNOfreq),(Mx_mean_allNOfreq(subset_allNOfreq)-Mx_steady)/f_wb_steady^2/Lwing^5,'-r','linewidth',1)
% % plot(rollaccel_norm*mod_value_freq(subset_freq),(Mx_mean_freq(subset_freq)-Mx_steady)/f_wb_steady^2/Lwing^5,'-g','linewidth',1)
% plot(rollaccel_norm*mod_value_stroke(subset_stroke),(Mx_mean_stroke(subset_stroke)-Mx_steady)/f_wb_steady^2/Lwing^5,'-b','linewidth',1)
% plot(rollaccel_norm*mod_value_pitch(subset_pitch),(Mx_mean_pitch(subset_pitch)-Mx_steady)/f_wb_steady^2/Lwing^5,'-c','linewidth',1)
% plot(rollaccel_norm*mod_value_dev(subset_dev),(Mx_mean_dev(subset_dev)-Mx_steady)/f_wb_steady^2/Lwing^5,'-m','linewidth',1)
% plot(rollaccel_norm*mod_value_allNOfreq(subset_allNOfreq),(Mx_mean_sum(subset_allNOfreq)-Mx_steady)/f_wb_steady^2/Lwing^5,'--','color',[.5 .5 .5],'linewidth',1)

subplot(3,3,4)
hold on
plot(rollaccel_norm*mod_value_allNOfreq(subset_allNOfreq),(Mx_mean_allNOfreq(subset_allNOfreq)-Mx_mean_allNOfreq(mod_value_allNOfreq==0))/f_wb_steady^2/Lwing^5,'-r','linewidth',1)
% plot(rollaccel_norm*mod_value_freq(subset_freq),(Mx_mean_freq(subset_freq)-Mx_mean_freq(mod_value_freq==0))/f_wb_steady^2/Lwing^5,'-g','linewidth',1)
plot(rollaccel_norm*mod_value_stroke(subset_stroke),(Mx_mean_stroke(subset_stroke)-Mx_mean_stroke(mod_value_stroke==0))/f_wb_steady^2/Lwing^5,'-b','linewidth',1)
plot(rollaccel_norm*mod_value_pitch(subset_pitch),(Mx_mean_pitch(subset_pitch)-Mx_mean_pitch(mod_value_pitch==0))/f_wb_steady^2/Lwing^5,'-c','linewidth',1)
plot(rollaccel_norm*mod_value_dev(subset_dev),(Mx_mean_dev(subset_dev)-Mx_mean_dev(mod_value_dev==0))/f_wb_steady^2/Lwing^5,'-m','linewidth',1)
plot(rollaccel_norm*mod_value_allNOfreq(subset_allNOfreq),(Mx_mean_sum(subset_allNOfreq)-Mx_mean_sum(mod_value_allNOfreq==0))/f_wb_steady^2/Lwing^5,'--','color',[.5 .5 .5],'linewidth',1)

xlabel('a roll norm','fontsize',10) 
ylabel('Mroll norm','fontsize',10)

% trendline
aRoll_norm = rollaccel_norm*mod_value_allNOfreq(subset_allNOfreq);
Mroll_norm = (Mx_mean_allNOfreq(subset_allNOfreq)-Mx_mean_allNOfreq(mod_value_allNOfreq==0))/f_wb_steady^2/Lwing^5;
aMroll_fit = polyfit(aRoll_norm,Mroll_norm,1);
Mroll_fit = polyval(aMroll_fit,[x_min x_max]);
plot([x_min x_max],Mroll_fit,'--k')

% axis equal
axis square
axis([x_min x_max y_min y_max])
set(gca,'XTick',[0:dx:x_max],'fontsize',8) 
set(gca,'YTick',[0:dy:y_max],'fontsize',8) 

%% Myaw
% load(yawDB_mean)
% 
% x_min = -.5;
% x_max = max(yawaccels);
% dx = max(yawaccels)/2;
% y_min = -.2;
% y_max = .8;
% dy = .4;
% 
% dmod = mod_value_all(2)-mod_value_all(1);
% subset_allNOfreq = find(yawaccel_norm*mod_value_allNOfreq>=x_min-dmod & yawaccel_norm*mod_value_allNOfreq<=x_max+dmod );
% % subset_freq = find(yawaccel_norm*mod_value_freq>=x_min-dmod & yawaccel_norm*mod_value_freq<=x_max+dmod );
% subset_stroke = find(yawaccel_norm*mod_value_stroke>=x_min-dmod & yawaccel_norm*mod_value_stroke<=x_max+dmod );
% subset_pitch = find(yawaccel_norm*mod_value_pitch>=x_min-dmod & yawaccel_norm*mod_value_pitch<=x_max+dmod );
% subset_dev = find(yawaccel_norm*mod_value_dev>=x_min-dmod & yawaccel_norm*mod_value_dev<=x_max+dmod );
% 
% % subplot(3,3,5)
% % hold on
% % plot(yawaccel_norm*mod_value_allNOfreq(subset_allNOfreq),(Mz_mean_allNOfreq(subset_allNOfreq)-Mz_steady)/f_wb_steady^2/Lwing^5,'-r','linewidth',1)
% % % plot(yawaccel_norm*mod_value_freq(subset_freq),(Mz_mean_freq(subset_freq)-Mz_steady)/f_wb_steady^2/Lwing^5,'-g','linewidth',1)
% % plot(yawaccel_norm*mod_value_stroke(subset_stroke),(Mz_mean_stroke(subset_stroke)-Mz_steady)/f_wb_steady^2/Lwing^5,'-b','linewidth',1)
% % plot(yawaccel_norm*mod_value_pitch(subset_pitch),(Mz_mean_pitch(subset_pitch)-Mz_steady)/f_wb_steady^2/Lwing^5,'-c','linewidth',1)
% % plot(yawaccel_norm*mod_value_dev(subset_dev),(Mz_mean_dev(subset_dev)-Mz_steady)/f_wb_steady^2/Lwing^5,'-m','linewidth',1)
% % plot(yawaccel_norm*mod_value_allNOfreq(subset_allNOfreq),(Mz_mean_sum(subset_allNOfreq)-Mz_steady)/f_wb_steady^2/Lwing^5,'--','color',[.5 .5 .5],'linewidth',1)
% 
% subplot(3,3,5)
% hold on
% plot(yawaccel_norm*mod_value_allNOfreq(subset_allNOfreq),-(Mz_mean_allNOfreq(subset_allNOfreq)-Mz_mean_allNOfreq(mod_value_allNOfreq==0))/f_wb_steady^2/Lwing^5,'-r','linewidth',1)
% % plot(yawaccel_norm*mod_value_freq(subset_freq),-(Mz_mean_freq(subset_freq)-Mz_mean_freq(mod_value_freq==0))/f_wb_steady^2/Lwing^5,'-g','linewidth',1)
% plot(yawaccel_norm*mod_value_stroke(subset_stroke),-(Mz_mean_stroke(subset_stroke)-Mz_mean_stroke(mod_value_stroke==0))/f_wb_steady^2/Lwing^5,'-b','linewidth',1)
% plot(yawaccel_norm*mod_value_pitch(subset_pitch),-(Mz_mean_pitch(subset_pitch)-Mz_mean_pitch(mod_value_pitch==0))/f_wb_steady^2/Lwing^5,'-c','linewidth',1)
% plot(yawaccel_norm*mod_value_dev(subset_dev),-(Mz_mean_dev(subset_dev)-Mz_mean_dev(mod_value_dev==0))/f_wb_steady^2/Lwing^5,'-m','linewidth',1)
% plot(yawaccel_norm*mod_value_allNOfreq(subset_allNOfreq),-(Mz_mean_sum(subset_allNOfreq)-Mz_mean_sum(mod_value_allNOfreq==0))/f_wb_steady^2/Lwing^5,'--','color',[.5 .5 .5],'linewidth',1)
% 
% xlabel('a yaw norm','fontsize',10) 
% ylabel('Myaw norm','fontsize',10)
% 
% % trendline
% aYaw_norm = yawaccel_norm*mod_value_allNOfreq(subset_allNOfreq);
% Myaw_norm = -(Mz_mean_allNOfreq(subset_allNOfreq)-Mz_mean_allNOfreq(mod_value_allNOfreq==0))/f_wb_steady^2/Lwing^5;
% aMyaw_fit = polyfit(aYaw_norm,Myaw_norm,1);
% Myaw_fit = polyval(aMyaw_fit,[x_min x_max]);
% plot([x_min x_max],Myaw_fit,'--k')
% 
% % axis equal
% axis square
% axis([x_min x_max y_min y_max])
% set(gca,'XTick',[0:dx:x_max],'fontsize',8) 
% set(gca,'YTick',[0:dy:y_max],'fontsize',8) 

%% Mpitch
load(pitchDB_mean)

x_min = -max(pitchaccels);
x_max = max(pitchaccels);
dx = max(pitchaccels);
y_min = -.4;
y_max = .4;
dy = .2;

dmod = mod_value_all(2)-mod_value_all(1);
subset_allNOfreq = find(pitchaccel_norm*mod_value_allNOfreq>=x_min-dmod & pitchaccel_norm*mod_value_allNOfreq<=x_max+dmod );
% subset_freq = find(pitchaccel_norm*mod_value_freq>=x_min-dmod & pitchaccel_norm*mod_value_freq<=x_max+dmod );
subset_stroke = find(pitchaccel_norm*mod_value_stroke>=x_min-dmod & pitchaccel_norm*mod_value_stroke<=x_max+dmod );
subset_pitch = find(pitchaccel_norm*mod_value_pitch>=x_min-dmod & pitchaccel_norm*mod_value_pitch<=x_max+dmod );
subset_dev = find(pitchaccel_norm*mod_value_dev>=x_min-dmod & pitchaccel_norm*mod_value_dev<=x_max+dmod );

% subplot(3,3,7)
% hold on
% plot(pitchaccel_norm*mod_value_allNOfreq(subset_allNOfreq),(My_mean_allNOfreq(subset_allNOfreq)-My_steady)/f_wb_steady^2/Lwing^5,'-r','linewidth',1)
% % plot(pitchaccel_norm*mod_value_freq(subset_freq),(My_mean_freq(subset_freq)-My_steady)/f_wb_steady^2/Lwing^5,'-g','linewidth',1)
% plot(pitchaccel_norm*mod_value_stroke(subset_stroke),(My_mean_stroke(subset_stroke)-My_steady)/f_wb_steady^2/Lwing^5,'-b','linewidth',1)
% plot(pitchaccel_norm*mod_value_pitch(subset_pitch),(My_mean_pitch(subset_pitch)-My_steady)/f_wb_steady^2/Lwing^5,'-c','linewidth',1)
% plot(pitchaccel_norm*mod_value_dev(subset_dev),(My_mean_dev(subset_dev)-My_steady)/f_wb_steady^2/Lwing^5,'-m','linewidth',1)
% plot(pitchaccel_norm*mod_value_allNOfreq(subset_allNOfreq),(My_mean_sum(subset_allNOfreq)-My_steady)/f_wb_steady^2/Lwing^5,'--','color',[.5 .5 .5],'linewidth',1)

subplot(3,3,7)
hold on
plot(pitchaccel_norm*mod_value_allNOfreq(subset_allNOfreq),(My_mean_allNOfreq(subset_allNOfreq)-My_mean_allNOfreq(mod_value_allNOfreq==0))/f_wb_steady^2/Lwing^5,'-r','linewidth',1)
% plot(pitchaccel_norm*mod_value_freq(subset_freq),(My_mean_freq(subset_freq)-My_mean_freq(mod_value_freq==0))/f_wb_steady^2/Lwing^5,'-g','linewidth',1)
plot(pitchaccel_norm*mod_value_stroke(subset_stroke),(My_mean_stroke(subset_stroke)-My_mean_stroke(mod_value_stroke==0))/f_wb_steady^2/Lwing^5,'-b','linewidth',1)
plot(pitchaccel_norm*mod_value_pitch(subset_pitch),(My_mean_pitch(subset_pitch)-My_mean_pitch(mod_value_pitch==0))/f_wb_steady^2/Lwing^5,'-c','linewidth',1)
plot(pitchaccel_norm*mod_value_dev(subset_dev),(My_mean_dev(subset_dev)-My_mean_dev(mod_value_dev==0))/f_wb_steady^2/Lwing^5,'-m','linewidth',1)
plot(pitchaccel_norm*mod_value_allNOfreq(subset_allNOfreq),(My_mean_sum(subset_allNOfreq)-My_mean_sum(mod_value_allNOfreq==0))/f_wb_steady^2/Lwing^5,'--','color',[.5 .5 .5],'linewidth',1)

xlabel('a pitch norm','fontsize',10) 
ylabel('Mpitch norm','fontsize',10)

% trendline
aPitch_norm = pitchaccel_norm*mod_value_allNOfreq(subset_allNOfreq);
Mpitch_norm = (My_mean_allNOfreq(subset_allNOfreq)-My_mean_allNOfreq(mod_value_allNOfreq==0))/f_wb_steady^2/Lwing^5;
aMpitch_fit = polyfit(aPitch_norm,Mpitch_norm,1);
Mpitch_fit = polyval(aMpitch_fit,[x_min x_max]);
plot([x_min x_max],Mpitch_fit,'--k')

% axis equal
% axis([x_min x_max y_min y_max])
% set(gca,'XTick',[0:dy:x_max],'fontsize',8) 
% set(gca,'YTick',[0:dy:y_max],'fontsize',8) 

% axis equal
axis square
axis([-x_max x_max -y_max y_max])
set(gca,'XTick',[-x_max:dx:x_max],'fontsize',8) 
set(gca,'YTick',[-x_max:dy:y_max],'fontsize',8) 

%% save plot
% mkdir('MSfigs')
% cd('MSfigs')
% saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch.fig')
% saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch.png')
% plot2svg(['MSfig_Fenhance_Mroll_Myaw_Mpitch.svg'])
% cd ..
% 
% %% save trenddata
% save('InertiaMoment_trends.mat','aRoll_norm','Mroll_norm','aMroll_fit','aYaw_norm','Myaw_norm','aMyaw_fit','aPitch_norm','Mpitch_norm','aMpitch_fit','f_wb_steady','Lwing')
