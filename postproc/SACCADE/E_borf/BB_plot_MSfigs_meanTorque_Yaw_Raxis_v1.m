% make borf plots mean F&M

clc
close all

%% Myaw
clear
load('MOD_norm_data.mat')
load(yawDB_mean)

x_min = -.02;
x_max = .08;
dx = x_max/4;
y_min = x_min;
y_max = x_max;
dy = dx;

dmod = mod_value_allNOfreq(2)-mod_value_allNOfreq(1);
subset_allNOfreq = find(Myaw_norm*mod_value_allNOfreq>=x_min-dmod & Myaw_norm*mod_value_allNOfreq<=x_max+dmod );
% subset_freq = find(Myaw_norm*mod_value_freq>=x_min-dmod & Myaw_norm*mod_value_freq<=x_max+dmod );
subset_stroke = find(Myaw_norm*mod_value_stroke>=x_min-dmod & Myaw_norm*mod_value_stroke<=x_max+dmod );
subset_pitch = find(Myaw_norm*mod_value_pitch>=x_min-dmod & Myaw_norm*mod_value_pitch<=x_max+dmod );
subset_dev = find(Myaw_norm*mod_value_dev>=x_min-dmod & Myaw_norm*mod_value_dev<=x_max+dmod );

% subplot(2,2,3)
% hold on
% plot(Myaw_norm*mod_value_allNOfreq(subset_allNOfreq),(Mz_mean_allNOfreq(subset_allNOfreq)-Mz_steady)/Mg_fly/Lwing,'-or','linewidth',1,'markersize',5,'markerfacecolor','r')
% % plot(Myaw_norm*mod_value_freq(subset_freq),(Mz_mean_freq(subset_freq)-Mz_steady)/Mg_fly/Lwing,'-og','linewidth',1,'markersize',5,'markerfacecolor','g')
% plot(Myaw_norm*mod_value_stroke(subset_stroke),(Mz_mean_stroke(subset_stroke)-Mz_steady)/Mg_fly/Lwing,'-ob','linewidth',1,'markersize',5,'markerfacecolor','b')
% plot(Myaw_norm*mod_value_pitch(subset_pitch),(Mz_mean_pitch(subset_pitch)-Mz_steady)/Mg_fly/Lwing,'-oc','linewidth',1,'markersize',5,'markerfacecolor','c')
% plot(Myaw_norm*mod_value_dev(subset_dev),(Mz_mean_dev(subset_dev)-Mz_steady)/Mg_fly/Lwing,'-om','linewidth',1,'markersize',5,'markerfacecolor','m')
% plot(Myaw_norm*mod_value_allNOfreq(subset_allNOfreq),(Mz_mean_sum(subset_allNOfreq)-Mz_steady)/Mg_fly/Lwing,'-o','color',[.5 .5 .5],'linewidth',1,'markersize',5,'markerfacecolor',[.5 .5 .5])

subplot(2,2,3)
hold on
plot(Myaw_norm*mod_value_allNOfreq(subset_allNOfreq),(Mz_mean_allNOfreq(subset_allNOfreq)-Mz_mean_allNOfreq(mod_value_allNOfreq==0))/Mg_fly/Lwing,'-or','linewidth',1,'markersize',5,'markerfacecolor','r')
% plot(Myaw_norm*mod_value_freq(subset_freq),(Mz_mean_freq(subset_freq)-Mz_mean_freq(mod_value_freq==0))/Mg_fly/Lwing,'-og','linewidth',1,'markersize',5,'markerfacecolor','g')
plot(Myaw_norm*mod_value_stroke(subset_stroke),(Mz_mean_stroke(subset_stroke)-Mz_mean_stroke(mod_value_stroke==0))/Mg_fly/Lwing,'-ob','linewidth',1,'markersize',5,'markerfacecolor','b')
plot(Myaw_norm*mod_value_pitch(subset_pitch),(Mz_mean_pitch(subset_pitch)-Mz_mean_pitch(mod_value_pitch==0))/Mg_fly/Lwing,'-oc','linewidth',1,'markersize',5,'markerfacecolor','c')
plot(Myaw_norm*mod_value_dev(subset_dev),(Mz_mean_dev(subset_dev)-Mz_mean_dev(mod_value_dev==0))/Mg_fly/Lwing,'-om','linewidth',1,'markersize',5,'markerfacecolor','m')
plot(Myaw_norm*mod_value_allNOfreq(subset_allNOfreq),(Mz_mean_sum(subset_allNOfreq)-Mz_mean_sum(mod_value_allNOfreq==0))/Mg_fly/Lwing,'-o','color',[.5 .5 .5],'linewidth',1,'markersize',5,'markerfacecolor',[.5 .5 .5])

% trendline

xlabel('Myaw/(m+a)g fly','fontsize',10) 
ylabel('Myaw/(m+a)g borf','fontsize',10)

% trendline
plot([x_min x_max],[x_min x_max],'-k')

% aYaw_norm = Myaw_norm*mod_value_allNOfreq(subset_allNOfreq);
% Myaw_norm = -(Mz_mean_allNOfreq(subset_allNOfreq)-Mz_mean_allNOfreq(mod_value_allNOfreq==0))/Mg_fly/Lwing;
% aMyaw_fit = polyfit(aYaw_norm,Myaw_norm,1);
% Myaw_fit = polyval(aMyaw_fit,[x_min x_max]);
% plot([x_min x_max],Myaw_fit,'-k')
% 
% % save trend data
% save('InertiaMoment_trend_yaw.mat','aYaw_norm','Myaw_norm','aMyaw_fit','f_wb_steady','Lwing')

axis equal
% axis square
axis([x_min x_max y_min y_max])
set(gca,'XTick',[0:dx:x_max],'fontsize',8) 
set(gca,'YTick',[0:dy:y_max],'fontsize',8) 


%% MaxisR
clear
load('MOD_norm_data.mat')
load(RaxisDB_mean)

x_min = -.02;
x_max = .06;
dx = x_max/3;
y_min = x_min;
y_max = x_max;
dy = dx;

dmod = mod_value_allNOfreq(2)-mod_value_allNOfreq(1);
subset_allNOfreq = find(M_R_norm*mod_value_allNOfreq>=x_min-dmod & M_R_norm*mod_value_allNOfreq<=x_max+dmod );
% subset_freq = find(M_R_norm*mod_value_freq>=x_min-dmod & M_R_norm*mod_value_freq<=x_max+dmod );
subset_stroke = find(M_R_norm*mod_value_stroke>=x_min-dmod & M_R_norm*mod_value_stroke<=x_max+dmod );
subset_pitch = find(M_R_norm*mod_value_pitch>=x_min-dmod & M_R_norm*mod_value_pitch<=x_max+dmod );
subset_dev = find(M_R_norm*mod_value_dev>=x_min-dmod & M_R_norm*mod_value_dev<=x_max+dmod );

% subplot(2,2,4)
% hold on
% plot(M_R_norm*mod_value_allNOfreq(subset_allNOfreq),(M_R_mean_allNOfreq(subset_allNOfreq)-M_R_steady)/Mg_fly/Lwing,'-or','linewidth',1,'markersize',5,'markerfacecolor','r')
% % plot(M_R_norm*mod_value_freq(subset_freq),(M_R_mean_freq(subset_freq)-M_R_steady)/Mg_fly/Lwing,'-og','linewidth',1,'markersize',5,'markerfacecolor','g')
% plot(M_R_norm*mod_value_stroke(subset_stroke),(M_R_mean_stroke(subset_stroke)-M_R_steady)/Mg_fly/Lwing,'-ob','linewidth',1,'markersize',5,'markerfacecolor','b')
% plot(M_R_norm*mod_value_pitch(subset_pitch),(M_R_mean_pitch(subset_pitch)-M_R_steady)/Mg_fly/Lwing,'-oc','linewidth',1,'markersize',5,'markerfacecolor','c')
% plot(M_R_norm*mod_value_dev(subset_dev),(M_R_mean_dev(subset_dev)-M_R_steady)/Mg_fly/Lwing,'-om','linewidth',1,'markersize',5,'markerfacecolor','m')
% plot(M_R_norm*mod_value_allNOfreq(subset_allNOfreq),(M_R_mean_sum(subset_allNOfreq)-M_R_steady)/Mg_fly/Lwing,'-o','color',[.5 .5 .5],'linewidth',1,'markersize',5,'markerfacecolor',[.5 .5 .5])

subplot(2,2,4)
hold on
plot(M_R_norm*mod_value_allNOfreq(subset_allNOfreq),(M_R_mean_allNOfreq(subset_allNOfreq)-M_R_mean_allNOfreq(mod_value_allNOfreq==0))/Mg_fly/Lwing,'-or','linewidth',1,'markersize',5,'markerfacecolor','r')
% plot(M_R_norm*mod_value_freq(subset_freq),(M_R_mean_freq(subset_freq)-M_R_mean_freq(mod_value_freq==0))/Mg_fly/Lwing,'-og','linewidth',1,'markersize',5,'markerfacecolor','g')
plot(M_R_norm*mod_value_stroke(subset_stroke),(M_R_mean_stroke(subset_stroke)-M_R_mean_stroke(mod_value_stroke==0))/Mg_fly/Lwing,'-ob','linewidth',1,'markersize',5,'markerfacecolor','b')
plot(M_R_norm*mod_value_pitch(subset_pitch),(M_R_mean_pitch(subset_pitch)-M_R_mean_pitch(mod_value_pitch==0))/Mg_fly/Lwing,'-oc','linewidth',1,'markersize',5,'markerfacecolor','c')
plot(M_R_norm*mod_value_dev(subset_dev),(M_R_mean_dev(subset_dev)-M_R_mean_dev(mod_value_dev==0))/Mg_fly/Lwing,'-om','linewidth',1,'markersize',5,'markerfacecolor','m')
plot(M_R_norm*mod_value_allNOfreq(subset_allNOfreq),(M_R_mean_sum(subset_allNOfreq)-M_R_mean_sum(mod_value_allNOfreq==0))/Mg_fly/Lwing,'-o','color',[.5 .5 .5],'linewidth',1,'markersize',5,'markerfacecolor',[.5 .5 .5])

xlabel('M_R/(m+a)g fly','fontsize',10) 
ylabel('M_R/(m+a)g borf','fontsize',10)

% trendline
plot([x_min x_max],[x_min x_max],'-k')

% aYaw_norm = M_R_norm*mod_value_allNOfreq(subset_allNOfreq);
% M_R_norm = -(M_R_mean_allNOfreq(subset_allNOfreq)-M_R_mean_allNOfreq(mod_value_allNOfreq==0))/Mg_fly/Lwing;
% aM_R_fit = polyfit(aYaw_norm,M_R_norm,1);
% M_R_fit = polyval(aM_R_fit,[x_min x_max]);
% plot([x_min x_max],M_R_fit,'-k')
% 
% % save trend data
% save('InertiaMoment_trend_yaw.mat','aYaw_norm','M_R_norm','aM_R_fit','f_wb_steady','Lwing')

axis equal
% axis square
axis([x_min x_max y_min y_max])
set(gca,'XTick',[0:dx:x_max],'fontsize',8) 
set(gca,'YTick',[0:dy:y_max],'fontsize',8) 

%% save plot
mkdir('MSfigs_NObutter')
cd('MSfigs_NObutter')
saveas(gca, 'MSfig_Myaw_MaxisR.fig')
saveas(gca, 'MSfig_Myaw_MaxisR.png')
plot2svg(['MSfig_Myaw_MaxisR.svg'])
cd ..

%% save trenddata
% save('InertiaMoment_trends.mat','aRoll_norm','Mroll_norm','aMroll_fit','aYaw_norm','Myaw_norm','aMyaw_fit','aPitch_norm','Mpitch_norm','aMpitch_fit','f_wb_steady','Lwing')
