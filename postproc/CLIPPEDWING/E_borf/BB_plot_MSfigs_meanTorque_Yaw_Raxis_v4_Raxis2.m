% make borf plots mean F&M kinematics components & cross talk

clc
close all

%% Myaw

clear
load('MOD_norm_data.mat')

if exist(yawDB_mean) == 2
    
load(yawDB_mean)

x_min = -.02;
x_max = max(yawTorques);
dx = -x_min;
y_min = x_min;
y_max = x_max;
dy = dx;

dmod = Myaw_norm*(mod_value_allNOfreq(2)-mod_value_allNOfreq(1));
subset_allNOfreq = find(Myaw_norm*mod_value_allNOfreq>=x_min-dmod & Myaw_norm*mod_value_allNOfreq<=x_max+dmod );
% subset_freq = find(Myaw_norm*mod_value_freq>=x_min-dmod & Myaw_norm*mod_value_freq<=x_max+dmod );
subset_stroke = find(Myaw_norm*mod_value_stroke>=x_min-dmod & Myaw_norm*mod_value_stroke<=x_max+dmod );
subset_pitch = find(Myaw_norm*mod_value_pitch>=x_min-dmod & Myaw_norm*mod_value_pitch<=x_max+dmod );
subset_dev = find(Myaw_norm*mod_value_dev>=x_min-dmod & Myaw_norm*mod_value_dev<=x_max+dmod );

% kin components
subplot(2,2,1)
hold on
plot(Myaw_norm*mod_value_allNOfreq(subset_allNOfreq),(Mz_mean_sum(subset_allNOfreq)-Mz_mean_sum(mod_value_allNOfreq==0))/Mg_fly/Lwing,'-o','color',[.5 .5 .5],'linewidth',1,'markersize',5,'markerfacecolor',[.5 .5 .5])
plot(Myaw_norm*mod_value_allNOfreq(subset_allNOfreq),(Mz_mean_allNOfreq(subset_allNOfreq)-Mz_mean_allNOfreq(mod_value_allNOfreq==0))/Mg_fly/Lwing,'-or','linewidth',1,'markersize',5,'markerfacecolor','r')
% plot(Myaw_norm*mod_value_freq(subset_freq),(Mz_mean_freq(subset_freq)-Mz_mean_freq(mod_value_allNOfreq==0))/Mg_fly/Lwing,'-oc','linewidth',1,'markersize',5,'markerfacecolor','c')
plot(Myaw_norm*mod_value_stroke(subset_stroke),(Mz_mean_stroke(subset_stroke)-Mz_mean_stroke(mod_value_allNOfreq==0))/Mg_fly/Lwing,'-ob','linewidth',1,'markersize',5,'markerfacecolor','b')
plot(Myaw_norm*mod_value_pitch(subset_pitch),(Mz_mean_pitch(subset_pitch)-Mz_mean_pitch(mod_value_allNOfreq==0))/Mg_fly/Lwing,'-og','linewidth',1,'markersize',5,'markerfacecolor','g')
plot(Myaw_norm*mod_value_dev(subset_dev),(Mz_mean_dev(subset_dev)-Mz_mean_dev(mod_value_allNOfreq==0))/Mg_fly/Lwing,'-om','linewidth',1,'markersize',5,'markerfacecolor','m')

% trendline
plot([x_min x_max],[x_min x_max],'-k')

xlabel('Tyaw/(g+a)ml fly','fontsize',10) 
ylabel('Tyaw/(g+a)ml borf','fontsize',10)
axis equal
% axis square
axis([x_min x_max y_min y_max])
set(gca,'XTick',[x_min:dx:x_max],'fontsize',8) 
set(gca,'YTick',[y_min:dy:y_max],'fontsize',8) 

% cross-talk
subplot(2,2,3)
hold on
plot(Myaw_norm*mod_value_allNOfreq(subset_allNOfreq),(Mz_mean_allNOfreq(subset_allNOfreq)-Mz_mean_allNOfreq(mod_value_allNOfreq==0))/Mg_fly/Lwing,'-o','color','k','linewidth',1,'markersize',5,'markerfacecolor','k')
% plot(Myaw_norm*mod_value_allNOfreq(subset_allNOfreq),(M_R_mean_allNOfreq(subset_allNOfreq)-M_R_mean_allNOfreq(mod_value_allNOfreq==0))/Mg_fly/Lwing,'-o','color',[0 .5 1],'linewidth',1,'markersize',5,'markerfacecolor',[0 .5 1])
% plot(Myaw_norm*mod_value_allNOfreq(subset_allNOfreq),(M_L_mean_allNOfreq(subset_allNOfreq)-M_L_mean_allNOfreq(mod_value_allNOfreq==0))/Mg_fly/Lwing,'-o','color',[.5 .5 .5],'linewidth',1,'markersize',5,'markerfacecolor',[.5 .5 .5])
plot(Myaw_norm*mod_value_allNOfreq(subset_allNOfreq),(Mx_mean_allNOfreq(subset_allNOfreq)-Mx_mean_allNOfreq(mod_value_allNOfreq==0))/Mg_fly/Lwing,'-o','color','r','linewidth',1,'markersize',5,'markerfacecolor','r')
plot(Myaw_norm*mod_value_allNOfreq(subset_allNOfreq),(My_mean_allNOfreq(subset_allNOfreq)-My_mean_allNOfreq(mod_value_allNOfreq==0))/Mg_fly/Lwing,'-o','color','b','linewidth',1,'markersize',5,'markerfacecolor','b')

% trendline
plot([2*x_min 2*x_max],[2*x_min 2*x_max],'-k')

xlabel('Tyaw/(g+a)ml fly','fontsize',10) 
ylabel('Tyaw/(g+a)ml borf','fontsize',10)
axis equal
% axis square
axis([x_min x_max y_min y_max])
set(gca,'XTick',[x_min:dx:x_max],'fontsize',8) 
set(gca,'YTick',[y_min:dy:y_max],'fontsize',8) 

end

%% MaxisR

clear
load('MOD_norm_data.mat')

if exist(RaxisDB_mean) == 2
    
load(RaxisDB_mean)

x_min = -max(RaxisTorques)
x_max = 0;
dx = x_max/2;
y_min = x_min;
y_max = x_max;
dy = dx;

dmod = M_R_norm*(mod_value_allNOfreq(2)-mod_value_allNOfreq(1));
% subset_allNOfreq = find(M_R_norm*mod_value_allNOfreq>=x_min-dmod & M_R_norm*mod_value_allNOfreq<=x_max+dmod );
% % subset_freq = find(M_R_norm*mod_value_freq>=x_min-dmod & M_R_norm*mod_value_freq<=x_max+dmod );
% subset_stroke = find(M_R_norm*mod_value_stroke>=x_min-dmod & M_R_norm*mod_value_stroke<=x_max+dmod );
% subset_pitch = find(M_R_norm*mod_value_pitch>=x_min-dmod & M_R_norm*mod_value_pitch<=x_max+dmod );
% subset_dev = find(M_R_norm*mod_value_dev>=x_min-dmod & M_R_norm*mod_value_dev<=x_max+dmod );

% no pos
subset_allNOfreq = find(M_R_norm*mod_value_allNOfreq>=x_min-dmod & M_R_norm*mod_value_allNOfreq<=x_max );
% subset_freq = find(M_R_norm*mod_value_freq>=x_min-dmod & M_R_norm*mod_value_freq<=x_max );
subset_stroke = find(M_R_norm*mod_value_stroke>=x_min-dmod & M_R_norm*mod_value_stroke<=x_max );
subset_pitch = find(M_R_norm*mod_value_pitch>=x_min-dmod & M_R_norm*mod_value_pitch<=x_max );
subset_dev = find(M_R_norm*mod_value_dev>=x_min-dmod & M_R_norm*mod_value_dev<=x_max );

% % no neg
% subset_allNOfreq = find(M_R_norm*mod_value_allNOfreq>=x_min & M_R_norm*mod_value_allNOfreq<=x_max+dmod );
% % subset_freq = find(M_R_norm*mod_value_freq>=x_min & M_R_norm*mod_value_freq<=x_max+dmod );
% subset_stroke = find(M_R_norm*mod_value_stroke>=x_min & M_R_norm*mod_value_stroke<=x_max+dmod );
% subset_pitch = find(M_R_norm*mod_value_pitch>=x_min & M_R_norm*mod_value_pitch<=x_max+dmod );
% subset_dev = find(M_R_norm*mod_value_dev>=x_min & M_R_norm*mod_value_dev<=x_max+dmod );

% kin components
subplot(2,2,2)
hold on
plot(M_R_norm*mod_value_allNOfreq(subset_allNOfreq),(M_R_mean_sum(subset_allNOfreq)-M_R_mean_sum(mod_value_allNOfreq==0))/Mg_fly/Lwing,'-o','color',[.5 .5 .5],'linewidth',1,'markersize',5,'markerfacecolor',[.5 .5 .5])
plot(M_R_norm*mod_value_allNOfreq(subset_allNOfreq),(M_R_mean_allNOfreq(subset_allNOfreq)-M_R_mean_allNOfreq(mod_value_allNOfreq==0))/Mg_fly/Lwing,'-or','linewidth',1,'markersize',5,'markerfacecolor','r')
% plot(M_R_norm*mod_value_freq(subset_freq),(M_R_mean_freq(subset_freq)-M_R_mean_freq(mod_value_freq==0))/Mg_fly/Lwing,'-oc','linewidth',1,'markersize',5,'markerfacecolor','c')
plot(M_R_norm*mod_value_stroke(subset_stroke),(M_R_mean_stroke(subset_stroke)-M_R_mean_stroke(mod_value_stroke==0))/Mg_fly/Lwing,'-ob','linewidth',1,'markersize',5,'markerfacecolor','b')
plot(M_R_norm*mod_value_pitch(subset_pitch),(M_R_mean_pitch(subset_pitch)-M_R_mean_pitch(mod_value_pitch==0))/Mg_fly/Lwing,'-og','linewidth',1,'markersize',5,'markerfacecolor','g')
plot(M_R_norm*mod_value_dev(subset_dev),(M_R_mean_dev(subset_dev)-M_R_mean_dev(mod_value_dev==0))/Mg_fly/Lwing,'-om','linewidth',1,'markersize',5,'markerfacecolor','m')

% trendline
plot([x_min -x_min],[x_min -x_min],'-k')

xlabel('T_R/(g+a)ml fly','fontsize',10) 
ylabel('T_R/(g+a)ml borf','fontsize',10)
axis equal
% axis square
axis([x_min -x_min y_min -y_min])
set(gca,'XTick',[x_min:dx:-x_min],'fontsize',8) 
set(gca,'YTick',[y_min:dy:-y_min],'fontsize',8) 

% cross-talk
subplot(2,2,4)
hold on
plot(M_R_norm*mod_value_allNOfreq(subset_allNOfreq),(Mz_mean_allNOfreq(subset_allNOfreq)-Mz_mean_allNOfreq(mod_value_allNOfreq==0))/Mg_fly/Lwing,'-o','color','k','linewidth',1,'markersize',5,'markerfacecolor','k')
plot(M_R_norm*mod_value_allNOfreq(subset_allNOfreq),(M_R_mean_allNOfreq(subset_allNOfreq)-M_R_mean_allNOfreq(mod_value_allNOfreq==0))/Mg_fly/Lwing,'-o','color',[0 .5 1],'linewidth',1,'markersize',5,'markerfacecolor',[0 .5 1])
plot(M_R_norm*mod_value_allNOfreq(subset_allNOfreq),(M_L_mean_allNOfreq(subset_allNOfreq)-M_L_mean_allNOfreq(mod_value_allNOfreq==0))/Mg_fly/Lwing,'-o','color',[.5 .5 .5],'linewidth',1,'markersize',5,'markerfacecolor',[.5 .5 .5])
% plot(M_R_norm*mod_value_allNOfreq(subset_allNOfreq),(Mx_mean_allNOfreq(subset_allNOfreq)-Mx_mean_allNOfreq(mod_value_allNOfreq==0))/Mg_fly/Lwing,'-o','color','r','linewidth',1,'markersize',5,'markerfacecolor','r')
% plot(M_R_norm*mod_value_allNOfreq(subset_allNOfreq),(My_mean_allNOfreq(subset_allNOfreq)-My_mean_allNOfreq(mod_value_allNOfreq==0))/Mg_fly/Lwing,'-o','color','b','linewidth',1,'markersize',5,'markerfacecolor','b')

% trendline
plot([x_min -x_min],[x_min -x_min],'-k')

xlabel('T_R/(g+a)ml fly','fontsize',10) 
ylabel('T_R/(g+a)ml borf','fontsize',10)
axis equal
% axis square
axis([x_min -x_min y_min -y_min])
set(gca,'XTick',[x_min:dx:-x_min],'fontsize',8) 
set(gca,'YTick',[y_min:dy:-y_min],'fontsize',8) 

end
%% save plot
mkdir('MSfigs_NObutter')
cd('MSfigs_NObutter')
saveas(gca, 'MSfig_Myaw_MaxisR_kincomp_crosstalk.fig')
saveas(gca, 'MSfig_Myaw_MaxisR_kincomp_crosstalk.png')
plot2svg(['MSfig_Myaw_MaxisR_kincomp_crosstalk.svg'])
cd ..

