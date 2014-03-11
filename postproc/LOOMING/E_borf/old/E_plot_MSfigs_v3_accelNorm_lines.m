% make borf plots mean F&M

clear
clc
close all

%% load norm data
load('norm_data.mat')

% FDB = 'borf_db_Fenhance_NOcali_means.mat'
% rollDB = 'borf_db_F_roll_LR_NOcali_means.mat'
% pitchDB = 'borf_db_PitchAccel_NOcali_means.mat'

FDB = 'borf_db_Fenhance_cali_means.mat'
rollDB = 'borf_db_F_roll_LR_cali_means.mat'
pitchDB = 'borf_db_PitchAccel_cali_means.mat'

%% Fenhance
load(FDB)

figure
subplot(3,3,1)
hold on

plot(Fenhance_norm*mod_value_all+1,F_mean_all/Mg_fly,'-r','linewidth',1)
plot(Fenhance_norm*mod_value_freq+1,F_mean_freq/Mg_fly,'-g','linewidth',1)
plot(Fenhance_norm*mod_value_stroke+1,F_mean_stroke/Mg_fly,'-b','linewidth',1)
plot(Fenhance_norm*mod_value_pitch+1,F_mean_pitch/Mg_fly,'-c','linewidth',1)
plot(Fenhance_norm*mod_value_dev+1,F_mean_dev/Mg_fly,'-m','linewidth',1)
% plot(Fenhance_norm*mod_value_all(1:MODmin),F_mean_sum/Mg_fly,'-k','linewidth',1)

xlabel('a+g/g','fontsize',10) 
ylabel('F/Mg','fontsize',10)

x_min = .75;
x_max = 2;
dx = .5;;
y_min = x_min;
y_max = x_max;
dy = dx;

% trendline
plot([x_min x_max],[x_min x_max],'--k')

axis equal
axis([x_min x_max y_min y_max])
set(gca,'XTick',[1:dx:x_max],'fontsize',8) 
set(gca,'YTick',[1:dy:y_max],'fontsize',8) 


%% Mroll
load(rollDB)

subplot(3,3,4)
hold on
plot(rollaccel_norm*mod_value_allNOfreq,(Mx_mean_allNOfreq-Mx_steady)/f_wb_steady^2/Lwing^5,'-r','linewidth',1)
% plot(rollaccel_norm*mod_value_freq,(Mx_mean_freq-Mx_steady)/f_wb_steady^2/Lwing^5,'-g','linewidth',1)
plot(rollaccel_norm*mod_value_stroke,(Mx_mean_stroke-Mx_steady)/f_wb_steady^2/Lwing^5,'-b','linewidth',1)
plot(rollaccel_norm*mod_value_pitch,(Mx_mean_pitch-Mx_steady)/f_wb_steady^2/Lwing^5,'-c','linewidth',1)
plot(rollaccel_norm*mod_value_dev,(Mx_mean_dev-Mx_steady)/f_wb_steady^2/Lwing^5,'-m','linewidth',1)
% plot(rollaccel_norm*mod_value_allNOfreq(1:MODmin),(Mx_mean_sum-Mx_steady)/f_wb_steady^2/Lwing^5,'-k','linewidth',1)

xlabel('a roll norm','fontsize',10) 
ylabel('Mroll norm','fontsize',10)

x_min = -1;
x_max = 5;
dx = 2.5;
y_min = -.25;
y_max = 1;
dy = .5;

% trendline
aRoll_norm = rollaccel_norm*mod_value_allNOfreq;
Mroll_norm = (Mx_mean_allNOfreq-Mx_steady)/f_wb_steady^2/Lwing^5;
aMroll_fit = polyfit(aRoll_norm,Mroll_norm,1);
Mroll_fit = polyval(aMroll_fit,[x_min x_max]);
plot([x_min x_max],Mroll_fit,'--k')

% axis equal
axis square
axis([x_min x_max y_min y_max])
set(gca,'XTick',[0:dx:x_max],'fontsize',8) 
set(gca,'YTick',[0:dy:y_max],'fontsize',8) 

%% Mpitch
load(pitchDB)

subplot(3,3,7)
hold on
plot(pitchaccel_norm*mod_value_allNOfreq,(My_mean_allNOfreq-My_steady)/f_wb_steady^2/Lwing^5,'-r','linewidth',1)
% plot(pitchaccel_norm*mod_value_freq,(My_mean_freq-My_steady)/f_wb_steady^2/Lwing^5,'-g','linewidth',1)
plot(pitchaccel_norm*mod_value_stroke,(My_mean_stroke-My_steady)/f_wb_steady^2/Lwing^5,'-b','linewidth',1)
plot(pitchaccel_norm*mod_value_pitch,(My_mean_pitch-My_steady)/f_wb_steady^2/Lwing^5,'-c','linewidth',1)
plot(pitchaccel_norm*mod_value_dev,(My_mean_dev-My_steady)/f_wb_steady^2/Lwing^5,'-m','linewidth',1)
% plot(pitchaccel_norm*mod_value_allNOfreq(1:MODmin),(My_mean_sum-My_steady),'-k','linewidth',1)

xlabel('a pitch norm','fontsize',10) 
ylabel('Mpitch norm','fontsize',10)

x_min = -3;
x_max = 3;
dx = 1;
y_min = -.4;
y_max = .4;
dy = .2;

% trendline
aPitch_norm = pitchaccel_norm*mod_value_allNOfreq;
Mpitch_norm = (My_mean_allNOfreq-My_steady)/f_wb_steady^2/Lwing^5;
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
saveas(gca, 'MSfig_Fenhance_Mroll_Mpitch.fig')
saveas(gca, 'MSfig_Fenhance_Mroll_Mpitch.png')
plot2svg(['MSfig_Fenhance_Mroll_Mpitch.svg'])
    
%% save trenddata
save('InertiaMoment_trends.mat','aRoll_norm','Mroll_norm','aMroll_fit','aPitch_norm','Mpitch_norm','aMpitch_fit','f_wb_steady','Lwing')
