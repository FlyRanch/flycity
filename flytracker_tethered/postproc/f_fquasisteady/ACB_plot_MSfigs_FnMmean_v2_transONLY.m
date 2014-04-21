% make borf plots mean F&M

clear
clc
close all

%% load MOD&norm data
load('MOD_norm_data_QSmodel.mat')

%% Fenhance
load(FDB)

x_min = .825;
x_max = 1 + max(Fenhanses);
dx = max(Fenhanses)/2;
y_min = x_min;
y_max = x_max;
dy = dx;

dmod = mod_value_Fenhance(2)-mod_value_Fenhance(1);
subset = find(Fenhance_norm*mod_value_Fenhance>=x_min-1-dmod & Fenhance_norm*mod_value_Fenhance<=x_max-1+dmod );

Ftot_trans_norm_Fenhance_sum_mean = Ftot_trans_norm_Fenhance_freq_mean(mod_value_Fenhance==0) +...
    (Ftot_trans_norm_Fenhance_freq_mean - Ftot_trans_norm_Fenhance_freq_mean(mod_value_Fenhance==0)) +...
    (Ftot_trans_norm_Fenhance_stroke_mean - Ftot_trans_norm_Fenhance_stroke_mean(mod_value_Fenhance==0)) +...
    (Ftot_trans_norm_Fenhance_dev_mean - Ftot_trans_norm_Fenhance_dev_mean(mod_value_Fenhance==0)) +...
    (Ftot_trans_norm_Fenhance_rot_mean - Ftot_trans_norm_Fenhance_rot_mean(mod_value_Fenhance==0));

Ftot_trans_norm_Fenhance_sum_mean = Ftot_trans_norm_Fenhance_freq_mean(mod_value_Fenhance==0) +...
    (Ftot_trans_norm_Fenhance_freq_mean - Ftot_trans_norm_Fenhance_freq_mean(mod_value_Fenhance==0)) +...
    (Ftot_trans_norm_Fenhance_stroke_mean - Ftot_trans_norm_Fenhance_stroke_mean(mod_value_Fenhance==0)) +...
    (Ftot_trans_norm_Fenhance_dev_mean - Ftot_trans_norm_Fenhance_dev_mean(mod_value_Fenhance==0)) +...
    (Ftot_trans_norm_Fenhance_rot_mean - Ftot_trans_norm_Fenhance_rot_mean(mod_value_Fenhance==0));

% figure
subplot(3,3,1)
hold on

% plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Ftot_trans_norm_Fenhance_all_mean(subset),'-r','linewidth',1)
% plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Ftot_trans_norm_Fenhance_freq_mean(subset),'-g','linewidth',1)
% plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Ftot_trans_norm_Fenhance_stroke_mean(subset),'-b','linewidth',1)
% plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Ftot_trans_norm_Fenhance_dev_mean(subset),'-m','linewidth',1)
% plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Ftot_trans_norm_Fenhance_rot_mean(subset),'-c','linewidth',1)
% plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Ftot_trans_norm_Fenhance_sum_mean(subset),'-','color',[.5 .5 .5],'linewidth',1)
% 
plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Ftot_trans_norm_Fenhance_all_mean(subset),'-r','linewidth',1)
plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Ftot_trans_norm_Fenhance_freq_mean(subset),'-g','linewidth',1)
plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Ftot_trans_norm_Fenhance_stroke_mean(subset),'-b','linewidth',1)
plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Ftot_trans_norm_Fenhance_dev_mean(subset),'-m','linewidth',1)
plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Ftot_trans_norm_Fenhance_rot_mean(subset),'-c','linewidth',1)
plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Ftot_trans_norm_Fenhance_sum_mean(subset),'--','color',[.5 .5 .5],'linewidth',1)

xlabel('a+g/g','fontsize',10) 
ylabel('F/Mg','fontsize',10)


% trendline
plot([x_min x_max],[x_min x_max],'-k')

axis equal
axis([x_min x_max y_min y_max])
set(gca,'XTick',[1:dx:x_max],'fontsize',8) 
set(gca,'YTick',[1:dy:y_max],'fontsize',8) 


%% Mroll
load(rollDB)

x_min = -1;
x_max = max(rollaccels);
dx = max(rollaccels)/2;
y_min = -.01;
y_max = .05;
dy = .025;

dmod = mod_value_RollAccel(2)-mod_value_RollAccel(1);
subset = find(rollaccel_norm*mod_value_RollAccel>=x_min-dmod & rollaccel_norm*mod_value_RollAccel<=x_max+dmod );

Mx_trans_norm_RollAccel_sum_mean = ...
    (Mx_trans_norm_RollAccel_stroke_mean - Mx_trans_norm_RollAccel_stroke_mean(mod_value_RollAccel==0)) +...
    (Mx_trans_norm_RollAccel_dev_mean - Mx_trans_norm_RollAccel_dev_mean(mod_value_RollAccel==0)) +...
    (Mx_trans_norm_RollAccel_rot_mean - Mx_trans_norm_RollAccel_rot_mean(mod_value_RollAccel==0));

Mx_trans_norm_RollAccel_sum_mean = ...
    (Mx_trans_norm_RollAccel_stroke_mean - Mx_trans_norm_RollAccel_stroke_mean(mod_value_RollAccel==0)) +...
    (Mx_trans_norm_RollAccel_dev_mean - Mx_trans_norm_RollAccel_dev_mean(mod_value_RollAccel==0)) +...
    (Mx_trans_norm_RollAccel_rot_mean - Mx_trans_norm_RollAccel_rot_mean(mod_value_RollAccel==0));

subplot(3,3,4)
hold on
plot(rollaccel_norm*mod_value_RollAccel(subset),-(Mx_trans_norm_RollAccel_allNOfreq_mean(subset)-Mx_trans_norm_RollAccel_allNOfreq_mean(mod_value_RollAccel==0)),'-r','linewidth',1)
% plot(rollaccel_norm*mod_value_RollAccel(subset),-(Mx_trans_norm_RollAccel_freq_mean(subset)-Mx_trans_norm_RollAccel_freq_mean(mod_value_RollAccel==0)),'-g','linewidth',1)
plot(rollaccel_norm*mod_value_RollAccel(subset),-(Mx_trans_norm_RollAccel_stroke_mean(subset)-Mx_trans_norm_RollAccel_stroke_mean(mod_value_RollAccel==0)),'-b','linewidth',1)
plot(rollaccel_norm*mod_value_RollAccel(subset),-(Mx_trans_norm_RollAccel_dev_mean(subset)-Mx_trans_norm_RollAccel_dev_mean(mod_value_RollAccel==0)),'-m','linewidth',1)
plot(rollaccel_norm*mod_value_RollAccel(subset),-(Mx_trans_norm_RollAccel_rot_mean(subset)-Mx_trans_norm_RollAccel_rot_mean(mod_value_RollAccel==0)),'-c','linewidth',1)
plot(rollaccel_norm*mod_value_RollAccel(subset),-(Mx_trans_norm_RollAccel_sum_mean(subset)-Mx_trans_norm_RollAccel_sum_mean(mod_value_RollAccel==0)),'--','color',[.5 .5 .5],'linewidth',1)

xlabel('a roll norm','fontsize',10) 
ylabel('Mroll norm','fontsize',10)

% trendline
aRoll_norm = rollaccel_norm*mod_value_RollAccel(subset);
Mroll_norm = -(Mx_trans_norm_RollAccel_allNOfreq_mean(subset)-Mx_trans_norm_RollAccel_allNOfreq_mean(mod_value_RollAccel==0));
aMroll_fit = polyfit(aRoll_norm,Mroll_norm',1);
Mroll_fit = polyval(aMroll_fit,[x_min x_max]);
plot([x_min x_max],Mroll_fit,'--k')

% axis equal
axis square
axis([x_min x_max y_min y_max])
set(gca,'XTick',[0:dx:x_max],'fontsize',8) 
set(gca,'YTick',[0:dy:y_max],'fontsize',8) 

%% Myaw
load(yawDB)

x_min = -.5;
x_max = max(yawaccels);
dx = max(yawaccels)/2;
y_min = -.01;
y_max = .02;
dy = .01;

dmod = mod_value_YawAccel(2)-mod_value_YawAccel(1);
subset = find(yawaccel_norm*mod_value_YawAccel>=x_min-dmod & yawaccel_norm*mod_value_YawAccel<=x_max+dmod );

Mz_trans_norm_YawAccel_sum_mean = ...
    (Mz_trans_norm_YawAccel_stroke_mean - Mz_trans_norm_YawAccel_stroke_mean(mod_value_YawAccel==0)) +...
    (Mz_trans_norm_YawAccel_dev_mean - Mz_trans_norm_YawAccel_dev_mean(mod_value_YawAccel==0)) +...
    (Mz_trans_norm_YawAccel_rot_mean - Mz_trans_norm_YawAccel_rot_mean(mod_value_YawAccel==0));

Mz_trans_norm_YawAccel_sum_mean = ...
    (Mz_trans_norm_YawAccel_stroke_mean - Mz_trans_norm_YawAccel_stroke_mean(mod_value_YawAccel==0)) +...
    (Mz_trans_norm_YawAccel_dev_mean - Mz_trans_norm_YawAccel_dev_mean(mod_value_YawAccel==0)) +...
    (Mz_trans_norm_YawAccel_rot_mean - Mz_trans_norm_YawAccel_rot_mean(mod_value_YawAccel==0));

subplot(3,3,5)
hold on
plot(yawaccel_norm*mod_value_YawAccel(subset),-(Mz_trans_norm_YawAccel_allNOfreq_mean(subset)-Mz_trans_norm_YawAccel_allNOfreq_mean(mod_value_YawAccel==0)),'-r','linewidth',1)
% plot(yawaccel_norm*mod_value_YawAccel(subset),-(Mz_trans_norm_YawAccel_freq_mean(subset)-Mz_trans_norm_YawAccel_freq_mean(mod_value_YawAccel==0)),'-g','linewidth',1)
plot(yawaccel_norm*mod_value_YawAccel(subset),-(Mz_trans_norm_YawAccel_stroke_mean(subset)-Mz_trans_norm_YawAccel_stroke_mean(mod_value_YawAccel==0)),'-b','linewidth',1)
plot(yawaccel_norm*mod_value_YawAccel(subset),-(Mz_trans_norm_YawAccel_dev_mean(subset)-Mz_trans_norm_YawAccel_dev_mean(mod_value_YawAccel==0)),'-m','linewidth',1)
plot(yawaccel_norm*mod_value_YawAccel(subset),-(Mz_trans_norm_YawAccel_rot_mean(subset)-Mz_trans_norm_YawAccel_rot_mean(mod_value_YawAccel==0)),'-c','linewidth',1)
plot(yawaccel_norm*mod_value_YawAccel(subset),-(Mz_trans_norm_YawAccel_sum_mean(subset)-Mz_trans_norm_YawAccel_sum_mean(mod_value_YawAccel==0)),'--','color',[.5 .5 .5],'linewidth',1)

xlabel('a yaw norm','fontsize',10) 
ylabel('Myaw norm','fontsize',10)

% trendline
aYaw_norm = yawaccel_norm*mod_value_YawAccel(subset);
Myaw_norm = -(Mz_trans_norm_YawAccel_allNOfreq_mean(subset)-Mz_trans_norm_YawAccel_allNOfreq_mean(mod_value_YawAccel==0));
aMyaw_fit = polyfit(aYaw_norm,Myaw_norm',1);
Myaw_fit = polyval(aMyaw_fit,[x_min x_max]);
plot([x_min x_max],Myaw_fit,'--k')

% axis equal
axis square
axis([x_min x_max y_min y_max])
set(gca,'XTick',[0:dx:x_max],'fontsize',8) 
set(gca,'YTick',[0:dy:y_max],'fontsize',8) 

%% Mpitch
load(pitchDB)

x_min = -max(pitchaccels);
x_max = max(pitchaccels);
dx = max(pitchaccels);
y_min = -.08;
y_max = .08;
dy = .04;

dmod = mod_value_PitchAccel(2)-mod_value_PitchAccel(1);
subset = find(pitchaccel_norm*mod_value_PitchAccel>=x_min-dmod & pitchaccel_norm*mod_value_PitchAccel<=x_max+dmod );

My_trans_norm_PitchAccel_sum_mean = ...
    (My_trans_norm_PitchAccel_stroke_mean - My_trans_norm_PitchAccel_stroke_mean(mod_value_PitchAccel==0)) +...
    (My_trans_norm_PitchAccel_dev_mean - My_trans_norm_PitchAccel_dev_mean(mod_value_PitchAccel==0)) +...
    (My_trans_norm_PitchAccel_rot_mean - My_trans_norm_PitchAccel_rot_mean(mod_value_PitchAccel==0));

My_trans_norm_PitchAccel_sum_mean = ...
    (My_trans_norm_PitchAccel_stroke_mean - My_trans_norm_PitchAccel_stroke_mean(mod_value_PitchAccel==0)) +...
    (My_trans_norm_PitchAccel_dev_mean - My_trans_norm_PitchAccel_dev_mean(mod_value_PitchAccel==0)) +...
    (My_trans_norm_PitchAccel_rot_mean - My_trans_norm_PitchAccel_rot_mean(mod_value_PitchAccel==0));

subplot(3,3,7)
hold on
plot(pitchaccel_norm*mod_value_PitchAccel(subset),(My_trans_norm_PitchAccel_allNOfreq_mean(subset)-My_trans_norm_PitchAccel_allNOfreq_mean(mod_value_PitchAccel==0)),'-r','linewidth',1)
% plot(pitchaccel_norm*mod_value_PitchAccel(subset),(My_trans_norm_PitchAccel_freq_mean(subset)-My_trans_norm_PitchAccel_freq_mean(mod_value_PitchAccel==0)),'-g','linewidth',1)
plot(pitchaccel_norm*mod_value_PitchAccel(subset),(My_trans_norm_PitchAccel_stroke_mean(subset)-My_trans_norm_PitchAccel_stroke_mean(mod_value_PitchAccel==0)),'-b','linewidth',1)
plot(pitchaccel_norm*mod_value_PitchAccel(subset),(My_trans_norm_PitchAccel_dev_mean(subset)-My_trans_norm_PitchAccel_dev_mean(mod_value_PitchAccel==0)),'-m','linewidth',1)
plot(pitchaccel_norm*mod_value_PitchAccel(subset),(My_trans_norm_PitchAccel_rot_mean(subset)-My_trans_norm_PitchAccel_rot_mean(mod_value_PitchAccel==0)),'-c','linewidth',1)
plot(pitchaccel_norm*mod_value_PitchAccel(subset),(My_trans_norm_PitchAccel_sum_mean(subset)-My_trans_norm_PitchAccel_sum_mean(mod_value_PitchAccel==0)),'--','color',[.5 .5 .5],'linewidth',1)

xlabel('a pitch norm','fontsize',10) 
ylabel('Mpitch norm','fontsize',10)

% trendline
aPitch_norm = pitchaccel_norm*mod_value_PitchAccel(subset);
Mpitch_norm = (My_trans_norm_PitchAccel_allNOfreq_mean(subset)-My_trans_norm_PitchAccel_allNOfreq_mean(mod_value_PitchAccel==0));
aMpitch_fit = polyfit(aPitch_norm,Mpitch_norm',1);
Mpitch_fit = polyval(aMpitch_fit,[x_min x_max]);
plot([x_min x_max],Mpitch_fit,'--k')

% axis equal
axis square
axis([x_min x_max y_min y_max])
set(gca,'XTick',[0:dx:x_max],'fontsize',8) 
set(gca,'YTick',[y_min:dy:y_max],'fontsize',8) 

%% save plot
mkdir('MSfigs')
cd('MSfigs')
saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch_QSmodel_transONLY.fig')
saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch_QSmodel_transONLY.png')
plot2svg(['MSfig_Fenhance_Mroll_Myaw_Mpitch_QSmodel_transONLY.svg'])
cd ..

%% save trenddata
save('InertiaMoment_trends_QSmodel_transONLY.mat','aRoll_norm','Mroll_norm','aMroll_fit','aYaw_norm','Myaw_norm','aMyaw_fit','aPitch_norm','Mpitch_norm','aMpitch_fit','f_wb_steady')
