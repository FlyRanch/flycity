% make borf plots mean F&M

clc
close all
figure

%% Fenhance
clear
% load('MOD_norm_data_QSmodel.mat')
load('norm_data_torque.mat')
load(FDB)

Fdeficit_on = 1;
Fdeficit_on = 0;

rot_on = 1;
% rot_on = 0;

Fdeficit = (Fenhance_norm*mod_value_Fenhance+1)./Ftot_trans_norm_Fenhance_all_mean';
Fdeficit_steady = Fdeficit(mod_value_Fenhance==0);
Fdeficit_mean = nanmean(Fdeficit);
% Fdeficit_compensation = Fdeficit_mean;
Fdeficit_compensation = Fdeficit_steady;

save('QSmodel_settings.mat','Fdeficit','Fdeficit_compensation','Fdeficit_on','rot_on')
if Fdeficit_on == 0
    Fdeficit_compensation = 1;
end

% plots
x_min = 1 + floor(10*min(mod_value_Fenhance)*Fenhance_norm)/10;
x_max = 1 + ceil(10*max(mod_value_Fenhance)*Fenhance_norm)/10;
dx = (x_max-1)/2;
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

Ftot_transNrot_norm_Fenhance_sum_mean = Ftot_transNrot_norm_Fenhance_freq_mean(mod_value_Fenhance==0) +...
    (Ftot_transNrot_norm_Fenhance_freq_mean - Ftot_transNrot_norm_Fenhance_freq_mean(mod_value_Fenhance==0)) +...
    (Ftot_transNrot_norm_Fenhance_stroke_mean - Ftot_transNrot_norm_Fenhance_stroke_mean(mod_value_Fenhance==0)) +...
    (Ftot_transNrot_norm_Fenhance_dev_mean - Ftot_transNrot_norm_Fenhance_dev_mean(mod_value_Fenhance==0)) +...
    (Ftot_transNrot_norm_Fenhance_rot_mean - Ftot_transNrot_norm_Fenhance_rot_mean(mod_value_Fenhance==0));

% kinematics components
figure(1)
subplot(2,3,1)
hold on

if rot_on == 0
    plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Fdeficit_compensation*Ftot_trans_norm_Fenhance_all_mean(subset),'-r.','linewidth',1)
    plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Fdeficit_compensation*Ftot_trans_norm_Fenhance_freq_mean(subset),'-g.','linewidth',1)
    plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Fdeficit_compensation*Ftot_trans_norm_Fenhance_stroke_mean(subset),'-b.','linewidth',1)
    plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Fdeficit_compensation*Ftot_trans_norm_Fenhance_dev_mean(subset),'-m.','linewidth',1)
    plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Fdeficit_compensation*Ftot_trans_norm_Fenhance_rot_mean(subset),'-c.','linewidth',1)
    plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Fdeficit_compensation*Ftot_trans_norm_Fenhance_sum_mean(subset),'-','color',[.5 .5 .5],'linewidth',1)
else
    plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Fdeficit_compensation*Ftot_transNrot_norm_Fenhance_all_mean(subset),'-r.','linewidth',1)
    plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Fdeficit_compensation*Ftot_transNrot_norm_Fenhance_freq_mean(subset),'-g.','linewidth',1)
    plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Fdeficit_compensation*Ftot_transNrot_norm_Fenhance_stroke_mean(subset),'-b.','linewidth',1)
    plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Fdeficit_compensation*Ftot_transNrot_norm_Fenhance_dev_mean(subset),'-m.','linewidth',1)
    plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Fdeficit_compensation*Ftot_transNrot_norm_Fenhance_rot_mean(subset),'-c.','linewidth',1)
    plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Fdeficit_compensation*Ftot_transNrot_norm_Fenhance_sum_mean(subset),'--','color',[.5 .5 .5],'linewidth',1)
end

xlabel('a+g/g fly','fontsize',10) 
ylabel('F/Mg qsmodel','fontsize',10)

% trendline
plot([x_min x_max],[x_min x_max],'-k')

axis equal
axis([x_min x_max y_min y_max])
set(gca,'XTick',[1:dx:x_max],'fontsize',8) 
set(gca,'YTick',[1:dy:y_max],'fontsize',8) 

% cross-talk
figure(2)
subplot(2,3,1)
hold on

if rot_on == 0
    plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Fdeficit_compensation*Fx_trans_norm_Fenhance_all_mean(subset),'-r.','linewidth',1)
    plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Fdeficit_compensation*Fy_trans_norm_Fenhance_all_mean(subset),'-g.','linewidth',1)
    plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Fdeficit_compensation*-Fz_trans_norm_Fenhance_all_mean(subset),'-b.','linewidth',1)
else
    plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Fdeficit_compensation*Fx_transNrot_norm_Fenhance_all_mean(subset),'-r.','linewidth',1)
    plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Fdeficit_compensation*Fy_transNrot_norm_Fenhance_all_mean(subset),'-g.','linewidth',1)
    plot(Fenhance_norm*mod_value_Fenhance(subset)+1,Fdeficit_compensation*-Fz_transNrot_norm_Fenhance_all_mean(subset),'-b.','linewidth',1)
end

xlabel('a+g/g fly','fontsize',10) 
ylabel('F/Mg qsmodel','fontsize',10)

% trendline
% plot([x_min x_max],[x_min x_max],'-k')

axis equal
axis([x_min x_max y_min y_max])
axis([x_min x_max 0 1.5])
set(gca,'XTick',[1:.5:x_max],'fontsize',8) 
set(gca,'YTick',[0:.5:y_max],'fontsize',8) 

%% Mroll
clear
% load('MOD_norm_data_QSmodel.mat')
load('norm_data_torque.mat')
load(rollDB)

load('QSmodel_settings.mat')
if Fdeficit_on == 0
    Fdeficit_compensation = 1;
end

x_min = floor(100*min(mod_value_RollTorque)*Mroll_norm)/100;
x_max = ceil(100*max(mod_value_RollTorque)*Mroll_norm)/100;
dx = (x_max)/2;
% x_min = -max(rollTorques)/2;
% x_max = max(rollTorques);
% dx = max(rollTorques)/2;
y_min = x_min;
y_max = x_max;
dy = dx;

dmod = mod_value_RollTorque(2)-mod_value_RollTorque(1);
subset = find(Mroll_norm*mod_value_RollTorque>=x_min-dmod & Mroll_norm*mod_value_RollTorque<=x_max+dmod );

Mx_trans_norm_RollTorque_sum_mean = ...
    ((Mx_trans_norm_RollTorque_stroke_mean - Mx_trans_norm_RollTorque_stroke_mean(mod_value_RollTorque==0)) +...
    (Mx_trans_norm_RollTorque_dev_mean - Mx_trans_norm_RollTorque_dev_mean(mod_value_RollTorque==0)) +...
    (Mx_trans_norm_RollTorque_rot_mean - Mx_trans_norm_RollTorque_rot_mean(mod_value_RollTorque==0)));

Mx_transNrot_norm_RollTorque_sum_mean = ...
    ((Mx_transNrot_norm_RollTorque_stroke_mean - Mx_transNrot_norm_RollTorque_stroke_mean(mod_value_RollTorque==0)) +...
    (Mx_transNrot_norm_RollTorque_dev_mean - Mx_transNrot_norm_RollTorque_dev_mean(mod_value_RollTorque==0)) +...
    (Mx_transNrot_norm_RollTorque_rot_mean - Mx_transNrot_norm_RollTorque_rot_mean(mod_value_RollTorque==0)));

% R&L axes
varname=dir('flyVars*')
varname = varname.name;
load(varname)

M_R_trans_norm_RollTorque_allNOfreq_mean =  Mx_trans_norm_RollTorque_allNOfreq_mean*cosd(rot_axis_angle) + My_trans_norm_RollTorque_allNOfreq_mean*sind(rot_axis_angle);
M_L_trans_norm_RollTorque_allNOfreq_mean = -Mx_trans_norm_RollTorque_allNOfreq_mean*sind(rot_axis_angle) + My_trans_norm_RollTorque_allNOfreq_mean*cosd(rot_axis_angle);

M_R_transNrot_norm_RollTorque_allNOfreq_mean =  Mx_transNrot_norm_RollTorque_allNOfreq_mean*cosd(rot_axis_angle) + My_transNrot_norm_RollTorque_allNOfreq_mean*sind(rot_axis_angle);
M_L_transNrot_norm_RollTorque_allNOfreq_mean = -Mx_transNrot_norm_RollTorque_allNOfreq_mean*sind(rot_axis_angle) + My_transNrot_norm_RollTorque_allNOfreq_mean*cosd(rot_axis_angle);

% kinematics components
figure(1)
subplot(2,3,2)
hold on

if rot_on == 0
    plot(Mroll_norm*mod_value_RollTorque(subset),-Fdeficit_compensation*(Mx_trans_norm_RollTorque_allNOfreq_mean(subset)-Mx_trans_norm_RollTorque_allNOfreq_mean(mod_value_RollTorque==0)),'-r.','linewidth',1)
    % plot(Mroll_norm*mod_value_RollTorque(subset),-Fdeficit_compensation*(Mx_trans_norm_RollTorque_freq_mean(subset)-Mx_trans_norm_RollTorque_freq_mean(mod_value_RollTorque==0)),'-g.','linewidth',1)
    plot(Mroll_norm*mod_value_RollTorque(subset),-Fdeficit_compensation*(Mx_trans_norm_RollTorque_stroke_mean(subset)-Mx_trans_norm_RollTorque_stroke_mean(mod_value_RollTorque==0)),'-b.','linewidth',1)
    plot(Mroll_norm*mod_value_RollTorque(subset),-Fdeficit_compensation*(Mx_trans_norm_RollTorque_dev_mean(subset)-Mx_trans_norm_RollTorque_dev_mean(mod_value_RollTorque==0)),'-m.','linewidth',1)
    plot(Mroll_norm*mod_value_RollTorque(subset),-Fdeficit_compensation*(Mx_trans_norm_RollTorque_rot_mean(subset)-Mx_trans_norm_RollTorque_rot_mean(mod_value_RollTorque==0)),'-c.','linewidth',1)
    plot(Mroll_norm*mod_value_RollTorque(subset),-Fdeficit_compensation*(Mx_trans_norm_RollTorque_sum_mean(subset)-Mx_trans_norm_RollTorque_sum_mean(mod_value_RollTorque==0)),'--','color',[.5 .5 .5],'linewidth',1)
else
    plot(Mroll_norm*mod_value_RollTorque(subset),-Fdeficit_compensation*(Mx_transNrot_norm_RollTorque_allNOfreq_mean(subset)-Mx_transNrot_norm_RollTorque_allNOfreq_mean(mod_value_RollTorque==0)),'-r.','linewidth',1)
    % plot(Mroll_norm*mod_value_RollTorque(subset),-Fdeficit_compensation*(Mx_transNrot_norm_RollTorque_freq_mean(subset)-Mx_transNrot_norm_RollTorque_freq_mean(mod_value_RollTorque==0)),'-g.','linewidth',1)
    plot(Mroll_norm*mod_value_RollTorque(subset),-Fdeficit_compensation*(Mx_transNrot_norm_RollTorque_stroke_mean(subset)-Mx_transNrot_norm_RollTorque_stroke_mean(mod_value_RollTorque==0)),'-b.','linewidth',1)
    plot(Mroll_norm*mod_value_RollTorque(subset),-Fdeficit_compensation*(Mx_transNrot_norm_RollTorque_dev_mean(subset)-Mx_transNrot_norm_RollTorque_dev_mean(mod_value_RollTorque==0)),'-m.','linewidth',1)
    plot(Mroll_norm*mod_value_RollTorque(subset),-Fdeficit_compensation*(Mx_transNrot_norm_RollTorque_rot_mean(subset)-Mx_transNrot_norm_RollTorque_rot_mean(mod_value_RollTorque==0)),'-c.','linewidth',1)
    plot(Mroll_norm*mod_value_RollTorque(subset),-Fdeficit_compensation*(Mx_transNrot_norm_RollTorque_sum_mean(subset)-Mx_transNrot_norm_RollTorque_sum_mean(mod_value_RollTorque==0)),'--','color',[.5 .5 .5],'linewidth',1)
end

xlabel('Mroll/(a+g)ml fly','fontsize',10) 
ylabel('Mroll/mgl qsmodel','fontsize',10)

% trendline
plot([x_min x_max],[x_min x_max],'-k')
% aRoll_norm = Mroll_norm*mod_value_RollTorque(subset);
% Mroll_norm = -(Mx_transNrot_norm_RollTorque_allNOfreq_mean(subset)-Mx_transNrot_norm_RollTorque_allNOfreq_mean(mod_value_RollTorque==0));
% aMroll_fit = polyfit(aRoll_norm,Mroll_norm',1);
% Mroll_fit = polyval(aMroll_fit,[x_min x_max]);
% plot([x_min x_max],Mroll_fit,'--k')

% axis equal
axis square
axis([x_min x_max y_min y_max])
set(gca,'XTick',[0:dx:x_max],'fontsize',8) 
set(gca,'YTick',[0:dy:y_max],'fontsize',8) 

% cross-talk
figure(2)
subplot(2,3,2)
hold on

if rot_on == 0
    plot(Mroll_norm*mod_value_RollTorque(subset),-Fdeficit_compensation*(Mx_trans_norm_RollTorque_allNOfreq_mean(subset)-Mx_trans_norm_RollTorque_allNOfreq_mean(mod_value_RollTorque==0)),'-r.','linewidth',1)
    plot(Mroll_norm*mod_value_RollTorque(subset), Fdeficit_compensation*(My_trans_norm_RollTorque_allNOfreq_mean(subset)-My_trans_norm_RollTorque_allNOfreq_mean(mod_value_RollTorque==0)),'-g.','linewidth',1)
    plot(Mroll_norm*mod_value_RollTorque(subset), Fdeficit_compensation*(Mz_trans_norm_RollTorque_allNOfreq_mean(subset)-Mz_trans_norm_RollTorque_allNOfreq_mean(mod_value_RollTorque==0)),'-b.','linewidth',1)
%     plot(Mroll_norm*mod_value_RollTorque(subset), Fdeficit_compensation*(M_R_trans_norm_RollTorque_allNOfreq_mean(subset)-M_R_trans_norm_RollTorque_allNOfreq_mean(mod_value_RollTorque==0)),'-c.','linewidth',1)
%     plot(Mroll_norm*mod_value_RollTorque(subset), Fdeficit_compensation*(M_L_trans_norm_RollTorque_allNOfreq_mean(subset)-M_L_trans_norm_RollTorque_allNOfreq_mean(mod_value_RollTorque==0)),'-m.','linewidth',1)
else
    plot(Mroll_norm*mod_value_RollTorque(subset),-Fdeficit_compensation*(Mx_transNrot_norm_RollTorque_allNOfreq_mean(subset)-Mx_transNrot_norm_RollTorque_allNOfreq_mean(mod_value_RollTorque==0)),'-r.','linewidth',1)
    plot(Mroll_norm*mod_value_RollTorque(subset), Fdeficit_compensation*(My_transNrot_norm_RollTorque_allNOfreq_mean(subset)-My_transNrot_norm_RollTorque_allNOfreq_mean(mod_value_RollTorque==0)),'-g.','linewidth',1)
    plot(Mroll_norm*mod_value_RollTorque(subset), Fdeficit_compensation*(Mz_transNrot_norm_RollTorque_allNOfreq_mean(subset)-Mz_transNrot_norm_RollTorque_allNOfreq_mean(mod_value_RollTorque==0)),'-b.','linewidth',1)
%     plot(Mroll_norm*mod_value_RollTorque(subset), Fdeficit_compensation*(M_R_transNrot_norm_RollTorque_allNOfreq_mean(subset)-M_R_transNrot_norm_RollTorque_allNOfreq_mean(mod_value_RollTorque==0)),'-c.','linewidth',1)
%     plot(Mroll_norm*mod_value_RollTorque(subset), Fdeficit_compensation*(M_L_transNrot_norm_RollTorque_allNOfreq_mean(subset)-M_L_transNrot_norm_RollTorque_allNOfreq_mean(mod_value_RollTorque==0)),'-m.','linewidth',1)
end

xlabel('Mroll/(a+g)ml fly','fontsize',10) 
ylabel('M/mgl qsmodel','fontsize',10)

% trendline
plot([x_min x_max],[x_min x_max],'-k')
% aRoll_norm = Mroll_norm*mod_value_RollTorque(subset);
% Mroll_norm = -(Mx_transNrot_norm_RollTorque_allNOfreq_mean(subset)-Mx_transNrot_norm_RollTorque_allNOfreq_mean(mod_value_RollTorque==0));
% aMroll_fit = polyfit(aRoll_norm,Mroll_norm',1);
% Mroll_fit = polyval(aMroll_fit,[x_min x_max]);
% plot([x_min x_max],Mroll_fit,'--k')

% axis equal
axis square
axis([x_min x_max y_min y_max])
set(gca,'XTick',[0:dx:x_max],'fontsize',8) 
set(gca,'YTick',[0:dy:y_max],'fontsize',8) 

%% Myaw
clear
% load('MOD_norm_data_QSmodel.mat')
load('norm_data_torque.mat')
load(yawDB)

load('QSmodel_settings.mat')
if Fdeficit_on == 0
    Fdeficit_compensation = 1;
end

x_min = floor(100*min(mod_value_YawTorque)*Myaw_norm)/100;
x_max = ceil(100*max(mod_value_YawTorque)*Myaw_norm)/100;
dx = (x_max)/2;
% x_min = -max(yawTorques)/2;
% x_max = max(yawTorques);
% dx = max(yawTorques)/2;
y_min = x_min;
y_max = x_max;
dy = dx;

dmod = mod_value_YawTorque(2)-mod_value_YawTorque(1);
subset = find(Myaw_norm*mod_value_YawTorque>=x_min-dmod & Myaw_norm*mod_value_YawTorque<=x_max+dmod );

Mz_trans_norm_YawTorque_sum_mean = ...
    (Mz_trans_norm_YawTorque_stroke_mean - Mz_trans_norm_YawTorque_stroke_mean(mod_value_YawTorque==0)) +...
    (Mz_trans_norm_YawTorque_dev_mean - Mz_trans_norm_YawTorque_dev_mean(mod_value_YawTorque==0)) +...
    (Mz_trans_norm_YawTorque_rot_mean - Mz_trans_norm_YawTorque_rot_mean(mod_value_YawTorque==0));

Mz_transNrot_norm_YawTorque_sum_mean = ...
    (Mz_transNrot_norm_YawTorque_stroke_mean - Mz_transNrot_norm_YawTorque_stroke_mean(mod_value_YawTorque==0)) +...
    (Mz_transNrot_norm_YawTorque_dev_mean - Mz_transNrot_norm_YawTorque_dev_mean(mod_value_YawTorque==0)) +...
    (Mz_transNrot_norm_YawTorque_rot_mean - Mz_transNrot_norm_YawTorque_rot_mean(mod_value_YawTorque==0));

% R&L axes
varname=dir('flyVars*')
varname = varname.name;
load(varname)

M_R_trans_norm_YawTorque_allNOfreq_mean =  Mx_trans_norm_YawTorque_allNOfreq_mean*cosd(rot_axis_angle) + My_trans_norm_YawTorque_allNOfreq_mean*sind(rot_axis_angle);
M_L_trans_norm_YawTorque_allNOfreq_mean = -Mx_trans_norm_YawTorque_allNOfreq_mean*sind(rot_axis_angle) + My_trans_norm_YawTorque_allNOfreq_mean*cosd(rot_axis_angle);

M_R_transNrot_norm_YawTorque_allNOfreq_mean =  Mx_transNrot_norm_YawTorque_allNOfreq_mean*cosd(rot_axis_angle) + My_transNrot_norm_YawTorque_allNOfreq_mean*sind(rot_axis_angle);
M_L_transNrot_norm_YawTorque_allNOfreq_mean = -Mx_transNrot_norm_YawTorque_allNOfreq_mean*sind(rot_axis_angle) + My_transNrot_norm_YawTorque_allNOfreq_mean*cosd(rot_axis_angle);

% kinematics components
figure(1)
subplot(2,3,4)
hold on

if rot_on == 0
    plot(Myaw_norm*mod_value_YawTorque(subset),Fdeficit_compensation*(Mz_trans_norm_YawTorque_allNOfreq_mean(subset)-Mz_trans_norm_YawTorque_allNOfreq_mean(mod_value_YawTorque==0)),'-r.','linewidth',1)
    % plot(Myaw_norm*mod_value_YawTorque(subset),Fdeficit_compensation*(Mz_trans_norm_YawTorque_freq_mean(subset)-Mz_trans_norm_YawTorque_freq_mean(mod_value_YawTorque==0)),'-g.','linewidth',1)
    plot(Myaw_norm*mod_value_YawTorque(subset),Fdeficit_compensation*(Mz_trans_norm_YawTorque_stroke_mean(subset)-Mz_trans_norm_YawTorque_stroke_mean(mod_value_YawTorque==0)),'-b.','linewidth',1)
    plot(Myaw_norm*mod_value_YawTorque(subset),Fdeficit_compensation*(Mz_trans_norm_YawTorque_dev_mean(subset)-Mz_trans_norm_YawTorque_dev_mean(mod_value_YawTorque==0)),'-m.','linewidth',1)
    plot(Myaw_norm*mod_value_YawTorque(subset),Fdeficit_compensation*(Mz_trans_norm_YawTorque_rot_mean(subset)-Mz_trans_norm_YawTorque_rot_mean(mod_value_YawTorque==0)),'-c.','linewidth',1)
    plot(Myaw_norm*mod_value_YawTorque(subset),Fdeficit_compensation*(Mz_trans_norm_YawTorque_sum_mean(subset)-Mz_trans_norm_YawTorque_sum_mean(mod_value_YawTorque==0)),'--','color',[.5 .5 .5],'linewidth',1)
else
    plot(Myaw_norm*mod_value_YawTorque(subset),Fdeficit_compensation*(Mz_transNrot_norm_YawTorque_allNOfreq_mean(subset)-Mz_transNrot_norm_YawTorque_allNOfreq_mean(mod_value_YawTorque==0)),'-r.','linewidth',1)
    % plot(Myaw_norm*mod_value_YawTorque(subset),Fdeficit_compensation*(Mz_transNrot_norm_YawTorque_freq_mean(subset)-Mz_transNrot_norm_YawTorque_freq_mean(mod_value_YawTorque==0)),'-g.','linewidth',1)
    plot(Myaw_norm*mod_value_YawTorque(subset),Fdeficit_compensation*(Mz_transNrot_norm_YawTorque_stroke_mean(subset)-Mz_transNrot_norm_YawTorque_stroke_mean(mod_value_YawTorque==0)),'-b.','linewidth',1)
    plot(Myaw_norm*mod_value_YawTorque(subset),Fdeficit_compensation*(Mz_transNrot_norm_YawTorque_dev_mean(subset)-Mz_transNrot_norm_YawTorque_dev_mean(mod_value_YawTorque==0)),'-m.','linewidth',1)
    plot(Myaw_norm*mod_value_YawTorque(subset),Fdeficit_compensation*(Mz_transNrot_norm_YawTorque_rot_mean(subset)-Mz_transNrot_norm_YawTorque_rot_mean(mod_value_YawTorque==0)),'-c.','linewidth',1)
    plot(Myaw_norm*mod_value_YawTorque(subset),Fdeficit_compensation*(Mz_transNrot_norm_YawTorque_sum_mean(subset)-Mz_transNrot_norm_YawTorque_sum_mean(mod_value_YawTorque==0)),'--','color',[.5 .5 .5],'linewidth',1)
end

xlabel('Myaw/(a+g)ml fly','fontsize',10) 
ylabel('Myaw/mgl qsmodel','fontsize',10)

% trendline
plot([x_min x_max],[x_min x_max],'-k')
% aYaw_norm = Myaw_norm*mod_value_YawTorque(subset);
% Myaw_norm = -(Mz_transNrot_norm_YawTorque_allNOfreq_mean(subset)-Mz_transNrot_norm_YawTorque_allNOfreq_mean(mod_value_YawTorque==0));
% aMyaw_fit = polyfit(aYaw_norm,Myaw_norm',1);
% Myaw_fit = polyval(aMyaw_fit,[x_min x_max]);
% plot([x_min x_max],Myaw_fit,'--k')

% axis equal
axis square
axis([x_min x_max y_min y_max])
set(gca,'XTick',[0:dx:x_max],'fontsize',8) 
set(gca,'YTick',[0:dy:y_max],'fontsize',8) 

% cross-talk roll&pitch
figure(2)
subplot(2,3,4)
hold on

if rot_on == 0
    plot(Myaw_norm*mod_value_YawTorque(subset),-Fdeficit_compensation*(Mx_trans_norm_YawTorque_allNOfreq_mean(subset)-Mx_trans_norm_YawTorque_allNOfreq_mean(mod_value_YawTorque==0)),'-r.','linewidth',1)
    plot(Myaw_norm*mod_value_YawTorque(subset), Fdeficit_compensation*(My_trans_norm_YawTorque_allNOfreq_mean(subset)-My_trans_norm_YawTorque_allNOfreq_mean(mod_value_YawTorque==0)),'-g.','linewidth',1)
    plot(Myaw_norm*mod_value_YawTorque(subset), Fdeficit_compensation*(Mz_trans_norm_YawTorque_allNOfreq_mean(subset)-Mz_trans_norm_YawTorque_allNOfreq_mean(mod_value_YawTorque==0)),'-b.','linewidth',1)
%     plot(Myaw_norm*mod_value_YawTorque(subset), Fdeficit_compensation*(M_R_trans_norm_YawTorque_allNOfreq_mean(subset)-M_R_trans_norm_YawTorque_allNOfreq_mean(mod_value_YawTorque==0)),'-c.','linewidth',1)
%     plot(Myaw_norm*mod_value_YawTorque(subset), Fdeficit_compensation*(M_L_trans_norm_YawTorque_allNOfreq_mean(subset)-M_L_trans_norm_YawTorque_allNOfreq_mean(mod_value_YawTorque==0)),'-m.','linewidth',1)
else
    plot(Myaw_norm*mod_value_YawTorque(subset),-Fdeficit_compensation*(Mx_transNrot_norm_YawTorque_allNOfreq_mean(subset)-Mx_transNrot_norm_YawTorque_allNOfreq_mean(mod_value_YawTorque==0)),'-r.','linewidth',1)
    plot(Myaw_norm*mod_value_YawTorque(subset), Fdeficit_compensation*(My_transNrot_norm_YawTorque_allNOfreq_mean(subset)-My_transNrot_norm_YawTorque_allNOfreq_mean(mod_value_YawTorque==0)),'-g.','linewidth',1)
    plot(Myaw_norm*mod_value_YawTorque(subset), Fdeficit_compensation*(Mz_transNrot_norm_YawTorque_allNOfreq_mean(subset)-Mz_transNrot_norm_YawTorque_allNOfreq_mean(mod_value_YawTorque==0)),'-b.','linewidth',1)
%     plot(Myaw_norm*mod_value_YawTorque(subset), Fdeficit_compensation*(M_R_transNrot_norm_YawTorque_allNOfreq_mean(subset)-M_R_transNrot_norm_YawTorque_allNOfreq_mean(mod_value_YawTorque==0)),'-c.','linewidth',1)
%     plot(Myaw_norm*mod_value_YawTorque(subset), Fdeficit_compensation*(M_L_transNrot_norm_YawTorque_allNOfreq_mean(subset)-M_L_transNrot_norm_YawTorque_allNOfreq_mean(mod_value_YawTorque==0)),'-m.','linewidth',1)
end

xlabel('Myaw/(a+g)ml fly','fontsize',10) 
ylabel('M/mgl qsmodel','fontsize',10)

% trendline
plot([x_min x_max],[x_min x_max],'-k')
% aYaw_norm = Myaw_norm*mod_value_YawTorque(subset);
% Myaw_norm = -(Mz_transNrot_norm_YawTorque_allNOfreq_mean(subset)-Mz_transNrot_norm_YawTorque_allNOfreq_mean(mod_value_YawTorque==0));
% aMyaw_fit = polyfit(aYaw_norm,Myaw_norm',1);
% Myaw_fit = polyval(aMyaw_fit,[x_min x_max]);
% plot([x_min x_max],Myaw_fit,'--k')

% axis equal
axis square
axis([x_min x_max y_min y_max])
set(gca,'XTick',[0:dx:x_max],'fontsize',8) 
set(gca,'YTick',[0:dy:y_max],'fontsize',8) 

% cross-talk axisL&R
figure(2)
subplot(2,3,6)
hold on

if rot_on == 0
%     plot(Myaw_norm*mod_value_YawTorque(subset),-Fdeficit_compensation*(Mx_trans_norm_YawTorque_allNOfreq_mean(subset)-Mx_trans_norm_YawTorque_allNOfreq_mean(mod_value_YawTorque==0)),'-r.','linewidth',1)
%     plot(Myaw_norm*mod_value_YawTorque(subset), Fdeficit_compensation*(My_trans_norm_YawTorque_allNOfreq_mean(subset)-My_trans_norm_YawTorque_allNOfreq_mean(mod_value_YawTorque==0)),'-g.','linewidth',1)
    plot(Myaw_norm*mod_value_YawTorque(subset), Fdeficit_compensation*(Mz_trans_norm_YawTorque_allNOfreq_mean(subset)-Mz_trans_norm_YawTorque_allNOfreq_mean(mod_value_YawTorque==0)),'-b.','linewidth',1)
    plot(Myaw_norm*mod_value_YawTorque(subset), Fdeficit_compensation*(M_R_trans_norm_YawTorque_allNOfreq_mean(subset)-M_R_trans_norm_YawTorque_allNOfreq_mean(mod_value_YawTorque==0)),'-c.','linewidth',1)
    plot(Myaw_norm*mod_value_YawTorque(subset), Fdeficit_compensation*(M_L_trans_norm_YawTorque_allNOfreq_mean(subset)-M_L_trans_norm_YawTorque_allNOfreq_mean(mod_value_YawTorque==0)),'-m.','linewidth',1)
else
%     plot(Myaw_norm*mod_value_YawTorque(subset),-Fdeficit_compensation*(Mx_transNrot_norm_YawTorque_allNOfreq_mean(subset)-Mx_transNrot_norm_YawTorque_allNOfreq_mean(mod_value_YawTorque==0)),'-r.','linewidth',1)
%     plot(Myaw_norm*mod_value_YawTorque(subset), Fdeficit_compensation*(My_transNrot_norm_YawTorque_allNOfreq_mean(subset)-My_transNrot_norm_YawTorque_allNOfreq_mean(mod_value_YawTorque==0)),'-g.','linewidth',1)
    plot(Myaw_norm*mod_value_YawTorque(subset), Fdeficit_compensation*(Mz_transNrot_norm_YawTorque_allNOfreq_mean(subset)-Mz_transNrot_norm_YawTorque_allNOfreq_mean(mod_value_YawTorque==0)),'-b.','linewidth',1)
    plot(Myaw_norm*mod_value_YawTorque(subset), Fdeficit_compensation*(M_R_transNrot_norm_YawTorque_allNOfreq_mean(subset)-M_R_transNrot_norm_YawTorque_allNOfreq_mean(mod_value_YawTorque==0)),'-c.','linewidth',1)
    plot(Myaw_norm*mod_value_YawTorque(subset), Fdeficit_compensation*(M_L_transNrot_norm_YawTorque_allNOfreq_mean(subset)-M_L_transNrot_norm_YawTorque_allNOfreq_mean(mod_value_YawTorque==0)),'-m.','linewidth',1)
end

xlabel('Myaw/(a+g)ml fly','fontsize',10) 
ylabel('M/mgl qsmodel','fontsize',10)

% trendline
plot([x_min x_max],[x_min x_max],'-k')
% aYaw_norm = Myaw_norm*mod_value_YawTorque(subset);
% Myaw_norm = -(Mz_transNrot_norm_YawTorque_allNOfreq_mean(subset)-Mz_transNrot_norm_YawTorque_allNOfreq_mean(mod_value_YawTorque==0));
% aMyaw_fit = polyfit(aYaw_norm,Myaw_norm',1);
% Myaw_fit = polyval(aMyaw_fit,[x_min x_max]);
% plot([x_min x_max],Myaw_fit,'--k')

% axis equal
axis square
axis([x_min x_max y_min y_max])
set(gca,'XTick',[0:dx:x_max],'fontsize',8) 
set(gca,'YTick',[0:dy:y_max],'fontsize',8) 

%% Mpitch
clear
% load('MOD_norm_data_QSmodel.mat')
load('norm_data_torque.mat')
load(pitchDB)

load('QSmodel_settings.mat')
if Fdeficit_on == 0
    Fdeficit_compensation = 1;
end

x_min = floor(100*min(mod_value_PitchTorque)*Mpitch_norm)/100;
x_max = ceil(100*max(mod_value_PitchTorque)*Mpitch_norm)/100;
dx = (x_max)/2;
% x_min = -max(yawTorques)/2;
% x_max = max(yawTorques);
% dx = max(yawTorques)/2;
y_min = x_min;
y_max = x_max;
dy = dx;

dmod = mod_value_PitchTorque(2)-mod_value_PitchTorque(1);
subset = find(Mpitch_norm*mod_value_PitchTorque>=x_min-dmod & Mpitch_norm*mod_value_PitchTorque<=x_max+dmod );

My_trans_norm_PitchTorque_sum_mean = ...
    (My_trans_norm_PitchTorque_stroke_mean - My_trans_norm_PitchTorque_stroke_mean(mod_value_PitchTorque==0)) +...
    (My_trans_norm_PitchTorque_dev_mean - My_trans_norm_PitchTorque_dev_mean(mod_value_PitchTorque==0)) +...
    (My_trans_norm_PitchTorque_rot_mean - My_trans_norm_PitchTorque_rot_mean(mod_value_PitchTorque==0));

My_transNrot_norm_PitchTorque_sum_mean = ...
    (My_transNrot_norm_PitchTorque_stroke_mean - My_transNrot_norm_PitchTorque_stroke_mean(mod_value_PitchTorque==0)) +...
    (My_transNrot_norm_PitchTorque_dev_mean - My_transNrot_norm_PitchTorque_dev_mean(mod_value_PitchTorque==0)) +...
    (My_transNrot_norm_PitchTorque_rot_mean - My_transNrot_norm_PitchTorque_rot_mean(mod_value_PitchTorque==0));

% R&L axes
varname=dir('flyVars*')
varname = varname.name;
load(varname)

M_R_trans_norm_PitchTorque_allNOfreq_mean =  Mx_trans_norm_PitchTorque_allNOfreq_mean*cosd(rot_axis_angle) + My_trans_norm_PitchTorque_allNOfreq_mean*sind(rot_axis_angle);
M_L_trans_norm_PitchTorque_allNOfreq_mean = -Mx_trans_norm_PitchTorque_allNOfreq_mean*sind(rot_axis_angle) + My_trans_norm_PitchTorque_allNOfreq_mean*cosd(rot_axis_angle);

M_R_transNrot_norm_PitchTorque_allNOfreq_mean =  Mx_transNrot_norm_PitchTorque_allNOfreq_mean*cosd(rot_axis_angle) + My_transNrot_norm_PitchTorque_allNOfreq_mean*sind(rot_axis_angle);
M_L_transNrot_norm_PitchTorque_allNOfreq_mean = -Mx_transNrot_norm_PitchTorque_allNOfreq_mean*sind(rot_axis_angle) + My_transNrot_norm_PitchTorque_allNOfreq_mean*cosd(rot_axis_angle);

% kinematics components
figure(1)
subplot(2,3,5)
hold on

if rot_on == 0
    plot(Mpitch_norm*mod_value_PitchTorque(subset),Fdeficit_compensation*(My_trans_norm_PitchTorque_allNOfreq_mean(subset)-My_trans_norm_PitchTorque_allNOfreq_mean(mod_value_PitchTorque==0)),'-r.','linewidth',1)
    % plot(Mpitch_norm*mod_value_PitchTorque(subset),Fdeficit_compensation*(My_trans_norm_PitchTorque_freq_mean(subset)-My_trans_norm_PitchTorque_freq_mean(mod_value_PitchTorque==0)),'-g.','linewidth',1)
    plot(Mpitch_norm*mod_value_PitchTorque(subset),Fdeficit_compensation*(My_trans_norm_PitchTorque_stroke_mean(subset)-My_trans_norm_PitchTorque_stroke_mean(mod_value_PitchTorque==0)),'-b.','linewidth',1)
    plot(Mpitch_norm*mod_value_PitchTorque(subset),Fdeficit_compensation*(My_trans_norm_PitchTorque_dev_mean(subset)-My_trans_norm_PitchTorque_dev_mean(mod_value_PitchTorque==0)),'-m.','linewidth',1)
    plot(Mpitch_norm*mod_value_PitchTorque(subset),Fdeficit_compensation*(My_trans_norm_PitchTorque_rot_mean(subset)-My_trans_norm_PitchTorque_rot_mean(mod_value_PitchTorque==0)),'-c.','linewidth',1)
    plot(Mpitch_norm*mod_value_PitchTorque(subset),Fdeficit_compensation*(My_trans_norm_PitchTorque_sum_mean(subset)-My_trans_norm_PitchTorque_sum_mean(mod_value_PitchTorque==0)),'--','color',[.5 .5 .5],'linewidth',1)
else
    plot(Mpitch_norm*mod_value_PitchTorque(subset),Fdeficit_compensation*(My_transNrot_norm_PitchTorque_allNOfreq_mean(subset)-My_transNrot_norm_PitchTorque_allNOfreq_mean(mod_value_PitchTorque==0)),'-r.','linewidth',1)
    % plot(Mpitch_norm*mod_value_PitchTorque(subset),Fdeficit_compensation*(My_transNrot_norm_PitchTorque_freq_mean(subset)-My_transNrot_norm_PitchTorque_freq_mean(mod_value_PitchTorque==0)),'-g.','linewidth',1)
    plot(Mpitch_norm*mod_value_PitchTorque(subset),Fdeficit_compensation*(My_transNrot_norm_PitchTorque_stroke_mean(subset)-My_transNrot_norm_PitchTorque_stroke_mean(mod_value_PitchTorque==0)),'-b.','linewidth',1)
    plot(Mpitch_norm*mod_value_PitchTorque(subset),Fdeficit_compensation*(My_transNrot_norm_PitchTorque_dev_mean(subset)-My_transNrot_norm_PitchTorque_dev_mean(mod_value_PitchTorque==0)),'-m.','linewidth',1)
    plot(Mpitch_norm*mod_value_PitchTorque(subset),Fdeficit_compensation*(My_transNrot_norm_PitchTorque_rot_mean(subset)-My_transNrot_norm_PitchTorque_rot_mean(mod_value_PitchTorque==0)),'-c.','linewidth',1)
    plot(Mpitch_norm*mod_value_PitchTorque(subset),Fdeficit_compensation*(My_transNrot_norm_PitchTorque_sum_mean(subset)-My_transNrot_norm_PitchTorque_sum_mean(mod_value_PitchTorque==0)),'--','color',[.5 .5 .5],'linewidth',1)
end

xlabel('Mpitch/(a+g)ml fly','fontsize',10) 
ylabel('Mpitch/mgl qsmodel','fontsize',10)

% trendline
plot([x_min x_max],[x_min x_max],'-k')
% aPitch_norm = Mpitch_norm*mod_value_PitchTorque(subset);
% Mpitch_norm = (My_transNrot_norm_PitchTorque_allNOfreq_mean(subset)-My_transNrot_norm_PitchTorque_allNOfreq_mean(mod_value_PitchTorque==0));
% aMpitch_fit = polyfit(aPitch_norm,Mpitch_norm',1);
% Mpitch_fit = polyval(aMpitch_fit,[x_min x_max]);
% plot([x_min x_max],Mpitch_fit,'--k')

% axis equal
axis square
axis([x_min x_max y_min y_max])
set(gca,'XTick',[x_min:dx:x_max],'fontsize',8) 
set(gca,'YTick',[y_min:dy:y_max],'fontsize',8) 

% cross-talk
figure(2)
subplot(2,3,5)
hold on

if rot_on == 0
    plot(Mpitch_norm*mod_value_PitchTorque(subset),-Fdeficit_compensation*(Mx_trans_norm_PitchTorque_allNOfreq_mean(subset)-Mx_trans_norm_PitchTorque_allNOfreq_mean(mod_value_PitchTorque==0)),'-r.','linewidth',1)
    plot(Mpitch_norm*mod_value_PitchTorque(subset), Fdeficit_compensation*(My_trans_norm_PitchTorque_allNOfreq_mean(subset)-My_trans_norm_PitchTorque_allNOfreq_mean(mod_value_PitchTorque==0)),'-g.','linewidth',1)
    plot(Mpitch_norm*mod_value_PitchTorque(subset), Fdeficit_compensation*(Mz_trans_norm_PitchTorque_allNOfreq_mean(subset)-Mz_trans_norm_PitchTorque_allNOfreq_mean(mod_value_PitchTorque==0)),'-b.','linewidth',1)
%     plot(Mpitch_norm*mod_value_PitchTorque(subset), Fdeficit_compensation*(M_R_trans_norm_PitchTorque_allNOfreq_mean(subset)-M_R_trans_norm_PitchTorque_allNOfreq_mean(mod_value_PitchTorque==0)),'-c.','linewidth',1)
%     plot(Mpitch_norm*mod_value_PitchTorque(subset), Fdeficit_compensation*(M_L_trans_norm_PitchTorque_allNOfreq_mean(subset)-M_L_trans_norm_PitchTorque_allNOfreq_mean(mod_value_PitchTorque==0)),'-m.','linewidth',1)
else
    plot(Mpitch_norm*mod_value_PitchTorque(subset),-Fdeficit_compensation*(Mx_transNrot_norm_PitchTorque_allNOfreq_mean(subset)-Mx_transNrot_norm_PitchTorque_allNOfreq_mean(mod_value_PitchTorque==0)),'-r.','linewidth',1)
    plot(Mpitch_norm*mod_value_PitchTorque(subset), Fdeficit_compensation*(My_transNrot_norm_PitchTorque_allNOfreq_mean(subset)-My_transNrot_norm_PitchTorque_allNOfreq_mean(mod_value_PitchTorque==0)),'-g.','linewidth',1)
    plot(Mpitch_norm*mod_value_PitchTorque(subset), Fdeficit_compensation*(Mz_transNrot_norm_PitchTorque_allNOfreq_mean(subset)-Mz_transNrot_norm_PitchTorque_allNOfreq_mean(mod_value_PitchTorque==0)),'-b.','linewidth',1)
%     plot(Mpitch_norm*mod_value_PitchTorque(subset), Fdeficit_compensation*(M_R_transNrot_norm_PitchTorque_allNOfreq_mean(subset)-M_R_transNrot_norm_PitchTorque_allNOfreq_mean(mod_value_PitchTorque==0)),'-c.','linewidth',1)
%     plot(Mpitch_norm*mod_value_PitchTorque(subset), Fdeficit_compensation*(M_L_transNrot_norm_PitchTorque_allNOfreq_mean(subset)-M_L_transNrot_norm_PitchTorque_allNOfreq_mean(mod_value_PitchTorque==0)),'-m.','linewidth',1)
end

xlabel('Mpitch/(a+g)ml fly','fontsize',10) 
ylabel('M/mgl qsmodel','fontsize',10)

% trendline
plot([x_min x_max],[x_min x_max],'-k')
% aPitch_norm = Mpitch_norm*mod_value_PitchTorque(subset);
% Mpitch_norm = (My_transNrot_norm_PitchTorque_allNOfreq_mean(subset)-My_transNrot_norm_PitchTorque_allNOfreq_mean(mod_value_PitchTorque==0));
% aMpitch_fit = polyfit(aPitch_norm,Mpitch_norm',1);
% Mpitch_fit = polyval(aMpitch_fit,[x_min x_max]);
% plot([x_min x_max],Mpitch_fit,'--k')

% axis equal
axis square
axis([x_min x_max y_min y_max])
set(gca,'XTick',[x_min:dx:x_max],'fontsize',8) 
set(gca,'YTick',[y_min:dy:y_max],'fontsize',8) 

%% MaxisR
clear
% load('MOD_norm_data_QSmodel.mat')
load('norm_data_torque.mat')
load(RaxisDB)

varname=dir('flyVars*')
varname = varname.name;
load(varname)

load('QSmodel_settings.mat')
if Fdeficit_on == 0
    Fdeficit_compensation = 1;
end

x_min = floor(100*min(mod_value_TorqueAxisR)*M_R_norm)/100;
x_max = ceil(100*max(mod_value_TorqueAxisR)*M_R_norm)/100;
dx = (x_max)/2;
% x_min = -max(yawTorques)/2;
% x_max = max(yawTorques);
% dx = max(yawTorques)/2;
y_min = x_min;
y_max = x_max;
dy = dx;

dmod = mod_value_TorqueAxisR(2)-mod_value_TorqueAxisR(1);
subset = find(M_R_norm*mod_value_TorqueAxisR>=x_min-dmod & M_R_norm*mod_value_TorqueAxisR<=x_max+dmod );

Mx_trans_norm_TorqueAxisR_sum_mean = ...
    (Mx_trans_norm_TorqueAxisR_stroke_mean - Mx_trans_norm_TorqueAxisR_stroke_mean(mod_value_TorqueAxisR==0)) +...
    (Mx_trans_norm_TorqueAxisR_dev_mean - Mx_trans_norm_TorqueAxisR_dev_mean(mod_value_TorqueAxisR==0)) +...
    (Mx_trans_norm_TorqueAxisR_rot_mean - Mx_trans_norm_TorqueAxisR_rot_mean(mod_value_TorqueAxisR==0));

Mx_transNrot_norm_TorqueAxisR_sum_mean = ...
    (Mx_transNrot_norm_TorqueAxisR_stroke_mean - Mx_transNrot_norm_TorqueAxisR_stroke_mean(mod_value_TorqueAxisR==0)) +...
    (Mx_transNrot_norm_TorqueAxisR_dev_mean - Mx_transNrot_norm_TorqueAxisR_dev_mean(mod_value_TorqueAxisR==0)) +...
    (Mx_transNrot_norm_TorqueAxisR_rot_mean - Mx_transNrot_norm_TorqueAxisR_rot_mean(mod_value_TorqueAxisR==0));

My_trans_norm_TorqueAxisR_sum_mean = ...
    (My_trans_norm_TorqueAxisR_stroke_mean - My_trans_norm_TorqueAxisR_stroke_mean(mod_value_TorqueAxisR==0)) +...
    (My_trans_norm_TorqueAxisR_dev_mean - My_trans_norm_TorqueAxisR_dev_mean(mod_value_TorqueAxisR==0)) +...
    (My_trans_norm_TorqueAxisR_rot_mean - My_trans_norm_TorqueAxisR_rot_mean(mod_value_TorqueAxisR==0));

My_transNrot_norm_TorqueAxisR_sum_mean = ...
    (My_transNrot_norm_TorqueAxisR_stroke_mean - My_transNrot_norm_TorqueAxisR_stroke_mean(mod_value_TorqueAxisR==0)) +...
    (My_transNrot_norm_TorqueAxisR_dev_mean - My_transNrot_norm_TorqueAxisR_dev_mean(mod_value_TorqueAxisR==0)) +...
    (My_transNrot_norm_TorqueAxisR_rot_mean - My_transNrot_norm_TorqueAxisR_rot_mean(mod_value_TorqueAxisR==0));

%% R&L axes
% allNOfreq
M_R_trans_norm_TorqueAxisR_allNOfreq_mean =  Mx_trans_norm_TorqueAxisR_allNOfreq_mean*cosd(rot_axis_angle) + My_trans_norm_TorqueAxisR_allNOfreq_mean*sind(rot_axis_angle);
M_L_trans_norm_TorqueAxisR_allNOfreq_mean = -Mx_trans_norm_TorqueAxisR_allNOfreq_mean*sind(rot_axis_angle) + My_trans_norm_TorqueAxisR_allNOfreq_mean*cosd(rot_axis_angle);

M_R_transNrot_norm_TorqueAxisR_allNOfreq_mean =  Mx_transNrot_norm_TorqueAxisR_allNOfreq_mean*cosd(rot_axis_angle) + My_transNrot_norm_TorqueAxisR_allNOfreq_mean*sind(rot_axis_angle);
M_L_transNrot_norm_TorqueAxisR_allNOfreq_mean = -Mx_transNrot_norm_TorqueAxisR_allNOfreq_mean*sind(rot_axis_angle) + My_transNrot_norm_TorqueAxisR_allNOfreq_mean*cosd(rot_axis_angle);

M_strkp_trans_norm_TorqueAxisR_allNOfreq_mean     = sqrt(Mx_trans_norm_TorqueAxisR_allNOfreq_mean    .^2 + My_trans_norm_TorqueAxisR_allNOfreq_mean    .^2);
M_strkp_transNrot_norm_TorqueAxisR_allNOfreq_mean = sqrt(Mx_transNrot_norm_TorqueAxisR_allNOfreq_mean.^2 + My_transNrot_norm_TorqueAxisR_allNOfreq_mean.^2);

M_strkpAxis_trans_norm_TorqueAxisR_allNOfreq_mean     = atand(My_trans_norm_TorqueAxisR_allNOfreq_mean     ./ Mx_trans_norm_TorqueAxisR_allNOfreq_mean);
M_strkpAxis_transNrot_norm_TorqueAxisR_allNOfreq_mean = atand(My_transNrot_norm_TorqueAxisR_allNOfreq_mean ./ Mx_transNrot_norm_TorqueAxisR_allNOfreq_mean);

% stroke
M_R_trans_norm_TorqueAxisR_stroke_mean =  Mx_trans_norm_TorqueAxisR_stroke_mean*cosd(rot_axis_angle) + My_trans_norm_TorqueAxisR_stroke_mean*sind(rot_axis_angle);
M_L_trans_norm_TorqueAxisR_stroke_mean = -Mx_trans_norm_TorqueAxisR_stroke_mean*sind(rot_axis_angle) + My_trans_norm_TorqueAxisR_stroke_mean*cosd(rot_axis_angle);

M_R_transNrot_norm_TorqueAxisR_stroke_mean =  Mx_transNrot_norm_TorqueAxisR_stroke_mean*cosd(rot_axis_angle) + My_transNrot_norm_TorqueAxisR_stroke_mean*sind(rot_axis_angle);
M_L_transNrot_norm_TorqueAxisR_stroke_mean = -Mx_transNrot_norm_TorqueAxisR_stroke_mean*sind(rot_axis_angle) + My_transNrot_norm_TorqueAxisR_stroke_mean*cosd(rot_axis_angle);

M_strkp_trans_norm_TorqueAxisR_stroke_mean     = sqrt(Mx_trans_norm_TorqueAxisR_stroke_mean    .^2 + My_trans_norm_TorqueAxisR_stroke_mean    .^2);
M_strkp_transNrot_norm_TorqueAxisR_stroke_mean = sqrt(Mx_transNrot_norm_TorqueAxisR_stroke_mean.^2 + My_transNrot_norm_TorqueAxisR_stroke_mean.^2);

M_strkpAxis_trans_norm_TorqueAxisR_stroke_mean     = atand(My_trans_norm_TorqueAxisR_stroke_mean     ./ Mx_trans_norm_TorqueAxisR_stroke_mean);
M_strkpAxis_transNrot_norm_TorqueAxisR_stroke_mean = atand(My_transNrot_norm_TorqueAxisR_stroke_mean ./ Mx_transNrot_norm_TorqueAxisR_stroke_mean);

% dev
M_R_trans_norm_TorqueAxisR_dev_mean =  Mx_trans_norm_TorqueAxisR_dev_mean*cosd(rot_axis_angle) + My_trans_norm_TorqueAxisR_dev_mean*sind(rot_axis_angle);
M_L_trans_norm_TorqueAxisR_dev_mean = -Mx_trans_norm_TorqueAxisR_dev_mean*sind(rot_axis_angle) + My_trans_norm_TorqueAxisR_dev_mean*cosd(rot_axis_angle);

M_R_transNrot_norm_TorqueAxisR_dev_mean =  Mx_transNrot_norm_TorqueAxisR_dev_mean*cosd(rot_axis_angle) + My_transNrot_norm_TorqueAxisR_dev_mean*sind(rot_axis_angle);
M_L_transNrot_norm_TorqueAxisR_dev_mean = -Mx_transNrot_norm_TorqueAxisR_dev_mean*sind(rot_axis_angle) + My_transNrot_norm_TorqueAxisR_dev_mean*cosd(rot_axis_angle);

M_strkp_trans_norm_TorqueAxisR_dev_mean     = sqrt(Mx_trans_norm_TorqueAxisR_dev_mean    .^2 + My_trans_norm_TorqueAxisR_dev_mean    .^2);
M_strkp_transNrot_norm_TorqueAxisR_dev_mean = sqrt(Mx_transNrot_norm_TorqueAxisR_dev_mean.^2 + My_transNrot_norm_TorqueAxisR_dev_mean.^2);

M_strkpAxis_trans_norm_TorqueAxisR_dev_mean     = atand(My_trans_norm_TorqueAxisR_dev_mean     ./ Mx_trans_norm_TorqueAxisR_dev_mean);
M_strkpAxis_transNrot_norm_TorqueAxisR_dev_mean = atand(My_transNrot_norm_TorqueAxisR_dev_mean ./ Mx_transNrot_norm_TorqueAxisR_dev_mean);

% rot
M_R_trans_norm_TorqueAxisR_rot_mean =  Mx_trans_norm_TorqueAxisR_rot_mean*cosd(rot_axis_angle) + My_trans_norm_TorqueAxisR_rot_mean*sind(rot_axis_angle);
M_L_trans_norm_TorqueAxisR_rot_mean = -Mx_trans_norm_TorqueAxisR_rot_mean*sind(rot_axis_angle) + My_trans_norm_TorqueAxisR_rot_mean*cosd(rot_axis_angle);

M_R_transNrot_norm_TorqueAxisR_rot_mean =  Mx_transNrot_norm_TorqueAxisR_rot_mean*cosd(rot_axis_angle) + My_transNrot_norm_TorqueAxisR_rot_mean*sind(rot_axis_angle);
M_L_transNrot_norm_TorqueAxisR_rot_mean = -Mx_transNrot_norm_TorqueAxisR_rot_mean*sind(rot_axis_angle) + My_transNrot_norm_TorqueAxisR_rot_mean*cosd(rot_axis_angle);

M_strkp_trans_norm_TorqueAxisR_rot_mean     = sqrt(Mx_trans_norm_TorqueAxisR_rot_mean    .^2 + My_trans_norm_TorqueAxisR_rot_mean    .^2);
M_strkp_transNrot_norm_TorqueAxisR_rot_mean = sqrt(Mx_transNrot_norm_TorqueAxisR_rot_mean.^2 + My_transNrot_norm_TorqueAxisR_rot_mean.^2);

M_strkpAxis_trans_norm_TorqueAxisR_rot_mean     = atand(My_trans_norm_TorqueAxisR_rot_mean     ./ Mx_trans_norm_TorqueAxisR_rot_mean);
M_strkpAxis_transNrot_norm_TorqueAxisR_rot_mean = atand(My_transNrot_norm_TorqueAxisR_rot_mean ./ Mx_transNrot_norm_TorqueAxisR_rot_mean);

% sum
M_R_trans_norm_TorqueAxisR_sum_mean =  Mx_trans_norm_TorqueAxisR_sum_mean*cosd(rot_axis_angle) + My_trans_norm_TorqueAxisR_sum_mean*sind(rot_axis_angle);
M_L_trans_norm_TorqueAxisR_sum_mean = -Mx_trans_norm_TorqueAxisR_sum_mean*sind(rot_axis_angle) + My_trans_norm_TorqueAxisR_sum_mean*cosd(rot_axis_angle);

M_R_transNrot_norm_TorqueAxisR_sum_mean =  Mx_transNrot_norm_TorqueAxisR_sum_mean*cosd(rot_axis_angle) + My_transNrot_norm_TorqueAxisR_sum_mean*sind(rot_axis_angle);
M_L_transNrot_norm_TorqueAxisR_sum_mean = -Mx_transNrot_norm_TorqueAxisR_sum_mean*sind(rot_axis_angle) + My_transNrot_norm_TorqueAxisR_sum_mean*cosd(rot_axis_angle);

M_strkp_trans_norm_TorqueAxisR_sum_mean     = sqrt(Mx_trans_norm_TorqueAxisR_sum_mean    .^2 + My_trans_norm_TorqueAxisR_sum_mean    .^2);
M_strkp_transNrot_norm_TorqueAxisR_sum_mean = sqrt(Mx_transNrot_norm_TorqueAxisR_sum_mean.^2 + My_transNrot_norm_TorqueAxisR_sum_mean.^2);

M_strkpAxis_trans_norm_TorqueAxisR_sum_mean     = atand(My_trans_norm_TorqueAxisR_sum_mean     ./ Mx_trans_norm_TorqueAxisR_sum_mean);
M_strkpAxis_transNrot_norm_TorqueAxisR_sum_mean = atand(My_transNrot_norm_TorqueAxisR_sum_mean ./ Mx_transNrot_norm_TorqueAxisR_sum_mean);

%% kinematics components
figure(1)

% M_R
subplot(2,3,3)
hold on

if rot_on == 0
    plot(M_R_norm*mod_value_TorqueAxisR(subset),Fdeficit_compensation*(M_R_trans_norm_TorqueAxisR_allNOfreq_mean(subset)-M_R_trans_norm_TorqueAxisR_allNOfreq_mean(mod_value_TorqueAxisR==0)),'-r.','linewidth',1)
    % plot(M_R_norm*mod_value_TorqueAxisR(subset),Fdeficit_compensation*(M_R_trans_norm_TorqueAxisR_freq_mean(subset)-M_R_trans_norm_TorqueAxisR_freq_mean(mod_value_TorqueAxisR==0)),'-g.','linewidth',1)
    plot(M_R_norm*mod_value_TorqueAxisR(subset),Fdeficit_compensation*(M_R_trans_norm_TorqueAxisR_stroke_mean(subset)-M_R_trans_norm_TorqueAxisR_stroke_mean(mod_value_TorqueAxisR==0)),'-b.','linewidth',1)
    plot(M_R_norm*mod_value_TorqueAxisR(subset),Fdeficit_compensation*(M_R_trans_norm_TorqueAxisR_dev_mean(subset)-M_R_trans_norm_TorqueAxisR_dev_mean(mod_value_TorqueAxisR==0)),'-m.','linewidth',1)
    plot(M_R_norm*mod_value_TorqueAxisR(subset),Fdeficit_compensation*(M_R_trans_norm_TorqueAxisR_rot_mean(subset)-M_R_trans_norm_TorqueAxisR_rot_mean(mod_value_TorqueAxisR==0)),'-c.','linewidth',1)
    plot(M_R_norm*mod_value_TorqueAxisR(subset),Fdeficit_compensation*(M_R_trans_norm_TorqueAxisR_sum_mean(subset)-M_R_trans_norm_TorqueAxisR_sum_mean(mod_value_TorqueAxisR==0)),'--','color',[.5 .5 .5],'linewidth',1)
else
    plot(M_R_norm*mod_value_TorqueAxisR(subset),Fdeficit_compensation*(M_R_transNrot_norm_TorqueAxisR_allNOfreq_mean(subset)-M_R_transNrot_norm_TorqueAxisR_allNOfreq_mean(mod_value_TorqueAxisR==0)),'-r.','linewidth',1)
    % plot(M_R_norm*mod_value_TorqueAxisR(subset),Fdeficit_compensation*(M_R_transNrot_norm_TorqueAxisR_freq_mean(subset)-M_R_transNrot_norm_TorqueAxisR_freq_mean(mod_value_TorqueAxisR==0)),'-g.','linewidth',1)
    plot(M_R_norm*mod_value_TorqueAxisR(subset),Fdeficit_compensation*(M_R_transNrot_norm_TorqueAxisR_stroke_mean(subset)-M_R_transNrot_norm_TorqueAxisR_stroke_mean(mod_value_TorqueAxisR==0)),'-b.','linewidth',1)
    plot(M_R_norm*mod_value_TorqueAxisR(subset),Fdeficit_compensation*(M_R_transNrot_norm_TorqueAxisR_dev_mean(subset)-M_R_transNrot_norm_TorqueAxisR_dev_mean(mod_value_TorqueAxisR==0)),'-m.','linewidth',1)
    plot(M_R_norm*mod_value_TorqueAxisR(subset),Fdeficit_compensation*(M_R_transNrot_norm_TorqueAxisR_rot_mean(subset)-M_R_transNrot_norm_TorqueAxisR_rot_mean(mod_value_TorqueAxisR==0)),'-c.','linewidth',1)
    plot(M_R_norm*mod_value_TorqueAxisR(subset),Fdeficit_compensation*(M_R_transNrot_norm_TorqueAxisR_sum_mean(subset)-M_R_transNrot_norm_TorqueAxisR_sum_mean(mod_value_TorqueAxisR==0)),'--','color',[.5 .5 .5],'linewidth',1)
end

xlabel('M_R/(a+g)ml fly','fontsize',10) 
ylabel('M_R/mgl qsmodel','fontsize',10)

% trendline
plot([x_min x_max],[x_min x_max],'-k')
% aPitch_norm = M_R_norm*mod_value_TorqueAxisR(subset);
% M_R_norm = (M_R_transNrot_norm_TorqueAxisR_allNOfreq_mean(subset)-M_R_transNrot_norm_TorqueAxisR_allNOfreq_mean(mod_value_TorqueAxisR==0));
% aM_R_fit = polyfit(aPitch_norm,M_R_norm',1);
% M_R_fit = polyval(aM_R_fit,[x_min x_max]);
% plot([x_min x_max],M_R_fit,'--k')

% axis equal
axis square
axis([x_min x_max y_min y_max])
set(gca,'XTick',[x_min:dx:x_max],'fontsize',8) 
set(gca,'YTick',[y_min:dy:y_max],'fontsize',8) 

% M_L
subplot(2,3,6)
hold on

if rot_on == 0
    plot(M_R_norm*mod_value_TorqueAxisR(subset),Fdeficit_compensation*(M_L_trans_norm_TorqueAxisR_allNOfreq_mean(subset)-M_L_trans_norm_TorqueAxisR_allNOfreq_mean(mod_value_TorqueAxisR==0)),'-r.','linewidth',1)
    % plot(M_R_norm*mod_value_TorqueAxisR(subset),Fdeficit_compensation*(M_L_trans_norm_TorqueAxisR_freq_mean(subset)-M_L_trans_norm_TorqueAxisR_freq_mean(mod_value_TorqueAxisR==0)),'-g.','linewidth',1)
    plot(M_R_norm*mod_value_TorqueAxisR(subset),Fdeficit_compensation*(M_L_trans_norm_TorqueAxisR_stroke_mean(subset)-M_L_trans_norm_TorqueAxisR_stroke_mean(mod_value_TorqueAxisR==0)),'-b.','linewidth',1)
    plot(M_R_norm*mod_value_TorqueAxisR(subset),Fdeficit_compensation*(M_L_trans_norm_TorqueAxisR_dev_mean(subset)-M_L_trans_norm_TorqueAxisR_dev_mean(mod_value_TorqueAxisR==0)),'-m.','linewidth',1)
    plot(M_R_norm*mod_value_TorqueAxisR(subset),Fdeficit_compensation*(M_L_trans_norm_TorqueAxisR_rot_mean(subset)-M_L_trans_norm_TorqueAxisR_rot_mean(mod_value_TorqueAxisR==0)),'-c.','linewidth',1)
    plot(M_R_norm*mod_value_TorqueAxisR(subset),Fdeficit_compensation*(M_L_trans_norm_TorqueAxisR_sum_mean(subset)-M_L_trans_norm_TorqueAxisR_sum_mean(mod_value_TorqueAxisR==0)),'--','color',[.5 .5 .5],'linewidth',1)
else
    plot(M_R_norm*mod_value_TorqueAxisR(subset),Fdeficit_compensation*(M_L_transNrot_norm_TorqueAxisR_allNOfreq_mean(subset)-M_L_transNrot_norm_TorqueAxisR_allNOfreq_mean(mod_value_TorqueAxisR==0)),'-r.','linewidth',1)
    % plot(M_R_norm*mod_value_TorqueAxisR(subset),Fdeficit_compensation*(M_L_transNrot_norm_TorqueAxisR_freq_mean(subset)-M_L_transNrot_norm_TorqueAxisR_freq_mean(mod_value_TorqueAxisR==0)),'-g.','linewidth',1)
    plot(M_R_norm*mod_value_TorqueAxisR(subset),Fdeficit_compensation*(M_L_transNrot_norm_TorqueAxisR_stroke_mean(subset)-M_L_transNrot_norm_TorqueAxisR_stroke_mean(mod_value_TorqueAxisR==0)),'-b.','linewidth',1)
    plot(M_R_norm*mod_value_TorqueAxisR(subset),Fdeficit_compensation*(M_L_transNrot_norm_TorqueAxisR_dev_mean(subset)-M_L_transNrot_norm_TorqueAxisR_dev_mean(mod_value_TorqueAxisR==0)),'-m.','linewidth',1)
    plot(M_R_norm*mod_value_TorqueAxisR(subset),Fdeficit_compensation*(M_L_transNrot_norm_TorqueAxisR_rot_mean(subset)-M_L_transNrot_norm_TorqueAxisR_rot_mean(mod_value_TorqueAxisR==0)),'-c.','linewidth',1)
    plot(M_R_norm*mod_value_TorqueAxisR(subset),Fdeficit_compensation*(M_L_transNrot_norm_TorqueAxisR_sum_mean(subset)-M_L_transNrot_norm_TorqueAxisR_sum_mean(mod_value_TorqueAxisR==0)),'--','color',[.5 .5 .5],'linewidth',1)
end

xlabel('M_R/(a+g)ml fly','fontsize',10) 
ylabel('M_L/mgl qsmodel','fontsize',10)

% % trendline
% plot([x_min x_max],[x_min x_max],'-k')
% aPitch_norm = M_R_norm*mod_value_TorqueAxisR(subset);
% M_R_norm = (M_L_transNrot_norm_TorqueAxisR_allNOfreq_mean(subset)-M_L_transNrot_norm_TorqueAxisR_allNOfreq_mean(mod_value_TorqueAxisR==0));
% aM_L_fit = polyfit(aPitch_norm,M_R_norm',1);
% M_L_fit = polyval(aM_L_fit,[x_min x_max]);
% plot([x_min x_max],M_L_fit,'--k')

% axis equal
axis square
axis([x_min x_max y_min y_max])
set(gca,'XTick',[x_min:dx:x_max],'fontsize',8) 
set(gca,'YTick',[y_min:dy:y_max],'fontsize',8) 

%% cross-talk
figure(2)
subplot(2,3,3)
hold on

if rot_on == 0
%     plot(M_R_norm*mod_value_TorqueAxisR(subset), Fdeficit_compensation*(Mx_trans_norm_TorqueAxisR_allNOfreq_mean(subset)-Mx_trans_norm_TorqueAxisR_allNOfreq_mean(mod_value_TorqueAxisR==0)),'-r.','linewidth',1)
%     plot(M_R_norm*mod_value_TorqueAxisR(subset), Fdeficit_compensation*(My_trans_norm_TorqueAxisR_allNOfreq_mean(subset)-My_trans_norm_TorqueAxisR_allNOfreq_mean(mod_value_TorqueAxisR==0)),'-g.','linewidth',1)
    plot(M_R_norm*mod_value_TorqueAxisR(subset), Fdeficit_compensation*(Mz_trans_norm_TorqueAxisR_allNOfreq_mean(subset)-Mz_trans_norm_TorqueAxisR_allNOfreq_mean(mod_value_TorqueAxisR==0)),'-b.','linewidth',1)
    plot(M_R_norm*mod_value_TorqueAxisR(subset), Fdeficit_compensation*(M_R_trans_norm_TorqueAxisR_allNOfreq_mean(subset)-M_R_trans_norm_TorqueAxisR_allNOfreq_mean(mod_value_TorqueAxisR==0)),'-c.','linewidth',1)
    plot(M_R_norm*mod_value_TorqueAxisR(subset), Fdeficit_compensation*(M_L_trans_norm_TorqueAxisR_allNOfreq_mean(subset)-M_L_trans_norm_TorqueAxisR_allNOfreq_mean(mod_value_TorqueAxisR==0)),'-m.','linewidth',1)
else
%     plot(M_R_norm*mod_value_TorqueAxisR(subset), Fdeficit_compensation*(Mx_transNrot_norm_TorqueAxisR_allNOfreq_mean(subset)-Mx_transNrot_norm_TorqueAxisR_allNOfreq_mean(mod_value_TorqueAxisR==0)),'-r.','linewidth',1)
%     plot(M_R_norm*mod_value_TorqueAxisR(subset), Fdeficit_compensation*(My_transNrot_norm_TorqueAxisR_allNOfreq_mean(subset)-My_transNrot_norm_TorqueAxisR_allNOfreq_mean(mod_value_TorqueAxisR==0)),'-g.','linewidth',1)
    plot(M_R_norm*mod_value_TorqueAxisR(subset), Fdeficit_compensation*(Mz_transNrot_norm_TorqueAxisR_allNOfreq_mean(subset)-Mz_transNrot_norm_TorqueAxisR_allNOfreq_mean(mod_value_TorqueAxisR==0)),'-b.','linewidth',1)
    plot(M_R_norm*mod_value_TorqueAxisR(subset), Fdeficit_compensation*(M_R_transNrot_norm_TorqueAxisR_allNOfreq_mean(subset)-M_R_transNrot_norm_TorqueAxisR_allNOfreq_mean(mod_value_TorqueAxisR==0)),'-c.','linewidth',1)
    plot(M_R_norm*mod_value_TorqueAxisR(subset), Fdeficit_compensation*(M_L_transNrot_norm_TorqueAxisR_allNOfreq_mean(subset)-M_L_transNrot_norm_TorqueAxisR_allNOfreq_mean(mod_value_TorqueAxisR==0)),'-m.','linewidth',1)
end

xlabel('M_R/(a+g)ml fly','fontsize',10) 
ylabel('M/mgl qsmodel','fontsize',10)

% trendline
plot([x_min x_max],[x_min x_max],'-k')
% aPitch_norm = M_R_norm*mod_value_TorqueAxisR(subset);
% M_R_norm = (M_R_transNrot_norm_TorqueAxisR_allNOfreq_mean(subset)-M_R_transNrot_norm_TorqueAxisR_allNOfreq_mean(mod_value_TorqueAxisR==0));
% aM_R_fit = polyfit(aPitch_norm,M_R_norm',1);
% M_R_fit = polyval(aM_R_fit,[x_min x_max]);
% plot([x_min x_max],M_R_fit,'--k')

% axis equal
axis square
axis([x_min x_max y_min y_max])
set(gca,'XTick',[x_min:dx:x_max],'fontsize',8) 
set(gca,'YTick',[y_min:dy:y_max],'fontsize',8) 

%% save plot
mkdir('MSfigs_QSmodel')
cd('MSfigs_QSmodel')

% kinematics components
figure(1)
subplot(2,3,4)
legend('all','stroke','dev','rot','sum','location','NW')

if Fdeficit_on == 0 && rot_on == 0
    saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch_MaxisR_QSmodel_transONLY.fig')
    saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch_MaxisR_QSmodel_transONLY.png')
    plot2svg(['MSfig_Fenhance_Mroll_Myaw_Mpitch_MaxisR_QSmodel_transONLY.svg'])
    
elseif Fdeficit_on == 0 && rot_on == 1
    saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch_MaxisR_QSmodel_transNrot.fig')
    saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch_MaxisR_QSmodel_transNrot.png')
    plot2svg(['MSfig_Fenhance_Mroll_Myaw_Mpitch_MaxisR_QSmodel_transNrot.svg'])
    
elseif Fdeficit_on == 1 && rot_on == 0
    saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch_MaxisR_QSmodel_Fcompensation_transONLY.fig')
    saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch_MaxisR_QSmodel_Fcompensation_transONLY.png')
    plot2svg(['MSfig_Fenhance_Mroll_Myaw_Mpitch_MaxisR_QSmodel_Fcompensation_transONLY.svg'])
    
elseif Fdeficit_on == 1 && rot_on == 1
    saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch_MaxisR_QSmodel_Fcompensation_transNrot.fig')
    saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch_MaxisR_QSmodel_Fcompensation_transNrot.png')
    plot2svg(['MSfig_Fenhance_Mroll_Myaw_Mpitch_MaxisR_QSmodel_Fcompensation_transNrot.svg'])
end
    
% cross-talk
figure(2)
subplot(2,3,1)
legend('Fx','Fy','Fz','location','NW')

subplot(2,3,4)
legend('Mroll','Mpitch','Myaw','location','NW')

subplot(2,3,6)
legend('Myaw','M_R','M_L','location','NW')

if Fdeficit_on == 0 && rot_on == 0
    saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch_MaxisR_QSmodel_crosstalk_transONLY.fig')
    saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch_MaxisR_QSmodel_crosstalk_transONLY.png')
    plot2svg(['MSfig_Fenhance_Mroll_Myaw_Mpitch_MaxisR_QSmodel_crosstalk_transONLY.svg'])
    
elseif Fdeficit_on == 0 && rot_on == 1
    saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch_MaxisR_QSmodel_crosstalk_transNrot.fig')
    saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch_MaxisR_QSmodel_crosstalk_transNrot.png')
    plot2svg(['MSfig_Fenhance_Mroll_Myaw_Mpitch_MaxisR_QSmodel_crosstalk_transNrot.svg'])
    
elseif Fdeficit_on == 1 && rot_on == 0
    saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch_MaxisR_QSmodel_crosstalk_Fcompensation_transONLY.fig')
    saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch_MaxisR_QSmodel_crosstalk_Fcompensation_transONLY.png')
    plot2svg(['MSfig_Fenhance_Mroll_Myaw_Mpitch_MaxisR_QSmodel_crosstalk_Fcompensation_transONLY.svg'])
    
elseif Fdeficit_on == 1 && rot_on == 1
    saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch_MaxisR_QSmodel_crosstalk_Fcompensation_transNrot.fig')
    saveas(gca, 'MSfig_Fenhance_Mroll_Myaw_Mpitch_MaxisR_QSmodel_crosstalk_Fcompensation_transNrot.png')
    plot2svg(['MSfig_Fenhance_Mroll_Myaw_Mpitch_MaxisR_QSmodel_crosstalk_Fcompensation_transNrot.svg'])
end
cd ..

%% delete settings file
delete('QSmodel_settings.mat')

%% save trenddata
% save('InertiaMoment_trends_QSmodel_transNrot.mat','aRoll_norm','Mroll_norm','aMroll_fit','aYaw_norm','Myaw_norm','aMyaw_fit','aPitch_norm','Mpitch_norm','aMpitch_fit','f_wb_steady')
