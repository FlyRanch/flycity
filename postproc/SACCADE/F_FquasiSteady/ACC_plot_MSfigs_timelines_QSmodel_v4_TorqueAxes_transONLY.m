% make borf plots F&M timelines_QSmodel

clc
close all


%% Fenhance
clear
% load('MOD_norm_data_QSmodel.mat')
load('norm_data_torque.mat')
load(FDB)

% R&L axes
varname=dir('flyVars*')
varname = varname.name;
load(varname)

M_R_trans_norm_Fenhance_all =  Mx_trans_norm_Fenhance_all*cosd(rot_axis_angle) + My_trans_norm_Fenhance_all*sind(rot_axis_angle);
M_L_trans_norm_Fenhance_all = -Mx_trans_norm_Fenhance_all*sind(rot_axis_angle) + My_trans_norm_Fenhance_all*cosd(rot_axis_angle);

M_R_transNrot_norm_Fenhance_all =  Mx_transNrot_norm_Fenhance_all*cosd(rot_axis_angle) + My_transNrot_norm_Fenhance_all*sind(rot_axis_angle);
M_L_transNrot_norm_Fenhance_all = -Mx_transNrot_norm_Fenhance_all*sind(rot_axis_angle) + My_transNrot_norm_Fenhance_all*cosd(rot_axis_angle);

% mod_values = FenhMods;
% color_map = [0 0 0; 0 .5 1; 0 1 1];
mod_values = [min(FenhMods) max(FenhMods)];
color_map = [0 0 0; 0 1 1];

figure
for i = 1:length(mod_values)
    mod_now = mod_values(i);
    
    mod_diff = abs(mod_value_Fenhance - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % Fx ALL
    val_now = Fx_trans_norm_Fenhance_all(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,1)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fx/Mg')
    
    % Fy ALL
    val_now = Fy_trans_norm_Fenhance_all(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,4)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fy/Mg')
    
    % Fz ALL
    val_now = Fz_trans_norm_Fenhance_all(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,7)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fz/Mg')
    
    % Mx ALL
    val_now = Mx_trans_norm_Fenhance_all(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,2)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mx/Mg')
    
    % My ALL
    val_now = My_trans_norm_Fenhance_all(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,5)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('My/Mg')
    
    % Mz ALL
    val_now = Mz_trans_norm_Fenhance_all(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,8)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mz/Mg')
    
    % M_R ALL
    val_now = M_R_trans_norm_Fenhance_all(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,3)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M_R/Mg')

    % M_L ALL
    val_now = M_L_trans_norm_Fenhance_all(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,6)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M_L/Mg')

    % Fz ALL & mean
    val_now = -Fz_trans_norm_Fenhance_all(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    val_mean = mean(val_now);
    
    subplot(3,3,9)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('-Fz/Mg')
    plot([0 1],[val_mean val_mean],'--','color',color_map(i,:),'linewidth',1)
    val_pp = csaps(t_now,val_now,.99999)
%     fnplt(val_pp);
        ylim([-.5 4])
    ylabel('Fz/Mg')
end

mkdir('MSfigs_QSmodel')
cd('MSfigs_QSmodel')
saveas(gca, 'MSfig_Fenhance_timelines_QSmodel_transONLY.fig')
saveas(gca, 'MSfig_Fenhance_timelines_QSmodel_transONLY.png')
plot2svg(['MSfig_Fenhance_timelines_QSmodel_transONLY.svg'])
cd ..

%% Mroll
clear
% load('MOD_norm_data_QSmodel.mat')
load('norm_data_torque.mat')
load(rollDB)

% R&L axes
varname=dir('flyVars*')
varname = varname.name;
load(varname)

M_R_trans_norm_RollTorque_allNOfreq =  Mx_trans_norm_RollTorque_allNOfreq*cosd(rot_axis_angle) + My_trans_norm_RollTorque_allNOfreq*sind(rot_axis_angle);
M_L_trans_norm_RollTorque_allNOfreq = -Mx_trans_norm_RollTorque_allNOfreq*sind(rot_axis_angle) + My_trans_norm_RollTorque_allNOfreq*cosd(rot_axis_angle);

M_R_transNrot_norm_RollTorque_allNOfreq =  Mx_transNrot_norm_RollTorque_allNOfreq*cosd(rot_axis_angle) + My_transNrot_norm_RollTorque_allNOfreq*sind(rot_axis_angle);
M_L_transNrot_norm_RollTorque_allNOfreq = -Mx_transNrot_norm_RollTorque_allNOfreq*sind(rot_axis_angle) + My_transNrot_norm_RollTorque_allNOfreq*cosd(rot_axis_angle);

% steady Mx
mod_now = 0;
mod_diff = abs(mod_value_RollTorque - mod_now);
n = find(mod_diff==min(mod_diff));
val_steady = Mx_trans_norm_RollTorque_allNOfreq(:,n);

% mod plots
mod_values = RollMods;
% color_map = [0 0 0; 0 0 1; 0 .5 1];
% mod_values = [min(RollMods) max(RollMods)];
color_map = [0 0 0; 0 1 1];

figure
for i = 1:length(mod_values)
    
    mod_now = mod_values(i);
    
    mod_diff = abs(mod_value_RollTorque - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % Fx ALL
    val_now = Fx_trans_norm_RollTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,1)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fx/Mg')
    
    % Fy ALL
    val_now = Fy_trans_norm_RollTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,4)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fy/Mg')
    
    % Fz ALL
    val_now = Fz_trans_norm_RollTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,7)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fz/Mg')
    
    % Mx ALL
    val_now = Mx_trans_norm_RollTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,2)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mx/Mg')
    
    % My ALL
    val_now = My_trans_norm_RollTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,5)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('My/Mg')
    
    % Mz ALL
    val_now = Mz_trans_norm_RollTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,8)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mz/Mg')
    
    % M_R ALL
    val_now = M_R_trans_norm_RollTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,3)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M_R/Mg')

    % M_L ALL
    val_now = M_L_trans_norm_RollTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,6)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M_L/Mg')

    % Mx ALL & mean
    val_now = -Mx_trans_norm_RollTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    val_mean = mean(val_now);
    
    subplot(3,3,9)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mx/Mg')
    plot([0 1],[val_mean val_mean],'--','color',color_map(i,:),'linewidth',1)
    val_pp = csaps(t_now,val_now,.99999)
%     fnplt(val_pp);
    ylim([-.1 .3])
    ylabel('Mx/Mg')    
end

mkdir('MSfigs_QSmodel')
cd('MSfigs_QSmodel')
saveas(gca, 'MSfig_Mroll_timelines_QSmodel_transONLY.fig')
saveas(gca, 'MSfig_Mroll_timelines_QSmodel_transONLY.png')
plot2svg(['MSfig_Mroll_timelines_QSmodel_transONLY.svg'])
cd ..

%% Myaw
clear
% load('MOD_norm_data_QSmodel.mat')
load('norm_data_torque.mat')
load(yawDB)

% R&L axes
varname=dir('flyVars*')
varname = varname.name;
load(varname)

M_R_trans_norm_YawTorque_allNOfreq =  Mx_trans_norm_YawTorque_allNOfreq*cosd(rot_axis_angle) + My_trans_norm_YawTorque_allNOfreq*sind(rot_axis_angle);
M_L_trans_norm_YawTorque_allNOfreq = -Mx_trans_norm_YawTorque_allNOfreq*sind(rot_axis_angle) + My_trans_norm_YawTorque_allNOfreq*cosd(rot_axis_angle);

M_R_transNrot_norm_YawTorque_allNOfreq =  Mx_transNrot_norm_YawTorque_allNOfreq*cosd(rot_axis_angle) + My_transNrot_norm_YawTorque_allNOfreq*sind(rot_axis_angle);
M_L_transNrot_norm_YawTorque_allNOfreq = -Mx_transNrot_norm_YawTorque_allNOfreq*sind(rot_axis_angle) + My_transNrot_norm_YawTorque_allNOfreq*cosd(rot_axis_angle);

% steady Mz
mod_now = 0;
mod_diff = abs(mod_value_YawTorque - mod_now);
n = find(mod_diff==min(mod_diff));
val_steady = Mz_trans_norm_YawTorque_allNOfreq(:,n);

% mod plots
mod_values = YawMods;
% color_map = [0 0 0; 0 0 1; 0 .5 1];
% mod_values = [min(YawMods) max(YawMods)];
color_map = [0 0 0; 0 1 1];

figure
for i = 1:length(mod_values)
    
    mod_now = mod_values(i);
    
    mod_diff = abs(mod_value_YawTorque - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % Fx ALL
    val_now = Fx_trans_norm_YawTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,1)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fx/Mg')
    
    % Fy ALL
    val_now = Fy_trans_norm_YawTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,4)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fy/Mg')
    
    % Fz ALL
    val_now = Fz_trans_norm_YawTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,7)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fz/Mg')
    
    % Mx ALL
    val_now = Mx_trans_norm_YawTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,2)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mx/Mg')
    
    % My ALL
    val_now = My_trans_norm_YawTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,5)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('My/Mg')
    
    % Mz ALL
    val_now = Mz_trans_norm_YawTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,8)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mz/Mg')
    
    % M_R ALL
    val_now = M_R_trans_norm_YawTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,3)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M_R/Mg')

    % M_L ALL
    val_now = M_L_trans_norm_YawTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,6)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M_L/Mg')

    % Mz ALL & mean
    val_now = Mz_trans_norm_YawTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    val_mean = mean(val_now);
    
    subplot(3,3,9)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mz/Mg')
    plot([0 1],[val_mean val_mean],'--','color',color_map(i,:),'linewidth',1)
    val_pp = csaps(t_now,val_now,.99999)
%     fnplt(val_pp);
    ylim([-.1 .3])
    ylabel('Mz/Mg')    
end

mkdir('MSfigs_QSmodel')
cd('MSfigs_QSmodel')
saveas(gca, 'MSfig_Myaw_timelines_QSmodel_transONLY.fig')
saveas(gca, 'MSfig_Myaw_timelines_QSmodel_transONLY.png')
plot2svg(['MSfig_Myaw_timelines_QSmodel_transONLY.svg'])
cd ..

%% Mpitch
clear
% load('MOD_norm_data_QSmodel.mat')
load('norm_data_torque.mat')
load(pitchDB)

% R&L axes
varname=dir('flyVars*')
varname = varname.name;
load(varname)

M_R_trans_norm_PitchTorque_allNOfreq =  Mx_trans_norm_PitchTorque_allNOfreq*cosd(rot_axis_angle) + My_trans_norm_PitchTorque_allNOfreq*sind(rot_axis_angle);
M_L_trans_norm_PitchTorque_allNOfreq = -Mx_trans_norm_PitchTorque_allNOfreq*sind(rot_axis_angle) + My_trans_norm_PitchTorque_allNOfreq*cosd(rot_axis_angle);

M_R_transNrot_norm_PitchTorque_allNOfreq =  Mx_transNrot_norm_PitchTorque_allNOfreq*cosd(rot_axis_angle) + My_transNrot_norm_PitchTorque_allNOfreq*sind(rot_axis_angle);
M_L_transNrot_norm_PitchTorque_allNOfreq = -Mx_transNrot_norm_PitchTorque_allNOfreq*sind(rot_axis_angle) + My_transNrot_norm_PitchTorque_allNOfreq*cosd(rot_axis_angle);

% steady My
mod_now = 0;
mod_diff = abs(mod_value_PitchTorque - mod_now);
n = find(mod_diff==min(mod_diff));
val_steady = My_trans_norm_PitchTorque_allNOfreq(:,n);

% mod plots
mod_values = PitchMods;
% color_map = [0 0 0; 0 0 1; 0 .5 1];
mod_values = [-max(PitchMods) min(PitchMods) max(PitchMods)];
color_map = [1 0 0; 0 0 0; 0 1 1];

figure
for i = 1:length(mod_values)
    
    mod_now = mod_values(i);
    
    mod_diff = abs(mod_value_PitchTorque - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % Fx ALL
    val_now = Fx_trans_norm_PitchTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,1)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fx/Mg')
    
    % Fy ALL
    val_now = Fy_trans_norm_PitchTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,4)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fy/Mg')
    
    % Fz ALL
    val_now = Fz_trans_norm_PitchTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,7)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fz/Mg')
    
    % Mx ALL
    val_now = Mx_trans_norm_PitchTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,2)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mx/Mg')
    
    % My ALL
    val_now = My_trans_norm_PitchTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,5)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('My/Mg')
    
    % Mz ALL
    val_now = Mz_trans_norm_PitchTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,8)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mz/Mg')
    
    % M_R ALL
    val_now = M_R_trans_norm_PitchTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,3)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M_R/Mg')

    % M_L ALL
    val_now = M_L_trans_norm_PitchTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,6)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M_L/Mg')
    
    % -My ALL & mean
    val_now = My_trans_norm_PitchTorque_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    val_mean = mean(val_now);
    
    subplot(3,3,9)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('My/Mg')
    plot([0 1],[val_mean val_mean],'--','color',color_map(i,:),'linewidth',1)
    val_pp = csaps(t_now,val_now,.99999)
%     fnplt(val_pp);
    ylim([-1 1])
    ylabel('My/Mg')
end

mkdir('MSfigs_QSmodel')
cd('MSfigs_QSmodel')
saveas(gca, 'MSfig_Mpitch_timelines_QSmodel_transONLY.fig')
saveas(gca, 'MSfig_Mpitch_timelines_QSmodel_transONLY.png')
plot2svg(['MSfig_Mpitch_timelines_QSmodel_transONLY.svg'])
cd ..

%% M_R
clear
% load('MOD_norm_data_QSmodel.mat')
load('norm_data_torque.mat')
load(RaxisDB)

% R&L axes
varname=dir('flyVars*')
varname = varname.name;
load(varname)

M_R_trans_norm_TorqueAxisR_allNOfreq =  Mx_trans_norm_TorqueAxisR_allNOfreq*cosd(rot_axis_angle) + My_trans_norm_TorqueAxisR_allNOfreq*sind(rot_axis_angle);
M_L_trans_norm_TorqueAxisR_allNOfreq = -Mx_trans_norm_TorqueAxisR_allNOfreq*sind(rot_axis_angle) + My_trans_norm_TorqueAxisR_allNOfreq*cosd(rot_axis_angle);

M_R_transNrot_norm_TorqueAxisR_allNOfreq =  Mx_transNrot_norm_TorqueAxisR_allNOfreq*cosd(rot_axis_angle) + My_transNrot_norm_TorqueAxisR_allNOfreq*sind(rot_axis_angle);
M_L_transNrot_norm_TorqueAxisR_allNOfreq = -Mx_transNrot_norm_TorqueAxisR_allNOfreq*sind(rot_axis_angle) + My_transNrot_norm_TorqueAxisR_allNOfreq*cosd(rot_axis_angle);

% steady My
mod_now = 0;
mod_diff = abs(mod_value_TorqueAxisR - mod_now);
n = find(mod_diff==min(mod_diff));
val_steady = M_R_trans_norm_TorqueAxisR_allNOfreq(:,n);

% mod plots
mod_values = RaxisMods;
% color_map = [0 0 0; 0 0 1; 0 .5 1];
mod_values = [-max(RaxisMods) min(RaxisMods) max(RaxisMods)];
color_map = [1 0 0; 0 0 0; 0 1 1];

figure
for i = 1:length(mod_values)
    
    mod_now = mod_values(i);
    
    mod_diff = abs(mod_value_TorqueAxisR - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % Fx ALL
    val_now = Fx_trans_norm_TorqueAxisR_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,1)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fx/Mg')
    
    % Fy ALL
    val_now = Fy_trans_norm_TorqueAxisR_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,4)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fy/Mg')
    
    % Fz ALL
    val_now = Fz_trans_norm_TorqueAxisR_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,7)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fz/Mg')
    
    % Mx ALL
    val_now = Mx_trans_norm_TorqueAxisR_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,2)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mx/Mg')
    
    % My ALL
    val_now = My_trans_norm_TorqueAxisR_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,5)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('My/Mg')
    
    % Mz ALL
    val_now = Mz_trans_norm_TorqueAxisR_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,8)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mz/Mg')
    
    % M_R ALL
    val_now = M_R_trans_norm_TorqueAxisR_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,3)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M_R/Mg')

    % M_L ALL
    val_now = M_L_trans_norm_TorqueAxisR_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,6)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M_L/Mg')
    
    % M_R ALL & mean
    val_now = M_R_trans_norm_TorqueAxisR_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    val_mean = mean(val_now);
    
    subplot(3,3,9)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M_R/Mg')
    plot([0 1],[val_mean val_mean],'--','color',color_map(i,:),'linewidth',1)
    val_pp = csaps(t_now,val_now,.99999)
%     fnplt(val_pp);
    ylim([-1 1])
    ylabel('M_R/Mg')
end

mkdir('MSfigs_QSmodel')
cd('MSfigs_QSmodel')
saveas(gca, 'MSfig_M_R_timelines_QSmodel_transONLY.fig')
saveas(gca, 'MSfig_M_R_timelines_QSmodel_transONLY.png')
plot2svg(['MSfig_M_R_timelines_QSmodel_transONLY.svg'])
cd ..


