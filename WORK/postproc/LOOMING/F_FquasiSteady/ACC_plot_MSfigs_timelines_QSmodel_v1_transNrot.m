% make borf plots F&M timelines_QSmodel

clc
close all


%% Fenhance
clear
load('MOD_norm_data_QSmodel.mat')
load(FDB)

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
    val_now = Fx_transNrot_norm_Fenhance_all(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,1)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fx/Mg')
    
    % Fy ALL
    val_now = Fy_transNrot_norm_Fenhance_all(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,4)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fy/Mg')
    
    % Fz ALL
    val_now = Fz_transNrot_norm_Fenhance_all(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,7)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fz/Mg')
    
    % Mx ALL
    val_now = Mx_transNrot_norm_Fenhance_all(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,2)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mx/Mg')
    
    % My ALL
    val_now = My_transNrot_norm_Fenhance_all(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,5)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('My/Mg')
    
    % Mz ALL
    val_now = Mz_transNrot_norm_Fenhance_all(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,8)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mz/Mg')
    
    % Fnorm ALL
    val_now = sqrt(Fx_transNrot_norm_Fenhance_all(:,n).^2 + Fy_transNrot_norm_Fenhance_all(:,n).^2 + Fz_transNrot_norm_Fenhance_all(:,n).^2);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,3)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('F/Mg')

    % Mnorm ALL
    val_now = sqrt(Mx_transNrot_norm_Fenhance_all(:,n).^2 + My_transNrot_norm_Fenhance_all(:,n).^2 + Mz_transNrot_norm_Fenhance_all(:,n).^2);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,6)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M/Mg')
    
    % Fz ALL & mean
    val_now = -Fz_transNrot_norm_Fenhance_all(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    val_mean = mean(val_now);
    
    subplot(3,3,9)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('-Fz/Mg')
    plot([0 1],[val_mean val_mean],'--','color',color_map(i,:),'linewidth',1)
    val_pp = csaps(t_now,val_now,.99999)
%     fnplt(val_pp);
        ylim([-.5 5])
    ylabel('Fz/Mg')
end

mkdir('MSfigs')
cd('MSfigs')
saveas(gca, 'MSfig_Fenhance_timelines_QSmodel_transNrot.fig')
saveas(gca, 'MSfig_Fenhance_timelines_QSmodel_transNrot.png')
plot2svg(['MSfig_Fenhance_timelines_QSmodel_transNrot.svg'])
cd ..

%% Mroll
clear
load('MOD_norm_data_QSmodel.mat')
load(rollDB)

% steady Mx
mod_now = 0;
mod_diff = abs(mod_value_RollAccel - mod_now);
n = find(mod_diff==min(mod_diff));
val_steady = Mx_transNrot_norm_RollAccel_allNOfreq(:,n);

% mod plots
mod_values = RollMods;
% color_map = [0 0 0; 0 0 1; 0 .5 1];
% mod_values = [min(RollMods) max(RollMods)];
color_map = [0 0 0; 0 1 1];

figure
for i = 1:length(mod_values)
    
    mod_now = mod_values(i);
    
    mod_diff = abs(mod_value_RollAccel - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % Fx ALL
    val_now = Fx_transNrot_norm_RollAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,1)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fx/Mg')
    
    % Fy ALL
    val_now = Fy_transNrot_norm_RollAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,4)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fy/Mg')
    
    % Fz ALL
    val_now = Fz_transNrot_norm_RollAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,7)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fz/Mg')
    
    % Mx ALL
    val_now = Mx_transNrot_norm_RollAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,2)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mx/Mg')
    
    % My ALL
    val_now = My_transNrot_norm_RollAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,5)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('My/Mg')
    
    % Mz ALL
    val_now = Mz_transNrot_norm_RollAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,8)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mz/Mg')
    
    % Fnorm ALL
    val_now = sqrt(Fx_transNrot_norm_RollAccel_allNOfreq(:,n).^2 + Fy_transNrot_norm_RollAccel_allNOfreq(:,n).^2 + Fz_transNrot_norm_RollAccel_allNOfreq(:,n).^2);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,3)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('F/Mg')

    % Mnorm ALL
    val_now = sqrt(Mx_transNrot_norm_RollAccel_allNOfreq(:,n).^2 + My_transNrot_norm_RollAccel_allNOfreq(:,n).^2 + Mz_transNrot_norm_RollAccel_allNOfreq(:,n).^2);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,6)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M*')
    
    % Mx ALL & mean
    val_now = -Mx_transNrot_norm_RollAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    val_mean = mean(val_now);
    
    subplot(3,3,9)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mx')
    plot([0 1],[val_mean val_mean],'--','color',color_map(i,:),'linewidth',1)
    val_pp = csaps(t_now,val_now,.99999)
%     fnplt(val_pp);
    ylim([-.1 .3])
    ylabel('Mx')    
end

mkdir('MSfigs')
cd('MSfigs')
saveas(gca, 'MSfig_Mroll_timelines_QSmodel_transNrot.fig')
saveas(gca, 'MSfig_Mroll_timelines_QSmodel_transNrot.png')
plot2svg(['MSfig_Mroll_timelines_QSmodel_transNrot.svg'])
cd ..

%% Myaw
clear
load('MOD_norm_data_QSmodel.mat')
load(yawDB)

% steady Mz
mod_now = 0;
mod_diff = abs(mod_value_YawAccel - mod_now);
n = find(mod_diff==min(mod_diff));
val_steady = Mz_transNrot_norm_YawAccel_allNOfreq(:,n);

% mod plots
mod_values = YawMods;
% color_map = [0 0 0; 0 0 1; 0 .5 1];
% mod_values = [min(YawMods) max(YawMods)];
color_map = [0 0 0; 0 1 1];

figure
for i = 1:length(mod_values)
    
    mod_now = mod_values(i);
    
    mod_diff = abs(mod_value_YawAccel - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % Fx ALL
    val_now = Fx_transNrot_norm_YawAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,1)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fx/Mg')
    
    % Fy ALL
    val_now = Fy_transNrot_norm_YawAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,4)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fy/Mg')
    
    % Fz ALL
    val_now = Fz_transNrot_norm_YawAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,7)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fz/Mg')
    
    % Mx ALL
    val_now = Mx_transNrot_norm_YawAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,2)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mx/Mg')
    
    % My ALL
    val_now = My_transNrot_norm_YawAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,5)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('My/Mg')
    
    % Mz ALL
    val_now = Mz_transNrot_norm_YawAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,8)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mz/Mg')
    
    % Fnorm ALL
    val_now = sqrt(Fx_transNrot_norm_YawAccel_allNOfreq(:,n).^2 + Fy_transNrot_norm_YawAccel_allNOfreq(:,n).^2 + Fz_transNrot_norm_YawAccel_allNOfreq(:,n).^2);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,3)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('F/Mg')

    % Mnorm ALL
    val_now = sqrt(Mx_transNrot_norm_YawAccel_allNOfreq(:,n).^2 + My_transNrot_norm_YawAccel_allNOfreq(:,n).^2 + Mz_transNrot_norm_YawAccel_allNOfreq(:,n).^2);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,6)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M*')
    
    % -Mz ALL & mean
    val_now = -Mz_transNrot_norm_YawAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    val_mean = mean(val_now);
    
    subplot(3,3,9)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mz')
    plot([0 1],[val_mean val_mean],'--','color',color_map(i,:),'linewidth',1)
    val_pp = csaps(t_now,val_now,.99999)
%     fnplt(val_pp);
    ylim([-.15 .2])
    ylabel('-Mz')    
end

mkdir('MSfigs')
cd('MSfigs')
saveas(gca, 'MSfig_Myaw_timelines_QSmodel_transNrot.fig')
saveas(gca, 'MSfig_Myaw_timelines_QSmodel_transNrot.png')
plot2svg(['MSfig_Myaw_timelines_QSmodel_transNrot.svg'])
cd ..

%% Mpitch
clear
load('MOD_norm_data_QSmodel.mat')
load(pitchDB)

% steady My
mod_now = 0;
mod_diff = abs(mod_value_PitchAccel - mod_now);
n = find(mod_diff==min(mod_diff));
val_steady = My_transNrot_norm_PitchAccel_allNOfreq(:,n);

% mod plots
mod_values = PitchMods;
% color_map = [0 0 0; 0 0 1; 0 .5 1];
% mod_values = [min(PitchMods) max(PitchMods)];
color_map = [1 0 0; 0 0 0; 0 1 1];

figure
for i = 1:length(mod_values)
    
    mod_now = mod_values(i);
    
    mod_diff = abs(mod_value_PitchAccel - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % Fx ALL
    val_now = Fx_transNrot_norm_PitchAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,1)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fx/Mg')
    
    % Fy ALL
    val_now = Fy_transNrot_norm_PitchAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,4)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fy/Mg')
    
    % Fz ALL
    val_now = Fz_transNrot_norm_PitchAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,7)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fz/Mg')
    
    % Mx ALL
    val_now = Mx_transNrot_norm_PitchAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,2)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mx/Mg')
    
    % My ALL
    val_now = My_transNrot_norm_PitchAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,5)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('My/Mg')
    
    % Mz ALL
    val_now = Mz_transNrot_norm_PitchAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,8)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mz/Mg')
    
    % Fnorm ALL
    val_now = sqrt(Fx_transNrot_norm_PitchAccel_allNOfreq(:,n).^2 + Fy_transNrot_norm_PitchAccel_allNOfreq(:,n).^2 + Fz_transNrot_norm_PitchAccel_allNOfreq(:,n).^2);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,3)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('F/Mg')

    % Mnorm ALL
    val_now = sqrt(Mx_transNrot_norm_PitchAccel_allNOfreq(:,n).^2 + My_transNrot_norm_PitchAccel_allNOfreq(:,n).^2 + Mz_transNrot_norm_PitchAccel_allNOfreq(:,n).^2);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,6)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M*')
    
    % -My ALL & mean
    val_now = My_transNrot_norm_PitchAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    val_mean = mean(val_now);
    
    subplot(3,3,9)
    plot(t_now,val_now,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('My')
    plot([0 1],[val_mean val_mean],'--','color',color_map(i,:),'linewidth',1)
    val_pp = csaps(t_now,val_now,.99999)
%     fnplt(val_pp);
    ylim([-1 1])
    ylabel('-My')
end

mkdir('MSfigs')
cd('MSfigs')
saveas(gca, 'MSfig_Mpitch_timelines_QSmodel_transNrot.fig')
saveas(gca, 'MSfig_Mpitch_timelines_QSmodel_transNrot.png')
plot2svg(['MSfig_Mpitch_timelines_QSmodel_transNrot.svg'])
cd ..


