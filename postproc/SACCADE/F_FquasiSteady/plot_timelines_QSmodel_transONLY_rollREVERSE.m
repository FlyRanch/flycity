% make borf plots F&M timelines_QSmodel

% clear
% clc
% close all

%% load MOD&norm data
load('MOD_norm_data_QSmodel.mat')

%% Fenhance
load(FDB)

% mod_values = FenhMods;
% color_map = [.5 .5 .5; 0 .5 1; 0 .5 1];
mod_values = [min(FenhMods) max(FenhMods)];
color_map = [.5 .5 .5; 0 .5 1];

figure(1)
for i = 1:length(mod_values)
    mod_now = mod_values(i);
    
    mod_diff = abs(mod_value_Fenhance - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % Fx ALL
    val_now = Fx_trans_norm_Fenhance_all(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,1)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fx/Mg')
    
    % Fy ALL
    val_now = Fy_trans_norm_Fenhance_all(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,4)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fy/Mg')
    
    % Fz ALL
    val_now = Fz_trans_norm_Fenhance_all(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,7)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fz/Mg')
    
    % Mx ALL
    val_now = Mx_trans_norm_Fenhance_all(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,2)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mx/Mg')
    
    % My ALL
    val_now = My_trans_norm_Fenhance_all(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,5)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('My/Mg')
    
    % Mz ALL
    val_now = Mz_trans_norm_Fenhance_all(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,8)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mz/Mg')
    
    % Fnorm ALL
    val_now = sqrt(Fx_trans_norm_Fenhance_all(:,n).^2 + Fy_trans_norm_Fenhance_all(:,n).^2 + Fz_trans_norm_Fenhance_all(:,n).^2);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,3)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('F/Mg')

    % Mnorm ALL
    val_now = sqrt(Mx_trans_norm_Fenhance_all(:,n).^2 + My_trans_norm_Fenhance_all(:,n).^2 + Mz_trans_norm_Fenhance_all(:,n).^2);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,6)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M/Mg')
    
    % Fz ALL & mean
    val_now = -Fz_trans_norm_Fenhance_all(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    val_mean = mean(val_now);
    
    subplot(3,3,9)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('-Fz/Mg')
%     plot([0 1],[val_mean val_mean],'--','color',color_map(i,:),'linewidth',1)
%     val_pp = csaps(t_now,val_now,.99999)
%     fnplt(val_pp);
        ylim([-.5 3.5])
    ylabel('Fz/Mg')
end

% mkdir('MSfigs')
% cd('MSfigs')
% saveas(gca, 'MSfig_Fenhance_timelines_BORFnQSmodel_transONLY.fig')
% saveas(gca, 'MSfig_Fenhance_timelines_BORFnQSmodel_transONLY.png')
% plot2svg(['MSfig_Fenhance_timelines_BORFnQSmodel_transONLY.svg'])
% cd ..

%% Mroll
load(rollDB)

% steady Mx
mod_now = 0;
mod_diff = abs(mod_value_RollAccel - mod_now);
n = find(mod_diff==min(mod_diff));
val_steady = Mx_trans_norm_RollAccel_allNOfreq(:,n);

% mod plots
mod_values = RollMods;
% color_map = [.5 .5 .5; 0 0 1; 0 .5 1];
% mod_values = [min(RollMods) max(RollMods)];
color_map = [.5 .5 .5; 0 .5 1];

figure(2)
for i = 1:length(mod_values)
    
    mod_now = mod_values(i);
    
    mod_diff = abs(mod_value_RollAccel - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % Fx ALL
    val_now = Fx_trans_norm_RollAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,1)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fx/Mg')
    
    % Fy ALL
    val_now = -Fy_trans_norm_RollAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,4)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fy/Mg')
    
    % Fz ALL
    val_now = Fz_trans_norm_RollAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,7)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fz/Mg')
    
    % Mx ALL
    val_now = -Mx_trans_norm_RollAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,2)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mx/Mg')
    
    % My ALL
    val_now = My_trans_norm_RollAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,5)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('My/Mg')
    
    % Mz ALL
    val_now = -Mz_trans_norm_RollAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,8)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mz/Mg')
    
    % Fnorm ALL
    val_now = sqrt(Fx_trans_norm_RollAccel_allNOfreq(:,n).^2 + Fy_trans_norm_RollAccel_allNOfreq(:,n).^2 + Fz_trans_norm_RollAccel_allNOfreq(:,n).^2);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,3)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('F/Mg')

    % Mnorm ALL
    val_now = sqrt(Mx_trans_norm_RollAccel_allNOfreq(:,n).^2 + My_trans_norm_RollAccel_allNOfreq(:,n).^2 + Mz_trans_norm_RollAccel_allNOfreq(:,n).^2);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,6)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M*')
    
    % Mx ALL & mean
    val_now = -Mx_trans_norm_RollAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    val_mean = mean(val_now);
    
    subplot(3,3,9)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mx')
%     plot([0 1],[val_mean val_mean],'--','color',color_map(i,:),'linewidth',1)
%     val_pp = csaps(t_now,val_now,.99999)
%     fnplt(val_pp);
    ylim([-1 2])
    ylabel('Mx')    
end

% mkdir('MSfigs')
% cd('MSfigs')
% saveas(gca, 'MSfig_Mroll_timelines_BORFnQSmodel_transONLY.fig')
% saveas(gca, 'MSfig_Mroll_timelines_BORFnQSmodel_transONLY.png')
% plot2svg(['MSfig_Mroll_timelines_BORFnQSmodel_transONLY.svg'])
% cd ..

%% Myaw
load(yawDB)

% steady Mz
mod_now = 0;
mod_diff = abs(mod_value_YawAccel - mod_now);
n = find(mod_diff==min(mod_diff));
val_steady = Mz_trans_norm_YawAccel_allNOfreq(:,n);

% mod plots
mod_values = YawMods;
% color_map = [.5 .5 .5; 0 0 1; 0 .5 1];
% mod_values = [min(YawMods) max(YawMods)];
color_map = [.5 .5 .5; 0 .5 1];

figure(3)
for i = 1:length(mod_values)
    
    mod_now = mod_values(i);
    
    mod_diff = abs(mod_value_YawAccel - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % Fx ALL
    val_now = Fx_trans_norm_YawAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,1)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fx/Mg')
    
    % Fy ALL
    val_now = Fy_trans_norm_YawAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,4)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fy/Mg')
    
    % Fz ALL
    val_now = Fz_trans_norm_YawAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,7)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fz/Mg')
    
    % Mx ALL
    val_now = Mx_trans_norm_YawAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,2)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mx/Mg')
    
    % My ALL
    val_now = My_trans_norm_YawAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,5)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('My/Mg')
    
    % Mz ALL
    val_now = Mz_trans_norm_YawAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,8)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mz/Mg')
    
    % Fnorm ALL
    val_now = sqrt(Fx_trans_norm_YawAccel_allNOfreq(:,n).^2 + Fy_trans_norm_YawAccel_allNOfreq(:,n).^2 + Fz_trans_norm_YawAccel_allNOfreq(:,n).^2);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,3)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('F/Mg')

    % Mnorm ALL
    val_now = sqrt(Mx_trans_norm_YawAccel_allNOfreq(:,n).^2 + My_trans_norm_YawAccel_allNOfreq(:,n).^2 + Mz_trans_norm_YawAccel_allNOfreq(:,n).^2);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,6)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M*')
    
    % -Mz ALL & mean
    val_now = -Mz_trans_norm_YawAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    val_mean = mean(val_now);
    
    subplot(3,3,9)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mz')
%     plot([0 1],[val_mean val_mean],'--','color',color_map(i,:),'linewidth',1)
%     val_pp = csaps(t_now,val_now,.99999)
%     fnplt(val_pp);
    ylim([-.5 1.5])
    ylabel('-Mz')    
end
% 
% mkdir('MSfigs')
% cd('MSfigs')
% saveas(gca, 'MSfig_Myaw_timelines_BORFnQSmodel_transONLY.fig')
% saveas(gca, 'MSfig_Myaw_timelines_BORFnQSmodel_transONLY.png')
% plot2svg(['MSfig_Myaw_timelines_BORFnQSmodel_transONLY.svg'])
% cd ..

%% Mpitch
load(pitchDB)

% steady My
mod_now = 0;
mod_diff = abs(mod_value_PitchAccel - mod_now);
n = find(mod_diff==min(mod_diff));
val_steady = My_trans_norm_PitchAccel_allNOfreq(:,n);

% mod plots
mod_values = PitchMods;
% color_map = [.5 .5 .5; 0 0 1; 0 .5 1];
% mod_values = [min(PitchMods) max(PitchMods)];
color_map = [1 .5 0; .5 .5 .5; 0 .5 1];

figure(4)
for i = 1:length(mod_values)
    
    mod_now = mod_values(i);
    
    mod_diff = abs(mod_value_PitchAccel - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % Fx ALL
    val_now = Fx_trans_norm_PitchAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,1)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fx/Mg')
    
    % Fy ALL
    val_now = Fy_trans_norm_PitchAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,4)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fy/Mg')
    
    % Fz ALL
    val_now = Fz_trans_norm_PitchAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,7)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fz/Mg')
    
    % Mx ALL
    val_now = Mx_trans_norm_PitchAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,2)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mx/Mg')
    
    % My ALL
    val_now = My_trans_norm_PitchAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,5)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('My/Mg')
    
    % Mz ALL
    val_now = Mz_trans_norm_PitchAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,8)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mz/Mg')
    
    % Fnorm ALL
    val_now = sqrt(Fx_trans_norm_PitchAccel_allNOfreq(:,n).^2 + Fy_trans_norm_PitchAccel_allNOfreq(:,n).^2 + Fz_trans_norm_PitchAccel_allNOfreq(:,n).^2);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,3)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('F/Mg')

    % Mnorm ALL
    val_now = sqrt(Mx_trans_norm_PitchAccel_allNOfreq(:,n).^2 + My_trans_norm_PitchAccel_allNOfreq(:,n).^2 + Mz_trans_norm_PitchAccel_allNOfreq(:,n).^2);
    t_now = [0:1/(length(val_now)-1):1];
    
    subplot(3,3,6)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M*')
    
    % -My ALL & mean
    val_now = My_trans_norm_PitchAccel_allNOfreq(:,n);
    t_now = [0:1/(length(val_now)-1):1];
    val_mean = mean(val_now);
    
    subplot(3,3,9)
    plot(t_now,val_now,'--','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('My')
%     plot([0 1],[val_mean val_mean],'--','color',color_map(i,:),'linewidth',1)
%     val_pp = csaps(t_now,val_now,.99999)
%     fnplt(val_pp);
    ylim([-8 8])
    ylabel('-My')
end
% 
% mkdir('MSfigs')
% cd('MSfigs')
% saveas(gca, 'MSfig_Mpitch_timelines_BORFnQSmodel_transONLY.fig')
% saveas(gca, 'MSfig_Mpitch_timelines_BORFnQSmodel_transONLY.png')
% plot2svg(['MSfig_Mpitch_timelines_BORFnQSmodel_transONLY.svg'])
% cd ..


