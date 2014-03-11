% make borf plots F&M timelines

clear
clc
close all

%% load MOD&norm data
load('MOD_norm_data.mat')

%% Fenhance
load(FDB_all)

% mod_values = FenhMods;
% color_map = [0 0 0; 0 .5 1; 0 1 1];
mod_values = [min(FenhMods) max(FenhMods)];
color_map = [0 0 0; 0 1 1];

%%%%% !!!!!!!!!!! NEW WRONG MODS !!!!!!!!!!!!!!
% FenhMods = [0 .7]/Fenhance_norm
% mod_values = FenhMods;
% color_map = [0 0 0; 0 0 1];

figure
F_all = sqrt(Fx_all.^2 + Fz_all.^2 + Fz_all.^2);
M_all = sqrt(Mx_all.^2 + Mz_all.^2 + Mz_all.^2);
for i = 1:length(mod_values)
    mod_now = mod_values(i);
    
    mod_diff = abs(mod_value_all - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % Fx ALL
    val_now(:,:) = Fx_all(:,:,n)/Mg_fly;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,1)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fx/Mg')
    
    % Fy ALL
    val_now(:,:) = Fy_all(:,:,n)/Mg_fly;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,4)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fy/Mg')
    
    % Fz ALL
    val_now(:,:) = Fz_all(:,:,n)/Mg_fly;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,7)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fz/Mg')
    
    % Mx ALL
    val_now(:,:) = Mx_all(:,:,n)/f_wb_steady^2/Lwing^5;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,2)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mx')
    
    % My ALL
    val_now(:,:) = My_all(:,:,n)/f_wb_steady^2/Lwing^5;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,5)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('My')
    
    % Mz ALL
    val_now(:,:) = Mz_all(:,:,n)/f_wb_steady^2/Lwing^5;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,8)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mz')
    
    % Fnorm ALL
    val_now(:,:) = F_all(:,:,n)/Mg_fly;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,3)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('F/Mg')
    
    % Mnorm ALL
    val_now(:,:) = M_all(:,:,n)/f_wb_steady^2/Lwing^5;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,6)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M')
    
    % -Fz ALL
    val_now(:,:) = -Fz_all(:,:,n)/Mg_fly;
    val_mean = nanmean(val_now,2);
    val_mean_mean = nanmean(val_mean);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,9)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    plot([0 1],[val_mean_mean val_mean_mean],'--','color',color_map(i,:),'linewidth',1)
    val_pp = csaps(t_now,val_mean,.9999)
    fnplt(val_pp);
        ylim([-1 4])
    ylabel('-Fz/Mg')
end

mkdir('MSfigs')
cd('MSfigs')
saveas(gca, 'MSfig_Fenhance_timelines.fig')
saveas(gca, 'MSfig_Fenhance_timelines.png')
plot2svg(['MSfig_Fenhance_timelines.svg'])
cd ..

%% Mroll
load(rollDB_all)

% steady Mx
mod_now = 0;
mod_diff = abs(mod_value_allNOfreq - mod_now);
n = find(mod_diff==min(mod_diff));
val_now(:,:) = Mx_allNOfreq(:,:,n)/f_wb_steady^2/Lwing^5;
val_steady = nanmean(val_now,2);

% mod plots
mod_values = RollMods;
% color_map = [0 0 0; 0 0 1; 0 .5 1];
% mod_values = [min(RollMods) max(RollMods)];
color_map = [0 0 0; 0 1 1];

figure
F_allNOfreq = sqrt(Fx_allNOfreq.^2 + Fz_allNOfreq.^2 + Fz_allNOfreq.^2);
M_allNOfreq = sqrt(Mx_allNOfreq.^2 + Mz_allNOfreq.^2 + Mz_allNOfreq.^2);
for i = 1:length(mod_values)
    mod_now = mod_values(i);
    
    mod_diff = abs(mod_value_allNOfreq - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % Fx ALL
    val_now(:,:) = Fx_allNOfreq(:,:,n);
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,1)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fx/Mg')
    
    % Fy ALL
    val_now(:,:) = Fy_allNOfreq(:,:,n)/Mg_fly;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,4)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fy/Mg')
    
    % Fz ALL
    val_now(:,:) = Fz_allNOfreq(:,:,n)/Mg_fly;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,7)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fz/Mg')
    
    % Mx ALL
    val_now(:,:) = Mx_allNOfreq(:,:,n)/f_wb_steady^2/Lwing^5;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,2)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mx')
    
    % My ALL
    val_now(:,:) = My_allNOfreq(:,:,n)/f_wb_steady^2/Lwing^5;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,5)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('My')
    
    % Mz ALL
    val_now(:,:) = Mz_allNOfreq(:,:,n)/f_wb_steady^2/Lwing^5;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,8)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mz')
    
    % Fnorm ALL
    val_now(:,:) = F_allNOfreq(:,:,n)/Mg_fly;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,3)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('F/Mg')
    
    % Mnorm ALL
    val_now(:,:) = M_allNOfreq(:,:,n)/f_wb_steady^2/Lwing^5;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,6)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M')
    
    % Mx ALL
    val_now(:,:) = Mx_allNOfreq(:,:,n)/f_wb_steady^2/Lwing^5;
    val_mean = nanmean(val_now,2);
    val_mean_mean = nanmean(val_mean);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,9)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    plot([0 1],[val_mean_mean val_mean_mean],'--','color',color_map(i,:),'linewidth',1)
    val_pp = csaps(t_now,val_mean,.9999)
    fnplt(val_pp);
    ylabel('Mx')
    ylim([-1 4])
end

mkdir('MSfigs')
cd('MSfigs')
saveas(gca, 'MSfig_Mroll_timelines.fig')
saveas(gca, 'MSfig_Mroll_timelines.png')
plot2svg(['MSfig_Mroll_timelines.svg'])
cd ..

figure
for i = 1:length(mod_values)
    mod_now = mod_values(i);
    
    mod_diff = abs(mod_value_allNOfreq - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % Mx-Mx_steady ALL
    val_now(:,:) = Mx_allNOfreq(:,:,n)/f_wb_steady^2/Lwing^5;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
%     subplot(3,3,9)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    val_pp = csaps(t_now,val_mean,.9999)
    fnplt(val_pp);
    hold on
    ylabel('Mx')
end

%% Mpitch
load(pitchDB_all)

% steady Mx
mod_now = 0;
mod_diff = abs(mod_value_allNOfreq - mod_now);
n = find(mod_diff==min(mod_diff));
val_now(:,:) = My_allNOfreq(:,:,n)/f_wb_steady^2/Lwing^5;
val_steady = nanmean(val_now(:));

% mod plots
mod_values = PitchMods;
color_map = [1 0 0; 0 0 0; 0 1 1];

figure
F_allNOfreq = sqrt(Fx_allNOfreq.^2 + Fz_allNOfreq.^2 + Fz_allNOfreq.^2);
M_allNOfreq = sqrt(Mx_allNOfreq.^2 + Mz_allNOfreq.^2 + Mz_allNOfreq.^2);
for i = 1:length(mod_values)
    mod_now = mod_values(i);
    
    mod_diff = abs(mod_value_allNOfreq - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % Fx ALL
    val_now(:,:) = Fx_allNOfreq(:,:,n)/Mg_fly;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,1)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fx/Mg')
    
    % Fy ALL
    val_now(:,:) = Fy_allNOfreq(:,:,n)/Mg_fly;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,4)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fy/Mg')
    
    % Fz ALL
    val_now(:,:) = Fz_allNOfreq(:,:,n)/Mg_fly;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,7)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fz/Mg')
    
    % Mx ALL
    val_now(:,:) = Mx_allNOfreq(:,:,n)/f_wb_steady^2/Lwing^5;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,2)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mx')
    
    % My ALL
    val_now(:,:) = My_allNOfreq(:,:,n)/f_wb_steady^2/Lwing^5;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,5)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('My')
    
    % Mz ALL
    val_now(:,:) = Mz_allNOfreq(:,:,n)/f_wb_steady^2/Lwing^5;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,8)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mz')
    
    % Fnorm ALL
    val_now(:,:) = F_allNOfreq(:,:,n)/Mg_fly;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,3)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('F/Mg')
    
    % Mnorm ALL
    val_now(:,:) = M_allNOfreq(:,:,n)/f_wb_steady^2/Lwing^5;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,6)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('M')
    
        
    % My-My_steady_mean ALL
    val_now(:,:) = My_allNOfreq(:,:,n)/f_wb_steady^2/Lwing^5;
    val_mean = nanmean(val_now,2) - val_steady;
    val_mean_mean = nanmean(val_mean);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    subplot(3,3,9)
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    plot([0 1],[val_mean_mean val_mean_mean],'--','color',color_map(i,:),'linewidth',1)
    val_pp = csaps(t_now,val_mean,.9999)
    fnplt(val_pp);
    ylabel('Mx-steady')
    ylim([-8 8])
end

mkdir('MSfigs')
cd('MSfigs')
saveas(gca, 'MSfig_Mpitch_timelines.fig')
saveas(gca, 'MSfig_Mpitch_timelines.png')
plot2svg(['MSfig_Mpitch_timelines.svg'])    
cd ..

