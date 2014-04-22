% make borf plots F&M timelines

clear
clc
close all

%% load norm data
load('MOD_norm_data.mat')

%% load norm data
% load('norm_data.mat')
% 
% FenhMods = Fenhanses/Fenhance_norm
% RollMods = rollaccels/rollaccel_norm
% % PitchMods = pitchaccels/pitchaccel_norm
% PitchMods = [-max(pitchaccels) 0 max(pitchaccels)]'/pitchaccel_norm
% 
% FDB = 'borf_db_Fenhance_NOcali_alldata.mat'
% rollDB = 'borf_db_F_roll_LR_NOcali_alldata.mat'
% pitchDB = 'borf_db_PitchAccel_NOcali_alldata.mat'
% 
% % FDB = 'borf_db_Fenhance_cali_alldata.mat'
% % rollDB = 'borf_db_F_roll_LR_cali_alldata.mat'
% % pitchDB = 'borf_db_PitchAccel_cali_alldata.mat'

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
color_map = [0 0 0; 0 0 1; 0 .5 1];

% rollAccelNorm 0to5
RollMods = [0:5]/rollaccel_norm;
mod_values = RollMods;
color_map = jet(length(mod_values))

figure
F_allNOfreq = sqrt(Fx_allNOfreq.^2 + Fz_allNOfreq.^2 + Fz_allNOfreq.^2);
M_allNOfreq = sqrt(Mx_allNOfreq.^2 + Mz_allNOfreq.^2 + Mz_allNOfreq.^2);
for i = 1:length(mod_values)
    mod_now = mod_values(i);
    
    mod_diff = abs(mod_value_allNOfreq - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    
    % Mx ALL
    val_now(:,:) = Mx_allNOfreq(:,:,n)/f_wb_steady^2/Lwing^5;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Mx')
    
end

mkdir('MSfigs')
cd('MSfigs')
saveas(gca, 'MSfig_Mroll_timelines_0to5.fig')
saveas(gca, 'MSfig_Mroll_timelines_0to5.png')
plot2svg(['MSfig_Mroll_timelines_0to5.svg'])
cd ..
