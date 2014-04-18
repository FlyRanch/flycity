% make borf plots F&M timelines

clear
clc
close all

%% load norm data
load('MOD_norm_data.mat')
% FDB_mean = 'borf_db_Fenhance_NOcali_INCvel_means.mat'
% FDB_all = 'borf_db_Fenhance_NOcali_INCvel_alldata.mat'

%% load norm data
% load('norm_data.mat')
% 
% FenhMods = Fenhanses/Fenhance_norm
% RollMods = rollaccels/rollaccel_norm
% % PitchMods = pitchaccels/pitchaccel_norm
% PitchMods = [-max(pitchaccels) 0 max(pitchaccels)]'/pitchaccel_norm
% 
% FDB = 'borf_db_Fenhance_NOcali_alldata.mat'
% rollDB = 'borf_db_F_roll_LR_NOcali_means.mat'
% pitchDB = 'borf_db_PitchAccel_NOcali_means.mat'
% 
% FDB = 'borf_db_Fenhance_cali_alldata.mat'
% rollDB = 'borf_db_F_roll_LR_cali_alldata.mat'
% pitchDB = 'borf_db_PitchAccel_cali_alldata.mat'

%% Fenhance
load(FDB_all)

mod_values = FenhMods;

%%%%% !!!!!!!!!!! NEW WRONG MODS !!!!!!!!!!!!!!
FenhMods = [0:.2:1]/Fenhance_norm
mod_values = FenhMods;
color_map = jet(length(FenhMods));

figure
F_all = sqrt(Fx_all.^2 + Fz_all.^2 + Fz_all.^2);
M_all = sqrt(Mx_all.^2 + Mz_all.^2 + Mz_all.^2);
for i = 1:length(mod_values)
    mod_now = mod_values(i);
    
    mod_diff = abs(mod_value_all - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % Fz ALL
    val_now(:,:) = Fz_all(:,:,n)/Mg_fly;
    val_mean = nanmean(val_now,2);
    val_std = nanstd(val_now')';
    
    val_mean = val_mean(isnan(val_mean)==0);
    val_std = val_std(isnan(val_std)==0);
    t_now = [0:1/(length(val_mean)-1):1];
    
    plot(t_now,val_mean,'-','color',color_map(i,:),'linewidth',1)
    hold on
    ylabel('Fz/Mg')
end

mkdir('MSfigs')
cd('MSfigs')
saveas(gca, 'MSfig_Fenhance_timelines_1to2g.fig')
saveas(gca, 'MSfig_Fenhance_timelines_1to2g.png')
plot2svg(['MSfig_Fenhance_timelines_1to2g.svg'])

% saveas(gca, 'MSfig_Fenhance_timelines_1to2g_INCvel.fig')
% saveas(gca, 'MSfig_Fenhance_timelines_1to2g_INCvel.png')
% plot2svg(['MSfig_Fenhance_timelines_1to2g_INCvel.svg'])
cd ..
