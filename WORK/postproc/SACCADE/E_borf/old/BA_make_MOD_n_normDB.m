% make normDB from WBmod files

clear
clc
close all

if exist('norm_data.mat') ~= 2

    name = dir('WBdataset_all_steadyNmods*')

    name = name.name;
    load(name)

    % save norm data
    save('norm_data.mat','rollaccel_norm','pitchaccel_norm','yawaccel_norm','Fenhance_norm',...
        'f_wb_steady')
end


%% make MOD&normDB
clear
clc
close all

butter_cut = 7; %Hz
butter_n = 5; % order

cali_on = 1
% cali_on = 0

% load norm data
load('norm_data.mat')
% 
% % mod locations
% Fenhanses = [0 .6]'
% rollaccels = [0 3]'
% % pitchaccels = [0 1.5 3]'
% pitchaccels = [0 2]'
% yawaccels = [0 2]'

FenhMods = Fenhanses/Fenhance_norm
RollMods = rollaccels/rollaccel_norm
% PitchMods = pitchaccels/pitchaccel_norm
PitchMods = [-max(pitchaccels) 0 max(pitchaccels)]'/pitchaccel_norm
YawMods = yawaccels/yawaccel_norm

if cali_on == 1

    FDB_mean = 'borf_db_Fenhance_cali_means.mat'
    rollDB_mean = 'borf_db_F_roll_LR_cali_means.mat'
    pitchDB_mean = 'borf_db_PitchAccel_cali_means.mat'
    yawDB_mean = 'borf_db_F_yaw_LR_cali_means.mat'

    FDB_all = 'borf_db_Fenhance_cali_alldata.mat'
    rollDB_all = 'borf_db_F_roll_LR_cali_alldata.mat'
    pitchDB_all = 'borf_db_PitchAccel_cali_alldata.mat'
    yawDB_all = 'borf_db_F_yaw_LR_cali_alldata.mat'

    rollDB_LnRsym = 'borf_db_F_roll_LR_cali_alldata_LnRsym.mat'
    yawDB_LnRsym = 'borf_db_F_yaw_LR_cali_alldata_LnRsym.mat'

else

    FDB_mean = 'borf_db_Fenhance_NOcali_means.mat'
    rollDB_mean = 'borf_db_F_roll_LR_NOcali_means.mat'
    pitchDB_mean = 'borf_db_PitchAccel_NOcali_means.mat'
    yawDB_mean = 'borf_db_F_yaw_LR_NOcali_means.mat'
    
    FDB_all = 'borf_db_Fenhance_NOcali_alldata.mat'
    rollDB_all = 'borf_db_F_roll_LR_NOcali_alldata.mat'
    pitchDB_all = 'borf_db_PitchAccel_NOcali_alldata.mat'
    yawDB_all = 'borf_db_F_yaw_LR_NOcali_alldata.mat'
    
    rollDB_LnRsym = 'borf_db_F_roll_LR_NOcali_alldata_LnRsym.mat'
    yawDB_LnRsym = 'borf_db_F_yaw_LR_NOcali_alldata_LnRsym.mat'

end

save('MOD_norm_data.mat')
