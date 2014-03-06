% make borf dataset Fenhance

clear
clc
close all
warning off

% save_name = 'F_roll_DLR';
save_name = 'F_roll_LR';

%% cali data
% addpath('/home/florian/DATA/flighttracker/looming_FlyCity_borf')
% cali_file = 'cali_matrix_interp_08132013.mat'
% load(cali_file)

%% settings
seq_length_max = 15000;
Nwb_total = 7;
Nwb = 4;
wb_start = 4;
wb_stop = wb_start + Nwb-1;

% robo data
Lrobo = .23;
crobo = .065;
ARrobo = Lrobo/crobo;

% fly data
Lwing = .00299; %m
ARwing_fly  = 3.1579; % wing aspect ratio of fly (based on Elzinga's R/c)

% Mfly = 113.9246*Lwing^3.008; % kg
% % M = 1e-6*.1078*(Lwing*1e3)^3.008; % kg
% 
% % 10 fresh fem & 10 fresh male, by MD
% Mfly = 2.3e-6;
% 
% % 7 fresh fem & 7 fresh male, by FTM
% Mfly = 2.3e-6;

% 8 fem & 7 male, starved for 24h, by FTM. 50%male/female
Mfly = 1.8e-6;
Mg_fly = Mfly*9.81;

% fluid data
v_air = 1.568e-5; %kinematic viscocity in m^2/s
v_oil = 115.0*1e-6; %kinematic viscocity in m^2/s
rho_air = 1.22521; %kg/m^3
rho_oil = 880.0; %kg/m^3

%% FTM scale
% cwing = Lwing/ARrobo;
% f_robo2fly = (v_air/v_oil) * (Lrobo/Lwing)^2;
% F_robo2fly = (rho_air/rho_oil) * (v_air/v_oil)^2;
% M_robo2fly = F_robo2fly * (Lwing/Lrobo);


%% Elzinga scaling
% l_scale_up = crobo/c_fly;
c_fly = Lwing/ARwing_fly; %mean chord of actual fly in m
f_robo2fly = (v_air*crobo*Lrobo)/(c_fly*Lwing*v_oil);
F_robo2fly = (rho_air/rho_oil)*((c_fly/crobo)^4)*((f_robo2fly)^2);
M_robo2fly = (rho_air/rho_oil)*((c_fly/crobo)^5)*((f_robo2fly)^2);

%% colormap
imax = 21;

% black2blue
colormap_black2blue(1:imax,1) = 0;
colormap_black2blue(1:imax,2) = [0:.5/(imax-1):.5];
colormap_black2blue(1:imax,3) = [0:1/(imax-1):1];

%% make db Fenhance

%% LEFT SYM
side = 'left';

% file_name = [save_name,'_all_sym_',side,'_'];
% makeDB_allMODs_NOcali_WBs

file_name = [save_name,'_all_const_freq_sym_',side,'_'];
makeDB_allNOfreqMODs_NOcali_WBs

% file_name = [save_name,'_freq_sym_',side,'_'];
% makeDB_freqMODs_NOcali_WBs

file_name = [save_name,'_stroke_sym_',side,'_'];
makeDB_strokeMODs_NOcali_WBs

file_name = [save_name,'_rot_sym_',side,'_'];
makeDB_pitchMODs_NOcali_WBs

file_name = [save_name,'_dev_sym_',side,'_'];
makeDB_devMODs_NOcali_WBs


%% store data

% allNOfreq
Fx_mean_allNOfreq_left = Fx_mean_allNOfreq;
Fy_mean_allNOfreq_left = Fy_mean_allNOfreq;
Fz_mean_allNOfreq_left = Fz_mean_allNOfreq;
Mx_mean_allNOfreq_left = Mx_mean_allNOfreq;
My_mean_allNOfreq_left = My_mean_allNOfreq;
Mz_mean_allNOfreq_left = Mz_mean_allNOfreq;

Fx_allNOfreq_left = Fx_allNOfreq;
Fy_allNOfreq_left = Fy_allNOfreq;
Fz_allNOfreq_left = Fz_allNOfreq;
Mx_allNOfreq_left = Mx_allNOfreq;
My_allNOfreq_left = My_allNOfreq;
Mz_allNOfreq_left = Mz_allNOfreq;

t_allNOfreq_left = t_allNOfreq;
mod_value_allNOfreq_left = mod_value_allNOfreq;

stroke_L_allNOfreq_left = stroke_L_allNOfreq;
pitch_L_allNOfreq_left = pitch_L_allNOfreq;
dev_L_allNOfreq_left = dev_L_allNOfreq;

stroke_R_allNOfreq_left = stroke_R_allNOfreq;
pitch_R_allNOfreq_left = pitch_R_allNOfreq;
dev_R_allNOfreq_left = dev_R_allNOfreq;

% stroke
Fx_mean_stroke_left = Fx_mean_stroke;
Fy_mean_stroke_left = Fy_mean_stroke;
Fz_mean_stroke_left = Fz_mean_stroke;
Mx_mean_stroke_left = Mx_mean_stroke;
My_mean_stroke_left = My_mean_stroke;
Mz_mean_stroke_left = Mz_mean_stroke;

Fx_stroke_left = Fx_stroke;
Fy_stroke_left = Fy_stroke;
Fz_stroke_left = Fz_stroke;
Mx_stroke_left = Mx_stroke;
My_stroke_left = My_stroke;
Mz_stroke_left = Mz_stroke;

t_stroke_left = t_stroke;
mod_value_stroke_left = mod_value_stroke;

stroke_L_stroke_left = stroke_L_stroke;
pitch_L_stroke_left = pitch_L_stroke;
dev_L_stroke_left = dev_L_stroke;

stroke_R_stroke_left = stroke_R_stroke;
pitch_R_stroke_left = pitch_R_stroke;
dev_R_stroke_left = dev_R_stroke;

% pitch
Fx_mean_pitch_left = Fx_mean_pitch;
Fy_mean_pitch_left = Fy_mean_pitch;
Fz_mean_pitch_left = Fz_mean_pitch;
Mx_mean_pitch_left = Mx_mean_pitch;
My_mean_pitch_left = My_mean_pitch;
Mz_mean_pitch_left = Mz_mean_pitch;

Fx_pitch_left = Fx_pitch;
Fy_pitch_left = Fy_pitch;
Fz_pitch_left = Fz_pitch;
Mx_pitch_left = Mx_pitch;
My_pitch_left = My_pitch;
Mz_pitch_left = Mz_pitch;

t_pitch_left = t_pitch;
mod_value_pitch_left = mod_value_pitch;

stroke_L_pitch_left = stroke_L_pitch;
pitch_L_pitch_left = pitch_L_pitch;
dev_L_pitch_left = dev_L_pitch;

stroke_R_pitch_left = stroke_R_pitch;
pitch_R_pitch_left = pitch_R_pitch;
dev_R_pitch_left = dev_R_pitch;

% dev
Fx_mean_dev_left = Fx_mean_dev;
Fy_mean_dev_left = Fy_mean_dev;
Fz_mean_dev_left = Fz_mean_dev;
Mx_mean_dev_left = Mx_mean_dev;
My_mean_dev_left = My_mean_dev;
Mz_mean_dev_left = Mz_mean_dev;

Fx_dev_left = Fx_dev;
Fy_dev_left = Fy_dev;
Fz_dev_left = Fz_dev;
Mx_dev_left = Mx_dev;
My_dev_left = My_dev;
Mz_dev_left = Mz_dev;

t_dev_left = t_dev;
mod_value_dev_left = mod_value_dev;

stroke_L_dev_left = stroke_L_dev;
pitch_L_dev_left = pitch_L_dev;
dev_L_dev_left = dev_L_dev;

stroke_R_dev_left = stroke_R_dev;
pitch_R_dev_left = pitch_R_dev;
dev_R_dev_left = dev_R_dev;

%% NO mirror
side = 'right';

% file_name = [save_name,'_all_sym_',side,'_'];
% makeDB_allMODs_NOcali_WBs

file_name = [save_name,'_all_const_freq_sym_',side,'_'];
makeDB_allNOfreqMODs_NOcali_WBs

% file_name = [save_name,'_freq_sym_',side,'_'];
% makeDB_freqMODs_NOcali_WBs

file_name = [save_name,'_stroke_sym_',side,'_'];
makeDB_strokeMODs_NOcali_WBs

file_name = [save_name,'_rot_sym_',side,'_'];
makeDB_pitchMODs_NOcali_WBs

file_name = [save_name,'_dev_sym_',side,'_'];
makeDB_devMODs_NOcali_WBs

%% store data

% allNOfreq
Fx_mean_allNOfreq_right = Fx_mean_allNOfreq;
Fy_mean_allNOfreq_right = Fy_mean_allNOfreq;
Fz_mean_allNOfreq_right = Fz_mean_allNOfreq;
Mx_mean_allNOfreq_right = Mx_mean_allNOfreq;
My_mean_allNOfreq_right = My_mean_allNOfreq;
Mz_mean_allNOfreq_right = Mz_mean_allNOfreq;

Fx_allNOfreq_right = Fx_allNOfreq;
Fy_allNOfreq_right = Fy_allNOfreq;
Fz_allNOfreq_right = Fz_allNOfreq;
Mx_allNOfreq_right = Mx_allNOfreq;
My_allNOfreq_right = My_allNOfreq;
Mz_allNOfreq_right = Mz_allNOfreq;

t_allNOfreq_right = t_allNOfreq;
mod_value_allNOfreq_right = mod_value_allNOfreq;

stroke_L_allNOfreq_right = stroke_L_allNOfreq;
pitch_L_allNOfreq_right = pitch_L_allNOfreq;
dev_L_allNOfreq_right = dev_L_allNOfreq;

stroke_R_allNOfreq_right = stroke_R_allNOfreq;
pitch_R_allNOfreq_right = pitch_R_allNOfreq;
dev_R_allNOfreq_right = dev_R_allNOfreq;

% stroke
Fx_mean_stroke_right = Fx_mean_stroke;
Fy_mean_stroke_right = Fy_mean_stroke;
Fz_mean_stroke_right = Fz_mean_stroke;
Mx_mean_stroke_right = Mx_mean_stroke;
My_mean_stroke_right = My_mean_stroke;
Mz_mean_stroke_right = Mz_mean_stroke;

Fx_stroke_right = Fx_stroke;
Fy_stroke_right = Fy_stroke;
Fz_stroke_right = Fz_stroke;
Mx_stroke_right = Mx_stroke;
My_stroke_right = My_stroke;
Mz_stroke_right = Mz_stroke;

t_stroke_right = t_stroke;
mod_value_stroke_right = mod_value_stroke;

stroke_L_stroke_right = stroke_L_stroke;
pitch_L_stroke_right = pitch_L_stroke;
dev_L_stroke_right = dev_L_stroke;

stroke_R_stroke_right = stroke_R_stroke;
pitch_R_stroke_right = pitch_R_stroke;
dev_R_stroke_right = dev_R_stroke;

% pitch
Fx_mean_pitch_right = Fx_mean_pitch;
Fy_mean_pitch_right = Fy_mean_pitch;
Fz_mean_pitch_right = Fz_mean_pitch;
Mx_mean_pitch_right = Mx_mean_pitch;
My_mean_pitch_right = My_mean_pitch;
Mz_mean_pitch_right = Mz_mean_pitch;

Fx_pitch_right = Fx_pitch;
Fy_pitch_right = Fy_pitch;
Fz_pitch_right = Fz_pitch;
Mx_pitch_right = Mx_pitch;
My_pitch_right = My_pitch;
Mz_pitch_right = Mz_pitch;

t_pitch_right = t_pitch;
mod_value_pitch_right = mod_value_pitch;

stroke_L_pitch_right = stroke_L_pitch;
pitch_L_pitch_right = pitch_L_pitch;
dev_L_pitch_right = dev_L_pitch;

stroke_R_pitch_right = stroke_R_pitch;
pitch_R_pitch_right = pitch_R_pitch;
dev_R_pitch_right = dev_R_pitch;

% dev
Fx_mean_dev_right = Fx_mean_dev;
Fy_mean_dev_right = Fy_mean_dev;
Fz_mean_dev_right = Fz_mean_dev;
Mx_mean_dev_right = Mx_mean_dev;
My_mean_dev_right = My_mean_dev;
Mz_mean_dev_right = Mz_mean_dev;

Fx_dev_right = Fx_dev;
Fy_dev_right = Fy_dev;
Fz_dev_right = Fz_dev;
Mx_dev_right = Mx_dev;
My_dev_right = My_dev;
Mz_dev_right = Mz_dev;

t_dev_right = t_dev;
mod_value_dev_right = mod_value_dev;

stroke_L_dev_right = stroke_L_dev;
pitch_L_dev_right = pitch_L_dev;
dev_L_dev_right = dev_L_dev;

stroke_R_dev_right = stroke_R_dev;
pitch_R_dev_right = pitch_R_dev;
dev_R_dev_right = dev_R_dev;

%% merge data

% allNOfreq
diff = abs(mod_value_allNOfreq_right - mod_value_allNOfreq_left);
if nanmax(diff(:)) == 0
    
    Fx_mean_allNOfreq = (Fx_mean_allNOfreq_right + Fx_mean_allNOfreq_left)/2;
    Fy_mean_allNOfreq = (Fy_mean_allNOfreq_right + Fy_mean_allNOfreq_left)/2;
    Fz_mean_allNOfreq = (Fz_mean_allNOfreq_right + Fz_mean_allNOfreq_left)/2;

    Mx_mean_allNOfreq = (Mx_mean_allNOfreq_right + Mx_mean_allNOfreq_left)/2;
    My_mean_allNOfreq = (My_mean_allNOfreq_right + My_mean_allNOfreq_left)/2;
    Mz_mean_allNOfreq = (Mz_mean_allNOfreq_right + Mz_mean_allNOfreq_left)/2;
else
    error = 'allNOfreq'
    pause
end

% stroke
diff = abs(mod_value_stroke_right - mod_value_stroke_left);
if nanmax(diff(:)) == 0
    
    Fx_mean_stroke = (Fx_mean_stroke_right + Fx_mean_stroke_left)/2;
    Fy_mean_stroke = (Fy_mean_stroke_right + Fy_mean_stroke_left)/2;
    Fz_mean_stroke = (Fz_mean_stroke_right + Fz_mean_stroke_left)/2;

    Mx_mean_stroke = (Mx_mean_stroke_right + Mx_mean_stroke_left)/2;
    My_mean_stroke = (My_mean_stroke_right + My_mean_stroke_left)/2;
    Mz_mean_stroke = (Mz_mean_stroke_right + Mz_mean_stroke_left)/2;
else
    error = 'stroke'
    pause
end

% pitch
diff = abs(mod_value_pitch_right - mod_value_pitch_left);
if nanmax(diff(:)) == 0
    
    Fx_mean_pitch = (Fx_mean_pitch_right + Fx_mean_pitch_left)/2;
    Fy_mean_pitch = (Fy_mean_pitch_right + Fy_mean_pitch_left)/2;
    Fz_mean_pitch = (Fz_mean_pitch_right + Fz_mean_pitch_left)/2;

    Mx_mean_pitch = (Mx_mean_pitch_right + Mx_mean_pitch_left)/2;
    My_mean_pitch = (My_mean_pitch_right + My_mean_pitch_left)/2;
    Mz_mean_pitch = (Mz_mean_pitch_right + Mz_mean_pitch_left)/2;
else
    error = 'pitch'
    pause
end

% dev
diff = abs(mod_value_dev_right - mod_value_dev_left);
if nanmax(diff(:)) == 0
    
    Fx_mean_dev = (Fx_mean_dev_right + Fx_mean_dev_left)/2;
    Fy_mean_dev = (Fy_mean_dev_right + Fy_mean_dev_left)/2;
    Fz_mean_dev = (Fz_mean_dev_right + Fz_mean_dev_left)/2;

    Mx_mean_dev = (Mx_mean_dev_right + Mx_mean_dev_left)/2;
    My_mean_dev = (My_mean_dev_right + My_mean_dev_left)/2;
    Mz_mean_dev = (Mz_mean_dev_right + Mz_mean_dev_left)/2;
else
    error = 'dev'
    pause
end

stroke_L_allNOfreq = stroke_L_allNOfreq_left;
pitch_L_allNOfreq = pitch_L_allNOfreq_left;
dev_L_allNOfreq = dev_L_allNOfreq_left;

stroke_L_stroke = stroke_L_stroke_left;
pitch_L_stroke = pitch_L_stroke_left;
dev_L_stroke = dev_L_stroke_left;

stroke_L_pitch = stroke_L_pitch_left;
pitch_L_pitch = pitch_L_pitch_left;
dev_L_pitch = dev_L_pitch_left;

stroke_L_dev = stroke_L_dev_left;
pitch_L_dev = pitch_L_dev_left;
dev_L_dev = dev_L_dev_left;

%% postprocess
MODmin = min([length(Mx_mean_allNOfreq) length(Mx_mean_stroke) length(Mx_mean_pitch) length(Mx_mean_dev)]);
Mx_steady = Mx_mean_allNOfreq(mod_value_allNOfreq == 0);
Mx_mean_sum = Mx_steady + (Mx_mean_stroke(1:MODmin)-Mx_steady) + (Mx_mean_pitch(1:MODmin)-Mx_steady) + (Mx_mean_dev(1:MODmin)-Mx_steady);
Mx_mean_sum_nopitch = Mx_steady + (Mx_mean_stroke(1:MODmin)-Mx_steady) + (Mx_mean_dev(1:MODmin)-Mx_steady);

%% plot data
cd ..
mkdir('figs_NOcali')
cd('figs_NOcali')
% plot_RollAccelNorm_LR_MODdata_mirror
plot_RollAccelNorm_LR_MODdata_LnRsym
% plot_RollAccelNorm_LR_MODdata

% save data
cd ..
clear exp_type ft kine mod_val t list i...
    t_now Fx_now Fy_now Fz_now Mx_now My_now Mz_now...
          Fx_robo Fy_robo Fz_robo Mx_robo My_robo Mz_robo...
          Fx_fly Fy_fly Fz_fly Mx_fly My_fly Mz_fly...
          Fx_cali_L Fy_cali_L Fz_cali_L Mx_cali_L My_cali_L Mz_cali_L...
          Fx_cali_R Fy_cali_R Fz_cali_R Mx_cali_R My_cali_R Mz_cali_R...
          Fx_cali_L_array Fy_cali_L_array Fz_cali_L_array Mx_cali_L_array My_cali_L_array Mz_cali_L_array...
          Fx_cali_R_array Fy_cali_R_array Fz_cali_R_array Mx_cali_R_array My_cali_R_array Mz_cali_R_array...
    stroke_L_now pitch_L_now dev_L_now stroke_R_now pitch_R_now dev_R_now...
    N_wb_now n_start n_stop frobo_now f_now

save(['borf_db_',save_name,'_NOcali_alldata_LnRsym.mat'])
% save(['borf_db_',save_name,'_NOcali_alldata_sym_',side,'.mat'])

clear Fx_all Fy_all Fz_all Mx_all My_all Mz_all...
    Fx_allNOfreq Fy_allNOfreq Fz_allNOfreq Mx_allNOfreq My_allNOfreq Mz_allNOfreq...
    Fx_freq Fy_freq Fz_freq Mx_freq My_freq Mz_freq...
    Fx_stroke Fy_stroke Fz_stroke Mx_stroke My_stroke Mz_stroke...
    Fx_pitch Fy_pitch Fz_pitch Mx_pitch My_pitch Mz_pitch...
    Fx_dev Fy_dev Fz_dev Mx_dev My_dev Mz_dev...
    stroke_L_all stroke_R_all pitch_L_all pitch_R_all dev_L_all dev_R_all...
    stroke_L_allNOfreq stroke_R_allNOfreq pitch_L_allNOfreq pitch_R_allNOfreq dev_L_allNOfreq dev_R_allNOfreq...
    stroke_L_freq stroke_R_freq pitch_L_freq pitch_R_freq dev_L_freq dev_R_freq...
    stroke_L_stroke stroke_R_stroke pitch_L_stroke pitch_R_stroke dev_L_stroke dev_R_stroke...
    stroke_L_pitch stroke_R_pitch pitch_L_pitch pitch_R_pitch dev_L_pitch dev_R_pitch...
    stroke_L_dev stroke_R_dev pitch_L_dev pitch_R_dev dev_L_dev dev_R_dev...
    t_all t_allNOfreq t_stroke t_pitch t_dev t_freq...
    Fx_all_right Fy_all_right Fz_all_right Mx_all_right My_all_right Mz_all_right...
    Fx_allNOfreq_right Fy_allNOfreq_right Fz_allNOfreq_right Mx_allNOfreq_right My_allNOfreq_right Mz_allNOfreq_right...
    Fx_freq_right Fy_freq_right Fz_freq_right Mx_freq_right My_freq_right Mz_freq_right...
    Fx_stroke_right Fy_stroke_right Fz_stroke_right Mx_stroke_right My_stroke_right Mz_stroke_right...
    Fx_pitch_right Fy_pitch_right Fz_pitch_right Mx_pitch_right My_pitch_right Mz_pitch_right...
    Fx_dev_right Fy_dev_right Fz_dev_right Mx_dev_right My_dev_right Mz_dev_right...
    stroke_L_all_right stroke_R_all_right pitch_L_all_right pitch_R_all_right dev_L_all_right dev_R_all_right...
    stroke_L_allNOfreq_right stroke_R_allNOfreq_right pitch_L_allNOfreq_right pitch_R_allNOfreq_right dev_L_allNOfreq_right dev_R_allNOfreq_right...
    stroke_L_freq_right stroke_R_freq_right pitch_L_freq_right pitch_R_freq_right dev_L_freq_right dev_R_freq_right...
    stroke_L_stroke_right stroke_R_stroke_right pitch_L_stroke_right pitch_R_stroke_right dev_L_stroke_right dev_R_stroke_right...
    stroke_L_pitch_right stroke_R_pitch_right pitch_L_pitch_right pitch_R_pitch_right dev_L_pitch_right dev_R_pitch_right...
    stroke_L_dev_right stroke_R_dev_right pitch_L_dev_right pitch_R_dev_right dev_L_dev_right dev_R_dev_right...
    t_all_right t_allNOfreq_right t_stroke_right t_pitch_right t_dev_right t_freq_right...
    Fx_all_left Fy_all_left Fz_all_left Mx_all_left My_all_left Mz_all_left...
    Fx_allNOfreq_left Fy_allNOfreq_left Fz_allNOfreq_left Mx_allNOfreq_left My_allNOfreq_left Mz_allNOfreq_left...
    Fx_freq_left Fy_freq_left Fz_freq_left Mx_freq_left My_freq_left Mz_freq_left...
    Fx_stroke_left Fy_stroke_left Fz_stroke_left Mx_stroke_left My_stroke_left Mz_stroke_left...
    Fx_pitch_left Fy_pitch_left Fz_pitch_left Mx_pitch_left My_pitch_left Mz_pitch_left...
    Fx_dev_left Fy_dev_left Fz_dev_left Mx_dev_left My_dev_left Mz_dev_left...
    stroke_L_all_left stroke_R_all_left pitch_L_all_left pitch_R_all_left dev_L_all_left dev_R_all_left...
    stroke_L_allNOfreq_left stroke_R_allNOfreq_left pitch_L_allNOfreq_left pitch_R_allNOfreq_left dev_L_allNOfreq_left dev_R_allNOfreq_left...
    stroke_L_freq_left stroke_R_freq_left pitch_L_freq_left pitch_R_freq_left dev_L_freq_left dev_R_freq_left...
    stroke_L_stroke_left stroke_R_stroke_left pitch_L_stroke_left pitch_R_stroke_left dev_L_stroke_left dev_R_stroke_left...
    stroke_L_pitch_left stroke_R_pitch_left pitch_L_pitch_left pitch_R_pitch_left dev_L_pitch_left dev_R_pitch_left...
    stroke_L_dev_left stroke_R_dev_left pitch_L_dev_left pitch_R_dev_left dev_L_dev_left dev_R_dev_left...
    t_all_left t_allNOfreq_left t_stroke_left t_pitch_left t_dev_left t_freq_left

save(['borf_db_',save_name,'_NOcali_means_LnRsym.mat'])
% save(['borf_db_',save_name,'_NOcali_means_sym_',side,'.mat'])
