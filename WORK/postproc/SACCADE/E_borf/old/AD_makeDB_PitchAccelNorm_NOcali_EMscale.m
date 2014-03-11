% make borf dataset Fenhance

clear
clc
close all
warning off

save_name = 'PitchAccel';

%% cali data
% addpath('/home/florian/DATA/flighttracker/looming_FlyCity_borf')
% cali_file = 'cali_matrix_08132013.mat'
% load(cali_file)

%% settings
butter_cut = 10;
butter_n = 5;

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

%% make db PitchAccel

% file_name = [save_name,'_all_'];
% makeDB_allMODs_NOcali_WBs

file_name = [save_name,'_all_const_freq_'];
makeDB_allNOfreqMODs_NOcali_WBs

% file_name = [save_name,'_freq_'];
% makeDB_freqMODs_NOcali_WBs

file_name = [save_name,'_stroke_'];
makeDB_strokeMODs_NOcali_WBs

file_name = [save_name,'_rot_'];
makeDB_pitchMODs_NOcali_WBs

file_name = [save_name,'_dev_'];
makeDB_devMODs_NOcali_WBs


%% postprocess
% MODmin = min([length(My_mean_freq) length(My_mean_stroke) length(My_mean_pitch) length(My_mean_dev)]);
% My_steady = My_mean_all(find(mod_value_all == 0));
% My_mean_sum = My_steady + (My_mean_freq(1:MODmin)-My_steady) + (My_mean_stroke(1:MODmin)-My_steady) + (My_mean_pitch(1:MODmin)-My_steady) + (My_mean_dev(1:MODmin)-My_steady);
% My_mean_sum_nopitch = My_steady + (My_mean_freq(1:MODmin)-My_steady) + (My_mean_stroke(1:MODmin)-My_steady) + (My_mean_dev(1:MODmin)-My_steady);

MODmin = min([length(My_mean_stroke) length(My_mean_pitch) length(My_mean_dev)]);
My_steady = My_mean_allNOfreq(find(mod_value_allNOfreq == 0));
My_mean_sum = My_steady + (My_mean_stroke(1:MODmin)-My_steady) + (My_mean_pitch(1:MODmin)-My_steady) + (My_mean_dev(1:MODmin)-My_steady);
My_mean_sum_nopitch = My_steady + (My_mean_stroke(1:MODmin)-My_steady) + (My_mean_dev(1:MODmin)-My_steady);

%% plot data
cd ..
mkdir('figs_NOcali')
cd('figs_NOcali')
plot_PitchAccelNorm_MODdata

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

save(['borf_db_',save_name,'_NOcali_alldata.mat'])

clear Fx_all Fy_all Fz_all Mx_all My_all Mz_all...
    Fx_freq Fy_freq Fz_freq Mx_freq My_freq Mz_freq...
    Fx_stroke Fy_stroke Fz_stroke Mx_stroke My_stroke Mz_stroke...
    Fx_pitch Fy_pitch Fz_pitch Mx_pitch My_pitch Mz_pitch...
    Fx_dev Fy_dev Fz_dev Mx_dev My_dev Mz_dev...
    stroke_L_all stroke_R_all pitch_L_all pitch_R_all dev_L_all dev_R_all...
    stroke_L_freq stroke_R_freq pitch_L_freq pitch_R_freq dev_L_freq dev_R_freq...
    stroke_L_stroke stroke_R_stroke pitch_L_stroke pitch_R_stroke dev_L_stroke dev_R_stroke...
    stroke_L_pitch stroke_R_pitch pitch_L_pitch pitch_R_pitch dev_L_pitch dev_R_pitch...
    stroke_L_dev stroke_R_dev pitch_L_dev pitch_R_dev dev_L_dev dev_R_dev...
    t_all t_stroke t_pitch t_dev t_freq
save(['borf_db_',save_name,'_NOcali_means.mat'])
