% make borf dataset F_roll_LR

clear
clc
close all

save_name = 'F_roll_DLR';
% save_name = 'F_roll_LR';

% settings
seq_length_max = 140000;
Nwb = 7;
wb_start = 3;
wb_stop = 6;

% robo data
Lrobo = .23;
crobo = .065;
ARrobo = Lrobo/crobo;

% fly data
Lwing = .00299; %m
cwing = Lwing/ARrobo;

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

% robo2fly
f_robo2fly = (v_air/v_oil) * (Lrobo/Lwing)^2;
F_robo2fly = (rho_air/rho_oil) * (v_air/v_oil)^2;
M_robo2fly = F_robo2fly * (Lwing/Lrobo);

%% Elzinga scaling
ARwing_fly  = 3.1579; % wing aspect ratio of fly (based on Elzinga's R/c)
c_fly  = Lwing/ARwing_fly; %mean chord of actual fly in m

l_scale_up = crobo/c_fly;
f_robo2fly = (v_air*crobo*Lrobo)/(c_fly*Lwing*v_oil);
F_robo2fly = (rho_air/rho_oil)*((c_fly/crobo)^4)*((f_robo2fly)^2);
M_robo2fly = (rho_air/rho_oil)*((c_fly/crobo)^5)*((f_robo2fly)^2);

%% NORM data
Nm2mNmm = 1e6;

%% colormap
imax = 21;

% black2blue
colormap_black2blue(1:imax,1) = 0;
colormap_black2blue(1:imax,2) = [0:.5/(imax-1):.5];
colormap_black2blue(1:imax,3) = [0:1/(imax-1):1];


%% make db F_roll_LR
makeDB_allMODs

makeDB_allNOfreqMODs

makeDB_freqMODs

makeDB_strokeMODs

makeDB_pitchMODs

makeDB_devMODs

MODmin = min([length(Mx_mean_freq) length(Mx_mean_stroke) length(Mx_mean_pitch) length(Mx_mean_dev)]);
Mx_steady = Mx_mean_all(find(mod_value_all == 0));
Mx_mean_sum = Mx_steady + (Mx_mean_freq(1:MODmin)-Mx_steady) + (Mx_mean_stroke(1:MODmin)-Mx_steady) + (Mx_mean_pitch(1:MODmin)-Mx_steady) + (Mx_mean_dev(1:MODmin)-Mx_steady);
Mx_mean_sum_nopitch = Mx_steady + (Mx_mean_freq(1:MODmin)-Mx_steady) + (Mx_mean_stroke(1:MODmin)-Mx_steady) + (Mx_mean_dev(1:MODmin)-Mx_steady);

%% plot data
cd ..
mkdir('figs')
cd('figs')
plot_RollAccel_LR_MODdata

% save data
cd ..
clear exp_type ft kine mod_val t list i t_now Fx_now Fy_now Fz_now Mx_now My_now Mz_now stroke_L_now pitch_L_now dev_L_now stroke_R_now pitch_R_now dev_R_now N_wb_now n_start n_stop frobo_now f_now
save(['borf_db_',save_name,'_alldata.mat'])

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
save(['borf_db_',save_name,'_means.mat'])
