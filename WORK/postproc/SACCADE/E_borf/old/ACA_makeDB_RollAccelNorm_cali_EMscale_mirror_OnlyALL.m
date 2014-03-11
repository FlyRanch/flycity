% make borf dataset Fenhance

clear
clc
close all
warning off

% save_name = 'F_roll_DLR';
save_name = 'F_roll_LR';
mirror_nr = 1;

%% cali data
addpath('/home/florian/DATA/flighttracker/looming_FlyCity_borf')
% cali_file = 'cali_matrix_interp_08132013.mat'
% cali_file = 'interp_cali_matrix_09292013.mat'
cali_file = 'interp_cali_matrix_10302013.mat'
load(cali_file)

%% settings
% butter_cut = 10;
% butter_n = 5;

seq_length_max = 15000;
Nwb_total = 20;
Nwb = 16;
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

%% mirrored FIRST
mirror_nr = 1;

% file_name = [save_name,'_all_mirror',num2str(mirror_nr),'_'];
% makeDB_allMODs_cali_interp

file_name = [save_name,'_all_const_freq_mirror',num2str(mirror_nr),'_'];
makeDB_allNOfreqMODs_cali_interp

% file_name = [save_name,'_freq_mirror',num2str(mirror_nr),'_'];
% makeDB_freqMODs_cali_interp

% file_name = [save_name,'_stroke_mirror',num2str(mirror_nr),'_'];
% makeDB_strokeMODs_cali_interp
% 
% file_name = [save_name,'_rot_mirror',num2str(mirror_nr),'_'];
% makeDB_pitchMODs_cali_interp
% 
% file_name = [save_name,'_dev_mirror',num2str(mirror_nr),'_'];
% makeDB_devMODs_cali_interp


%% store data

% allNOfreq
Fx_mean_allNOfreq_mirror1 = Fx_mean_allNOfreq;
Fy_mean_allNOfreq_mirror1 = Fy_mean_allNOfreq;
Fz_mean_allNOfreq_mirror1 = Fz_mean_allNOfreq;
Mx_mean_allNOfreq_mirror1 = Mx_mean_allNOfreq;
My_mean_allNOfreq_mirror1 = My_mean_allNOfreq;
Mz_mean_allNOfreq_mirror1 = Mz_mean_allNOfreq;

Fx_allNOfreq_mirror1 = Fx_allNOfreq;
Fy_allNOfreq_mirror1 = Fy_allNOfreq;
Fz_allNOfreq_mirror1 = Fz_allNOfreq;
Mx_allNOfreq_mirror1 = Mx_allNOfreq;
My_allNOfreq_mirror1 = My_allNOfreq;
Mz_allNOfreq_mirror1 = Mz_allNOfreq;

t_allNOfreq_mirror1 = t_allNOfreq;
mod_value_allNOfreq_mirror1 = mod_value_allNOfreq;

stroke_L_allNOfreq_mirror1 = stroke_L_allNOfreq;
pitch_L_allNOfreq_mirror1 = pitch_L_allNOfreq;
dev_L_allNOfreq_mirror1 = dev_L_allNOfreq;

stroke_R_allNOfreq_mirror1 = stroke_R_allNOfreq;
pitch_R_allNOfreq_mirror1 = pitch_R_allNOfreq;
dev_R_allNOfreq_mirror1 = dev_R_allNOfreq;

% % stroke
% Fx_mean_stroke_mirror1 = Fx_mean_stroke;
% Fy_mean_stroke_mirror1 = Fy_mean_stroke;
% Fz_mean_stroke_mirror1 = Fz_mean_stroke;
% Mx_mean_stroke_mirror1 = Mx_mean_stroke;
% My_mean_stroke_mirror1 = My_mean_stroke;
% Mz_mean_stroke_mirror1 = Mz_mean_stroke;
% 
% Fx_stroke_mirror1 = Fx_stroke;
% Fy_stroke_mirror1 = Fy_stroke;
% Fz_stroke_mirror1 = Fz_stroke;
% Mx_stroke_mirror1 = Mx_stroke;
% My_stroke_mirror1 = My_stroke;
% Mz_stroke_mirror1 = Mz_stroke;
% 
% t_stroke_mirror1 = t_stroke;
% mod_value_stroke_mirror1 = mod_value_stroke;
% 
% stroke_L_stroke_mirror1 = stroke_L_stroke;
% pitch_L_stroke_mirror1 = pitch_L_stroke;
% dev_L_stroke_mirror1 = dev_L_stroke;
% 
% stroke_R_stroke_mirror1 = stroke_R_stroke;
% pitch_R_stroke_mirror1 = pitch_R_stroke;
% dev_R_stroke_mirror1 = dev_R_stroke;
% 
% % pitch
% Fx_mean_pitch_mirror1 = Fx_mean_pitch;
% Fy_mean_pitch_mirror1 = Fy_mean_pitch;
% Fz_mean_pitch_mirror1 = Fz_mean_pitch;
% Mx_mean_pitch_mirror1 = Mx_mean_pitch;
% My_mean_pitch_mirror1 = My_mean_pitch;
% Mz_mean_pitch_mirror1 = Mz_mean_pitch;
% 
% Fx_pitch_mirror1 = Fx_pitch;
% Fy_pitch_mirror1 = Fy_pitch;
% Fz_pitch_mirror1 = Fz_pitch;
% Mx_pitch_mirror1 = Mx_pitch;
% My_pitch_mirror1 = My_pitch;
% Mz_pitch_mirror1 = Mz_pitch;
% 
% t_pitch_mirror1 = t_pitch;
% mod_value_pitch_mirror1 = mod_value_pitch;
% 
% stroke_L_pitch_mirror1 = stroke_L_pitch;
% pitch_L_pitch_mirror1 = pitch_L_pitch;
% dev_L_pitch_mirror1 = dev_L_pitch;
% 
% stroke_R_pitch_mirror1 = stroke_R_pitch;
% pitch_R_pitch_mirror1 = pitch_R_pitch;
% dev_R_pitch_mirror1 = dev_R_pitch;
% 
% % dev
% Fx_mean_dev_mirror1 = Fx_mean_dev;
% Fy_mean_dev_mirror1 = Fy_mean_dev;
% Fz_mean_dev_mirror1 = Fz_mean_dev;
% Mx_mean_dev_mirror1 = Mx_mean_dev;
% My_mean_dev_mirror1 = My_mean_dev;
% Mz_mean_dev_mirror1 = Mz_mean_dev;
% 
% Fx_dev_mirror1 = Fx_dev;
% Fy_dev_mirror1 = Fy_dev;
% Fz_dev_mirror1 = Fz_dev;
% Mx_dev_mirror1 = Mx_dev;
% My_dev_mirror1 = My_dev;
% Mz_dev_mirror1 = Mz_dev;
% 
% t_dev_mirror1 = t_dev;
% mod_value_dev_mirror1 = mod_value_dev;
% 
% stroke_L_dev_mirror1 = stroke_L_dev;
% pitch_L_dev_mirror1 = pitch_L_dev;
% dev_L_dev_mirror1 = dev_L_dev;
% 
% stroke_R_dev_mirror1 = stroke_R_dev;
% pitch_R_dev_mirror1 = pitch_R_dev;
% dev_R_dev_mirror1 = dev_R_dev;

%% NO mirror
mirror_nr = 0;

% file_name = [save_name,'_all_mirror',num2str(mirror_nr),'_'];
% makeDB_allMODs_cali_interp

file_name = [save_name,'_all_const_freq_mirror',num2str(mirror_nr),'_'];
makeDB_allNOfreqMODs_cali_interp

% file_name = [save_name,'_freq_mirror',num2str(mirror_nr),'_'];
% makeDB_freqMODs_cali_interp

% file_name = [save_name,'_stroke_mirror',num2str(mirror_nr),'_'];
% makeDB_strokeMODs_cali_interp
% 
% file_name = [save_name,'_rot_mirror',num2str(mirror_nr),'_'];
% makeDB_pitchMODs_cali_interp
% 
% file_name = [save_name,'_dev_mirror',num2str(mirror_nr),'_'];
% makeDB_devMODs_cali_interp

%% store data

% allNOfreq
Fx_mean_allNOfreq_mirror0 = Fx_mean_allNOfreq;
Fy_mean_allNOfreq_mirror0 = Fy_mean_allNOfreq;
Fz_mean_allNOfreq_mirror0 = Fz_mean_allNOfreq;
Mx_mean_allNOfreq_mirror0 = Mx_mean_allNOfreq;
My_mean_allNOfreq_mirror0 = My_mean_allNOfreq;
Mz_mean_allNOfreq_mirror0 = Mz_mean_allNOfreq;

Fx_allNOfreq_mirror0 = Fx_allNOfreq;
Fy_allNOfreq_mirror0 = Fy_allNOfreq;
Fz_allNOfreq_mirror0 = Fz_allNOfreq;
Mx_allNOfreq_mirror0 = Mx_allNOfreq;
My_allNOfreq_mirror0 = My_allNOfreq;
Mz_allNOfreq_mirror0 = Mz_allNOfreq;

t_allNOfreq_mirror0 = t_allNOfreq;
mod_value_allNOfreq_mirror0 = mod_value_allNOfreq;

stroke_L_allNOfreq_mirror0 = stroke_L_allNOfreq;
pitch_L_allNOfreq_mirror0 = pitch_L_allNOfreq;
dev_L_allNOfreq_mirror0 = dev_L_allNOfreq;

stroke_R_allNOfreq_mirror0 = stroke_R_allNOfreq;
pitch_R_allNOfreq_mirror0 = pitch_R_allNOfreq;
dev_R_allNOfreq_mirror0 = dev_R_allNOfreq;

% % stroke
% Fx_mean_stroke_mirror0 = Fx_mean_stroke;
% Fy_mean_stroke_mirror0 = Fy_mean_stroke;
% Fz_mean_stroke_mirror0 = Fz_mean_stroke;
% Mx_mean_stroke_mirror0 = Mx_mean_stroke;
% My_mean_stroke_mirror0 = My_mean_stroke;
% Mz_mean_stroke_mirror0 = Mz_mean_stroke;
% 
% Fx_stroke_mirror0 = Fx_stroke;
% Fy_stroke_mirror0 = Fy_stroke;
% Fz_stroke_mirror0 = Fz_stroke;
% Mx_stroke_mirror0 = Mx_stroke;
% My_stroke_mirror0 = My_stroke;
% Mz_stroke_mirror0 = Mz_stroke;
% 
% t_stroke_mirror0 = t_stroke;
% mod_value_stroke_mirror0 = mod_value_stroke;
% 
% stroke_L_stroke_mirror0 = stroke_L_stroke;
% pitch_L_stroke_mirror0 = pitch_L_stroke;
% dev_L_stroke_mirror0 = dev_L_stroke;
% 
% stroke_R_stroke_mirror0 = stroke_R_stroke;
% pitch_R_stroke_mirror0 = pitch_R_stroke;
% dev_R_stroke_mirror0 = dev_R_stroke;
% 
% % pitch
% Fx_mean_pitch_mirror0 = Fx_mean_pitch;
% Fy_mean_pitch_mirror0 = Fy_mean_pitch;
% Fz_mean_pitch_mirror0 = Fz_mean_pitch;
% Mx_mean_pitch_mirror0 = Mx_mean_pitch;
% My_mean_pitch_mirror0 = My_mean_pitch;
% Mz_mean_pitch_mirror0 = Mz_mean_pitch;
% 
% Fx_pitch_mirror0 = Fx_pitch;
% Fy_pitch_mirror0 = Fy_pitch;
% Fz_pitch_mirror0 = Fz_pitch;
% Mx_pitch_mirror0 = Mx_pitch;
% My_pitch_mirror0 = My_pitch;
% Mz_pitch_mirror0 = Mz_pitch;
% 
% t_pitch_mirror0 = t_pitch;
% mod_value_pitch_mirror0 = mod_value_pitch;
% 
% stroke_L_pitch_mirror0 = stroke_L_pitch;
% pitch_L_pitch_mirror0 = pitch_L_pitch;
% dev_L_pitch_mirror0 = dev_L_pitch;
% 
% stroke_R_pitch_mirror0 = stroke_R_pitch;
% pitch_R_pitch_mirror0 = pitch_R_pitch;
% dev_R_pitch_mirror0 = dev_R_pitch;
% 
% % dev
% Fx_mean_dev_mirror0 = Fx_mean_dev;
% Fy_mean_dev_mirror0 = Fy_mean_dev;
% Fz_mean_dev_mirror0 = Fz_mean_dev;
% Mx_mean_dev_mirror0 = Mx_mean_dev;
% My_mean_dev_mirror0 = My_mean_dev;
% Mz_mean_dev_mirror0 = Mz_mean_dev;
% 
% Fx_dev_mirror0 = Fx_dev;
% Fy_dev_mirror0 = Fy_dev;
% Fz_dev_mirror0 = Fz_dev;
% Mx_dev_mirror0 = Mx_dev;
% My_dev_mirror0 = My_dev;
% Mz_dev_mirror0 = Mz_dev;
% 
% t_dev_mirror0 = t_dev;
% mod_value_dev_mirror0 = mod_value_dev;
% 
% stroke_L_dev_mirror0 = stroke_L_dev;
% pitch_L_dev_mirror0 = pitch_L_dev;
% dev_L_dev_mirror0 = dev_L_dev;
% 
% stroke_R_dev_mirror0 = stroke_R_dev;
% pitch_R_dev_mirror0 = pitch_R_dev;
% dev_R_dev_mirror0 = dev_R_dev;

%% merge data

% allNOfreq
diff = abs(stroke_L_allNOfreq_mirror0 - stroke_R_allNOfreq_mirror1);
if nanmax(diff(:)) == 0
    
    Fx_mean_allNOfreq = (Fx_mean_allNOfreq_mirror0 + Fx_mean_allNOfreq_mirror1)/2;
    Fy_mean_allNOfreq = (Fy_mean_allNOfreq_mirror0 - Fy_mean_allNOfreq_mirror1)/2;
    Fz_mean_allNOfreq = (Fz_mean_allNOfreq_mirror0 + Fz_mean_allNOfreq_mirror1)/2;

    Mx_mean_allNOfreq = (Mx_mean_allNOfreq_mirror0 - Mx_mean_allNOfreq_mirror1)/2;
    My_mean_allNOfreq = (My_mean_allNOfreq_mirror0 + My_mean_allNOfreq_mirror1)/2;
    Mz_mean_allNOfreq = (Mz_mean_allNOfreq_mirror0 - Mz_mean_allNOfreq_mirror1)/2;
    
    Fx_allNOfreq = (Fx_allNOfreq_mirror0 + Fx_allNOfreq_mirror1)/2;
    Fy_allNOfreq = (Fy_allNOfreq_mirror0 - Fy_allNOfreq_mirror1)/2;
    Fz_allNOfreq = (Fz_allNOfreq_mirror0 + Fz_allNOfreq_mirror1)/2;

    Mx_allNOfreq = (Mx_allNOfreq_mirror0 - Mx_allNOfreq_mirror1)/2;
    My_allNOfreq = (My_allNOfreq_mirror0 + My_allNOfreq_mirror1)/2;
    Mz_allNOfreq = (Mz_allNOfreq_mirror0 - Mz_allNOfreq_mirror1)/2;
else
    error = 'allNOfreq'
    pause
end

% % stroke
% diff = abs(stroke_L_stroke_mirror0 - stroke_R_stroke_mirror1);
% if nanmax(diff(:)) == 0
%     
%     Fx_mean_stroke = (Fx_mean_stroke_mirror0 + Fx_mean_stroke_mirror1)/2;
%     Fy_mean_stroke = (Fy_mean_stroke_mirror0 - Fy_mean_stroke_mirror1)/2;
%     Fz_mean_stroke = (Fz_mean_stroke_mirror0 + Fz_mean_stroke_mirror1)/2;
% 
%     Mx_mean_stroke = (Mx_mean_stroke_mirror0 - Mx_mean_stroke_mirror1)/2;
%     My_mean_stroke = (My_mean_stroke_mirror0 + My_mean_stroke_mirror1)/2;
%     Mz_mean_stroke = (Mz_mean_stroke_mirror0 - Mz_mean_stroke_mirror1)/2;
%     
%     Fx_stroke = (Fx_stroke_mirror0 + Fx_stroke_mirror1)/2;
%     Fy_stroke = (Fy_stroke_mirror0 - Fy_stroke_mirror1)/2;
%     Fz_stroke = (Fz_stroke_mirror0 + Fz_stroke_mirror1)/2;
% 
%     Mx_stroke = (Mx_stroke_mirror0 - Mx_stroke_mirror1)/2;
%     My_stroke = (My_stroke_mirror0 + My_stroke_mirror1)/2;
%     Mz_stroke = (Mz_stroke_mirror0 - Mz_stroke_mirror1)/2;
% else
%     error = 'stroke'
%     pause
% end
% 
% % pitch
% diff = abs(pitch_L_pitch_mirror0 - pitch_R_pitch_mirror1);
% if nanmax(diff(:)) == 0
%     
%     Fx_mean_pitch = (Fx_mean_pitch_mirror0 + Fx_mean_pitch_mirror1)/2;
%     Fy_mean_pitch = (Fy_mean_pitch_mirror0 - Fy_mean_pitch_mirror1)/2;
%     Fz_mean_pitch = (Fz_mean_pitch_mirror0 + Fz_mean_pitch_mirror1)/2;
% 
%     Mx_mean_pitch = (Mx_mean_pitch_mirror0 - Mx_mean_pitch_mirror1)/2;
%     My_mean_pitch = (My_mean_pitch_mirror0 + My_mean_pitch_mirror1)/2;
%     Mz_mean_pitch = (Mz_mean_pitch_mirror0 - Mz_mean_pitch_mirror1)/2;
%     
%     Fx_pitch = (Fx_pitch_mirror0 + Fx_pitch_mirror1)/2;
%     Fy_pitch = (Fy_pitch_mirror0 - Fy_pitch_mirror1)/2;
%     Fz_pitch = (Fz_pitch_mirror0 + Fz_pitch_mirror1)/2;
% 
%     Mx_pitch = (Mx_pitch_mirror0 - Mx_pitch_mirror1)/2;
%     My_pitch = (My_pitch_mirror0 + My_pitch_mirror1)/2;
%     Mz_pitch = (Mz_pitch_mirror0 - Mz_pitch_mirror1)/2;
% else
%     error = 'pitch'
%     pause
% end
% 
% % dev
% diff = abs(dev_L_dev_mirror0 - dev_R_dev_mirror1);
% if nanmax(diff(:)) == 0
%     
%     Fx_mean_dev = (Fx_mean_dev_mirror0 + Fx_mean_dev_mirror1)/2;
%     Fy_mean_dev = (Fy_mean_dev_mirror0 - Fy_mean_dev_mirror1)/2;
%     Fz_mean_dev = (Fz_mean_dev_mirror0 + Fz_mean_dev_mirror1)/2;
% 
%     Mx_mean_dev = (Mx_mean_dev_mirror0 - Mx_mean_dev_mirror1)/2;
%     My_mean_dev = (My_mean_dev_mirror0 + My_mean_dev_mirror1)/2;
%     Mz_mean_dev = (Mz_mean_dev_mirror0 - Mz_mean_dev_mirror1)/2;
%     
%     Fx_dev = (Fx_dev_mirror0 + Fx_dev_mirror1)/2;
%     Fy_dev = (Fy_dev_mirror0 - Fy_dev_mirror1)/2;
%     Fz_dev = (Fz_dev_mirror0 + Fz_dev_mirror1)/2;
% 
%     Mx_dev = (Mx_dev_mirror0 - Mx_dev_mirror1)/2;
%     My_dev = (My_dev_mirror0 + My_dev_mirror1)/2;
%     Mz_dev = (Mz_dev_mirror0 - Mz_dev_mirror1)/2;
% else
%     error = 'dev'
%     pause
% end

%% postprocess
% MODmin = min([length(Mx_mean_allNOfreq) length(Mx_mean_stroke) length(Mx_mean_pitch) length(Mx_mean_dev)]);
% Mx_steady = Mx_mean_allNOfreq(mod_value_allNOfreq == 0);
% Mx_mean_sum = Mx_steady + (Mx_mean_stroke(1:MODmin)-Mx_steady) + (Mx_mean_pitch(1:MODmin)-Mx_steady) + (Mx_mean_dev(1:MODmin)-Mx_steady);
% Mx_mean_sum_nopitch = Mx_steady + (Mx_mean_stroke(1:MODmin)-Mx_steady) + (Mx_mean_dev(1:MODmin)-Mx_steady);

%% plot data
cd ..
mkdir('figs_cali')
cd('figs_cali')
% plot_RollAccelNorm_LR_MODdata_mirror
% plot_RollAccelNorm_LR_MODdata
figure
hold on
plot(mod_value_allNOfreq,Mx_mean_allNOfreq,'.b')
plot(mod_value_allNOfreq,My_mean_allNOfreq,'.g')
plot(mod_value_allNOfreq,Mz_mean_allNOfreq,'.r')
legend('Mroll','Mpitch','Myaw','location','nw')
grid on
xlabel('MOD value')
ylabel('M')
saveas(gca, [save_name,'_MODvalue_allNOfreqMOD.png'])
cd ..

% save data
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

save(['borf_db_',save_name,'_cali_alldata.mat'])
% save(['borf_db_',save_name,'_cali_alldata_mirror',num2str(mirror_nr),'.mat'])

clear cali_L cali_R str_interp dev_interp rot_interp...
    Fx_all Fy_all Fz_all Mx_all My_all Mz_all...
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
    Fx_all_mirror0 Fy_all_mirror0 Fz_all_mirror0 Mx_all_mirror0 My_all_mirror0 Mz_all_mirror0...
    Fx_allNOfreq_mirror0 Fy_allNOfreq_mirror0 Fz_allNOfreq_mirror0 Mx_allNOfreq_mirror0 My_allNOfreq_mirror0 Mz_allNOfreq_mirror0...
    Fx_freq_mirror0 Fy_freq_mirror0 Fz_freq_mirror0 Mx_freq_mirror0 My_freq_mirror0 Mz_freq_mirror0...
    Fx_stroke_mirror0 Fy_stroke_mirror0 Fz_stroke_mirror0 Mx_stroke_mirror0 My_stroke_mirror0 Mz_stroke_mirror0...
    Fx_pitch_mirror0 Fy_pitch_mirror0 Fz_pitch_mirror0 Mx_pitch_mirror0 My_pitch_mirror0 Mz_pitch_mirror0...
    Fx_dev_mirror0 Fy_dev_mirror0 Fz_dev_mirror0 Mx_dev_mirror0 My_dev_mirror0 Mz_dev_mirror0...
    stroke_L_all_mirror0 stroke_R_all_mirror0 pitch_L_all_mirror0 pitch_R_all_mirror0 dev_L_all_mirror0 dev_R_all_mirror0...
    stroke_L_allNOfreq_mirror0 stroke_R_allNOfreq_mirror0 pitch_L_allNOfreq_mirror0 pitch_R_allNOfreq_mirror0 dev_L_allNOfreq_mirror0 dev_R_allNOfreq_mirror0...
    stroke_L_freq_mirror0 stroke_R_freq_mirror0 pitch_L_freq_mirror0 pitch_R_freq_mirror0 dev_L_freq_mirror0 dev_R_freq_mirror0...
    stroke_L_stroke_mirror0 stroke_R_stroke_mirror0 pitch_L_stroke_mirror0 pitch_R_stroke_mirror0 dev_L_stroke_mirror0 dev_R_stroke_mirror0...
    stroke_L_pitch_mirror0 stroke_R_pitch_mirror0 pitch_L_pitch_mirror0 pitch_R_pitch_mirror0 dev_L_pitch_mirror0 dev_R_pitch_mirror0...
    stroke_L_dev_mirror0 stroke_R_dev_mirror0 pitch_L_dev_mirror0 pitch_R_dev_mirror0 dev_L_dev_mirror0 dev_R_dev_mirror0...
    t_all_mirror0 t_allNOfreq_mirror0 t_stroke_mirror0 t_pitch_mirror0 t_dev_mirror0 t_freq_mirror0...
    Fx_all_mirror1 Fy_all_mirror1 Fz_all_mirror1 Mx_all_mirror1 My_all_mirror1 Mz_all_mirror1...
    Fx_allNOfreq_mirror1 Fy_allNOfreq_mirror1 Fz_allNOfreq_mirror1 Mx_allNOfreq_mirror1 My_allNOfreq_mirror1 Mz_allNOfreq_mirror1...
    Fx_freq_mirror1 Fy_freq_mirror1 Fz_freq_mirror1 Mx_freq_mirror1 My_freq_mirror1 Mz_freq_mirror1...
    Fx_stroke_mirror1 Fy_stroke_mirror1 Fz_stroke_mirror1 Mx_stroke_mirror1 My_stroke_mirror1 Mz_stroke_mirror1...
    Fx_pitch_mirror1 Fy_pitch_mirror1 Fz_pitch_mirror1 Mx_pitch_mirror1 My_pitch_mirror1 Mz_pitch_mirror1...
    Fx_dev_mirror1 Fy_dev_mirror1 Fz_dev_mirror1 Mx_dev_mirror1 My_dev_mirror1 Mz_dev_mirror1...
    stroke_L_all_mirror1 stroke_R_all_mirror1 pitch_L_all_mirror1 pitch_R_all_mirror1 dev_L_all_mirror1 dev_R_all_mirror1...
    stroke_L_allNOfreq_mirror1 stroke_R_allNOfreq_mirror1 pitch_L_allNOfreq_mirror1 pitch_R_allNOfreq_mirror1 dev_L_allNOfreq_mirror1 dev_R_allNOfreq_mirror1...
    stroke_L_freq_mirror1 stroke_R_freq_mirror1 pitch_L_freq_mirror1 pitch_R_freq_mirror1 dev_L_freq_mirror1 dev_R_freq_mirror1...
    stroke_L_stroke_mirror1 stroke_R_stroke_mirror1 pitch_L_stroke_mirror1 pitch_R_stroke_mirror1 dev_L_stroke_mirror1 dev_R_stroke_mirror1...
    stroke_L_pitch_mirror1 stroke_R_pitch_mirror1 pitch_L_pitch_mirror1 pitch_R_pitch_mirror1 dev_L_pitch_mirror1 dev_R_pitch_mirror1...
    stroke_L_dev_mirror1 stroke_R_dev_mirror1 pitch_L_dev_mirror1 pitch_R_dev_mirror1 dev_L_dev_mirror1 dev_R_dev_mirror1...
    t_all_mirror1 t_allNOfreq_mirror1 t_stroke_mirror1 t_pitch_mirror1 t_dev_mirror1 t_freq_mirror1

save(['borf_db_',save_name,'_cali_means.mat'])
% save(['borf_db_',save_name,'_cali_means_mirror',num2str(mirror_nr),'.mat'])
