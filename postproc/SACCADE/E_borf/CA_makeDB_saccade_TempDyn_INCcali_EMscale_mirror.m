% make borf dataset Fenhance

clear
clc
close all
warning off


%% fly variables
var_file = dir('flyVar*')
load(var_file.name)

%% cali data
% cali_file = dir('cali_matrix*');
cali_file = dir('interp_cali_matrix*');
cali_file = cali_file.name
load(cali_file)

%% borf data
load_dir = dir('Saccade_Temporal_0*');
load_dir = load_dir.name
load_name = 'F_sacc_ts';
save_name =  load_name;

cd(load_dir)

%% settings
butter_cut = 10;
butter_n = 5;

Nwb = 27;
WB_start = 6;
WB_stop = 1;
WB0 = 4;

% robo data
Lrobo = .23;
crobo = .065;
ARrobo = Lrobo/crobo;

% fluid data
v_air = 1.568e-5; %kinematic viscocity in m^2/s
v_oil = 115.0*1e-6; %kinematic viscocity in m^2/s
% rho_air = 1.22521; %kg/m^3
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
% imax = 21;
% 
% % black2blue
% colormap_black2blue(1:imax,1) = 0;
% colormap_black2blue(1:imax,2) = [0:.5/(imax-1):.5];
% colormap_black2blue(1:imax,3) = [0:1/(imax-1):1];

%% make db Fenhance

%% mirrored FIRST
mirror_nr = 1;

file_name = [save_name,'_mirror',num2str(mirror_nr),'_'];

calc_WBstart % determine WB start only once

makeDB_TempDyn_cali_interp

%% store data & remove pre&post WBs
Fx_mean_mean_mirror1 = Fx_mean_mean(WB_start:Nwb-WB_stop);
Fy_mean_mean_mirror1 = Fy_mean_mean(WB_start:Nwb-WB_stop);
Fz_mean_mean_mirror1 = Fz_mean_mean(WB_start:Nwb-WB_stop);
Mx_mean_mean_mirror1 = Mx_mean_mean(WB_start:Nwb-WB_stop);
My_mean_mean_mirror1 = My_mean_mean(WB_start:Nwb-WB_stop);
Mz_mean_mean_mirror1 = Mz_mean_mean(WB_start:Nwb-WB_stop);

Fx_mean_all_mirror1 = Fx_mean_all(wb_start(WB_start):wb_start(Nwb-WB_stop));
Fy_mean_all_mirror1 = Fy_mean_all(wb_start(WB_start):wb_start(Nwb-WB_stop));
Fz_mean_all_mirror1 = Fz_mean_all(wb_start(WB_start):wb_start(Nwb-WB_stop));
Mx_mean_all_mirror1 = Mx_mean_all(wb_start(WB_start):wb_start(Nwb-WB_stop));
My_mean_all_mirror1 = My_mean_all(wb_start(WB_start):wb_start(Nwb-WB_stop));
Mz_mean_all_mirror1 = Mz_mean_all(wb_start(WB_start):wb_start(Nwb-WB_stop));

Fx_all_mirror1 = Fx_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);
Fy_all_mirror1 = Fy_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);
Fz_all_mirror1 = Fz_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);
Mx_all_mirror1 = Mx_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);
My_all_mirror1 = My_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);
Mz_all_mirror1 = Mz_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);

t_all_mirror1 = t_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);

stroke_L_all_mirror1 = stroke_L_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);
pitch_L_all_mirror1 = pitch_L_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);
dev_L_all_mirror1 = dev_L_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);

stroke_R_all_mirror1 = stroke_R_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);
pitch_R_all_mirror1 = pitch_R_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);
dev_R_all_mirror1 = dev_R_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);

%% NO mirror
mirror_nr = 0;

file_name = [save_name,'_mirror',num2str(mirror_nr),'_'];
makeDB_TempDyn_cali_interp

%% store data & remove pre&post WBs
Fx_mean_mean_mirror0 = Fx_mean_mean(WB_start:Nwb-WB_stop);
Fy_mean_mean_mirror0 = Fy_mean_mean(WB_start:Nwb-WB_stop);
Fz_mean_mean_mirror0 = Fz_mean_mean(WB_start:Nwb-WB_stop);
Mx_mean_mean_mirror0 = Mx_mean_mean(WB_start:Nwb-WB_stop);
My_mean_mean_mirror0 = My_mean_mean(WB_start:Nwb-WB_stop);
Mz_mean_mean_mirror0 = Mz_mean_mean(WB_start:Nwb-WB_stop);

Fx_mean_all_mirror0 = Fx_mean_all(wb_start(WB_start):wb_start(Nwb-WB_stop));
Fy_mean_all_mirror0 = Fy_mean_all(wb_start(WB_start):wb_start(Nwb-WB_stop));
Fz_mean_all_mirror0 = Fz_mean_all(wb_start(WB_start):wb_start(Nwb-WB_stop));
Mx_mean_all_mirror0 = Mx_mean_all(wb_start(WB_start):wb_start(Nwb-WB_stop));
My_mean_all_mirror0 = My_mean_all(wb_start(WB_start):wb_start(Nwb-WB_stop));
Mz_mean_all_mirror0 = Mz_mean_all(wb_start(WB_start):wb_start(Nwb-WB_stop));

Fx_all_mirror0 = Fx_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);
Fy_all_mirror0 = Fy_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);
Fz_all_mirror0 = Fz_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);
Mx_all_mirror0 = Mx_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);
My_all_mirror0 = My_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);
Mz_all_mirror0 = Mz_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);

t_all_mirror0 = t_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);

stroke_L_all_mirror0 = stroke_L_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);
pitch_L_all_mirror0 = pitch_L_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);
dev_L_all_mirror0 = dev_L_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);

stroke_R_all_mirror0 = stroke_R_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);
pitch_R_all_mirror0 = pitch_R_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);
dev_R_all_mirror0 = dev_R_all(wb_start(WB_start):wb_start(Nwb-WB_stop),:);

% time
t_mean = t_mean(WB_start:Nwb-WB_stop);

%% merge data
diff = abs(stroke_L_all_mirror0 - stroke_R_all_mirror1);
if nanmax(diff(:)) == 0
    
    %% wb mean
    Fx_mean_mean = (Fx_mean_mean_mirror0 + Fx_mean_mean_mirror1)/2;
    Fy_mean_mean = (Fy_mean_mean_mirror0 - Fy_mean_mean_mirror1)/2;
    Fz_mean_mean = (Fz_mean_mean_mirror0 + Fz_mean_mean_mirror1)/2;

    Mx_mean_mean = (Mx_mean_mean_mirror0 - Mx_mean_mean_mirror1)/2;
    My_mean_mean = (My_mean_mean_mirror0 + My_mean_mean_mirror1)/2;
    Mz_mean_mean = (Mz_mean_mean_mirror0 - Mz_mean_mean_mirror1)/2;
    
    M_R_mean_mean = Mx_mean_mean * cosd(rot_axis_angle) + My_mean_mean * sind(rot_axis_angle);
    M_L_mean_mean = Mx_mean_mean * sind(rot_axis_angle) - My_mean_mean * cosd(rot_axis_angle);
    
    %% mean wb
    Fx_mean_all = (Fx_mean_all_mirror0 + Fx_mean_all_mirror1)/2;
    Fy_mean_all = (Fy_mean_all_mirror0 - Fy_mean_all_mirror1)/2;
    Fz_mean_all = (Fz_mean_all_mirror0 + Fz_mean_all_mirror1)/2;

    Mx_mean_all = (Mx_mean_all_mirror0 - Mx_mean_all_mirror1)/2;
    My_mean_all = (My_mean_all_mirror0 + My_mean_all_mirror1)/2;
    Mz_mean_all = (Mz_mean_all_mirror0 - Mz_mean_all_mirror1)/2;
    
    M_R_mean_all = Mx_mean_all * cosd(rot_axis_angle) + My_mean_all * sind(rot_axis_angle);
    M_L_mean_all = Mx_mean_all * sind(rot_axis_angle) - My_mean_all * cosd(rot_axis_angle);
    
    %% all trials
    Fx_all = (Fx_all_mirror0 + Fx_all_mirror1)/2;
    Fy_all = (Fy_all_mirror0 - Fy_all_mirror1)/2;
    Fz_all = (Fz_all_mirror0 + Fz_all_mirror1)/2;

    Mx_all = (Mx_all_mirror0 - Mx_all_mirror1)/2;
    My_all = (My_all_mirror0 + My_all_mirror1)/2;
    Mz_all = (Mz_all_mirror0 - Mz_all_mirror1)/2;
    
    M_R_all = Mx_all * cosd(rot_axis_angle) + My_all * sind(rot_axis_angle);
    M_L_all = Mx_all * sind(rot_axis_angle) - My_all * cosd(rot_axis_angle);
else
    error = 'all'
    pause
end

%% postprocess: wingbeat 1 to WB0: My=0
Mx_WBstart = nanmean(Mx_mean_mean(1:WB0));
My_WBstart = nanmean(My_mean_mean(1:WB0));
Mz_WBstart = nanmean(Mz_mean_mean(1:WB0));

Mx_mean_mean_minWBstart = Mx_mean_mean - Mx_WBstart;
Mx_mean_all_minWBstart = Mx_mean_all - Mx_WBstart;
Mx_all_minWBstart = Mx_all - Mx_WBstart;

My_mean_mean_minWBstart = My_mean_mean - My_WBstart;
My_mean_all_minWBstart = My_mean_all - My_WBstart;
My_all_minWBstart = My_all - My_WBstart;

Mz_mean_mean_minWBstart = Mz_mean_mean - Mz_WBstart;
Mz_mean_all_minWBstart = Mz_mean_all - Mz_WBstart;
Mz_all_minWBstart = Mz_all - Mz_WBstart;

M_R_mean_mean_minWBstart = Mx_mean_mean * cosd(rot_axis_angle) + My_mean_mean_minWBstart * sind(rot_axis_angle);
M_L_mean_mean_minWBstart = Mx_mean_mean * sind(rot_axis_angle) - My_mean_mean_minWBstart * cosd(rot_axis_angle);

M_R_mean_all_minWBstart = Mx_mean_all * cosd(rot_axis_angle) + My_mean_all_minWBstart * sind(rot_axis_angle);
M_L_mean_all_minWBstart = Mx_mean_all * sind(rot_axis_angle) - My_mean_all_minWBstart * cosd(rot_axis_angle);

M_R_all_minWBstart = Mx_all * cosd(rot_axis_angle) + My_all_minWBstart * sind(rot_axis_angle);
M_L_all_minWBstart = Mx_all * sind(rot_axis_angle) - My_all_minWBstart * cosd(rot_axis_angle);

%% normalize
Fx_mean_mean_norm = Fx_mean_mean / Mg_fly;
Fy_mean_mean_norm = Fy_mean_mean / Mg_fly;
Fz_mean_mean_norm = Fz_mean_mean / Mg_fly;

Mx_mean_mean_norm = Mx_mean_mean / Mg_fly / Lwing;
My_mean_mean_norm = My_mean_mean / Mg_fly / Lwing;
Mz_mean_mean_norm = Mz_mean_mean / Mg_fly / Lwing;

M_R_mean_mean_norm = M_R_mean_mean / Mg_fly / Lwing;
M_L_mean_mean_norm = M_L_mean_mean / Mg_fly / Lwing;

My_mean_mean_minWBstart_norm = My_mean_mean_minWBstart / Mg_fly / Lwing;
M_R_mean_mean_minWBstart_norm = M_R_mean_mean_minWBstart / Mg_fly / Lwing;
M_L_mean_mean_minWBstart_norm = M_L_mean_mean_minWBstart / Mg_fly / Lwing;

Fx_mean_all_norm = Fx_mean_all / Mg_fly;
Fy_mean_all_norm = Fy_mean_all / Mg_fly;
Fz_mean_all_norm = Fz_mean_all / Mg_fly;

Mx_mean_all_norm = Mx_mean_all / Mg_fly / Lwing;
My_mean_all_norm = My_mean_all / Mg_fly / Lwing;
Mz_mean_all_norm = Mz_mean_all / Mg_fly / Lwing;

M_R_mean_all_norm = M_R_mean_all / Mg_fly / Lwing;
M_L_mean_all_norm = M_L_mean_all / Mg_fly / Lwing;

My_mean_all_minWBstart_norm = My_mean_all_minWBstart / Mg_fly / Lwing;
M_R_mean_all_minWBstart_norm = M_R_mean_all_minWBstart / Mg_fly / Lwing;
M_L_mean_all_minWBstart_norm = M_L_mean_all_minWBstart / Mg_fly / Lwing;

Fx_all_norm = Fx_all / Mg_fly;
Fy_all_norm = Fy_all / Mg_fly;
Fz_all_norm = Fz_all / Mg_fly;

Mx_all_norm = Mx_all / Mg_fly / Lwing;
My_all_norm = My_all / Mg_fly / Lwing;
Mz_all_norm = Mz_all / Mg_fly / Lwing;

M_R_all_norm = M_R_all / Mg_fly / Lwing;
M_L_all_norm = M_L_all / Mg_fly / Lwing;

My_all_minWBstart_norm = My_all_minWBstart / Mg_fly / Lwing;
M_R_all_minWBstart_norm = M_R_all_minWBstart / Mg_fly / Lwing;
M_L_all_minWBstart_norm = M_L_all_minWBstart / Mg_fly / Lwing;

%% plot data
cd ..
mkdir('figs_TempDyn_cali')
cd('figs_TempDyn_cali')
plot_ForceTorque_TempDyn

%% save data
cd ..
clear ft kine t list i j file_now stroke_now step_interp...
    t_now Fx_now Fy_now Fz_now Mx_now My_now Mz_now...
          Fx_robo_now Fy_robo_now Fz_robo_now Mx_robo_now My_robo_now Mz_robo_now...
          Fx_fly Fy_fly Fz_fly Mx_fly My_fly Mz_fly...
          ...
          t_all_mirror0 t_all_mirror1...
          Fx_all_mirror0 Fy_all_mirror0 Fz_all_mirror0 Mx_all_mirror0 My_all_mirror0 Mz_all_mirror0...
          Fx_all_mirror1 Fy_all_mirror1 Fz_all_mirror1 Mx_all_mirror1 My_all_mirror1 Mz_all_mirror1...
          Fx_mean_all_mirror0 Fy_mean_all_mirror0 Fz_mean_all_mirror0 Mx_mean_all_mirror0 My_mean_all_mirror0 Mz_mean_all_mirror0...
          Fx_mean_all_mirror1 Fy_mean_all_mirror1 Fz_mean_all_mirror1 Mx_mean_all_mirror1 My_mean_all_mirror1 Mz_mean_all_mirror1...
          Fx_mean_mean_mirror0 Fy_mean_mean_mirror0 Fz_mean_mean_mirror0 Mx_mean_mean_mirror0 My_mean_mean_mirror0 Mz_mean_mean_mirror0...
          Fx_mean_mean_mirror1 Fy_mean_mean_mirror1 Fz_mean_mean_mirror1 Mx_mean_mean_mirror1 My_mean_mean_mirror1 Mz_mean_mean_mirror1...
          ...
          Fx_cali_L Fy_cali_L Fz_cali_L Mx_cali_L My_cali_L Mz_cali_L...
          Fx_cali_R Fy_cali_R Fz_cali_R Mx_cali_R My_cali_R Mz_cali_R...
          Fx_cali_L_array Fy_cali_L_array Fz_cali_L_array Mx_cali_L_array My_cali_L_array Mz_cali_L_array...
          Fx_cali_R_array Fy_cali_R_array Fz_cali_R_array Mx_cali_R_array My_cali_R_array Mz_cali_R_array...
          ...
          str_interp dev_interp rot_interp str_steps_interp dev_steps_interp rot_steps_interp...
    stroke_L_now pitch_L_now dev_L_now stroke_R_now pitch_R_now dev_R_now...
    stroke_L_all_mirror0 pitch_L_all_mirror0 dev_L_all_mirror0 stroke_R_all_mirror0 pitch_R_all_mirror0 dev_R_all_mirror0...
    stroke_L_all_mirror1 pitch_L_all_mirror1 dev_L_all_mirror1 stroke_R_all_mirror1 pitch_R_all_mirror1 dev_R_all_mirror1

save(['borf_db_',save_name,'_cali_alldata.mat'])
% save(['borf_db_',save_name,'_cali_alldata_mirror',num2str(mirror_nr),'.mat'])
