
clear
clc
close all
warning off

%% const
plot_on = 0;
plot_on = 1;

plot_fit = 0;
% plot_fit = 1;

%% colormap for Aratio: blue to white to red
cmap_RdBu =cbrewer('div','RdBu',100);
cmap_BuRd = flipud(cmap_RdBu);

% black at edges of cmap
cmap_BuRd_edge = cmap_BuRd;
cmap_BuRd_edge(1,:) = [0 0 0];
cmap_BuRd_edge(end,:) = [0 0 0];

cmap_Aratio = cmap_BuRd;

%% colormap for A+: reds
% cmap_Reds =cbrewer('seq','Reds',100);
% cmap_AdAi = cmap_Reds;
% 
% cmap_YlOrRd =cbrewer('seq','YlOrRd',100);
% cmap_AdAi = cmap_YlOrRd;

cmap_hot =colormap(hot(100));
cmap_hot = flipud(cmap_hot);
cmap_AdAi = cmap_hot;

% %% colormap for A+: pink to white to green
% cmap_BrBG =cbrewer('div','BrBG',100);
% cmap_BGBr = flipud(cmap_BrBG);
% 
% % black at edges of cmap
% cmap_BGBr_edge = cmap_BGBr;
% cmap_BGBr_edge(1,:) = [0 0 0];
% cmap_BGBr_edge(end,:) = [0 0 0];
% 
% cmap_AdAi = cmap_BGBr;

%% colormap for Forces: purple to white to green
cmap_PRGn =cbrewer('div','PRGn',100);
cmap_GnPR = flipud(cmap_PRGn);

% black at edges of cmap
cmap_GnPR_edge = cmap_GnPR;
cmap_GnPR_edge(1,:) = [0 0 0];
cmap_GnPR_edge(end,:) = [0 0 0];

cmap_F = cmap_GnPR;
cmap_F_neg = cmap_PRGn;

%% colormap for Torques: brown to white to turquise
cmap_PiYG =cbrewer('div','PiYG',100);
cmap_YGPi = flipud(cmap_PiYG);

% black at edges of cmap
cmap_YGPi_edge = cmap_YGPi;
cmap_YGPi_edge(1,:) = [0 0 0];
cmap_YGPi_edge(end,:) = [0 0 0];

cmap_T = cmap_YGPi;

%% load data
% fly data
var_file = dir('flyVar*');
load(var_file.name)

% steady Wb F&M&CoM data
load('steadyWB_FnMnCoM_data.mat')

% Fenhance_by_freqMod_4SteadyWB
load('Fenhance_by_freqMod_4SteadyWB.mat')

%% cali CF for all (others have ERROR: Mx is assymetric)
% cali_template = 'cali_matrix_interp_*.mat';
cali_file = 'cali_matrix_interp_L_Full_R_Full_CF.mat';
load(cali_file)

% robo data
Lrobo = .23;
crobo = .065;
ARrobo = Lrobo/crobo;

% fluid data
v_air = 1.568e-5; %kinematic viscocity in m^2/s
v_oil = 115.0*1e-6; %kinematic viscocity in m^2/s
% rho_air = 1.22521; %kg/m^3
rho_oil = 880.0; %kg/m^3

%% Elzinga scaling
% l_scale_up = crobo/c_fly;
c_fly = Lwing/ARwing_fly; %mean chord of actual fly in m
f_robo2fly = (v_air*crobo*Lrobo)/(c_fly*Lwing*v_oil);
F_robo2fly = (rho_air/rho_oil)*((c_fly/crobo)^4)*((f_robo2fly)^2);
M_robo2fly = (rho_air/rho_oil)*((c_fly/crobo)^5)*((f_robo2fly)^2);

%% cut wing geometry (FROM BW IMAGES OF SOLIDWORKS PROJECTIONS)
geom_file = dir('BorfMorphCutDatabase*');
load(geom_file.name)
BorfMorphCutData.cut_ratio = BorfMorphCutData.cut_perc/100;

%% run data
Nwb = 10;
wb_start = 2;
wb_stop = 10;
m=0;
p=0;

%% intact wing amp decrease
file_template = '*reduce*.mat';
Aincrease = 0;

cut_type_now = 0;  % no cut
cut_perc_now = 100;
dirs = dir('L_Full_R_Full_Acrylic*');

for d = 1:length(dirs)
    dir_now = dirs(d).name;
    
    if dirs(d).isdir == 1
        cd(dir_now)

        file_names = dir(file_template);

        for i = 1:length(file_names)
            file_now = file_names(i).name;
            load(file_now)

            m=m+1;
            cut_ratio_now  = cut_perc_now/100;
            extract_clipped_FnMnMorph_AmpIncrease_INCcali_incFcut_tempDyn

        end
        cd ..
        
    end
end

%% intact wing amp increase
file_template = '*increase*.mat';
Aincrease = 1;

cut_type_now = 0;  % no cut
cut_perc_now = 100;
dirs = dir('L_Full_R_Full_Acrylic*');

for d = 1:length(dirs)
    dir_now = dirs(d).name;
    
    if dirs(d).isdir == 1
        cd(dir_now)

        file_names = dir(file_template);

        for i = 1:length(file_names)
            file_now = file_names(i).name;
            load(file_now)

            m=m+1;
            cut_ratio_now  = cut_perc_now/100;
            extract_clipped_FnMnMorph_AmpIncrease_INCcali_incFcut_tempDyn

        end
        cd ..
        
    end
end

%% Distal clipped wing
file_template = '*increase*.mat';
Aincrease = 1;

cut_type_now = 1;  % TIP
dirs = dir('L_Full_R_Distal*');

for d = 1:length(dirs)
    dir_now = dirs(d).name;
    
    if dirs(d).isdir == 1
        cd(dir_now)

        %% clipped amp increase
        file_names = dir(file_template);

        for i = 1:length(file_names)
            file_now = file_names(i).name;
            load(file_now)

            m=m+1;
            cut_perc_now = str2num(dir_now(end-1:end));
            cut_ratio_now  = cut_perc_now/100;
            extract_clipped_FnMnMorph_AmpIncrease_INCcali_incFcut_tempDyn
        end
        
        cd ..
    end
end


%% Trailing edge clipped wing
file_template = '*increase*.mat';
Aincrease = 1;

cut_type_now = 2;  % TE
dirs = dir('L_Full_R_Trailing*');

for d = 1:length(dirs)
    dir_now = dirs(d).name;
    
    if dirs(d).isdir == 1
        cd(dir_now)

        %% clipped amp increase
        file_names = dir(file_template);

        for i = 1:length(file_names)
            file_now = file_names(i).name;
            load(file_now)

            m=m+1;
            cut_perc_now = str2num(dir_now(end-1:end));
            cut_ratio_now  = cut_perc_now/100;
            extract_clipped_FnMnMorph_AmpIncrease_INCcali_incFcut_tempDyn
        end
        
        cd ..
    end
end

Time_norm = [0:Nwb/(size(time_fly,1)-1):Nwb]';

%% save data
Stroke_Amplitude_ratio = Amp_ratio;
WingLength_ratio = l_ratio;
WingArea_ratio = A_ratio;

save('FRUITFLY_WINGDAMAGE_ROBOTICFLY_DATABASE','Time_norm',...
    'Fx_norm','Fy_norm','Fz_norm',...
    'Mx_norm','My_norm','Mz_norm',...
    'Stroke_IntactWing','Deviation_IntactWing','Rotation_IntactWing',...
    'Stroke_CutWing','Deviation_CutWing','Rotation_CutWing',...
    ...
    'Stroke_Amplitude_ratio',...
    'cut_ratio',...
    'cut_perc',...
    'cut_type',...
    ...
    'WingLength_ratio',...
    'WingArea_ratio',...
    'S1_ratio',...
    'S2_ratio',...
    'S3_ratio')

