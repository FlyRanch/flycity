
clear
clc
close all
warning off

%% const
% fly data
var_file = dir('flyVar*');
load(var_file.name)

% steadyWB F&M&CoM data
load('steadyWB_FnMnCoM_data.mat')

% Fenhance_by_freqMod_4SteadyWB
load('Fenhance_by_freqMod_4SteadyWB.mat')

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

%% cali CF for all (others have ERROR: Mx is assymetric)
% cali_template = 'cali_matrix_interp_*.mat';
cali_file = 'cali_matrix_interp_L_Full_R_Full_CF.mat';
load(cali_file)

%% load & save intact
cut_type_now = 0;  % no cut

% dirs = dir('L_Full_R_Full*')
% dirs = dir('L_Full_R_Full_CF*')
dirs = dir('L_Full_R_Full_Acrylic*');

% file_template = 'F_*_0.mat';
% file_template = 'F_*increase_0.mat';
file_template = 'F_*reduce_0.mat';

for i = 1:length(dirs)
    dir_now = dirs(i);
    if dir_now.isdir == 1
        cd(dir_now.name)
        
        file_name = dir(file_template);
        load(file_name(1).name)
        
        m=m+1;
        cut_perc_now = 100;
        
        cut_ratio(m,1) = cut_perc_now/100;
        cut_perc(m,1) = cut_perc_now;
        cut_type(m,1) = cut_type_now;

        n_total = length(t);
        n_wb = n_total/Nwb;
        n_start = round(wb_start*n_wb+1);
        n_stop = round(wb_stop*n_wb);
        
        % borf geometry
        n_geom = find(BorfMorphCutData.cut_type == cut_type_now);
        
        cut_perc_geom(m,1) = BorfMorphCutData.cut_perc(n_geom);
        cut_ratio_geom(m,1) = BorfMorphCutData.cut_ratio(n_geom);
        cut_type_geom(m,1) = BorfMorphCutData.cut_type(n_geom);
        
        l_ratio(m,1) = BorfMorphCutData.WingLength_ratio(n_geom);
        CoA_ratio(m,1) = BorfMorphCutData.CoA_ratio(n_geom);
        A_ratio(m,1) = BorfMorphCutData.WingArea_ratio(n_geom);
        S1_ratio(m,1) =BorfMorphCutData.FirstMoment_ratio(n_geom);
        S2_ratio(m,1) =BorfMorphCutData.SecondMoment_ratio(n_geom);
        S3_ratio(m,1) =BorfMorphCutData.ThirdMoment_ratio(n_geom);
        
        CoA_normL(m,1) = BorfMorphCutData.CoA_norm(n_geom);
        A_normL(m,1) = BorfMorphCutData.WingArea_norm(n_geom);
        S1_normL(m,1) =BorfMorphCutData.FirstMoment_norm(n_geom);
        S2_normL(m,1) =BorfMorphCutData.SecondMoment_norm(n_geom);
        S3_normL(m,1) =BorfMorphCutData.ThirdMoment_norm(n_geom);
        
        l_normA(m,1) = BorfMorphCutData.WingLength_normA(n_geom);
        CoA_normA(m,1) = BorfMorphCutData.CoA_normA(n_geom);
        S1_normA(m,1) =BorfMorphCutData.FirstMoment_normA(n_geom);
        S2_normA(m,1) =BorfMorphCutData.SecondMoment_normA(n_geom);
        S3_normA(m,1) =BorfMorphCutData.ThirdMoment_normA(n_geom);

        extract_clipped_FnM_INCcali
        
        cd ..

    end
end

%% load & save distal data
cut_type_now = 1;  % TIP

dirs = dir('L_Full_R_Distal*');

file_template = 'F_*increase_0.mat';

for i = 1:length(dirs)
    dir_now = dirs(i);
    if dir_now.isdir == 1
        cd(dir_now.name)
        
        file_name = dir(file_template);
        load(file_name.name)
        
        m=m+1;
        cut_perc_now = str2num(dir_now.name(end-1:end));
        
        cut_ratio(m,1) = cut_perc_now/100;
        cut_perc(m,1) = cut_perc_now;
        cut_type(m,1) = cut_type_now;

        n_total = length(t);
        n_wb = n_total/Nwb;
        n_start = round(wb_start*n_wb+1);
        n_stop = round(wb_stop*n_wb);

        % borf geometry
        n_geom = find(BorfMorphCutData.cut_type == cut_type_now & BorfMorphCutData.cut_ratio == cut_perc(m,1)/100);
        
        cut_perc_geom(m,1) = BorfMorphCutData.cut_perc(n_geom);
        cut_ratio_geom(m,1) = BorfMorphCutData.cut_ratio(n_geom);
        cut_type_geom(m,1) = BorfMorphCutData.cut_type(n_geom);
        
        l_ratio(m,1) = BorfMorphCutData.WingLength_ratio(n_geom);
        CoA_ratio(m,1) = BorfMorphCutData.CoA_ratio(n_geom);
        A_ratio(m,1) = BorfMorphCutData.WingArea_ratio(n_geom);
        S1_ratio(m,1) =BorfMorphCutData.FirstMoment_ratio(n_geom);
        S2_ratio(m,1) =BorfMorphCutData.SecondMoment_ratio(n_geom);
        S3_ratio(m,1) =BorfMorphCutData.ThirdMoment_ratio(n_geom);
        
        CoA_normL(m,1) = BorfMorphCutData.CoA_norm(n_geom);
        A_normL(m,1) = BorfMorphCutData.WingArea_norm(n_geom);
        S1_normL(m,1) =BorfMorphCutData.FirstMoment_norm(n_geom);
        S2_normL(m,1) =BorfMorphCutData.SecondMoment_norm(n_geom);
        S3_normL(m,1) =BorfMorphCutData.ThirdMoment_norm(n_geom);
        
        l_normA(m,1) = BorfMorphCutData.WingLength_normA(n_geom);
        CoA_normA(m,1) = BorfMorphCutData.CoA_normA(n_geom);
        S1_normA(m,1) =BorfMorphCutData.FirstMoment_normA(n_geom);
        S2_normA(m,1) =BorfMorphCutData.SecondMoment_normA(n_geom);
        S3_normA(m,1) =BorfMorphCutData.ThirdMoment_normA(n_geom);

        extract_clipped_FnM_INCcali
        
        cd ..

    end
end

%% load & save trailing edge data
cut_type_now = 2;  % TE

dirs = dir('L_Full_R_Trailing*');

file_template = 'F_*increase_0.mat';

for i = 1:length(dirs)
    dir_now = dirs(i);
    if dir_now.isdir == 1
        cd(dir_now.name)
        
        file_name = dir(file_template);
        load(file_name.name)
        
        m=m+1;
        cut_perc_now = str2num(dir_now.name(end-1:end));
        
        cut_ratio(m,1) = cut_perc_now/100;
        cut_perc(m,1) = cut_perc_now;
        cut_type(m,1) = cut_type_now;

        n_total = length(t);
        n_wb = n_total/Nwb;
        n_start = round(wb_start*n_wb+1);
        n_stop = round(wb_stop*n_wb);

        % borf geometry
        n_geom = find(BorfMorphCutData.cut_type == cut_type_now & BorfMorphCutData.cut_ratio == cut_perc(m,1)/100);
        
        cut_perc_geom(m,1) = BorfMorphCutData.cut_perc(n_geom);
        cut_ratio_geom(m,1) = BorfMorphCutData.cut_ratio(n_geom);
        cut_type_geom(m,1) = BorfMorphCutData.cut_type(n_geom);
        
        l_ratio(m,1) = BorfMorphCutData.WingLength_ratio(n_geom);
        CoA_ratio(m,1) = BorfMorphCutData.CoA_ratio(n_geom);
        A_ratio(m,1) = BorfMorphCutData.WingArea_ratio(n_geom);
        S1_ratio(m,1) =BorfMorphCutData.FirstMoment_ratio(n_geom);
        S2_ratio(m,1) =BorfMorphCutData.SecondMoment_ratio(n_geom);
        S3_ratio(m,1) =BorfMorphCutData.ThirdMoment_ratio(n_geom);
        
        CoA_normL(m,1) = BorfMorphCutData.CoA_norm(n_geom);
        A_normL(m,1) = BorfMorphCutData.WingArea_norm(n_geom);
        S1_normL(m,1) =BorfMorphCutData.FirstMoment_norm(n_geom);
        S2_normL(m,1) =BorfMorphCutData.SecondMoment_norm(n_geom);
        S3_normL(m,1) =BorfMorphCutData.ThirdMoment_norm(n_geom);
        
        l_normA(m,1) = BorfMorphCutData.WingLength_normA(n_geom);
        CoA_normA(m,1) = BorfMorphCutData.CoA_normA(n_geom);
        S1_normA(m,1) =BorfMorphCutData.FirstMoment_normA(n_geom);
        S2_normA(m,1) =BorfMorphCutData.SecondMoment_normA(n_geom);
        S3_normA(m,1) =BorfMorphCutData.ThirdMoment_normA(n_geom);
        
        extract_clipped_FnM_INCcali
        
        cd ..

    end
end

%% calc My at CoM;
My_norm_CoM = My_norm - d_norm_steady*Fz_norm;

%% calc Forces & Torques at wingbeat frequency of cut wing flies
Fx_norm_freqMod = Fx_norm*Fnorm_clip_steady;
Fy_norm_freqMod = Fy_norm*Fnorm_clip_steady;
Fz_norm_freqMod = Fz_norm*Fnorm_clip_steady;

Mx_norm_freqMod = Mx_norm*Fnorm_clip_steady;
My_norm_freqMod = My_norm*Fnorm_clip_steady;
Mz_norm_freqMod = Mz_norm*Fnorm_clip_steady;

My_norm_freqMod_CoM = My_norm_CoM*Fnorm_clip_steady;

%% calc linear fits
calcLinearFits_FvsS2nMvsS3_steadyFreq
calcLinearFits_FvsS2nMvsS3_FreqModClippedFly

%% plot F-S2 & M-S3 (My@"CoM") for steady freq & cut fly freq
mkdir('figures_cutWing_robofly')
cd('figures_cutWing_robofly')

plot_FvsS2nMvsS3_FreqSteady_n_FreqModClippedFly

saveas(gcf,['FvsS2_MvsS3_robofly_CutWing_steadyWB_INCcali_MyAtCoM_freqSteadyNclippedFly.fig'])
saveas(gcf,['FvsS2_MvsS3_robofly_CutWing_steadyWB_INCcali_MyAtCoM_freqSteadyNclippedFly.png'])
saveas(gcf,['FvsS2_MvsS3_robofly_CutWing_steadyWB_INCcali_MyAtCoM_freqSteadyNclippedFly.svg'])
% plot2svg(['FvsS2_MvsS3_robofly_CutWing_steadyWB_INCcali_MyAtCoM_freqSteadyNclippedFly.svg'])

cd ..

%% save data
save('roboflyDB_CutWing_steadyWB_INCcaliCF_INCclippedWingFreq',...
    'freq_clip_mean',...
    'freq_steady',...
    'time',...
    'time_norm',...
    'fx_norm',...
    'fy_norm',...
    'fz_norm',...
    'mx_norm',...
    'my_norm',...
    'mz_norm',...
    'Fx_norm',...
    'Fy_norm',...
    'Fz_norm',...
    'Mx_norm',...
    'My_norm',...
    'Mz_norm',...
    'My_norm_CoM',...
    'Fx_norm_freqMod',...
    'Fy_norm_freqMod',...
    'Fz_norm_freqMod',...
    'Mx_norm_freqMod',...
    'My_norm_freqMod',...
    'Mz_norm_freqMod',...
    'My_norm_freqMod_CoM',...
    'cut_ratio',...
    'cut_perc',...
    'cut_type',...
    'cut_perc_geom',...
    'cut_ratio_geom',...
    'cut_type_geom',...
    'l_ratio',...
    'CoA_ratio',...
    'A_ratio',...
    'S1_ratio',...
    'S2_ratio',...
    'S3_ratio',...
    'CoA_normL',...
    'A_normL',...
    'S1_normL',...
    'S2_normL',...
    'S3_normL',...
    'CoA_normA',...
    'l_normA',...
    'S1_normA',...
    'S2_normA',...
    'S3_normA',...
    'Fx_S2_fit',...
    'Fx_S2_fit_error',...
    'Fy_S2_fit',...
    'Fy_S2_fit_error',...
    'Fz_S2_fit',...
    'Fz_S2_fit_error',...
    'Mx_S3_fit',...
    'Mx_S3_fit_error',...
    'My_S3_fit',...
    'My_S3_fit_error',...
    'Mz_S3_fit',...
    'Mz_S3_fit_error',...
    'MxMinSteady_S3_fit',...
    'MxMinSteady_S3_fit_error',...
    'MyMinSteady_S3_fit',...
    'MyMinSteady_S3_fit_error',...
    'MzMinSteady_S3_fit',...
    'MzMinSteady_S3_fit_error',...
    'My_CoM_S3_fit',...
    'My_CoM_S3_fit_error',...
    'MyMinSteady_CoM_S3_fit',...
    'MyMinSteady_CoM_S3_fit_error',...
    'Fx_S2_fit_freqMod',...
    'Fx_S2_fit_error_freqMod',...
    'Fy_S2_fit_freqMod',...
    'Fy_S2_fit_error_freqMod',...
    'Fz_S2_fit_freqMod',...
    'Fz_S2_fit_error_freqMod',...
    'Mx_S3_fit_freqMod',...
    'Mx_S3_fit_error_freqMod',...
    'My_S3_fit_freqMod',...
    'My_S3_fit_error_freqMod',...
    'Mz_S3_fit_freqMod',...
    'Mz_S3_fit_error_freqMod',...
    'MxMinSteady_S3_fit_freqMod',...
    'MxMinSteady_S3_fit_error_freqMod',...
    'MyMinSteady_S3_fit_freqMod',...
    'MyMinSteady_S3_fit_error_freqMod',...
    'MzMinSteady_S3_fit_freqMod',...
    'MzMinSteady_S3_fit_error_freqMod',...
    'My_CoM_S3_fit_freqMod',...
    'My_CoM_S3_fit_error_freqMod',...
    'MyMinSteady_CoM_S3_fit_freqMod',...
    'MyMinSteady_CoM_S3_fit_error_freqMod');





