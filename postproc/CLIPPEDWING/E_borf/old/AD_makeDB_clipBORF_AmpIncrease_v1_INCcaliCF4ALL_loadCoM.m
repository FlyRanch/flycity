
clear
clc
close all

%% const

% plot_on = 0;
plot_on = 1;

% fly data
var_file = dir('flyVar*');
load(var_file.name)

% steady Wb F&M&CoM data
load('steadyWB_FnMnCoM_data.mat')

% cali data
cali_file = 'cali_matrix_interp_Full_CF.mat';
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

file_template = '*increase*.mat';

%% intact wing
cut_type_now = 0;  % no cut
cut_perc_now = 100;
dirs = dir('L_Full_R_Full_Acrylic*');

for d = 1:length(dirs)
    dir_now = dirs(d).name;
    
    if dirs(d).isdir == 1
        cd(dir_now)

        %% clipped amp increase
        file_names = dir(file_template);

        for i = 1:length(file_names)
            file_now = file_names(i).name
            load(file_now)

            m=m+1;
            cut_ratio_now  = cut_perc_now/100;
            extract_clipped_FnMnMorph_AmpIncrease_INCcali

        end
        
        %% calc linear fits for F&M versus Stroke Amplitude
        p = p+1;
        n_now = find(cut_type==cut_type_now & cut_ratio==cut_ratio_now);
        calc_fits_clipped_FnM_vs_StrokeAmp
        
        cd ..

        %% plot F-Amp & M-Amp (My@"CoG")
        if plot_on == 1
            mkdir('FnMvsStrokeAmplitude_plots')
            cd('FnMvsStrokeAmplitude_plots')
            plot_clipped_FnM_vs_StrokeAmp

            saveas(gcf,[dir_now,'.fig'])
            saveas(gcf,[dir_now,'.png'])
            plot2svg([dir_now,'.svg'])
            
            cd ..
        end
    end
end

%% Distal clipped wing
cut_type_now = 1;  % TIP
dirs = dir('L_Full_R_Distal*');

for d = 1:length(dirs)
    dir_now = dirs(d).name;
    
    if dirs(d).isdir == 1
        cd(dir_now)

        %% clipped amp increase
        file_names = dir(file_template);

        for i = 1:length(file_names)
            file_now = file_names(i).name
            load(file_now)

            m=m+1;
            cut_perc_now = str2num(dir_now(end-1:end));
            cut_ratio_now  = cut_perc_now/100;
            extract_clipped_FnMnMorph_AmpIncrease_INCcali
        end
        
        %% calc linear fits for F&M versus Stroke Amplitude
        p = p+1;
        n_now = find(cut_type==cut_type_now & cut_ratio==cut_ratio_now);
        calc_fits_clipped_FnM_vs_StrokeAmp
        
        cd ..

        %% plot F-Amp & M-Amp (My@"CoG")
        if plot_on == 1
            mkdir('FnMvsStrokeAmplitude_plots')
            cd('FnMvsStrokeAmplitude_plots')
            plot_clipped_FnM_vs_StrokeAmp

            saveas(gcf,[dir_now,'.fig'])
            saveas(gcf,[dir_now,'.png'])
            plot2svg([dir_now,'.svg'])
            
            cd ..
        end
    end
end


%% Trailing edge clipped wing
cut_type_now = 2;  % TE
dirs = dir('L_Full_R_Trailing*');

for d = 1:length(dirs)
    dir_now = dirs(d).name;
    
    if dirs(d).isdir == 1
        cd(dir_now)

        %% clipped amp increase
        file_names = dir(file_template);

        for i = 1:length(file_names)
            file_now = file_names(i).name
            load(file_now)

            m=m+1;
            cut_perc_now = str2num(dir_now(end-1:end));
            cut_ratio_now  = cut_perc_now/100;
            extract_clipped_FnMnMorph_AmpIncrease_INCcali
        end
        
        %% calc linear fits for F&M versus Stroke Amplitude
        p = p+1;
        n_now = find(cut_type==cut_type_now & cut_ratio==cut_ratio_now);
        calc_fits_clipped_FnM_vs_StrokeAmp
        
        cd ..

        %% plot F-Amp & M-Amp (My@"CoG")
        if plot_on == 1
            mkdir('FnMvsStrokeAmplitude_plots')
            cd('FnMvsStrokeAmplitude_plots')
            plot_clipped_FnM_vs_StrokeAmp

            saveas(gcf,[dir_now,'.fig'])
            saveas(gcf,[dir_now,'.png'])
            plot2svg([dir_now,'.svg'])
            
            cd ..
        end
    end
end



%% save data
save('roboflyDB_FnM_vs_StrokeAmplitude_CutWing_INCcali',...
    'cut_ratio_fit',...
    'cut_type_fit',...
    'l_ratio_fit',...
    'CoA_ratio_fit',...
    'A_ratio_fit',...
    'S1_ratio_fit',...
    'S2_ratio_fit',...
    'S3_ratio_fit',...
    ...
    'Amp_ratio',...
    'AmpStroke_intact',...
    'AmpStroke_cut',...
    ...
    'cut_ratio',...
    'cut_perc',...
    'cut_type',...
    'cut_perc_geom',...
    'cut_ratio_geom',...
    'cut_type_geom',...
    ...
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
    ...
    'Fx_norm',...
    'Fy_norm',...
    'Fz_norm',...
    'Mx_norm',...
    'My_norm',...
    'Mz_norm',...
    ...
    'Fy_norm_MinSteady',...
    'Mx_norm_MinSteady',...
    'My_norm_MinSteady',...
    'Mz_norm_MinSteady',...
    'My_norm_CoM',...
    ...
    'Fx_Amp_fit',...
    'Fx_Amp_fit_error',...
    'Fy_Amp_fit',...
    'Fy_Amp_fit_error',...
    'Fz_Amp_fit',...
    'Fz_Amp_fit_error',...
    'Fy_MinSteady_Amp_fit',...
    'Fy_MinSteady_Amp_fit_error',...
    ...
    'Mx_Amp_fit',...
    'Mx_Amp_fit_error',...
    'My_Amp_fit',...
    'My_Amp_fit_error',...
    'Mz_Amp_fit',...
    'Mz_Amp_fit_error',...
    ...
    'Mx_MinSteady_Amp_fit',...
    'Mx_MinSteady_Amp_fit_error',...
    'My_MinSteady_Amp_fit',...
    'My_MinSteady_Amp_fit_error',...
    'Mz_MinSteady_Amp_fit',...
    'Mz_MinSteady_Amp_fit_error',...
    ...
    'My_CoM_Amp_fit',...
    'My_CoM_Amp_fit_error',...
    'My_MinSteady_Amp_fit2',...
    'My_MinSteady_Amp_fit2_error',...
    'My_CoM_Amp_fit2',...
    'My_CoM_Amp_fit2_error');

