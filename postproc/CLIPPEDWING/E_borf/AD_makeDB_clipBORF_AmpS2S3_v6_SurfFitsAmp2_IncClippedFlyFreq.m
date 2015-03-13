
clear
clc
close all

%% const
fit_type = 'poly12';

% colormap: blue to white to red
cmap_surf=cbrewer('div','RdBu',100);
cmap_surf = flipud(cmap_surf);

% % black at edges
% cmap_surf(1,:) = [0 0 0];
% cmap_surf(end,:) = [0 0 0];

% % colormap:jet
% cmap_surf = jet(100);
% cmap = flipud(cmap_surf);

plot_on = 0;
plot_on = 1;

plot_fit = 0;
% plot_fit = 1;

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
            file_now = file_names(i).name;
            load(file_now)

            m=m+1;
            cut_ratio_now  = cut_perc_now/100;
            extract_clipped_FnMnMorph_AmpIncrease_INCcali_incFcut

        end
        cd ..
        
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
            file_now = file_names(i).name;
            load(file_now)

            m=m+1;
            cut_perc_now = str2num(dir_now(end-1:end));
            cut_ratio_now  = cut_perc_now/100;
            extract_clipped_FnMnMorph_AmpIncrease_INCcali_incFcut
        end
        
        cd ..
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
            file_now = file_names(i).name;
            load(file_now)

            m=m+1;
            cut_perc_now = str2num(dir_now(end-1:end));
            cut_ratio_now  = cut_perc_now/100;
            extract_clipped_FnMnMorph_AmpIncrease_INCcali_incFcut
        end
        
        cd ..
    end
end

%% calc Forces & Torques at wingbeat frequency of cut wing flies
Fx_norm_freqMod = Fx_norm*Fnorm_clip_steady;
Fy_norm_freqMod = Fy_norm*Fnorm_clip_steady;
Fz_norm_freqMod = Fz_norm*Fnorm_clip_steady;

Mx_norm_freqMod = Mx_norm*Fnorm_clip_steady;
My_norm_freqMod = My_norm*Fnorm_clip_steady;
Mz_norm_freqMod = Mz_norm*Fnorm_clip_steady;

% min zero offset
Fy_norm_MinSteady_freqMod = Fy_norm_MinSteady*Fnorm_clip_steady;
Mx_norm_MinSteady_freqMod = Mx_norm_MinSteady*Fnorm_clip_steady;
My_norm_MinSteady_freqMod = My_norm_MinSteady*Fnorm_clip_steady;
Mz_norm_MinSteady_freqMod = Mz_norm_MinSteady*Fnorm_clip_steady;

% @CoM
My_norm_CoM_freqMod = My_norm_CoM*Fnorm_clip_steady;

% single wing
Fx_norm_cut_freqMod = Fx_norm_cut*Fnorm_clip_steady;
Fy_norm_cut_freqMod = Fy_norm_cut*Fnorm_clip_steady;
Fz_norm_cut_freqMod = Fz_norm_cut*Fnorm_clip_steady;

Fy_norm_cut_MinSteady_freqMod = Fy_norm_cut_MinSteady*Fnorm_clip_steady;

%% calc surface fits F&Ms vs Amp vs S2nS3 @ steady wb frequency
% F&M vs Amp: parabola & F&M vs S2&S3: linear

% F-Aratio-S2ratio SurfFit
[Fx_Amp_S2_SurfFit, Fx_Amp_S2_SurfFit_error] = createSurfaceFit(Amp_ratio, S2_ratio ,Fx_norm ,fit_type ,plot_fit);
[Fy_Amp_S2_SurfFit, Fy_Amp_S2_SurfFit_error] = createSurfaceFit(Amp_ratio, S2_ratio ,Fy_norm ,fit_type ,plot_fit);
[Fz_Amp_S2_SurfFit, Fz_Amp_S2_SurfFit_error] = createSurfaceFit(Amp_ratio, S2_ratio ,Fz_norm ,fit_type ,plot_fit);

[Fy_MinSteady_Amp_S2_SurfFit, Fy_MinSteady_Amp_S2_SurfFit_error] = createSurfaceFit(Amp_ratio, S2_ratio ,Fy_norm_MinSteady ,fit_type ,plot_fit);

% M-Aratio-S3ratio SurfFit
[Mx_Amp_S3_SurfFit, Mx_Amp_S3_SurfFit_error] = createSurfaceFit(Amp_ratio, S3_ratio ,Mx_norm ,fit_type ,plot_fit);
[My_Amp_S3_SurfFit, My_Amp_S3_SurfFit_error] = createSurfaceFit(Amp_ratio, S3_ratio ,My_norm ,fit_type ,plot_fit);
[Mz_Amp_S3_SurfFit, Mz_Amp_S3_SurfFit_error] = createSurfaceFit(Amp_ratio, S3_ratio ,Mz_norm ,fit_type ,plot_fit);

[Mx_MinSteady_Amp_S3_SurfFit, Mx_MinSteady_Amp_S3_SurfFit_error] = createSurfaceFit(Amp_ratio, S3_ratio ,Mx_norm_MinSteady ,fit_type ,plot_fit);
[My_MinSteady_Amp_S3_SurfFit, My_MinSteady_Amp_S3_SurfFit_error] = createSurfaceFit(Amp_ratio, S3_ratio ,My_norm_MinSteady ,fit_type ,plot_fit);
[Mz_MinSteady_Amp_S3_SurfFit, Mz_MinSteady_Amp_S3_SurfFit_error] = createSurfaceFit(Amp_ratio, S3_ratio ,Mz_norm_MinSteady ,fit_type ,plot_fit);

[My_CoM_Amp_S3_SurfFit, My_CoM_Amp_S3_SurfFit_error] = createSurfaceFit(Amp_ratio, S3_ratio ,My_norm_CoM ,fit_type ,plot_fit);

%% calc surface fits F&Ms vs Amp vs S2nS3 @ clipped fly wb frequency
% F&M vs Amp: parabola & F&M vs S2&S3: linear

% F-Aratio-S2ratio SurfFit
[Fx_Amp_S2_SurfFit_freqMod, Fx_Amp_S2_SurfFit_error_freqMod] = createSurfaceFit(Amp_ratio, S2_ratio ,Fx_norm_freqMod ,fit_type ,plot_fit);
[Fy_Amp_S2_SurfFit_freqMod, Fy_Amp_S2_SurfFit_error_freqMod] = createSurfaceFit(Amp_ratio, S2_ratio ,Fy_norm_freqMod ,fit_type ,plot_fit);
[Fz_Amp_S2_SurfFit_freqMod, Fz_Amp_S2_SurfFit_error_freqMod] = createSurfaceFit(Amp_ratio, S2_ratio ,Fz_norm_freqMod ,fit_type ,plot_fit);

[Fy_MinSteady_Amp_S2_SurfFit_freqMod, Fy_MinSteady_Amp_S2_SurfFit_error_freqMod] = createSurfaceFit(Amp_ratio, S2_ratio ,Fy_norm_MinSteady_freqMod ,fit_type ,plot_fit);

% M-Aratio-S3ratio SurfFit
[Mx_Amp_S3_SurfFit_freqMod, Mx_Amp_S3_SurfFit_error_freqMod] = createSurfaceFit(Amp_ratio, S3_ratio ,Mx_norm_freqMod ,fit_type ,plot_fit);
[My_Amp_S3_SurfFit_freqMod, My_Amp_S3_SurfFit_error_freqMod] = createSurfaceFit(Amp_ratio, S3_ratio ,My_norm_freqMod ,fit_type ,plot_fit);
[Mz_Amp_S3_SurfFit_freqMod, Mz_Amp_S3_SurfFit_error_freqMod] = createSurfaceFit(Amp_ratio, S3_ratio ,Mz_norm_freqMod ,fit_type ,plot_fit);

[Mx_MinSteady_Amp_S3_SurfFit_freqMod, Mx_MinSteady_Amp_S3_SurfFit_error_freqMod] = createSurfaceFit(Amp_ratio, S3_ratio ,Mx_norm_MinSteady_freqMod ,fit_type ,plot_fit);
[My_MinSteady_Amp_S3_SurfFit_freqMod, My_MinSteady_Amp_S3_SurfFit_error_freqMod] = createSurfaceFit(Amp_ratio, S3_ratio ,My_norm_MinSteady_freqMod ,fit_type ,plot_fit);
[Mz_MinSteady_Amp_S3_SurfFit_freqMod, Mz_MinSteady_Amp_S3_SurfFit_error_freqMod] = createSurfaceFit(Amp_ratio, S3_ratio ,Mz_norm_MinSteady_freqMod ,fit_type ,plot_fit);

[My_CoM_Amp_S3_SurfFit_freqMod, My_CoM_Amp_S3_SurfFit_error_freqMod] = createSurfaceFit(Amp_ratio, S3_ratio ,My_norm_CoM_freqMod ,fit_type ,plot_fit);


%% plot surface fits F&Ms vs Amp vs S2nS3
if plot_on == 1

    mkdir('figures_cutWing_robofly')
    cd('figures_cutWing_robofly')
    
    % @ steady wb frequency
    close all
    plot_clipped_FnM_vs_StrokeAmp_vs_S2nS3_SurfFits

    figure(1)
    saveas(gcf,'ForcesVsStrokeAmpVsS2_SurfFits_freqSteadyWB.fig')
    saveas(gcf,'ForcesVsStrokeAmpVsS2_SurfFits_freqSteadyWB.png')
    plot2svg('ForcesVsStrokeAmpVsS2_SurfFits_freqSteadyWB.svg')

    figure(2)
    saveas(gcf,'TorqueVsStrokeAmpVsS3_SurfFits_minSteady_MyAtCoM_freqSteadyWB.fig')
    saveas(gcf,'TorqueVsStrokeAmpVsS3_SurfFits_minSteady_MyAtCoM_freqSteadyWB.png')
    plot2svg('TorqueVsStrokeAmpVsS3_SurfFits_minSteady_MyAtCoM_freqSteadyWB.svg')

    figure(3)
    saveas(gcf,'FzMx_vs_StrokeAmp_vs_S2S3_SurfFits_MxminSteady_freqSteadyWB.fig')
    saveas(gcf,'FzMx_vs_StrokeAmp_vs_S2S3_SurfFits_MxminSteady_freqSteadyWB.png')
    plot2svg('FzMx_vs_StrokeAmp_vs_S2S3_SurfFits_MxminSteady_freqSteadyWB.svg')
    
    % @ clipped fly wb frequency
    close all
    plot_clipped_FnM_vs_AmpnS2nS3_SurfFits_ClipFlyWBfreq
    
    figure(1)
    saveas(gcf,'ForcesVsStrokeAmpVsS2_SurfFits_ClipFlyWBfreq.fig')
    saveas(gcf,'ForcesVsStrokeAmpVsS2_SurfFits_ClipFlyWBfreq.png')
    plot2svg('ForcesVsStrokeAmpVsS2_SurfFits_ClipFlyWBfreq.svg')

    figure(2)
    saveas(gcf,'TorqueVsStrokeAmpVsS3_SurfFits_minSteady_MyAtCoM_ClipFlyWBfreq.fig')
    saveas(gcf,'TorqueVsStrokeAmpVsS3_SurfFits_minSteady_MyAtCoM_ClipFlyWBfreq.png')
    plot2svg('TorqueVsStrokeAmpVsS3_SurfFits_minSteady_MyAtCoM_ClipFlyWBfreq.svg')

    figure(3)
    saveas(gcf,'FzMx_vs_StrokeAmp_vs_S2S3_SurfFits_MxminSteady_ClipFlyWBfreq.fig')
    saveas(gcf,'FzMx_vs_StrokeAmp_vs_S2S3_SurfFits_MxminSteady_ClipFlyWBfreq.png')
    plot2svg('FzMx_vs_StrokeAmp_vs_S2S3_SurfFits_MxminSteady_ClipFlyWBfreq.svg')

    cd ..
end



%% save data
save('roboflyDB_CutWing_FnM_vs_StrokeAmplitude_vs_S2nS3_INCcaliCF',...
    'fit_type',...
    ...
    'Fx_Amp_S2_SurfFit' ,'Fx_Amp_S2_SurfFit_error',...
    'Fy_Amp_S2_SurfFit' ,'Fy_Amp_S2_SurfFit_error',...
    'Fz_Amp_S2_SurfFit' ,'Fz_Amp_S2_SurfFit_error',...
    'Fy_MinSteady_Amp_S2_SurfFit' ,'Fy_MinSteady_Amp_S2_SurfFit_error',...
    'Mx_Amp_S3_SurfFit' ,'Mx_Amp_S3_SurfFit_error',...
    'My_Amp_S3_SurfFit' ,'My_Amp_S3_SurfFit_error',...
    'Mz_Amp_S3_SurfFit' ,'Mz_Amp_S3_SurfFit_error',...
    'Mx_MinSteady_Amp_S3_SurfFit' ,'Mx_MinSteady_Amp_S3_SurfFit_error',...
    'My_MinSteady_Amp_S3_SurfFit' ,'My_MinSteady_Amp_S3_SurfFit_error',...
    'Mz_MinSteady_Amp_S3_SurfFit' ,'Mz_MinSteady_Amp_S3_SurfFit_error',...
    'My_CoM_Amp_S3_SurfFit' ,'My_CoM_Amp_S3_SurfFit_error',...
    ...
    'Fx_Amp_S2_SurfFit_freqMod' ,'Fx_Amp_S2_SurfFit_error_freqMod',...
    'Fy_Amp_S2_SurfFit_freqMod' ,'Fy_Amp_S2_SurfFit_error_freqMod',...
    'Fz_Amp_S2_SurfFit_freqMod' ,'Fz_Amp_S2_SurfFit_error_freqMod',...
    'Fy_MinSteady_Amp_S2_SurfFit_freqMod' ,'Fy_MinSteady_Amp_S2_SurfFit_error_freqMod',...
    'Mx_Amp_S3_SurfFit_freqMod' ,'Mx_Amp_S3_SurfFit_error_freqMod',...
    'My_Amp_S3_SurfFit_freqMod' ,'My_Amp_S3_SurfFit_error_freqMod',...
    'Mz_Amp_S3_SurfFit_freqMod' ,'Mz_Amp_S3_SurfFit_error_freqMod',...
    'Mx_MinSteady_Amp_S3_SurfFit_freqMod' ,'Mx_MinSteady_Amp_S3_SurfFit_error_freqMod',...
    'My_MinSteady_Amp_S3_SurfFit_freqMod' ,'My_MinSteady_Amp_S3_SurfFit_error_freqMod',...
    'Mz_MinSteady_Amp_S3_SurfFit_freqMod' ,'Mz_MinSteady_Amp_S3_SurfFit_error_freqMod',...
    'My_CoM_Amp_S3_SurfFit_freqMod' ,'My_CoM_Amp_S3_SurfFit_error_freqMod',...
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
    'Fx_norm_cut',...
    'Fy_norm_cut',...
    'Fz_norm_cut',...
    'Fy_norm_cut_MinSteady');

