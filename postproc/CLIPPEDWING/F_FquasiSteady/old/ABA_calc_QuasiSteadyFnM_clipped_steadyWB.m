clear;
clc;
close all
warning off

load('bodyNwingModel_4qsModel.mat')
load('WBdataset_steady_1603WBs.mat')
% 
% loadname=dir('WBmod_S2S3AmpRatioFunc_*')
% loadname = loadname.name;
% load(loadname)
% 

plot_on = 1
% plot_on = 0

%% constants
nr_sect = settings.nr_chord_sect;
nr_timepoints = settings.nr_timepoints;

Mg_fly = body_model.Mg_fly;
l_wing = wing_model.length/1000;
    
t_norm = 0:1/(nr_timepoints-1):1;

%% steady WB
f_steady = f_wb_steady_meanCIstd(1);

% fouriers coeffs
stroke_coeffs_steady = stroke_steady_fourier_coeffs_binmean;
dev_coeffs_steady = dev_steady_fourier_coeffs_binmean;
rot_coeffs_steady = pitch_steady_fourier_coeffs_binmean;

[stroke_steady] = deg2rad(calc_val_fourier_series_4thN8th_order(t_norm,stroke_coeffs_steady,0));
[dev_steady] = deg2rad(calc_val_fourier_series_4thN8th_order(t_norm,dev_coeffs_steady,0));
[rot_steady] = deg2rad(calc_val_fourier_series_4thN8th_order(t_norm,rot_coeffs_steady,0));

%% WB kin now
freq = f_steady;

stroke_L = stroke_steady;
dev_L = dev_steady;
rot_L = rot_steady;

stroke_R = stroke_steady;
dev_R = dev_steady;
rot_R = rot_steady;

% transl+rot forces
rot_on=1;
[ FM_strkpln, Vel_wingtip ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_on )

%% normalize
FM_norm = FM_strkpln;

FM_norm.FM(:,1:3) = FM_norm.FM(:,1:3) / Mg_fly;
FM_norm.FM(:,4:6) = FM_norm.FM(:,1:3) / Mg_fly / l_wing;

Fz_mean = nanmean(FM_norm.FM(:,3))

% if plot_on == 1
%     mkdir('QSmodel_SteadyWB')
%     cd('QSmodel_SteadyWB')
% end
% 

plot(FM_norm.FM(:,3))

