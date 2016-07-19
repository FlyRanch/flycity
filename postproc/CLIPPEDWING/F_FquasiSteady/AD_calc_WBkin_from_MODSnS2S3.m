clear
clc

S2ratio = .5;
S3ratio = .5;

nr_timepoints = 200;
t_norm = 0:1/(nr_timepoints-1):1;

% save('fruitfly_damagedwing_kinMODS.mat',...
%     'solAdAiRatio','S2','S3',...
%     'f_wb_steady_meanCIstd','freq_S2S3AmpRatioFunc_steadyWBs_asympFit10',...
%     'stroke_steady_fourier_coeffs_binmean','dev_steady_fourier_coeffs_binmean','pitch_steady_fourier_coeffs_binmean',...
%     'strokeMOD_L_S2S3AmpRatioFunc_fourier_coeffs_binmean','devMOD_L_S2S3AmpRatioFunc_fourier_coeffs_binmean','pitchMOD_L_S2S3AmpRatioFunc_fourier_coeffs_binmean',...
%     'strokeMOD_R_S2S3AmpRatioFunc_fourier_coeffs_binmean','devMOD_R_S2S3AmpRatioFunc_fourier_coeffs_binmean','pitchMOD_R_S2S3AmpRatioFunc_fourier_coeffs_binmean')
    
load('fruitfly_damagedwing_kinMODS.mat')

% save('fruitfly_damagedwing_tipcut_params.mat','span_ratio','S2ratio','S3ratio','S2S3AmpRatioFunc_now','freq_now')
load('fruitfly_damagedwing_tipcut_params.mat')

%% steady WB
f_steady = f_wb_steady_meanCIstd(1);

stroke_coeffs_steady = stroke_steady_fourier_coeffs_binmean;
dev_coeffs_steady = dev_steady_fourier_coeffs_binmean;
rot_coeffs_steady = pitch_steady_fourier_coeffs_binmean;

[stroke_steady] = (calc_val_fourier_series_4thN8th_order(t_norm,stroke_coeffs_steady,0))';
[dev_steady] = (calc_val_fourier_series_4thN8th_order(t_norm,dev_coeffs_steady,0))';
[rot_steady] = (calc_val_fourier_series_4thN8th_order(t_norm,rot_coeffs_steady,0))';

sol = subs(solAdAiRatio,S2,1);
sol = subs(sol,S3,1);
S2S3AmpRatioFunc_NONclipped = eval(sol);

%% damaged fly WB kin
% MODS
strokeMOD_coeffs_intact = strokeMOD_L_S2S3AmpRatioFunc_fourier_coeffs_binmean;
devMOD_coeffs_intact = devMOD_L_S2S3AmpRatioFunc_fourier_coeffs_binmean;
rotMOD_coeffs_intact = pitchMOD_L_S2S3AmpRatioFunc_fourier_coeffs_binmean;

[strokeMOD_intact] = (calc_val_fourier_series_4thN8th_order(t_norm,strokeMOD_coeffs_intact,0))';
[devMOD_intact] = (calc_val_fourier_series_4thN8th_order(t_norm,devMOD_coeffs_intact,0))';
[rotMOD_intact] = (calc_val_fourier_series_4thN8th_order(t_norm,rotMOD_coeffs_intact,0))';

strokeMOD_coeffs_damaged = strokeMOD_R_S2S3AmpRatioFunc_fourier_coeffs_binmean;
devMOD_coeffs_damaged = devMOD_R_S2S3AmpRatioFunc_fourier_coeffs_binmean;
rotMOD_coeffs_damaged = pitchMOD_R_S2S3AmpRatioFunc_fourier_coeffs_binmean;

[strokeMOD_damaged] = (calc_val_fourier_series_4thN8th_order(t_norm,strokeMOD_coeffs_damaged,0))';
[devMOD_damaged] = (calc_val_fourier_series_4thN8th_order(t_norm,devMOD_coeffs_damaged,0))';
[rotMOD_damaged] = (calc_val_fourier_series_4thN8th_order(t_norm,rotMOD_coeffs_damaged,0))';

% MOD scaling based on S2ratio & S3ratio
% sol = subs(solAdAiRatio,S2,S2ratio);
% sol = subs(sol,S3,S3ratio);
% S2S3AmpRatioFunc_now = eval(sol);

% damaged fly WB kin
% freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit10,S2S3AmpRatioFunc_now);

stroke_intact_now   = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * strokeMOD_intact + stroke_steady;    
dev_intact_now      = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * devMOD_intact    + dev_steady;    
rot_intact_now      = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * rotMOD_intact    + rot_steady;    

stroke_damaged_now   = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * strokeMOD_damaged + stroke_steady;    
dev_damaged_now      = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * devMOD_damaged    + dev_steady;    
rot_damaged_now      = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * rotMOD_damaged    + rot_steady;    
