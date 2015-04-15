clear;
clc;
% close all
warning off

load('bodyNwingModel_4qsModel.mat')

loadname=dir('WBdataset_steadyNclipMods_S2S3AmpRatioFunc*')
loadname = loadname.name;
load(loadname)

plot_on = 1
% plot_on = 0

rot_on=1;
rot_on=0;

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

[stroke_steady] = (calc_val_fourier_series_4thN8th_order(t_norm,stroke_coeffs_steady,0))';
[dev_steady] = (calc_val_fourier_series_4thN8th_order(t_norm,dev_coeffs_steady,0))';
[rot_steady] = (calc_val_fourier_series_4thN8th_order(t_norm,rot_coeffs_steady,0))';

%% damaged fly WB kin MODs
f_clipped_fly = f_wb_S2S3AmpRatioFunc_meanCIstd(1);
freqMOD = freqMOD_wb_S2S3AmpRatioFunc_meanCIstd(1);

% fouriers coeffs
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

%% S2 & S3 & S2S3AmpRatioFunc of intact wing
h_sect = wing_model.length/nr_sect;
y_sect = wing_model.y_sect_R(:,2);
c_sect = wing_model.chords_R';

S2_sect_element = 1/12 .* c_sect .* h_sect.^3;
S3_sect_element = 1/36 .* c_sect .* h_sect.^4;

S2_sect = c_sect .* h_sect .* y_sect.^2;
S3_sect = c_sect .* h_sect .* y_sect.^3;

S2_intact = sum(S2_sect);
S3_intact = sum(S3_sect);

sol = subs(solAdAiRatio,S2,1);
sol = subs(sol,S3,1);
S2S3AmpRatioFunc_NONclipped = eval(sol);

%% wing clipped & kin MODs
sect_min = round(3*nr_sect/4);
sect_max = nr_sect;

%% loop with different cuts
for i = 1:(sect_max-sect_min+1)

    N_cut = sect_max -i +1;
    
    % S2 & S3 of damaged wing
    y_sect = wing_model.y_sect_R(1:N_cut,2);
    c_sect = wing_model.chords_R(1:N_cut)';
    
    S2_sect = c_sect .* h_sect .* y_sect.^2;
    S3_sect = c_sect .* h_sect .* y_sect.^3;

    S2_damaged = sum(S2_sect);
    S3_damaged = sum(S3_sect);
    
    S2ratio = S2_damaged/S2_intact;
    S3ratio = S3_damaged/S3_intact;

    %% WB kin MOD based on stroke amplitude RATIO
    sol = subs(solAdAiRatio,S2,S2ratio);
    sol = subs(sol,S3,S3ratio);
    S2S3AmpRatioFunc_now = eval(sol);
    
    freq_now   = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * freqMOD + f_steady;    
    
    stroke_intact_now   = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * strokeMOD_intact + stroke_steady;    
    dev_intact_now      = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * devMOD_intact    + dev_steady;    
    rot_intact_now      = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * rotMOD_intact    + rot_steady;    

    stroke_damaged_now   = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * strokeMOD_damaged + stroke_steady;    
    dev_damaged_now      = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * devMOD_damaged    + dev_steady;    
    rot_damaged_now      = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * rotMOD_damaged    + rot_steady;    
    
    %% WB kin ALL MODs
%     freq = f_clipped_fly;
    freq = freq_now;

    stroke_L = deg2rad(stroke_intact_now);
    dev_L = deg2rad(dev_intact_now);
    rot_L = deg2rad(rot_intact_now);

    stroke_R = deg2rad(stroke_damaged_now);
    dev_R = deg2rad(dev_damaged_now);
    rot_R = deg2rad(rot_damaged_now);

    % qs forces & torques
    [ FM_strkpln, Vel_wingtip ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_on );

    Fx_intact_all = FM_strkpln.FM_L.Fx_strkpln_L;
    Fy_intact_all = FM_strkpln.FM_L.Fy_strkpln_L;
    Fz_intact_all = FM_strkpln.FM_L.Fz_strkpln_L;
    
    Mx_intact_all = FM_strkpln.FM_L.Mx_strkpln_L;
    My_intact_all = FM_strkpln.FM_L.My_strkpln_L;
    Mz_intact_all = FM_strkpln.FM_L.Mz_strkpln_L;
    
    Fx_damaged_all = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_damaged_all = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_damaged_all = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_damaged_all = FM_strkpln.FM_R.Mx_strkpln_R;
    My_damaged_all = FM_strkpln.FM_R.My_strkpln_R;
    Mz_damaged_all = FM_strkpln.FM_R.Mz_strkpln_R;
    
    Fx_NONdamaged_all = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_NONdamaged_all = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_NONdamaged_all = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_NONdamaged_all = FM_strkpln.FM_R.Mx_strkpln_R;
    My_NONdamaged_all = FM_strkpln.FM_R.My_strkpln_R;
    Mz_NONdamaged_all = FM_strkpln.FM_R.Mz_strkpln_R;
    
    % remove forces and torques at removed part of wing (n>N_cut)
    Fx_damaged_all(N_cut+1:end,:) = nan;
    Fy_damaged_all(N_cut+1:end,:) = nan;
    Fz_damaged_all(N_cut+1:end,:) = nan;
    
    Mx_damaged_all(N_cut+1:end,:) = nan;
    My_damaged_all(N_cut+1:end,:) = nan;
    Mz_damaged_all(N_cut+1:end,:) = nan;

    % total force & torque
    Fx_damaged_total_all(:,i) = (nansum(Fx_damaged_all,1) + nansum(Fx_intact_all,1))' / Mg_fly;
    Fy_damaged_total_all(:,i) = (nansum(Fy_damaged_all,1) + nansum(Fy_intact_all,1))' / Mg_fly;
    Fz_damaged_total_all(:,i) = (nansum(Fz_damaged_all,1) + nansum(Fz_intact_all,1))' / Mg_fly;
    
    Mx_damaged_total_all(:,i) = (nansum(Mx_damaged_all,1) + nansum(Mx_intact_all,1))' / Mg_fly / l_wing;
    My_damaged_total_all(:,i) = (nansum(My_damaged_all,1) + nansum(My_intact_all,1))' / Mg_fly / l_wing;
    Mz_damaged_total_all(:,i) = (nansum(Mz_damaged_all,1) + nansum(Mz_intact_all,1))' / Mg_fly / l_wing;
    
    Fx_NONdamaged_total_all(:,i) = (nansum(Fx_NONdamaged_all,1) + nansum(Fx_intact_all,1))' / Mg_fly;
    Fy_NONdamaged_total_all(:,i) = (nansum(Fy_NONdamaged_all,1) + nansum(Fy_intact_all,1))' / Mg_fly;
    Fz_NONdamaged_total_all(:,i) = (nansum(Fz_NONdamaged_all,1) + nansum(Fz_intact_all,1))' / Mg_fly;
    
    Mx_NONdamaged_total_all(:,i) = (nansum(Mx_NONdamaged_all,1) + nansum(Mx_intact_all,1))' / Mg_fly / l_wing;
    My_NONdamaged_total_all(:,i) = (nansum(My_NONdamaged_all,1) + nansum(My_intact_all,1))' / Mg_fly / l_wing;
    Mz_NONdamaged_total_all(:,i) = (nansum(Mz_NONdamaged_all,1) + nansum(Mz_intact_all,1))' / Mg_fly / l_wing;
    
    %% WB kin FREQ MODs
%     freq = f_clipped_fly;
    freq = freq_steady;

    stroke_L = deg2rad(stroke_steady);
    dev_L = deg2rad(dev_steady);
    rot_L = deg2rad(rot_steady);

    stroke_R = deg2rad(stroke_steady);
    dev_R = deg2rad(dev_steady);
    rot_R = deg2rad(rot_steady);

    % qs forces & torques
    [ FM_strkpln, Vel_wingtip ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_on );

    Fx_intact_all = FM_strkpln.FM_L.Fx_strkpln_L;
    Fy_intact_all = FM_strkpln.FM_L.Fy_strkpln_L;
    Fz_intact_all = FM_strkpln.FM_L.Fz_strkpln_L;
    
    Mx_intact_all = FM_strkpln.FM_L.Mx_strkpln_L;
    My_intact_all = FM_strkpln.FM_L.My_strkpln_L;
    Mz_intact_all = FM_strkpln.FM_L.Mz_strkpln_L;
    
    Fx_damaged_all = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_damaged_all = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_damaged_all = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_damaged_all = FM_strkpln.FM_R.Mx_strkpln_R;
    My_damaged_all = FM_strkpln.FM_R.My_strkpln_R;
    Mz_damaged_all = FM_strkpln.FM_R.Mz_strkpln_R;
    
    Fx_NONdamaged_all = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_NONdamaged_all = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_NONdamaged_all = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_NONdamaged_all = FM_strkpln.FM_R.Mx_strkpln_R;
    My_NONdamaged_all = FM_strkpln.FM_R.My_strkpln_R;
    Mz_NONdamaged_all = FM_strkpln.FM_R.Mz_strkpln_R;
    
    % remove forces and torques at removed part of wing (n>N_cut)
    Fx_damaged_all(N_cut+1:end,:) = nan;
    Fy_damaged_all(N_cut+1:end,:) = nan;
    Fz_damaged_all(N_cut+1:end,:) = nan;
    
    Mx_damaged_all(N_cut+1:end,:) = nan;
    My_damaged_all(N_cut+1:end,:) = nan;
    Mz_damaged_all(N_cut+1:end,:) = nan;

    % total force & torque
    Fx_damaged_total_all(:,i) = (nansum(Fx_damaged_all,1) + nansum(Fx_intact_all,1))' / Mg_fly;
    Fy_damaged_total_all(:,i) = (nansum(Fy_damaged_all,1) + nansum(Fy_intact_all,1))' / Mg_fly;
    Fz_damaged_total_all(:,i) = (nansum(Fz_damaged_all,1) + nansum(Fz_intact_all,1))' / Mg_fly;
    
    Mx_damaged_total_all(:,i) = (nansum(Mx_damaged_all,1) + nansum(Mx_intact_all,1))' / Mg_fly / l_wing;
    My_damaged_total_all(:,i) = (nansum(My_damaged_all,1) + nansum(My_intact_all,1))' / Mg_fly / l_wing;
    Mz_damaged_total_all(:,i) = (nansum(Mz_damaged_all,1) + nansum(Mz_intact_all,1))' / Mg_fly / l_wing;
    
    Fx_NONdamaged_total_all(:,i) = (nansum(Fx_NONdamaged_all,1) + nansum(Fx_intact_all,1))' / Mg_fly;
    Fy_NONdamaged_total_all(:,i) = (nansum(Fy_NONdamaged_all,1) + nansum(Fy_intact_all,1))' / Mg_fly;
    Fz_NONdamaged_total_all(:,i) = (nansum(Fz_NONdamaged_all,1) + nansum(Fz_intact_all,1))' / Mg_fly;
    
    Mx_NONdamaged_total_all(:,i) = (nansum(Mx_NONdamaged_all,1) + nansum(Mx_intact_all,1))' / Mg_fly / l_wing;
    My_NONdamaged_total_all(:,i) = (nansum(My_NONdamaged_all,1) + nansum(My_intact_all,1))' / Mg_fly / l_wing;
    Mz_NONdamaged_total_all(:,i) = (nansum(Mz_NONdamaged_all,1) + nansum(Mz_intact_all,1))' / Mg_fly / l_wing;
    
    %% WB kin ALL MODs
%     freq = f_clipped_fly;
    freq = freq_steady;

    stroke_L = deg2rad(stroke_intact_now);
    dev_L = deg2rad(dev_steady);
    rot_L = deg2rad(rot_steady);

    stroke_R = deg2rad(stroke_damaged_now);
    dev_R = deg2rad(dev_steady);
    rot_R = deg2rad(rot_steady);

    % qs forces & torques
    [ FM_strkpln, Vel_wingtip ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_on );

    Fx_intact_all = FM_strkpln.FM_L.Fx_strkpln_L;
    Fy_intact_all = FM_strkpln.FM_L.Fy_strkpln_L;
    Fz_intact_all = FM_strkpln.FM_L.Fz_strkpln_L;
    
    Mx_intact_all = FM_strkpln.FM_L.Mx_strkpln_L;
    My_intact_all = FM_strkpln.FM_L.My_strkpln_L;
    Mz_intact_all = FM_strkpln.FM_L.Mz_strkpln_L;
    
    Fx_damaged_all = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_damaged_all = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_damaged_all = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_damaged_all = FM_strkpln.FM_R.Mx_strkpln_R;
    My_damaged_all = FM_strkpln.FM_R.My_strkpln_R;
    Mz_damaged_all = FM_strkpln.FM_R.Mz_strkpln_R;
    
    Fx_NONdamaged_all = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_NONdamaged_all = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_NONdamaged_all = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_NONdamaged_all = FM_strkpln.FM_R.Mx_strkpln_R;
    My_NONdamaged_all = FM_strkpln.FM_R.My_strkpln_R;
    Mz_NONdamaged_all = FM_strkpln.FM_R.Mz_strkpln_R;
    
    % remove forces and torques at removed part of wing (n>N_cut)
    Fx_damaged_all(N_cut+1:end,:) = nan;
    Fy_damaged_all(N_cut+1:end,:) = nan;
    Fz_damaged_all(N_cut+1:end,:) = nan;
    
    Mx_damaged_all(N_cut+1:end,:) = nan;
    My_damaged_all(N_cut+1:end,:) = nan;
    Mz_damaged_all(N_cut+1:end,:) = nan;

    % total force & torque
    Fx_damaged_total_all(:,i) = (nansum(Fx_damaged_all,1) + nansum(Fx_intact_all,1))' / Mg_fly;
    Fy_damaged_total_all(:,i) = (nansum(Fy_damaged_all,1) + nansum(Fy_intact_all,1))' / Mg_fly;
    Fz_damaged_total_all(:,i) = (nansum(Fz_damaged_all,1) + nansum(Fz_intact_all,1))' / Mg_fly;
    
    Mx_damaged_total_all(:,i) = (nansum(Mx_damaged_all,1) + nansum(Mx_intact_all,1))' / Mg_fly / l_wing;
    My_damaged_total_all(:,i) = (nansum(My_damaged_all,1) + nansum(My_intact_all,1))' / Mg_fly / l_wing;
    Mz_damaged_total_all(:,i) = (nansum(Mz_damaged_all,1) + nansum(Mz_intact_all,1))' / Mg_fly / l_wing;
    
    Fx_NONdamaged_total_all(:,i) = (nansum(Fx_NONdamaged_all,1) + nansum(Fx_intact_all,1))' / Mg_fly;
    Fy_NONdamaged_total_all(:,i) = (nansum(Fy_NONdamaged_all,1) + nansum(Fy_intact_all,1))' / Mg_fly;
    Fz_NONdamaged_total_all(:,i) = (nansum(Fz_NONdamaged_all,1) + nansum(Fz_intact_all,1))' / Mg_fly;
    
    Mx_NONdamaged_total_all(:,i) = (nansum(Mx_NONdamaged_all,1) + nansum(Mx_intact_all,1))' / Mg_fly / l_wing;
    My_NONdamaged_total_all(:,i) = (nansum(My_NONdamaged_all,1) + nansum(My_intact_all,1))' / Mg_fly / l_wing;
    Mz_NONdamaged_total_all(:,i) = (nansum(Mz_NONdamaged_all,1) + nansum(Mz_intact_all,1))' / Mg_fly / l_wing;
    
    %% WB kin ALL MODs
%     freq = f_clipped_fly;
    freq = freq_steady;

    stroke_L = deg2rad(stroke_steady);
    dev_L = deg2rad(dev_intact_now);
    rot_L = deg2rad(rot_steady);

    stroke_R = deg2rad(stroke_steady);
    dev_R = deg2rad(dev_damaged_now);
    rot_R = deg2rad(rot_steady);

    % qs forces & torques
    [ FM_strkpln, Vel_wingtip ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_on );

    Fx_intact_all = FM_strkpln.FM_L.Fx_strkpln_L;
    Fy_intact_all = FM_strkpln.FM_L.Fy_strkpln_L;
    Fz_intact_all = FM_strkpln.FM_L.Fz_strkpln_L;
    
    Mx_intact_all = FM_strkpln.FM_L.Mx_strkpln_L;
    My_intact_all = FM_strkpln.FM_L.My_strkpln_L;
    Mz_intact_all = FM_strkpln.FM_L.Mz_strkpln_L;
    
    Fx_damaged_all = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_damaged_all = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_damaged_all = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_damaged_all = FM_strkpln.FM_R.Mx_strkpln_R;
    My_damaged_all = FM_strkpln.FM_R.My_strkpln_R;
    Mz_damaged_all = FM_strkpln.FM_R.Mz_strkpln_R;
    
    Fx_NONdamaged_all = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_NONdamaged_all = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_NONdamaged_all = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_NONdamaged_all = FM_strkpln.FM_R.Mx_strkpln_R;
    My_NONdamaged_all = FM_strkpln.FM_R.My_strkpln_R;
    Mz_NONdamaged_all = FM_strkpln.FM_R.Mz_strkpln_R;
    
    % remove forces and torques at removed part of wing (n>N_cut)
    Fx_damaged_all(N_cut+1:end,:) = nan;
    Fy_damaged_all(N_cut+1:end,:) = nan;
    Fz_damaged_all(N_cut+1:end,:) = nan;
    
    Mx_damaged_all(N_cut+1:end,:) = nan;
    My_damaged_all(N_cut+1:end,:) = nan;
    Mz_damaged_all(N_cut+1:end,:) = nan;

    % total force & torque
    Fx_damaged_total_all(:,i) = (nansum(Fx_damaged_all,1) + nansum(Fx_intact_all,1))' / Mg_fly;
    Fy_damaged_total_all(:,i) = (nansum(Fy_damaged_all,1) + nansum(Fy_intact_all,1))' / Mg_fly;
    Fz_damaged_total_all(:,i) = (nansum(Fz_damaged_all,1) + nansum(Fz_intact_all,1))' / Mg_fly;
    
    Mx_damaged_total_all(:,i) = (nansum(Mx_damaged_all,1) + nansum(Mx_intact_all,1))' / Mg_fly / l_wing;
    My_damaged_total_all(:,i) = (nansum(My_damaged_all,1) + nansum(My_intact_all,1))' / Mg_fly / l_wing;
    Mz_damaged_total_all(:,i) = (nansum(Mz_damaged_all,1) + nansum(Mz_intact_all,1))' / Mg_fly / l_wing;
    
    Fx_NONdamaged_total_all(:,i) = (nansum(Fx_NONdamaged_all,1) + nansum(Fx_intact_all,1))' / Mg_fly;
    Fy_NONdamaged_total_all(:,i) = (nansum(Fy_NONdamaged_all,1) + nansum(Fy_intact_all,1))' / Mg_fly;
    Fz_NONdamaged_total_all(:,i) = (nansum(Fz_NONdamaged_all,1) + nansum(Fz_intact_all,1))' / Mg_fly;
    
    Mx_NONdamaged_total_all(:,i) = (nansum(Mx_NONdamaged_all,1) + nansum(Mx_intact_all,1))' / Mg_fly / l_wing;
    My_NONdamaged_total_all(:,i) = (nansum(My_NONdamaged_all,1) + nansum(My_intact_all,1))' / Mg_fly / l_wing;
    Mz_NONdamaged_total_all(:,i) = (nansum(Mz_NONdamaged_all,1) + nansum(Mz_intact_all,1))' / Mg_fly / l_wing;
    
    %% WB kin ALL MODs
%     freq = f_clipped_fly;
    freq = freq_steady;

    stroke_L = deg2rad(stroke_steady);
    dev_L = deg2rad(dev_steady);
    rot_L = deg2rad(rot_intact_now);

    stroke_R = deg2rad(stroke_steady);
    dev_R = deg2rad(dev_steady);
    rot_R = deg2rad(rot_damaged_now);

    % qs forces & torques
    [ FM_strkpln, Vel_wingtip ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_on );

    Fx_intact_all = FM_strkpln.FM_L.Fx_strkpln_L;
    Fy_intact_all = FM_strkpln.FM_L.Fy_strkpln_L;
    Fz_intact_all = FM_strkpln.FM_L.Fz_strkpln_L;
    
    Mx_intact_all = FM_strkpln.FM_L.Mx_strkpln_L;
    My_intact_all = FM_strkpln.FM_L.My_strkpln_L;
    Mz_intact_all = FM_strkpln.FM_L.Mz_strkpln_L;
    
    Fx_damaged_all = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_damaged_all = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_damaged_all = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_damaged_all = FM_strkpln.FM_R.Mx_strkpln_R;
    My_damaged_all = FM_strkpln.FM_R.My_strkpln_R;
    Mz_damaged_all = FM_strkpln.FM_R.Mz_strkpln_R;
    
    Fx_NONdamaged_all = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_NONdamaged_all = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_NONdamaged_all = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_NONdamaged_all = FM_strkpln.FM_R.Mx_strkpln_R;
    My_NONdamaged_all = FM_strkpln.FM_R.My_strkpln_R;
    Mz_NONdamaged_all = FM_strkpln.FM_R.Mz_strkpln_R;
    
    % remove forces and torques at removed part of wing (n>N_cut)
    Fx_damaged_all(N_cut+1:end,:) = nan;
    Fy_damaged_all(N_cut+1:end,:) = nan;
    Fz_damaged_all(N_cut+1:end,:) = nan;
    
    Mx_damaged_all(N_cut+1:end,:) = nan;
    My_damaged_all(N_cut+1:end,:) = nan;
    Mz_damaged_all(N_cut+1:end,:) = nan;

    % total force & torque
    Fx_damaged_total_all(:,i) = (nansum(Fx_damaged_all,1) + nansum(Fx_intact_all,1))' / Mg_fly;
    Fy_damaged_total_all(:,i) = (nansum(Fy_damaged_all,1) + nansum(Fy_intact_all,1))' / Mg_fly;
    Fz_damaged_total_all(:,i) = (nansum(Fz_damaged_all,1) + nansum(Fz_intact_all,1))' / Mg_fly;
    
    Mx_damaged_total_all(:,i) = (nansum(Mx_damaged_all,1) + nansum(Mx_intact_all,1))' / Mg_fly / l_wing;
    My_damaged_total_all(:,i) = (nansum(My_damaged_all,1) + nansum(My_intact_all,1))' / Mg_fly / l_wing;
    Mz_damaged_total_all(:,i) = (nansum(Mz_damaged_all,1) + nansum(Mz_intact_all,1))' / Mg_fly / l_wing;
    
    Fx_NONdamaged_total_all(:,i) = (nansum(Fx_NONdamaged_all,1) + nansum(Fx_intact_all,1))' / Mg_fly;
    Fy_NONdamaged_total_all(:,i) = (nansum(Fy_NONdamaged_all,1) + nansum(Fy_intact_all,1))' / Mg_fly;
    Fz_NONdamaged_total_all(:,i) = (nansum(Fz_NONdamaged_all,1) + nansum(Fz_intact_all,1))' / Mg_fly;
    
    Mx_NONdamaged_total_all(:,i) = (nansum(Mx_NONdamaged_all,1) + nansum(Mx_intact_all,1))' / Mg_fly / l_wing;
    My_NONdamaged_total_all(:,i) = (nansum(My_NONdamaged_all,1) + nansum(My_intact_all,1))' / Mg_fly / l_wing;
    Mz_NONdamaged_total_all(:,i) = (nansum(Mz_NONdamaged_all,1) + nansum(Mz_intact_all,1))' / Mg_fly / l_wing;
    
    %% store data
    S2ratios(i,1) = S2ratio;
    S3ratios(i,1) = S3ratio;
    S2S3AmpRatioFuncs(i,1) = S2S3AmpRatioFunc_now;
end


    %% wingbeat average forces & torques
    Fx_damaged_mean_all = nanmean(Fx_damaged_total_all);
    Fy_damaged_mean_all = nanmean(Fy_damaged_total_all);
    Fz_damaged_mean_all = nanmean(Fz_damaged_total_all);
    
    Mx_damaged_mean_all = nanmean(Mx_damaged_total_all);
    My_damaged_mean_all = nanmean(My_damaged_total_all);
    Mz_damaged_mean_all = nanmean(Mz_damaged_total_all);
    
    Fx_NONdamaged_mean_all = nanmean(Fx_NONdamaged_total_all);
    Fy_NONdamaged_mean_all = nanmean(Fy_NONdamaged_total_all);
    Fz_NONdamaged_mean_all = nanmean(Fz_NONdamaged_total_all);
    
    Mx_NONdamaged_mean_all = nanmean(Mx_NONdamaged_total_all);
    My_NONdamaged_mean_all = nanmean(My_NONdamaged_total_all);
    Mz_NONdamaged_mean_all = nanmean(Mz_NONdamaged_total_all);

    figure
    subplot(2,1,1)
    hold on
    plot(S2ratios,Fx_damaged_mean_all,'o-b')
    plot(S2ratios,Fy_damaged_mean_all,'o-g')
    plot(S2ratios,Fz_damaged_mean_all,'o-r')
    
    plot(S2ratios,Fx_NONdamaged_mean_all,'*-b')
    plot(S2ratios,Fy_NONdamaged_mean_all,'*-g')
    plot(S2ratios,Fz_NONdamaged_mean_all,'*-r')
    
    subplot(2,1,2)
    hold on
    plot(S2ratios,Mx_damaged_mean_all,'o-b')
    plot(S2ratios,My_damaged_mean_all,'o-g')
    plot(S2ratios,Mz_damaged_mean_all,'o-r')
    
    plot(S2ratios,Mx_NONdamaged_mean_all,'*-b')
    plot(S2ratios,My_NONdamaged_mean_all,'*-g')
    plot(S2ratios,Mz_NONdamaged_mean_all,'*-r')
    
    
    

