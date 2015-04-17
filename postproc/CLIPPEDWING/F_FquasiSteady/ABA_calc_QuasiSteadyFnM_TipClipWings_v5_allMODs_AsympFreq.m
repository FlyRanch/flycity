clear;
clc;
close all
warning off

load('bodyNwingModel_4qsModel.mat')

loadname=dir('WBdataset_steadyNclipMods_S2S3AmpRatioFunc.mat')
loadname = loadname.name;
load(loadname)

plot_on = 1
% plot_on = 0

rot_on=1;
% rot_on=0;

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

for freq_asymFitNr = 10


%% loop with different cuts
for i = 1:(sect_max-sect_min+1)

    N_cut = sect_max -i +1
    
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
    
    % asymptotic freq fit
    if freq_asymFitNr == 1
        freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit1,S2S3AmpRatioFunc_now);
    elseif freq_asymFitNr == 2
        freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit2,S2S3AmpRatioFunc_now);
    elseif freq_asymFitNr == 3
        freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit3,S2S3AmpRatioFunc_now);
    elseif freq_asymFitNr == 4
        freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit4,S2S3AmpRatioFunc_now);
    elseif freq_asymFitNr == 5
        freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit5,S2S3AmpRatioFunc_now);
    elseif freq_asymFitNr == 9
        freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit9,S2S3AmpRatioFunc_now);
    elseif freq_asymFitNr == 10
        freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit10,S2S3AmpRatioFunc_now);
    elseif freq_asymFitNr == 11
        freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit11,S2S3AmpRatioFunc_now);
    elseif freq_asymFitNr == 15
        freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit15,S2S3AmpRatioFunc_now);
    end        
        
        
%     freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_parabFit,S2S3AmpRatioFunc_now);
%     freq_now            = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * freqMOD + f_steady;    
    
    stroke_intact_now   = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * strokeMOD_intact + stroke_steady;    
    dev_intact_now      = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * devMOD_intact    + dev_steady;    
    rot_intact_now      = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * rotMOD_intact    + rot_steady;    

    stroke_damaged_now   = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * strokeMOD_damaged + stroke_steady;    
    dev_damaged_now      = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * devMOD_damaged    + dev_steady;    
    rot_damaged_now      = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * rotMOD_damaged    + rot_steady;    
    
    %% store data
    S2ratios(i,1) = S2ratio;
    S3ratios(i,1) = S3ratio;
    S2S3AmpRatioFuncs(i,1) = S2S3AmpRatioFunc_now;
    
    %% WB kin NO MODs
    freq = freq_steady;
%     freq = f_clipped_fly;
%     freq = freq_now;

    stroke_L = deg2rad(stroke_steady);
    dev_L = deg2rad(dev_steady);
    rot_L = deg2rad(rot_steady);

    stroke_R = deg2rad(stroke_steady);
    dev_R = deg2rad(dev_steady);
    rot_R = deg2rad(rot_steady);

    % qs forces & torques
    [ FM_strkpln, ~ ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_on );

    Fx_intact_steady = FM_strkpln.FM_L.Fx_strkpln_L;
    Fy_intact_steady = FM_strkpln.FM_L.Fy_strkpln_L;
    Fz_intact_steady = FM_strkpln.FM_L.Fz_strkpln_L;
    
    Mx_intact_steady = FM_strkpln.FM_L.Mx_strkpln_L;
    My_intact_steady = FM_strkpln.FM_L.My_strkpln_L;
    Mz_intact_steady = FM_strkpln.FM_L.Mz_strkpln_L;
    
    Fx_damaged_steady = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_damaged_steady = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_damaged_steady = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_damaged_steady = FM_strkpln.FM_R.Mx_strkpln_R;
    My_damaged_steady = FM_strkpln.FM_R.My_strkpln_R;
    Mz_damaged_steady = FM_strkpln.FM_R.Mz_strkpln_R;
    
    Fx_NONdamaged_steady = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_NONdamaged_steady = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_NONdamaged_steady = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_NONdamaged_steady = FM_strkpln.FM_R.Mx_strkpln_R;
    My_NONdamaged_steady = FM_strkpln.FM_R.My_strkpln_R;
    Mz_NONdamaged_steady = FM_strkpln.FM_R.Mz_strkpln_R;
    
    % remove forces and torques at removed part of wing (n>N_cut)
    Fx_damaged_steady(N_cut+1:end,:) = nan;
    Fy_damaged_steady(N_cut+1:end,:) = nan;
    Fz_damaged_steady(N_cut+1:end,:) = nan;
    
    Mx_damaged_steady(N_cut+1:end,:) = nan;
    My_damaged_steady(N_cut+1:end,:) = nan;
    Mz_damaged_steady(N_cut+1:end,:) = nan;

    % total force & torque
    Fx_damaged_total_steady(:,i) = (nansum(Fx_damaged_steady,1) + nansum(Fx_intact_steady,1))' / Mg_fly;
    Fy_damaged_total_steady(:,i) = (nansum(Fy_damaged_steady,1) + nansum(Fy_intact_steady,1))' / Mg_fly;
    Fz_damaged_total_steady(:,i) = (nansum(Fz_damaged_steady,1) + nansum(Fz_intact_steady,1))' / Mg_fly;
    
    Mx_damaged_total_steady(:,i) = (nansum(Mx_damaged_steady,1) + nansum(Mx_intact_steady,1))' / Mg_fly / (l_wing*1000);
    My_damaged_total_steady(:,i) = (nansum(My_damaged_steady,1) + nansum(My_intact_steady,1))' / Mg_fly / (l_wing*1000);
    Mz_damaged_total_steady(:,i) = (nansum(Mz_damaged_steady,1) + nansum(Mz_intact_steady,1))' / Mg_fly / (l_wing*1000);
    
    Fx_NONdamaged_total_steady(:,i) = (nansum(Fx_NONdamaged_steady,1) + nansum(Fx_intact_steady,1))' / Mg_fly;
    Fy_NONdamaged_total_steady(:,i) = (nansum(Fy_NONdamaged_steady,1) + nansum(Fy_intact_steady,1))' / Mg_fly;
    Fz_NONdamaged_total_steady(:,i) = (nansum(Fz_NONdamaged_steady,1) + nansum(Fz_intact_steady,1))' / Mg_fly;
    
    Mx_NONdamaged_total_steady(:,i) = (nansum(Mx_NONdamaged_steady,1) + nansum(Mx_intact_steady,1))' / Mg_fly / (l_wing*1000);
    My_NONdamaged_total_steady(:,i) = (nansum(My_NONdamaged_steady,1) + nansum(My_intact_steady,1))' / Mg_fly / (l_wing*1000);
    Mz_NONdamaged_total_steady(:,i) = (nansum(Mz_NONdamaged_steady,1) + nansum(Mz_intact_steady,1))' / Mg_fly / (l_wing*1000);
    
    %% WB kin ALL MODs
%     freq = freq_steady;
%     freq = f_clipped_fly;
    freq = freq_now;

    stroke_L = deg2rad(stroke_intact_now);
    dev_L = deg2rad(dev_intact_now);
    rot_L = deg2rad(rot_intact_now);

    stroke_R = deg2rad(stroke_damaged_now);
    dev_R = deg2rad(dev_damaged_now);
    rot_R = deg2rad(rot_damaged_now);

    % qs forces & torques
    [ FM_strkpln, ~ ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_on );

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
    
    Mx_damaged_total_all(:,i) = (nansum(Mx_damaged_all,1) + nansum(Mx_intact_all,1))' / Mg_fly / (l_wing*1000);
    My_damaged_total_all(:,i) = (nansum(My_damaged_all,1) + nansum(My_intact_all,1))' / Mg_fly / (l_wing*1000);
    Mz_damaged_total_all(:,i) = (nansum(Mz_damaged_all,1) + nansum(Mz_intact_all,1))' / Mg_fly / (l_wing*1000);
    
    Fx_NONdamaged_total_all(:,i) = (nansum(Fx_NONdamaged_all,1) + nansum(Fx_intact_all,1))' / Mg_fly;
    Fy_NONdamaged_total_all(:,i) = (nansum(Fy_NONdamaged_all,1) + nansum(Fy_intact_all,1))' / Mg_fly;
    Fz_NONdamaged_total_all(:,i) = (nansum(Fz_NONdamaged_all,1) + nansum(Fz_intact_all,1))' / Mg_fly;
    
    Mx_NONdamaged_total_all(:,i) = (nansum(Mx_NONdamaged_all,1) + nansum(Mx_intact_all,1))' / Mg_fly / (l_wing*1000);
    My_NONdamaged_total_all(:,i) = (nansum(My_NONdamaged_all,1) + nansum(My_intact_all,1))' / Mg_fly / (l_wing*1000);
    Mz_NONdamaged_total_all(:,i) = (nansum(Mz_NONdamaged_all,1) + nansum(Mz_intact_all,1))' / Mg_fly / (l_wing*1000);
    
    %% WB kin FREQ MODs
%     freq = freq_steady;
%     freq = f_clipped_fly;
    freq = freq_now;

    stroke_L = deg2rad(stroke_steady);
    dev_L = deg2rad(dev_steady);
    rot_L = deg2rad(rot_steady);

    stroke_R = deg2rad(stroke_steady);
    dev_R = deg2rad(dev_steady);
    rot_R = deg2rad(rot_steady);

    % qs forces & torques
    [ FM_strkpln, ~ ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_on );

    Fx_intact_freqMOD = FM_strkpln.FM_L.Fx_strkpln_L;
    Fy_intact_freqMOD = FM_strkpln.FM_L.Fy_strkpln_L;
    Fz_intact_freqMOD = FM_strkpln.FM_L.Fz_strkpln_L;
    
    Mx_intact_freqMOD = FM_strkpln.FM_L.Mx_strkpln_L;
    My_intact_freqMOD = FM_strkpln.FM_L.My_strkpln_L;
    Mz_intact_freqMOD = FM_strkpln.FM_L.Mz_strkpln_L;
    
    Fx_damaged_freqMOD = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_damaged_freqMOD = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_damaged_freqMOD = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_damaged_freqMOD = FM_strkpln.FM_R.Mx_strkpln_R;
    My_damaged_freqMOD = FM_strkpln.FM_R.My_strkpln_R;
    Mz_damaged_freqMOD = FM_strkpln.FM_R.Mz_strkpln_R;
    
    Fx_NONdamaged_freqMOD = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_NONdamaged_freqMOD = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_NONdamaged_freqMOD = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_NONdamaged_freqMOD = FM_strkpln.FM_R.Mx_strkpln_R;
    My_NONdamaged_freqMOD = FM_strkpln.FM_R.My_strkpln_R;
    Mz_NONdamaged_freqMOD = FM_strkpln.FM_R.Mz_strkpln_R;
    
    % remove forces and torques at removed part of wing (n>N_cut)
    Fx_damaged_freqMOD(N_cut+1:end,:) = nan;
    Fy_damaged_freqMOD(N_cut+1:end,:) = nan;
    Fz_damaged_freqMOD(N_cut+1:end,:) = nan;
    
    Mx_damaged_freqMOD(N_cut+1:end,:) = nan;
    My_damaged_freqMOD(N_cut+1:end,:) = nan;
    Mz_damaged_freqMOD(N_cut+1:end,:) = nan;

    % total force & torque
    Fx_damaged_total_freqMOD(:,i) = (nansum(Fx_damaged_freqMOD,1) + nansum(Fx_intact_freqMOD,1))' / Mg_fly;
    Fy_damaged_total_freqMOD(:,i) = (nansum(Fy_damaged_freqMOD,1) + nansum(Fy_intact_freqMOD,1))' / Mg_fly;
    Fz_damaged_total_freqMOD(:,i) = (nansum(Fz_damaged_freqMOD,1) + nansum(Fz_intact_freqMOD,1))' / Mg_fly;
    
    Mx_damaged_total_freqMOD(:,i) = (nansum(Mx_damaged_freqMOD,1) + nansum(Mx_intact_freqMOD,1))' / Mg_fly / (l_wing*1000);
    My_damaged_total_freqMOD(:,i) = (nansum(My_damaged_freqMOD,1) + nansum(My_intact_freqMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_damaged_total_freqMOD(:,i) = (nansum(Mz_damaged_freqMOD,1) + nansum(Mz_intact_freqMOD,1))' / Mg_fly / (l_wing*1000);
    
    Fx_NONdamaged_total_freqMOD(:,i) = (nansum(Fx_NONdamaged_freqMOD,1) + nansum(Fx_intact_freqMOD,1))' / Mg_fly;
    Fy_NONdamaged_total_freqMOD(:,i) = (nansum(Fy_NONdamaged_freqMOD,1) + nansum(Fy_intact_freqMOD,1))' / Mg_fly;
    Fz_NONdamaged_total_freqMOD(:,i) = (nansum(Fz_NONdamaged_freqMOD,1) + nansum(Fz_intact_freqMOD,1))' / Mg_fly;
    
    Mx_NONdamaged_total_freqMOD(:,i) = (nansum(Mx_NONdamaged_freqMOD,1) + nansum(Mx_intact_freqMOD,1))' / Mg_fly / (l_wing*1000);
    My_NONdamaged_total_freqMOD(:,i) = (nansum(My_NONdamaged_freqMOD,1) + nansum(My_intact_freqMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_NONdamaged_total_freqMOD(:,i) = (nansum(Mz_NONdamaged_freqMOD,1) + nansum(Mz_intact_freqMOD,1))' / Mg_fly / (l_wing*1000);
    
    %% WB kin stroke MODs
    freq = freq_steady;
%     freq = f_clipped_fly;
%     freq = freq_now;

    stroke_L = deg2rad(stroke_intact_now);
    dev_L = deg2rad(dev_steady);
    rot_L = deg2rad(rot_steady);

    stroke_R = deg2rad(stroke_damaged_now);
    dev_R = deg2rad(dev_steady);
    rot_R = deg2rad(rot_steady);

    % qs forces & torques
    [ FM_strkpln, ~ ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_on );

    Fx_intact_strokeMOD = FM_strkpln.FM_L.Fx_strkpln_L;
    Fy_intact_strokeMOD = FM_strkpln.FM_L.Fy_strkpln_L;
    Fz_intact_strokeMOD = FM_strkpln.FM_L.Fz_strkpln_L;
    
    Mx_intact_strokeMOD = FM_strkpln.FM_L.Mx_strkpln_L;
    My_intact_strokeMOD = FM_strkpln.FM_L.My_strkpln_L;
    Mz_intact_strokeMOD = FM_strkpln.FM_L.Mz_strkpln_L;
    
    Fx_damaged_strokeMOD = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_damaged_strokeMOD = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_damaged_strokeMOD = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_damaged_strokeMOD = FM_strkpln.FM_R.Mx_strkpln_R;
    My_damaged_strokeMOD = FM_strkpln.FM_R.My_strkpln_R;
    Mz_damaged_strokeMOD = FM_strkpln.FM_R.Mz_strkpln_R;
    
    Fx_NONdamaged_strokeMOD = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_NONdamaged_strokeMOD = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_NONdamaged_strokeMOD = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_NONdamaged_strokeMOD = FM_strkpln.FM_R.Mx_strkpln_R;
    My_NONdamaged_strokeMOD = FM_strkpln.FM_R.My_strkpln_R;
    Mz_NONdamaged_strokeMOD = FM_strkpln.FM_R.Mz_strkpln_R;
    
    % remove forces and torques at removed part of wing (n>N_cut)
    Fx_damaged_strokeMOD(N_cut+1:end,:) = nan;
    Fy_damaged_strokeMOD(N_cut+1:end,:) = nan;
    Fz_damaged_strokeMOD(N_cut+1:end,:) = nan;
    
    Mx_damaged_strokeMOD(N_cut+1:end,:) = nan;
    My_damaged_strokeMOD(N_cut+1:end,:) = nan;
    Mz_damaged_strokeMOD(N_cut+1:end,:) = nan;

    % total force & torque
    Fx_damaged_total_strokeMOD(:,i) = (nansum(Fx_damaged_strokeMOD,1) + nansum(Fx_intact_strokeMOD,1))' / Mg_fly;
    Fy_damaged_total_strokeMOD(:,i) = (nansum(Fy_damaged_strokeMOD,1) + nansum(Fy_intact_strokeMOD,1))' / Mg_fly;
    Fz_damaged_total_strokeMOD(:,i) = (nansum(Fz_damaged_strokeMOD,1) + nansum(Fz_intact_strokeMOD,1))' / Mg_fly;
    
    Mx_damaged_total_strokeMOD(:,i) = (nansum(Mx_damaged_strokeMOD,1) + nansum(Mx_intact_strokeMOD,1))' / Mg_fly / (l_wing*1000);
    My_damaged_total_strokeMOD(:,i) = (nansum(My_damaged_strokeMOD,1) + nansum(My_intact_strokeMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_damaged_total_strokeMOD(:,i) = (nansum(Mz_damaged_strokeMOD,1) + nansum(Mz_intact_strokeMOD,1))' / Mg_fly / (l_wing*1000);
    
    Fx_NONdamaged_total_strokeMOD(:,i) = (nansum(Fx_NONdamaged_strokeMOD,1) + nansum(Fx_intact_strokeMOD,1))' / Mg_fly;
    Fy_NONdamaged_total_strokeMOD(:,i) = (nansum(Fy_NONdamaged_strokeMOD,1) + nansum(Fy_intact_strokeMOD,1))' / Mg_fly;
    Fz_NONdamaged_total_strokeMOD(:,i) = (nansum(Fz_NONdamaged_strokeMOD,1) + nansum(Fz_intact_strokeMOD,1))' / Mg_fly;
    
    Mx_NONdamaged_total_strokeMOD(:,i) = (nansum(Mx_NONdamaged_strokeMOD,1) + nansum(Mx_intact_strokeMOD,1))' / Mg_fly / (l_wing*1000);
    My_NONdamaged_total_strokeMOD(:,i) = (nansum(My_NONdamaged_strokeMOD,1) + nansum(My_intact_strokeMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_NONdamaged_total_strokeMOD(:,i) = (nansum(Mz_NONdamaged_strokeMOD,1) + nansum(Mz_intact_strokeMOD,1))' / Mg_fly / (l_wing*1000);
    
    %% WB kin dev MOD
    freq = freq_steady;
%     freq = f_clipped_fly;
%     freq = freq_now;

    stroke_L = deg2rad(stroke_steady);
    dev_L = deg2rad(dev_intact_now);
    rot_L = deg2rad(rot_steady);

    stroke_R = deg2rad(stroke_steady);
    dev_R = deg2rad(dev_damaged_now);
    rot_R = deg2rad(rot_steady);

    % qs forces & torques
    [ FM_strkpln, ~ ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_on );

    Fx_intact_devMOD = FM_strkpln.FM_L.Fx_strkpln_L;
    Fy_intact_devMOD = FM_strkpln.FM_L.Fy_strkpln_L;
    Fz_intact_devMOD = FM_strkpln.FM_L.Fz_strkpln_L;
    
    Mx_intact_devMOD = FM_strkpln.FM_L.Mx_strkpln_L;
    My_intact_devMOD = FM_strkpln.FM_L.My_strkpln_L;
    Mz_intact_devMOD = FM_strkpln.FM_L.Mz_strkpln_L;
    
    Fx_damaged_devMOD = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_damaged_devMOD = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_damaged_devMOD = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_damaged_devMOD = FM_strkpln.FM_R.Mx_strkpln_R;
    My_damaged_devMOD = FM_strkpln.FM_R.My_strkpln_R;
    Mz_damaged_devMOD = FM_strkpln.FM_R.Mz_strkpln_R;
    
    Fx_NONdamaged_devMOD = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_NONdamaged_devMOD = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_NONdamaged_devMOD = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_NONdamaged_devMOD = FM_strkpln.FM_R.Mx_strkpln_R;
    My_NONdamaged_devMOD = FM_strkpln.FM_R.My_strkpln_R;
    Mz_NONdamaged_devMOD = FM_strkpln.FM_R.Mz_strkpln_R;
    
    % remove forces and torques at removed part of wing (n>N_cut)
    Fx_damaged_devMOD(N_cut+1:end,:) = nan;
    Fy_damaged_devMOD(N_cut+1:end,:) = nan;
    Fz_damaged_devMOD(N_cut+1:end,:) = nan;
    
    Mx_damaged_devMOD(N_cut+1:end,:) = nan;
    My_damaged_devMOD(N_cut+1:end,:) = nan;
    Mz_damaged_devMOD(N_cut+1:end,:) = nan;

    % total force & torque
    Fx_damaged_total_devMOD(:,i) = (nansum(Fx_damaged_devMOD,1) + nansum(Fx_intact_devMOD,1))' / Mg_fly;
    Fy_damaged_total_devMOD(:,i) = (nansum(Fy_damaged_devMOD,1) + nansum(Fy_intact_devMOD,1))' / Mg_fly;
    Fz_damaged_total_devMOD(:,i) = (nansum(Fz_damaged_devMOD,1) + nansum(Fz_intact_devMOD,1))' / Mg_fly;
    
    Mx_damaged_total_devMOD(:,i) = (nansum(Mx_damaged_devMOD,1) + nansum(Mx_intact_devMOD,1))' / Mg_fly / (l_wing*1000);
    My_damaged_total_devMOD(:,i) = (nansum(My_damaged_devMOD,1) + nansum(My_intact_devMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_damaged_total_devMOD(:,i) = (nansum(Mz_damaged_devMOD,1) + nansum(Mz_intact_devMOD,1))' / Mg_fly / (l_wing*1000);
    
    Fx_NONdamaged_total_devMOD(:,i) = (nansum(Fx_NONdamaged_devMOD,1) + nansum(Fx_intact_devMOD,1))' / Mg_fly;
    Fy_NONdamaged_total_devMOD(:,i) = (nansum(Fy_NONdamaged_devMOD,1) + nansum(Fy_intact_devMOD,1))' / Mg_fly;
    Fz_NONdamaged_total_devMOD(:,i) = (nansum(Fz_NONdamaged_devMOD,1) + nansum(Fz_intact_devMOD,1))' / Mg_fly;
    
    Mx_NONdamaged_total_devMOD(:,i) = (nansum(Mx_NONdamaged_devMOD,1) + nansum(Mx_intact_devMOD,1))' / Mg_fly / (l_wing*1000);
    My_NONdamaged_total_devMOD(:,i) = (nansum(My_NONdamaged_devMOD,1) + nansum(My_intact_devMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_NONdamaged_total_devMOD(:,i) = (nansum(Mz_NONdamaged_devMOD,1) + nansum(Mz_intact_devMOD,1))' / Mg_fly / (l_wing*1000);
    
    %% WB kin rot MOD
    freq = freq_steady;
%     freq = f_clipped_fly;
%     freq = freq_now;

    stroke_L = deg2rad(stroke_steady);
    dev_L = deg2rad(dev_steady);
    rot_L = deg2rad(rot_intact_now);

    stroke_R = deg2rad(stroke_steady);
    dev_R = deg2rad(dev_steady);
    rot_R = deg2rad(rot_damaged_now);

    % qs forces & torques
    [ FM_strkpln, ~ ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_on );

    Fx_intact_rotMOD = FM_strkpln.FM_L.Fx_strkpln_L;
    Fy_intact_rotMOD = FM_strkpln.FM_L.Fy_strkpln_L;
    Fz_intact_rotMOD = FM_strkpln.FM_L.Fz_strkpln_L;
    
    Mx_intact_rotMOD = FM_strkpln.FM_L.Mx_strkpln_L;
    My_intact_rotMOD = FM_strkpln.FM_L.My_strkpln_L;
    Mz_intact_rotMOD = FM_strkpln.FM_L.Mz_strkpln_L;
    
    Fx_damaged_rotMOD = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_damaged_rotMOD = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_damaged_rotMOD = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_damaged_rotMOD = FM_strkpln.FM_R.Mx_strkpln_R;
    My_damaged_rotMOD = FM_strkpln.FM_R.My_strkpln_R;
    Mz_damaged_rotMOD = FM_strkpln.FM_R.Mz_strkpln_R;
    
    Fx_NONdamaged_rotMOD = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_NONdamaged_rotMOD = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_NONdamaged_rotMOD = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_NONdamaged_rotMOD = FM_strkpln.FM_R.Mx_strkpln_R;
    My_NONdamaged_rotMOD = FM_strkpln.FM_R.My_strkpln_R;
    Mz_NONdamaged_rotMOD = FM_strkpln.FM_R.Mz_strkpln_R;
    
    % remove forces and torques at removed part of wing (n>N_cut)
    Fx_damaged_rotMOD(N_cut+1:end,:) = nan;
    Fy_damaged_rotMOD(N_cut+1:end,:) = nan;
    Fz_damaged_rotMOD(N_cut+1:end,:) = nan;
    
    Mx_damaged_rotMOD(N_cut+1:end,:) = nan;
    My_damaged_rotMOD(N_cut+1:end,:) = nan;
    Mz_damaged_rotMOD(N_cut+1:end,:) = nan;

    % total force & torque
    Fx_damaged_total_rotMOD(:,i) = (nansum(Fx_damaged_rotMOD,1) + nansum(Fx_intact_rotMOD,1))' / Mg_fly;
    Fy_damaged_total_rotMOD(:,i) = (nansum(Fy_damaged_rotMOD,1) + nansum(Fy_intact_rotMOD,1))' / Mg_fly;
    Fz_damaged_total_rotMOD(:,i) = (nansum(Fz_damaged_rotMOD,1) + nansum(Fz_intact_rotMOD,1))' / Mg_fly;
    
    Mx_damaged_total_rotMOD(:,i) = (nansum(Mx_damaged_rotMOD,1) + nansum(Mx_intact_rotMOD,1))' / Mg_fly / (l_wing*1000);
    My_damaged_total_rotMOD(:,i) = (nansum(My_damaged_rotMOD,1) + nansum(My_intact_rotMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_damaged_total_rotMOD(:,i) = (nansum(Mz_damaged_rotMOD,1) + nansum(Mz_intact_rotMOD,1))' / Mg_fly / (l_wing*1000);
    
    Fx_NONdamaged_total_rotMOD(:,i) = (nansum(Fx_NONdamaged_rotMOD,1) + nansum(Fx_intact_rotMOD,1))' / Mg_fly;
    Fy_NONdamaged_total_rotMOD(:,i) = (nansum(Fy_NONdamaged_rotMOD,1) + nansum(Fy_intact_rotMOD,1))' / Mg_fly;
    Fz_NONdamaged_total_rotMOD(:,i) = (nansum(Fz_NONdamaged_rotMOD,1) + nansum(Fz_intact_rotMOD,1))' / Mg_fly;
    
    Mx_NONdamaged_total_rotMOD(:,i) = (nansum(Mx_NONdamaged_rotMOD,1) + nansum(Mx_intact_rotMOD,1))' / Mg_fly / (l_wing*1000);
    My_NONdamaged_total_rotMOD(:,i) = (nansum(My_NONdamaged_rotMOD,1) + nansum(My_intact_rotMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_NONdamaged_total_rotMOD(:,i) = (nansum(Mz_NONdamaged_rotMOD,1) + nansum(Mz_intact_rotMOD,1))' / Mg_fly / (l_wing*1000);
    
end

    %% wingbeat average forces & torques: steady wb
    Fx_damaged_mean_steady = nanmean(Fx_damaged_total_steady);
    Fy_damaged_mean_steady = nanmean(Fy_damaged_total_steady);
    Fz_damaged_mean_steady = nanmean(Fz_damaged_total_steady);
    
    Mx_damaged_mean_steady = nanmean(Mx_damaged_total_steady);
    My_damaged_mean_steady = nanmean(My_damaged_total_steady);
    Mz_damaged_mean_steady = nanmean(Mz_damaged_total_steady);
    
    Fx_NONdamaged_mean_steady = nanmean(Fx_NONdamaged_total_steady);
    Fy_NONdamaged_mean_steady = nanmean(Fy_NONdamaged_total_steady);
    Fz_NONdamaged_mean_steady = nanmean(Fz_NONdamaged_total_steady);
    
    Mx_NONdamaged_mean_steady = nanmean(Mx_NONdamaged_total_steady);
    My_NONdamaged_mean_steady = nanmean(My_NONdamaged_total_steady);
    Mz_NONdamaged_mean_steady = nanmean(Mz_NONdamaged_total_steady);

    %% wingbeat average forces & torques: ALL MODs
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

    %% wingbeat average forces & torques: freq MOD
    Fx_damaged_mean_freqMOD = nanmean(Fx_damaged_total_freqMOD);
    Fy_damaged_mean_freqMOD = nanmean(Fy_damaged_total_freqMOD);
    Fz_damaged_mean_freqMOD = nanmean(Fz_damaged_total_freqMOD);
    
    Mx_damaged_mean_freqMOD = nanmean(Mx_damaged_total_freqMOD);
    My_damaged_mean_freqMOD = nanmean(My_damaged_total_freqMOD);
    Mz_damaged_mean_freqMOD = nanmean(Mz_damaged_total_freqMOD);
    
    Fx_NONdamaged_mean_freqMOD = nanmean(Fx_NONdamaged_total_freqMOD);
    Fy_NONdamaged_mean_freqMOD = nanmean(Fy_NONdamaged_total_freqMOD);
    Fz_NONdamaged_mean_freqMOD = nanmean(Fz_NONdamaged_total_freqMOD);
    
    Mx_NONdamaged_mean_freqMOD = nanmean(Mx_NONdamaged_total_freqMOD);
    My_NONdamaged_mean_freqMOD = nanmean(My_NONdamaged_total_freqMOD);
    Mz_NONdamaged_mean_freqMOD = nanmean(Mz_NONdamaged_total_freqMOD);

    %% wingbeat average forces & torques: stroke MOD
    Fx_damaged_mean_strokeMOD = nanmean(Fx_damaged_total_strokeMOD);
    Fy_damaged_mean_strokeMOD = nanmean(Fy_damaged_total_strokeMOD);
    Fz_damaged_mean_strokeMOD = nanmean(Fz_damaged_total_strokeMOD);
    
    Mx_damaged_mean_strokeMOD = nanmean(Mx_damaged_total_strokeMOD);
    My_damaged_mean_strokeMOD = nanmean(My_damaged_total_strokeMOD);
    Mz_damaged_mean_strokeMOD = nanmean(Mz_damaged_total_strokeMOD);
    
    Fx_NONdamaged_mean_strokeMOD = nanmean(Fx_NONdamaged_total_strokeMOD);
    Fy_NONdamaged_mean_strokeMOD = nanmean(Fy_NONdamaged_total_strokeMOD);
    Fz_NONdamaged_mean_strokeMOD = nanmean(Fz_NONdamaged_total_strokeMOD);
    
    Mx_NONdamaged_mean_strokeMOD = nanmean(Mx_NONdamaged_total_strokeMOD);
    My_NONdamaged_mean_strokeMOD = nanmean(My_NONdamaged_total_strokeMOD);
    Mz_NONdamaged_mean_strokeMOD = nanmean(Mz_NONdamaged_total_strokeMOD);

    %% wingbeat average forces & torques: dev MOD
    Fx_damaged_mean_devMOD = nanmean(Fx_damaged_total_devMOD);
    Fy_damaged_mean_devMOD = nanmean(Fy_damaged_total_devMOD);
    Fz_damaged_mean_devMOD = nanmean(Fz_damaged_total_devMOD);
    
    Mx_damaged_mean_devMOD = nanmean(Mx_damaged_total_devMOD);
    My_damaged_mean_devMOD = nanmean(My_damaged_total_devMOD);
    Mz_damaged_mean_devMOD = nanmean(Mz_damaged_total_devMOD);
    
    Fx_NONdamaged_mean_devMOD = nanmean(Fx_NONdamaged_total_devMOD);
    Fy_NONdamaged_mean_devMOD = nanmean(Fy_NONdamaged_total_devMOD);
    Fz_NONdamaged_mean_devMOD = nanmean(Fz_NONdamaged_total_devMOD);
    
    Mx_NONdamaged_mean_devMOD = nanmean(Mx_NONdamaged_total_devMOD);
    My_NONdamaged_mean_devMOD = nanmean(My_NONdamaged_total_devMOD);
    Mz_NONdamaged_mean_devMOD = nanmean(Mz_NONdamaged_total_devMOD);

    %% wingbeat average forces & torques: rot MOD
    Fx_damaged_mean_rotMOD = nanmean(Fx_damaged_total_rotMOD);
    Fy_damaged_mean_rotMOD = nanmean(Fy_damaged_total_rotMOD);
    Fz_damaged_mean_rotMOD = nanmean(Fz_damaged_total_rotMOD);
    
    Mx_damaged_mean_rotMOD = nanmean(Mx_damaged_total_rotMOD);
    My_damaged_mean_rotMOD = nanmean(My_damaged_total_rotMOD);
    Mz_damaged_mean_rotMOD = nanmean(Mz_damaged_total_rotMOD);
    
    Fx_NONdamaged_mean_rotMOD = nanmean(Fx_NONdamaged_total_rotMOD);
    Fy_NONdamaged_mean_rotMOD = nanmean(Fy_NONdamaged_total_rotMOD);
    Fz_NONdamaged_mean_rotMOD = nanmean(Fz_NONdamaged_total_rotMOD);
    
    Mx_NONdamaged_mean_rotMOD = nanmean(Mx_NONdamaged_total_rotMOD);
    My_NONdamaged_mean_rotMOD = nanmean(My_NONdamaged_total_rotMOD);
    Mz_NONdamaged_mean_rotMOD = nanmean(Mz_NONdamaged_total_rotMOD);

    save(['allMODs_TipClip_freqAsym',num2str(freq_asymFitNr),'.mat'])
    
    %% plot
    
    mkdir('qsModel_FnM_TipCut')
    cd('qsModel_FnM_TipCut')

    % ALL MODs VS NO MODs
    figure
    subplot(1,2,1)
    hold on
    plot(S2ratios,Fx_damaged_mean_all,'o-k','markersize',10,'markerfacecolor','b')
    plot(S2ratios,Fy_damaged_mean_all,'o-k','markersize',10,'markerfacecolor','r')
    plot(S2ratios,Fz_damaged_mean_all,'o-k','markersize',10,'markerfacecolor','g')
    
    plot(S2ratios,Fx_damaged_mean_steady,'o-k','markersize',10,'markerfacecolor','c')
    plot(S2ratios,Fy_damaged_mean_steady,'o-k','markersize',10,'markerfacecolor','m')
    plot(S2ratios,Fz_damaged_mean_steady,'o-k','markersize',10,'markerfacecolor','y')
    
    legend('x-axis Mod','y-axis Mod','z-axis Mod','x-axis steady','y-axis steady','z-axis steady','location','E')
    xlabel('S2 ratio')
    ylabel('normalized forces F/mg')
    axis([0.5 1 -1.5 .25])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-1.5:.25:1)
    
    subplot(1,2,2)
    hold on
    plot(S2ratios,Mx_damaged_mean_all,'o-k','markersize',10,'markerfacecolor','b')
    plot(S2ratios,My_damaged_mean_all-My_damaged_mean_all(S2ratios==1),'o-k','markersize',10,'markerfacecolor','r')
    plot(S2ratios,Mz_damaged_mean_all,'o-k','markersize',10,'markerfacecolor','g')
    
    plot(S2ratios,Mx_damaged_mean_steady,'o-k','markersize',10,'markerfacecolor','c')
    plot(S2ratios,My_damaged_mean_steady-My_damaged_mean_all(S2ratios==1),'o-k','markersize',10,'markerfacecolor','m')
    plot(S2ratios,Mz_damaged_mean_steady,'o-k','markersize',10,'markerfacecolor','y')
    
    xlabel('S3 ratio')
    ylabel('normalized torques T/mgl')
    axis([0.5 1 -.15 .2])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-.15:.05:.2)    
    
    saveas(gca,['FnM_WBmod_TipClip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['FnM_WBmod_TipClip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['FnM_WBmod_TipClip_asympFit',num2str(freq_asymFitNr),'.svg'])

    % compare MODs Fz & Mx
    
    figure
    subplot(1,2,1)
    hold on
    plot(S2ratios,Fz_damaged_mean_steady,'o-k','markersize',10,'markerfacecolor','w')
    plot(S2ratios,Fz_damaged_mean_all,'o-k','markersize',10,'markerfacecolor',[.5 .5 .5])
    plot(S2ratios,Fz_damaged_mean_freqMOD,'o-k','markersize',10,'markerfacecolor','b')
    plot(S2ratios,Fz_damaged_mean_strokeMOD,'o-k','markersize',10,'markerfacecolor','g')
    plot(S2ratios,Fz_damaged_mean_devMOD,'o-k','markersize',10,'markerfacecolor','r')
    plot(S2ratios,Fz_damaged_mean_rotMOD,'o-k','markersize',10,'markerfacecolor','c')
    
    legend('steady','modified','frequency','stroke','deviation','rotation','location','SE')
    xlabel('S2 ratio')
    ylabel('normalized vertical force Fz/mg')
    axis([0.5 1 -1.5 0])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-1.5:.25:1)
    
    subplot(1,2,2)
    hold on
    plot(S2ratios,Mx_damaged_mean_steady,'o-k','markersize',10,'markerfacecolor','w')
    plot(S2ratios,Mx_damaged_mean_all,'o-k','markersize',10,'markerfacecolor',[.5 .5 .5])
    plot(S2ratios,Mx_damaged_mean_freqMOD,'o-k','markersize',10,'markerfacecolor','b')
    plot(S2ratios,Mx_damaged_mean_strokeMOD,'o-k','markersize',10,'markerfacecolor','g')
    plot(S2ratios,Mx_damaged_mean_devMOD,'o-k','markersize',10,'markerfacecolor','r')
    plot(S2ratios,Mx_damaged_mean_rotMOD,'o-k','markersize',10,'markerfacecolor','c')
    
    xlabel('S3 ratio')
    ylabel('normalized roll torque Tx/mgl')
    axis([0.5 1 0 .3])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-.2:.1:.3)    

    saveas(gca,['Fz_Mx_WBmodComponents_TipClip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['Fz_Mx_WBmodComponents_TipClip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['Fz_Mx_WBmodComponents_TipClip_asympFit',num2str(freq_asymFitNr),'.svg'])
    
    % compare MODs: Fy & Mz
    
    figure
    subplot(1,2,1)
    hold on
    plot(S2ratios,Fy_damaged_mean_steady,'o-k','markersize',10,'markerfacecolor','w')
    plot(S2ratios,Fy_damaged_mean_all,'o-k','markersize',10,'markerfacecolor',[.5 .5 .5])
    plot(S2ratios,Fy_damaged_mean_freqMOD,'o-k','markersize',10,'markerfacecolor','b')
    plot(S2ratios,Fy_damaged_mean_strokeMOD,'o-k','markersize',10,'markerfacecolor','g')
    plot(S2ratios,Fy_damaged_mean_devMOD,'o-k','markersize',10,'markerfacecolor','r')
    plot(S2ratios,Fy_damaged_mean_rotMOD,'o-k','markersize',10,'markerfacecolor','c')
    
    legend('steady','modified','frequency','stroke','deviation','rotation','location','SE')
    xlabel('S2 ratio')
    ylabel('normalized sideways force Fy/mg')
    axis([0.5 1 -.5 0])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-1.5:.25:1)
    
    subplot(1,2,2)
    hold on
    plot(S2ratios,Mz_damaged_mean_steady,'o-k','markersize',10,'markerfacecolor','w')
    plot(S2ratios,Mz_damaged_mean_all,'o-k','markersize',10,'markerfacecolor',[.5 .5 .5])
    plot(S2ratios,Mz_damaged_mean_freqMOD,'o-k','markersize',10,'markerfacecolor','b')
    plot(S2ratios,Mz_damaged_mean_strokeMOD,'o-k','markersize',10,'markerfacecolor','g')
    plot(S2ratios,Mz_damaged_mean_devMOD,'o-k','markersize',10,'markerfacecolor','r')
    plot(S2ratios,Mz_damaged_mean_rotMOD,'o-k','markersize',10,'markerfacecolor','c')
    
    xlabel('S3 ratio')
    ylabel('normalized yaw torque Tz/mgl')
    axis([0.5 1 -.2 .1])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-.2:.1:.3)    

    saveas(gca,['Fy_Mz_WBmodComponents_TipClip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['Fy_Mz_WBmodComponents_TipClip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['Fy_Mz_WBmodComponents_TipClip_asympFit',num2str(freq_asymFitNr),'.svg'])
    
    cd ..
    end