clear;
clc;
close all
warning off

load('bodyNwingModel_4qsModel.mat')

loadname=dir('WBdataset_steadyNclipMods_S2S3AmpRatioFunc.mat')
loadname = loadname.name;
load(loadname)

load('MODsNstats_bodyNfreq_indiv.mat')

plot_on = 1
% plot_on = 0

rot_on=1;
% rot_on=0;

%% constants
% freq & roll fit type
% fit_type = 1   % asymp10
% fit_type = 2   % spline999;
% fit_type = 3   % linear;
% fit_type = 4   % power2;

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
sect_min = 16;
sect_max = nr_sect;

for fit_type = 2
% for fit_type = 1:4

fit_type
% fit_type = 1   % asymp10
% fit_type = 2   % spline999;
% fit_type = 3   % linear;
% fit_type = 4   % power2;

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
    
%     % asymptotic freq fit
%     if freq_asymFitNr == 1
%         freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit1,S2S3AmpRatioFunc_now);
%     elseif freq_asymFitNr == 2
%         freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit2,S2S3AmpRatioFunc_now);
%     elseif freq_asymFitNr == 3
%         freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit3,S2S3AmpRatioFunc_now);
%     elseif freq_asymFitNr == 4
%         freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit4,S2S3AmpRatioFunc_now);
%     elseif freq_asymFitNr == 5
%         freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit5,S2S3AmpRatioFunc_now);
%     elseif freq_asymFitNr == 9
%         freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit9,S2S3AmpRatioFunc_now);
%     elseif freq_asymFitNr == 10
%         freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit10,S2S3AmpRatioFunc_now);
%     elseif freq_asymFitNr == 11
%         freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit11,S2S3AmpRatioFunc_now);
%     elseif freq_asymFitNr == 15
%         freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_asympFit15,S2S3AmpRatioFunc_now);
%     end        
%         
%         
%     freq_now            = feval(freq_S2S3AmpRatioFunc_steadyWBs_parabFit,S2S3AmpRatioFunc_now);
%     freq_now            = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * freqMOD + f_steady;    

    stroke_intact_now   = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * strokeMOD_intact + stroke_steady;    
    dev_intact_now      = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * devMOD_intact    + dev_steady;    
    rot_intact_now      = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * rotMOD_intact    + rot_steady;    

    stroke_damaged_now   = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * strokeMOD_damaged + stroke_steady;    
    dev_damaged_now      = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * devMOD_damaged    + dev_steady;    
    rot_damaged_now      = (S2S3AmpRatioFunc_now-S2S3AmpRatioFunc_NONclipped) * rotMOD_damaged    + rot_steady;    

    %% indiv-based freq & roll
    % asympFit10
    freq_asymp10           = feval(freq_S2S3AmpRatioFunc_indiv_asympFit10,S2S3AmpRatioFunc_now);
    roll_asymp10           = feval(roll_S2S3AmpRatioFunc_indiv_asympFit10,S2S3AmpRatioFunc_now);

    % spline999
    freq_spline999         = feval(freq_S2S3AmpRatioFunc_indiv_smooth999,S2S3AmpRatioFunc_now);
    roll_spline999         = feval(roll_S2S3AmpRatioFunc_indiv_smooth999,S2S3AmpRatioFunc_now);

    % linear
    freq_linear            = feval(freq_S2S3AmpRatioFunc_indiv_linear,S2S3AmpRatioFunc_now);
    roll_linear            = feval(roll_S2S3AmpRatioFunc_indiv_linear,S2S3AmpRatioFunc_now);

    % power2
    freq_power2            = feval(freq_S2S3AmpRatioFunc_indiv_power2,S2S3AmpRatioFunc_now);
    roll_power2            = feval(roll_S2S3AmpRatioFunc_indiv_power2,S2S3AmpRatioFunc_now);
    
    if fit_type == 1
        freq_now = freq_asymp10;
        roll_now = roll_asymp10;
    elseif fit_type == 2
        freq_now = freq_spline999;
        roll_now = roll_spline999;
    elseif fit_type == 3
        freq_now = freq_linear;
        roll_now = roll_linear;
    elseif fit_type == 4
        freq_now = freq_power2;
        roll_now = roll_power2;
    end
    
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

    % total force & torque damaged & intact wing
    Fx_damagedwing_steady(:,i) = (nansum(Fx_damaged_steady,1))' / Mg_fly;
    Fy_damagedwing_steady(:,i) = (nansum(Fy_damaged_steady,1))' / Mg_fly;
    Fz_damagedwing_steady(:,i) = (nansum(Fz_damaged_steady,1))' / Mg_fly;
    
    Mx_damagedwing_steady(:,i) = (nansum(Mx_damaged_steady,1))' / Mg_fly / (l_wing*1000);
    My_damagedwing_steady(:,i) = (nansum(My_damaged_steady,1))' / Mg_fly / (l_wing*1000);
    Mz_damagedwing_steady(:,i) = (nansum(Mz_damaged_steady,1))' / Mg_fly / (l_wing*1000);
    
    Fx_intactwing_steady(:,i) = (nansum(Fx_intact_steady,1))' / Mg_fly;
    Fy_intactwing_steady(:,i) = (nansum(Fy_intact_steady,1))' / Mg_fly;
    Fz_intactwing_steady(:,i) = (nansum(Fz_intact_steady,1))' / Mg_fly;
    
    Mx_intactwing_steady(:,i) = (nansum(Mx_intact_steady,1))' / Mg_fly / (l_wing*1000);
    My_intactwing_steady(:,i) = (nansum(My_intact_steady,1))' / Mg_fly / (l_wing*1000);
    Mz_intactwing_steady(:,i) = (nansum(Mz_intact_steady,1))' / Mg_fly / (l_wing*1000);

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
    
    %% remove forces and torques at removed part of wing (n>N_cut)
    Fx_damaged_all(N_cut+1:end,:) = nan;
    Fy_damaged_all(N_cut+1:end,:) = nan;
    Fz_damaged_all(N_cut+1:end,:) = nan;
    
    Mx_damaged_all(N_cut+1:end,:) = nan;
    My_damaged_all(N_cut+1:end,:) = nan;
    Mz_damaged_all(N_cut+1:end,:) = nan;

    %% store data
    % wingbeat kinematics
    roll_all(i,1) = roll_now;
    freq_all(i,1) = freq_now;

    stroke_intact_all(:,i) = (stroke_intact_now);
    dev_intact_all(:,i) = (dev_intact_now);
    rot_intact_all(:,i) = (rot_intact_now);

    stroke_damaged_all(:,i) = (stroke_damaged_now);
    dev_damaged_all(:,i) = (dev_damaged_now);
    rot_damaged_all(:,i) = (rot_damaged_now);
    
    % total force & torque damaged & intact wing
    Fx_damagedwing_all(:,i) = (nansum(Fx_damaged_all,1))' / Mg_fly;
    Fy_damagedwing_all(:,i) = (nansum(Fy_damaged_all,1))' / Mg_fly;
    Fz_damagedwing_all(:,i) = (nansum(Fz_damaged_all,1))' / Mg_fly;
    
    Mx_damagedwing_all(:,i) = (nansum(Mx_damaged_all,1))' / Mg_fly / (l_wing*1000);
    My_damagedwing_all(:,i) = (nansum(My_damaged_all,1))' / Mg_fly / (l_wing*1000);
    Mz_damagedwing_all(:,i) = (nansum(Mz_damaged_all,1))' / Mg_fly / (l_wing*1000);
    
    Fx_intactwing_all(:,i) = (nansum(Fx_intact_all,1))' / Mg_fly;
    Fy_intactwing_all(:,i) = (nansum(Fy_intact_all,1))' / Mg_fly;
    Fz_intactwing_all(:,i) = (nansum(Fz_intact_all,1))' / Mg_fly;
    
    Mx_intactwing_all(:,i) = (nansum(Mx_intact_all,1))' / Mg_fly / (l_wing*1000);
    My_intactwing_all(:,i) = (nansum(My_intact_all,1))' / Mg_fly / (l_wing*1000);
    Mz_intactwing_all(:,i) = (nansum(Mz_intact_all,1))' / Mg_fly / (l_wing*1000);

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

    % total force & torque damaged & intact wing
    Fx_damagedwing_freqMOD(:,i) = (nansum(Fx_damaged_freqMOD,1))' / Mg_fly;
    Fy_damagedwing_freqMOD(:,i) = (nansum(Fy_damaged_freqMOD,1))' / Mg_fly;
    Fz_damagedwing_freqMOD(:,i) = (nansum(Fz_damaged_freqMOD,1))' / Mg_fly;
    
    Mx_damagedwing_freqMOD(:,i) = (nansum(Mx_damaged_freqMOD,1))' / Mg_fly / (l_wing*1000);
    My_damagedwing_freqMOD(:,i) = (nansum(My_damaged_freqMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_damagedwing_freqMOD(:,i) = (nansum(Mz_damaged_freqMOD,1))' / Mg_fly / (l_wing*1000);
    
    Fx_intactwing_freqMOD(:,i) = (nansum(Fx_intact_freqMOD,1))' / Mg_fly;
    Fy_intactwing_freqMOD(:,i) = (nansum(Fy_intact_freqMOD,1))' / Mg_fly;
    Fz_intactwing_freqMOD(:,i) = (nansum(Fz_intact_freqMOD,1))' / Mg_fly;
    
    Mx_intactwing_freqMOD(:,i) = (nansum(Mx_intact_freqMOD,1))' / Mg_fly / (l_wing*1000);
    My_intactwing_freqMOD(:,i) = (nansum(My_intact_freqMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_intactwing_freqMOD(:,i) = (nansum(Mz_intact_freqMOD,1))' / Mg_fly / (l_wing*1000);

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

    % total force & torque damaged & intact wing
    Fx_damagedwing_strokeMOD(:,i) = (nansum(Fx_damaged_strokeMOD,1))' / Mg_fly;
    Fy_damagedwing_strokeMOD(:,i) = (nansum(Fy_damaged_strokeMOD,1))' / Mg_fly;
    Fz_damagedwing_strokeMOD(:,i) = (nansum(Fz_damaged_strokeMOD,1))' / Mg_fly;
    
    Mx_damagedwing_strokeMOD(:,i) = (nansum(Mx_damaged_strokeMOD,1))' / Mg_fly / (l_wing*1000);
    My_damagedwing_strokeMOD(:,i) = (nansum(My_damaged_strokeMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_damagedwing_strokeMOD(:,i) = (nansum(Mz_damaged_strokeMOD,1))' / Mg_fly / (l_wing*1000);
    
    Fx_intactwing_strokeMOD(:,i) = (nansum(Fx_intact_strokeMOD,1))' / Mg_fly;
    Fy_intactwing_strokeMOD(:,i) = (nansum(Fy_intact_strokeMOD,1))' / Mg_fly;
    Fz_intactwing_strokeMOD(:,i) = (nansum(Fz_intact_strokeMOD,1))' / Mg_fly;
    
    Mx_intactwing_strokeMOD(:,i) = (nansum(Mx_intact_strokeMOD,1))' / Mg_fly / (l_wing*1000);
    My_intactwing_strokeMOD(:,i) = (nansum(My_intact_strokeMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_intactwing_strokeMOD(:,i) = (nansum(Mz_intact_strokeMOD,1))' / Mg_fly / (l_wing*1000);

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

    % total force & torque damaged & intact wing
    Fx_damagedwing_devMOD(:,i) = (nansum(Fx_damaged_devMOD,1))' / Mg_fly;
    Fy_damagedwing_devMOD(:,i) = (nansum(Fy_damaged_devMOD,1))' / Mg_fly;
    Fz_damagedwing_devMOD(:,i) = (nansum(Fz_damaged_devMOD,1))' / Mg_fly;
    
    Mx_damagedwing_devMOD(:,i) = (nansum(Mx_damaged_devMOD,1))' / Mg_fly / (l_wing*1000);
    My_damagedwing_devMOD(:,i) = (nansum(My_damaged_devMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_damagedwing_devMOD(:,i) = (nansum(Mz_damaged_devMOD,1))' / Mg_fly / (l_wing*1000);
    
    Fx_intactwing_devMOD(:,i) = (nansum(Fx_intact_devMOD,1))' / Mg_fly;
    Fy_intactwing_devMOD(:,i) = (nansum(Fy_intact_devMOD,1))' / Mg_fly;
    Fz_intactwing_devMOD(:,i) = (nansum(Fz_intact_devMOD,1))' / Mg_fly;
    
    Mx_intactwing_devMOD(:,i) = (nansum(Mx_intact_devMOD,1))' / Mg_fly / (l_wing*1000);
    My_intactwing_devMOD(:,i) = (nansum(My_intact_devMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_intactwing_devMOD(:,i) = (nansum(Mz_intact_devMOD,1))' / Mg_fly / (l_wing*1000);

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

    % total force & torque damaged & intact wing
    Fx_damagedwing_rotMOD(:,i) = (nansum(Fx_damaged_rotMOD,1))' / Mg_fly;
    Fy_damagedwing_rotMOD(:,i) = (nansum(Fy_damaged_rotMOD,1))' / Mg_fly;
    Fz_damagedwing_rotMOD(:,i) = (nansum(Fz_damaged_rotMOD,1))' / Mg_fly;
    
    Mx_damagedwing_rotMOD(:,i) = (nansum(Mx_damaged_rotMOD,1))' / Mg_fly / (l_wing*1000);
    My_damagedwing_rotMOD(:,i) = (nansum(My_damaged_rotMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_damagedwing_rotMOD(:,i) = (nansum(Mz_damaged_rotMOD,1))' / Mg_fly / (l_wing*1000);
    
    Fx_intactwing_rotMOD(:,i) = (nansum(Fx_intact_rotMOD,1))' / Mg_fly;
    Fy_intactwing_rotMOD(:,i) = (nansum(Fy_intact_rotMOD,1))' / Mg_fly;
    Fz_intactwing_rotMOD(:,i) = (nansum(Fz_intact_rotMOD,1))' / Mg_fly;
    
    Mx_intactwing_rotMOD(:,i) = (nansum(Mx_intact_rotMOD,1))' / Mg_fly / (l_wing*1000);
    My_intactwing_rotMOD(:,i) = (nansum(My_intact_rotMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_intactwing_rotMOD(:,i) = (nansum(Mz_intact_rotMOD,1))' / Mg_fly / (l_wing*1000);

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
    
    %% WB kin freq & Rot MODs
%     freq = freq_steady;
%     freq = f_clipped_fly;
    freq = freq_now;

    stroke_L = deg2rad(stroke_steady);
    dev_L = deg2rad(dev_steady);
    rot_L = deg2rad(rot_intact_now);

    stroke_R = deg2rad(stroke_steady);
    dev_R = deg2rad(dev_steady);
    rot_R = deg2rad(rot_damaged_now);

    % qs forces & torques
    [ FM_strkpln, ~ ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_on );

    Fx_intact_freqRotMOD = FM_strkpln.FM_L.Fx_strkpln_L;
    Fy_intact_freqRotMOD = FM_strkpln.FM_L.Fy_strkpln_L;
    Fz_intact_freqRotMOD = FM_strkpln.FM_L.Fz_strkpln_L;
    
    Mx_intact_freqRotMOD = FM_strkpln.FM_L.Mx_strkpln_L;
    My_intact_freqRotMOD = FM_strkpln.FM_L.My_strkpln_L;
    Mz_intact_freqRotMOD = FM_strkpln.FM_L.Mz_strkpln_L;
    
    Fx_damaged_freqRotMOD = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_damaged_freqRotMOD = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_damaged_freqRotMOD = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_damaged_freqRotMOD = FM_strkpln.FM_R.Mx_strkpln_R;
    My_damaged_freqRotMOD = FM_strkpln.FM_R.My_strkpln_R;
    Mz_damaged_freqRotMOD = FM_strkpln.FM_R.Mz_strkpln_R;
    
    Fx_NONdamaged_freqRotMOD = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_NONdamaged_freqRotMOD = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_NONdamaged_freqRotMOD = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_NONdamaged_freqRotMOD = FM_strkpln.FM_R.Mx_strkpln_R;
    My_NONdamaged_freqRotMOD = FM_strkpln.FM_R.My_strkpln_R;
    Mz_NONdamaged_freqRotMOD = FM_strkpln.FM_R.Mz_strkpln_R;
    
    % remove forces and torques at removed part of wing (n>N_cut)
    Fx_damaged_freqRotMOD(N_cut+1:end,:) = nan;
    Fy_damaged_freqRotMOD(N_cut+1:end,:) = nan;
    Fz_damaged_freqRotMOD(N_cut+1:end,:) = nan;
    
    Mx_damaged_freqRotMOD(N_cut+1:end,:) = nan;
    My_damaged_freqRotMOD(N_cut+1:end,:) = nan;
    Mz_damaged_freqRotMOD(N_cut+1:end,:) = nan;

    % total force & torque damaged & intact wing
    Fx_damagedwing_freqRotMOD(:,i) = (nansum(Fx_damaged_freqRotMOD,1))' / Mg_fly;
    Fy_damagedwing_freqRotMOD(:,i) = (nansum(Fy_damaged_freqRotMOD,1))' / Mg_fly;
    Fz_damagedwing_freqRotMOD(:,i) = (nansum(Fz_damaged_freqRotMOD,1))' / Mg_fly;
    
    Mx_damagedwing_freqRotMOD(:,i) = (nansum(Mx_damaged_freqRotMOD,1))' / Mg_fly / (l_wing*1000);
    My_damagedwing_freqRotMOD(:,i) = (nansum(My_damaged_freqRotMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_damagedwing_freqRotMOD(:,i) = (nansum(Mz_damaged_freqRotMOD,1))' / Mg_fly / (l_wing*1000);
    
    Fx_intactwing_freqRotMOD(:,i) = (nansum(Fx_intact_freqRotMOD,1))' / Mg_fly;
    Fy_intactwing_freqRotMOD(:,i) = (nansum(Fy_intact_freqRotMOD,1))' / Mg_fly;
    Fz_intactwing_freqRotMOD(:,i) = (nansum(Fz_intact_freqRotMOD,1))' / Mg_fly;
    
    Mx_intactwing_freqRotMOD(:,i) = (nansum(Mx_intact_freqRotMOD,1))' / Mg_fly / (l_wing*1000);
    My_intactwing_freqRotMOD(:,i) = (nansum(My_intact_freqRotMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_intactwing_freqRotMOD(:,i) = (nansum(Mz_intact_freqRotMOD,1))' / Mg_fly / (l_wing*1000);

    % total force & torque
    Fx_damaged_total_freqRotMOD(:,i) = (nansum(Fx_damaged_freqRotMOD,1) + nansum(Fx_intact_freqRotMOD,1))' / Mg_fly;
    Fy_damaged_total_freqRotMOD(:,i) = (nansum(Fy_damaged_freqRotMOD,1) + nansum(Fy_intact_freqRotMOD,1))' / Mg_fly;
    Fz_damaged_total_freqRotMOD(:,i) = (nansum(Fz_damaged_freqRotMOD,1) + nansum(Fz_intact_freqRotMOD,1))' / Mg_fly;
    
    Mx_damaged_total_freqRotMOD(:,i) = (nansum(Mx_damaged_freqRotMOD,1) + nansum(Mx_intact_freqRotMOD,1))' / Mg_fly / (l_wing*1000);
    My_damaged_total_freqRotMOD(:,i) = (nansum(My_damaged_freqRotMOD,1) + nansum(My_intact_freqRotMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_damaged_total_freqRotMOD(:,i) = (nansum(Mz_damaged_freqRotMOD,1) + nansum(Mz_intact_freqRotMOD,1))' / Mg_fly / (l_wing*1000);
    
    Fx_NONdamaged_total_freqRotMOD(:,i) = (nansum(Fx_NONdamaged_freqRotMOD,1) + nansum(Fx_intact_freqRotMOD,1))' / Mg_fly;
    Fy_NONdamaged_total_freqRotMOD(:,i) = (nansum(Fy_NONdamaged_freqRotMOD,1) + nansum(Fy_intact_freqRotMOD,1))' / Mg_fly;
    Fz_NONdamaged_total_freqRotMOD(:,i) = (nansum(Fz_NONdamaged_freqRotMOD,1) + nansum(Fz_intact_freqRotMOD,1))' / Mg_fly;
    
    Mx_NONdamaged_total_freqRotMOD(:,i) = (nansum(Mx_NONdamaged_freqRotMOD,1) + nansum(Mx_intact_freqRotMOD,1))' / Mg_fly / (l_wing*1000);
    My_NONdamaged_total_freqRotMOD(:,i) = (nansum(My_NONdamaged_freqRotMOD,1) + nansum(My_intact_freqRotMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_NONdamaged_total_freqRotMOD(:,i) = (nansum(Mz_NONdamaged_freqRotMOD,1) + nansum(Mz_intact_freqRotMOD,1))' / Mg_fly / (l_wing*1000);
    
    %% WB kin freq & Dev MODs
%     freq = freq_steady;
%     freq = f_clipped_fly;
    freq = freq_now;

    stroke_L = deg2rad(stroke_steady);
    dev_L = deg2rad(dev_intact_now);
    rot_L = deg2rad(rot_steady);

    stroke_R = deg2rad(stroke_steady);
    dev_R = deg2rad(dev_damaged_now);
    rot_R = deg2rad(rot_steady);

    % qs forces & torques
    [ FM_strkpln, ~ ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_on );

    Fx_intact_freqDevMOD = FM_strkpln.FM_L.Fx_strkpln_L;
    Fy_intact_freqDevMOD = FM_strkpln.FM_L.Fy_strkpln_L;
    Fz_intact_freqDevMOD = FM_strkpln.FM_L.Fz_strkpln_L;
    
    Mx_intact_freqDevMOD = FM_strkpln.FM_L.Mx_strkpln_L;
    My_intact_freqDevMOD = FM_strkpln.FM_L.My_strkpln_L;
    Mz_intact_freqDevMOD = FM_strkpln.FM_L.Mz_strkpln_L;
    
    Fx_damaged_freqDevMOD = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_damaged_freqDevMOD = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_damaged_freqDevMOD = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_damaged_freqDevMOD = FM_strkpln.FM_R.Mx_strkpln_R;
    My_damaged_freqDevMOD = FM_strkpln.FM_R.My_strkpln_R;
    Mz_damaged_freqDevMOD = FM_strkpln.FM_R.Mz_strkpln_R;
    
    Fx_NONdamaged_freqDevMOD = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_NONdamaged_freqDevMOD = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_NONdamaged_freqDevMOD = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_NONdamaged_freqDevMOD = FM_strkpln.FM_R.Mx_strkpln_R;
    My_NONdamaged_freqDevMOD = FM_strkpln.FM_R.My_strkpln_R;
    Mz_NONdamaged_freqDevMOD = FM_strkpln.FM_R.Mz_strkpln_R;
    
    % remove forces and torques at removed part of wing (n>N_cut)
    Fx_damaged_freqDevMOD(N_cut+1:end,:) = nan;
    Fy_damaged_freqDevMOD(N_cut+1:end,:) = nan;
    Fz_damaged_freqDevMOD(N_cut+1:end,:) = nan;
    
    Mx_damaged_freqDevMOD(N_cut+1:end,:) = nan;
    My_damaged_freqDevMOD(N_cut+1:end,:) = nan;
    Mz_damaged_freqDevMOD(N_cut+1:end,:) = nan;

    % total force & torque damaged & intact wing
    Fx_damagedwing_freqDevMOD(:,i) = (nansum(Fx_damaged_freqDevMOD,1))' / Mg_fly;
    Fy_damagedwing_freqDevMOD(:,i) = (nansum(Fy_damaged_freqDevMOD,1))' / Mg_fly;
    Fz_damagedwing_freqDevMOD(:,i) = (nansum(Fz_damaged_freqDevMOD,1))' / Mg_fly;
    
    Mx_damagedwing_freqDevMOD(:,i) = (nansum(Mx_damaged_freqDevMOD,1))' / Mg_fly / (l_wing*1000);
    My_damagedwing_freqDevMOD(:,i) = (nansum(My_damaged_freqDevMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_damagedwing_freqDevMOD(:,i) = (nansum(Mz_damaged_freqDevMOD,1))' / Mg_fly / (l_wing*1000);
    
    Fx_intactwing_freqDevMOD(:,i) = (nansum(Fx_intact_freqDevMOD,1))' / Mg_fly;
    Fy_intactwing_freqDevMOD(:,i) = (nansum(Fy_intact_freqDevMOD,1))' / Mg_fly;
    Fz_intactwing_freqDevMOD(:,i) = (nansum(Fz_intact_freqDevMOD,1))' / Mg_fly;
    
    Mx_intactwing_freqDevMOD(:,i) = (nansum(Mx_intact_freqDevMOD,1))' / Mg_fly / (l_wing*1000);
    My_intactwing_freqDevMOD(:,i) = (nansum(My_intact_freqDevMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_intactwing_freqDevMOD(:,i) = (nansum(Mz_intact_freqDevMOD,1))' / Mg_fly / (l_wing*1000);

    % total force & torque
    Fx_damaged_total_freqDevMOD(:,i) = (nansum(Fx_damaged_freqDevMOD,1) + nansum(Fx_intact_freqDevMOD,1))' / Mg_fly;
    Fy_damaged_total_freqDevMOD(:,i) = (nansum(Fy_damaged_freqDevMOD,1) + nansum(Fy_intact_freqDevMOD,1))' / Mg_fly;
    Fz_damaged_total_freqDevMOD(:,i) = (nansum(Fz_damaged_freqDevMOD,1) + nansum(Fz_intact_freqDevMOD,1))' / Mg_fly;
    
    Mx_damaged_total_freqDevMOD(:,i) = (nansum(Mx_damaged_freqDevMOD,1) + nansum(Mx_intact_freqDevMOD,1))' / Mg_fly / (l_wing*1000);
    My_damaged_total_freqDevMOD(:,i) = (nansum(My_damaged_freqDevMOD,1) + nansum(My_intact_freqDevMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_damaged_total_freqDevMOD(:,i) = (nansum(Mz_damaged_freqDevMOD,1) + nansum(Mz_intact_freqDevMOD,1))' / Mg_fly / (l_wing*1000);
    
    Fx_NONdamaged_total_freqDevMOD(:,i) = (nansum(Fx_NONdamaged_freqDevMOD,1) + nansum(Fx_intact_freqDevMOD,1))' / Mg_fly;
    Fy_NONdamaged_total_freqDevMOD(:,i) = (nansum(Fy_NONdamaged_freqDevMOD,1) + nansum(Fy_intact_freqDevMOD,1))' / Mg_fly;
    Fz_NONdamaged_total_freqDevMOD(:,i) = (nansum(Fz_NONdamaged_freqDevMOD,1) + nansum(Fz_intact_freqDevMOD,1))' / Mg_fly;
    
    Mx_NONdamaged_total_freqDevMOD(:,i) = (nansum(Mx_NONdamaged_freqDevMOD,1) + nansum(Mx_intact_freqDevMOD,1))' / Mg_fly / (l_wing*1000);
    My_NONdamaged_total_freqDevMOD(:,i) = (nansum(My_NONdamaged_freqDevMOD,1) + nansum(My_intact_freqDevMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_NONdamaged_total_freqDevMOD(:,i) = (nansum(Mz_NONdamaged_freqDevMOD,1) + nansum(Mz_intact_freqDevMOD,1))' / Mg_fly / (l_wing*1000);
    
    %% WB kin freq & rot & dev MODs
%     freq = freq_steady;
%     freq = f_clipped_fly;
    freq = freq_now;

    stroke_L = deg2rad(stroke_steady);
    dev_L = deg2rad(dev_intact_now);
    rot_L = deg2rad(rot_intact_now);

    stroke_R = deg2rad(stroke_steady);
    dev_R = deg2rad(dev_damaged_now);
    rot_R = deg2rad(rot_damaged_now);

    % qs forces & torques
    [ FM_strkpln, ~ ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_on );

    Fx_intact_freqRotDevMOD = FM_strkpln.FM_L.Fx_strkpln_L;
    Fy_intact_freqRotDevMOD = FM_strkpln.FM_L.Fy_strkpln_L;
    Fz_intact_freqRotDevMOD = FM_strkpln.FM_L.Fz_strkpln_L;
    
    Mx_intact_freqRotDevMOD = FM_strkpln.FM_L.Mx_strkpln_L;
    My_intact_freqRotDevMOD = FM_strkpln.FM_L.My_strkpln_L;
    Mz_intact_freqRotDevMOD = FM_strkpln.FM_L.Mz_strkpln_L;
    
    Fx_damaged_freqRotDevMOD = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_damaged_freqRotDevMOD = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_damaged_freqRotDevMOD = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_damaged_freqRotDevMOD = FM_strkpln.FM_R.Mx_strkpln_R;
    My_damaged_freqRotDevMOD = FM_strkpln.FM_R.My_strkpln_R;
    Mz_damaged_freqRotDevMOD = FM_strkpln.FM_R.Mz_strkpln_R;
    
    Fx_NONdamaged_freqRotDevMOD = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_NONdamaged_freqRotDevMOD = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_NONdamaged_freqRotDevMOD = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_NONdamaged_freqRotDevMOD = FM_strkpln.FM_R.Mx_strkpln_R;
    My_NONdamaged_freqRotDevMOD = FM_strkpln.FM_R.My_strkpln_R;
    Mz_NONdamaged_freqRotDevMOD = FM_strkpln.FM_R.Mz_strkpln_R;
    
    % remove forces and torques at removed part of wing (n>N_cut)
    Fx_damaged_freqRotDevMOD(N_cut+1:end,:) = nan;
    Fy_damaged_freqRotDevMOD(N_cut+1:end,:) = nan;
    Fz_damaged_freqRotDevMOD(N_cut+1:end,:) = nan;
    
    Mx_damaged_freqRotDevMOD(N_cut+1:end,:) = nan;
    My_damaged_freqRotDevMOD(N_cut+1:end,:) = nan;
    Mz_damaged_freqRotDevMOD(N_cut+1:end,:) = nan;

    % total force & torque damaged & intact wing
    Fx_damagedwing_freqRotDevMOD(:,i) = (nansum(Fx_damaged_freqRotDevMOD,1))' / Mg_fly;
    Fy_damagedwing_freqRotDevMOD(:,i) = (nansum(Fy_damaged_freqRotDevMOD,1))' / Mg_fly;
    Fz_damagedwing_freqRotDevMOD(:,i) = (nansum(Fz_damaged_freqRotDevMOD,1))' / Mg_fly;
    
    Mx_damagedwing_freqRotDevMOD(:,i) = (nansum(Mx_damaged_freqRotDevMOD,1))' / Mg_fly / (l_wing*1000);
    My_damagedwing_freqRotDevMOD(:,i) = (nansum(My_damaged_freqRotDevMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_damagedwing_freqRotDevMOD(:,i) = (nansum(Mz_damaged_freqRotDevMOD,1))' / Mg_fly / (l_wing*1000);
    
    Fx_intactwing_freqRotDevMOD(:,i) = (nansum(Fx_intact_freqRotDevMOD,1))' / Mg_fly;
    Fy_intactwing_freqRotDevMOD(:,i) = (nansum(Fy_intact_freqRotDevMOD,1))' / Mg_fly;
    Fz_intactwing_freqRotDevMOD(:,i) = (nansum(Fz_intact_freqRotDevMOD,1))' / Mg_fly;
    
    Mx_intactwing_freqRotDevMOD(:,i) = (nansum(Mx_intact_freqRotDevMOD,1))' / Mg_fly / (l_wing*1000);
    My_intactwing_freqRotDevMOD(:,i) = (nansum(My_intact_freqRotDevMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_intactwing_freqRotDevMOD(:,i) = (nansum(Mz_intact_freqRotDevMOD,1))' / Mg_fly / (l_wing*1000);

    % total force & torque
    Fx_damaged_total_freqRotDevMOD(:,i) = (nansum(Fx_damaged_freqRotDevMOD,1) + nansum(Fx_intact_freqRotDevMOD,1))' / Mg_fly;
    Fy_damaged_total_freqRotDevMOD(:,i) = (nansum(Fy_damaged_freqRotDevMOD,1) + nansum(Fy_intact_freqRotDevMOD,1))' / Mg_fly;
    Fz_damaged_total_freqRotDevMOD(:,i) = (nansum(Fz_damaged_freqRotDevMOD,1) + nansum(Fz_intact_freqRotDevMOD,1))' / Mg_fly;
    
    Mx_damaged_total_freqRotDevMOD(:,i) = (nansum(Mx_damaged_freqRotDevMOD,1) + nansum(Mx_intact_freqRotDevMOD,1))' / Mg_fly / (l_wing*1000);
    My_damaged_total_freqRotDevMOD(:,i) = (nansum(My_damaged_freqRotDevMOD,1) + nansum(My_intact_freqRotDevMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_damaged_total_freqRotDevMOD(:,i) = (nansum(Mz_damaged_freqRotDevMOD,1) + nansum(Mz_intact_freqRotDevMOD,1))' / Mg_fly / (l_wing*1000);
    
    Fx_NONdamaged_total_freqRotDevMOD(:,i) = (nansum(Fx_NONdamaged_freqRotDevMOD,1) + nansum(Fx_intact_freqRotDevMOD,1))' / Mg_fly;
    Fy_NONdamaged_total_freqRotDevMOD(:,i) = (nansum(Fy_NONdamaged_freqRotDevMOD,1) + nansum(Fy_intact_freqRotDevMOD,1))' / Mg_fly;
    Fz_NONdamaged_total_freqRotDevMOD(:,i) = (nansum(Fz_NONdamaged_freqRotDevMOD,1) + nansum(Fz_intact_freqRotDevMOD,1))' / Mg_fly;
    
    Mx_NONdamaged_total_freqRotDevMOD(:,i) = (nansum(Mx_NONdamaged_freqRotDevMOD,1) + nansum(Mx_intact_freqRotDevMOD,1))' / Mg_fly / (l_wing*1000);
    My_NONdamaged_total_freqRotDevMOD(:,i) = (nansum(My_NONdamaged_freqRotDevMOD,1) + nansum(My_intact_freqRotDevMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_NONdamaged_total_freqRotDevMOD(:,i) = (nansum(Mz_NONdamaged_freqRotDevMOD,1) + nansum(Mz_intact_freqRotDevMOD,1))' / Mg_fly / (l_wing*1000);
    
    %% WB kin stroke & dev MODs
    freq = freq_steady;
%     freq = f_clipped_fly;
%     freq = freq_now;

    stroke_L = deg2rad(stroke_intact_now);
    dev_L = deg2rad(dev_intact_now);
    rot_L = deg2rad(rot_steady);

    stroke_R = deg2rad(stroke_damaged_now);
    dev_R = deg2rad(dev_damaged_now);
    rot_R = deg2rad(rot_steady);

    % qs forces & torques
    [ FM_strkpln, ~ ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_on );

    Fx_intact_strokeDevMOD = FM_strkpln.FM_L.Fx_strkpln_L;
    Fy_intact_strokeDevMOD = FM_strkpln.FM_L.Fy_strkpln_L;
    Fz_intact_strokeDevMOD = FM_strkpln.FM_L.Fz_strkpln_L;
    
    Mx_intact_strokeDevMOD = FM_strkpln.FM_L.Mx_strkpln_L;
    My_intact_strokeDevMOD = FM_strkpln.FM_L.My_strkpln_L;
    Mz_intact_strokeDevMOD = FM_strkpln.FM_L.Mz_strkpln_L;
    
    Fx_damaged_strokeDevMOD = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_damaged_strokeDevMOD = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_damaged_strokeDevMOD = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_damaged_strokeDevMOD = FM_strkpln.FM_R.Mx_strkpln_R;
    My_damaged_strokeDevMOD = FM_strkpln.FM_R.My_strkpln_R;
    Mz_damaged_strokeDevMOD = FM_strkpln.FM_R.Mz_strkpln_R;
    
    Fx_NONdamaged_strokeDevMOD = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_NONdamaged_strokeDevMOD = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_NONdamaged_strokeDevMOD = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_NONdamaged_strokeDevMOD = FM_strkpln.FM_R.Mx_strkpln_R;
    My_NONdamaged_strokeDevMOD = FM_strkpln.FM_R.My_strkpln_R;
    Mz_NONdamaged_strokeDevMOD = FM_strkpln.FM_R.Mz_strkpln_R;
    
    % remove forces and torques at removed part of wing (n>N_cut)
    Fx_damaged_strokeDevMOD(N_cut+1:end,:) = nan;
    Fy_damaged_strokeDevMOD(N_cut+1:end,:) = nan;
    Fz_damaged_strokeDevMOD(N_cut+1:end,:) = nan;
    
    Mx_damaged_strokeDevMOD(N_cut+1:end,:) = nan;
    My_damaged_strokeDevMOD(N_cut+1:end,:) = nan;
    Mz_damaged_strokeDevMOD(N_cut+1:end,:) = nan;

    % total force & torque damaged & intact wing
    Fx_damagedwing_strokeDevMOD(:,i) = (nansum(Fx_damaged_strokeDevMOD,1))' / Mg_fly;
    Fy_damagedwing_strokeDevMOD(:,i) = (nansum(Fy_damaged_strokeDevMOD,1))' / Mg_fly;
    Fz_damagedwing_strokeDevMOD(:,i) = (nansum(Fz_damaged_strokeDevMOD,1))' / Mg_fly;
    
    Mx_damagedwing_strokeDevMOD(:,i) = (nansum(Mx_damaged_strokeDevMOD,1))' / Mg_fly / (l_wing*1000);
    My_damagedwing_strokeDevMOD(:,i) = (nansum(My_damaged_strokeDevMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_damagedwing_strokeDevMOD(:,i) = (nansum(Mz_damaged_strokeDevMOD,1))' / Mg_fly / (l_wing*1000);
    
    Fx_intactwing_strokeDevMOD(:,i) = (nansum(Fx_intact_strokeDevMOD,1))' / Mg_fly;
    Fy_intactwing_strokeDevMOD(:,i) = (nansum(Fy_intact_strokeDevMOD,1))' / Mg_fly;
    Fz_intactwing_strokeDevMOD(:,i) = (nansum(Fz_intact_strokeDevMOD,1))' / Mg_fly;
    
    Mx_intactwing_strokeDevMOD(:,i) = (nansum(Mx_intact_strokeDevMOD,1))' / Mg_fly / (l_wing*1000);
    My_intactwing_strokeDevMOD(:,i) = (nansum(My_intact_strokeDevMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_intactwing_strokeDevMOD(:,i) = (nansum(Mz_intact_strokeDevMOD,1))' / Mg_fly / (l_wing*1000);

    % total force & torque
    Fx_damaged_total_strokeDevMOD(:,i) = (nansum(Fx_damaged_strokeDevMOD,1) + nansum(Fx_intact_strokeDevMOD,1))' / Mg_fly;
    Fy_damaged_total_strokeDevMOD(:,i) = (nansum(Fy_damaged_strokeDevMOD,1) + nansum(Fy_intact_strokeDevMOD,1))' / Mg_fly;
    Fz_damaged_total_strokeDevMOD(:,i) = (nansum(Fz_damaged_strokeDevMOD,1) + nansum(Fz_intact_strokeDevMOD,1))' / Mg_fly;
    
    Mx_damaged_total_strokeDevMOD(:,i) = (nansum(Mx_damaged_strokeDevMOD,1) + nansum(Mx_intact_strokeDevMOD,1))' / Mg_fly / (l_wing*1000);
    My_damaged_total_strokeDevMOD(:,i) = (nansum(My_damaged_strokeDevMOD,1) + nansum(My_intact_strokeDevMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_damaged_total_strokeDevMOD(:,i) = (nansum(Mz_damaged_strokeDevMOD,1) + nansum(Mz_intact_strokeDevMOD,1))' / Mg_fly / (l_wing*1000);
    
    Fx_NONdamaged_total_strokeDevMOD(:,i) = (nansum(Fx_NONdamaged_strokeDevMOD,1) + nansum(Fx_intact_strokeDevMOD,1))' / Mg_fly;
    Fy_NONdamaged_total_strokeDevMOD(:,i) = (nansum(Fy_NONdamaged_strokeDevMOD,1) + nansum(Fy_intact_strokeDevMOD,1))' / Mg_fly;
    Fz_NONdamaged_total_strokeDevMOD(:,i) = (nansum(Fz_NONdamaged_strokeDevMOD,1) + nansum(Fz_intact_strokeDevMOD,1))' / Mg_fly;
    
    Mx_NONdamaged_total_strokeDevMOD(:,i) = (nansum(Mx_NONdamaged_strokeDevMOD,1) + nansum(Mx_intact_strokeDevMOD,1))' / Mg_fly / (l_wing*1000);
    My_NONdamaged_total_strokeDevMOD(:,i) = (nansum(My_NONdamaged_strokeDevMOD,1) + nansum(My_intact_strokeDevMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_NONdamaged_total_strokeDevMOD(:,i) = (nansum(Mz_NONdamaged_strokeDevMOD,1) + nansum(Mz_intact_strokeDevMOD,1))' / Mg_fly / (l_wing*1000);
    
    %% WB kin stroke & Rot MODs
    freq = freq_steady;
%     freq = f_clipped_fly;
%     freq = freq_now;

    stroke_L = deg2rad(stroke_intact_now);
    dev_L = deg2rad(dev_steady);
    rot_L = deg2rad(rot_intact_now);

    stroke_R = deg2rad(stroke_damaged_now);
    dev_R = deg2rad(dev_steady);
    rot_R = deg2rad(rot_damaged_now);

    % qs forces & torques
    [ FM_strkpln, ~ ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_on );

    Fx_intact_strokeRotMOD = FM_strkpln.FM_L.Fx_strkpln_L;
    Fy_intact_strokeRotMOD = FM_strkpln.FM_L.Fy_strkpln_L;
    Fz_intact_strokeRotMOD = FM_strkpln.FM_L.Fz_strkpln_L;
    
    Mx_intact_strokeRotMOD = FM_strkpln.FM_L.Mx_strkpln_L;
    My_intact_strokeRotMOD = FM_strkpln.FM_L.My_strkpln_L;
    Mz_intact_strokeRotMOD = FM_strkpln.FM_L.Mz_strkpln_L;
    
    Fx_damaged_strokeRotMOD = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_damaged_strokeRotMOD = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_damaged_strokeRotMOD = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_damaged_strokeRotMOD = FM_strkpln.FM_R.Mx_strkpln_R;
    My_damaged_strokeRotMOD = FM_strkpln.FM_R.My_strkpln_R;
    Mz_damaged_strokeRotMOD = FM_strkpln.FM_R.Mz_strkpln_R;
    
    Fx_NONdamaged_strokeRotMOD = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_NONdamaged_strokeRotMOD = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_NONdamaged_strokeRotMOD = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_NONdamaged_strokeRotMOD = FM_strkpln.FM_R.Mx_strkpln_R;
    My_NONdamaged_strokeRotMOD = FM_strkpln.FM_R.My_strkpln_R;
    Mz_NONdamaged_strokeRotMOD = FM_strkpln.FM_R.Mz_strkpln_R;
    
    % remove forces and torques at removed part of wing (n>N_cut)
    Fx_damaged_strokeRotMOD(N_cut+1:end,:) = nan;
    Fy_damaged_strokeRotMOD(N_cut+1:end,:) = nan;
    Fz_damaged_strokeRotMOD(N_cut+1:end,:) = nan;
    
    Mx_damaged_strokeRotMOD(N_cut+1:end,:) = nan;
    My_damaged_strokeRotMOD(N_cut+1:end,:) = nan;
    Mz_damaged_strokeRotMOD(N_cut+1:end,:) = nan;

    % total force & torque damaged & intact wing
    Fx_damagedwing_strokeRotMOD(:,i) = (nansum(Fx_damaged_strokeRotMOD,1))' / Mg_fly;
    Fy_damagedwing_strokeRotMOD(:,i) = (nansum(Fy_damaged_strokeRotMOD,1))' / Mg_fly;
    Fz_damagedwing_strokeRotMOD(:,i) = (nansum(Fz_damaged_strokeRotMOD,1))' / Mg_fly;
    
    Mx_damagedwing_strokeRotMOD(:,i) = (nansum(Mx_damaged_strokeRotMOD,1))' / Mg_fly / (l_wing*1000);
    My_damagedwing_strokeRotMOD(:,i) = (nansum(My_damaged_strokeRotMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_damagedwing_strokeRotMOD(:,i) = (nansum(Mz_damaged_strokeRotMOD,1))' / Mg_fly / (l_wing*1000);
    
    Fx_intactwing_strokeRotMOD(:,i) = (nansum(Fx_intact_strokeRotMOD,1))' / Mg_fly;
    Fy_intactwing_strokeRotMOD(:,i) = (nansum(Fy_intact_strokeRotMOD,1))' / Mg_fly;
    Fz_intactwing_strokeRotMOD(:,i) = (nansum(Fz_intact_strokeRotMOD,1))' / Mg_fly;
    
    Mx_intactwing_strokeRotMOD(:,i) = (nansum(Mx_intact_strokeRotMOD,1))' / Mg_fly / (l_wing*1000);
    My_intactwing_strokeRotMOD(:,i) = (nansum(My_intact_strokeRotMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_intactwing_strokeRotMOD(:,i) = (nansum(Mz_intact_strokeRotMOD,1))' / Mg_fly / (l_wing*1000);

    % total force & torque
    Fx_damaged_total_strokeRotMOD(:,i) = (nansum(Fx_damaged_strokeRotMOD,1) + nansum(Fx_intact_strokeRotMOD,1))' / Mg_fly;
    Fy_damaged_total_strokeRotMOD(:,i) = (nansum(Fy_damaged_strokeRotMOD,1) + nansum(Fy_intact_strokeRotMOD,1))' / Mg_fly;
    Fz_damaged_total_strokeRotMOD(:,i) = (nansum(Fz_damaged_strokeRotMOD,1) + nansum(Fz_intact_strokeRotMOD,1))' / Mg_fly;
    
    Mx_damaged_total_strokeRotMOD(:,i) = (nansum(Mx_damaged_strokeRotMOD,1) + nansum(Mx_intact_strokeRotMOD,1))' / Mg_fly / (l_wing*1000);
    My_damaged_total_strokeRotMOD(:,i) = (nansum(My_damaged_strokeRotMOD,1) + nansum(My_intact_strokeRotMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_damaged_total_strokeRotMOD(:,i) = (nansum(Mz_damaged_strokeRotMOD,1) + nansum(Mz_intact_strokeRotMOD,1))' / Mg_fly / (l_wing*1000);
    
    Fx_NONdamaged_total_strokeRotMOD(:,i) = (nansum(Fx_NONdamaged_strokeRotMOD,1) + nansum(Fx_intact_strokeRotMOD,1))' / Mg_fly;
    Fy_NONdamaged_total_strokeRotMOD(:,i) = (nansum(Fy_NONdamaged_strokeRotMOD,1) + nansum(Fy_intact_strokeRotMOD,1))' / Mg_fly;
    Fz_NONdamaged_total_strokeRotMOD(:,i) = (nansum(Fz_NONdamaged_strokeRotMOD,1) + nansum(Fz_intact_strokeRotMOD,1))' / Mg_fly;
    
    Mx_NONdamaged_total_strokeRotMOD(:,i) = (nansum(Mx_NONdamaged_strokeRotMOD,1) + nansum(Mx_intact_strokeRotMOD,1))' / Mg_fly / (l_wing*1000);
    My_NONdamaged_total_strokeRotMOD(:,i) = (nansum(My_NONdamaged_strokeRotMOD,1) + nansum(My_intact_strokeRotMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_NONdamaged_total_strokeRotMOD(:,i) = (nansum(Mz_NONdamaged_strokeRotMOD,1) + nansum(Mz_intact_strokeRotMOD,1))' / Mg_fly / (l_wing*1000);
    
    %% WB kin stroke & dev & rot MODs
    freq = freq_steady;
%     freq = f_clipped_fly;
%     freq = freq_now;

    stroke_L = deg2rad(stroke_intact_now);
    dev_L    = deg2rad(dev_intact_now);
    rot_L    = deg2rad(rot_intact_now);

    stroke_R = deg2rad(stroke_damaged_now);
    dev_R    = deg2rad(dev_damaged_now);
    rot_R    = deg2rad(rot_damaged_now);

    % qs forces & torques
    [ FM_strkpln, ~ ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_on );

    Fx_intact_strokeDevRotMOD = FM_strkpln.FM_L.Fx_strkpln_L;
    Fy_intact_strokeDevRotMOD = FM_strkpln.FM_L.Fy_strkpln_L;
    Fz_intact_strokeDevRotMOD = FM_strkpln.FM_L.Fz_strkpln_L;
    
    Mx_intact_strokeDevRotMOD = FM_strkpln.FM_L.Mx_strkpln_L;
    My_intact_strokeDevRotMOD = FM_strkpln.FM_L.My_strkpln_L;
    Mz_intact_strokeDevRotMOD = FM_strkpln.FM_L.Mz_strkpln_L;
    
    Fx_damaged_strokeDevRotMOD = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_damaged_strokeDevRotMOD = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_damaged_strokeDevRotMOD = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_damaged_strokeDevRotMOD = FM_strkpln.FM_R.Mx_strkpln_R;
    My_damaged_strokeDevRotMOD = FM_strkpln.FM_R.My_strkpln_R;
    Mz_damaged_strokeDevRotMOD = FM_strkpln.FM_R.Mz_strkpln_R;
    
    Fx_NONdamaged_strokeDevRotMOD = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_NONdamaged_strokeDevRotMOD = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_NONdamaged_strokeDevRotMOD = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_NONdamaged_strokeDevRotMOD = FM_strkpln.FM_R.Mx_strkpln_R;
    My_NONdamaged_strokeDevRotMOD = FM_strkpln.FM_R.My_strkpln_R;
    Mz_NONdamaged_strokeDevRotMOD = FM_strkpln.FM_R.Mz_strkpln_R;
    
    % remove forces and torques at removed part of wing (n>N_cut)
    Fx_damaged_strokeDevRotMOD(N_cut+1:end,:) = nan;
    Fy_damaged_strokeDevRotMOD(N_cut+1:end,:) = nan;
    Fz_damaged_strokeDevRotMOD(N_cut+1:end,:) = nan;
    
    Mx_damaged_strokeDevRotMOD(N_cut+1:end,:) = nan;
    My_damaged_strokeDevRotMOD(N_cut+1:end,:) = nan;
    Mz_damaged_strokeDevRotMOD(N_cut+1:end,:) = nan;

    % total force & torque damaged & intact wing
    Fx_damagedwing_strokeDevRotMOD(:,i) = (nansum(Fx_damaged_strokeDevRotMOD,1))' / Mg_fly;
    Fy_damagedwing_strokeDevRotMOD(:,i) = (nansum(Fy_damaged_strokeDevRotMOD,1))' / Mg_fly;
    Fz_damagedwing_strokeDevRotMOD(:,i) = (nansum(Fz_damaged_strokeDevRotMOD,1))' / Mg_fly;
    
    Mx_damagedwing_strokeDevRotMOD(:,i) = (nansum(Mx_damaged_strokeDevRotMOD,1))' / Mg_fly / (l_wing*1000);
    My_damagedwing_strokeDevRotMOD(:,i) = (nansum(My_damaged_strokeDevRotMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_damagedwing_strokeDevRotMOD(:,i) = (nansum(Mz_damaged_strokeDevRotMOD,1))' / Mg_fly / (l_wing*1000);
    
    Fx_intactwing_strokeDevRotMOD(:,i) = (nansum(Fx_intact_strokeDevRotMOD,1))' / Mg_fly;
    Fy_intactwing_strokeDevRotMOD(:,i) = (nansum(Fy_intact_strokeDevRotMOD,1))' / Mg_fly;
    Fz_intactwing_strokeDevRotMOD(:,i) = (nansum(Fz_intact_strokeDevRotMOD,1))' / Mg_fly;
    
    Mx_intactwing_strokeDevRotMOD(:,i) = (nansum(Mx_intact_strokeDevRotMOD,1))' / Mg_fly / (l_wing*1000);
    My_intactwing_strokeDevRotMOD(:,i) = (nansum(My_intact_strokeDevRotMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_intactwing_strokeDevRotMOD(:,i) = (nansum(Mz_intact_strokeDevRotMOD,1))' / Mg_fly / (l_wing*1000);

    % total force & torque
    Fx_damaged_total_strokeDevRotMOD(:,i) = (nansum(Fx_damaged_strokeDevRotMOD,1) + nansum(Fx_intact_strokeDevRotMOD,1))' / Mg_fly;
    Fy_damaged_total_strokeDevRotMOD(:,i) = (nansum(Fy_damaged_strokeDevRotMOD,1) + nansum(Fy_intact_strokeDevRotMOD,1))' / Mg_fly;
    Fz_damaged_total_strokeDevRotMOD(:,i) = (nansum(Fz_damaged_strokeDevRotMOD,1) + nansum(Fz_intact_strokeDevRotMOD,1))' / Mg_fly;
    
    Mx_damaged_total_strokeDevRotMOD(:,i) = (nansum(Mx_damaged_strokeDevRotMOD,1) + nansum(Mx_intact_strokeDevRotMOD,1))' / Mg_fly / (l_wing*1000);
    My_damaged_total_strokeDevRotMOD(:,i) = (nansum(My_damaged_strokeDevRotMOD,1) + nansum(My_intact_strokeDevRotMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_damaged_total_strokeDevRotMOD(:,i) = (nansum(Mz_damaged_strokeDevRotMOD,1) + nansum(Mz_intact_strokeDevRotMOD,1))' / Mg_fly / (l_wing*1000);
    
    Fx_NONdamaged_total_strokeDevRotMOD(:,i) = (nansum(Fx_NONdamaged_strokeDevRotMOD,1) + nansum(Fx_intact_strokeDevRotMOD,1))' / Mg_fly;
    Fy_NONdamaged_total_strokeDevRotMOD(:,i) = (nansum(Fy_NONdamaged_strokeDevRotMOD,1) + nansum(Fy_intact_strokeDevRotMOD,1))' / Mg_fly;
    Fz_NONdamaged_total_strokeDevRotMOD(:,i) = (nansum(Fz_NONdamaged_strokeDevRotMOD,1) + nansum(Fz_intact_strokeDevRotMOD,1))' / Mg_fly;
    
    Mx_NONdamaged_total_strokeDevRotMOD(:,i) = (nansum(Mx_NONdamaged_strokeDevRotMOD,1) + nansum(Mx_intact_strokeDevRotMOD,1))' / Mg_fly / (l_wing*1000);
    My_NONdamaged_total_strokeDevRotMOD(:,i) = (nansum(My_NONdamaged_strokeDevRotMOD,1) + nansum(My_intact_strokeDevRotMOD,1))' / Mg_fly / (l_wing*1000);
    Mz_NONdamaged_total_strokeDevRotMOD(:,i) = (nansum(Mz_NONdamaged_strokeDevRotMOD,1) + nansum(Mz_intact_strokeDevRotMOD,1))' / Mg_fly / (l_wing*1000);
    
    
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

    %% wingbeat average forces & torques: freq & Rot MOD
    Fx_damaged_mean_freqRotMOD = nanmean(Fx_damaged_total_freqRotMOD);
    Fy_damaged_mean_freqRotMOD = nanmean(Fy_damaged_total_freqRotMOD);
    Fz_damaged_mean_freqRotMOD = nanmean(Fz_damaged_total_freqRotMOD);
    
    Mx_damaged_mean_freqRotMOD = nanmean(Mx_damaged_total_freqRotMOD);
    My_damaged_mean_freqRotMOD = nanmean(My_damaged_total_freqRotMOD);
    Mz_damaged_mean_freqRotMOD = nanmean(Mz_damaged_total_freqRotMOD);
    
    Fx_NONdamaged_mean_freqRotMOD = nanmean(Fx_NONdamaged_total_freqRotMOD);
    Fy_NONdamaged_mean_freqRotMOD = nanmean(Fy_NONdamaged_total_freqRotMOD);
    Fz_NONdamaged_mean_freqRotMOD = nanmean(Fz_NONdamaged_total_freqRotMOD);
    
    Mx_NONdamaged_mean_freqRotMOD = nanmean(Mx_NONdamaged_total_freqRotMOD);
    My_NONdamaged_mean_freqRotMOD = nanmean(My_NONdamaged_total_freqRotMOD);
    Mz_NONdamaged_mean_freqRotMOD = nanmean(Mz_NONdamaged_total_freqRotMOD);

    %% wingbeat average forces & torques: freq & Dev MOD
    Fx_damaged_mean_freqDevMOD = nanmean(Fx_damaged_total_freqDevMOD);
    Fy_damaged_mean_freqDevMOD = nanmean(Fy_damaged_total_freqDevMOD);
    Fz_damaged_mean_freqDevMOD = nanmean(Fz_damaged_total_freqDevMOD);
    
    Mx_damaged_mean_freqDevMOD = nanmean(Mx_damaged_total_freqDevMOD);
    My_damaged_mean_freqDevMOD = nanmean(My_damaged_total_freqDevMOD);
    Mz_damaged_mean_freqDevMOD = nanmean(Mz_damaged_total_freqDevMOD);
    
    Fx_NONdamaged_mean_freqDevMOD = nanmean(Fx_NONdamaged_total_freqDevMOD);
    Fy_NONdamaged_mean_freqDevMOD = nanmean(Fy_NONdamaged_total_freqDevMOD);
    Fz_NONdamaged_mean_freqDevMOD = nanmean(Fz_NONdamaged_total_freqDevMOD);
    
    Mx_NONdamaged_mean_freqDevMOD = nanmean(Mx_NONdamaged_total_freqDevMOD);
    My_NONdamaged_mean_freqDevMOD = nanmean(My_NONdamaged_total_freqDevMOD);
    Mz_NONdamaged_mean_freqDevMOD = nanmean(Mz_NONdamaged_total_freqDevMOD);

    %% wingbeat average forces & torques: freq & Rot & dev MOD
    Fx_damaged_mean_freqRotDevMOD = nanmean(Fx_damaged_total_freqRotDevMOD);
    Fy_damaged_mean_freqRotDevMOD = nanmean(Fy_damaged_total_freqRotDevMOD);
    Fz_damaged_mean_freqRotDevMOD = nanmean(Fz_damaged_total_freqRotDevMOD);
    
    Mx_damaged_mean_freqRotDevMOD = nanmean(Mx_damaged_total_freqRotDevMOD);
    My_damaged_mean_freqRotDevMOD = nanmean(My_damaged_total_freqRotDevMOD);
    Mz_damaged_mean_freqRotDevMOD = nanmean(Mz_damaged_total_freqRotDevMOD);
    
    Fx_NONdamaged_mean_freqRotDevMOD = nanmean(Fx_NONdamaged_total_freqRotDevMOD);
    Fy_NONdamaged_mean_freqRotDevMOD = nanmean(Fy_NONdamaged_total_freqRotDevMOD);
    Fz_NONdamaged_mean_freqRotDevMOD = nanmean(Fz_NONdamaged_total_freqRotDevMOD);
    
    Mx_NONdamaged_mean_freqRotDevMOD = nanmean(Mx_NONdamaged_total_freqRotDevMOD);
    My_NONdamaged_mean_freqRotDevMOD = nanmean(My_NONdamaged_total_freqRotDevMOD);
    Mz_NONdamaged_mean_freqRotDevMOD = nanmean(Mz_NONdamaged_total_freqRotDevMOD);

    %% wingbeat average forces & torques: stroke & dev MOD
    Fx_damaged_mean_strokeDevMOD = nanmean(Fx_damaged_total_strokeDevMOD);
    Fy_damaged_mean_strokeDevMOD = nanmean(Fy_damaged_total_strokeDevMOD);
    Fz_damaged_mean_strokeDevMOD = nanmean(Fz_damaged_total_strokeDevMOD);
    
    Mx_damaged_mean_strokeDevMOD = nanmean(Mx_damaged_total_strokeDevMOD);
    My_damaged_mean_strokeDevMOD = nanmean(My_damaged_total_strokeDevMOD);
    Mz_damaged_mean_strokeDevMOD = nanmean(Mz_damaged_total_strokeDevMOD);
    
    Fx_NONdamaged_mean_strokeDevMOD = nanmean(Fx_NONdamaged_total_strokeDevMOD);
    Fy_NONdamaged_mean_strokeDevMOD = nanmean(Fy_NONdamaged_total_strokeDevMOD);
    Fz_NONdamaged_mean_strokeDevMOD = nanmean(Fz_NONdamaged_total_strokeDevMOD);
    
    Mx_NONdamaged_mean_strokeDevMOD = nanmean(Mx_NONdamaged_total_strokeDevMOD);
    My_NONdamaged_mean_strokeDevMOD = nanmean(My_NONdamaged_total_strokeDevMOD);
    Mz_NONdamaged_mean_strokeDevMOD = nanmean(Mz_NONdamaged_total_strokeDevMOD);

    %% wingbeat average forces & torques: stroke & Rot MOD
    Fx_damaged_mean_strokeRotMOD = nanmean(Fx_damaged_total_strokeRotMOD);
    Fy_damaged_mean_strokeRotMOD = nanmean(Fy_damaged_total_strokeRotMOD);
    Fz_damaged_mean_strokeRotMOD = nanmean(Fz_damaged_total_strokeRotMOD);
    
    Mx_damaged_mean_strokeRotMOD = nanmean(Mx_damaged_total_strokeRotMOD);
    My_damaged_mean_strokeRotMOD = nanmean(My_damaged_total_strokeRotMOD);
    Mz_damaged_mean_strokeRotMOD = nanmean(Mz_damaged_total_strokeRotMOD);
    
    Fx_NONdamaged_mean_strokeRotMOD = nanmean(Fx_NONdamaged_total_strokeRotMOD);
    Fy_NONdamaged_mean_strokeRotMOD = nanmean(Fy_NONdamaged_total_strokeRotMOD);
    Fz_NONdamaged_mean_strokeRotMOD = nanmean(Fz_NONdamaged_total_strokeRotMOD);
    
    Mx_NONdamaged_mean_strokeRotMOD = nanmean(Mx_NONdamaged_total_strokeRotMOD);
    My_NONdamaged_mean_strokeRotMOD = nanmean(My_NONdamaged_total_strokeRotMOD);
    Mz_NONdamaged_mean_strokeRotMOD = nanmean(Mz_NONdamaged_total_strokeRotMOD);

    %% wingbeat average forces & torques: stroke & dev & rot MOD
    Fx_damaged_mean_strokeDevRotMOD = nanmean(Fx_damaged_total_strokeDevRotMOD);
    Fy_damaged_mean_strokeDevRotMOD = nanmean(Fy_damaged_total_strokeDevRotMOD);
    Fz_damaged_mean_strokeDevRotMOD = nanmean(Fz_damaged_total_strokeDevRotMOD);
    
    Mx_damaged_mean_strokeDevRotMOD = nanmean(Mx_damaged_total_strokeDevRotMOD);
    My_damaged_mean_strokeDevRotMOD = nanmean(My_damaged_total_strokeDevRotMOD);
    Mz_damaged_mean_strokeDevRotMOD = nanmean(Mz_damaged_total_strokeDevRotMOD);
    
    Fx_NONdamaged_mean_strokeDevRotMOD = nanmean(Fx_NONdamaged_total_strokeDevRotMOD);
    Fy_NONdamaged_mean_strokeDevRotMOD = nanmean(Fy_NONdamaged_total_strokeDevRotMOD);
    Fz_NONdamaged_mean_strokeDevRotMOD = nanmean(Fz_NONdamaged_total_strokeDevRotMOD);
    
    Mx_NONdamaged_mean_strokeDevRotMOD = nanmean(Mx_NONdamaged_total_strokeDevRotMOD);
    My_NONdamaged_mean_strokeDevRotMOD = nanmean(My_NONdamaged_total_strokeDevRotMOD);
    Mz_NONdamaged_mean_strokeDevRotMOD = nanmean(Mz_NONdamaged_total_strokeDevRotMOD);
    
    %% add roll
%     Fy&Fz
    % ALL+ROLL
    Fy_damaged_mean_all_rolled = -Fz_damaged_mean_all .* sind(roll_all') + Fy_damaged_mean_all .* cosd(roll_all');
    Fz_damaged_mean_all_rolled =  Fz_damaged_mean_all .* cosd(roll_all') + Fy_damaged_mean_all .* sind(roll_all');

    Fy_NONdamaged_mean_all_rolled = -Fz_NONdamaged_mean_all .* sind(roll_all') + Fy_NONdamaged_mean_all .* cosd(roll_all');
    Fz_NONdamaged_mean_all_rolled =  Fz_NONdamaged_mean_all .* cosd(roll_all') + Fy_NONdamaged_mean_all .* sind(roll_all');

    % STEADY+ROLL
    Fy_damaged_mean_steady_rolled = -Fz_damaged_mean_steady .* sind(roll_all') + Fy_damaged_mean_steady .* cosd(roll_all');
    Fz_damaged_mean_steady_rolled =  Fz_damaged_mean_steady .* cosd(roll_all') + Fy_damaged_mean_steady .* sind(roll_all');

    Fy_NONdamaged_mean_steady_rolled = -Fz_NONdamaged_mean_steady .* sind(roll_all') + Fy_NONdamaged_mean_steady .* cosd(roll_all');
    Fz_NONdamaged_mean_steady_rolled =  Fz_NONdamaged_mean_steady .* cosd(roll_all') + Fy_NONdamaged_mean_steady .* sind(roll_all');

    % STROKE+ROLL
    Fy_damaged_mean_strokeMOD_rolled = -Fz_damaged_mean_strokeMOD .* sind(roll_all') + Fy_damaged_mean_strokeMOD .* cosd(roll_all');
    Fz_damaged_mean_strokeMOD_rolled =  Fz_damaged_mean_strokeMOD .* cosd(roll_all') + Fy_damaged_mean_strokeMOD .* sind(roll_all');

    Fy_NONdamaged_mean_strokeMOD_rolled = -Fz_NONdamaged_mean_strokeMOD .* sind(roll_all') + Fy_NONdamaged_mean_strokeMOD .* cosd(roll_all');
    Fz_NONdamaged_mean_strokeMOD_rolled =  Fz_NONdamaged_mean_strokeMOD .* cosd(roll_all') + Fy_NONdamaged_mean_strokeMOD .* sind(roll_all');

    % STROKE+DEV+ROLL
    Fy_damaged_mean_strokeDevMOD_rolled = -Fz_damaged_mean_strokeDevMOD .* sind(roll_all') + Fy_damaged_mean_strokeDevMOD .* cosd(roll_all');
    Fz_damaged_mean_strokeDevMOD_rolled =  Fz_damaged_mean_strokeDevMOD .* cosd(roll_all') + Fy_damaged_mean_strokeDevMOD .* sind(roll_all');

    Fy_NONdamaged_mean_strokeDevMOD_rolled = -Fz_NONdamaged_mean_strokeDevMOD .* sind(roll_all') + Fy_NONdamaged_mean_strokeDevMOD .* cosd(roll_all');
    Fz_NONdamaged_mean_strokeDevMOD_rolled =  Fz_NONdamaged_mean_strokeDevMOD .* cosd(roll_all') + Fy_NONdamaged_mean_strokeDevMOD .* sind(roll_all');

    % STROKE+ROT+ROLL
    Fy_damaged_mean_strokeRotMOD_rolled = -Fz_damaged_mean_strokeRotMOD .* sind(roll_all') + Fy_damaged_mean_strokeRotMOD .* cosd(roll_all');
    Fz_damaged_mean_strokeRotMOD_rolled =  Fz_damaged_mean_strokeRotMOD .* cosd(roll_all') + Fy_damaged_mean_strokeRotMOD .* sind(roll_all');

    Fy_NONdamaged_mean_strokeRotMOD_rolled = -Fz_NONdamaged_mean_strokeRotMOD .* sind(roll_all') + Fy_NONdamaged_mean_strokeRotMOD .* cosd(roll_all');
    Fz_NONdamaged_mean_strokeRotMOD_rolled =  Fz_NONdamaged_mean_strokeRotMOD .* cosd(roll_all') + Fy_NONdamaged_mean_strokeRotMOD .* sind(roll_all');

    % STROKE+DEV+ROT+ROLL
    Fy_damaged_mean_strokeDevRotMOD_rolled = -Fz_damaged_mean_strokeDevRotMOD .* sind(roll_all') + Fy_damaged_mean_strokeDevRotMOD .* cosd(roll_all');
    Fz_damaged_mean_strokeDevRotMOD_rolled =  Fz_damaged_mean_strokeDevRotMOD .* cosd(roll_all') + Fy_damaged_mean_strokeDevRotMOD .* sind(roll_all');

    Fy_NONdamaged_mean_strokeDevRotMOD_rolled = -Fz_NONdamaged_mean_strokeDevRotMOD .* sind(roll_all') + Fy_NONdamaged_mean_strokeDevRotMOD .* cosd(roll_all');
    Fz_NONdamaged_mean_strokeDevRotMOD_rolled =  Fz_NONdamaged_mean_strokeDevRotMOD .* cosd(roll_all') + Fy_NONdamaged_mean_strokeDevRotMOD .* sind(roll_all');

%     My&Mz
    % ALL+ROLL
    My_damaged_mean_all_rolled = -Mz_damaged_mean_all .* sind(roll_all') + My_damaged_mean_all .* cosd(roll_all');
    Mz_damaged_mean_all_rolled =  Mz_damaged_mean_all .* cosd(roll_all') + (My_damaged_mean_all-My_damaged_mean_all(S3ratios==1)) .* sind(roll_all');

    My_NONdamaged_mean_all_rolled = -Mz_NONdamaged_mean_all .* sind(roll_all') + My_NONdamaged_mean_all .* cosd(roll_all');
    Mz_NONdamaged_mean_all_rolled =  Mz_NONdamaged_mean_all .* cosd(roll_all') + (My_NONdamaged_mean_all-My_damaged_mean_all(S3ratios==1)) .* sind(roll_all');

    % STEADY+ROLL
    My_damaged_mean_steady_rolled = -Mz_damaged_mean_steady .* sind(roll_all') + My_damaged_mean_steady .* cosd(roll_all');
    Mz_damaged_mean_steady_rolled =  Mz_damaged_mean_steady .* cosd(roll_all') + (My_damaged_mean_steady-My_damaged_mean_all(S3ratios==1)) .* sind(roll_all');

    My_NONdamaged_mean_steady_rolled = -Mz_NONdamaged_mean_steady .* sind(roll_all') + My_NONdamaged_mean_steady .* cosd(roll_all');
    Mz_NONdamaged_mean_steady_rolled =  Mz_NONdamaged_mean_steady .* cosd(roll_all') + (My_NONdamaged_mean_steady-My_damaged_mean_all(S3ratios==1)) .* sind(roll_all');

    % STROKE+ROLL
    My_damaged_mean_strokeMOD_rolled = -Mz_damaged_mean_strokeMOD .* sind(roll_all') + My_damaged_mean_strokeMOD .* cosd(roll_all');
    Mz_damaged_mean_strokeMOD_rolled =  Mz_damaged_mean_strokeMOD .* cosd(roll_all') + (My_damaged_mean_strokeMOD-My_damaged_mean_all(S3ratios==1)) .* sind(roll_all');

    My_NONdamaged_mean_strokeMOD_rolled = -Mz_NONdamaged_mean_strokeMOD .* sind(roll_all') + My_NONdamaged_mean_strokeMOD .* cosd(roll_all');
    Mz_NONdamaged_mean_strokeMOD_rolled =  Mz_NONdamaged_mean_strokeMOD .* cosd(roll_all') + (My_NONdamaged_mean_strokeMOD-My_damaged_mean_all(S3ratios==1)) .* sind(roll_all');

    % FREQ+ROLL
    My_damaged_mean_freqMOD_rolled = -Mz_damaged_mean_freqMOD .* sind(roll_all') + My_damaged_mean_freqMOD .* cosd(roll_all');
    Mz_damaged_mean_freqMOD_rolled =  Mz_damaged_mean_freqMOD .* cosd(roll_all') + (My_damaged_mean_freqMOD-My_damaged_mean_all(S3ratios==1)) .* sind(roll_all');

    My_NONdamaged_mean_freqMOD_rolled = -Mz_NONdamaged_mean_freqMOD .* sind(roll_all') + My_NONdamaged_mean_freqMOD .* cosd(roll_all');
    Mz_NONdamaged_mean_freqMOD_rolled =  Mz_NONdamaged_mean_freqMOD .* cosd(roll_all') + (My_NONdamaged_mean_freqMOD-My_damaged_mean_all(S3ratios==1)) .* sind(roll_all');

    % FREQ+ROLL+ROT
    My_damaged_mean_freqRotMOD_rolled = -Mz_damaged_mean_freqRotMOD .* sind(roll_all') + My_damaged_mean_freqRotMOD .* cosd(roll_all');
    Mz_damaged_mean_freqRotMOD_rolled =  Mz_damaged_mean_freqRotMOD .* cosd(roll_all') + (My_damaged_mean_freqRotMOD-My_damaged_mean_all(S3ratios==1)) .* sind(roll_all');

    My_NONdamaged_mean_freqRotMOD_rolled = -Mz_NONdamaged_mean_freqRotMOD .* sind(roll_all') + My_NONdamaged_mean_freqRotMOD .* cosd(roll_all');
    Mz_NONdamaged_mean_freqRotMOD_rolled =  Mz_NONdamaged_mean_freqRotMOD .* cosd(roll_all') + (My_NONdamaged_mean_freqRotMOD-My_damaged_mean_all(S3ratios==1)) .* sind(roll_all');

    % FREQ+ROLL+ROT+DEV
    My_damaged_mean_freqRotDevMOD_rolled = -Mz_damaged_mean_freqRotDevMOD .* sind(roll_all') + My_damaged_mean_freqRotDevMOD .* cosd(roll_all');
    Mz_damaged_mean_freqRotDevMOD_rolled =  Mz_damaged_mean_freqRotDevMOD .* cosd(roll_all') + (My_damaged_mean_freqRotDevMOD-My_damaged_mean_all(S3ratios==1)) .* sind(roll_all');

    My_NONdamaged_mean_freqRotDevMOD_rolled = -Mz_NONdamaged_mean_freqRotDevMOD .* sind(roll_all') + My_NONdamaged_mean_freqRotDevMOD .* cosd(roll_all');
    Mz_NONdamaged_mean_freqRotDevMOD_rolled =  Mz_NONdamaged_mean_freqRotDevMOD .* cosd(roll_all') + (My_NONdamaged_mean_freqRotDevMOD-My_damaged_mean_all(S3ratios==1)) .* sind(roll_all');

    % STROKE+DEV+ROLL
    My_damaged_mean_strokeDevMOD_rolled = -Mz_damaged_mean_strokeDevMOD .* sind(roll_all') + My_damaged_mean_strokeDevMOD .* cosd(roll_all');
    Mz_damaged_mean_strokeDevMOD_rolled =  Mz_damaged_mean_strokeDevMOD .* cosd(roll_all') + (My_damaged_mean_strokeDevMOD-My_damaged_mean_all(S3ratios==1)) .* sind(roll_all');

    My_NONdamaged_mean_strokeDevMOD_rolled = -Mz_NONdamaged_mean_strokeDevMOD .* sind(roll_all') + My_NONdamaged_mean_strokeDevMOD .* cosd(roll_all');
    Mz_NONdamaged_mean_strokeDevMOD_rolled =  Mz_NONdamaged_mean_strokeDevMOD .* cosd(roll_all') + (My_NONdamaged_mean_strokeDevMOD-My_damaged_mean_all(S3ratios==1)) .* sind(roll_all');

    % STROKE+ROT+ROLL
    My_damaged_mean_strokeRotMOD_rolled = -Mz_damaged_mean_strokeRotMOD .* sind(roll_all') + My_damaged_mean_strokeRotMOD .* cosd(roll_all');
    Mz_damaged_mean_strokeRotMOD_rolled =  Mz_damaged_mean_strokeRotMOD .* cosd(roll_all') + (My_damaged_mean_strokeRotMOD-My_damaged_mean_all(S3ratios==1)) .* sind(roll_all');

    My_NONdamaged_mean_strokeRotMOD_rolled = -Mz_NONdamaged_mean_strokeRotMOD .* sind(roll_all') + My_NONdamaged_mean_strokeRotMOD .* cosd(roll_all');
    Mz_NONdamaged_mean_strokeRotMOD_rolled =  Mz_NONdamaged_mean_strokeRotMOD .* cosd(roll_all') + (My_NONdamaged_mean_strokeRotMOD-My_damaged_mean_all(S3ratios==1)) .* sind(roll_all');

    % STROKE+DEV+ROT+ROLL
    My_damaged_mean_strokeDevRotMOD_rolled = -Mz_damaged_mean_strokeDevRotMOD .* sind(roll_all') + My_damaged_mean_strokeDevRotMOD .* cosd(roll_all');
    Mz_damaged_mean_strokeDevRotMOD_rolled =  Mz_damaged_mean_strokeDevRotMOD .* cosd(roll_all') + (My_damaged_mean_strokeDevRotMOD-My_damaged_mean_all(S3ratios==1)) .* sind(roll_all');

    My_NONdamaged_mean_strokeDevRotMOD_rolled = -Mz_NONdamaged_mean_strokeDevRotMOD .* sind(roll_all') + My_NONdamaged_mean_strokeDevRotMOD .* cosd(roll_all');
    Mz_NONdamaged_mean_strokeDevRotMOD_rolled =  Mz_NONdamaged_mean_strokeDevRotMOD .* cosd(roll_all') + (My_NONdamaged_mean_strokeDevRotMOD-My_damaged_mean_all(S3ratios==1)) .* sind(roll_all');

    %% save data
    if fit_type == 1
        save(['allMODs_TipClip_indivFreqNrollFit_Asym10_roton',num2str(rot_on),'.mat'])
        if plot_on == 1
            mkdir(['qsModel_FnM_indivFreqNrollFit_Asym10_roton',num2str(rot_on)])
            cd(['qsModel_FnM_indivFreqNrollFit_Asym10_roton',num2str(rot_on)])
        end
    elseif fit_type == 2
        save(['allMODs_TipClip_indivFreqNrollFit_spline999_roton',num2str(rot_on),'.mat'])
        if plot_on == 1
            mkdir(['qsModel_FnM_indivFreqNrollFit_spline999_roton',num2str(rot_on)])
            cd(['qsModel_FnM_indivFreqNrollFit_spline999_roton',num2str(rot_on)])
        end
    elseif fit_type == 3
        save(['allMODs_TipClip_indivFreqNrollFit_linear_roton',num2str(rot_on),'.mat'])
        if plot_on == 1
            mkdir(['qsModel_FnM_indivFreqNrollFit_linear_roton',num2str(rot_on)])
            cd(['qsModel_FnM_indivFreqNrollFit_linear_roton',num2str(rot_on)])
        end
    elseif fit_type == 4
        save(['allMODs_TipClip_indivFreqNrollFit_power2_roton',num2str(rot_on),'.mat'])
        if plot_on == 1
            mkdir(['qsModel_FnM_indivFreqNrollFit_power2_roton',num2str(rot_on)])
            cd(['qsModel_FnM_indivFreqNrollFit_power2_roton',num2str(rot_on)])
        end
    end
    
    %% plot
    if plot_on == 1

    % ALL MODs VS NO MODs
    figure
    subplot(1,2,1)
    hold on
%     plot(S2ratios,Fx_damaged_mean_all,'o-k','markersize',10,'markerfacecolor','b')
%     plot(S2ratios,Fy_damaged_mean_all,'o-k','markersize',10,'markerfacecolor','r')
%     plot(S2ratios,Fz_damaged_mean_all,'o-k','markersize',10,'markerfacecolor','g')
    
    plot(S2ratios,Fx_damaged_mean_all,'o-k','markersize',10,'markerfacecolor','b')
    plot(S2ratios,Fy_damaged_mean_all_rolled,'o-k','markersize',10,'markerfacecolor','r')
    plot(S2ratios,Fz_damaged_mean_all_rolled,'o-k','markersize',10,'markerfacecolor','g')
    
    plot(S2ratios,Fx_damaged_mean_steady,'o-k','markersize',10,'markerfacecolor','c')
    plot(S2ratios,Fy_damaged_mean_steady,'o-k','markersize',10,'markerfacecolor','m')
    plot(S2ratios,Fz_damaged_mean_steady,'o-k','markersize',10,'markerfacecolor','y')
    
    legend('x-axis Mod','y-axis Mod','z-axis Mod','x-axis steady','y-axis steady','z-axis steady','location','S')
    xlabel('S2 ratio')
    ylabel('normalized forces F/mg')
    axis([0.5 1 -2.5 .25])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-2.5:.5:1)
    
    subplot(1,2,2)
    hold on
    plot(S3ratios,Mx_damaged_mean_all,'o-k','markersize',10,'markerfacecolor','b')
    plot(S3ratios,My_damaged_mean_all_rolled-My_damaged_mean_all_rolled(S3ratios==1),'o-k','markersize',10,'markerfacecolor','r')
    plot(S3ratios,Mz_damaged_mean_all_rolled,'o-k','markersize',10,'markerfacecolor','g')

    plot(S3ratios,Mx_damaged_mean_steady,'o-k','markersize',10,'markerfacecolor','c')
    plot(S3ratios,My_damaged_mean_steady-My_damaged_mean_all(S3ratios==1),'o-k','markersize',10,'markerfacecolor','m')
    plot(S3ratios,Mz_damaged_mean_steady,'o-k','markersize',10,'markerfacecolor','y')
    
    xlabel('S3 ratio')
    ylabel('normalized torques T/mgl')
    axis([0.5 1 -.15 .2])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-.15:.05:.2)    
    
    saveas(gca,['FnM_WBmod_TipClip.fig'])
    saveas(gca,['FnM_WBmod_TipClip.png'])
    plot2svg(['FnM_WBmod_TipClip.svg'])

    % compare MODs Fz & Mx CUMULATIVE FREQ FIRST & ROLL LAST
    figure
    subplot(1,2,1)
    hold on
    plot(S2ratios,Fz_damaged_mean_steady,'o-k','markersize',10,'markerfacecolor','w')
    plot(S2ratios,Fz_damaged_mean_freqMOD,'o-k','markersize',10,'markerfacecolor','y')
    plot(S2ratios,Fz_damaged_mean_freqDevMOD,'o-k','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S2ratios,Fz_damaged_mean_freqRotDevMOD,'o-k','markersize',10,'markerfacecolor','r')
    plot(S2ratios,Fz_damaged_mean_all,'o-k','markersize',10,'markerfacecolor',[.5 0 0])
    plot(S2ratios,Fz_damaged_mean_all_rolled,'o-k','markersize',10,'markerfacecolor','k')
   
    legend('none','+freq','+dev','+rot','+stroke','+roll','location','NE')
    xlabel('S2 ratio')
    ylabel('normalized vertical force Fz/mg')
    axis([0.5 1 -1.5 0])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-1.5:.25:1)
    
    subplot(1,2,2)
    hold on
    plot(S3ratios,Mx_damaged_mean_steady,'o-k','markersize',10,'markerfacecolor','w')
    plot(S3ratios,Mx_damaged_mean_freqMOD,'o-k','markersize',10,'markerfacecolor','y')
    plot(S3ratios,Mx_damaged_mean_freqDevMOD,'o-k','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S3ratios,Mx_damaged_mean_freqRotDevMOD,'o-k','markersize',10,'markerfacecolor','r')
    plot(S3ratios,Mx_damaged_mean_all,'o-k','markersize',10,'markerfacecolor',[.5 0 0])
    plot(S3ratios,Mx_damaged_mean_all,'o-k','markersize',10,'markerfacecolor','k')
    
    xlabel('S3 ratio')
    ylabel('normalized roll torque Tx/mgl')
    axis([0.5 1 -.05 .25])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-.05:.05:.3)    

    saveas(gca,['Fz_Mx_WBmodCombi_freqDevRotStroke_TipClip.fig'])
    saveas(gca,['Fz_Mx_WBmodCombi_freqDevRotStroke_TipClip.png'])
    plot2svg(['Fz_Mx_WBmodCombi_freqDevRotStroke_TipClip.svg'])

    % compare MODs Fz & Mx CUMULATIVE STROKE FIRST & FREQ LAST
    figure
    subplot(1,2,1)
    hold on
%     plot(S2ratios,Fz_damaged_mean_steady,'o-k','markersize',10,'markerfacecolor','w')
%     plot(S2ratios,Fz_damaged_mean_strokeMOD,'o-k','markersize',10,'markerfacecolor','y')
%     plot(S2ratios,Fz_damaged_mean_strokeDevMOD,'o-k','markersize',10,'markerfacecolor',[1 .5 0])
%     plot(S2ratios,Fz_damaged_mean_strokeDevRotMOD,'o-k','markersize',10,'markerfacecolor','r')
%     plot(S2ratios,Fz_damaged_mean_all,'o-k','markersize',10,'markerfacecolor',[.5 0 0])
%     plot(S2ratios,Fz_damaged_mean_all_rolled,'o-k','markersize',10,'markerfacecolor','k')
    
    plot(S2ratios,Fz_damaged_mean_steady,'o-k','markersize',10,'markerfacecolor','w')
    plot(S2ratios,Fz_damaged_mean_strokeMOD,'o-k','markersize',10,'markerfacecolor','y')
    plot(S2ratios,Fz_damaged_mean_strokeRotMOD,'o-k','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S2ratios,Fz_damaged_mean_strokeDevRotMOD,'o-k','markersize',10,'markerfacecolor','r')
    plot(S2ratios,Fz_damaged_mean_strokeDevRotMOD_rolled,'o-k','markersize',10,'markerfacecolor',[.5 0 0])
    plot(S2ratios,Fz_damaged_mean_all_rolled,'o-k','markersize',10,'markerfacecolor','k')
    
    legend('none','+stroke','+rot','+dev','+roll','+freq','location','NE')
    xlabel('S2 ratio')
    ylabel('normalized vertical force Fz/mg')
    axis([0.5 1 -1.5 0])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-1.5:.25:1)
    
    subplot(1,2,2)
    hold on
    plot(S3ratios,Mx_damaged_mean_steady,'o-k','markersize',10,'markerfacecolor','w')
    plot(S3ratios,Mx_damaged_mean_strokeMOD,'o-k','markersize',10,'markerfacecolor','y')
    plot(S3ratios,Mx_damaged_mean_strokeRotMOD,'o-k','markersize',10,'markerfacecolor',[1 .5 0])
    plot(S3ratios,Mx_damaged_mean_strokeDevRotMOD,'o-k','markersize',10,'markerfacecolor','r')
    plot(S3ratios,Mx_damaged_mean_strokeDevRotMOD,'o-k','markersize',10,'markerfacecolor',[.5 0 0])
    plot(S3ratios,Mx_damaged_mean_all,'o-k','markersize',10,'markerfacecolor','k')
    
    xlabel('S3 ratio')
    ylabel('normalized roll torque Tx/mgl')
    axis([0.5 1 -.05 .25])
    set(gca,'xtick',0:.5:1)
    set(gca,'ytick',-.05:.05:.3)    

    saveas(gca,['Fz_Mx_WBmodCombi_strokeRotDevRollFreq_TipClip.fig'])
    saveas(gca,['Fz_Mx_WBmodCombi_strokeRotDevRollFreq_TipClip.png'])
    plot2svg(['Fz_Mx_WBmodCombi_strokeRotDevRollFreq_TipClip.svg'])
    
    cd ..
    end
end