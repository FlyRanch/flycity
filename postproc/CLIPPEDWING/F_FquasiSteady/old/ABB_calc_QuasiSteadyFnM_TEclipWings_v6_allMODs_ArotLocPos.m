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

rot_lift_on=1;
% rot_lift_on=0;

save_on = 1;
% save_on = 0;

all_means_on = 0;
steady_on = 0;
all_on = 0;
freq_on = 0;
stroke_on = 0;
dev_on = 0;
wingrot_on = 0;

all_means_on = 1;
steady_on = 1;
all_on = 1;
freq_on = 1;
stroke_on = 1;
dev_on = 1;
wingrot_on = 1;

% location of peak rotation angle
loc_set = 3
if loc_set == 1
    maxpeak_nr = 1;
    minpeak_nr = 2;
    
elseif loc_set == 2
    maxpeak_nr = 2;
    minpeak_nr = 3;

elseif loc_set == 3
    maxpeak_nr = 3;
    minpeak_nr = 1;
end

%% constants
nr_sect = settings.nr_chord_sect;
nr_timepoints = settings.nr_timepoints;

n_rot_lim1 = .2*nr_timepoints;
n_rot_lim2 = .75*nr_timepoints;
n_rot_lim3 = .5*nr_timepoints;

Mg_fly = body_model.Mg_fly;                             % [kg]
l_wing = wing_model.length*1e-3;                        % [m]
area_wing = wing_model.area*1e-6;                       % [m^2]
massperArea = wing_model.mass/(wing_model.area*1e-6);   % [kg/m^2]

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
h_sect = wing_model.length/nr_sect*1e-3;    % [m]
x_sect = wing_model.x_LE_R'*1e-3;           % [m]
y_sect = wing_model.y_sect_R(:,2)*1e-3;     % [m]
c_sect = wing_model.chords_R'*1e-3;         % [m]

c_max = max(c_sect);

S2_sect_element = 1/12 .* c_sect .* h_sect.^3;
S3_sect_element = 1/36 .* c_sect .* h_sect.^4;

S2_sect = c_sect .* h_sect .* y_sect.^2;
S3_sect = c_sect .* h_sect .* y_sect.^3;

S2_intact = nansum(S2_sect);
S3_intact = nansum(S3_sect);

sol = subs(solAdAiRatio,S2,1);
sol = subs(sol,S3,1);
S2S3AmpRatioFunc_NONclipped = eval(sol);

% wing section mass
mass_wing_sect_intact = massperArea * h_sect * c_sect;
mass_virtual_sect_intact = 1/4* pi * rho_air * h_sect * c_sect.^2;
mass_tot_sect_intact = mass_wing_sect_intact + mass_virtual_sect_intact;

% Iwing & Iaddedmass
Iwing_sect = 1/12 * massperArea * h_sect * c_sect.^3 + massperArea * h_sect * c_sect .* (1/2*c_sect - x_sect).^2;
Iaddedmass_sect = 1/64 * pi * rho_air * h_sect * c_sect.^4 + 1/4* pi * rho_air * h_sect * c_sect.^2 .* (1/2*c_sect - x_sect).^2;
Itot_sect = Iwing_sect + Iaddedmass_sect;

Iwing_intact = nansum(Iwing_sect);
Iaddedmass_intact = nansum(Iaddedmass_sect);
Itot_intact = nansum(Itot_sect);

%% wing clipped & kin MODs
TEclip_ratio_min = 0.45;
TEclip_Dratio = -0.05;
TEclip_ratios = [1:TEclip_Dratio:TEclip_ratio_min]';

TEclip_ratios = [1 .8 .65 .55 .475 .425]'

for freq_asymFitNr = 10


%% loop with different cuts
for i = 1:length(TEclip_ratios)

    %% S2 & S3 of damaged wing
    TEclip_ratio_now = TEclip_ratios(i)
    c_cut = TEclip_ratio_now * c_max;
    
    c_sect_cut = c_sect;
    c_sect_cut(c_sect_cut>c_cut) = c_cut;
    
    S2_sect = c_sect_cut .* h_sect .* y_sect.^2;
    S3_sect = c_sect_cut .* h_sect .* y_sect.^3;

    S2_damaged = nansum(S2_sect);
    S3_damaged = nansum(S3_sect);
    
    S2ratio = S2_damaged/S2_intact;
    S3ratio = S3_damaged/S3_intact;
    
    % wing section mass
    mass_wing_sect_damaged = massperArea * h_sect * c_sect;
    mass_virtual_sect_damaged = 1/4* pi * rho_air * h_sect * c_sect.^2;
    mass_tot_sect_damaged = mass_wing_sect_damaged + mass_virtual_sect_damaged;

    % Iwing & Iaddedmass
    Iwing_sect_cut = 1/12 * massperArea * h_sect * c_sect_cut.^3 + massperArea * h_sect * c_sect_cut .* (1/2*c_sect_cut - x_sect).^2;
    Iaddedmass_sect_cut = 1/64 * pi * rho_air * h_sect * c_sect_cut.^4 + 1/4* pi * rho_air * h_sect * c_sect_cut.^2 .* (1/2*c_sect_cut - x_sect).^2;
    Itot_sect_cut = Iwing_sect_cut + Iaddedmass_sect_cut;

    Iwing_damaged = nansum(Iwing_sect_cut);
    Iaddedmass_damaged = nansum(Iaddedmass_sect_cut);
    Itot_damaged = nansum(Itot_sect_cut);

%% write damaged wing data to wingmodel
    wing_model.chords_R = 1e3*c_sect_cut';    % [mm]
    
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
    
    Iwing(i,1) = Iwing_damaged;
    Iaddedmass(i,1) = Iaddedmass_damaged;
    Itotal(i,1) = Itot_damaged;
    
    %% WB kin NO MODs
    if steady_on == 1
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
    [ FM_strkpln, kine ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_lift_on );

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
    
    %% calc spanwise wing torque & rotation angle amplitude 
    chords_L        = 1e-3*wing_model.chords_L';      % [ m ]
    x_LE_L          = 1e-3*wing_model.x_LE_L';        % [ m ]
    chords_R        = 1e-3*wing_model.chords_R';      % [ m ]
    x_LE_R          = 1e-3*wing_model.x_LE_R';        % [ m ]

    rot_ddot_L = kine.eta_ddot_L;
    rot_ddot_R = kine.eta_ddot_R;
    
    Maero_sect_L = FM_strkpln.FM_L.My_L;
    Maero_sect_R = FM_strkpln.FM_R.My_R;
    
    clear ax_L ax_R az_L az_R Mmass_sect_L Mmass_sect_R
    for n = 1:length(rot_L)
        
        % torque arm
        x_cm_L = (x_LE_L-chords_L/2)*sin(rot_L(n));
        z_cm_L = (x_LE_L-chords_L/2)*cos(rot_L(n));

        x_cm_R = (x_LE_R-chords_R/2)*sin(rot_R(n));
        z_cm_R = (x_LE_R-chords_R/2)*cos(rot_R(n));

        ax_L(:,:) = kine.Udot_left(1,n,:)*1e-3;      % [ m/s^2 ]
        az_L(:,:) = kine.Udot_left(3,n,:)*1e-3;      % [ m/s^2 ]

        ax_R(:,:) = kine.Udot_right(1,n,:)*1e-3;     % [ m/s^2 ]
        az_R(:,:) = kine.Udot_right(3,n,:)*1e-3;     % [ m/s^2 ]
        
        Mmass_sect_L(:,n) = mass_tot_sect_intact  .* (x_cm_L.*az_L + z_cm_L.*ax_L);
        Mmass_sect_R(:,n) = mass_tot_sect_damaged .* (x_cm_R.*az_R + z_cm_R.*ax_R);
    end

    Mwing_mass_intactwing_steady(:,i)  = nansum(Mmass_sect_L)' / Mg_fly / l_wing;
    Mwing_mass_damagedwing_steady(:,i) = nansum(Mmass_sect_R)' / Mg_fly / l_wing;
    
    Mwing_inertia_intactwing_steady(:,i)  = Itot_intact * rot_ddot_L / Mg_fly / l_wing;
    Mwing_inertia_damagedwing_steady(:,i) = Itot_damaged * rot_ddot_R / Mg_fly / l_wing;
    
    Mwing_aero_intactwing_steady(:,i)  = nansum(Maero_sect_L)' / Mg_fly / (l_wing*1e3);
    Mwing_aero_damagedwing_steady(:,i) = nansum(Maero_sect_R)' / Mg_fly / (l_wing*1e3);
    
    Mwing_total_intactwing_steady(:,i)  = Mwing_mass_intactwing_steady(:,i) + Mwing_inertia_intactwing_steady(:,i) + Mwing_aero_intactwing_steady(:,i);
    Mwing_total_damagedwing_steady(:,i)  = Mwing_mass_damagedwing_steady(:,i) + Mwing_inertia_damagedwing_steady(:,i) + Mwing_aero_damagedwing_steady(:,i);
    
    %% torque at min & max rotation angle & rotation amplitude
%     [rot_maxpeaks_L,rot_maxpeaklocs_L] = findpeaks(rot_L);
%     n_rot_max_L = rot_maxpeaklocs_L(maxpeak_nr);
%     rot_max_L = rot_L(n_rot_max_L);
% 
%     [rot_minpeaks_L,rot_minpeaklocs_L] = findpeaks(-rot_L);
%     n_rot_min_L = rot_minpeaklocs_L(minpeak_nr);
%     rot_min_L = rot_L(n_rot_min_L);
% 
%     [rot_maxpeaks_R,rot_maxpeaklocs_R] = findpeaks(rot_R);
%     n_rot_max_R = rot_maxpeaklocs_R(maxpeak_nr);
%     rot_max_R = rot_R(n_rot_max_R);
% 
%     [rot_minpeaks_R,rot_minpeaklocs_R] = findpeaks(-rot_R);
%     n_rot_min_R = rot_minpeaklocs_R(minpeak_nr);
%     rot_min_R = rot_R(n_rot_min_R);
% 
%     % plot
% %     subplot(2,1,1)
% %     plot(t_norm,rot_L)
% %     hold on
% %     plot(t_norm(rot_maxpeaklocs_L),rot_maxpeaks_L,'or')
% %     plot(t_norm(rot_minpeaklocs_L),rot_L(rot_minpeaklocs_L),'og')
% %     
% %     plot(t_norm(n_rot_max_L),rot_max_L,'*r')
% %     plot(t_norm(n_rot_min_L),rot_min_L,'*g')
% % 
% %     subplot(2,1,2)
% %     plot(t_norm,rot_R)
% %     hold on
% %     plot(t_norm(rot_maxpeaklocs_R),rot_maxpeaks_R,'or')
% %     plot(t_norm(rot_minpeaklocs_R),rot_R(rot_minpeaklocs_R),'og')
% %     
% %     plot(t_norm(n_rot_max_R),rot_max_R,'*r')
% %     plot(t_norm(n_rot_min_R),rot_min_R,'*g')
%     
    %% store data
%     % rotation angle amplitude etc
%     Arot_intactwing_steady(i,1) = rad2deg(rot_max_L - rot_min_L);
%     Arot_damagedwing_steady(i,1) = rad2deg(rot_max_R - rot_min_R);
%     
%     dMwing_total_intactwing_steady(i,1) = Mwing_total_intactwing_steady(n_rot_max_L,i) - Mwing_total_intactwing_steady(n_rot_min_L,i);
%     dMwing_total_damagedwing_steady(i,1) = Mwing_total_intactwing_steady(n_rot_max_R,i) - Mwing_total_intactwing_steady(n_rot_min_R,i);

    % total force & torque damaged & intact wing
    Fx_damagedwing_steady(:,i) = (nansum(Fx_damaged_steady,1))' / Mg_fly;
    Fy_damagedwing_steady(:,i) = (nansum(Fy_damaged_steady,1))' / Mg_fly;
    Fz_damagedwing_steady(:,i) = (nansum(Fz_damaged_steady,1))' / Mg_fly;
    
    Mx_damagedwing_steady(:,i) = (nansum(Mx_damaged_steady,1))' / Mg_fly / (l_wing*1e3);
    My_damagedwing_steady(:,i) = (nansum(My_damaged_steady,1))' / Mg_fly / (l_wing*1e3);
    Mz_damagedwing_steady(:,i) = (nansum(Mz_damaged_steady,1))' / Mg_fly / (l_wing*1e3);
    
    Fx_intactwing_steady(:,i) = (nansum(Fx_intact_steady,1))' / Mg_fly;
    Fy_intactwing_steady(:,i) = (nansum(Fy_intact_steady,1))' / Mg_fly;
    Fz_intactwing_steady(:,i) = (nansum(Fz_intact_steady,1))' / Mg_fly;
    
    Mx_intactwing_steady(:,i) = (nansum(Mx_intact_steady,1))' / Mg_fly / (l_wing*1e3);
    My_intactwing_steady(:,i) = (nansum(My_intact_steady,1))' / Mg_fly / (l_wing*1e3);
    Mz_intactwing_steady(:,i) = (nansum(Mz_intact_steady,1))' / Mg_fly / (l_wing*1e3);

    % total force & torque
    Fx_damaged_total_steady(:,i) = (nansum(Fx_damaged_steady,1) + nansum(Fx_intact_steady,1))' / Mg_fly;
    Fy_damaged_total_steady(:,i) = (nansum(Fy_damaged_steady,1) + nansum(Fy_intact_steady,1))' / Mg_fly;
    Fz_damaged_total_steady(:,i) = (nansum(Fz_damaged_steady,1) + nansum(Fz_intact_steady,1))' / Mg_fly;
    
    Mx_damaged_total_steady(:,i) = (nansum(Mx_damaged_steady,1) + nansum(Mx_intact_steady,1))' / Mg_fly / (l_wing*1e3);
    My_damaged_total_steady(:,i) = (nansum(My_damaged_steady,1) + nansum(My_intact_steady,1))' / Mg_fly / (l_wing*1e3);
    Mz_damaged_total_steady(:,i) = (nansum(Mz_damaged_steady,1) + nansum(Mz_intact_steady,1))' / Mg_fly / (l_wing*1e3);
    
    Fx_NONdamaged_total_steady(:,i) = (nansum(Fx_NONdamaged_steady,1) + nansum(Fx_intact_steady,1))' / Mg_fly;
    Fy_NONdamaged_total_steady(:,i) = (nansum(Fy_NONdamaged_steady,1) + nansum(Fy_intact_steady,1))' / Mg_fly;
    Fz_NONdamaged_total_steady(:,i) = (nansum(Fz_NONdamaged_steady,1) + nansum(Fz_intact_steady,1))' / Mg_fly;
    
    Mx_NONdamaged_total_steady(:,i) = (nansum(Mx_NONdamaged_steady,1) + nansum(Mx_intact_steady,1))' / Mg_fly / (l_wing*1e3);
    My_NONdamaged_total_steady(:,i) = (nansum(My_NONdamaged_steady,1) + nansum(My_intact_steady,1))' / Mg_fly / (l_wing*1e3);
    Mz_NONdamaged_total_steady(:,i) = (nansum(Mz_NONdamaged_steady,1) + nansum(Mz_intact_steady,1))' / Mg_fly / (l_wing*1e3);
    end
    
    %% WB kin ALL MODs
    if all_on == 1
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
    [ FM_strkpln, kine ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_lift_on );

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
    
    %% calc spanwise wing torque & rotation angle amplitude 
    chords_L        = 1e-3*wing_model.chords_L';      % [ m ]
    x_LE_L          = 1e-3*wing_model.x_LE_L';        % [ m ]
    chords_R        = 1e-3*wing_model.chords_R';      % [ m ]
    x_LE_R          = 1e-3*wing_model.x_LE_R';        % [ m ]

    rot_ddot_L = kine.eta_ddot_L;
    rot_ddot_R = kine.eta_ddot_R;
    
    Maero_sect_L = FM_strkpln.FM_L.My_L;
    Maero_sect_R = FM_strkpln.FM_R.My_R;
    
    clear ax_L ax_R az_L az_R Mmass_sect_L Mmass_sect_R
    for n = 1:length(rot_L)
        
        % torque arm
        x_cm_L = (x_LE_L-chords_L/2)*sin(rot_L(n));
        z_cm_L = (x_LE_L-chords_L/2)*cos(rot_L(n));

        x_cm_R = (x_LE_R-chords_R/2)*sin(rot_R(n));
        z_cm_R = (x_LE_R-chords_R/2)*cos(rot_R(n));

        ax_L(:,:) = kine.Udot_left(1,n,:)*1e-3;      % [ m/s^2 ]
        az_L(:,:) = kine.Udot_left(3,n,:)*1e-3;      % [ m/s^2 ]

        ax_R(:,:) = kine.Udot_right(1,n,:)*1e-3;     % [ m/s^2 ]
        az_R(:,:) = kine.Udot_right(3,n,:)*1e-3;     % [ m/s^2 ]
        
        Mmass_sect_L(:,n) = mass_tot_sect_intact  .* (x_cm_L.*az_L + z_cm_L.*ax_L);
        Mmass_sect_R(:,n) = mass_tot_sect_damaged .* (x_cm_R.*az_R + z_cm_R.*ax_R);
    end

    Mwing_mass_intactwing_all(:,i)  = nansum(Mmass_sect_L)' / Mg_fly / l_wing;
    Mwing_mass_damagedwing_all(:,i) = nansum(Mmass_sect_R)' / Mg_fly / l_wing;
    
    Mwing_inertia_intactwing_all(:,i)  = Itot_intact * rot_ddot_L / Mg_fly / l_wing;
    Mwing_inertia_damagedwing_all(:,i) = Itot_damaged * rot_ddot_R / Mg_fly / l_wing;
    
    Mwing_aero_intactwing_all(:,i)  = nansum(Maero_sect_L)' / Mg_fly / (l_wing*1e3);
    Mwing_aero_damagedwing_all(:,i) = nansum(Maero_sect_R)' / Mg_fly / (l_wing*1e3);
    
    Mwing_total_intactwing_all(:,i)  = Mwing_mass_intactwing_all(:,i) + Mwing_inertia_intactwing_all(:,i) + Mwing_aero_intactwing_all(:,i);
    Mwing_total_damagedwing_all(:,i)  = Mwing_mass_damagedwing_all(:,i) + Mwing_inertia_damagedwing_all(:,i) + Mwing_aero_damagedwing_all(:,i);
    
    %% torque at min & max rotation angle & rotation amplitude

    [rot_maxpeaks_L,rot_maxpeaklocs_L] = findpeaks(rot_L);
    [rot_minpeaks_L,rot_minpeaklocs_L] = findpeaks(-rot_L);
    [rot_maxpeaks_R,rot_maxpeaklocs_R] = findpeaks(rot_R);
    [rot_minpeaks_R,rot_minpeaklocs_R] = findpeaks(-rot_R);

    if loc_set == 0
        rot_max_L = max(rot_L(1:n_rot_lim1));
        rot_min_L = min(rot_L(1:n_rot_lim2));

        n_rot_max_L = find(rot_L==rot_max_L);
        n_rot_min_L = find(rot_L==rot_min_L);

        rot_max_R = max(rot_R(1:n_rot_lim1));
        rot_min_R = min(rot_R(1:n_rot_lim2));

        n_rot_max_R = find(rot_R==rot_max_R);
        n_rot_min_R = find(rot_R==rot_min_R);
        
    else

    %     if maxpeak_nr>length(rot_maxpeaklocs_L)
    %         n_rot_max_L = nan;
    %         rot_max_L = nan;
        if maxpeak_nr == 3
            n_rot_max_L = rot_maxpeaklocs_L(end);
            rot_max_L = rot_L(n_rot_max_L);
            if n_rot_max_L < n_rot_lim3
                n_rot_max_L = nan;
                rot_max_L = nan;
            end
            
        else
            n_rot_max_L = rot_maxpeaklocs_L(maxpeak_nr);
            rot_max_L = rot_L(n_rot_max_L);
            if n_rot_max_L > n_rot_lim3
                n_rot_max_L = nan;
                rot_max_L = nan;
            end
        end


    %     if minpeak_nr>length(rot_minpeaklocs_L)
    %         n_rot_min_L = nan;
    %         rot_min_L = nan;
        if minpeak_nr == 3
            n_rot_min_L = rot_minpeaklocs_L(end);
            rot_min_L = rot_L(n_rot_min_L);
            if n_rot_min_L < n_rot_lim3
                n_rot_min_L = nan;
                rot_min_L = nan;
            end
            
        else
            n_rot_min_L = rot_minpeaklocs_L(minpeak_nr);
            rot_min_L = rot_L(n_rot_min_L);
            if n_rot_min_L > n_rot_lim3
                n_rot_min_L = nan;
                rot_min_L = nan;
            end
        end

    %     if maxpeak_nr>length(rot_maxpeaklocs_R)
    %         n_rot_max_R = nan;
    %         rot_max_R = nan;
        if maxpeak_nr == 3
            n_rot_max_R = rot_maxpeaklocs_R(end);
            rot_max_R = rot_R(n_rot_max_R);
            if n_rot_max_R < n_rot_lim3
                n_rot_max_R = nan;
                rot_max_R = nan;
            end
            
        else
            n_rot_max_R = rot_maxpeaklocs_R(maxpeak_nr);
            rot_max_R = rot_R(n_rot_max_R);
            if n_rot_max_R > n_rot_lim3
                n_rot_max_R = nan;
                rot_max_R = nan;
            end
        end


    %     if minpeak_nr>length(rot_minpeaklocs_R)
    %         n_rot_min_R = nan;
    %         rot_min_R = nan;
        if minpeak_nr == 3
            n_rot_min_R = rot_minpeaklocs_R(end);
            rot_min_R = rot_R(n_rot_min_R);
            if n_rot_min_R < n_rot_lim3
                n_rot_min_R = nan;
                rot_min_R = nan;
            end
            
        else
            n_rot_min_R = rot_minpeaklocs_R(minpeak_nr);
            rot_min_R = rot_R(n_rot_min_R);
            if n_rot_min_R > n_rot_lim3
                n_rot_min_R = nan;
                rot_min_R = nan;
            end
        end
    end

    % plot
    subplot(2,1,1)
    plot(t_norm,rot_L)
    hold on
    plot(t_norm(rot_maxpeaklocs_L),rot_maxpeaks_L,'or')
    plot(t_norm(rot_minpeaklocs_L),rot_L(rot_minpeaklocs_L),'og')
    
    if isnan(n_rot_max_L) == 0
        plot(t_norm(n_rot_max_L),rot_max_L,'*r')
    end
    if isnan(n_rot_min_L) == 0
        plot(t_norm(n_rot_min_L),rot_min_L,'*g')
    end
    
    subplot(2,1,2)
    plot(t_norm,rot_R)
    hold on
    plot(t_norm(rot_maxpeaklocs_R),rot_maxpeaks_R,'or')
    plot(t_norm(rot_minpeaklocs_R),rot_R(rot_minpeaklocs_R),'og')
    
    if isnan(n_rot_max_R) == 0
        plot(t_norm(n_rot_max_R),rot_max_R,'*r')
    end
    if isnan(n_rot_min_R) == 0
        plot(t_norm(n_rot_min_R),rot_min_R,'*g')
    end
    
    %% store data
    % rotation angle amplitude etc
    Arot_intactwing_all(i,1) = rad2deg(rot_max_L - rot_min_L);
    Arot_damagedwing_all(i,1) = rad2deg(rot_max_R - rot_min_R);
    
    if isnan(n_rot_max_L) == 0 & isnan(n_rot_min_L) == 0 
        dMwing_total_intactwing_all(i,1)   = Mwing_total_intactwing_all(n_rot_max_L,i)   - Mwing_total_intactwing_all(n_rot_min_L,i);
        dMwing_mass_intactwing_all(i,1)    = Mwing_mass_intactwing_all(n_rot_max_L,i)    - Mwing_mass_intactwing_all(n_rot_min_L,i);
        dMwing_inertia_intactwing_all(i,1) = Mwing_inertia_intactwing_all(n_rot_max_L,i) - Mwing_inertia_intactwing_all(n_rot_min_L,i);
        dMwing_aero_intactwing_all(i,1)    = Mwing_aero_intactwing_all(n_rot_max_L,i)    - Mwing_aero_intactwing_all(n_rot_min_L,i);
    else
        dMwing_total_intactwing_all(i,1)   = nan;
        dMwing_mass_intactwing_all(i,1)    = nan;
        dMwing_inertia_intactwing_all(i,1) = nan;
        dMwing_aero_intactwing_all(i,1)    = nan;
    end
    
    if isnan(n_rot_max_R) == 0 & isnan(n_rot_min_R) == 0 
        dMwing_total_damagedwing_all(i,1)   = Mwing_total_damagedwing_all(n_rot_max_R,i)   - Mwing_total_damagedwing_all(n_rot_min_R,i);
        dMwing_mass_damagedwing_all(i,1)    = Mwing_mass_damagedwing_all(n_rot_max_R,i)    - Mwing_mass_damagedwing_all(n_rot_min_R,i);
        dMwing_inertia_damagedwing_all(i,1) = Mwing_inertia_damagedwing_all(n_rot_max_R,i) - Mwing_inertia_damagedwing_all(n_rot_min_R,i);
        dMwing_aero_damagedwing_all(i,1)    = Mwing_aero_damagedwing_all(n_rot_max_R,i)    - Mwing_aero_damagedwing_all(n_rot_min_R,i);
    else
        dMwing_total_damagedwing_all(i,1)   = nan;
        dMwing_mass_damagedwing_all(i,1)    = nan;
        dMwing_inertia_damagedwing_all(i,1) = nan;
        dMwing_aero_damagedwing_all(i,1)    = nan;
    end

    n_rot_max_intactwing_all(i,:) = n_rot_max_L;
    n_rot_max_damagedwing_all(i,:) = n_rot_max_R;
    
    n_rot_min_intactwing_all(i,:) = n_rot_min_L;
    n_rot_min_damagedwing_all(i,:) = n_rot_min_R;
    
    rot_max_intactwing_all(i,:) = rot_max_L;
    rot_max_damagedwing_all(i,:) = rot_max_R;
    
    rot_min_intactwing_all(i,:) = rot_min_L;
    rot_min_damagedwing_all(i,:) = rot_min_R;
    
    rot_intactwing_all(i,:) = rot_L;
    rot_damagedwing_all(i,:) = rot_R;
    
    % wingbeat kinematics
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
    
    Mx_damagedwing_all(:,i) = (nansum(Mx_damaged_all,1))' / Mg_fly / (l_wing*1e3);
    My_damagedwing_all(:,i) = (nansum(My_damaged_all,1))' / Mg_fly / (l_wing*1e3);
    Mz_damagedwing_all(:,i) = (nansum(Mz_damaged_all,1))' / Mg_fly / (l_wing*1e3);
    
    Fx_intactwing_all(:,i) = (nansum(Fx_intact_all,1))' / Mg_fly;
    Fy_intactwing_all(:,i) = (nansum(Fy_intact_all,1))' / Mg_fly;
    Fz_intactwing_all(:,i) = (nansum(Fz_intact_all,1))' / Mg_fly;
    
    Mx_intactwing_all(:,i) = (nansum(Mx_intact_all,1))' / Mg_fly / (l_wing*1e3);
    My_intactwing_all(:,i) = (nansum(My_intact_all,1))' / Mg_fly / (l_wing*1e3);
    Mz_intactwing_all(:,i) = (nansum(Mz_intact_all,1))' / Mg_fly / (l_wing*1e3);

    % total force & torque
    Fx_damaged_total_all(:,i) = (nansum(Fx_damaged_all,1) + nansum(Fx_intact_all,1))' / Mg_fly;
    Fy_damaged_total_all(:,i) = (nansum(Fy_damaged_all,1) + nansum(Fy_intact_all,1))' / Mg_fly;
    Fz_damaged_total_all(:,i) = (nansum(Fz_damaged_all,1) + nansum(Fz_intact_all,1))' / Mg_fly;
    
    Mx_damaged_total_all(:,i) = (nansum(Mx_damaged_all,1) + nansum(Mx_intact_all,1))' / Mg_fly / (l_wing*1e3);
    My_damaged_total_all(:,i) = (nansum(My_damaged_all,1) + nansum(My_intact_all,1))' / Mg_fly / (l_wing*1e3);
    Mz_damaged_total_all(:,i) = (nansum(Mz_damaged_all,1) + nansum(Mz_intact_all,1))' / Mg_fly / (l_wing*1e3);
    
    Fx_NONdamaged_total_all(:,i) = (nansum(Fx_NONdamaged_all,1) + nansum(Fx_intact_all,1))' / Mg_fly;
    Fy_NONdamaged_total_all(:,i) = (nansum(Fy_NONdamaged_all,1) + nansum(Fy_intact_all,1))' / Mg_fly;
    Fz_NONdamaged_total_all(:,i) = (nansum(Fz_NONdamaged_all,1) + nansum(Fz_intact_all,1))' / Mg_fly;
    
    Mx_NONdamaged_total_all(:,i) = (nansum(Mx_NONdamaged_all,1) + nansum(Mx_intact_all,1))' / Mg_fly / (l_wing*1e3);
    My_NONdamaged_total_all(:,i) = (nansum(My_NONdamaged_all,1) + nansum(My_intact_all,1))' / Mg_fly / (l_wing*1e3);
    Mz_NONdamaged_total_all(:,i) = (nansum(Mz_NONdamaged_all,1) + nansum(Mz_intact_all,1))' / Mg_fly / (l_wing*1e3);
    end
    
    %% WB kin FREQ MODs
    if freq_on == 1
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
    [ FM_strkpln, kine ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_lift_on );

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
    
    %% calc spanwise wing torque & rotation angle amplitude 
    chords_L        = 1e-3*wing_model.chords_L';      % [ m ]
    x_LE_L          = 1e-3*wing_model.x_LE_L';        % [ m ]
    chords_R        = 1e-3*wing_model.chords_R';      % [ m ]
    x_LE_R          = 1e-3*wing_model.x_LE_R';        % [ m ]

    rot_ddot_L = kine.eta_ddot_L;
    rot_ddot_R = kine.eta_ddot_R;
    
    Maero_sect_L = FM_strkpln.FM_L.My_L;
    Maero_sect_R = FM_strkpln.FM_R.My_R;
    
    clear ax_L ax_R az_L az_R Mmass_sect_L Mmass_sect_R
    for n = 1:length(rot_L)
        
        % torque arm
        x_cm_L = (x_LE_L-chords_L/2)*sin(rot_L(n));
        z_cm_L = (x_LE_L-chords_L/2)*cos(rot_L(n));

        x_cm_R = (x_LE_R-chords_R/2)*sin(rot_R(n));
        z_cm_R = (x_LE_R-chords_R/2)*cos(rot_R(n));

        ax_L(:,:) = kine.Udot_left(1,n,:)*1e-3;      % [ m/s^2 ]
        az_L(:,:) = kine.Udot_left(3,n,:)*1e-3;      % [ m/s^2 ]

        ax_R(:,:) = kine.Udot_right(1,n,:)*1e-3;     % [ m/s^2 ]
        az_R(:,:) = kine.Udot_right(3,n,:)*1e-3;     % [ m/s^2 ]
        
        Mmass_sect_L(:,n) = mass_tot_sect_intact  .* (x_cm_L.*az_L + z_cm_L.*ax_L);
        Mmass_sect_R(:,n) = mass_tot_sect_damaged .* (x_cm_R.*az_R + z_cm_R.*ax_R);
    end

    Mwing_mass_intactwing_freqMOD(:,i)  = nansum(Mmass_sect_L)' / Mg_fly / l_wing;
    Mwing_mass_damagedwing_freqMOD(:,i) = nansum(Mmass_sect_R)' / Mg_fly / l_wing;
    
    Mwing_inertia_intactwing_freqMOD(:,i)  = Itot_intact * rot_ddot_L / Mg_fly / l_wing;
    Mwing_inertia_damagedwing_freqMOD(:,i) = Itot_damaged * rot_ddot_R / Mg_fly / l_wing;
    
    Mwing_aero_intactwing_freqMOD(:,i)  = nansum(Maero_sect_L)' / Mg_fly / (l_wing*1e3);
    Mwing_aero_damagedwing_freqMOD(:,i) = nansum(Maero_sect_R)' / Mg_fly / (l_wing*1e3);
    
    Mwing_total_intactwing_freqMOD(:,i)  = Mwing_mass_intactwing_freqMOD(:,i) + Mwing_inertia_intactwing_freqMOD(:,i) + Mwing_aero_intactwing_freqMOD(:,i);
    Mwing_total_damagedwing_freqMOD(:,i)  = Mwing_mass_damagedwing_freqMOD(:,i) + Mwing_inertia_damagedwing_freqMOD(:,i) + Mwing_aero_damagedwing_freqMOD(:,i);

    %% torque at min & max rotation angle & rotation amplitude
%     [rot_maxpeaks_L,rot_maxpeaklocs_L] = findpeaks(rot_L);
%     n_rot_max_L = rot_maxpeaklocs_L(maxpeak_nr);
%     rot_max_L = rot_L(n_rot_max_L);
% 
%     [rot_minpeaks_L,rot_minpeaklocs_L] = findpeaks(-rot_L);
%     n_rot_min_L = rot_minpeaklocs_L(minpeak_nr);
%     rot_min_L = rot_L(n_rot_min_L);
% 
%     [rot_maxpeaks_R,rot_maxpeaklocs_R] = findpeaks(rot_R);
%     n_rot_max_R = rot_maxpeaklocs_R(maxpeak_nr);
%     rot_max_R = rot_R(n_rot_max_R);
% 
%     [rot_minpeaks_R,rot_minpeaklocs_R] = findpeaks(-rot_R);
%     n_rot_min_R = rot_minpeaklocs_R(minpeak_nr);
%     rot_min_R = rot_R(n_rot_min_R);
% 
%     % plot
% %     subplot(2,1,1)
% %     plot(t_norm,rot_L)
% %     hold on
% %     plot(t_norm(rot_maxpeaklocs_L),rot_maxpeaks_L,'or')
% %     plot(t_norm(rot_minpeaklocs_L),rot_L(rot_minpeaklocs_L),'og')
% %     
% %     plot(t_norm(n_rot_max_L),rot_max_L,'*r')
% %     plot(t_norm(n_rot_min_L),rot_min_L,'*g')
% % 
% %     subplot(2,1,2)
% %     plot(t_norm,rot_R)
% %     hold on
% %     plot(t_norm(rot_maxpeaklocs_R),rot_maxpeaks_R,'or')
% %     plot(t_norm(rot_minpeaklocs_R),rot_R(rot_minpeaklocs_R),'og')
% %     
% %     plot(t_norm(n_rot_max_R),rot_max_R,'*r')
% %     plot(t_norm(n_rot_min_R),rot_min_R,'*g')
%     
    %% store data
%     % rotation angle amplitude etc
%     Arot_intactwing_freqMOD(i,1) = rad2deg(rot_max_L - rot_min_L);
%     Arot_damagedwing_freqMOD(i,1) = rad2deg(rot_max_R - rot_min_R);
%     
%     dMwing_total_intactwing_freqMOD(i,1) = Mwing_total_intactwing_freqMOD(n_rot_max_L,i) - Mwing_total_intactwing_freqMOD(n_rot_min_L,i);
%     dMwing_total_damagedwing_freqMOD(i,1) = Mwing_total_intactwing_freqMOD(n_rot_max_R,i) - Mwing_total_intactwing_freqMOD(n_rot_min_R,i);

    % total force & torque damaged & intact wing
    Fx_damagedwing_freqMOD(:,i) = (nansum(Fx_damaged_freqMOD,1))' / Mg_fly;
    Fy_damagedwing_freqMOD(:,i) = (nansum(Fy_damaged_freqMOD,1))' / Mg_fly;
    Fz_damagedwing_freqMOD(:,i) = (nansum(Fz_damaged_freqMOD,1))' / Mg_fly;
    
    Mx_damagedwing_freqMOD(:,i) = (nansum(Mx_damaged_freqMOD,1))' / Mg_fly / (l_wing*1e3);
    My_damagedwing_freqMOD(:,i) = (nansum(My_damaged_freqMOD,1))' / Mg_fly / (l_wing*1e3);
    Mz_damagedwing_freqMOD(:,i) = (nansum(Mz_damaged_freqMOD,1))' / Mg_fly / (l_wing*1e3);
    
    Fx_intactwing_freqMOD(:,i) = (nansum(Fx_intact_freqMOD,1))' / Mg_fly;
    Fy_intactwing_freqMOD(:,i) = (nansum(Fy_intact_freqMOD,1))' / Mg_fly;
    Fz_intactwing_freqMOD(:,i) = (nansum(Fz_intact_freqMOD,1))' / Mg_fly;
    
    Mx_intactwing_freqMOD(:,i) = (nansum(Mx_intact_freqMOD,1))' / Mg_fly / (l_wing*1e3);
    My_intactwing_freqMOD(:,i) = (nansum(My_intact_freqMOD,1))' / Mg_fly / (l_wing*1e3);
    Mz_intactwing_freqMOD(:,i) = (nansum(Mz_intact_freqMOD,1))' / Mg_fly / (l_wing*1e3);

    % total force & torque
    Fx_damaged_total_freqMOD(:,i) = (nansum(Fx_damaged_freqMOD,1) + nansum(Fx_intact_freqMOD,1))' / Mg_fly;
    Fy_damaged_total_freqMOD(:,i) = (nansum(Fy_damaged_freqMOD,1) + nansum(Fy_intact_freqMOD,1))' / Mg_fly;
    Fz_damaged_total_freqMOD(:,i) = (nansum(Fz_damaged_freqMOD,1) + nansum(Fz_intact_freqMOD,1))' / Mg_fly;
    
    Mx_damaged_total_freqMOD(:,i) = (nansum(Mx_damaged_freqMOD,1) + nansum(Mx_intact_freqMOD,1))' / Mg_fly / (l_wing*1e3);
    My_damaged_total_freqMOD(:,i) = (nansum(My_damaged_freqMOD,1) + nansum(My_intact_freqMOD,1))' / Mg_fly / (l_wing*1e3);
    Mz_damaged_total_freqMOD(:,i) = (nansum(Mz_damaged_freqMOD,1) + nansum(Mz_intact_freqMOD,1))' / Mg_fly / (l_wing*1e3);
    
    Fx_NONdamaged_total_freqMOD(:,i) = (nansum(Fx_NONdamaged_freqMOD,1) + nansum(Fx_intact_freqMOD,1))' / Mg_fly;
    Fy_NONdamaged_total_freqMOD(:,i) = (nansum(Fy_NONdamaged_freqMOD,1) + nansum(Fy_intact_freqMOD,1))' / Mg_fly;
    Fz_NONdamaged_total_freqMOD(:,i) = (nansum(Fz_NONdamaged_freqMOD,1) + nansum(Fz_intact_freqMOD,1))' / Mg_fly;
    
    Mx_NONdamaged_total_freqMOD(:,i) = (nansum(Mx_NONdamaged_freqMOD,1) + nansum(Mx_intact_freqMOD,1))' / Mg_fly / (l_wing*1e3);
    My_NONdamaged_total_freqMOD(:,i) = (nansum(My_NONdamaged_freqMOD,1) + nansum(My_intact_freqMOD,1))' / Mg_fly / (l_wing*1e3);
    Mz_NONdamaged_total_freqMOD(:,i) = (nansum(Mz_NONdamaged_freqMOD,1) + nansum(Mz_intact_freqMOD,1))' / Mg_fly / (l_wing*1e3);
    end
    
    %% WB kin stroke MODs
    if stroke_on == 1
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
    [ FM_strkpln, kine ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_lift_on );

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
    
    %% calc spanwise wing torque & rotation angle amplitude 
    chords_L        = 1e-3*wing_model.chords_L';      % [ m ]
    x_LE_L          = 1e-3*wing_model.x_LE_L';        % [ m ]
    chords_R        = 1e-3*wing_model.chords_R';      % [ m ]
    x_LE_R          = 1e-3*wing_model.x_LE_R';        % [ m ]

    rot_ddot_L = kine.eta_ddot_L;
    rot_ddot_R = kine.eta_ddot_R;
    
    Maero_sect_L = FM_strkpln.FM_L.My_L;
    Maero_sect_R = FM_strkpln.FM_R.My_R;
    
    clear ax_L ax_R az_L az_R Mmass_sect_L Mmass_sect_R
    for n = 1:length(rot_L)
        
        % torque arm
        x_cm_L = (x_LE_L-chords_L/2)*sin(rot_L(n));
        z_cm_L = (x_LE_L-chords_L/2)*cos(rot_L(n));

        x_cm_R = (x_LE_R-chords_R/2)*sin(rot_R(n));
        z_cm_R = (x_LE_R-chords_R/2)*cos(rot_R(n));

        ax_L(:,:) = kine.Udot_left(1,n,:)*1e-3;      % [ m/s^2 ]
        az_L(:,:) = kine.Udot_left(3,n,:)*1e-3;      % [ m/s^2 ]

        ax_R(:,:) = kine.Udot_right(1,n,:)*1e-3;     % [ m/s^2 ]
        az_R(:,:) = kine.Udot_right(3,n,:)*1e-3;     % [ m/s^2 ]
        
        Mmass_sect_L(:,n) = mass_tot_sect_intact  .* (x_cm_L.*az_L + z_cm_L.*ax_L);
        Mmass_sect_R(:,n) = mass_tot_sect_damaged .* (x_cm_R.*az_R + z_cm_R.*ax_R);
    end

    Mwing_mass_intactwing_strokeMOD(:,i)  = nansum(Mmass_sect_L)' / Mg_fly / l_wing;
    Mwing_mass_damagedwing_strokeMOD(:,i) = nansum(Mmass_sect_R)' / Mg_fly / l_wing;
    
    Mwing_inertia_intactwing_strokeMOD(:,i)  = Itot_intact * rot_ddot_L / Mg_fly / l_wing;
    Mwing_inertia_damagedwing_strokeMOD(:,i) = Itot_damaged * rot_ddot_R / Mg_fly / l_wing;
    
    Mwing_aero_intactwing_strokeMOD(:,i)  = nansum(Maero_sect_L)' / Mg_fly / (l_wing*1e3);
    Mwing_aero_damagedwing_strokeMOD(:,i) = nansum(Maero_sect_R)' / Mg_fly / (l_wing*1e3);
    
    Mwing_total_intactwing_strokeMOD(:,i)  = Mwing_mass_intactwing_strokeMOD(:,i) + Mwing_inertia_intactwing_strokeMOD(:,i) + Mwing_aero_intactwing_strokeMOD(:,i);
    Mwing_total_damagedwing_strokeMOD(:,i)  = Mwing_mass_damagedwing_strokeMOD(:,i) + Mwing_inertia_damagedwing_strokeMOD(:,i) + Mwing_aero_damagedwing_strokeMOD(:,i);

    %% torque at min & max rotation angle & rotation amplitude
%     [rot_maxpeaks_L,rot_maxpeaklocs_L] = findpeaks(rot_L);
%     n_rot_max_L = rot_maxpeaklocs_L(maxpeak_nr);
%     rot_max_L = rot_L(n_rot_max_L);
% 
%     [rot_minpeaks_L,rot_minpeaklocs_L] = findpeaks(-rot_L);
%     n_rot_min_L = rot_minpeaklocs_L(minpeak_nr);
%     rot_min_L = rot_L(n_rot_min_L);
% 
%     [rot_maxpeaks_R,rot_maxpeaklocs_R] = findpeaks(rot_R);
%     n_rot_max_R = rot_maxpeaklocs_R(maxpeak_nr);
%     rot_max_R = rot_R(n_rot_max_R);
% 
%     [rot_minpeaks_R,rot_minpeaklocs_R] = findpeaks(-rot_R);
%     n_rot_min_R = rot_minpeaklocs_R(minpeak_nr);
%     rot_min_R = rot_R(n_rot_min_R);
% 
%     % plot
% %     subplot(2,1,1)
% %     plot(t_norm,rot_L)
% %     hold on
% %     plot(t_norm(rot_maxpeaklocs_L),rot_maxpeaks_L,'or')
% %     plot(t_norm(rot_minpeaklocs_L),rot_L(rot_minpeaklocs_L),'og')
% %     
% %     plot(t_norm(n_rot_max_L),rot_max_L,'*r')
% %     plot(t_norm(n_rot_min_L),rot_min_L,'*g')
% % 
% %     subplot(2,1,2)
% %     plot(t_norm,rot_R)
% %     hold on
% %     plot(t_norm(rot_maxpeaklocs_R),rot_maxpeaks_R,'or')
% %     plot(t_norm(rot_minpeaklocs_R),rot_R(rot_minpeaklocs_R),'og')
% %     
% %     plot(t_norm(n_rot_max_R),rot_max_R,'*r')
% %     plot(t_norm(n_rot_min_R),rot_min_R,'*g')
    
    %% store data
%     % rotation angle amplitude etc
%     Arot_intactwing_strokeMOD(i,1) = rad2deg(rot_max_L - rot_min_L);
%     Arot_damagedwing_strokeMOD(i,1) = rad2deg(rot_max_R - rot_min_R);
%     
%     dMwing_total_intactwing_strokeMOD(i,1) = Mwing_total_intactwing_strokeMOD(n_rot_max_L,i) - Mwing_total_intactwing_strokeMOD(n_rot_min_L,i);
%     dMwing_total_damagedwing_strokeMOD(i,1) = Mwing_total_intactwing_strokeMOD(n_rot_max_R,i) - Mwing_total_intactwing_strokeMOD(n_rot_min_R,i);

    % total force & torque damaged & intact wing
    Fx_damagedwing_strokeMOD(:,i) = (nansum(Fx_damaged_strokeMOD,1))' / Mg_fly;
    Fy_damagedwing_strokeMOD(:,i) = (nansum(Fy_damaged_strokeMOD,1))' / Mg_fly;
    Fz_damagedwing_strokeMOD(:,i) = (nansum(Fz_damaged_strokeMOD,1))' / Mg_fly;
    
    Mx_damagedwing_strokeMOD(:,i) = (nansum(Mx_damaged_strokeMOD,1))' / Mg_fly / (l_wing*1e3);
    My_damagedwing_strokeMOD(:,i) = (nansum(My_damaged_strokeMOD,1))' / Mg_fly / (l_wing*1e3);
    Mz_damagedwing_strokeMOD(:,i) = (nansum(Mz_damaged_strokeMOD,1))' / Mg_fly / (l_wing*1e3);
    
    Fx_intactwing_strokeMOD(:,i) = (nansum(Fx_intact_strokeMOD,1))' / Mg_fly;
    Fy_intactwing_strokeMOD(:,i) = (nansum(Fy_intact_strokeMOD,1))' / Mg_fly;
    Fz_intactwing_strokeMOD(:,i) = (nansum(Fz_intact_strokeMOD,1))' / Mg_fly;
    
    Mx_intactwing_strokeMOD(:,i) = (nansum(Mx_intact_strokeMOD,1))' / Mg_fly / (l_wing*1e3);
    My_intactwing_strokeMOD(:,i) = (nansum(My_intact_strokeMOD,1))' / Mg_fly / (l_wing*1e3);
    Mz_intactwing_strokeMOD(:,i) = (nansum(Mz_intact_strokeMOD,1))' / Mg_fly / (l_wing*1e3);

    % total force & torque
    Fx_damaged_total_strokeMOD(:,i) = (nansum(Fx_damaged_strokeMOD,1) + nansum(Fx_intact_strokeMOD,1))' / Mg_fly;
    Fy_damaged_total_strokeMOD(:,i) = (nansum(Fy_damaged_strokeMOD,1) + nansum(Fy_intact_strokeMOD,1))' / Mg_fly;
    Fz_damaged_total_strokeMOD(:,i) = (nansum(Fz_damaged_strokeMOD,1) + nansum(Fz_intact_strokeMOD,1))' / Mg_fly;
    
    Mx_damaged_total_strokeMOD(:,i) = (nansum(Mx_damaged_strokeMOD,1) + nansum(Mx_intact_strokeMOD,1))' / Mg_fly / (l_wing*1e3);
    My_damaged_total_strokeMOD(:,i) = (nansum(My_damaged_strokeMOD,1) + nansum(My_intact_strokeMOD,1))' / Mg_fly / (l_wing*1e3);
    Mz_damaged_total_strokeMOD(:,i) = (nansum(Mz_damaged_strokeMOD,1) + nansum(Mz_intact_strokeMOD,1))' / Mg_fly / (l_wing*1e3);
    
    Fx_NONdamaged_total_strokeMOD(:,i) = (nansum(Fx_NONdamaged_strokeMOD,1) + nansum(Fx_intact_strokeMOD,1))' / Mg_fly;
    Fy_NONdamaged_total_strokeMOD(:,i) = (nansum(Fy_NONdamaged_strokeMOD,1) + nansum(Fy_intact_strokeMOD,1))' / Mg_fly;
    Fz_NONdamaged_total_strokeMOD(:,i) = (nansum(Fz_NONdamaged_strokeMOD,1) + nansum(Fz_intact_strokeMOD,1))' / Mg_fly;
    
    Mx_NONdamaged_total_strokeMOD(:,i) = (nansum(Mx_NONdamaged_strokeMOD,1) + nansum(Mx_intact_strokeMOD,1))' / Mg_fly / (l_wing*1e3);
    My_NONdamaged_total_strokeMOD(:,i) = (nansum(My_NONdamaged_strokeMOD,1) + nansum(My_intact_strokeMOD,1))' / Mg_fly / (l_wing*1e3);
    Mz_NONdamaged_total_strokeMOD(:,i) = (nansum(Mz_NONdamaged_strokeMOD,1) + nansum(Mz_intact_strokeMOD,1))' / Mg_fly / (l_wing*1e3);
    end
    
    %% WB kin dev MOD
    if dev_on == 1
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
    [ FM_strkpln, kine ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_lift_on );

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
    
    %% calc spanwise wing torque & rotation angle amplitude 
    chords_L        = 1e-3*wing_model.chords_L';      % [ m ]
    x_LE_L          = 1e-3*wing_model.x_LE_L';        % [ m ]
    chords_R        = 1e-3*wing_model.chords_R';      % [ m ]
    x_LE_R          = 1e-3*wing_model.x_LE_R';        % [ m ]

    rot_ddot_L = kine.eta_ddot_L;
    rot_ddot_R = kine.eta_ddot_R;
    
    Maero_sect_L = FM_strkpln.FM_L.My_L;
    Maero_sect_R = FM_strkpln.FM_R.My_R;
    
    clear ax_L ax_R az_L az_R Mmass_sect_L Mmass_sect_R
    for n = 1:length(rot_L)
        
        % torque arm
        x_cm_L = (x_LE_L-chords_L/2)*sin(rot_L(n));
        z_cm_L = (x_LE_L-chords_L/2)*cos(rot_L(n));

        x_cm_R = (x_LE_R-chords_R/2)*sin(rot_R(n));
        z_cm_R = (x_LE_R-chords_R/2)*cos(rot_R(n));

        ax_L(:,:) = kine.Udot_left(1,n,:)*1e-3;      % [ m/s^2 ]
        az_L(:,:) = kine.Udot_left(3,n,:)*1e-3;      % [ m/s^2 ]

        ax_R(:,:) = kine.Udot_right(1,n,:)*1e-3;     % [ m/s^2 ]
        az_R(:,:) = kine.Udot_right(3,n,:)*1e-3;     % [ m/s^2 ]
        
        Mmass_sect_L(:,n) = mass_tot_sect_intact  .* (x_cm_L.*az_L + z_cm_L.*ax_L);
        Mmass_sect_R(:,n) = mass_tot_sect_damaged .* (x_cm_R.*az_R + z_cm_R.*ax_R);
    end

    Mwing_mass_intactwing_devMOD(:,i)  = nansum(Mmass_sect_L)' / Mg_fly / l_wing;
    Mwing_mass_damagedwing_devMOD(:,i) = nansum(Mmass_sect_R)' / Mg_fly / l_wing;
    
    Mwing_inertia_intactwing_devMOD(:,i)  = Itot_intact * rot_ddot_L / Mg_fly / l_wing;
    Mwing_inertia_damagedwing_devMOD(:,i) = Itot_damaged * rot_ddot_R / Mg_fly / l_wing;
    
    Mwing_aero_intactwing_devMOD(:,i)  = nansum(Maero_sect_L)' / Mg_fly / (l_wing*1e3);
    Mwing_aero_damagedwing_devMOD(:,i) = nansum(Maero_sect_R)' / Mg_fly / (l_wing*1e3);
    
    Mwing_total_intactwing_devMOD(:,i)  = Mwing_mass_intactwing_devMOD(:,i) + Mwing_inertia_intactwing_devMOD(:,i) + Mwing_aero_intactwing_devMOD(:,i);
    Mwing_total_damagedwing_devMOD(:,i)  = Mwing_mass_damagedwing_devMOD(:,i) + Mwing_inertia_damagedwing_devMOD(:,i) + Mwing_aero_damagedwing_devMOD(:,i);

    %% torque at min & max rotation angle & rotation amplitude
%     [rot_maxpeaks_L,rot_maxpeaklocs_L] = findpeaks(rot_L);
%     n_rot_max_L = rot_maxpeaklocs_L(maxpeak_nr);
%     rot_max_L = rot_L(n_rot_max_L);
% 
%     [rot_minpeaks_L,rot_minpeaklocs_L] = findpeaks(-rot_L);
%     n_rot_min_L = rot_minpeaklocs_L(minpeak_nr);
%     rot_min_L = rot_L(n_rot_min_L);
% 
%     [rot_maxpeaks_R,rot_maxpeaklocs_R] = findpeaks(rot_R);
%     n_rot_max_R = rot_maxpeaklocs_R(maxpeak_nr);
%     rot_max_R = rot_R(n_rot_max_R);
% 
%     [rot_minpeaks_R,rot_minpeaklocs_R] = findpeaks(-rot_R);
%     n_rot_min_R = rot_minpeaklocs_R(minpeak_nr);
%     rot_min_R = rot_R(n_rot_min_R);
% 
%     % plot
% %     subplot(2,1,1)
% %     plot(t_norm,rot_L)
% %     hold on
% %     plot(t_norm(rot_maxpeaklocs_L),rot_maxpeaks_L,'or')
% %     plot(t_norm(rot_minpeaklocs_L),rot_L(rot_minpeaklocs_L),'og')
% %     
% %     plot(t_norm(n_rot_max_L),rot_max_L,'*r')
% %     plot(t_norm(n_rot_min_L),rot_min_L,'*g')
% % 
% %     subplot(2,1,2)
% %     plot(t_norm,rot_R)
% %     hold on
% %     plot(t_norm(rot_maxpeaklocs_R),rot_maxpeaks_R,'or')
% %     plot(t_norm(rot_minpeaklocs_R),rot_R(rot_minpeaklocs_R),'og')
% %     
% %     plot(t_norm(n_rot_max_R),rot_max_R,'*r')
% %     plot(t_norm(n_rot_min_R),rot_min_R,'*g')
    
    %% store data
%     % rotation angle amplitude etc
%     Arot_intactwing_devMOD(i,1) = rad2deg(rot_max_L - rot_min_L);
%     Arot_damagedwing_devMOD(i,1) = rad2deg(rot_max_R - rot_min_R);
%     
%     dMwing_total_intactwing_devMOD(i,1) = Mwing_total_intactwing_devMOD(n_rot_max_L,i) - Mwing_total_intactwing_devMOD(n_rot_min_L,i);
%     dMwing_total_damagedwing_devMOD(i,1) = Mwing_total_intactwing_devMOD(n_rot_max_R,i) - Mwing_total_intactwing_devMOD(n_rot_min_R,i);

    % total force & torque damaged & intact wing
    Fx_damagedwing_devMOD(:,i) = (nansum(Fx_damaged_devMOD,1))' / Mg_fly;
    Fy_damagedwing_devMOD(:,i) = (nansum(Fy_damaged_devMOD,1))' / Mg_fly;
    Fz_damagedwing_devMOD(:,i) = (nansum(Fz_damaged_devMOD,1))' / Mg_fly;
    
    Mx_damagedwing_devMOD(:,i) = (nansum(Mx_damaged_devMOD,1))' / Mg_fly / (l_wing*1e3);
    My_damagedwing_devMOD(:,i) = (nansum(My_damaged_devMOD,1))' / Mg_fly / (l_wing*1e3);
    Mz_damagedwing_devMOD(:,i) = (nansum(Mz_damaged_devMOD,1))' / Mg_fly / (l_wing*1e3);
    
    Fx_intactwing_devMOD(:,i) = (nansum(Fx_intact_devMOD,1))' / Mg_fly;
    Fy_intactwing_devMOD(:,i) = (nansum(Fy_intact_devMOD,1))' / Mg_fly;
    Fz_intactwing_devMOD(:,i) = (nansum(Fz_intact_devMOD,1))' / Mg_fly;
    
    Mx_intactwing_devMOD(:,i) = (nansum(Mx_intact_devMOD,1))' / Mg_fly / (l_wing*1e3);
    My_intactwing_devMOD(:,i) = (nansum(My_intact_devMOD,1))' / Mg_fly / (l_wing*1e3);
    Mz_intactwing_devMOD(:,i) = (nansum(Mz_intact_devMOD,1))' / Mg_fly / (l_wing*1e3);

    % total force & torque
    Fx_damaged_total_devMOD(:,i) = (nansum(Fx_damaged_devMOD,1) + nansum(Fx_intact_devMOD,1))' / Mg_fly;
    Fy_damaged_total_devMOD(:,i) = (nansum(Fy_damaged_devMOD,1) + nansum(Fy_intact_devMOD,1))' / Mg_fly;
    Fz_damaged_total_devMOD(:,i) = (nansum(Fz_damaged_devMOD,1) + nansum(Fz_intact_devMOD,1))' / Mg_fly;
    
    Mx_damaged_total_devMOD(:,i) = (nansum(Mx_damaged_devMOD,1) + nansum(Mx_intact_devMOD,1))' / Mg_fly / (l_wing*1e3);
    My_damaged_total_devMOD(:,i) = (nansum(My_damaged_devMOD,1) + nansum(My_intact_devMOD,1))' / Mg_fly / (l_wing*1e3);
    Mz_damaged_total_devMOD(:,i) = (nansum(Mz_damaged_devMOD,1) + nansum(Mz_intact_devMOD,1))' / Mg_fly / (l_wing*1e3);
    
    Fx_NONdamaged_total_devMOD(:,i) = (nansum(Fx_NONdamaged_devMOD,1) + nansum(Fx_intact_devMOD,1))' / Mg_fly;
    Fy_NONdamaged_total_devMOD(:,i) = (nansum(Fy_NONdamaged_devMOD,1) + nansum(Fy_intact_devMOD,1))' / Mg_fly;
    Fz_NONdamaged_total_devMOD(:,i) = (nansum(Fz_NONdamaged_devMOD,1) + nansum(Fz_intact_devMOD,1))' / Mg_fly;
    
    Mx_NONdamaged_total_devMOD(:,i) = (nansum(Mx_NONdamaged_devMOD,1) + nansum(Mx_intact_devMOD,1))' / Mg_fly / (l_wing*1e3);
    My_NONdamaged_total_devMOD(:,i) = (nansum(My_NONdamaged_devMOD,1) + nansum(My_intact_devMOD,1))' / Mg_fly / (l_wing*1e3);
    Mz_NONdamaged_total_devMOD(:,i) = (nansum(Mz_NONdamaged_devMOD,1) + nansum(Mz_intact_devMOD,1))' / Mg_fly / (l_wing*1e3);
    end
    
    %% WB kin rot MOD
    if wingrot_on == 1
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
    [ FM_strkpln, kine ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_lift_on );

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
    
    %% calc spanwise wing torque & rotation angle amplitude 
    chords_L        = 1e-3*wing_model.chords_L';      % [ m ]
    x_LE_L          = 1e-3*wing_model.x_LE_L';        % [ m ]
    chords_R        = 1e-3*wing_model.chords_R';      % [ m ]
    x_LE_R          = 1e-3*wing_model.x_LE_R';        % [ m ]

    rot_ddot_L = kine.eta_ddot_L;
    rot_ddot_R = kine.eta_ddot_R;
    
    Maero_sect_L = FM_strkpln.FM_L.My_L;
    Maero_sect_R = FM_strkpln.FM_R.My_R;
    
    clear ax_L ax_R az_L az_R Mmass_sect_L Mmass_sect_R
    for n = 1:length(rot_L)
        
        % torque arm
        x_cm_L = (x_LE_L-chords_L/2)*sin(rot_L(n));
        z_cm_L = (x_LE_L-chords_L/2)*cos(rot_L(n));

        x_cm_R = (x_LE_R-chords_R/2)*sin(rot_R(n));
        z_cm_R = (x_LE_R-chords_R/2)*cos(rot_R(n));

        ax_L(:,:) = kine.Udot_left(1,n,:)*1e-3;      % [ m/s^2 ]
        az_L(:,:) = kine.Udot_left(3,n,:)*1e-3;      % [ m/s^2 ]

        ax_R(:,:) = kine.Udot_right(1,n,:)*1e-3;     % [ m/s^2 ]
        az_R(:,:) = kine.Udot_right(3,n,:)*1e-3;     % [ m/s^2 ]
        
        Mmass_sect_L(:,n) = mass_tot_sect_intact  .* (x_cm_L.*az_L + z_cm_L.*ax_L);
        Mmass_sect_R(:,n) = mass_tot_sect_damaged .* (x_cm_R.*az_R + z_cm_R.*ax_R);
    end

    Mwing_mass_intactwing_rotMOD(:,i)  = nansum(Mmass_sect_L)' / Mg_fly / l_wing;
    Mwing_mass_damagedwing_rotMOD(:,i) = nansum(Mmass_sect_R)' / Mg_fly / l_wing;
    
    Mwing_inertia_intactwing_rotMOD(:,i)  = Itot_intact * rot_ddot_L / Mg_fly / l_wing;
    Mwing_inertia_damagedwing_rotMOD(:,i) = Itot_damaged * rot_ddot_R / Mg_fly / l_wing;
    
    Mwing_aero_intactwing_rotMOD(:,i)  = nansum(Maero_sect_L)' / Mg_fly / (l_wing*1e3);
    Mwing_aero_damagedwing_rotMOD(:,i) = nansum(Maero_sect_R)' / Mg_fly / (l_wing*1e3);
    
    Mwing_total_intactwing_rotMOD(:,i)  = Mwing_mass_intactwing_rotMOD(:,i) + Mwing_inertia_intactwing_rotMOD(:,i) + Mwing_aero_intactwing_rotMOD(:,i);
    Mwing_total_damagedwing_rotMOD(:,i)  = Mwing_mass_damagedwing_rotMOD(:,i) + Mwing_inertia_damagedwing_rotMOD(:,i) + Mwing_aero_damagedwing_rotMOD(:,i);

    %% torque at min & max rotation angle & rotation amplitude
%     [rot_maxpeaks_L,rot_maxpeaklocs_L] = findpeaks(rot_L);
%     n_rot_max_L = rot_maxpeaklocs_L(maxpeak_nr);
%     rot_max_L = rot_L(n_rot_max_L);
% 
%     [rot_minpeaks_L,rot_minpeaklocs_L] = findpeaks(-rot_L);
%     n_rot_min_L = rot_minpeaklocs_L(minpeak_nr);
%     rot_min_L = rot_L(n_rot_min_L);
% 
%     [rot_maxpeaks_R,rot_maxpeaklocs_R] = findpeaks(rot_R);
%     n_rot_max_R = rot_maxpeaklocs_R(maxpeak_nr);
%     rot_max_R = rot_R(n_rot_max_R);
% 
%     [rot_minpeaks_R,rot_minpeaklocs_R] = findpeaks(-rot_R);
%     n_rot_min_R = rot_minpeaklocs_R(minpeak_nr);
%     rot_min_R = rot_R(n_rot_min_R);
% 
%     % plot
% %     subplot(2,1,1)
% %     plot(t_norm,rot_L)
% %     hold on
% %     plot(t_norm(rot_maxpeaklocs_L),rot_maxpeaks_L,'or')
% %     plot(t_norm(rot_minpeaklocs_L),rot_L(rot_minpeaklocs_L),'og')
% %     
% %     plot(t_norm(n_rot_max_L),rot_max_L,'*r')
% %     plot(t_norm(n_rot_min_L),rot_min_L,'*g')
% % 
% %     subplot(2,1,2)
% %     plot(t_norm,rot_R)
% %     hold on
% %     plot(t_norm(rot_maxpeaklocs_R),rot_maxpeaks_R,'or')
% %     plot(t_norm(rot_minpeaklocs_R),rot_R(rot_minpeaklocs_R),'og')
% %     
% %     plot(t_norm(n_rot_max_R),rot_max_R,'*r')
% %     plot(t_norm(n_rot_min_R),rot_min_R,'*g')
%     
    %% store data
%     % rotation angle amplitude etc
%     Arot_intactwing_rotMOD(i,1) = rad2deg(rot_max_L - rot_min_L);
%     Arot_damagedwing_rotMOD(i,1) = rad2deg(rot_max_R - rot_min_R);
%     
%     dMwing_total_intactwing_rotMOD(i,1) = Mwing_total_intactwing_rotMOD(n_rot_max_L,i) - Mwing_total_intactwing_rotMOD(n_rot_min_L,i);
%     dMwing_total_damagedwing_rotMOD(i,1) = Mwing_total_intactwing_rotMOD(n_rot_max_R,i) - Mwing_total_intactwing_rotMOD(n_rot_min_R,i);

    % total force & torque damaged & intact wing
    Fx_damagedwing_rotMOD(:,i) = (nansum(Fx_damaged_rotMOD,1))' / Mg_fly;
    Fy_damagedwing_rotMOD(:,i) = (nansum(Fy_damaged_rotMOD,1))' / Mg_fly;
    Fz_damagedwing_rotMOD(:,i) = (nansum(Fz_damaged_rotMOD,1))' / Mg_fly;
    
    Mx_damagedwing_rotMOD(:,i) = (nansum(Mx_damaged_rotMOD,1))' / Mg_fly / (l_wing*1e3);
    My_damagedwing_rotMOD(:,i) = (nansum(My_damaged_rotMOD,1))' / Mg_fly / (l_wing*1e3);
    Mz_damagedwing_rotMOD(:,i) = (nansum(Mz_damaged_rotMOD,1))' / Mg_fly / (l_wing*1e3);
    
    Fx_intactwing_rotMOD(:,i) = (nansum(Fx_intact_rotMOD,1))' / Mg_fly;
    Fy_intactwing_rotMOD(:,i) = (nansum(Fy_intact_rotMOD,1))' / Mg_fly;
    Fz_intactwing_rotMOD(:,i) = (nansum(Fz_intact_rotMOD,1))' / Mg_fly;
    
    Mx_intactwing_rotMOD(:,i) = (nansum(Mx_intact_rotMOD,1))' / Mg_fly / (l_wing*1e3);
    My_intactwing_rotMOD(:,i) = (nansum(My_intact_rotMOD,1))' / Mg_fly / (l_wing*1e3);
    Mz_intactwing_rotMOD(:,i) = (nansum(Mz_intact_rotMOD,1))' / Mg_fly / (l_wing*1e3);

    % total force & torque
    Fx_damaged_total_rotMOD(:,i) = (nansum(Fx_damaged_rotMOD,1) + nansum(Fx_intact_rotMOD,1))' / Mg_fly;
    Fy_damaged_total_rotMOD(:,i) = (nansum(Fy_damaged_rotMOD,1) + nansum(Fy_intact_rotMOD,1))' / Mg_fly;
    Fz_damaged_total_rotMOD(:,i) = (nansum(Fz_damaged_rotMOD,1) + nansum(Fz_intact_rotMOD,1))' / Mg_fly;
    
    Mx_damaged_total_rotMOD(:,i) = (nansum(Mx_damaged_rotMOD,1) + nansum(Mx_intact_rotMOD,1))' / Mg_fly / (l_wing*1e3);
    My_damaged_total_rotMOD(:,i) = (nansum(My_damaged_rotMOD,1) + nansum(My_intact_rotMOD,1))' / Mg_fly / (l_wing*1e3);
    Mz_damaged_total_rotMOD(:,i) = (nansum(Mz_damaged_rotMOD,1) + nansum(Mz_intact_rotMOD,1))' / Mg_fly / (l_wing*1e3);
    
    Fx_NONdamaged_total_rotMOD(:,i) = (nansum(Fx_NONdamaged_rotMOD,1) + nansum(Fx_intact_rotMOD,1))' / Mg_fly;
    Fy_NONdamaged_total_rotMOD(:,i) = (nansum(Fy_NONdamaged_rotMOD,1) + nansum(Fy_intact_rotMOD,1))' / Mg_fly;
    Fz_NONdamaged_total_rotMOD(:,i) = (nansum(Fz_NONdamaged_rotMOD,1) + nansum(Fz_intact_rotMOD,1))' / Mg_fly;
    
    Mx_NONdamaged_total_rotMOD(:,i) = (nansum(Mx_NONdamaged_rotMOD,1) + nansum(Mx_intact_rotMOD,1))' / Mg_fly / (l_wing*1e3);
    My_NONdamaged_total_rotMOD(:,i) = (nansum(My_NONdamaged_rotMOD,1) + nansum(My_intact_rotMOD,1))' / Mg_fly / (l_wing*1e3);
    Mz_NONdamaged_total_rotMOD(:,i) = (nansum(Mz_NONdamaged_rotMOD,1) + nansum(Mz_intact_rotMOD,1))' / Mg_fly / (l_wing*1e3);
    end
end

if all_means_on == 1
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
else
    
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
end

%% change min-max into max-min
if loc_set == 3
%     Arot_intactwing_steady = -Arot_intactwing_steady;
%     Arot_damagedwing_steady = -Arot_damagedwing_steady;
%     
%     dMwing_total_intactwing_steady = -dMwing_total_intactwing_steady;
%     dMwing_total_damagedwing_steady = -dMwing_total_damagedwing_steady;
    
    Arot_intactwing_all = -Arot_intactwing_all;
    Arot_damagedwing_all = -Arot_damagedwing_all;
    
    dMwing_total_intactwing_all = -dMwing_total_intactwing_all;
    dMwing_total_damagedwing_all = -dMwing_total_damagedwing_all;
    
%     Arot_intactwing_freqMOD = -Arot_intactwing_freqMOD;
%     Arot_damagedwing_freqMOD = -Arot_damagedwing_freqMOD;
%     
%     dMwing_total_intactwing_freqMOD = -dMwing_total_intactwing_freqMOD;
%     dMwing_total_damagedwing_freqMOD = -dMwing_total_damagedwing_freqMOD;
%     
%     Arot_intactwing_strokeMOD = -Arot_intactwing_strokeMOD;
%     Arot_damagedwing_strokeMOD = -Arot_damagedwing_strokeMOD;
%     
%     dMwing_total_intactwing_strokeMOD = -dMwing_total_intactwing_strokeMOD;
%     dMwing_total_damagedwing_strokeMOD = -dMwing_total_damagedwing_strokeMOD;
%     
%     Arot_intactwing_devMOD = -Arot_intactwing_devMOD;
%     Arot_damagedwing_devMOD = -Arot_damagedwing_devMOD;
%     
%     dMwing_total_intactwing_devMOD = -dMwing_total_intactwing_devMOD;
%     dMwing_total_damagedwing_devMOD = -dMwing_total_damagedwing_devMOD;
%     
%     Arot_intactwing_rotMOD = -Arot_intactwing_rotMOD;
%     Arot_damagedwing_rotMOD = -Arot_damagedwing_rotMOD;
%     
%     dMwing_total_intactwing_rotMOD = -dMwing_total_intactwing_rotMOD;
%     dMwing_total_damagedwing_rotMOD = -dMwing_total_damagedwing_rotMOD;
    
end

%% save data & plot
    save(['allMODs_TEclip_freqAsym',num2str(freq_asymFitNr),'_peakloc',num2str(loc_set),'.mat'])
    
    mkdir('qsModel_FnM_TEcut')
    cd('qsModel_FnM_TEcut')
    
    % rotation angle throughout wingbeat
    subplot(2,1,1)
    title('intact wing')
    
    subplot(2,1,2)
    title('damaged wing')
    
    saveas(gca,['rotationangles_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.fig'])
    saveas(gca,['rotationangles_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.png'])
    plot2svg(['rotationangles_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.svg'])
    
%% wing rotation spring system test
figure
plot(Arot_intactwing_all,dMwing_total_intactwing_all,'k-s','markerfacecolor','m','markersize',10)
hold on
plot(Arot_damagedwing_all,dMwing_total_damagedwing_all,'k-d','markerfacecolor','r','markersize',10)
xlabel('wing rotation amplitude')
ylabel('spanwise wing torque amplitude')
legend('intact wing','damaged wing')

saveas(gca,['Arot_vs_DMwing_TIPclip_asympFit',num2str(freq_asymFitNr),'_peakloc',num2str(loc_set),'.fig'])
saveas(gca,['Arot_vs_DMwing_TIPclip_asympFit',num2str(freq_asymFitNr),'_peakloc',num2str(loc_set),'.png'])
plot2svg(['Arot_vs_DMwing_TIPclip_asympFit',num2str(freq_asymFitNr),'_peakloc',num2str(loc_set),'.svg'])
    
%% forces & torques plots
if all_means_on == 1
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
    
    saveas(gca,['FnM_WBmod_TEclip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['FnM_WBmod_TEclip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['FnM_WBmod_TEclip_asympFit',num2str(freq_asymFitNr),'.svg'])

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

    saveas(gca,['Fz_Mx_WBmodComponents_TEclip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['Fz_Mx_WBmodComponents_TEclip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['Fz_Mx_WBmodComponents_TEclip_asympFit',num2str(freq_asymFitNr),'.svg'])
    
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

    saveas(gca,['Fy_Mz_WBmodComponents_TEclip_asympFit',num2str(freq_asymFitNr),'.fig'])
    saveas(gca,['Fy_Mz_WBmodComponents_TEclip_asympFit',num2str(freq_asymFitNr),'.png'])
    plot2svg(['Fy_Mz_WBmodComponents_TEclip_asympFit',num2str(freq_asymFitNr),'.svg'])
end
    cd ..
    end