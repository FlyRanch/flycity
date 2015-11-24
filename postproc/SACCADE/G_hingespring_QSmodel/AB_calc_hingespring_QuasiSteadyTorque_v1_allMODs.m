clear;
clc;
close all
warning off

load('bodyNwingModel_4qsModel.mat')

loadname_steady=dir('WBdataset_steady_1603WBs.mat')
loadname_steady = loadname_steady.name;
load(loadname_steady)

loadname=dir('WBmod_torquebased_YawTorque_457WBs.mat')
loadname = loadname.name;
load(loadname)

plot_on = 1
% plot_on = 0

% qs model INC rotational lift
rot_on=1;

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

%% wing kin MODs
YawMods_min = 0;
YawMods_max = 1.5;
dYawMods = .1;

YawMods = [YawMods_min:dYawMods:YawMods_max]';

%% constants
nr_sect = settings.nr_chord_sect;
nr_timepoints = settings.nr_timepoints;

nr_timepoints = 500;
nr_timepoints = 5000;

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

%% rwd fly WB kin MODs
f_YawTorque_fly = f_wb_YawTorque_meanCIstd(1);
freqMOD = freqMOD_wb_YawTorque_meanCIstd(1);

% fouriers coeffs
strokeMOD_coeffs_fwd = strokeMOD_fwd_YawTorque_fourier_coeffs_binmean;
devMOD_coeffs_fwd = devMOD_fwd_YawTorque_fourier_coeffs_binmean;
rotMOD_coeffs_fwd = pitchMOD_fwd_YawTorque_fourier_coeffs_binmean;

[strokeMOD_fwd] = (calc_val_fourier_series_4thN8th_order(t_norm,strokeMOD_coeffs_fwd,0))';
[devMOD_fwd] = (calc_val_fourier_series_4thN8th_order(t_norm,devMOD_coeffs_fwd,0))';
[rotMOD_fwd] = (calc_val_fourier_series_4thN8th_order(t_norm,rotMOD_coeffs_fwd,0))';

strokeMOD_coeffs_rwd = strokeMOD_rwd_YawTorque_fourier_coeffs_binmean;
devMOD_coeffs_rwd = devMOD_rwd_YawTorque_fourier_coeffs_binmean;
rotMOD_coeffs_rwd = pitchMOD_rwd_YawTorque_fourier_coeffs_binmean;

[strokeMOD_rwd] = (calc_val_fourier_series_4thN8th_order(t_norm,strokeMOD_coeffs_rwd,0))';
[devMOD_rwd] = (calc_val_fourier_series_4thN8th_order(t_norm,devMOD_coeffs_rwd,0))';
[rotMOD_rwd] = (calc_val_fourier_series_4thN8th_order(t_norm,rotMOD_coeffs_rwd,0))';

%% mass & spanwise inertial moment of wings
rho_air = 1.225;                            % [kg/m^3]

h_sect = wing_model.length/nr_sect*1e-3;    % [m]
x_sect = wing_model.x_LE_R'*1e-3;           % [m]
y_sect = wing_model.y_sect_R(:,2)*1e-3;     % [m]
c_sect = wing_model.chords_R'*1e-3;         % [m]

c_max = max(c_sect);

% wing section mass
mass_wing_sect = massperArea * h_sect * c_sect;
mass_virtual_sect = 1/4* pi * rho_air * h_sect * c_sect.^2;
mass_tot_sect = mass_wing_sect + mass_virtual_sect;

% Iwing & Iaddedmass
Iwing_sect = 1/12 * massperArea * h_sect * c_sect.^3 + massperArea * h_sect * c_sect .* (1/2*c_sect - x_sect).^2;
Iaddedmass_sect = 1/64 * pi * rho_air * h_sect * c_sect.^4 + 1/4* pi * rho_air * h_sect * c_sect.^2 .* (1/2*c_sect - x_sect).^2;
Itot_sect = Iwing_sect + Iaddedmass_sect;

Iwing = nansum(Iwing_sect);
Iaddedmass = nansum(Iaddedmass_sect);
Itot = nansum(Itot_sect);

%% loop with different mods
for i = 1:length(YawMods)

    YawMod_now = YawMods(i)
    YawTorque_now = YawMod_now * Myaw_norm;
        
    freq_now            = YawMod_now * freqMOD + f_steady;    
    
    stroke_fwd_now   = YawMod_now * strokeMOD_fwd + stroke_steady;    
    dev_fwd_now      = YawMod_now * devMOD_fwd    + dev_steady;    
    rot_fwd_now      = YawMod_now * rotMOD_fwd    + rot_steady;    

    stroke_rwd_now   = YawMod_now * strokeMOD_rwd + stroke_steady;    
    dev_rwd_now      = YawMod_now * devMOD_rwd    + dev_steady;    
    rot_rwd_now      = YawMod_now * rotMOD_rwd    + rot_steady;    
    
    %% WB kin ALL MODs
    freq = freq_now;

    stroke_L = deg2rad(stroke_fwd_now);
    dev_L = deg2rad(dev_fwd_now);
    rot_L = deg2rad(rot_fwd_now);

    stroke_R = deg2rad(stroke_rwd_now);
    dev_R = deg2rad(dev_rwd_now);
    rot_R = deg2rad(rot_rwd_now);

    % qs forces & torques
    [ FM_strkpln, kine ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev_L, rot_L, stroke_L, dev_R, rot_R, stroke_R, freq, rot_on );

    Fx_L_all = FM_strkpln.FM_L.Fx_strkpln_L;
    Fy_L_all = FM_strkpln.FM_L.Fy_strkpln_L;
    Fz_L_all = FM_strkpln.FM_L.Fz_strkpln_L;
    
    Mx_L_all = FM_strkpln.FM_L.Mx_strkpln_L;
    My_L_all = FM_strkpln.FM_L.My_strkpln_L;
    Mz_L_all = FM_strkpln.FM_L.Mz_strkpln_L;
    
    Fx_R_all = FM_strkpln.FM_R.Fx_strkpln_R;
    Fy_R_all = FM_strkpln.FM_R.Fy_strkpln_R;
    Fz_R_all = FM_strkpln.FM_R.Fz_strkpln_R;
    
    Mx_R_all = FM_strkpln.FM_R.Mx_strkpln_R;
    My_R_all = FM_strkpln.FM_R.My_strkpln_R;
    Mz_R_all = FM_strkpln.FM_R.Mz_strkpln_R;
    
    %% calc spanwise wing torque & rotation angle amplitude 

    rot_ddot_L = kine.eta_ddot_L;
    rot_ddot_R = kine.eta_ddot_R;
    
    Maero_sect_L = FM_strkpln.FM_L.My_L;
    Maero_sect_R = FM_strkpln.FM_R.My_R;
    
    clear ax_L ax_R az_L az_R Mmass_sect_L Mmass_sect_R
    for n = 1:length(rot_L)
        
        % torque arm
        x_cm_L = (x_sect-c_sect/2)*sin(rot_L(n));
        z_cm_L = (x_sect-c_sect/2)*cos(rot_L(n));

        x_cm_R = (x_sect-c_sect/2)*sin(rot_R(n));
        z_cm_R = (x_sect-c_sect/2)*cos(rot_R(n));

        ax_L(:,:) = kine.Udot_left(1,n,:)*1e-3;      % [ m/s^2 ]
        az_L(:,:) = kine.Udot_left(3,n,:)*1e-3;      % [ m/s^2 ]

        ax_R(:,:) = kine.Udot_right(1,n,:)*1e-3;     % [ m/s^2 ]
        az_R(:,:) = kine.Udot_right(3,n,:)*1e-3;     % [ m/s^2 ]
        
        Mmass_sect_L(:,n) = mass_tot_sect  .* (x_cm_L.*az_L + z_cm_L.*ax_L);
        Mmass_sect_R(:,n) = mass_tot_sect .* (x_cm_R.*az_R + z_cm_R.*ax_R);
    end

    Mwing_mass_fwd_all(:,i)  = nansum(Mmass_sect_L)' / Mg_fly / l_wing;
    Mwing_mass_rwd_all(:,i) = nansum(Mmass_sect_R)' / Mg_fly / l_wing;
    
    Mwing_inertia_fwd_all(:,i)  = Itot * rot_ddot_L / Mg_fly / l_wing;
    Mwing_inertia_rwd_all(:,i) = Itot * rot_ddot_R / Mg_fly / l_wing;
    
    Mwing_aero_fwd_all(:,i)  = nansum(Maero_sect_L)' / Mg_fly / (l_wing*1e3);
    Mwing_aero_rwd_all(:,i) = nansum(Maero_sect_R)' / Mg_fly / (l_wing*1e3);
    
    Mwing_total_fwd_all(:,i)  = Mwing_mass_fwd_all(:,i) + Mwing_inertia_fwd_all(:,i) + Mwing_aero_fwd_all(:,i);
    Mwing_total_rwd_all(:,i)  = Mwing_mass_rwd_all(:,i) + Mwing_inertia_rwd_all(:,i) + Mwing_aero_rwd_all(:,i);
    
    %% torque at min & max rotation angle
    [rot_maxpeaks_L,rot_maxpeaklocs_L] = findpeaks(rot_L);
    n_rot_max_L = rot_maxpeaklocs_L(maxpeak_nr);
    rot_max_L = rot_L(n_rot_max_L);

    [rot_minpeaks_L,rot_minpeaklocs_L] = findpeaks(-rot_L);
    n_rot_min_L = rot_minpeaklocs_L(minpeak_nr);
    rot_min_L = rot_L(n_rot_min_L);

    [rot_maxpeaks_R,rot_maxpeaklocs_R] = findpeaks(rot_R);
    n_rot_max_R = rot_maxpeaklocs_R(maxpeak_nr);
    rot_max_R = rot_R(n_rot_max_R);

    [rot_minpeaks_R,rot_minpeaklocs_R] = findpeaks(-rot_R);
    n_rot_min_R = rot_minpeaklocs_R(minpeak_nr);
    rot_min_R = rot_R(n_rot_min_R);

    % plot
    subplot(2,1,1)
    plot(t_norm,rot_L)
    hold on
    plot(t_norm(rot_maxpeaklocs_L),rot_maxpeaks_L,'or')
    plot(t_norm(rot_minpeaklocs_L),rot_L(rot_minpeaklocs_L),'og')
    
    plot(t_norm(n_rot_max_L),rot_max_L,'*r')
    plot(t_norm(n_rot_min_L),rot_min_L,'*g')

    subplot(2,1,2)
    plot(t_norm,rot_R)
    hold on
    plot(t_norm(rot_maxpeaklocs_R),rot_maxpeaks_R,'or')
    plot(t_norm(rot_minpeaklocs_R),rot_R(rot_minpeaklocs_R),'og')
    
    plot(t_norm(n_rot_max_R),rot_max_R,'*r')
    plot(t_norm(n_rot_min_R),rot_min_R,'*g')

    %% store data
    % rotation angle amplitude etc
    Arot_fwd_all(i,1) = rad2deg(rot_max_L - rot_min_L);
    Arot_rwd_all(i,1) = rad2deg(rot_max_R - rot_min_R);

    n_rot_max_fwd_all(i,:) = n_rot_max_L;
    n_rot_max_rwd_all(i,:) = n_rot_max_R;
    
    n_rot_min_fwd_all(i,:) = n_rot_min_L;
    n_rot_min_rwd_all(i,:) = n_rot_min_R;
    
    rot_max_fwd_all(i,:) = rot_max_L;
    rot_max_rwd_all(i,:) = rot_max_R;
    
    rot_min_fwd_all(i,:) = rot_min_L;
    rot_min_rwd_all(i,:) = rot_min_R;
    
    % spanwise torques @ rot max & rot min
    Mwing_total_fwd_atRotmax(i,1) = Mwing_total_fwd_all(n_rot_max_L,i);
    Mwing_total_fwd_atRotmin(i,1) = Mwing_total_fwd_all(n_rot_min_L,i);
    
    Mwing_total_rwd_atRotmax(i,1) = Mwing_total_rwd_all(n_rot_max_R,i);
    Mwing_total_rwd_atRotmin(i,1) = Mwing_total_rwd_all(n_rot_min_R,i);

    Mwing_mass_fwd_atRotmax(i,1) = Mwing_mass_fwd_all(n_rot_max_L,i);
    Mwing_mass_fwd_atRotmin(i,1) = Mwing_mass_fwd_all(n_rot_min_L,i);
    
    Mwing_mass_rwd_atRotmax(i,1) = Mwing_mass_rwd_all(n_rot_max_R,i);
    Mwing_mass_rwd_atRotmin(i,1) = Mwing_mass_rwd_all(n_rot_min_R,i);

    Mwing_inertia_fwd_atRotmax(i,1) = Mwing_inertia_fwd_all(n_rot_max_L,i);
    Mwing_inertia_fwd_atRotmin(i,1) = Mwing_inertia_fwd_all(n_rot_min_L,i);
    
    Mwing_inertia_rwd_atRotmax(i,1) = Mwing_inertia_rwd_all(n_rot_max_R,i);
    Mwing_inertia_rwd_atRotmin(i,1) = Mwing_inertia_rwd_all(n_rot_min_R,i);

    Mwing_aero_fwd_atRotmax(i,1) = Mwing_aero_fwd_all(n_rot_max_L,i);
    Mwing_aero_fwd_atRotmin(i,1) = Mwing_aero_fwd_all(n_rot_min_L,i);
    
    Mwing_aero_rwd_atRotmax(i,1) = Mwing_aero_rwd_all(n_rot_max_R,i);
    Mwing_aero_rwd_atRotmin(i,1) = Mwing_aero_rwd_all(n_rot_min_R,i);

    % spanwise torque change
    dMwing_total_fwd_all(i,1) = (Mwing_total_fwd_all(n_rot_max_L,i) - Mwing_total_fwd_all(n_rot_min_L,i));
    dMwing_total_rwd_all(i,1) = (Mwing_total_rwd_all(n_rot_max_R,i) - Mwing_total_rwd_all(n_rot_min_R,i));

    dMwing_mass_fwd_all(i,1) = (Mwing_mass_fwd_all(n_rot_max_L,i) - Mwing_mass_fwd_all(n_rot_min_L,i));
    dMwing_mass_rwd_all(i,1) = (Mwing_mass_rwd_all(n_rot_max_R,i) - Mwing_mass_rwd_all(n_rot_min_R,i));

    dMwing_inertia_fwd_all(i,1) = (Mwing_inertia_fwd_all(n_rot_max_L,i) - Mwing_inertia_fwd_all(n_rot_min_L,i));
    dMwing_inertia_rwd_all(i,1) = (Mwing_inertia_rwd_all(n_rot_max_R,i) - Mwing_inertia_rwd_all(n_rot_min_R,i));

    dMwing_aero_fwd_all(i,1) = (Mwing_aero_fwd_all(n_rot_max_L,i) - Mwing_aero_fwd_all(n_rot_min_L,i));
    dMwing_aero_rwd_all(i,1) = (Mwing_aero_rwd_all(n_rot_max_R,i) - Mwing_aero_rwd_all(n_rot_min_R,i));
    
    % wingbeat kinematics
    freq_all(i,1) = freq_now;

    stroke_fwd_all(:,i) = (stroke_fwd_now)';
    dev_fwd_all(:,i) = (dev_fwd_now)';
    rot_fwd_all(:,i) = (rot_fwd_now)';

    stroke_rwd_all(:,i) = (stroke_rwd_now)';
    dev_rwd_all(:,i) = (dev_rwd_now)';
    rot_rwd_all(:,i) = (rot_rwd_now)';
    
    % body forces & torques
    
    % input yaw torque
    YawTorques_all(i,1) = YawTorque_now;
    
    % total force & torque rwd & fwd wing
    Fx_rwd_all(:,i) = (nansum(Fx_R_all,1))' / Mg_fly;
    Fy_rwd_all(:,i) = (nansum(Fy_R_all,1))' / Mg_fly;
    Fz_rwd_all(:,i) = (nansum(Fz_R_all,1))' / Mg_fly;
    
    Mx_rwd_all(:,i) = (nansum(Mx_R_all,1))' / Mg_fly / (l_wing*1e3);
    My_rwd_all(:,i) = (nansum(My_R_all,1))' / Mg_fly / (l_wing*1e3);
    Mz_rwd_all(:,i) = (nansum(Mz_R_all,1))' / Mg_fly / (l_wing*1e3);
    
    Fx_fwd_all(:,i) = (nansum(Fx_L_all,1))' / Mg_fly;
    Fy_fwd_all(:,i) = (nansum(Fy_L_all,1))' / Mg_fly;
    Fz_fwd_all(:,i) = (nansum(Fz_L_all,1))' / Mg_fly;
    
    Mx_fwd_all(:,i) = (nansum(Mx_L_all,1))' / Mg_fly / (l_wing*1e3);
    My_fwd_all(:,i) = (nansum(My_L_all,1))' / Mg_fly / (l_wing*1e3);
    Mz_fwd_all(:,i) = (nansum(Mz_L_all,1))' / Mg_fly / (l_wing*1e3);

    % total force & torque
    Fx_total_all(:,i) = (nansum(Fx_R_all,1) + nansum(Fx_L_all,1))' / Mg_fly;
    Fy_total_all(:,i) = (nansum(Fy_R_all,1) + nansum(Fy_L_all,1))' / Mg_fly;
    Fz_total_all(:,i) = (nansum(Fz_R_all,1) + nansum(Fz_L_all,1))' / Mg_fly;
    
    Mx_total_all(:,i) = (nansum(Mx_R_all,1) + nansum(Mx_L_all,1))' / Mg_fly / (l_wing*1e3);
    My_total_all(:,i) = (nansum(My_R_all,1) + nansum(My_L_all,1))' / Mg_fly / (l_wing*1e3);
    Mz_total_all(:,i) = (nansum(Mz_R_all,1) + nansum(Mz_L_all,1))' / Mg_fly / (l_wing*1e3);
    
end

    %% wingbeat average forces & torques: ALL MODs
    Fx_mean_all = nanmean(Fx_total_all);
    Fy_mean_all = nanmean(Fy_total_all);
    Fz_mean_all = nanmean(Fz_total_all);
    
    Mx_mean_all = nanmean(Mx_total_all);
    My_mean_all = nanmean(My_total_all);
    Mz_mean_all = nanmean(Mz_total_all);
    
%% neg Arot & dM
if loc_set == 3
    Arot_fwd_all = -Arot_fwd_all;
    Arot_rwd_all = -Arot_rwd_all;
    
    dMwing_total_fwd_all = -dMwing_total_fwd_all;
    dMwing_total_rwd_all = -dMwing_total_rwd_all;
    
    dMwing_inertia_fwd_all = -dMwing_inertia_fwd_all;
    dMwing_inertia_rwd_all = -dMwing_inertia_rwd_all;
    
    dMwing_mass_fwd_all = -dMwing_mass_fwd_all;
    dMwing_mass_rwd_all = -dMwing_mass_rwd_all;
    
    dMwing_aero_fwd_all = -dMwing_aero_fwd_all;
    dMwing_aero_rwd_all = -dMwing_aero_rwd_all;
end

    %% calc spring constants
    
    %dMtotal
    k_ArotdMtotal_fwd_coeffs_norm_deg = polyfit(Arot_fwd_all,dMwing_total_fwd_all,1);
    k_ArotdMtotal_rwd_coeffs_norm_deg = polyfit(Arot_rwd_all,dMwing_total_rwd_all,1);

    k_ArotdMtotal_fwd_coeffs_norm_rad = polyfit(pi/180*(Arot_fwd_all),dMwing_total_fwd_all,1);
    k_ArotdMtotal_rwd_coeffs_norm_rad = polyfit(pi/180*(Arot_rwd_all),dMwing_total_rwd_all,1);

    k_ArotdMtotal_fwd_norm_deg = k_ArotdMtotal_fwd_coeffs_norm_deg(1);
    k_ArotdMtotal_rwd_norm_deg = k_ArotdMtotal_rwd_coeffs_norm_deg(1);

    k_ArotdMtotal_fwd_norm_rad = k_ArotdMtotal_fwd_coeffs_norm_rad(1);
    k_ArotdMtotal_rwd_norm_rad = k_ArotdMtotal_rwd_coeffs_norm_rad(1);
    
    %dMmass
    k_ArotdMmass_fwd_coeffs_norm_deg = polyfit(Arot_fwd_all,dMwing_mass_fwd_all,1);
    k_ArotdMmass_rwd_coeffs_norm_deg = polyfit(Arot_rwd_all,dMwing_mass_rwd_all,1);

    k_ArotdMmass_fwd_coeffs_norm_rad = polyfit(pi/180*(Arot_fwd_all),dMwing_mass_fwd_all,1);
    k_ArotdMmass_rwd_coeffs_norm_rad = polyfit(pi/180*(Arot_rwd_all),dMwing_mass_rwd_all,1);

    k_ArotdMmass_fwd_norm_deg = k_ArotdMmass_fwd_coeffs_norm_deg(1);
    k_ArotdMmass_rwd_norm_deg = k_ArotdMmass_rwd_coeffs_norm_deg(1);

    k_ArotdMmass_fwd_norm_rad = k_ArotdMmass_fwd_coeffs_norm_rad(1);
    k_ArotdMmass_rwd_norm_rad = k_ArotdMmass_rwd_coeffs_norm_rad(1);
    
    %dMinertia
    k_ArotdMinertia_fwd_coeffs_norm_deg = polyfit(Arot_fwd_all,dMwing_inertia_fwd_all,1);
    k_ArotdMinertia_rwd_coeffs_norm_deg = polyfit(Arot_rwd_all,dMwing_inertia_rwd_all,1);

    k_ArotdMinertia_fwd_coeffs_norm_rad = polyfit(pi/180*(Arot_fwd_all),dMwing_inertia_fwd_all,1);
    k_ArotdMinertia_rwd_coeffs_norm_rad = polyfit(pi/180*(Arot_rwd_all),dMwing_inertia_rwd_all,1);

    k_ArotdMinertia_fwd_norm_deg = k_ArotdMinertia_fwd_coeffs_norm_deg(1);
    k_ArotdMinertia_rwd_norm_deg = k_ArotdMinertia_rwd_coeffs_norm_deg(1);

    k_ArotdMinertia_fwd_norm_rad = k_ArotdMinertia_fwd_coeffs_norm_rad(1);
    k_ArotdMinertia_rwd_norm_rad = k_ArotdMinertia_rwd_coeffs_norm_rad(1);
    
    %dMaero
    k_ArotdMaero_fwd_coeffs_norm_deg = polyfit(Arot_fwd_all,dMwing_aero_fwd_all,1);
    k_ArotdMaero_rwd_coeffs_norm_deg = polyfit(Arot_rwd_all,dMwing_aero_rwd_all,1);

    k_ArotdMaero_fwd_coeffs_norm_rad = polyfit(pi/180*(Arot_fwd_all),dMwing_aero_fwd_all,1);
    k_ArotdMaero_rwd_coeffs_norm_rad = polyfit(pi/180*(Arot_rwd_all),dMwing_aero_rwd_all,1);

    k_ArotdMaero_fwd_norm_deg = k_ArotdMaero_fwd_coeffs_norm_deg(1);
    k_ArotdMaero_rwd_norm_deg = k_ArotdMaero_rwd_coeffs_norm_deg(1);

    k_ArotdMaero_fwd_norm_rad = k_ArotdMaero_fwd_coeffs_norm_rad(1);
    k_ArotdMaero_rwd_norm_rad = k_ArotdMaero_rwd_coeffs_norm_rad(1);
    
%% calc rot0 as function of Tyaw
rot0_max_fwd_all_rad = rot_max_fwd_all + Mwing_total_fwd_atRotmax/k_ArotdMtotal_fwd_norm_rad;
rot0_min_fwd_all_rad = rot_min_fwd_all + Mwing_total_fwd_atRotmin/k_ArotdMtotal_fwd_norm_rad;

rot0_max_rwd_all_rad = rot_max_rwd_all + Mwing_total_rwd_atRotmax/k_ArotdMtotal_rwd_norm_rad;
rot0_min_rwd_all_rad = rot_min_rwd_all + Mwing_total_rwd_atRotmin/k_ArotdMtotal_rwd_norm_rad;

rot_max_fwd_all_deg = rad2deg(rot_max_fwd_all);
rot_min_fwd_all_deg = rad2deg(rot_min_fwd_all);

rot_max_rwd_all_deg = rad2deg(rot_max_rwd_all);
rot_min_rwd_all_deg = rad2deg(rot_min_rwd_all);

rot0_max_fwd_all_deg = rot_max_fwd_all_deg + Mwing_total_fwd_atRotmax/k_ArotdMtotal_fwd_norm_deg;
rot0_min_fwd_all_deg = rot_min_fwd_all_deg + Mwing_total_fwd_atRotmin/k_ArotdMtotal_fwd_norm_deg;

rot0_max_rwd_all_deg = rot_max_rwd_all_deg + Mwing_total_rwd_atRotmax/k_ArotdMtotal_rwd_norm_deg;
rot0_min_rwd_all_deg = rot_min_rwd_all_deg + Mwing_total_rwd_atRotmin/k_ArotdMtotal_rwd_norm_deg;

% trendlines
rot0Tyaw_coeff_fwd_rotmax_deg = polyfit([YawTorques_all],[rot0_max_fwd_all_deg],1);
rot0Tyaw_coeff_rwd_rotmax_deg = polyfit([YawTorques_all],[rot0_max_rwd_all_deg],1);

rot0Tyaw_coeff_fwd_rotmin_deg = polyfit([YawTorques_all],[rot0_min_fwd_all_deg],1);
rot0Tyaw_coeff_rwd_rotmin_deg = polyfit([YawTorques_all],[rot0_min_rwd_all_deg],1);
    
rot0Tyaw_coeff_fwd_rotmax_rad = polyfit([YawTorques_all],[rot0_max_fwd_all_rad],1);
rot0Tyaw_coeff_rwd_rotmax_rad = polyfit([YawTorques_all],[rot0_max_rwd_all_rad],1);

rot0Tyaw_coeff_fwd_rotmin_rad = polyfit([YawTorques_all],[rot0_min_fwd_all_rad],1);
rot0Tyaw_coeff_rwd_rotmin_rad = polyfit([YawTorques_all],[rot0_min_rwd_all_rad],1);
    
    %% save data & plots
    save(['hingespring_Arot_vs_dTspanwise_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.mat'])
    
    mkdir('qsModel_FnM_n_Arot_VS_dTspanwise')
    cd('qsModel_FnM_n_Arot_VS_dTspanwise')
    
    subplot(2,1,1)
    title('forwards moving wing')
    
    subplot(2,1,2)
    title('backwards moving wing')
    
    saveas(gca,['rotationangles_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.fig'])
    saveas(gca,['rotationangles_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.png'])
    plot2svg(['rotationangles_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.svg'])
    
%% plot rotation angles throughout wingbeat
figure
hold on
plot(t_norm,rot_fwd_all(:,1)-90,'-k','linewidth',2)
plot(t_norm,rot_fwd_all(:,end)-90,'-r','linewidth',2)
plot(t_norm,rot_rwd_all(:,end)-90,'-b','linewidth',2)

plot(t_norm(n_rot_max_fwd_all(1)),rot_fwd_all(n_rot_max_fwd_all(1),1)-90,'ok','linewidth',2)
plot(t_norm(n_rot_min_fwd_all(1)),rot_fwd_all(n_rot_min_fwd_all(1),1)-90,'ok','linewidth',2)

plot(t_norm(n_rot_max_fwd_all(end)),rot_fwd_all(n_rot_max_fwd_all(end),end)-90,'or','linewidth',2)
plot(t_norm(n_rot_min_fwd_all(end)),rot_fwd_all(n_rot_min_fwd_all(end),end)-90,'or','linewidth',2)

plot(t_norm(n_rot_max_rwd_all(end)),rot_rwd_all(n_rot_max_rwd_all(end),end)-90,'ob','linewidth',2)
plot(t_norm(n_rot_min_rwd_all(end)),rot_rwd_all(n_rot_min_rwd_all(end),end)-90,'ob','linewidth',2)

legend('steady','fwd','rwd')
axis tight
axis square
axis([0,1,-60,60])
set(gca,'xtick',[0 1])
set(gca,'ytick',[-60 0 60])
xlabel('wingbeat cycle')
ylabel('wing rotation angle [deg]')

saveas(gca,['rot_vs_WBcycle_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.fig'])
saveas(gca,['rot_vs_WBcycle_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.png'])
plot2svg(['rot_vs_WBcycle_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.svg'])

%% plot rotation angle amplitude VS change in spanwise torque
cmap_hot =colormap(hot(100));
cmap_hot = flipud(cmap_hot);
cmap_Tyaw = cmap_hot;

if loc_set == 1
    Arot_min = 128;
    Arot_max = 135;

    dM_min = .45;
    dM_max = .475;

    rot0max_min = 40;
    rot0max_max = 50;

    rot0min_min = 140;
    rot0min_max = 150;
    
elseif loc_set == 2
    Arot_min = 97.5;
    Arot_max = 100;

    dM_min = .138;
    dM_max = .155;

    rot0max_min = 40;
    rot0max_max = 50;

    rot0min_min = 140;
    rot0min_max = 150;
    
elseif loc_set == 3
    Arot_min = 92;
    Arot_max = 97;

    dM_min = -.02;
    dM_max = .06;

    rot0max_min = 40;
    rot0max_max = 50;

    rot0min_min = 140;
    rot0min_max = 150;
end

YawTorque_min = 0;
YawTorque_max = .1;


%% total torque
figure
hold on

    for i = 1:length(YawTorques_all)
        color_nr = round(99/(YawTorque_max-YawTorque_min)*(YawTorques_all(i)-YawTorque_min)+1);
        if color_nr<1
            color_nr=1;
        elseif color_nr>size(cmap_Tyaw,1)
            color_nr=size(cmap_Tyaw,1)
        end

        plot(Arot_fwd_all(i),dMwing_total_fwd_all(i),'dk-','markerfacecolor',cmap_Tyaw(color_nr,:),'markersize',10)
        plot(Arot_rwd_all(i),dMwing_total_rwd_all(i),'o-','markerfacecolor',cmap_Tyaw(color_nr,:),'markersize',10,'color',[.5 .5 .5])
    end

%     plot(Arot_fwd_all,dMwing_total_fwd_all,'-k')
%     plot(Arot_rwd_all,dMwing_total_rwd_all,'-','color',[.5 .5 .5])

    plot([min(Arot_fwd_all) max(Arot_fwd_all)],polyval(k_ArotdMtotal_fwd_coeffs_norm_deg,[min(Arot_fwd_all) max(Arot_fwd_all)]),'-k','linewidth',2)
    plot([min(Arot_rwd_all) max(Arot_rwd_all)],polyval(k_ArotdMtotal_rwd_coeffs_norm_deg,[min(Arot_rwd_all) max(Arot_rwd_all)]),'-k','linewidth',2,'color',[.5 .5 .5])

    legend('fwd','rwd')
    axis tight
    axis square
    axis([Arot_min,Arot_max,dM_min,dM_max])
    set(gca,'xtick',Arot_min:(Arot_max-Arot_min):Arot_max)
    set(gca,'ytick',dM_min:(dM_max-dM_min):dM_max)
    xlabel('rotation angle amplitude [deg]')
    ylabel('normalized spanwise torque change')
    colormap(cmap_Tyaw)
    caxis([YawTorque_min YawTorque_max])
    h = colorbar('location','northoutside'); 
    title(h,'Yaw Torque')
    set(h,'xtick',YawTorque_min:(YawTorque_max-YawTorque_min)/2:YawTorque_max)
    
    saveas(gca,['Arot_vs_dTspanwise_vs_YawTorque_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.fig'])
    saveas(gca,['Arot_vs_dTspanwise_vs_YawTorque_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.png'])
    plot2svg(['Arot_vs_dTspanwise_vs_YawTorque_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.svg'])

%% sep torques
figure
subplot(2,2,1)
hold on
subplot(2,2,2)
hold on
subplot(2,2,3)
hold on
subplot(2,2,4)
hold on

    for i = 1:length(YawTorques_all)
        color_nr = round(99/(YawTorque_max-YawTorque_min)*(YawTorques_all(i)-YawTorque_min)+1);
        if color_nr<1
            color_nr=1;
        elseif color_nr>size(cmap_Tyaw,1)
            color_nr=size(cmap_Tyaw,1)
        end

        subplot(2,2,1)
        plot(Arot_fwd_all(i),dMwing_total_fwd_all(i),'dk-','markerfacecolor',cmap_Tyaw(color_nr,:),'markersize',10)
        plot(Arot_rwd_all(i),dMwing_total_rwd_all(i),'o-','markerfacecolor',cmap_Tyaw(color_nr,:),'markersize',10,'color',[.5 .5 .5])

        subplot(2,2,2)
        plot(Arot_fwd_all(i),dMwing_inertia_fwd_all(i),'dk-','markerfacecolor',cmap_Tyaw(color_nr,:),'markersize',10)
        plot(Arot_rwd_all(i),dMwing_inertia_rwd_all(i),'o-','markerfacecolor',cmap_Tyaw(color_nr,:),'markersize',10,'color',[.5 .5 .5])

        subplot(2,2,3)
        plot(Arot_fwd_all(i),dMwing_mass_fwd_all(i),'dk-','markerfacecolor',cmap_Tyaw(color_nr,:),'markersize',10)
        plot(Arot_rwd_all(i),dMwing_mass_rwd_all(i),'o-','markerfacecolor',cmap_Tyaw(color_nr,:),'markersize',10,'color',[.5 .5 .5])

        subplot(2,2,4)
        plot(Arot_fwd_all(i),dMwing_aero_fwd_all(i),'dk-','markerfacecolor',cmap_Tyaw(color_nr,:),'markersize',10)
        plot(Arot_rwd_all(i),dMwing_aero_rwd_all(i),'o-','markerfacecolor',cmap_Tyaw(color_nr,:),'markersize',10,'color',[.5 .5 .5])
    end

    subplot(2,2,1)

%     plot(Arot_fwd_all,dMwing_total_fwd_all,'-k')
%     plot(Arot_rwd_all,dMwing_total_rwd_all,'-','color',[.5 .5 .5])

    plot([min(Arot_fwd_all) max(Arot_fwd_all)],polyval(k_ArotdMtotal_fwd_coeffs_norm_deg,[min(Arot_fwd_all) max(Arot_fwd_all)]),'-k','linewidth',2)
    plot([min(Arot_rwd_all) max(Arot_rwd_all)],polyval(k_ArotdMtotal_rwd_coeffs_norm_deg,[min(Arot_rwd_all) max(Arot_rwd_all)]),'-k','linewidth',2,'color',[.5 .5 .5])

    axis tight
    axis square
    dM_min_sub = -.02;
    axis([Arot_min,Arot_max,dM_min_sub,dM_min_sub+dM_max-dM_min])
    set(gca,'xtick',Arot_min:(Arot_max-Arot_min):Arot_max)
    set(gca,'ytick',[dM_min_sub dM_min_sub+dM_max-dM_min])
    xlabel('angle amplitude [deg]')
    ylabel('torque change')
    title('total torque')

    subplot(2,2,2)

%     plot(Arot_fwd_all,dMwing_inertia_fwd_all,'-k')
%     plot(Arot_rwd_all,dMwing_inertia_rwd_all,'-','color',[.5 .5 .5])

    plot([min(Arot_fwd_all) max(Arot_fwd_all)],polyval(k_ArotdMinertia_fwd_coeffs_norm_deg,[min(Arot_fwd_all) max(Arot_fwd_all)]),'-k','linewidth',2)
    plot([min(Arot_rwd_all) max(Arot_rwd_all)],polyval(k_ArotdMinertia_rwd_coeffs_norm_deg,[min(Arot_rwd_all) max(Arot_rwd_all)]),'-k','linewidth',2,'color',[.5 .5 .5])

    axis tight
    axis square
    dM_min_sub = -.12;
    axis([Arot_min,Arot_max,dM_min_sub,dM_min_sub+dM_max-dM_min])
    set(gca,'xtick',Arot_min:(Arot_max-Arot_min):Arot_max)
    set(gca,'ytick',[dM_min_sub dM_min_sub+dM_max-dM_min])
    xlabel('angle amplitude [deg]')
    ylabel('torque change')
    title('intertial torque')
    
    subplot(2,2,3)

%     plot(Arot_fwd_all,dMwing_mass_fwd_all,'-k')
%     plot(Arot_rwd_all,dMwing_mass_rwd_all,'-','color',[.5 .5 .5])

    plot([min(Arot_fwd_all) max(Arot_fwd_all)],polyval(k_ArotdMmass_fwd_coeffs_norm_deg,[min(Arot_fwd_all) max(Arot_fwd_all)]),'-k','linewidth',2)
    plot([min(Arot_rwd_all) max(Arot_rwd_all)],polyval(k_ArotdMmass_rwd_coeffs_norm_deg,[min(Arot_rwd_all) max(Arot_rwd_all)]),'-k','linewidth',2,'color',[.5 .5 .5])

    axis tight
    axis square
    dM_min_sub = -.02;
    axis([Arot_min,Arot_max,dM_min_sub,dM_min_sub+dM_max-dM_min])
    set(gca,'xtick',Arot_min:(Arot_max-Arot_min):Arot_max)
    set(gca,'ytick',[dM_min_sub dM_min_sub+dM_max-dM_min])
    xlabel('angle amplitude [deg]')
    ylabel('torque change')
    title('mass torque')

    subplot(2,2,4)

%     plot(Arot_fwd_all,dMwing_aero_fwd_all,'-k')
%     plot(Arot_rwd_all,dMwing_aero_rwd_all,'-','color',[.5 .5 .5])

    plot([min(Arot_fwd_all) max(Arot_fwd_all)],polyval(k_ArotdMaero_fwd_coeffs_norm_deg,[min(Arot_fwd_all) max(Arot_fwd_all)]),'-k','linewidth',2)
    plot([min(Arot_rwd_all) max(Arot_rwd_all)],polyval(k_ArotdMaero_rwd_coeffs_norm_deg,[min(Arot_rwd_all) max(Arot_rwd_all)]),'-k','linewidth',2,'color',[.5 .5 .5])

    axis tight
    axis square
    dM_min_sub = .09;
    axis([Arot_min,Arot_max,dM_min_sub,dM_min_sub+dM_max-dM_min])
    set(gca,'xtick',Arot_min:(Arot_max-Arot_min):Arot_max)
    set(gca,'ytick',[dM_min_sub dM_min_sub+dM_max-dM_min])
    xlabel('angle amplitude [deg]')
    ylabel('torque change')
    title('aerodynamic torque')
    
    saveas(gca,['Arot_vs_dTspanwise_vs_YawTorque_components_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.fig'])
    saveas(gca,['Arot_vs_dTspanwise_vs_YawTorque_components_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.png'])
    plot2svg(['Arot_vs_dTspanwise_vs_YawTorque_components_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.svg'])
    
%% plot rot0 as function of Tyaw
figure

% rot0 @ rot_max
subplot(2,2,1)
hold on
plot(YawTorques_all,rot0_max_fwd_all_deg,'dk','markerfacecolor','r','markersize',10)
plot(YawTorques_all,rot0_max_rwd_all_deg,'ok','markerfacecolor','b','markersize',10)

plot([min(YawTorques_all) max(YawTorques_all)],polyval(rot0Tyaw_coeff_fwd_rotmax_deg,[min(YawTorques_all) max(YawTorques_all)]),'-k','linewidth',2)
plot([min(YawTorques_all) max(YawTorques_all)],polyval(rot0Tyaw_coeff_rwd_rotmax_deg,[min(YawTorques_all) max(YawTorques_all)]),'-b','color',[.5 .5 .5],'linewidth',2)

legend('fwd@ds','rwd@ds')
axis tight
% axis square
axis([YawTorque_min,YawTorque_max,rot0max_min,rot0max_max])
set(gca,'xtick',YawTorque_min:(YawTorque_max-YawTorque_min):YawTorque_max)
set(gca,'ytick',rot0max_min:(rot0max_max-rot0max_min):rot0max_max)
xlabel('normalized yaw torque')
ylabel('rot0 [deg]')
title('rot0 @ mid downstroke')

% rot0 @ rot_min
subplot(2,2,2)
hold on
plot(YawTorques_all,rot0_min_fwd_all_deg,'dk','markerfacecolor','r','markersize',10)
plot(YawTorques_all,rot0_min_rwd_all_deg,'ok','markerfacecolor','b','markersize',10)

plot([min(YawTorques_all) max(YawTorques_all)],polyval(rot0Tyaw_coeff_fwd_rotmin_deg,[min(YawTorques_all) max(YawTorques_all)]),'-k','linewidth',2)
plot([min(YawTorques_all) max(YawTorques_all)],polyval(rot0Tyaw_coeff_rwd_rotmin_deg,[min(YawTorques_all) max(YawTorques_all)]),'-','color',[.5 .5 .5],'linewidth',2)

legend('fwd@us','rwd@us')
axis tight
% axis square
axis([YawTorque_min,YawTorque_max,rot0min_min,rot0min_max])
set(gca,'xtick',YawTorque_min:(YawTorque_max-YawTorque_min):YawTorque_max)
set(gca,'ytick',rot0min_min:(rot0min_max-rot0min_min):rot0min_max)
xlabel('normalized yaw torque')
ylabel('rot0 [deg]')
title('rot0 @ mid upstroke')

saveas(gca,['rot0_vs_YawTorque_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.fig'])
saveas(gca,['rot0_vs_YawTorque_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.png'])
plot2svg(['rot0_vs_YawTorque_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.svg'])

%% plot ALL MODs
%     figure
%     subplot(1,2,1)
%     hold on
%     plot(YawTorques_all,Fx_mean_all,'o-k','markersize',10,'markerfacecolor','b')
%     plot(YawTorques_all,Fy_mean_all,'o-k','markersize',10,'markerfacecolor','r')
%     plot(YawTorques_all,Fz_mean_all,'o-k','markersize',10,'markerfacecolor','g')
%     
%     legend('x-axis Mod','y-axis Mod','z-axis Mod','location','E')
%     xlabel('Yaw Torque')
%     ylabel('normalized forces F/mg')
%     axis([YawTorque_min YawTorque_max -1 .25])
%     set(gca,'xtick',[YawTorque_min YawTorque_max])
%     set(gca,'ytick',[-1 0 .25])
%     
%     subplot(1,2,2)
%     hold on
%     plot(YawTorques_all,Mx_mean_all,'o-k','markersize',10,'markerfacecolor','b')
%     plot(YawTorques_all,My_mean_all-My_mean_all(YawTorques_all==0),'o-k','markersize',10,'markerfacecolor','r')
%     plot(YawTorques_all,Mz_mean_all,'o-k','markersize',10,'markerfacecolor','g')
%     
%     xlabel('Yaw Torque')
%     ylabel('normalized torques T/mgl')
%     axis([YawTorque_min YawTorque_max -.01 .03])
%     set(gca,'xtick',[YawTorque_min YawTorque_max])
%     set(gca,'ytick',[-.01 0 .03])    
%     
%     saveas(gca,['FnM_YawMod.fig'])
%     saveas(gca,['FnM_YawMod.png'])
%     plot2svg(['FnM_YawMod.svg'])

%% plot separate MODs Fz & Mx
%     figure
%     subplot(1,2,1)
%     hold on
%     plot(YawTorques_all,Fz_mean_sum,'o-k','markersize',10,'markerfacecolor','w')
%     plot(YawTorques_all,Fz_mean_all,'o-k','markersize',10,'markerfacecolor',[.5 .5 .5])
%     plot(YawTorques_all,Fz_mean_freqMOD,'o-k','markersize',10,'markerfacecolor','b')
%     plot(YawTorques_all,Fz_mean_strokeMOD,'o-k','markersize',10,'markerfacecolor','g')
%     plot(YawTorques_all,Fz_mean_devMOD,'o-k','markersize',10,'markerfacecolor','r')
%     plot(YawTorques_all,Fz_mean_rotMOD,'o-k','markersize',10,'markerfacecolor','c')
%     
%     legend('sum','modified','frequency','stroke','deviation','rotation','location','SE')
%     xlabel('Yaw Torque')
%     ylabel('normalized vertical force Fz/mg')
%     axis([YawTorque_min YawTorque_max -1 .25])
%     set(gca,'xtick',[YawTorque_min YawTorque_max])
%     set(gca,'ytick',[-1 0 .25])
%     
%     subplot(1,2,2)
%     hold on
%     plot(YawTorques_all,Mx_mean_sum,'o-k','markersize',10,'markerfacecolor','w')
%     plot(YawTorques_all,Mx_mean_all,'o-k','markersize',10,'markerfacecolor',[.5 .5 .5])
%     plot(YawTorques_all,Mx_mean_freqMOD,'o-k','markersize',10,'markerfacecolor','b')
%     plot(YawTorques_all,Mx_mean_strokeMOD,'o-k','markersize',10,'markerfacecolor','g')
%     plot(YawTorques_all,Mx_mean_devMOD,'o-k','markersize',10,'markerfacecolor','r')
%     plot(YawTorques_all,Mx_mean_rotMOD,'o-k','markersize',10,'markerfacecolor','c')
%     
%     xlabel('Yaw Torque')
%     ylabel('normalized roll torque Tx/mgl')
%     axis([YawTorque_min YawTorque_max -.01 .03])
%     set(gca,'xtick',[YawTorque_min YawTorque_max])
%     set(gca,'ytick',[-.01 0 .03])    
% 
%     saveas(gca,['FnM_YawModComponents.fig'])
%     saveas(gca,['FnM_YawModComponents.png'])
%     plot2svg(['FnM_YawModComponents.svg'])
    
    cd ..
