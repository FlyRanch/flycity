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
% rot_on=0;

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

% nr_timepoints = 500;
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
    rot_L_wingref = kine.eta_L;
    rot_R_wingref = kine.eta_R;

    rot_ddot_L = kine.eta_ddot_L;
    rot_ddot_R = kine.eta_ddot_R;
    
    Maero_sect_L = FM_strkpln.FM_L.My_L;
    Maero_sect_R = FM_strkpln.FM_R.My_R;
    
    clear ax_L ax_R az_L az_R Mmass_sect_L Mmass_sect_R
    for n = 1:length(rot_L_wingref)
        
        % torque arm in wing ref frame
        x_cm_L = (x_sect-c_sect/2);
        x_cm_R = (x_sect-c_sect/2);

        az_L(:,:) = kine.Udot_left(3,n,:)*1e-3;      % [ m/s^2 ]
        az_R(:,:) = kine.Udot_right(3,n,:)*1e-3;     % [ m/s^2 ]
        
        Mmass_sect_L(:,n) = mass_tot_sect .* x_cm_L .*az_L;
        Mmass_sect_R(:,n) = mass_tot_sect .* x_cm_R .*az_R;
    end

    Mwing_mass_fwd_all(:,i)  = nansum(Mmass_sect_L)' / Mg_fly / l_wing;
    Mwing_mass_rwd_all(:,i)  = nansum(Mmass_sect_R)' / Mg_fly / l_wing;
    
    Mwing_inertia_fwd_all(:,i)  = rot_ddot_L * Itot / Mg_fly / l_wing;
    Mwing_inertia_rwd_all(:,i)  = rot_ddot_R * Itot / Mg_fly / l_wing;
    
    Mwing_aero_fwd_all(:,i)  = nansum(Maero_sect_L)' / Mg_fly / (l_wing*1e3);
    Mwing_aero_rwd_all(:,i)  = nansum(Maero_sect_R)' / Mg_fly / (l_wing*1e3);
    
    %% sum of torques equals angular acceleration: I*rot_ddot = Taero + Tmass + Tspring
    Mwing_spring_fwd_all(:,i)  = Mwing_inertia_fwd_all(:,i) - Mwing_mass_fwd_all(:,i) - Mwing_aero_fwd_all(:,i);
    Mwing_spring_rwd_all(:,i)  = Mwing_inertia_rwd_all(:,i) - Mwing_mass_rwd_all(:,i) - Mwing_aero_rwd_all(:,i);
    
    %% torque at min & max rotation angle
    [rot_maxpeaks_L,rot_maxpeaklocs_L] = findpeaks(rot_L);
    if loc_set == 3
        n_rot_us_L = rot_maxpeaklocs_L(maxpeak_nr);
        rot_us_L = rot_L(n_rot_us_L);
    else
        n_rot_ds_L = rot_maxpeaklocs_L(maxpeak_nr);
        rot_ds_L = rot_L(n_rot_ds_L);
    end
    
    [rot_minpeaks_L,rot_minpeaklocs_L] = findpeaks(-rot_L);
    if loc_set == 3
        n_rot_ds_L = rot_minpeaklocs_L(minpeak_nr);
        rot_ds_L = rot_L(n_rot_ds_L);
    else
        n_rot_us_L = rot_minpeaklocs_L(minpeak_nr);
        rot_us_L = rot_L(n_rot_us_L);
    end

    [rot_maxpeaks_R,rot_maxpeaklocs_R] = findpeaks(rot_R);
    if loc_set == 3
        n_rot_us_R = rot_maxpeaklocs_R(maxpeak_nr);
        rot_us_R = rot_R(n_rot_us_R);
    else
        n_rot_ds_R = rot_maxpeaklocs_R(maxpeak_nr);
        rot_ds_R = rot_R(n_rot_ds_R);
    end
    
    [rot_minpeaks_R,rot_minpeaklocs_R] = findpeaks(-rot_R);
    if loc_set == 3
        n_rot_ds_R = rot_minpeaklocs_R(minpeak_nr);
        rot_ds_R = rot_R(n_rot_ds_R);
    else
        n_rot_us_R = rot_minpeaklocs_R(minpeak_nr);
        rot_us_R = rot_R(n_rot_us_R);
    end

    % plot
    subplot(2,1,1)
    plot(t_norm,rot_L)
    hold on
    plot(t_norm(rot_maxpeaklocs_L),rot_L(rot_maxpeaklocs_L),'or')
    plot(t_norm(rot_minpeaklocs_L),rot_L(rot_minpeaklocs_L),'og')
    
    plot(t_norm(n_rot_ds_L),rot_ds_L,'*c')
    plot(t_norm(n_rot_us_L),rot_us_L,'*m')

    subplot(2,1,2)
    plot(t_norm,rot_R)
    hold on
    plot(t_norm(rot_maxpeaklocs_R),rot_R(rot_maxpeaklocs_R),'or')
    plot(t_norm(rot_minpeaklocs_R),rot_R(rot_minpeaklocs_R),'og')
    
    plot(t_norm(n_rot_ds_R),rot_ds_R,'*c')
    plot(t_norm(n_rot_us_R),rot_us_R,'*m')

    %% store data
    % rot@ds & rot@us
    n_rot_ds_fwd_all(i,:) = n_rot_ds_L;
    n_rot_ds_rwd_all(i,:) = n_rot_ds_R;
    
    n_rot_us_fwd_all(i,:) = n_rot_us_L;
    n_rot_us_rwd_all(i,:) = n_rot_us_R;
    
    rot_ds_fwd_all(i,:) = rot_ds_L;
    rot_ds_rwd_all(i,:) = rot_ds_R;
    
    rot_us_fwd_all(i,:) = rot_us_L;
    rot_us_rwd_all(i,:) = rot_us_R;
    
    % spanwise torques @ rot ds & rot us
    Mwing_spring_fwd_atRotds(i,1) = Mwing_spring_fwd_all(n_rot_ds_L,i);
    Mwing_spring_fwd_atRotus(i,1) = Mwing_spring_fwd_all(n_rot_us_L,i);
    
    Mwing_spring_rwd_atRotds(i,1) = Mwing_spring_rwd_all(n_rot_ds_R,i);
    Mwing_spring_rwd_atRotus(i,1) = Mwing_spring_rwd_all(n_rot_us_R,i);

    Mwing_mass_fwd_atRotds(i,1) = Mwing_mass_fwd_all(n_rot_ds_L,i);
    Mwing_mass_fwd_atRotus(i,1) = Mwing_mass_fwd_all(n_rot_us_L,i);
    
    Mwing_mass_rwd_atRotds(i,1) = Mwing_mass_rwd_all(n_rot_ds_R,i);
    Mwing_mass_rwd_atRotus(i,1) = Mwing_mass_rwd_all(n_rot_us_R,i);

    Mwing_inertia_fwd_atRotds(i,1) = Mwing_inertia_fwd_all(n_rot_ds_L,i);
    Mwing_inertia_fwd_atRotus(i,1) = Mwing_inertia_fwd_all(n_rot_us_L,i);
    
    Mwing_inertia_rwd_atRotds(i,1) = Mwing_inertia_rwd_all(n_rot_ds_R,i);
    Mwing_inertia_rwd_atRotus(i,1) = Mwing_inertia_rwd_all(n_rot_us_R,i);

    Mwing_aero_fwd_atRotds(i,1) = Mwing_aero_fwd_all(n_rot_ds_L,i);
    Mwing_aero_fwd_atRotus(i,1) = Mwing_aero_fwd_all(n_rot_us_L,i);
    
    Mwing_aero_rwd_atRotds(i,1) = Mwing_aero_rwd_all(n_rot_ds_R,i);
    Mwing_aero_rwd_atRotus(i,1) = Mwing_aero_rwd_all(n_rot_us_R,i);
    
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
    
    %% rotations & torques in 2D wing ref frame: NEGATIVE 3D WING REF FRAME
    rot_ds_fwd_all_deg = rad2deg(rot_ds_fwd_all) -90;
    rot_us_fwd_all_deg = rad2deg(rot_us_fwd_all) -90;

    rot_ds_rwd_all_deg = rad2deg(rot_ds_rwd_all) -90;
    rot_us_rwd_all_deg = rad2deg(rot_us_rwd_all) -90;
    
    Mwing_spring_fwd_all_2Dref = -Mwing_spring_fwd_all;
    Mwing_spring_rwd_all_2Dref = -Mwing_spring_rwd_all;
    
    Mwing_spring_fwd_atRotds_2Dref = -Mwing_spring_fwd_atRotds;
    Mwing_spring_rwd_atRotds_2Dref = -Mwing_spring_rwd_atRotds;
    
    Mwing_spring_fwd_atRotus_2Dref = -Mwing_spring_fwd_atRotus;
    Mwing_spring_rwd_atRotus_2Dref = -Mwing_spring_rwd_atRotus;
    
    %% calc winghinge rotational spring coeffs
    dMwing_spring_fwd = Mwing_spring_fwd_atRotds_2Dref - Mwing_spring_fwd_atRotus_2Dref;
    dMwing_spring_rwd = Mwing_spring_rwd_atRotds_2Dref - Mwing_spring_rwd_atRotus_2Dref;
    
    Arot_fwd_deg = rot_ds_fwd_all_deg - rot_us_fwd_all_deg;
    Arot_rwd_deg = rot_ds_rwd_all_deg - rot_us_rwd_all_deg;
    
    Arot_fwd_rad = rot_ds_fwd_all - rot_us_fwd_all;
    Arot_rwd_rad = rot_ds_rwd_all - rot_us_rwd_all;
    
    % spring constants
    k_fwd_norm_deg = -dMwing_spring_fwd ./ Arot_fwd_deg;
    k_rwd_norm_deg = -dMwing_spring_rwd ./ Arot_rwd_deg;
    
    k_fwd_norm_rad = -dMwing_spring_fwd ./ Arot_fwd_rad;
    k_rwd_norm_rad = -dMwing_spring_rwd ./ Arot_rwd_rad;
    
    % spring setpoint
    rot0_fwd_norm_deg = rot_ds_fwd_all_deg + Mwing_spring_fwd_atRotds_2Dref./k_fwd_norm_deg;
    rot0_rwd_norm_deg = rot_ds_rwd_all_deg + Mwing_spring_rwd_atRotds_2Dref./k_rwd_norm_deg;

    rot0_fwd_norm_rad = rot_ds_fwd_all + Mwing_spring_fwd_atRotds_2Dref./k_fwd_norm_rad;
    rot0_rwd_norm_rad = rot_ds_rwd_all + Mwing_spring_rwd_atRotds_2Dref./k_rwd_norm_rad;
    
    %% save data & plots
    save(['hingespring_k_alpha0_Arot_dTspanwise_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.mat'])
    
    mkdir('qsModel_FnM_k_alpha0_n_rotVsTspanwise')
    cd('qsModel_FnM_k_alpha0_n_rotVsTspanwise')
    
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

plot(t_norm(n_rot_ds_fwd_all(1)),rot_fwd_all(n_rot_ds_fwd_all(1),1)-90,'ok','linewidth',2)
plot(t_norm(n_rot_us_fwd_all(1)),rot_fwd_all(n_rot_us_fwd_all(1),1)-90,'ok','linewidth',2)

plot(t_norm(n_rot_ds_fwd_all(end)),rot_fwd_all(n_rot_ds_fwd_all(end),end)-90,'or','linewidth',2)
plot(t_norm(n_rot_us_fwd_all(end)),rot_fwd_all(n_rot_us_fwd_all(end),end)-90,'or','linewidth',2)

plot(t_norm(n_rot_ds_rwd_all(end)),rot_rwd_all(n_rot_ds_rwd_all(end),end)-90,'ob','linewidth',2)
plot(t_norm(n_rot_us_rwd_all(end)),rot_rwd_all(n_rot_us_rwd_all(end),end)-90,'ob','linewidth',2)

legend('steady','fwd','rwd')
axis tight
% axis square
axis([0,1,-75,75])
set(gca,'xtick',[0 1])
set(gca,'ytick',[-60 0 60])
xlabel('wingbeat cycle')
ylabel('wing rotation angle [deg]')

saveas(gca,['rot_vs_WBcycle_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.fig'])
saveas(gca,['rot_vs_WBcycle_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.png'])
plot2svg(['rot_vs_WBcycle_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.svg'])

%% plot spanwise hinge torque throughout wingbeat
figure
hold on
plot(t_norm,Mwing_spring_fwd_all_2Dref(:,1),'-k','linewidth',2)
plot(t_norm,Mwing_spring_fwd_all_2Dref(:,end),'-r','linewidth',2)
plot(t_norm,Mwing_spring_rwd_all_2Dref(:,end),'-b','linewidth',2)

plot(t_norm(n_rot_ds_fwd_all(1)),Mwing_spring_fwd_all_2Dref(n_rot_ds_fwd_all(1),1),'ok','linewidth',2)
plot(t_norm(n_rot_us_fwd_all(1)),Mwing_spring_fwd_all_2Dref(n_rot_us_fwd_all(1),1),'ok','linewidth',2)

plot(t_norm(n_rot_ds_fwd_all(end)),Mwing_spring_fwd_all_2Dref(n_rot_ds_fwd_all(end),end),'or','linewidth',2)
plot(t_norm(n_rot_us_fwd_all(end)),Mwing_spring_fwd_all_2Dref(n_rot_us_fwd_all(end),end),'or','linewidth',2)

plot(t_norm(n_rot_ds_rwd_all(end)),Mwing_spring_rwd_all_2Dref(n_rot_ds_rwd_all(end),end),'ob','linewidth',2)
plot(t_norm(n_rot_us_rwd_all(end)),Mwing_spring_rwd_all_2Dref(n_rot_us_rwd_all(end),end),'ob','linewidth',2)

legend('steady','fwd','rwd')
grid on
axis tight
% axis square
axis([0,1,-.25,.25])
set(gca,'xtick',[0 1])
set(gca,'ytick',[-.25 0 .25])
xlabel('wingbeat cycle')
ylabel('normalized spanwise hinge torque')

saveas(gca,['Twing_vs_WBcycle_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.fig'])
saveas(gca,['Twing_vs_WBcycle_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.png'])
plot2svg(['Twing_vs_WBcycle_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.svg'])

%% plot spanwise hinge torque versus rotation angles
figure
hold on
plot(rot_fwd_all(:,1)-90,Mwing_spring_fwd_all_2Dref(:,1),'-k','linewidth',2)
plot(rot_fwd_all(:,end)-90,Mwing_spring_fwd_all_2Dref(:,end),'-r','linewidth',2)
plot(rot_rwd_all(:,end)-90,Mwing_spring_rwd_all_2Dref(:,end),'-b','linewidth',2)

plot(rot_fwd_all(n_rot_ds_fwd_all(1),1)-90,Mwing_spring_fwd_all_2Dref(n_rot_ds_fwd_all(1),1),'ok','linewidth',2)
plot(rot_fwd_all(n_rot_us_fwd_all(1),1)-90,Mwing_spring_fwd_all_2Dref(n_rot_us_fwd_all(1),1),'ok','linewidth',2)

plot(rot_fwd_all(n_rot_ds_fwd_all(end),end)-90,Mwing_spring_fwd_all_2Dref(n_rot_ds_fwd_all(end),end),'or','linewidth',2)
plot(rot_fwd_all(n_rot_us_fwd_all(end),end)-90,Mwing_spring_fwd_all_2Dref(n_rot_us_fwd_all(end),end),'or','linewidth',2)

plot(rot_rwd_all(n_rot_ds_rwd_all(end),end)-90,Mwing_spring_rwd_all_2Dref(n_rot_ds_rwd_all(end),end),'ob','linewidth',2)
plot(rot_rwd_all(n_rot_us_rwd_all(end),end)-90,Mwing_spring_rwd_all_2Dref(n_rot_us_rwd_all(end),end),'ob','linewidth',2)

legend('steady','fwd','rwd')
grid on
axis tight
% axis square
axis([-75,75,-.25,.25])
set(gca,'xtick',[-75 0 75])
set(gca,'ytick',[-.25 0 .25])
xlabel('wing rotation angle [deg]')
ylabel('normalized spanwise hinge torque')

saveas(gca,['Twing_vs_rot_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.fig'])
saveas(gca,['Twing_vs_rot_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.png'])
plot2svg(['Twing_vs_rot_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.svg'])

%% plot rotation angles & spanwise torque throughout wingbeat
figure

% rotation angles
subplot(2,2,1)
hold on
plot(t_norm,rot_fwd_all(:,1)-90,'-k','linewidth',2)
plot(t_norm,rot_fwd_all(:,end)-90,'-r','linewidth',2)
plot(t_norm,rot_rwd_all(:,end)-90,'-b','linewidth',2)

plot(t_norm(n_rot_ds_fwd_all(1)),rot_fwd_all(n_rot_ds_fwd_all(1),1)-90,'ok','linewidth',2)
plot(t_norm(n_rot_us_fwd_all(1)),rot_fwd_all(n_rot_us_fwd_all(1),1)-90,'ok','linewidth',2)

plot(t_norm(n_rot_ds_fwd_all(end)),rot_fwd_all(n_rot_ds_fwd_all(end),end)-90,'or','linewidth',2)
plot(t_norm(n_rot_us_fwd_all(end)),rot_fwd_all(n_rot_us_fwd_all(end),end)-90,'or','linewidth',2)

plot(t_norm(n_rot_ds_rwd_all(end)),rot_rwd_all(n_rot_ds_rwd_all(end),end)-90,'ob','linewidth',2)
plot(t_norm(n_rot_us_rwd_all(end)),rot_rwd_all(n_rot_us_rwd_all(end),end)-90,'ob','linewidth',2)

legend('steady','fwd','rwd')
axis tight
% axis square
axis([0,1,-60,60])
set(gca,'xtick',[0 1])
set(gca,'ytick',[-60 0 60])
xlabel('wingbeat cycle')
ylabel('wing rotation angle [deg]')

%  spanwise hinge torque throughout wingbeat
subplot(2,2,2)
hold on
plot(t_norm,Mwing_spring_fwd_all_2Dref(:,1),'-k','linewidth',2)
plot(t_norm,Mwing_spring_fwd_all_2Dref(:,end),'-r','linewidth',2)
plot(t_norm,Mwing_spring_rwd_all_2Dref(:,end),'-b','linewidth',2)

plot(t_norm(n_rot_ds_fwd_all(1)),Mwing_spring_fwd_all_2Dref(n_rot_ds_fwd_all(1),1),'ok','linewidth',2)
plot(t_norm(n_rot_us_fwd_all(1)),Mwing_spring_fwd_all_2Dref(n_rot_us_fwd_all(1),1),'ok','linewidth',2)

plot(t_norm(n_rot_ds_fwd_all(end)),Mwing_spring_fwd_all_2Dref(n_rot_ds_fwd_all(end),end),'or','linewidth',2)
plot(t_norm(n_rot_us_fwd_all(end)),Mwing_spring_fwd_all_2Dref(n_rot_us_fwd_all(end),end),'or','linewidth',2)

plot(t_norm(n_rot_ds_rwd_all(end)),Mwing_spring_rwd_all_2Dref(n_rot_ds_rwd_all(end),end),'ob','linewidth',2)
plot(t_norm(n_rot_us_rwd_all(end)),Mwing_spring_rwd_all_2Dref(n_rot_us_rwd_all(end),end),'ob','linewidth',2)

legend('steady','fwd','rwd')
grid on
axis tight
% axis square
axis([0,1,-.25,.25])
set(gca,'xtick',[0 1])
set(gca,'ytick',[-.25 0 .25])
xlabel('wingbeat cycle')
ylabel('normalized spanwise hinge torque')

saveas(gca,['rot_n_Twing_vs_WBcycle_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.fig'])
saveas(gca,['rot_n_Twing_vs_WBcycle_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.png'])
plot2svg(['rot_n_Twing_vs_WBcycle_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.svg'])

%% plot rotation angle amplitude VS change in spanwise torque
cmap_hot =colormap(hot(100));
cmap_hot = flipud(cmap_hot);
cmap_Tyaw = cmap_hot;

if loc_set == 1
    rot_ds_min = 60;
    rot_ds_max = 65;
    
    M_ds_min = .06;
    M_ds_max = .12;

    rot_us_min = -71;
    rot_us_max = -64;
    
    M_us_min = -.17;
    M_us_max = -.09;

    rot0_min = -5;
    rot0_max =  15;

    k_min = -3.8e-3;
    k_max = -3.2e-3;
    
elseif loc_set == 2
    rot_ds_min = 52;
    rot_ds_max = 57;
    
    M_ds_min = -.030;
    M_ds_max = -.000;

    rot_us_min = -47;
    rot_us_max = -42;
    
    M_us_min = .040;
    M_us_max = .070;

    rot0_min = -5;
    rot0_max =  15;

    k_min = -1.6e-3;
    k_max = -1.3e-3;
    
elseif loc_set == 3
    rot_ds_min = 50;
    rot_ds_max = 55;
    
    M_ds_min = -.11;
    M_ds_max = -.08;

    rot_us_min = -47;
    rot_us_max = -40;
    
    M_us_min = .09;
    M_us_max = .12;

    rot0_min = -180;
    rot0_max =  180;

    k_min = -6e-4;
    k_max =  2e-4;
end

YawTorque_min = 0;
YawTorque_max = .1;


%% rot vs Tspanwise
figure

% @ downstroke
subplot(2,2,3)
hold on

    for i = 1:length(YawTorques_all)
        color_nr = round(99/(YawTorque_max-YawTorque_min)*(YawTorques_all(i)-YawTorque_min)+1);
        if color_nr<1
            color_nr=1;
        elseif color_nr>size(cmap_Tyaw,1)
            color_nr=size(cmap_Tyaw,1)
        end

        plot(rot_ds_fwd_all_deg(i),Mwing_spring_fwd_atRotds(i),'dk-','markerfacecolor',cmap_Tyaw(color_nr,:),'markersize',10)
        plot(rot_ds_rwd_all_deg(i),Mwing_spring_rwd_atRotds(i),'o-','markerfacecolor',cmap_Tyaw(color_nr,:),'markersize',10,'color',[.5 .5 .5])
    end

    linfit_rotMtotal_ds_fwd_norm_deg = polyfit(rot_ds_fwd_all_deg,Mwing_spring_fwd_atRotds,1);
    linfit_rotMtotal_ds_rwd_norm_deg = polyfit(rot_ds_rwd_all_deg,Mwing_spring_rwd_atRotds,1);

    plot([min(rot_ds_fwd_all_deg) max(rot_ds_fwd_all_deg)],polyval(linfit_rotMtotal_ds_fwd_norm_deg,[min(rot_ds_fwd_all_deg) max(rot_ds_fwd_all_deg)]),'-k','linewidth',2)
    plot([min(rot_ds_rwd_all_deg) max(rot_ds_rwd_all_deg)],polyval(linfit_rotMtotal_ds_rwd_norm_deg,[min(rot_ds_rwd_all_deg) max(rot_ds_rwd_all_deg)]),'-k','linewidth',2,'color',[.5 .5 .5])

    legend('fwd@ds','rwd@ds','location','southeast')
    axis tight
%     axis square
    axis([rot_ds_min,rot_ds_max,M_ds_min,M_ds_max])
    set(gca,'xtick',rot_ds_min:(rot_ds_max-rot_ds_min):rot_ds_max)
    set(gca,'ytick',M_ds_min:(M_ds_max-M_ds_min):M_ds_max)
    xlabel('rotation angle [deg]')
    ylabel('normalized spanwise torque')
%     colormap(cmap_Tyaw)
%     caxis([YawTorque_min YawTorque_max])
%     h = colorbar('location','northoutside'); 
%     title(h,'Normalized Yaw Torque')
%     set(h,'xtick',YawTorque_min:(YawTorque_max-YawTorque_min)/2:YawTorque_max)
title('downstroke')
    
% @ upstroke
subplot(2,2,4)
hold on

    for i = 1:length(YawTorques_all)
        color_nr = round(99/(YawTorque_max-YawTorque_min)*(YawTorques_all(i)-YawTorque_min)+1);
        if color_nr<1
            color_nr=1;
        elseif color_nr>size(cmap_Tyaw,1)
            color_nr=size(cmap_Tyaw,1)
        end

        plot(rot_us_fwd_all_deg(i),Mwing_spring_fwd_atRotus(i),'dk-','markerfacecolor',cmap_Tyaw(color_nr,:),'markersize',10)
        plot(rot_us_rwd_all_deg(i),Mwing_spring_rwd_atRotus(i),'o-','markerfacecolor',cmap_Tyaw(color_nr,:),'markersize',10,'color',[.5 .5 .5])
    end

    linfit_rotMtotal_us_fwd_norm_deg = polyfit(rot_us_fwd_all_deg,Mwing_spring_fwd_atRotus,1);
    linfit_rotMtotal_us_rwd_norm_deg = polyfit(rot_us_rwd_all_deg,Mwing_spring_rwd_atRotus,1);

    plot([min(rot_us_fwd_all_deg) max(rot_us_fwd_all_deg)],polyval(linfit_rotMtotal_us_fwd_norm_deg,[min(rot_us_fwd_all_deg) max(rot_us_fwd_all_deg)]),'-k','linewidth',2)
    plot([min(rot_us_rwd_all_deg) max(rot_us_rwd_all_deg)],polyval(linfit_rotMtotal_us_rwd_norm_deg,[min(rot_us_rwd_all_deg) max(rot_us_rwd_all_deg)]),'-k','linewidth',2,'color',[.5 .5 .5])

    legend('fwd@us','rwd@us','location','northwest')
    axis tight
%     axis square
    axis([rot_us_min,rot_us_max,M_us_min,M_us_max])
    set(gca,'xtick',rot_us_min:(rot_us_max-rot_us_min):rot_us_max)
    set(gca,'ytick',M_us_min:(M_us_max-M_us_min):M_us_max)
    xlabel('rotation angle [deg]')
    ylabel('normalized spanwise torque')
%     colormap(cmap_Tyaw)
%     caxis([YawTorque_min YawTorque_max])
%     h = colorbar('location','northoutside'); 
%     title(h,'Normalized Yaw Torque')
%     set(h,'xtick',YawTorque_min:(YawTorque_max-YawTorque_min)/2:YawTorque_max)
title('upstroke')
    
    subplot(2,2,1)
    colormap(cmap_Tyaw)
    caxis([YawTorque_min YawTorque_max])
    h = colorbar('location','northoutside'); 
    title(h,'Normalized Yaw Torque')
    set(h,'xtick',YawTorque_min:(YawTorque_max-YawTorque_min)/2:YawTorque_max)
    axis off
    
    subplot(2,2,2)
    colormap(cmap_Tyaw)
    caxis([YawTorque_min YawTorque_max])
    h = colorbar('location','northoutside'); 
    title(h,'Normalized Yaw Torque')
    set(h,'xtick',YawTorque_min:(YawTorque_max-YawTorque_min)/2:YawTorque_max)
    axis off
    
    saveas(gca,['rot_vs_Tspanwise_vs_YawTorque_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.fig'])
    saveas(gca,['rot_vs_Tspanwise_vs_YawTorque_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.png'])
    plot2svg(['rot_vs_Tspanwise_vs_YawTorque_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.svg'])

%% plot k & rot0 as function of Tyaw
figure

% k
subplot(2,2,1)
hold on
plot(YawTorques_all,k_fwd_norm_deg,'dk','markerfacecolor','r','markersize',10)
plot(YawTorques_all,k_rwd_norm_deg,'ok','markerfacecolor','b','markersize',10)

linfit_kTyaw_fwd_norm_deg = polyfit(YawTorques_all,k_fwd_norm_deg,1)
linfit_kTyaw_rwd_norm_deg = polyfit(YawTorques_all,k_rwd_norm_deg,1)

plot([min(YawTorques_all) max(YawTorques_all)],polyval(linfit_kTyaw_fwd_norm_deg,[min(YawTorques_all) max(YawTorques_all)]),'-k','linewidth',2)
plot([min(YawTorques_all) max(YawTorques_all)],polyval(linfit_kTyaw_rwd_norm_deg,[min(YawTorques_all) max(YawTorques_all)]),'-b','color',[.5 .5 .5],'linewidth',2)

legend('fwd@ds','rwd@ds')
axis tight
% % axis square
% axis([YawTorque_min,YawTorque_max,k_min,k_max])
% set(gca,'xtick',YawTorque_min:(YawTorque_max-YawTorque_min):YawTorque_max)
% set(gca,'ytick',k_min:(k_max-k_min):k_max)
xlabel('normalized yaw torque')
ylabel('spring stiffness')
title('k')

% rot0
subplot(2,2,2)
hold on
plot(YawTorques_all,rot0_fwd_norm_deg,'dk','markerfacecolor','r','markersize',10)
plot(YawTorques_all,rot0_rwd_norm_deg,'ok','markerfacecolor','b','markersize',10)

linfit_rot0Tyaw_fwd_norm_deg = polyfit(YawTorques_all,rot0_fwd_norm_deg,1)
linfit_rot0Tyaw_rwd_norm_deg = polyfit(YawTorques_all,rot0_rwd_norm_deg,1)

plot([min(YawTorques_all) max(YawTorques_all)],polyval(linfit_rot0Tyaw_fwd_norm_deg,[min(YawTorques_all) max(YawTorques_all)]),'-k','linewidth',2)
plot([min(YawTorques_all) max(YawTorques_all)],polyval(linfit_rot0Tyaw_rwd_norm_deg,[min(YawTorques_all) max(YawTorques_all)]),'-b','color',[.5 .5 .5],'linewidth',2)

% legend('fwd@ds','rwd@ds')
axis tight
% % axis square
% axis([YawTorque_min,YawTorque_max,rot0_min,rot0_max])
% set(gca,'xtick',YawTorque_min:(YawTorque_max-YawTorque_min):YawTorque_max)
% set(gca,'ytick',rot0_min:(rot0_max-rot0_min):rot0_max)
xlabel('normalized yaw torque')
ylabel('rot0 [deg]')
title('rot0')

saveas(gca,['k_n_rot0_vs_YawTorque_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.fig'])
saveas(gca,['k_n_rot0_vs_YawTorque_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.png'])
plot2svg(['k_n_rot0_vs_YawTorque_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.svg'])

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
