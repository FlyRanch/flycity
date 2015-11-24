clear;
clc;
close all
warning off

load('bodyNwingModel_4qsModel.mat')

loadname_steady=dir('WBdataset_steady_1603WBs.mat')
loadname_steady = loadname_steady.name;
load(loadname_steady)

plot_on = 1
% plot_on = 0

% qs model INC rotational lift
rot_on=1;
% rot_on=0;

% location of peak rotation angle
loc_set = 2
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
axis_loc = [1];

nr_sect = settings.nr_chord_sect;
nr_timepoints = settings.nr_timepoints;

% nr_timepoints = 500;
% nr_timepoints = 5000;

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

%% loop with different rotation axes
for i = 1:length(axis_loc)

    %% WB kin ALL MODs
    freq = f_steady;
    stroke = deg2rad(stroke_steady);
    dev = deg2rad(dev_steady);
    rot = deg2rad(rot_steady);

    % qs forces & torques
    [ FM_strkpln, kine ] = quasi_steady_FnMnWingkin_atTimeNspanSections( settings, body_model, wing_model, dev, rot, stroke, dev, rot, stroke, freq, rot_on );
    
    %% calc spanwise wing torque & rotation angle amplitude 
    rot_dot = kine.eta_dot_R;
    rot_ddot = kine.eta_ddot_R;
    Maero_sect = FM_strkpln.FM_R.My_R;
    
    clear az Mmass_sect
    for n = 1:length(rot_ddot)
        
        % torque arm in wing ref frame
        x_cm = (x_sect-c_sect/2);
        az(:,:) = kine.Udot_left(3,n,:)*1e-3;      % [ m/s^2 ]
        
        Mmass_sect(:,n) = mass_tot_sect .* x_cm .*az;
    end

    Mwing_mass_all(:,i)  = nansum(Mmass_sect)' / Mg_fly / l_wing;
    Mwing_inertia_all(:,i)  = rot_ddot * Itot / Mg_fly / l_wing;
    Mwing_aero_all(:,i)  = nansum(Maero_sect)' / Mg_fly / (l_wing*1e3);
    
    %% sum of torques equals angular acceleration: I*rot_ddot + C*rot_dot + k*(rot-rot0) = Taero + Tmass
    % hing torque = - wing torque
    Mwing_springdamper_all(:,i)  = Mwing_inertia_all(:,i) - Mwing_mass_all(:,i) - Mwing_aero_all(:,i);
    Mwing_inertiamass_all(:,i)   = Mwing_inertia_all(:,i) - Mwing_mass_all(:,i);
    
    %% torque at min & max rotation angle
    [rot_maxpeaks,rot_maxpeaklocs] = findpeaks(rot);
    if loc_set == 3
        n_rot_us = rot_maxpeaklocs(maxpeak_nr);
        rot_us = rot(n_rot_us);
    else
        n_rot_ds = rot_maxpeaklocs(maxpeak_nr);
        rot_ds = rot(n_rot_ds);
    end
    
    [rot_minpeaks,rot_minpeaklocs] = findpeaks(-rot);
    if loc_set == 3
        n_rot_ds = rot_minpeaklocs(minpeak_nr);
        rot_ds = rot(n_rot_ds);
    else
        n_rot_us = rot_minpeaklocs(minpeak_nr);
        rot_us = rot(n_rot_us);
    end

    %% store data
    % rot@ds & rot@us
    n_rot_ds_all(i,:) = n_rot_ds;
    n_rot_us_all(i,:) = n_rot_us;
    rot_ds_all(i,:) = rot_ds;
    rot_us_all(i,:) = rot_us;
    
    % spanwise torques @ rot ds & rot us
    Mwing_springdamper_atRotds(i,1) = Mwing_springdamper_all(n_rot_ds,i);
    Mwing_springdamper_atRotus(i,1) = Mwing_springdamper_all(n_rot_us,i);

    Mwing_inertiamass_atRotds(i,1) = Mwing_inertiamass_all(n_rot_ds,i);
    Mwing_inertiamass_atRotus(i,1) = Mwing_inertiamass_all(n_rot_us,i);

    Mwing_inertia_atRotds(i,1) = Mwing_inertia_all(n_rot_ds,i);
    Mwing_inertia_atRotus(i,1) = Mwing_inertia_all(n_rot_us,i);

    Mwing_mass_atRotds(i,1) = Mwing_mass_all(n_rot_ds,i);
    Mwing_mass_atRotus(i,1) = Mwing_mass_all(n_rot_us,i);

    Mwing_aero_atRotds(i,1) = Mwing_aero_all(n_rot_ds,i);
    Mwing_aero_atRotus(i,1) = Mwing_aero_all(n_rot_us,i);
end
    
    %% rotations & torques in 2D wing ref frame: NEGATIVE 3D WING REF FRAME
    rot_dot_deg = -rad2deg(rot_dot);
    
    rot_ds_all_deg = rad2deg(rot_ds_all) -90;
    rot_us_all_deg = rad2deg(rot_us_all) -90;
    
    Mwing_springdamper_all_2Dref = -Mwing_springdamper_all;
    Mwing_springdamper_atRotds_2Dref = -Mwing_springdamper_atRotds;
    Mwing_springdamper_atRotus_2Dref = -Mwing_springdamper_atRotus;
    
    Mwing_inertiamass_all_2Dref = -Mwing_inertiamass_all;
    Mwing_inertiamass_atRotds_2Dref = -Mwing_inertiamass_atRotds;
    Mwing_inertiamass_atRotus_2Dref = -Mwing_inertiamass_atRotus;
    
    Mwing_inertia_all_2Dref = -Mwing_inertia_all;
    Mwing_inertia_atRotds_2Dref = -Mwing_inertia_atRotds;
    Mwing_inertia_atRotus_2Dref = -Mwing_inertia_atRotus;

    Mwing_mass_all_2Dref = -Mwing_mass_all;
    Mwing_mass_atRotds_2Dref = -Mwing_mass_atRotds;
    Mwing_mass_atRotus_2Dref = -Mwing_mass_atRotus;
    
    Mwing_aero_all_2Dref = -Mwing_aero_all;
    Mwing_aero_atRotds_2Dref = -Mwing_aero_atRotds;
    Mwing_aero_atRotus_2Dref = -Mwing_aero_atRotus;
    
    %% calc winghinge rotational springdamper coeffs
    p = polyfitn([rot_dot_deg,rot_steady],Mwing_springdamper_all_2Dref,1);
    C_norm_deg_model = p.Coefficients(1);
    k_norm_deg_model = p.Coefficients(2);
    rot0_norm_deg_model = -p.Coefficients(3)/p.Coefficients(2);
    Mwing_springdamper_all_2Dref_model = C_norm_deg_model*rot_dot_deg + k_norm_deg_model*(rot_steady-rot0_norm_deg_model);

%     figure
%     plot(rot_steady,Mwing_springdamper_all_2Dref)
%     hold on
%     plot(rot_steady,Mwing_springdamper_all_2Dref_model_total,'r')

    % at zero rot rate (damping zero)
    dMwing_springdamper = Mwing_springdamper_atRotds_2Dref - Mwing_springdamper_atRotus_2Dref;
    
    Arot_deg = rot_ds_all_deg - rot_us_all_deg;
    Arot_rad = rot_ds_all - rot_us_all;
    
    % springdamper constants
    k_norm_deg = -dMwing_springdamper ./ Arot_deg;
    k_norm_rad = -dMwing_springdamper ./ Arot_rad;
    
    % springdamper setpoint
    rot0_norm_deg = rot_ds_all_deg + Mwing_springdamper_atRotds_2Dref./k_norm_deg;
    rot0_norm_rad = rot_ds_all + Mwing_springdamper_atRotds_2Dref./k_norm_rad;
    
    %% save data & plots
    save(['hingespringdamper_k_alpha0_STEADY_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.mat'])
    
    mkdir('qsModel_FnM_k_alpha0_n_rotVsTspanwise')
    cd('qsModel_FnM_k_alpha0_n_rotVsTspanwise')
    
%% plot rotation angles throughout wingbeat
figure
subplot(2,2,1)
hold on
plot(t_norm,rot_steady-90,'-k','linewidth',2)
plot(t_norm(n_rot_ds_all),rot_steady(n_rot_ds_all,1)-90,'ok','linewidth',2)
plot(t_norm(n_rot_us_all),rot_steady(n_rot_us_all,1)-90,'ok','linewidth',2)

axis tight
% axis square
axis([0,1,-75,75])
set(gca,'xtick',[0 1])
set(gca,'ytick',[-75 0 75])
xlabel('wingbeat cycle')
ylabel('wing rotation angle [deg]')
grid on

subplot(2,2,3)
hold on
plot(t_norm,rot_dot_deg,'-k','linewidth',2)
plot(t_norm(n_rot_ds_all),rot_dot_deg(n_rot_ds_all,1)-90,'ok','linewidth',2)
plot(t_norm(n_rot_us_all),rot_dot_deg(n_rot_us_all,1)-90,'ok','linewidth',2)

axis tight
% axis square
axis([0,1,-2e5,2e5])
set(gca,'xtick',[0 1])
set(gca,'ytick',[-2e5 0 2e5])
xlabel('wingbeat cycle')
ylabel('angular rate [deg/s]')
grid on

% plot spanwise hinge torque versus rotation angles
subplot(2,2,2)
hold on
plot(rot_steady-90,Mwing_springdamper_all_2Dref,'-k*','linewidth',2)
plot(rot_steady-90,Mwing_springdamper_all_2Dref_model,'-r*','linewidth',2)

plot(rot_steady(n_rot_ds_all(1),1)-90,Mwing_springdamper_all_2Dref(n_rot_ds_all(1),1),'ok','linewidth',2)
plot(rot_steady(n_rot_us_all(1),1)-90,Mwing_springdamper_all_2Dref(n_rot_us_all(1),1),'ok','linewidth',2)

grid on
axis tight
% axis square
axis([-75,75,-.25,.25])
set(gca,'xtick',[-75 0 75])
set(gca,'ytick',[-.25 0 .25])
xlabel('wing rotation angle [deg]')
ylabel('normalized spanwise hinge torque')

% plot spanwise hinge torque versus rotation rate
subplot(2,2,4)
hold on
plot(rot_dot_deg-90,Mwing_springdamper_all_2Dref,'-k','linewidth',2)
plot(rot_dot_deg-90,Mwing_springdamper_all_2Dref_model,'-r','linewidth',2)

plot(rot_dot_deg(n_rot_ds_all(1),1),Mwing_springdamper_all_2Dref(n_rot_ds_all(1),1),'ok','linewidth',2)
plot(rot_dot_deg(n_rot_us_all(1),1),Mwing_springdamper_all_2Dref(n_rot_us_all(1),1),'ok','linewidth',2)

grid on
axis tight
% axis square
axis([-2e5,2e5,-.25,.25])
set(gca,'xtick',[-2e5 0 2e5])
set(gca,'ytick',[-.25 0 .25])
xlabel('wing rotation angle [deg]')
ylabel('normalized spanwise hinge torque')

saveas(gca,['rotNrotrate_VS_Twing_STEADY_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.fig'])
saveas(gca,['rotNrotrate_VS_Twing_STEADY_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.png'])
plot2svg(['rotNrotrate_VS_Twing_STEADY_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.svg'])



%% spanwise hinge torque throughout wingbeat
figure
hold on
plot(t_norm,Mwing_springdamper_all_2Dref,'-k','linewidth',2)
plot(t_norm,Mwing_springdamper_all_2Dref_model,'--k','linewidth',2)
plot(t_norm,Mwing_inertiamass_all_2Dref,'-b','linewidth',2)
plot(t_norm,Mwing_aero_all_2Dref,'-r','linewidth',2)

plot(t_norm(n_rot_ds_all(1)),Mwing_springdamper_all_2Dref(n_rot_ds_all(1),1),'ok','linewidth',2)
plot(t_norm(n_rot_us_all(1)),Mwing_springdamper_all_2Dref(n_rot_us_all(1),1),'ok','linewidth',2)

grid on
axis tight
legend('damped spring','wing inertia','aerodynamic')
% axis square
axis([0,1,-.25,.25])
set(gca,'xtick',[0 1])
set(gca,'ytick',[-.25 0 .25])
xlabel('wingbeat cycle')
ylabel('normalized spanwise hinge torque')

saveas(gca,['Twing_vs_WBcycle_STEADY_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.fig'])
saveas(gca,['Twing_vs_WBcycle_STEADY_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.png'])
plot2svg(['Twing_vs_WBcycle_STEADY_peakloc',num2str(loc_set),'_steps',num2str(nr_timepoints),'.svg'])
