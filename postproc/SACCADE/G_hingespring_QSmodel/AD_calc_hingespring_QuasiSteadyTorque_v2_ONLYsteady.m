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
nr_sect = settings.nr_chord_sect;
nr_timepoints = settings.nr_timepoints;

% nr_timepoints = 500;
% nr_timepoints = 5000;

Mg_fly = body_model.Mg_fly;                             % [kg]
l_wing = wing_model.length*1e-3;                        % [m]
area_wing = wing_model.area*1e-6;                       % [m^2]
massperArea = wing_model.mass/(wing_model.area*1e-6);   % [kg/m^2]


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

%% steady WB

% freq & time scales
f_steady = f_wb_steady_meanCIstd(1);

t_norm = 0:1/(nr_timepoints-1):1;
t = t_norm/f_steady;
dt = (1/f_steady)/(nr_timepoints-1);  

% kinematics
stroke_coeffs_steady = stroke_steady_fourier_coeffs_binmean;
dev_coeffs_steady = dev_steady_fourier_coeffs_binmean;
rot_coeffs_steady = pitch_steady_fourier_coeffs_binmean;

[stroke_steady] = (calc_val_fourier_series_4thN8th_order(t_norm,stroke_coeffs_steady,0))';
[dev_steady] = (calc_val_fourier_series_4thN8th_order(t_norm,dev_coeffs_steady,0))';
[rot_steady] = (calc_val_fourier_series_4thN8th_order(t_norm,rot_coeffs_steady,0))';

[stroke_dot_steady] = (calc_val_dot_fourier_series_4thN8th_order(t_norm,stroke_coeffs_steady,0))'*f_steady;
[dev_dot_steady] = (calc_val_dot_fourier_series_4thN8th_order(t_norm,dev_coeffs_steady,0))'*f_steady;
[rot_dot_steady] = (calc_val_dot_fourier_series_4thN8th_order(t_norm,rot_coeffs_steady,0))'*f_steady;

[stroke_dot_dot_steady] = (calc_val_dot_dot_fourier_series_4thN8th_order(t_norm,stroke_coeffs_steady,0))'*f_steady^2;
[dev_dot_dot_steady] = (calc_val_dot_dot_fourier_series_4thN8th_order(t_norm,dev_coeffs_steady,0))'*f_steady^2;
[rot_dot_dot_steady] = (calc_val_dot_dot_fourier_series_4thN8th_order(t_norm,rot_coeffs_steady,0))'*f_steady^2;
% 
% stroke_dot_steady2 = gradient(stroke_steady,dt);
% stroke_dot_dot_steady2 = gradient(stroke_dot_steady,dt);
% 
% dev_dot_steady2 = gradient(dev_steady,dt);
% dev_dot_dot_steady2 = gradient(dev_dot_steady,dt);
% 
% rot_dot_steady2 = gradient(rot_steady,dt);
% rot_dot_dot_steady2 = gradient(rot_dot_steady,dt);

%% qs forces & torques
wing_kin.freq = f_steady;

wing_kin.stroke_L = deg2rad(stroke_steady);
wing_kin.rot_L = deg2rad(rot_steady);
wing_kin.dev_L = deg2rad(dev_steady);

wing_kin.stroke_dot_L = deg2rad(stroke_dot_steady);
wing_kin.rot_dot_L = deg2rad(rot_dot_steady);
wing_kin.dev_dot_L = deg2rad(dev_dot_steady);

wing_kin.stroke_dot_dot_L = deg2rad(stroke_dot_dot_steady);
wing_kin.rot_dot_dot_L = deg2rad(rot_dot_dot_steady);
wing_kin.dev_dot_dot_L = deg2rad(dev_dot_dot_steady);

wing_kin.stroke_R = deg2rad(stroke_steady);
wing_kin.rot_R = deg2rad(rot_steady);
wing_kin.dev_R = deg2rad(dev_steady);

wing_kin.stroke_dot_R = deg2rad(stroke_dot_steady);
wing_kin.rot_dot_R = deg2rad(rot_dot_steady);
wing_kin.dev_dot_R = deg2rad(dev_dot_steady);

wing_kin.stroke_dot_dot_R = deg2rad(stroke_dot_dot_steady);
wing_kin.rot_dot_dot_R = deg2rad(rot_dot_dot_steady);
wing_kin.dev_dot_dot_R = deg2rad(dev_dot_dot_steady);

[ FM_strkpln, kine ] = quasi_steady_FnMnWingkin_atTimeNspanSections_fourier( settings, body_model, wing_model, wing_kin, rot_on );
    
    %% calc spanwise wing torque in wing ref frame (spanwise y-axis positive towards tip): ds:pitch up pos

    % inertial torque
    rot_dot_dot_rad_wingref = kine.eta_ddot_R;
    Mwing_inertia  = rot_dot_dot_rad_wingref * Itot / Mg_fly / l_wing;

    % aero torque
    Maero_sect = FM_strkpln.FM_R.My_R;
    Mwing_aero  = nansum(Maero_sect)' / Mg_fly / (l_wing*1e3);
    
    % mass torque
    clear az Mmass_sect
    az(:,:) = kine.Udot_left(3,:,:)*1e-3;      % [ m/s^2 ]
    x_cm = (x_sect-c_sect/2);
    for n = 1:size(az,1)
        Mmass_sect(:,n) = mass_tot_sect .* x_cm .*az(n,:)';
    end
    Mwing_mass  = nansum(Mmass_sect)' / Mg_fly / l_wing;
    
    %% sum of torques equals angular acceleration: I*rot_ddot + C*rot_dot + k*(rot-rot0) = Taero + Tmass
    % hing torque = - wing torque
    Mwing_springdamper  = -Mwing_inertia + Mwing_mass + Mwing_aero;
    Mwing_inertiamass   = -Mwing_inertia + Mwing_mass;
    
    figure
    hold on
    plot(t_norm,Mwing_springdamper,'k')
    plot(t_norm,-Mwing_inertia,'b')
    plot(t_norm,+Mwing_aero,'r')
    plot(t_norm,+Mwing_mass,'g')
    grid on
    
    %% torque at min & max rotation angle
    [rot_maxpeaks,rot_maxpeaklocs] = findpeaks(rot_steady);
    if loc_set == 3
        n_rot_us = rot_maxpeaklocs(maxpeak_nr);
        rot_us = rot_steady(n_rot_us);
    else
        n_rot_ds = rot_maxpeaklocs(maxpeak_nr);
        rot_ds = rot_steady(n_rot_ds);
    end
    
    [rot_minpeaks,rot_minpeaklocs] = findpeaks(-rot_steady);
    if loc_set == 3
        n_rot_ds = rot_minpeaklocs(minpeak_nr);
        rot_ds = rot_steady(n_rot_ds);
    else
        n_rot_us = rot_minpeaklocs(minpeak_nr);
        rot_us = rot_steady(n_rot_us);
    end

    %% store data
    % spanwise torques @ rot ds & rot us
    Mwing_springdamper_atRotds = Mwing_springdamper(n_rot_ds);
    Mwing_springdamper_atRotus = Mwing_springdamper(n_rot_us);

    Mwing_inertiamass_atRotds = Mwing_inertiamass(n_rot_ds);
    Mwing_inertiamass_atRotus = Mwing_inertiamass(n_rot_us);

    Mwing_inertia_atRotds = Mwing_inertia(n_rot_ds);
    Mwing_inertia_atRotus = Mwing_inertia(n_rot_us);

    Mwing_mass_atRotds = Mwing_mass(n_rot_ds);
    Mwing_mass_atRotus = Mwing_mass(n_rot_us);

    Mwing_aero_atRotds = Mwing_aero(n_rot_ds);
    Mwing_aero_atRotus = Mwing_aero(n_rot_us);
    
    %% torques in 2D wing ref frame: NEGATIVE 3D WING REF FRAME
    Mwing_springdamper_2Dref = -Mwing_springdamper;
    Mwing_springdamper_atRotds_2Dref = -Mwing_springdamper_atRotds;
    Mwing_springdamper_atRotus_2Dref = -Mwing_springdamper_atRotus;
    
    Mwing_inertiamass_2Dref = -Mwing_inertiamass;
    Mwing_inertiamass_atRotds_2Dref = -Mwing_inertiamass_atRotds;
    Mwing_inertiamass_atRotus_2Dref = -Mwing_inertiamass_atRotus;
    
    Mwing_inertia_2Dref = -Mwing_inertia;
    Mwing_inertia_atRotds_2Dref = -Mwing_inertia_atRotds;
    Mwing_inertia_atRotus_2Dref = -Mwing_inertia_atRotus;

    Mwing_mass_2Dref = -Mwing_mass;
    Mwing_mass_atRotds_2Dref = -Mwing_mass_atRotds;
    Mwing_mass_atRotus_2Dref = -Mwing_mass_atRotus;
    
    Mwing_aero_2Dref = -Mwing_aero;
    Mwing_aero_atRotds_2Dref = -Mwing_aero_atRotds;
    Mwing_aero_atRotus_2Dref = -Mwing_aero_atRotus;
    
    %% calc winghinge rotational springdamper coeffs
    p = polyfitn([rot_dot_steady,rot_steady],Mwing_springdamper_2Dref,1);
    C_norm_deg_model = p.Coefficients(1);
    k_norm_deg_model = p.Coefficients(2);
    rot0_norm_deg_model = -p.Coefficients(3)/p.Coefficients(2);
    Mwing_springdamper_2Dref_model = C_norm_deg_model*rot_dot_steady + k_norm_deg_model*(rot_steady-rot0_norm_deg_model);
    Mwing_spring_2Dref_model = k_norm_deg_model*(rot_steady-rot0_norm_deg_model);
    Mwing_damper_2Dref_model = C_norm_deg_model*rot_dot_steady;

    figure
    hold on
    plot(t_norm,Mwing_springdamper_2Dref)
    plot(t_norm,Mwing_springdamper_2Dref_model,'r')
    plot(t_norm,Mwing_spring_2Dref_model,'c')
    plot(t_norm,Mwing_damper_2Dref_model,'m')
    
    figure
    hold on
    plot(rot_steady,Mwing_springdamper_2Dref)
    plot(rot_steady,Mwing_springdamper_2Dref_model,'r')

    %% save data & plots
    save(['hingespringdamper_coeffs_STEADY_steps',num2str(nr_timepoints),'.mat'])
    
    mkdir('qsModel_FnM_k_alpha0_n_rotVsTspanwise')
    cd('qsModel_FnM_k_alpha0_n_rotVsTspanwise')
  
    
%% plot kine
figure
subplot(3,1,1)
hold on
plot(t_norm,stroke_steady)
plot(t_norm,dev_steady,'r')
plot(t_norm,rot_steady-90,'g')

axis([0,1,-75,75])
set(gca,'xtick',[0 1])
set(gca,'ytick',[-75 0 75])
xlabel('wingbeat cycle')
ylabel('angles [deg]')
grid on

subplot(3,1,2)
hold on
plot(t_norm,stroke_dot_steady)
plot(t_norm,dev_dot_steady,'r')
plot(t_norm,rot_dot_steady,'g')

axis([0,1,-1.5e5,1.5e5])
set(gca,'xtick',[0 1])
set(gca,'ytick',[-1.5e5 0 1.5e5])
xlabel('wingbeat cycle')
ylabel('angle rate [deg/s]')
grid on

subplot(3,1,3)
hold on
plot(t_norm,stroke_dot_dot_steady)
plot(t_norm,dev_dot_dot_steady,'r')
plot(t_norm,rot_dot_dot_steady,'g')

axis([0,1,-8e8,8e8])
set(gca,'xtick',[0 1])
set(gca,'ytick',[-8e8 0 8e8])
xlabel('wingbeat cycle')
ylabel('angle accel [deg/s^2]')
grid on

saveas(gca,['wingkin_steps',num2str(nr_timepoints),'.fig'])
saveas(gca,['wingkin_steps',num2str(nr_timepoints),'.png'])
plot2svg(['wingkin_steps',num2str(nr_timepoints),'.svg'])

%% plot rotation angles throughout wingbeat
figure
subplot(2,2,1)
hold on
subplot(2,2,2)
hold on
subplot(2,2,3)
hold on
subplot(2,2,4)
hold on

for i = 1:length(t_norm)-1
    
    plot(t_norm(i:i+1),rot_steady(i:i+1),'-k','linewidth',2)
end
    
plot(t_norm(n_rot_ds),rot_steady(n_rot_ds,1)-90,'or','linewidth',2)
plot(t_norm(n_rot_us),rot_steady(n_rot_us,1)-90,'ob','linewidth',2)

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
plot(t_norm,rot_dot_steady,'-k','linewidth',2)
plot(t_norm(n_rot_ds),rot_dot_steady(n_rot_ds,1)-90,'ok','linewidth',2)
plot(t_norm(n_rot_us),rot_dot_steady(n_rot_us,1)-90,'ok','linewidth',2)

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
plot(rot_steady-90,Mwing_springdamper_2Dref,'-k','linewidth',2)
plot(rot_steady-90,Mwing_springdamper_2Dref_model,'-r','linewidth',2)

plot(rot_steady(n_rot_ds(1),1)-90,Mwing_springdamper_2Dref(n_rot_ds(1),1),'ok','linewidth',2)
plot(rot_steady(n_rot_us(1),1)-90,Mwing_springdamper_2Dref(n_rot_us(1),1),'ok','linewidth',2)

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
plot(rot_dot_steady-90,Mwing_springdamper_2Dref,'-k','linewidth',2)
plot(rot_dot_steady-90,Mwing_springdamper_2Dref_model,'-r','linewidth',2)

plot(rot_dot_steady(n_rot_ds(1),1),Mwing_springdamper_2Dref(n_rot_ds(1),1),'ok','linewidth',2)
plot(rot_dot_steady(n_rot_us(1),1),Mwing_springdamper_2Dref(n_rot_us(1),1),'ok','linewidth',2)

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
plot(t_norm,Mwing_springdamper_2Dref,'-k','linewidth',2)
plot(t_norm,Mwing_springdamper_2Dref_model,'--k','linewidth',2)
plot(t_norm,Mwing_inertiamass_2Dref,'-b','linewidth',2)
plot(t_norm,Mwing_aero_2Dref,'-r','linewidth',2)

plot(t_norm(n_rot_ds(1)),Mwing_springdamper_2Dref(n_rot_ds(1),1),'ok','linewidth',2)
plot(t_norm(n_rot_us(1)),Mwing_springdamper_2Dref(n_rot_us(1),1),'ok','linewidth',2)

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
