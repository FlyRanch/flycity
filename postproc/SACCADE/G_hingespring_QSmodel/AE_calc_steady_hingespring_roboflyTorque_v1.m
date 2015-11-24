clear;
clc;
% close all
warning off

load('bodyNwingModel_4qsModel.mat')

loadname_kin_steady=dir('WBdataset_steady_1603WBs.mat')
loadname_kin_steady = loadname_kin_steady.name;
load(loadname_kin_steady)

loadname_robofly_steady=dir('force_torque_steady_robofly.mat')
loadname_robofly_steady = loadname_robofly_steady.name;
load(loadname_robofly_steady)

%% constants
rot_on =1;

nr_sect = settings.nr_chord_sect;
nr_timepoints = settings.nr_timepoints;

% nr_timepoints = 500;
% nr_timepoints = 5000;

% cyclic colormap
jet_wrap = vertcat(jet(nr_timepoints/2),flipud(jet(nr_timepoints/2)));
jet_map = jet(nr_timepoints);

%% wing params
rho_air = 1.225;                            % [kg/m^3]

Mg_fly = body_model.Mg_fly;                             % [kg]
l_wing = wing_model.length*1e-3;                        % [m]
area_wing = wing_model.area*1e-6;                       % [m^2]
massperArea = wing_model.mass/(wing_model.area*1e-6);   % [kg/m^2]

h_sect = wing_model.length/nr_sect*1e-3;    % [m]
x_sect = wing_model.x_LE_R'*1e-3;           % [m]
y_sect = wing_model.y_sect_R(:,2)*1e-3;     % [m]
c_sect = wing_model.chords_R'*1e-3;         % [m]

c_mean = mean(c_sect);                      % [m]
c_max = max(c_sect);

%% steady WB
% freq & time scales
f_steady = f_wb_steady_meanCIstd(1);

t_norm = 0:1/(nr_timepoints-1):1;
t_norm = t_norm';
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

%% robofly forces & torques
% interpolate
t_steady_interp = t_steady;
t_steady_interp(end) = t(end);
Fx_norm_robofly_world = interp1(t_steady_interp,Fx_norm_steady,t);
Fy_norm_robofly_world = interp1(t_steady_interp,Fy_norm_steady,t);
Fz_norm_robofly_world = interp1(t_steady_interp,Fz_norm_steady,t);

% rotate F around stroke angle, deviation angle & rotation angle
for i = 1:length(t)
    stroke_now = stroke_steady(i);
    dev_now = dev_steady(i);
    rot_now = rot_steady(i);
    F_now = [Fx_norm_robofly_world(i); Fy_norm_robofly_world(i); Fy_norm_robofly_world(i)];
    
    F_wing_now =    inv([cosd(rot_now) 0 sind(rot_now); 0 1 0; -sind(rot_now) 0 cosd(rot_now)]) *...
                    inv([1 0 0; 0 cosd(dev_now) -sind(dev_now); 0 sind(dev_now) cosd(dev_now)]) *...
                    inv([cosd(stroke_now) -sind(stroke_now) 0; sind(stroke_now) cosd(stroke_now) 0; 0 0 1]) *...
                    F_now;

    Fx_norm_robofly(i,1) = F_wing_now(1);
    Fy_norm_robofly(i,1) = F_wing_now(2);
    Fz_norm_robofly(i,1) = F_wing_now(3);
end

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

[ FM_qs_strkpln, kine ] = quasi_steady_FnMnWingkin_atTimeNspanSections_fourier( settings, body_model, wing_model, wing_kin, rot_on );

%% aero torque

% % weighted mean aero center (ac) @ 1/4 cord
% x_ca = (x_sect-c_sect/4);
% x_ca_norm = x_ca / l_wing;
% x_ca_norm_mean = mean(x_ca_norm);
% x_ca_norm_weighted = mean(x_ca .* c_sect) / mean(c_sect);

% weighted mean aero center
alfa = kine.alfa_R(1,:)';
for i=1:length(rot_steady)
    c_ratio_now = 0.82*abs(alfa(i))/pi+0.05;
    
    % aero center - rot axis
    x_ca_norm_now = (c_sect*c_ratio_now-x_sect)  / l_wing;
    x_ca_norm_weighted(i,1) = mean(x_ca_norm_now .* c_sect) / mean(c_sect);
    
    % aero center - LEADING EDGE
    x_caLE_norm_now = (c_sect*c_ratio_now)  / l_wing;
    x_caLE_norm_weighted(i,1) = mean(x_caLE_norm_now .* c_sect) / mean(c_sect);
end

% robofly aero torque
My_norm_robofly = Fz_norm_robofly .* x_ca_norm_weighted;
My_norm_robofly_LE = Fz_norm_robofly .* x_caLE_norm_weighted;

% qs aero torque
Fz_norm_qs_world    = FM_qs_strkpln.FM(:,3)/ Mg_fly;
Fz_norm_qs          = nansum(FM_qs_strkpln.FM_R.Fz_R)'/ Mg_fly;
My_norm_qs          = nansum(FM_qs_strkpln.FM_R.My_R)'  / Mg_fly / l_wing / 1e3;

x_ca_norm_qs        = My_norm_qs./Fz_norm_qs;
x_ca_norm_mean_qs   = mean(x_ca_norm_qs);

My_norm_qs2         = Fz_norm_qs .* x_ca_norm_weighted;
My_norm_qs2_LE      = Fz_norm_qs .* x_caLE_norm_weighted;

% !!! use qs model torque !!!
Mwing_aero_norm = My_norm_robofly;
Mwing_aero_norm_LE = My_norm_robofly_LE;

%% mass torque
% wing section mass
mass_wing_sect = massperArea * h_sect * c_sect;
mass_virtual_sect = 1/4* pi * rho_air * h_sect * c_sect.^2;
mass_tot_sect = mass_wing_sect + mass_virtual_sect;

% wing section torque
clear az Mmass_sect
az(:,:) = kine.Udot_left(3,:,:)*1e-3;      % [ m/s^2 ]
x_cm = (x_sect-c_sect/2);
x_cm_LE = (-c_sect/2);
for n = 1:size(az,1)
    Mmass_sect(:,n) = mass_tot_sect .* x_cm .*az(n,:)';
    Mmass_sect_LE(:,n) = mass_tot_sect .* x_cm_LE .*az(n,:)';
end
Mwing_mass_norm  = nansum(Mmass_sect)' / Mg_fly / l_wing;
Mwing_mass_norm_LE  = nansum(Mmass_sect_LE)' / Mg_fly / l_wing;

%% inertial torque
% Iwing
Iwing_sect = 1/12 * massperArea * h_sect * c_sect.^3 + massperArea * h_sect * c_sect .* (1/2*c_sect - x_sect).^2;
Itot = nansum(Iwing_sect);

Iwing_sect_LE = 1/12 * massperArea * h_sect * c_sect.^3 + massperArea * h_sect * c_sect .* (1/2*c_sect).^2;
Itot_LE = nansum(Iwing_sect_LE);

% %Iwing & Iaddedmass
% Iwing_sect = 1/12 * massperArea * h_sect * c_sect.^3 + massperArea * h_sect * c_sect .* (1/2*c_sect - x_sect).^2;
% Iaddedmass_sect = 1/64 * pi * rho_air * h_sect * c_sect.^4 + 1/4* pi * rho_air * h_sect * c_sect.^2 .* (1/2*c_sect - x_sect).^2;
% Itot_sect = Iwing_sect + Iaddedmass_sect;
% 
% Iwing = nansum(Iwing_sect);
% Iaddedmass = nansum(Iaddedmass_sect);
% Itot = nansum(Itot_sect);

% inertial torque (!!! in wing reference frame !!!)
rot_dot_dot_rad_wingref = kine.eta_ddot_R;
% rot_dot_dot_rad_wingref = -deg2rad(rot_dot_dot_steady);
Mwing_inertia_norm  = rot_dot_dot_rad_wingref' * Itot / Mg_fly / l_wing;
Mwing_inertia_norm_LE  = rot_dot_dot_rad_wingref' * Itot_LE / Mg_fly / l_wing;

%% sum of torques equals angular acceleration: I*rot_ddot + C*rot_dot + k*(rot-rot0) = Taero + Tmass
% hing torque = - wing torque
Mwing_springdamper_norm     = -Mwing_inertia_norm + Mwing_mass_norm + Mwing_aero_norm;
Mwing_inertiamass_norm      = -Mwing_inertia_norm + Mwing_mass_norm;

Mwing_inertiamass_norm_LE   = -Mwing_inertia_norm_LE + Mwing_mass_norm_LE;
Mwing_springdamper_norm_LE  = -Mwing_inertia_norm_LE + Mwing_mass_norm_LE + Mwing_aero_norm_LE;

%% torques in 2D wing ref frame: NEGATIVE 3D WING REF FRAME
Mwing_springdamper_2Dref = -Mwing_springdamper_norm;
Mwing_aero_2Dref = -Mwing_aero_norm;
Mwing_inertiamass_2Dref = -Mwing_inertiamass_norm;
Mwing_inertia_2Dref = -Mwing_inertia_norm;
Mwing_mass_2Dref = -Mwing_mass_norm;

Mwing_springdamper_2Dref_LE = -Mwing_springdamper_norm_LE;
Mwing_aero_2Dref_LE = -Mwing_aero_norm_LE;
Mwing_inertiamass_2Dref_LE = -Mwing_inertiamass_norm_LE;
Mwing_inertia_2Dref_LE = -Mwing_inertia_norm_LE;
Mwing_mass_2Dref_LE = -Mwing_mass_norm_LE;

%% calc winghinge rotational springdamper coeffs
p = polyfitn([rot_dot_steady',rot_steady'],Mwing_springdamper_2Dref,1);
C_norm_deg_model = p.Coefficients(1);
k_norm_deg_model = p.Coefficients(2);
rot0_norm_deg_model = -p.Coefficients(3)/p.Coefficients(2);
Mwing_springdamper_2Dref_model = C_norm_deg_model*rot_dot_steady + k_norm_deg_model*(rot_steady-rot0_norm_deg_model);
Mwing_spring_2Dref_model = k_norm_deg_model*(rot_steady-rot0_norm_deg_model);
Mwing_damper_2Dref_model = C_norm_deg_model*rot_dot_steady;

p_LE = polyfitn([rot_dot_steady',rot_steady'],Mwing_springdamper_2Dref_LE,1);
C_norm_deg_model_LE = p_LE.Coefficients(1);
k_norm_deg_model_LE = p_LE.Coefficients(2);
rot0_norm_deg_model_LE = -p_LE.Coefficients(3)/p_LE.Coefficients(2);
Mwing_springdamper_2Dref_model_LE = C_norm_deg_model_LE*rot_dot_steady + k_norm_deg_model_LE*(rot_steady-rot0_norm_deg_model_LE);
Mwing_spring_2Dref_model_LE = k_norm_deg_model_LE*(rot_steady-rot0_norm_deg_model_LE);
Mwing_damper_2Dref_model_LE = C_norm_deg_model_LE*rot_dot_steady;

%% save data & plots
save(['hingespringdamper_coeffs_STEADY_steps',num2str(nr_timepoints),'.mat'])

mkdir('qsModel_FnM_k_alpha0_n_rotVsTspanwise')
cd('qsModel_FnM_k_alpha0_n_rotVsTspanwise')

%% plot aero forces & torques    
figure
subplot(2,2,1)
hold on
plot(t_norm,Fz_norm_robofly_world)
plot(t_norm,Fz_norm_qs_world,'r')

subplot(2,2,2)
hold on
plot(t_norm,Fz_norm_robofly)
plot(t_norm,Fz_norm_qs,'r')

subplot(2,2,3)
hold on
plot(t_norm,My_norm_robofly)
plot(t_norm,My_norm_qs,'r')
plot(t_norm,My_norm_qs2,'g')

plot(t_norm,My_norm_robofly_LE,'c')
plot(t_norm,My_norm_qs2_LE,'m')

subplot(2,2,4)
hold on
plot(t_norm,x_ca_norm_weighted)
plot(t_norm,x_caLE_norm_weighted,'c')

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

%% plot rotation angles & torques throughout wingbeat
figure
subplot(2,3,1)
hold on
subplot(2,3,2)
hold on
subplot(2,3,3)
hold on
subplot(2,3,4)
hold on
subplot(2,3,5)
hold on
subplot(2,3,6)
hold on

for i = 1:length(t_norm)-1
    subplot(2,3,1)
    plot(t_norm(i:i+1),rot_steady(i:i+1),'-','color',jet_map(i,:),'linewidth',2)

    subplot(2,3,4)
    plot(t_norm(i:i+1),rot_dot_steady(i:i+1),'-','color',jet_map(i,:),'linewidth',2)

    subplot(2,3,2)
    plot(t_norm(i:i+1),Mwing_springdamper_2Dref(i:i+1),'-','color',jet_map(i,:),'linewidth',2)

    subplot(2,3,5)
    plot(t_norm(i:i+1),Mwing_springdamper_2Dref_LE(i:i+1),'-','color',jet_map(i,:),'linewidth',2)

    subplot(2,3,3)
    plot(rot_steady(i:i+1),Mwing_springdamper_2Dref(i:i+1),'-','color',jet_map(i,:),'linewidth',2)

    subplot(2,3,6)
    plot(rot_steady(i:i+1),Mwing_springdamper_2Dref_LE(i:i+1),'-','color',jet_map(i,:),'linewidth',2)
end
    subplot(2,3,2)
    plot(t_norm,Mwing_springdamper_2Dref_model,'k','linewidth',2)
    
    subplot(2,3,5)
    plot(t_norm,Mwing_springdamper_2Dref_model_LE,'k','linewidth',2)

    subplot(2,3,3)
    plot(rot_steady,Mwing_springdamper_2Dref_model,'k','linewidth',2)
    
    subplot(2,3,6)
    plot(rot_steady,Mwing_springdamper_2Dref_model_LE,'k','linewidth',2)

subplot(2,3,1)
axis tight
% axis square
axis([0,1,0,180])
set(gca,'xtick',[0 1])
set(gca,'ytick',[0 90 180])
xlabel('wingbeat cycle')
ylabel('wing rotation angle [deg]')
grid on

subplot(2,3,4)
axis tight
% axis square
axis([0,1,-2e5,2e5])
set(gca,'xtick',[0 1])
set(gca,'ytick',[-2e5 0 2e5])
xlabel('wingbeat cycle')
ylabel('angular rate [deg/s]')
grid on

subplot(2,3,2)
axis tight
% axis square
axis([0,1,-.25,.25])
set(gca,'xtick',[0 90 180])
set(gca,'ytick',[-.25 0 .25])
xlabel('wing rotation angle [deg]')
ylabel('normalized hinge torque')
grid on

subplot(2,3,5)
axis tight
% axis square
axis([0,1,-.25,.25])
set(gca,'xtick',[0 1])
set(gca,'ytick',[-.25 0 .25])
xlabel('wingbeat cycle')
ylabel('normalized LE torque')
grid on

subplot(2,3,3)
axis tight
% axis square
axis([0,180,-.25,.25])
set(gca,'xtick',[0 1])
set(gca,'ytick',[-.25 0 .25])
xlabel('wingbeat cycle')
ylabel('normalized hinge torque')
grid on

subplot(2,3,6)
axis tight
% axis square
axis([0,180,-.25,.25])
set(gca,'xtick',[0 90 180])
set(gca,'ytick',[-.25 0 .25])
xlabel('wing rotation angle [deg]')
ylabel('normalized LE torque')
grid on

saveas(gca,['rotNrotrate_VS_Twing_hingeAndLErotation_steps',num2str(nr_timepoints),'.fig'])
saveas(gca,['rotNrotrate_VS_Twing_hingeAndLErotation_steps',num2str(nr_timepoints),'.png'])
plot2svg(['rotNrotrate_VS_Twing_hingeAndLErotation_steps',num2str(nr_timepoints),'.svg'])

%% spanwise hinge torque throughout wingbeat
figure
subplot(2,1,1)
hold on
plot(t_norm,Mwing_springdamper_2Dref,'-k','linewidth',2)
plot(t_norm,-Mwing_inertia_2Dref,'-b','linewidth',2)
plot(t_norm,Mwing_mass_2Dref,'-g','linewidth',2)
plot(t_norm,Mwing_aero_2Dref,'-r','linewidth',2)

grid on
axis tight
legend('damped spring','wing inertia','mass','aerodynamic','location','northwest')
% axis square
axis([0,1,-.25,.25])
set(gca,'xtick',[0 1])
set(gca,'ytick',[-.25 0 .25])
xlabel('wingbeat cycle')
ylabel('normalized spanwise hinge torque')

subplot(2,1,2)
hold on
plot(t_norm,Mwing_springdamper_2Dref_LE,'-k','linewidth',2)
plot(t_norm,-Mwing_inertia_2Dref_LE,'-b','linewidth',2)
plot(t_norm,Mwing_mass_2Dref_LE,'-g','linewidth',2)
plot(t_norm,Mwing_aero_2Dref_LE,'-r','linewidth',2)

grid on
axis tight
legend('damped spring','wing inertia','aerodynamic','model')
% axis square
axis([0,1,-.25,.25])
set(gca,'xtick',[0 1])
set(gca,'ytick',[-.25 0 .25])
xlabel('wingbeat cycle')
ylabel('normalized spanwise LE torque')

saveas(gca,['Twing_vs_WBcycle_hingeAndLErotation_steps',num2str(nr_timepoints),'.fig'])
saveas(gca,['Twing_vs_WBcycle_hingeAndLErotation_steps',num2str(nr_timepoints),'.png'])
plot2svg(['Twing_vs_WBcycle_hingeAndLErotation_steps',num2str(nr_timepoints),'.svg'])

cd ..
