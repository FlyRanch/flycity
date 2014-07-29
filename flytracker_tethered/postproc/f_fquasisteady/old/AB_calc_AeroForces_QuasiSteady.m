% clear;
clc;
close all
warning off

addpath('/home/florian/Dropbox/WORK/flytracker/flytracker');
addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/mex/');
addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/core/');
addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/results/');

% load('WBdataset_all_steadyNmods_accelNorm.mat')
% load('WBdataset_steady_1603WBs.mat')
% load('BodyWingModel_WingSections_TotalMean.mat')
% load('borf_db_Fenhance_NOcali_alldata.mat')

%% constants
Mfly = const.Mfly;
Mg_fly = const.Mg_fly;
stroke_angle = const.stroke_angle;

wing_length = wing_model.length;
nr_sect = wing_model.nr_sect;
nr_timepoints = wing_model.nr_timepoints;

rot_lift_on = 1;
rot_lift_off = 0;

%% wingbeat kin

% Strokeplane2body rotation
R_strk = [cosd(stroke_angle) 0 sind(stroke_angle); 0 1 0; -sind(stroke_angle) 0 cosd(stroke_angle)];

% body non-moving/rotating
u_strk = [0 0 0];
w_strk = [0 0 0];

% calc steady WB kin
freq = nanmean(f_wb_steady);
stroke_coeffs = stroke_steady_fourier_coeffs_binmean;
dev_coeffs = dev_steady_fourier_coeffs_binmean;
rot_coeffs = pitch_steady_fourier_coeffs_binmean;

[stroke_L,stroke_R,stroke_dot_L,stroke_dot_R,stroke_dot_dot_L,stroke_dot_dot_R,...
    dev_L,dev_R,dev_dot_L,dev_dot_R,dev_dot_dot_L,dev_dot_dot_R,...
    rot_L,rot_R,rot_dot_L,rot_dot_R,rot_dot_dot_L,rot_dot_dot_R,t,dt,t_norm,dt_norm] = ...
    calc_WBkin_fouriercoeffs(nr_timepoints,freq,stroke_coeffs,dev_coeffs,rot_coeffs);

%% WB data
wb_loc(1) = 1;
wb_loc(2) = length(stroke_L);
down_loc(1) = 1;
down_loc(2) = find(stroke_L == max(stroke_L));
up_loc(1) = find(stroke_L == max(stroke_L))+1;
up_loc(2) = length(stroke_L);

[kine] = calc_angular_velocities_INCdots(stroke_L,dev_L,rot_L,stroke_R,dev_R,rot_R,stroke_dot_L,dev_dot_L,rot_dot_L,stroke_dot_R,dev_dot_R,rot_dot_R,stroke_dot_dot_L,dev_dot_dot_L,rot_dot_dot_L,stroke_dot_dot_R,dev_dot_dot_R,rot_dot_dot_R,freq,R_strk);

kine.R_strk = R_strk;
kine.u_strk(:,nr_timepoints) = u_strk;
kine.w_strk(:,nr_timepoints) = w_strk;

wb.wb_loc          = wb_loc;
wb.down_loc        = down_loc;
wb.up_loc          = up_loc;
wb.dt = dt;

%% NO rotational lift
[ FM_strkpln, FM_L, FM_R ,U_left, U_right, alfa_L, alfa_R, alfa_dot_L, alfa_dot_R ] =...
    Aerodynamic_forces( kine, body_model, wing_model, wb, rot_lift_off);

Fx_trans = FM_strkpln(1,:)./Mg_fly;
Fy_trans = FM_strkpln(2,:)./Mg_fly;
Fz_trans = FM_strkpln(3,:)./Mg_fly;

Mx_trans = FM_strkpln(4,:)./freq^2./(1e-3*wing_length)^5*1e-3;
My_trans = FM_strkpln(5,:)./freq^2./(1e-3*wing_length)^5*1e-3;
Mz_trans = FM_strkpln(6,:)./freq^2./(1e-3*wing_length)^5*1e-3;

F_trans = sqrt(Fx_trans.^2 + Fy_trans.^2 + Fz_trans.^2);

figure
subplot(3,2,1)
hold on
plot(t_norm,Fx_trans)

subplot(3,2,3)
hold on
plot(t_norm,Fy_trans)

subplot(3,2,5)
hold on
plot(t_norm,Fz_trans)

subplot(3,2,2)
hold on
plot(t_norm,Mx_trans)

subplot(3,2,4)
hold on
plot(t_norm,My_trans)

subplot(3,2,6)
hold on
plot(t_norm,Mz_trans)

%% INC rotational lift
[ FM_strkpln, FM_L, FM_R ,U_left, U_right, alfa_L, alfa_R, alfa_dot_L, alfa_dot_R ] = Aerodynamic_forces( kine, body_model, wing_model, wb, rot_lift_on);

Fx_transNrot = FM_strkpln(1,:)./Mg_fly;
Fy_transNrot = FM_strkpln(2,:)./Mg_fly;
Fz_transNrot = FM_strkpln(3,:)./Mg_fly;

Mx_transNrot = FM_strkpln(4,:)./freq^2./(1e-3*wing_length)^5*1e-3;
My_transNrot = FM_strkpln(5,:)./freq^2./(1e-3*wing_length)^5*1e-3;
Mz_transNrot = FM_strkpln(6,:)./freq^2./(1e-3*wing_length)^5*1e-3;

F_transNrot = sqrt(Fx_transNrot.^2 + Fy_transNrot.^2 + Fz_transNrot.^2);

subplot(3,2,1)
hold on
plot(t_norm,Fx_transNrot,'--r')
grid on

subplot(3,2,3)
hold on
plot(t_norm,Fy_transNrot,'--r')
grid on

subplot(3,2,5)
hold on
plot(t_norm,Fz_transNrot,'--r')
grid on

subplot(3,2,2)
hold on
plot(t_norm,Mx_transNrot,'--r')
grid on

subplot(3,2,4)
hold on
plot(t_norm,My_transNrot,'--r')
grid on

subplot(3,2,6)
hold on
plot(t_norm,Mz_transNrot,'--r')
grid on

%% borf data
loc = find(mod_value_all==0)

Fx_borf = Fx_all(:,3,loc)./Mg_fly;
Fy_borf = Fy_all(:,3,loc)./Mg_fly;
Fz_borf = Fz_all(:,3,loc)./Mg_fly;

Mx_borf = Mx_all(:,3,loc)./freq^2./(1e-3*wing_length)^5;
My_borf = My_all(:,3,loc)./freq^2./(1e-3*wing_length)^5;
Mz_borf = Mz_all(:,3,loc)./freq^2./(1e-3*wing_length)^5;

Fx_borf = Fx_borf(isnan(Fx_borf)==0);
Fy_borf = Fy_borf(isnan(Fy_borf)==0);
Fz_borf = Fz_borf(isnan(Fz_borf)==0);

Mx_borf = Mx_borf(isnan(Mx_borf)==0);
My_borf = My_borf(isnan(My_borf)==0);
Mz_borf = Mz_borf(isnan(Mz_borf)==0);

t_borf = [0:1/(length(Fx_borf)-1):1];

subplot(3,2,1)
plot(t_borf,Fx_borf,'c')
xlabel('t*','fontsize',10) 
ylabel('Fx','fontsize',10)

subplot(3,2,3)
plot(t_borf,-Fy_borf,'c')
xlabel('t*','fontsize',10) 
ylabel('Fy','fontsize',10)

subplot(3,2,5)
plot(t_borf,-Fz_borf,'c')
xlabel('t*','fontsize',10) 
ylabel('Fz','fontsize',10)

subplot(3,2,2)
plot(t_borf,Mx_borf,'c')
xlabel('t*','fontsize',10) 
ylabel('Mx','fontsize',10)

subplot(3,2,4)
plot(t_borf,-My_borf,'c')
xlabel('t*','fontsize',10) 
ylabel('My','fontsize',10)

subplot(3,2,6)
plot(t_borf,-Mz_borf,'c')
xlabel('t*','fontsize',10) 
ylabel('Mz','fontsize',10)

saveas(gca,'FnM_quasisteadyNborf_steadyWB.png')

%% plot Ftotal

figure
plot(t_norm,F_trans)
hold on
plot(t_norm,F_transNrot,'--r')
grid on

F_borf = sqrt(Fx_borf.^2 + Fy_borf.^2 + Fz_borf.^2);
plot(t_borf,F_borf,'c')
xlabel('t*','fontsize',10) 
ylabel('Ftotal','fontsize',10)

saveas(gca,'Ftotal_quasisteadyNborf_steadyWB.png')

%% plot alpha & U
figure
subplot(2,2,1)
plot(t_norm,alfa_L')
hold on
plot(t_norm,alfa_R')
xlabel('t*','fontsize',10) 
ylabel('AoA','fontsize',10)

subplot(2,2,2)
ux_L(:,:) = U_left(1,:,:);
ux_R(:,:) = U_right(1,:,:);
plot(t_norm,ux_L)
hold on
plot(t_norm,ux_R)
xlabel('t*','fontsize',10) 
ylabel('Ux','fontsize',10)

subplot(2,2,3)
uz_L(:,:) = U_left(3,:,:);
uz_R(:,:) = U_right(3,:,:);
plot(t_norm,uz_L)
hold on
plot(t_norm,uz_R)
xlabel('t*','fontsize',10) 
ylabel('Uz','fontsize',10)

subplot(2,2,4)
U_L = sqrt(ux_L.^2 + uz_L.^2);
U_R = sqrt(ux_R.^2 + uz_R.^2);
plot(t_norm,U_L)
hold on
plot(t_norm,U_R)
xlabel('t*','fontsize',10) 
ylabel('Utot','fontsize',10)

saveas(gca,'AoAnU_quasisteadyModel_steadyWB.png')



