% clear;
clc;
close all
warning off

addpath('/home/florian/Dropbox/WORK/flytracker/flytracker');
addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/mex/');
addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/core/');
addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/results/');

load('BodyWingModel_WingSections_TotalMean.mat')
load('WBdataset_steady_1603WBs.mat')
load('WBmod_accelbased_RollAccel_582WBs.mat')
load('norm_data.mat')

plot_on = 1
plot_on = 0
if plot_on == 1
    mkdir('ModRollAccel')
    cd('ModRollAccel')
end

%% constants
Mfly = const.Mfly;
Mg_fly = const.Mg_fly;
stroke_angle = const.stroke_angle;

wing_length = wing_model.length;
nr_sect = wing_model.nr_sect;
nr_timepoints = wing_model.nr_timepoints;

rot_lift_on = 1;
rot_lift_off = 0;

R_strk = body_model.R_strk;

% steady WB fouriers
stroke_coeffs_steady = stroke_steady_fourier_coeffs_binmean;
dev_coeffs_steady = dev_steady_fourier_coeffs_binmean;
rot_coeffs_steady = pitch_steady_fourier_coeffs_binmean;

%% body dynamics
% NO body vel&rot
body_model.u_strk(:,nr_timepoints) = [0 0 0];
body_model.w_strk(:,nr_timepoints) = [0 0 0];

%% MOD values
mod_min = -.5;
mod_max = 2.5;
dmod = .125;
mod_value = [mod_min:dmod:mod_max]';

for i = 1:length(mod_value)
    mod_value_now = mod_value(i);
    
    %% ALLnoFREQ MODS
    % set freq FIRST (same for MODnSteady)
    freq = f_wb_steady;

    % WBmod fouriers
    stroke_coeffs_MOD_L = strokeMOD_down_RollAccel_fourier_coeffs_binmean;
    dev_coeffs_MOD_L = devMOD_down_RollAccel_fourier_coeffs_binmean;
    rot_coeffs_MOD_L = pitchMOD_down_RollAccel_fourier_coeffs_binmean;
    
    stroke_coeffs_MOD_R = strokeMOD_up_RollAccel_fourier_coeffs_binmean;
    dev_coeffs_MOD_R = devMOD_up_RollAccel_fourier_coeffs_binmean;
    rot_coeffs_MOD_R = pitchMOD_up_RollAccel_fourier_coeffs_binmean;
    
    % calc FnM SYMMETRIC WB
    [FM_trans,FM_trans_norm,FM_transNrot,FM_transNrot_norm,...
        FM_L_trans,FM_L_trans_norm,FM_L_transNrot,FM_L_transNrot_norm,...
        FM_R_trans,FM_R_trans_norm,FM_R_transNrot,FM_R_transNrot_norm] =...
        calc_QSforcesNmoments_steadyNmodNbodyKin_LnR(mod_value_now,const,wing_model,body_model,freq,...
        stroke_coeffs_steady,dev_coeffs_steady,rot_coeffs_steady,...
        stroke_coeffs_MOD_L,dev_coeffs_MOD_L,rot_coeffs_MOD_L,...
        stroke_coeffs_MOD_R,dev_coeffs_MOD_R,rot_coeffs_MOD_R,plot_on);
    
    Fx_trans_norm_RollAccel_allNOfreq(:,i) = FM_trans_norm(1,:);
    Fy_trans_norm_RollAccel_allNOfreq(:,i) = FM_trans_norm(2,:);
    Fz_trans_norm_RollAccel_allNOfreq(:,i) = FM_trans_norm(3,:);
    
    Mx_trans_norm_RollAccel_allNOfreq(:,i) = FM_trans_norm(4,:);
    My_trans_norm_RollAccel_allNOfreq(:,i) = FM_trans_norm(5,:);
    Mz_trans_norm_RollAccel_allNOfreq(:,i) = FM_trans_norm(6,:);
    
    Fx_transNrot_norm_RollAccel_allNOfreq(:,i) = FM_transNrot_norm(1,:);
    Fy_transNrot_norm_RollAccel_allNOfreq(:,i) = FM_transNrot_norm(2,:);
    Fz_transNrot_norm_RollAccel_allNOfreq(:,i) = FM_transNrot_norm(3,:);
    
    Mx_transNrot_norm_RollAccel_allNOfreq(:,i) = FM_transNrot_norm(4,:);
    My_transNrot_norm_RollAccel_allNOfreq(:,i) = FM_transNrot_norm(5,:);
    Mz_transNrot_norm_RollAccel_allNOfreq(:,i) = FM_transNrot_norm(6,:);

    %% STROKE
    % set freq FIRST (same for MODnSteady)
    freq = f_wb_steady;

    % WBmod fouriers
    stroke_coeffs_MOD_L = strokeMOD_down_RollAccel_fourier_coeffs_binmean;
    dev_coeffs_MOD_L = [];
    rot_coeffs_MOD_L = [];
    
    stroke_coeffs_MOD_R = strokeMOD_up_RollAccel_fourier_coeffs_binmean;
    dev_coeffs_MOD_R = [];
    rot_coeffs_MOD_R = [];
    
    % calc FnM SYMMETRIC WB
    [FM_trans,FM_trans_norm,FM_transNrot,FM_transNrot_norm,...
        FM_L_trans,FM_L_trans_norm,FM_L_transNrot,FM_L_transNrot_norm,...
        FM_R_trans,FM_R_trans_norm,FM_R_transNrot,FM_R_transNrot_norm] =...
        calc_QSforcesNmoments_steadyNmodNbodyKin_LnR(mod_value_now,const,wing_model,body_model,freq,...
        stroke_coeffs_steady,dev_coeffs_steady,rot_coeffs_steady,...
        stroke_coeffs_MOD_L,dev_coeffs_MOD_L,rot_coeffs_MOD_L,...
        stroke_coeffs_MOD_R,dev_coeffs_MOD_R,rot_coeffs_MOD_R,0);
    
    Fx_trans_norm_RollAccel_stroke(:,i) = FM_trans_norm(1,:);
    Fy_trans_norm_RollAccel_stroke(:,i) = FM_trans_norm(2,:);
    Fz_trans_norm_RollAccel_stroke(:,i) = FM_trans_norm(3,:);
    
    Mx_trans_norm_RollAccel_stroke(:,i) = FM_trans_norm(4,:);
    My_trans_norm_RollAccel_stroke(:,i) = FM_trans_norm(5,:);
    Mz_trans_norm_RollAccel_stroke(:,i) = FM_trans_norm(6,:);
    
    Fx_transNrot_norm_RollAccel_stroke(:,i) = FM_transNrot_norm(1,:);
    Fy_transNrot_norm_RollAccel_stroke(:,i) = FM_transNrot_norm(2,:);
    Fz_transNrot_norm_RollAccel_stroke(:,i) = FM_transNrot_norm(3,:);
    
    Mx_transNrot_norm_RollAccel_stroke(:,i) = FM_transNrot_norm(4,:);
    My_transNrot_norm_RollAccel_stroke(:,i) = FM_transNrot_norm(5,:);
    Mz_transNrot_norm_RollAccel_stroke(:,i) = FM_transNrot_norm(6,:);
    
    %% DEVIATION
    % set freq FIRST (same for MODnSteady)
    freq = f_wb_steady;

    % WBmod fouriers
    stroke_coeffs_MOD_L = [];
    dev_coeffs_MOD_L = devMOD_down_RollAccel_fourier_coeffs_binmean;
    rot_coeffs_MOD_L = [];
    
    stroke_coeffs_MOD_R = [];
    dev_coeffs_MOD_R = devMOD_up_RollAccel_fourier_coeffs_binmean;
    rot_coeffs_MOD_R = [];
    
    % calc FnM SYMMETRIC WB
    [FM_trans,FM_trans_norm,FM_transNrot,FM_transNrot_norm,...
        FM_L_trans,FM_L_trans_norm,FM_L_transNrot,FM_L_transNrot_norm,...
        FM_R_trans,FM_R_trans_norm,FM_R_transNrot,FM_R_transNrot_norm] =...
        calc_QSforcesNmoments_steadyNmodNbodyKin_LnR(mod_value_now,const,wing_model,body_model,freq,...
        stroke_coeffs_steady,dev_coeffs_steady,rot_coeffs_steady,...
        stroke_coeffs_MOD_L,dev_coeffs_MOD_L,rot_coeffs_MOD_L,...
        stroke_coeffs_MOD_R,dev_coeffs_MOD_R,rot_coeffs_MOD_R,0);
    
    Fx_trans_norm_RollAccel_dev(:,i) = FM_trans_norm(1,:);
    Fy_trans_norm_RollAccel_dev(:,i) = FM_trans_norm(2,:);
    Fz_trans_norm_RollAccel_dev(:,i) = FM_trans_norm(3,:);
    
    Mx_trans_norm_RollAccel_dev(:,i) = FM_trans_norm(4,:);
    My_trans_norm_RollAccel_dev(:,i) = FM_trans_norm(5,:);
    Mz_trans_norm_RollAccel_dev(:,i) = FM_trans_norm(6,:);
    
    Fx_transNrot_norm_RollAccel_dev(:,i) = FM_transNrot_norm(1,:);
    Fy_transNrot_norm_RollAccel_dev(:,i) = FM_transNrot_norm(2,:);
    Fz_transNrot_norm_RollAccel_dev(:,i) = FM_transNrot_norm(3,:);
    
    Mx_transNrot_norm_RollAccel_dev(:,i) = FM_transNrot_norm(4,:);
    My_transNrot_norm_RollAccel_dev(:,i) = FM_transNrot_norm(5,:);
    Mz_transNrot_norm_RollAccel_dev(:,i) = FM_transNrot_norm(6,:);
    
    %% ROTATION
    % set freq FIRST (same for MODnSteady)
    freq = f_wb_steady;

    % WBmod fouriers
    stroke_coeffs_MOD_L = [];
    dev_coeffs_MOD_L = [];
    rot_coeffs_MOD_L = pitchMOD_down_RollAccel_fourier_coeffs_binmean;
    
    stroke_coeffs_MOD_R = [];
    dev_coeffs_MOD_R = [];
    rot_coeffs_MOD_R = pitchMOD_up_RollAccel_fourier_coeffs_binmean;
    
    % calc FnM SYMMETRIC WB
    [FM_trans,FM_trans_norm,FM_transNrot,FM_transNrot_norm,...
        FM_L_trans,FM_L_trans_norm,FM_L_transNrot,FM_L_transNrot_norm,...
        FM_R_trans,FM_R_trans_norm,FM_R_transNrot,FM_R_transNrot_norm] =...
        calc_QSforcesNmoments_steadyNmodNbodyKin_LnR(mod_value_now,const,wing_model,body_model,freq,...
        stroke_coeffs_steady,dev_coeffs_steady,rot_coeffs_steady,...
        stroke_coeffs_MOD_L,dev_coeffs_MOD_L,rot_coeffs_MOD_L,...
        stroke_coeffs_MOD_R,dev_coeffs_MOD_R,rot_coeffs_MOD_R,0);
    
    Fx_trans_norm_RollAccel_rot(:,i) = FM_trans_norm(1,:);
    Fy_trans_norm_RollAccel_rot(:,i) = FM_trans_norm(2,:);
    Fz_trans_norm_RollAccel_rot(:,i) = FM_trans_norm(3,:);
    
    Mx_trans_norm_RollAccel_rot(:,i) = FM_trans_norm(4,:);
    My_trans_norm_RollAccel_rot(:,i) = FM_trans_norm(5,:);
    Mz_trans_norm_RollAccel_rot(:,i) = FM_trans_norm(6,:);
    
    Fx_transNrot_norm_RollAccel_rot(:,i) = FM_transNrot_norm(1,:);
    Fy_transNrot_norm_RollAccel_rot(:,i) = FM_transNrot_norm(2,:);
    Fz_transNrot_norm_RollAccel_rot(:,i) = FM_transNrot_norm(3,:);
    
    Mx_transNrot_norm_RollAccel_rot(:,i) = FM_transNrot_norm(4,:);
    My_transNrot_norm_RollAccel_rot(:,i) = FM_transNrot_norm(5,:);
    Mz_transNrot_norm_RollAccel_rot(:,i) = FM_transNrot_norm(6,:);
end

if plot_on == 1
    cd ..
end

mod_value_RollAccel = mod_value;

% mean values
Fx_trans_norm_RollAccel_allNOfreq_mean = mean(Fx_trans_norm_RollAccel_allNOfreq);
Fy_trans_norm_RollAccel_allNOfreq_mean = mean(Fy_trans_norm_RollAccel_allNOfreq);
Fz_trans_norm_RollAccel_allNOfreq_mean = mean(Fz_trans_norm_RollAccel_allNOfreq);

Mx_trans_norm_RollAccel_allNOfreq_mean = mean(Mx_trans_norm_RollAccel_allNOfreq);
My_trans_norm_RollAccel_allNOfreq_mean = mean(My_trans_norm_RollAccel_allNOfreq);
Mz_trans_norm_RollAccel_allNOfreq_mean = mean(Mz_trans_norm_RollAccel_allNOfreq);

Fx_transNrot_norm_RollAccel_allNOfreq_mean = mean(Fx_transNrot_norm_RollAccel_allNOfreq);
Fy_transNrot_norm_RollAccel_allNOfreq_mean = mean(Fy_transNrot_norm_RollAccel_allNOfreq);
Fz_transNrot_norm_RollAccel_allNOfreq_mean = mean(Fz_transNrot_norm_RollAccel_allNOfreq);

Mx_transNrot_norm_RollAccel_allNOfreq_mean = mean(Mx_transNrot_norm_RollAccel_allNOfreq);
My_transNrot_norm_RollAccel_allNOfreq_mean = mean(My_transNrot_norm_RollAccel_allNOfreq);
Mz_transNrot_norm_RollAccel_allNOfreq_mean = mean(Mz_transNrot_norm_RollAccel_allNOfreq);

Fx_trans_norm_RollAccel_stroke_mean = mean(Fx_trans_norm_RollAccel_stroke);
Fy_trans_norm_RollAccel_stroke_mean = mean(Fy_trans_norm_RollAccel_stroke);
Fz_trans_norm_RollAccel_stroke_mean = mean(Fz_trans_norm_RollAccel_stroke);

Mx_trans_norm_RollAccel_stroke_mean = mean(Mx_trans_norm_RollAccel_stroke);
My_trans_norm_RollAccel_stroke_mean = mean(My_trans_norm_RollAccel_stroke);
Mz_trans_norm_RollAccel_stroke_mean = mean(Mz_trans_norm_RollAccel_stroke);

Fx_transNrot_norm_RollAccel_stroke_mean = mean(Fx_transNrot_norm_RollAccel_stroke);
Fy_transNrot_norm_RollAccel_stroke_mean = mean(Fy_transNrot_norm_RollAccel_stroke);
Fz_transNrot_norm_RollAccel_stroke_mean = mean(Fz_transNrot_norm_RollAccel_stroke);

Mx_transNrot_norm_RollAccel_stroke_mean = mean(Mx_transNrot_norm_RollAccel_stroke);
My_transNrot_norm_RollAccel_stroke_mean = mean(My_transNrot_norm_RollAccel_stroke);
Mz_transNrot_norm_RollAccel_stroke_mean = mean(Mz_transNrot_norm_RollAccel_stroke);

Fx_trans_norm_RollAccel_dev_mean = mean(Fx_trans_norm_RollAccel_dev);
Fy_trans_norm_RollAccel_dev_mean = mean(Fy_trans_norm_RollAccel_dev);
Fz_trans_norm_RollAccel_dev_mean = mean(Fz_trans_norm_RollAccel_dev);

Mx_trans_norm_RollAccel_dev_mean = mean(Mx_trans_norm_RollAccel_dev);
My_trans_norm_RollAccel_dev_mean = mean(My_trans_norm_RollAccel_dev);
Mz_trans_norm_RollAccel_dev_mean = mean(Mz_trans_norm_RollAccel_dev);

Fx_transNrot_norm_RollAccel_dev_mean = mean(Fx_transNrot_norm_RollAccel_dev);
Fy_transNrot_norm_RollAccel_dev_mean = mean(Fy_transNrot_norm_RollAccel_dev);
Fz_transNrot_norm_RollAccel_dev_mean = mean(Fz_transNrot_norm_RollAccel_dev);

Mx_transNrot_norm_RollAccel_dev_mean = mean(Mx_transNrot_norm_RollAccel_dev);
My_transNrot_norm_RollAccel_dev_mean = mean(My_transNrot_norm_RollAccel_dev);
Mz_transNrot_norm_RollAccel_dev_mean = mean(Mz_transNrot_norm_RollAccel_dev);

Fx_trans_norm_RollAccel_rot_mean = mean(Fx_trans_norm_RollAccel_rot);
Fy_trans_norm_RollAccel_rot_mean = mean(Fy_trans_norm_RollAccel_rot);
Fz_trans_norm_RollAccel_rot_mean = mean(Fz_trans_norm_RollAccel_rot);

Mx_trans_norm_RollAccel_rot_mean = mean(Mx_trans_norm_RollAccel_rot);
My_trans_norm_RollAccel_rot_mean = mean(My_trans_norm_RollAccel_rot);
Mz_trans_norm_RollAccel_rot_mean = mean(Mz_trans_norm_RollAccel_rot);

Fx_transNrot_norm_RollAccel_rot_mean = mean(Fx_transNrot_norm_RollAccel_rot);
Fy_transNrot_norm_RollAccel_rot_mean = mean(Fy_transNrot_norm_RollAccel_rot);
Fz_transNrot_norm_RollAccel_rot_mean = mean(Fz_transNrot_norm_RollAccel_rot);

Mx_transNrot_norm_RollAccel_rot_mean = mean(Mx_transNrot_norm_RollAccel_rot);
My_transNrot_norm_RollAccel_rot_mean = mean(My_transNrot_norm_RollAccel_rot);
Mz_transNrot_norm_RollAccel_rot_mean = mean(Mz_transNrot_norm_RollAccel_rot);

Ftot_trans_norm_RollAccel_allNOfreq_mean = sqrt(Fx_trans_norm_RollAccel_allNOfreq_mean.^2 + Fy_trans_norm_RollAccel_allNOfreq_mean.^2 + Fz_trans_norm_RollAccel_allNOfreq_mean.^2);
Mtot_trans_norm_RollAccel_allNOfreq_mean = sqrt(Mx_trans_norm_RollAccel_allNOfreq_mean.^2 + My_trans_norm_RollAccel_allNOfreq_mean.^2 + Mz_trans_norm_RollAccel_allNOfreq_mean.^2);
Ftot_transNrot_norm_RollAccel_allNOfreq_mean = sqrt(Fx_transNrot_norm_RollAccel_allNOfreq_mean.^2 + Fy_transNrot_norm_RollAccel_allNOfreq_mean.^2 + Fz_transNrot_norm_RollAccel_allNOfreq_mean.^2);
Mtot_transNrot_norm_RollAccel_allNOfreq_mean = sqrt(Mx_transNrot_norm_RollAccel_allNOfreq_mean.^2 + My_transNrot_norm_RollAccel_allNOfreq_mean.^2 + Mz_transNrot_norm_RollAccel_allNOfreq_mean.^2);

Ftot_trans_norm_RollAccel_stroke_mean = sqrt(Fx_trans_norm_RollAccel_stroke_mean.^2 + Fy_trans_norm_RollAccel_stroke_mean.^2 + Fz_trans_norm_RollAccel_stroke_mean.^2);
Mtot_trans_norm_RollAccel_stroke_mean = sqrt(Mx_trans_norm_RollAccel_stroke_mean.^2 + My_trans_norm_RollAccel_stroke_mean.^2 + Mz_trans_norm_RollAccel_stroke_mean.^2);
Ftot_transNrot_norm_RollAccel_stroke_mean = sqrt(Fx_transNrot_norm_RollAccel_stroke_mean.^2 + Fy_transNrot_norm_RollAccel_stroke_mean.^2 + Fz_transNrot_norm_RollAccel_stroke_mean.^2);
Mtot_transNrot_norm_RollAccel_stroke_mean = sqrt(Mx_transNrot_norm_RollAccel_stroke_mean.^2 + My_transNrot_norm_RollAccel_stroke_mean.^2 + Mz_transNrot_norm_RollAccel_stroke_mean.^2);

Ftot_trans_norm_RollAccel_dev_mean = sqrt(Fx_trans_norm_RollAccel_dev_mean.^2 + Fy_trans_norm_RollAccel_dev_mean.^2 + Fz_trans_norm_RollAccel_dev_mean.^2);
Mtot_trans_norm_RollAccel_dev_mean = sqrt(Mx_trans_norm_RollAccel_dev_mean.^2 + My_trans_norm_RollAccel_dev_mean.^2 + Mz_trans_norm_RollAccel_dev_mean.^2);
Ftot_transNrot_norm_RollAccel_dev_mean = sqrt(Fx_transNrot_norm_RollAccel_dev_mean.^2 + Fy_transNrot_norm_RollAccel_dev_mean.^2 + Fz_transNrot_norm_RollAccel_dev_mean.^2);
Mtot_transNrot_norm_RollAccel_dev_mean = sqrt(Mx_transNrot_norm_RollAccel_dev_mean.^2 + My_transNrot_norm_RollAccel_dev_mean.^2 + Mz_transNrot_norm_RollAccel_dev_mean.^2);

Ftot_trans_norm_RollAccel_rot_mean = sqrt(Fx_trans_norm_RollAccel_rot_mean.^2 + Fy_trans_norm_RollAccel_rot_mean.^2 + Fz_trans_norm_RollAccel_rot_mean.^2);
Mtot_trans_norm_RollAccel_rot_mean = sqrt(Mx_trans_norm_RollAccel_rot_mean.^2 + My_trans_norm_RollAccel_rot_mean.^2 + Mz_trans_norm_RollAccel_rot_mean.^2);
Ftot_transNrot_norm_RollAccel_rot_mean = sqrt(Fx_transNrot_norm_RollAccel_rot_mean.^2 + Fy_transNrot_norm_RollAccel_rot_mean.^2 + Fz_transNrot_norm_RollAccel_rot_mean.^2);
Mtot_transNrot_norm_RollAccel_rot_mean = sqrt(Mx_transNrot_norm_RollAccel_rot_mean.^2 + My_transNrot_norm_RollAccel_rot_mean.^2 + Mz_transNrot_norm_RollAccel_rot_mean.^2);

save('FnMqs_RollAccel.mat','mod_value_RollAccel',...
    'Fx_trans_norm_RollAccel_allNOfreq','Fy_trans_norm_RollAccel_allNOfreq','Fz_trans_norm_RollAccel_allNOfreq',...
    'Mx_trans_norm_RollAccel_allNOfreq','My_trans_norm_RollAccel_allNOfreq','Mz_trans_norm_RollAccel_allNOfreq',...
    'Fx_transNrot_norm_RollAccel_allNOfreq','Fy_transNrot_norm_RollAccel_allNOfreq','Fz_transNrot_norm_RollAccel_allNOfreq',...
    'Mx_transNrot_norm_RollAccel_allNOfreq','My_transNrot_norm_RollAccel_allNOfreq','Mz_transNrot_norm_RollAccel_allNOfreq',...
    'Fx_trans_norm_RollAccel_stroke','Fy_trans_norm_RollAccel_stroke','Fz_trans_norm_RollAccel_stroke',...
    'Mx_trans_norm_RollAccel_stroke','My_trans_norm_RollAccel_stroke','Mz_trans_norm_RollAccel_stroke',...
    'Fx_transNrot_norm_RollAccel_stroke','Fy_transNrot_norm_RollAccel_stroke','Fz_transNrot_norm_RollAccel_stroke',...
    'Mx_transNrot_norm_RollAccel_stroke','My_transNrot_norm_RollAccel_stroke','Mz_transNrot_norm_RollAccel_stroke',...
    'Fx_trans_norm_RollAccel_dev','Fy_trans_norm_RollAccel_dev','Fz_trans_norm_RollAccel_dev',...
    'Mx_trans_norm_RollAccel_dev','My_trans_norm_RollAccel_dev','Mz_trans_norm_RollAccel_dev',...
    'Fx_transNrot_norm_RollAccel_dev','Fy_transNrot_norm_RollAccel_dev','Fz_transNrot_norm_RollAccel_dev',...
    'Mx_transNrot_norm_RollAccel_dev','My_transNrot_norm_RollAccel_dev','Mz_transNrot_norm_RollAccel_dev',...
    'Fx_trans_norm_RollAccel_rot','Fy_trans_norm_RollAccel_rot','Fz_trans_norm_RollAccel_rot',...
    'Mx_trans_norm_RollAccel_rot','My_trans_norm_RollAccel_rot','Mz_trans_norm_RollAccel_rot',...
    'Fx_transNrot_norm_RollAccel_rot','Fy_transNrot_norm_RollAccel_rot','Fz_transNrot_norm_RollAccel_rot',...
    'Mx_transNrot_norm_RollAccel_rot','My_transNrot_norm_RollAccel_rot','Mz_transNrot_norm_RollAccel_rot',...
    'Fx_trans_norm_RollAccel_allNOfreq_mean','Fy_trans_norm_RollAccel_allNOfreq_mean','Fz_trans_norm_RollAccel_allNOfreq_mean',...
    'Mx_trans_norm_RollAccel_allNOfreq_mean','My_trans_norm_RollAccel_allNOfreq_mean','Mz_trans_norm_RollAccel_allNOfreq_mean',...
    'Fx_transNrot_norm_RollAccel_allNOfreq_mean','Fy_transNrot_norm_RollAccel_allNOfreq_mean','Fz_transNrot_norm_RollAccel_allNOfreq_mean',...
    'Mx_transNrot_norm_RollAccel_allNOfreq_mean','My_transNrot_norm_RollAccel_allNOfreq_mean','Mz_transNrot_norm_RollAccel_allNOfreq_mean',...
    'Fx_trans_norm_RollAccel_stroke_mean','Fy_trans_norm_RollAccel_stroke_mean','Fz_trans_norm_RollAccel_stroke_mean',...
    'Mx_trans_norm_RollAccel_stroke_mean','My_trans_norm_RollAccel_stroke_mean','Mz_trans_norm_RollAccel_stroke_mean',...
    'Fx_transNrot_norm_RollAccel_stroke_mean','Fy_transNrot_norm_RollAccel_stroke_mean','Fz_transNrot_norm_RollAccel_stroke_mean',...
    'Mx_transNrot_norm_RollAccel_stroke_mean','My_transNrot_norm_RollAccel_stroke_mean','Mz_transNrot_norm_RollAccel_stroke_mean',...
    'Fx_trans_norm_RollAccel_dev_mean','Fy_trans_norm_RollAccel_dev_mean','Fz_trans_norm_RollAccel_dev_mean',...
    'Mx_trans_norm_RollAccel_dev_mean','My_trans_norm_RollAccel_dev_mean','Mz_trans_norm_RollAccel_dev_mean',...
    'Fx_transNrot_norm_RollAccel_dev_mean','Fy_transNrot_norm_RollAccel_dev_mean','Fz_transNrot_norm_RollAccel_dev_mean',...
    'Mx_transNrot_norm_RollAccel_dev_mean','My_transNrot_norm_RollAccel_dev_mean','Mz_transNrot_norm_RollAccel_dev_mean',...
    'Fx_trans_norm_RollAccel_rot_mean','Fy_trans_norm_RollAccel_rot_mean','Fz_trans_norm_RollAccel_rot_mean',...
    'Mx_trans_norm_RollAccel_rot_mean','My_trans_norm_RollAccel_rot_mean','Mz_trans_norm_RollAccel_rot_mean',...
    'Fx_transNrot_norm_RollAccel_rot_mean','Fy_transNrot_norm_RollAccel_rot_mean','Fz_transNrot_norm_RollAccel_rot_mean',...
    'Mx_transNrot_norm_RollAccel_rot_mean','My_transNrot_norm_RollAccel_rot_mean','Mz_transNrot_norm_RollAccel_rot_mean',...
    'Ftot_trans_norm_RollAccel_allNOfreq_mean',...
    'Mtot_trans_norm_RollAccel_allNOfreq_mean',...
    'Ftot_transNrot_norm_RollAccel_allNOfreq_mean',...
    'Mtot_transNrot_norm_RollAccel_allNOfreq_mean',...
    'Ftot_trans_norm_RollAccel_stroke_mean',...
    'Mtot_trans_norm_RollAccel_stroke_mean',...
    'Ftot_transNrot_norm_RollAccel_stroke_mean',...
    'Mtot_transNrot_norm_RollAccel_stroke_mean',...
    'Ftot_trans_norm_RollAccel_dev_mean',...
    'Mtot_trans_norm_RollAccel_dev_mean',...
    'Ftot_transNrot_norm_RollAccel_dev_mean',...
    'Mtot_transNrot_norm_RollAccel_dev_mean',...
    'Ftot_trans_norm_RollAccel_rot_mean',...
    'Mtot_trans_norm_RollAccel_rot_mean',...
    'Ftot_transNrot_norm_RollAccel_rot_mean',...
    'Mtot_transNrot_norm_RollAccel_rot_mean')

