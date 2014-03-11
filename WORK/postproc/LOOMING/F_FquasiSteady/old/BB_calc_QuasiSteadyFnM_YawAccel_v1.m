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
        calc_QSforcesNmoments_steadyNmodNbodyKin_LnR(const,wing_model,body_model,freq,...
        stroke_coeffs_steady,dev_coeffs_steady,rot_coeffs_steady,...
        stroke_coeffs_MOD_L,dev_coeffs_MOD_L,rot_coeffs_MOD_L,...
        stroke_coeffs_MOD_R,dev_coeffs_MOD_R,rot_coeffs_MOD_R);
    
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
        calc_QSforcesNmoments_steadyNmodNbodyKin_LnR(const,wing_model,body_model,freq,...
        stroke_coeffs_steady,dev_coeffs_steady,rot_coeffs_steady,...
        stroke_coeffs_MOD_L,dev_coeffs_MOD_L,rot_coeffs_MOD_L,...
        stroke_coeffs_MOD_R,dev_coeffs_MOD_R,rot_coeffs_MOD_R);
    
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
        calc_QSforcesNmoments_steadyNmodNbodyKin_LnR(const,wing_model,body_model,freq,...
        stroke_coeffs_steady,dev_coeffs_steady,rot_coeffs_steady,...
        stroke_coeffs_MOD_L,dev_coeffs_MOD_L,rot_coeffs_MOD_L,...
        stroke_coeffs_MOD_R,dev_coeffs_MOD_R,rot_coeffs_MOD_R);
    
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
        calc_QSforcesNmoments_steadyNmodNbodyKin_LnR(const,wing_model,body_model,freq,...
        stroke_coeffs_steady,dev_coeffs_steady,rot_coeffs_steady,...
        stroke_coeffs_MOD_L,dev_coeffs_MOD_L,rot_coeffs_MOD_L,...
        stroke_coeffs_MOD_R,dev_coeffs_MOD_R,rot_coeffs_MOD_R);
    
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

mod_value_RollAccel = mod_value;

Ftot_trans_norm_RollAccel_allNOfreq = sqrt(Fx_trans_norm_RollAccel_allNOfreq.^2 + Fy_trans_norm_RollAccel_allNOfreq.^2 + Fz_trans_norm_RollAccel_allNOfreq.^2);
Mtot_trans_norm_RollAccel_allNOfreq = sqrt(Mx_trans_norm_RollAccel_allNOfreq.^2 + My_trans_norm_RollAccel_allNOfreq.^2 + Mz_trans_norm_RollAccel_allNOfreq.^2);
Ftot_transNrot_norm_RollAccel_allNOfreq = sqrt(Fx_transNrot_norm_RollAccel_allNOfreq.^2 + Fy_transNrot_norm_RollAccel_allNOfreq.^2 + Fz_transNrot_norm_RollAccel_allNOfreq.^2);
Mtot_transNrot_norm_RollAccel_allNOfreq = sqrt(Mx_transNrot_norm_RollAccel_allNOfreq.^2 + My_transNrot_norm_RollAccel_allNOfreq.^2 + Mz_transNrot_norm_RollAccel_allNOfreq.^2);

Ftot_trans_norm_RollAccel_stroke = sqrt(Fx_trans_norm_RollAccel_stroke.^2 + Fy_trans_norm_RollAccel_stroke.^2 + Fz_trans_norm_RollAccel_stroke.^2);
Mtot_trans_norm_RollAccel_stroke = sqrt(Mx_trans_norm_RollAccel_stroke.^2 + My_trans_norm_RollAccel_stroke.^2 + Mz_trans_norm_RollAccel_stroke.^2);
Ftot_transNrot_norm_RollAccel_stroke = sqrt(Fx_transNrot_norm_RollAccel_stroke.^2 + Fy_transNrot_norm_RollAccel_stroke.^2 + Fz_transNrot_norm_RollAccel_stroke.^2);
Mtot_transNrot_norm_RollAccel_stroke = sqrt(Mx_transNrot_norm_RollAccel_stroke.^2 + My_transNrot_norm_RollAccel_stroke.^2 + Mz_transNrot_norm_RollAccel_stroke.^2);

Ftot_trans_norm_RollAccel_dev = sqrt(Fx_trans_norm_RollAccel_dev.^2 + Fy_trans_norm_RollAccel_dev.^2 + Fz_trans_norm_RollAccel_dev.^2);
Mtot_trans_norm_RollAccel_dev = sqrt(Mx_trans_norm_RollAccel_dev.^2 + My_trans_norm_RollAccel_dev.^2 + Mz_trans_norm_RollAccel_dev.^2);
Ftot_transNrot_norm_RollAccel_dev = sqrt(Fx_transNrot_norm_RollAccel_dev.^2 + Fy_transNrot_norm_RollAccel_dev.^2 + Fz_transNrot_norm_RollAccel_dev.^2);
Mtot_transNrot_norm_RollAccel_dev = sqrt(Mx_transNrot_norm_RollAccel_dev.^2 + My_transNrot_norm_RollAccel_dev.^2 + Mz_transNrot_norm_RollAccel_dev.^2);

Ftot_trans_norm_RollAccel_rot = sqrt(Fx_trans_norm_RollAccel_rot.^2 + Fy_trans_norm_RollAccel_rot.^2 + Fz_trans_norm_RollAccel_rot.^2);
Mtot_trans_norm_RollAccel_rot = sqrt(Mx_trans_norm_RollAccel_rot.^2 + My_trans_norm_RollAccel_rot.^2 + Mz_trans_norm_RollAccel_rot.^2);
Ftot_transNrot_norm_RollAccel_rot = sqrt(Fx_transNrot_norm_RollAccel_rot.^2 + Fy_transNrot_norm_RollAccel_rot.^2 + Fz_transNrot_norm_RollAccel_rot.^2);
Mtot_transNrot_norm_RollAccel_rot = sqrt(Mx_transNrot_norm_RollAccel_rot.^2 + My_transNrot_norm_RollAccel_rot.^2 + Mz_transNrot_norm_RollAccel_rot.^2);

save('FnMqs_RollAccel.mat','mod_value_RollAccel',...
    'Ftot_trans_norm_RollAccel_allNOfreq','Mtot_trans_norm_RollAccel_allNOfreq',...
    'Ftot_transNrot_norm_RollAccel_allNOfreq','Mtot_transNrot_norm_RollAccel_allNOfreq',...
    'Fx_trans_norm_RollAccel_allNOfreq','Fy_trans_norm_RollAccel_allNOfreq','Fz_trans_norm_RollAccel_allNOfreq',...
    'Mx_trans_norm_RollAccel_allNOfreq','My_trans_norm_RollAccel_allNOfreq','Mz_trans_norm_RollAccel_allNOfreq',...
    'Fx_transNrot_norm_RollAccel_allNOfreq','Fy_transNrot_norm_RollAccel_allNOfreq','Fz_transNrot_norm_RollAccel_allNOfreq',...
    'Mx_transNrot_norm_RollAccel_allNOfreq','My_transNrot_norm_RollAccel_allNOfreq','Mz_transNrot_norm_RollAccel_allNOfreq',...
    'Ftot_trans_norm_RollAccel_stroke','Mtot_trans_norm_RollAccel_stroke',...
    'Ftot_transNrot_norm_RollAccel_stroke','Mtot_transNrot_norm_RollAccel_stroke',...
    'Fx_trans_norm_RollAccel_stroke','Fy_trans_norm_RollAccel_stroke','Fz_trans_norm_RollAccel_stroke',...
    'Mx_trans_norm_RollAccel_stroke','My_trans_norm_RollAccel_stroke','Mz_trans_norm_RollAccel_stroke',...
    'Fx_transNrot_norm_RollAccel_stroke','Fy_transNrot_norm_RollAccel_stroke','Fz_transNrot_norm_RollAccel_stroke',...
    'Mx_transNrot_norm_RollAccel_stroke','My_transNrot_norm_RollAccel_stroke','Mz_transNrot_norm_RollAccel_stroke',...
    'Ftot_trans_norm_RollAccel_dev','Mtot_trans_norm_RollAccel_dev',...
    'Ftot_transNrot_norm_RollAccel_dev','Mtot_transNrot_norm_RollAccel_dev',...
    'Fx_trans_norm_RollAccel_dev','Fy_trans_norm_RollAccel_dev','Fz_trans_norm_RollAccel_dev',...
    'Mx_trans_norm_RollAccel_dev','My_trans_norm_RollAccel_dev','Mz_trans_norm_RollAccel_dev',...
    'Fx_transNrot_norm_RollAccel_dev','Fy_transNrot_norm_RollAccel_dev','Fz_transNrot_norm_RollAccel_dev',...
    'Mx_transNrot_norm_RollAccel_dev','My_transNrot_norm_RollAccel_dev','Mz_transNrot_norm_RollAccel_dev',...
    'Ftot_trans_norm_RollAccel_rot','Mtot_trans_norm_RollAccel_rot',...
    'Ftot_transNrot_norm_RollAccel_rot','Mtot_transNrot_norm_RollAccel_rot',...
    'Fx_trans_norm_RollAccel_rot','Fy_trans_norm_RollAccel_rot','Fz_trans_norm_RollAccel_rot',...
    'Mx_trans_norm_RollAccel_rot','My_trans_norm_RollAccel_rot','Mz_trans_norm_RollAccel_rot',...
    'Fx_transNrot_norm_RollAccel_rot','Fy_transNrot_norm_RollAccel_rot','Fz_transNrot_norm_RollAccel_rot',...
    'Mx_transNrot_norm_RollAccel_rot','My_transNrot_norm_RollAccel_rot','Mz_transNrot_norm_RollAccel_rot')


