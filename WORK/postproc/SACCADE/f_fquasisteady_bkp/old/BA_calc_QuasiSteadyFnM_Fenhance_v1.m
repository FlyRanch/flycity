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
load('WBmod_accelbased_Fenhanced_719WBs.mat')
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
    
    %% ALL MODS
    % set freq FIRST (same for MODnSteady)
    freq = f_wb_steady + mod_value_now*freqMOD_wb_Fenhance_meanCIstd(1);

    % WBmod fouriers
    stroke_coeffs_MOD = strokeMOD_Fenhance_fourier_coeffs_binmean;
    dev_coeffs_MOD = devMOD_Fenhance_fourier_coeffs_binmean;
    rot_coeffs_MOD = pitchMOD_Fenhance_fourier_coeffs_binmean;
    
    % calc FnM SYMMETRIC WB
    [FM_trans,FM_trans_norm,FM_transNrot,FM_transNrot_norm,...
        FM_L_trans,FM_L_trans_norm,FM_L_transNrot,FM_L_transNrot_norm,...
        FM_R_trans,FM_R_trans_norm,FM_R_transNrot,FM_R_transNrot_norm] =...
        calc_QSforcesNmoments_steadyNmodNbodyKin_symm(const,wing_model,body_model,freq,...
        stroke_coeffs_steady,dev_coeffs_steady,rot_coeffs_steady,...
        stroke_coeffs_MOD,dev_coeffs_MOD,rot_coeffs_MOD);
    
    Fx_trans_norm_Fenhance_all(:,i) = FM_trans_norm(1,:);
    Fy_trans_norm_Fenhance_all(:,i) = FM_trans_norm(2,:);
    Fz_trans_norm_Fenhance_all(:,i) = FM_trans_norm(3,:);
    
    Mx_trans_norm_Fenhance_all(:,i) = FM_trans_norm(4,:);
    My_trans_norm_Fenhance_all(:,i) = FM_trans_norm(5,:);
    Mz_trans_norm_Fenhance_all(:,i) = FM_trans_norm(6,:);
    
    Fx_transNrot_norm_Fenhance_all(:,i) = FM_transNrot_norm(1,:);
    Fy_transNrot_norm_Fenhance_all(:,i) = FM_transNrot_norm(2,:);
    Fz_transNrot_norm_Fenhance_all(:,i) = FM_transNrot_norm(3,:);
    
    Mx_transNrot_norm_Fenhance_all(:,i) = FM_transNrot_norm(4,:);
    My_transNrot_norm_Fenhance_all(:,i) = FM_transNrot_norm(5,:);
    Mz_transNrot_norm_Fenhance_all(:,i) = FM_transNrot_norm(6,:);

    %% FREQ
    % set freq FIRST (same for MODnSteady)
    freq = f_wb_steady + mod_value_now*freqMOD_wb_Fenhance_meanCIstd(1);

    % WBmod fouriers
    stroke_coeffs_MOD = [];
    dev_coeffs_MOD = [];
    rot_coeffs_MOD = [];
    
    % calc FnM SYMMETRIC WB
    [FM_trans,FM_trans_norm,FM_transNrot,FM_transNrot_norm,...
        FM_L_trans,FM_L_trans_norm,FM_L_transNrot,FM_L_transNrot_norm,...
        FM_R_trans,FM_R_trans_norm,FM_R_transNrot,FM_R_transNrot_norm] =...
        calc_QSforcesNmoments_steadyNmodNbodyKin_symm(const,wing_model,body_model,freq,...
        stroke_coeffs_steady,dev_coeffs_steady,rot_coeffs_steady,...
        stroke_coeffs_MOD,dev_coeffs_MOD,rot_coeffs_MOD);
    
    Fx_trans_norm_Fenhance_freq(:,i) = FM_trans_norm(1,:);
    Fy_trans_norm_Fenhance_freq(:,i) = FM_trans_norm(2,:);
    Fz_trans_norm_Fenhance_freq(:,i) = FM_trans_norm(3,:);
    
    Mx_trans_norm_Fenhance_freq(:,i) = FM_trans_norm(4,:);
    My_trans_norm_Fenhance_freq(:,i) = FM_trans_norm(5,:);
    Mz_trans_norm_Fenhance_freq(:,i) = FM_trans_norm(6,:);
    
    Fx_transNrot_norm_Fenhance_freq(:,i) = FM_transNrot_norm(1,:);
    Fy_transNrot_norm_Fenhance_freq(:,i) = FM_transNrot_norm(2,:);
    Fz_transNrot_norm_Fenhance_freq(:,i) = FM_transNrot_norm(3,:);
    
    Mx_transNrot_norm_Fenhance_freq(:,i) = FM_transNrot_norm(4,:);
    My_transNrot_norm_Fenhance_freq(:,i) = FM_transNrot_norm(5,:);
    Mz_transNrot_norm_Fenhance_freq(:,i) = FM_transNrot_norm(6,:);
    
    %% STROKE
    % set freq FIRST (same for MODnSteady)
    freq = f_wb_steady;

    % WBmod fouriers
    stroke_coeffs_MOD = strokeMOD_Fenhance_fourier_coeffs_binmean;
    dev_coeffs_MOD = [];
    rot_coeffs_MOD = [];
    
    % calc FnM SYMMETRIC WB
    [FM_trans,FM_trans_norm,FM_transNrot,FM_transNrot_norm,...
        FM_L_trans,FM_L_trans_norm,FM_L_transNrot,FM_L_transNrot_norm,...
        FM_R_trans,FM_R_trans_norm,FM_R_transNrot,FM_R_transNrot_norm] =...
        calc_QSforcesNmoments_steadyNmodNbodyKin_symm(const,wing_model,body_model,freq,...
        stroke_coeffs_steady,dev_coeffs_steady,rot_coeffs_steady,...
        stroke_coeffs_MOD,dev_coeffs_MOD,rot_coeffs_MOD);
    
    Fx_trans_norm_Fenhance_stroke(:,i) = FM_trans_norm(1,:);
    Fy_trans_norm_Fenhance_stroke(:,i) = FM_trans_norm(2,:);
    Fz_trans_norm_Fenhance_stroke(:,i) = FM_trans_norm(3,:);
    
    Mx_trans_norm_Fenhance_stroke(:,i) = FM_trans_norm(4,:);
    My_trans_norm_Fenhance_stroke(:,i) = FM_trans_norm(5,:);
    Mz_trans_norm_Fenhance_stroke(:,i) = FM_trans_norm(6,:);
    
    Fx_transNrot_norm_Fenhance_stroke(:,i) = FM_transNrot_norm(1,:);
    Fy_transNrot_norm_Fenhance_stroke(:,i) = FM_transNrot_norm(2,:);
    Fz_transNrot_norm_Fenhance_stroke(:,i) = FM_transNrot_norm(3,:);
    
    Mx_transNrot_norm_Fenhance_stroke(:,i) = FM_transNrot_norm(4,:);
    My_transNrot_norm_Fenhance_stroke(:,i) = FM_transNrot_norm(5,:);
    Mz_transNrot_norm_Fenhance_stroke(:,i) = FM_transNrot_norm(6,:);
    
    %% DEVIATION
    % set freq FIRST (same for MODnSteady)
    freq = f_wb_steady;

    % WBmod fouriers
    stroke_coeffs_MOD = [];
    dev_coeffs_MOD = devMOD_Fenhance_fourier_coeffs_binmean;
    rot_coeffs_MOD = [];
    
    % calc FnM SYMMETRIC WB
    [FM_trans,FM_trans_norm,FM_transNrot,FM_transNrot_norm,...
        FM_L_trans,FM_L_trans_norm,FM_L_transNrot,FM_L_transNrot_norm,...
        FM_R_trans,FM_R_trans_norm,FM_R_transNrot,FM_R_transNrot_norm] =...
        calc_QSforcesNmoments_steadyNmodNbodyKin_symm(const,wing_model,body_model,freq,...
        stroke_coeffs_steady,dev_coeffs_steady,rot_coeffs_steady,...
        stroke_coeffs_MOD,dev_coeffs_MOD,rot_coeffs_MOD);
    
    Fx_trans_norm_Fenhance_dev(:,i) = FM_trans_norm(1,:);
    Fy_trans_norm_Fenhance_dev(:,i) = FM_trans_norm(2,:);
    Fz_trans_norm_Fenhance_dev(:,i) = FM_trans_norm(3,:);
    
    Mx_trans_norm_Fenhance_dev(:,i) = FM_trans_norm(4,:);
    My_trans_norm_Fenhance_dev(:,i) = FM_trans_norm(5,:);
    Mz_trans_norm_Fenhance_dev(:,i) = FM_trans_norm(6,:);
    
    Fx_transNrot_norm_Fenhance_dev(:,i) = FM_transNrot_norm(1,:);
    Fy_transNrot_norm_Fenhance_dev(:,i) = FM_transNrot_norm(2,:);
    Fz_transNrot_norm_Fenhance_dev(:,i) = FM_transNrot_norm(3,:);
    
    Mx_transNrot_norm_Fenhance_dev(:,i) = FM_transNrot_norm(4,:);
    My_transNrot_norm_Fenhance_dev(:,i) = FM_transNrot_norm(5,:);
    Mz_transNrot_norm_Fenhance_dev(:,i) = FM_transNrot_norm(6,:);
    
    %% ROTATION
    % set freq FIRST (same for MODnSteady)
    freq = f_wb_steady;

    % WBmod fouriers
    stroke_coeffs_MOD = [];
    dev_coeffs_MOD = [];
    rot_coeffs_MOD = pitchMOD_Fenhance_fourier_coeffs_binmean;
    
    % calc FnM SYMMETRIC WB
    [FM_trans,FM_trans_norm,FM_transNrot,FM_transNrot_norm,...
        FM_L_trans,FM_L_trans_norm,FM_L_transNrot,FM_L_transNrot_norm,...
        FM_R_trans,FM_R_trans_norm,FM_R_transNrot,FM_R_transNrot_norm] =...
        calc_QSforcesNmoments_steadyNmodNbodyKin_symm(const,wing_model,body_model,freq,...
        stroke_coeffs_steady,dev_coeffs_steady,rot_coeffs_steady,...
        stroke_coeffs_MOD,dev_coeffs_MOD,rot_coeffs_MOD);
    
    Fx_trans_norm_Fenhance_rot(:,i) = FM_trans_norm(1,:);
    Fy_trans_norm_Fenhance_rot(:,i) = FM_trans_norm(2,:);
    Fz_trans_norm_Fenhance_rot(:,i) = FM_trans_norm(3,:);
    
    Mx_trans_norm_Fenhance_rot(:,i) = FM_trans_norm(4,:);
    My_trans_norm_Fenhance_rot(:,i) = FM_trans_norm(5,:);
    Mz_trans_norm_Fenhance_rot(:,i) = FM_trans_norm(6,:);
    
    Fx_transNrot_norm_Fenhance_rot(:,i) = FM_transNrot_norm(1,:);
    Fy_transNrot_norm_Fenhance_rot(:,i) = FM_transNrot_norm(2,:);
    Fz_transNrot_norm_Fenhance_rot(:,i) = FM_transNrot_norm(3,:);
    
    Mx_transNrot_norm_Fenhance_rot(:,i) = FM_transNrot_norm(4,:);
    My_transNrot_norm_Fenhance_rot(:,i) = FM_transNrot_norm(5,:);
    Mz_transNrot_norm_Fenhance_rot(:,i) = FM_transNrot_norm(6,:);
end

mod_value_Fenhance = mod_value;

Ftot_trans_norm_Fenhance_all = sqrt(Fx_trans_norm_Fenhance_all.^2 + Fy_trans_norm_Fenhance_all.^2 + Fz_trans_norm_Fenhance_all.^2);
Mtot_trans_norm_Fenhance_all = sqrt(Mx_trans_norm_Fenhance_all.^2 + My_trans_norm_Fenhance_all.^2 + Mz_trans_norm_Fenhance_all.^2);
Ftot_transNrot_norm_Fenhance_all = sqrt(Fx_transNrot_norm_Fenhance_all.^2 + Fy_transNrot_norm_Fenhance_all.^2 + Fz_transNrot_norm_Fenhance_all.^2);
Mtot_transNrot_norm_Fenhance_all = sqrt(Mx_transNrot_norm_Fenhance_all.^2 + My_transNrot_norm_Fenhance_all.^2 + Mz_transNrot_norm_Fenhance_all.^2);

Ftot_trans_norm_Fenhance_freq = sqrt(Fx_trans_norm_Fenhance_freq.^2 + Fy_trans_norm_Fenhance_freq.^2 + Fz_trans_norm_Fenhance_freq.^2);
Mtot_trans_norm_Fenhance_freq = sqrt(Mx_trans_norm_Fenhance_freq.^2 + My_trans_norm_Fenhance_freq.^2 + Mz_trans_norm_Fenhance_freq.^2);
Ftot_transNrot_norm_Fenhance_freq = sqrt(Fx_transNrot_norm_Fenhance_freq.^2 + Fy_transNrot_norm_Fenhance_freq.^2 + Fz_transNrot_norm_Fenhance_freq.^2);
Mtot_transNrot_norm_Fenhance_freq = sqrt(Mx_transNrot_norm_Fenhance_freq.^2 + My_transNrot_norm_Fenhance_freq.^2 + Mz_transNrot_norm_Fenhance_freq.^2);

Ftot_trans_norm_Fenhance_stroke = sqrt(Fx_trans_norm_Fenhance_stroke.^2 + Fy_trans_norm_Fenhance_stroke.^2 + Fz_trans_norm_Fenhance_stroke.^2);
Mtot_trans_norm_Fenhance_stroke = sqrt(Mx_trans_norm_Fenhance_stroke.^2 + My_trans_norm_Fenhance_stroke.^2 + Mz_trans_norm_Fenhance_stroke.^2);
Ftot_transNrot_norm_Fenhance_stroke = sqrt(Fx_transNrot_norm_Fenhance_stroke.^2 + Fy_transNrot_norm_Fenhance_stroke.^2 + Fz_transNrot_norm_Fenhance_stroke.^2);
Mtot_transNrot_norm_Fenhance_stroke = sqrt(Mx_transNrot_norm_Fenhance_stroke.^2 + My_transNrot_norm_Fenhance_stroke.^2 + Mz_transNrot_norm_Fenhance_stroke.^2);

Ftot_trans_norm_Fenhance_dev = sqrt(Fx_trans_norm_Fenhance_dev.^2 + Fy_trans_norm_Fenhance_dev.^2 + Fz_trans_norm_Fenhance_dev.^2);
Mtot_trans_norm_Fenhance_dev = sqrt(Mx_trans_norm_Fenhance_dev.^2 + My_trans_norm_Fenhance_dev.^2 + Mz_trans_norm_Fenhance_dev.^2);
Ftot_transNrot_norm_Fenhance_dev = sqrt(Fx_transNrot_norm_Fenhance_dev.^2 + Fy_transNrot_norm_Fenhance_dev.^2 + Fz_transNrot_norm_Fenhance_dev.^2);
Mtot_transNrot_norm_Fenhance_dev = sqrt(Mx_transNrot_norm_Fenhance_dev.^2 + My_transNrot_norm_Fenhance_dev.^2 + Mz_transNrot_norm_Fenhance_dev.^2);

Ftot_trans_norm_Fenhance_rot = sqrt(Fx_trans_norm_Fenhance_rot.^2 + Fy_trans_norm_Fenhance_rot.^2 + Fz_trans_norm_Fenhance_rot.^2);
Mtot_trans_norm_Fenhance_rot = sqrt(Mx_trans_norm_Fenhance_rot.^2 + My_trans_norm_Fenhance_rot.^2 + Mz_trans_norm_Fenhance_rot.^2);
Ftot_transNrot_norm_Fenhance_rot = sqrt(Fx_transNrot_norm_Fenhance_rot.^2 + Fy_transNrot_norm_Fenhance_rot.^2 + Fz_transNrot_norm_Fenhance_rot.^2);
Mtot_transNrot_norm_Fenhance_rot = sqrt(Mx_transNrot_norm_Fenhance_rot.^2 + My_transNrot_norm_Fenhance_rot.^2 + Mz_transNrot_norm_Fenhance_rot.^2);

save('FnMqs_Fenhance.mat','mod_value_Fenhance',...
    'Ftot_trans_norm_Fenhance_all','Mtot_trans_norm_Fenhance_all',...
    'Ftot_transNrot_norm_Fenhance_all','Mtot_transNrot_norm_Fenhance_all',...
    'Fx_trans_norm_Fenhance_all','Fy_trans_norm_Fenhance_all','Fz_trans_norm_Fenhance_all',...
    'Mx_trans_norm_Fenhance_all','My_trans_norm_Fenhance_all','Mz_trans_norm_Fenhance_all',...
    'Fx_transNrot_norm_Fenhance_all','Fy_transNrot_norm_Fenhance_all','Fz_transNrot_norm_Fenhance_all',...
    'Mx_transNrot_norm_Fenhance_all','My_transNrot_norm_Fenhance_all','Mz_transNrot_norm_Fenhance_all',...
    'Ftot_trans_norm_Fenhance_freq','Mtot_trans_norm_Fenhance_freq',...
    'Ftot_transNrot_norm_Fenhance_freq','Mtot_transNrot_norm_Fenhance_freq',...
    'Fx_trans_norm_Fenhance_freq','Fy_trans_norm_Fenhance_freq','Fz_trans_norm_Fenhance_freq',...
    'Mx_trans_norm_Fenhance_freq','My_trans_norm_Fenhance_freq','Mz_trans_norm_Fenhance_freq',...
    'Fx_transNrot_norm_Fenhance_freq','Fy_transNrot_norm_Fenhance_freq','Fz_transNrot_norm_Fenhance_freq',...
    'Mx_transNrot_norm_Fenhance_freq','My_transNrot_norm_Fenhance_freq','Mz_transNrot_norm_Fenhance_freq',...
    'Ftot_trans_norm_Fenhance_stroke','Mtot_trans_norm_Fenhance_stroke',...
    'Ftot_transNrot_norm_Fenhance_stroke','Mtot_transNrot_norm_Fenhance_stroke',...
    'Fx_trans_norm_Fenhance_stroke','Fy_trans_norm_Fenhance_stroke','Fz_trans_norm_Fenhance_stroke',...
    'Mx_trans_norm_Fenhance_stroke','My_trans_norm_Fenhance_stroke','Mz_trans_norm_Fenhance_stroke',...
    'Fx_transNrot_norm_Fenhance_stroke','Fy_transNrot_norm_Fenhance_stroke','Fz_transNrot_norm_Fenhance_stroke',...
    'Mx_transNrot_norm_Fenhance_stroke','My_transNrot_norm_Fenhance_stroke','Mz_transNrot_norm_Fenhance_stroke',...
    'Ftot_trans_norm_Fenhance_dev','Mtot_trans_norm_Fenhance_dev',...
    'Ftot_transNrot_norm_Fenhance_dev','Mtot_transNrot_norm_Fenhance_dev',...
    'Fx_trans_norm_Fenhance_dev','Fy_trans_norm_Fenhance_dev','Fz_trans_norm_Fenhance_dev',...
    'Mx_trans_norm_Fenhance_dev','My_trans_norm_Fenhance_dev','Mz_trans_norm_Fenhance_dev',...
    'Fx_transNrot_norm_Fenhance_dev','Fy_transNrot_norm_Fenhance_dev','Fz_transNrot_norm_Fenhance_dev',...
    'Mx_transNrot_norm_Fenhance_dev','My_transNrot_norm_Fenhance_dev','Mz_transNrot_norm_Fenhance_dev',...
    'Ftot_trans_norm_Fenhance_rot','Mtot_trans_norm_Fenhance_rot',...
    'Ftot_transNrot_norm_Fenhance_rot','Mtot_transNrot_norm_Fenhance_rot',...
    'Fx_trans_norm_Fenhance_rot','Fy_trans_norm_Fenhance_rot','Fz_trans_norm_Fenhance_rot',...
    'Mx_trans_norm_Fenhance_rot','My_trans_norm_Fenhance_rot','Mz_trans_norm_Fenhance_rot',...
    'Fx_transNrot_norm_Fenhance_rot','Fy_transNrot_norm_Fenhance_rot','Fz_transNrot_norm_Fenhance_rot',...
    'Mx_transNrot_norm_Fenhance_rot','My_transNrot_norm_Fenhance_rot','Mz_transNrot_norm_Fenhance_rot')


