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

mod_min = 0;
mod_max = 1;
dmod = .5;

mod_value = [mod_min:dmod:mod_max]';

for i = 1:length(mod_value)
    mod_value_now = mod_value(i);
    
    %% ALL MODS
    % set freq FIRST (same for MODnSteady)
    freq = f_wb_steady + mod_value_now*freqMOD_wb_Fenhance_meanCIstd(1);

    % WBmod fouriers
    stroke_coeffs_MOD_L = strokeMOD_Fenhance_fourier_coeffs_binmean;
    dev_coeffs_MOD_L = devMOD_Fenhance_fourier_coeffs_binmean;
    rot_coeffs_MOD_L = pitchMOD_Fenhance_fourier_coeffs_binmean;
    
    stroke_coeffs_MOD_R = strokeMOD_Fenhance_fourier_coeffs_binmean;
    dev_coeffs_MOD_R = devMOD_Fenhance_fourier_coeffs_binmean;
    rot_coeffs_MOD_R = pitchMOD_Fenhance_fourier_coeffs_binmean;
    
    % calc FnM SYMMETRIC & MOD WB
    [FM_trans,FM_trans_norm,FM_transNrot,FM_transNrot_norm,...
        FM_L_trans,FM_L_trans_norm,FM_L_transNrot,FM_L_transNrot_norm,...
        FM_R_trans,FM_R_trans_norm,FM_R_transNrot,FM_R_transNrot_norm] =...
        calc_QSforcesNmoments_TESTnModNbodyKin_LnR_TEST(mod_value_now,const,wing_model,body_model,freq,...
        stroke_coeffs_steady,dev_coeffs_steady,rot_coeffs_steady,...
        stroke_coeffs_MOD_L,dev_coeffs_MOD_L,rot_coeffs_MOD_L,...
        stroke_coeffs_MOD_R,dev_coeffs_MOD_R,rot_coeffs_MOD_R);
    
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

end

mod_value_Fenhance = mod_value;

% mean values
Fx_trans_norm_Fenhance_all_mean = mean(Fx_trans_norm_Fenhance_all);
Fy_trans_norm_Fenhance_all_mean = mean(Fy_trans_norm_Fenhance_all);
Fz_trans_norm_Fenhance_all_mean = mean(Fz_trans_norm_Fenhance_all);

Mx_trans_norm_Fenhance_all_mean = mean(Mx_trans_norm_Fenhance_all);
My_trans_norm_Fenhance_all_mean = mean(My_trans_norm_Fenhance_all);
Mz_trans_norm_Fenhance_all_mean = mean(Mz_trans_norm_Fenhance_all);

Fx_transNrot_norm_Fenhance_all_mean = mean(Fx_transNrot_norm_Fenhance_all);
Fy_transNrot_norm_Fenhance_all_mean = mean(Fy_transNrot_norm_Fenhance_all);
Fz_transNrot_norm_Fenhance_all_mean = mean(Fz_transNrot_norm_Fenhance_all);

Mx_transNrot_norm_Fenhance_all_mean = mean(Mx_transNrot_norm_Fenhance_all);
My_transNrot_norm_Fenhance_all_mean = mean(My_transNrot_norm_Fenhance_all);
Mz_transNrot_norm_Fenhance_all_mean = mean(Mz_transNrot_norm_Fenhance_all);


save('FnMqs_Fenhance_TEST.mat','mod_value_Fenhance',...
    'Fx_trans_norm_Fenhance_all','Fy_trans_norm_Fenhance_all','Fz_trans_norm_Fenhance_all',...
    'Mx_trans_norm_Fenhance_all','My_trans_norm_Fenhance_all','Mz_trans_norm_Fenhance_all',...
    'Fx_transNrot_norm_Fenhance_all','Fy_transNrot_norm_Fenhance_all','Fz_transNrot_norm_Fenhance_all',...
    'Mx_transNrot_norm_Fenhance_all','My_transNrot_norm_Fenhance_all','Mz_transNrot_norm_Fenhance_all')

