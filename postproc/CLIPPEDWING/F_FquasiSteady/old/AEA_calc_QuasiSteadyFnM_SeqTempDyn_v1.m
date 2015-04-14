clear;
clc;
close all
warning off

addpath('/home/florian/Dropbox/WORK/flytracker/flytracker');
addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/mex/');
addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/core/');
addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/results/');

load('BodyWingModel_WingSections_TotalMean.mat')

loadname=dir('WBdataset_temporal_dynamics_*')
loadname = loadname.name;
load(loadname)

load('norm_data_torque.mat')

plot_on = 1
% plot_on = 0
if plot_on == 1
    mkdir('MSfigs_QSmodel_TempDynamics')
    cd('MSfigs_QSmodel_TempDynamics')
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

nr_timepoints = length(t_norm_wb_seq_bins_mean_all);

%% body dynamics
% NO body vel&rot
body_model.u_strk(:,nr_timepoints) = [0 0 0];
body_model.w_strk(:,nr_timepoints) = [0 0 0];

R_strk = body_model.R_strk;
u_strk = body_model.u_strk;
w_strk = body_model.w_strk;

%% WBkin & F&M
    calc_QSforcesNmoments_LnR_timeseries

    Fx_trans_norm_timeseries = FM_trans_norm(1,:);
    Fy_trans_norm_timeseries = FM_trans_norm(2,:);
    Fz_trans_norm_timeseries = FM_trans_norm(3,:);
    
    Mx_trans_norm_timeseries = FM_trans_norm(4,:);
    My_trans_norm_timeseries = FM_trans_norm(5,:);
    Mz_trans_norm_timeseries = FM_trans_norm(6,:);
    
    Fx_transNrot_norm_timeseries = FM_transNrot_norm(1,:);
    Fy_transNrot_norm_timeseries = FM_transNrot_norm(2,:);
    Fz_transNrot_norm_timeseries = FM_transNrot_norm(3,:);
    
    Mx_transNrot_norm_timeseries = FM_transNrot_norm(4,:);
    My_transNrot_norm_timeseries = FM_transNrot_norm(5,:);
    Mz_transNrot_norm_timeseries = FM_transNrot_norm(6,:);

    Fx_trans_norm_timeseries_WBmean = FM_trans_norm_WBmean(:,1);
    Fy_trans_norm_timeseries_WBmean = FM_trans_norm_WBmean(:,2);
    Fz_trans_norm_timeseries_WBmean = FM_trans_norm_WBmean(:,3);
    
    Mx_trans_norm_timeseries_WBmean = FM_trans_norm_WBmean(:,4);
    My_trans_norm_timeseries_WBmean = FM_trans_norm_WBmean(:,5);
    Mz_trans_norm_timeseries_WBmean = FM_trans_norm_WBmean(:,6);
    
    Fx_transNrot_norm_timeseries_WBmean = FM_transNrot_norm_WBmean(:,1);
    Fy_transNrot_norm_timeseries_WBmean = FM_transNrot_norm_WBmean(:,2);
    Fz_transNrot_norm_timeseries_WBmean = FM_transNrot_norm_WBmean(:,3);
    
    Mx_transNrot_norm_timeseries_WBmean = FM_transNrot_norm_WBmean(:,4);
    My_transNrot_norm_timeseries_WBmean = FM_transNrot_norm_WBmean(:,5);
    Mz_transNrot_norm_timeseries_WBmean = FM_transNrot_norm_WBmean(:,6);
    
    %% set pitch torque at start to zero
    My_trans_norm_timeseries_WBmean_NOcomp = My_trans_norm_timeseries_WBmean;
    My_transNrot_norm_timeseries_WBmean_NOcomp = My_transNrot_norm_timeseries_WBmean;
    
    My0 = My_trans_norm_timeseries_WBmean(isnan(My_trans_norm_timeseries_WBmean)==0);
    My0 = My0(1);
    My_trans_norm_timeseries_WBmean = My_trans_norm_timeseries_WBmean - My0;
    
    My0 = My_transNrot_norm_timeseries_WBmean(isnan(My_transNrot_norm_timeseries_WBmean)==0);
    My0 = My0(1);
    My_transNrot_norm_timeseries_WBmean = My_transNrot_norm_timeseries_WBmean - My0;
    
    %% axisRnL
    F_R_trans_norm_timeseries_WBmean = Fx_trans_norm_timeseries_WBmean * cosd(rot_axis_angle) + ...
        Fy_trans_norm_timeseries_WBmean * sind(rot_axis_angle);
    F_L_trans_norm_timeseries_WBmean = -Fx_trans_norm_timeseries_WBmean * sind(rot_axis_angle) + ...
        Fy_trans_norm_timeseries_WBmean * cosd(rot_axis_angle);

    F_R_transNrot_norm_timeseries_WBmean = Fx_transNrot_norm_timeseries_WBmean * cosd(rot_axis_angle) + ...
        Fy_transNrot_norm_timeseries_WBmean * sind(rot_axis_angle);
    F_L_transNrot_norm_timeseries_WBmean = -Fx_transNrot_norm_timeseries_WBmean * sind(rot_axis_angle) + ...
        Fy_transNrot_norm_timeseries_WBmean * cosd(rot_axis_angle);

    F_R_trans_norm_timeseries = Fx_trans_norm_timeseries * cosd(rot_axis_angle) + ...
        Fy_trans_norm_timeseries * sind(rot_axis_angle);
    F_L_trans_norm_timeseries = -Fx_trans_norm_timeseries * sind(rot_axis_angle) + ...
        Fy_trans_norm_timeseries * cosd(rot_axis_angle);

    F_R_transNrot_norm_timeseries = Fx_transNrot_norm_timeseries * cosd(rot_axis_angle) + ...
        Fy_transNrot_norm_timeseries * sind(rot_axis_angle);
    F_L_transNrot_norm_timeseries = -Fx_transNrot_norm_timeseries * sind(rot_axis_angle) + ...
        Fy_transNrot_norm_timeseries * cosd(rot_axis_angle);
    
    M_R_trans_norm_timeseries_WBmean = Mx_trans_norm_timeseries_WBmean * cosd(rot_axis_angle) + ...
        My_trans_norm_timeseries_WBmean * sind(rot_axis_angle);
    M_L_trans_norm_timeseries_WBmean = -Mx_trans_norm_timeseries_WBmean * sind(rot_axis_angle) + ...
        My_trans_norm_timeseries_WBmean * cosd(rot_axis_angle);

    M_R_transNrot_norm_timeseries_WBmean = Mx_transNrot_norm_timeseries_WBmean * cosd(rot_axis_angle) + ...
        My_transNrot_norm_timeseries_WBmean * sind(rot_axis_angle);
    M_L_transNrot_norm_timeseries_WBmean = -Mx_transNrot_norm_timeseries_WBmean * sind(rot_axis_angle) + ...
        My_transNrot_norm_timeseries_WBmean * cosd(rot_axis_angle);

    M_R_trans_norm_timeseries = Mx_trans_norm_timeseries * cosd(rot_axis_angle) + ...
        My_trans_norm_timeseries * sind(rot_axis_angle);
    M_L_trans_norm_timeseries = -Mx_trans_norm_timeseries * sind(rot_axis_angle) + ...
        My_trans_norm_timeseries * cosd(rot_axis_angle);

    M_R_transNrot_norm_timeseries = Mx_transNrot_norm_timeseries * cosd(rot_axis_angle) + ...
        My_transNrot_norm_timeseries * sind(rot_axis_angle);
    M_L_transNrot_norm_timeseries = -Mx_transNrot_norm_timeseries * sind(rot_axis_angle) + ...
        My_transNrot_norm_timeseries * cosd(rot_axis_angle);
    
if plot_on == 1
    cd ..
end

%% save data
save('FnMqs_timeseries.mat',...
    't_wb_seq_bins_mean_all','t_wb_seq_mean_all',...
    'Fx_trans_norm_timeseries','Fy_trans_norm_timeseries','Fz_trans_norm_timeseries',...
    'Mx_trans_norm_timeseries','My_trans_norm_timeseries','Mz_trans_norm_timeseries',...
    'Fx_transNrot_norm_timeseries','Fy_transNrot_norm_timeseries','Fz_transNrot_norm_timeseries',...
    'Mx_transNrot_norm_timeseries','My_transNrot_norm_timeseries','Mz_transNrot_norm_timeseries',...
    'Fx_trans_norm_timeseries_WBmean','Fy_trans_norm_timeseries_WBmean','Fz_trans_norm_timeseries_WBmean',...
    'Mx_trans_norm_timeseries_WBmean','My_trans_norm_timeseries_WBmean','Mz_trans_norm_timeseries_WBmean',...
    'Fx_transNrot_norm_timeseries_WBmean','Fy_transNrot_norm_timeseries_WBmean','Fz_transNrot_norm_timeseries_WBmean',...
    'Mx_transNrot_norm_timeseries_WBmean','My_transNrot_norm_timeseries_WBmean','Mz_transNrot_norm_timeseries_WBmean',...
    'F_R_trans_norm_timeseries','F_L_trans_norm_timeseries',...
    'M_R_trans_norm_timeseries','M_L_trans_norm_timeseries',...
    'F_R_transNrot_norm_timeseries','F_L_transNrot_norm_timeseries',...
    'M_R_transNrot_norm_timeseries','M_L_transNrot_norm_timeseries',...
    'F_R_trans_norm_timeseries_WBmean','F_L_trans_norm_timeseries_WBmean',...
    'M_R_trans_norm_timeseries_WBmean','M_L_trans_norm_timeseries_WBmean',...
    'F_R_transNrot_norm_timeseries_WBmean','F_L_transNrot_norm_timeseries_WBmean',...
    'M_R_transNrot_norm_timeseries_WBmean','M_L_transNrot_norm_timeseries_WBmean',...
    'My_trans_norm_timeseries_WBmean_NOcomp','My_transNrot_norm_timeseries_WBmean_NOcomp');





















