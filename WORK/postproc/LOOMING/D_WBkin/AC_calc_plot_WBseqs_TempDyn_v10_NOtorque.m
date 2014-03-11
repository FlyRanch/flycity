clc
clear
close all

addpath('/home/florian/Dropbox/WORK/toolbox')
addpath('/home/florian/Dropbox/WORK/toolbox/plot2svg')
addpath('/home/florian/Dropbox/WORK/toolbox/CircStat')
addpath('/home/florian/Dropbox/WORK/toolbox/flytracker')

%% load data
loadname=dir('WBdataset_all_*')
loadname=loadname.name;
load(loadname)

% steady_name=dir('WBdataset_steady_*')
% steady_name=steady_name.name;
% load(steady_name)

const_name=dir('flyVars*')
const_name=const_name.name;
load(const_name)

load('clusters.mat')

%% settings
n_bins = size(stroke_wb_L_bins,1);
wb_max = 30;
wb_min4mean = max(seq_nr)/3; % only mean if half of the sequences are providing data
wb_min4mean = 10; % only mean if half of the sequences are providing data
wb_min4mean_clusters = 5;

f_wb = (f_wb_L + f_wb_R)./2;

define_vars_timeline

%% align
for seq_now = 1:max(seq_nr)
    
    frame_nr_2resp_now = frame_nr_2resp(seq_nr==seq_now);
    frame_nr_2resp_now_pre = frame_nr_2resp_now(frame_nr_2resp_now<0);
    frame_nr_2resp_now_post = frame_nr_2resp_now(frame_nr_2resp_now>=0);
    
    %% wb kin
    f_wb_now = f_wb(seq_nr==seq_now);
    f_wb_now_pre = f_wb_now(frame_nr_2resp_now<0);
    f_wb_now_post = f_wb_now(frame_nr_2resp_now>=0);
    
    stroke_wb_L_bins_now = stroke_wb_L_bins(:,seq_nr==seq_now);
    stroke_wb_L_bins_now_pre = stroke_wb_L_bins_now(:,frame_nr_2resp_now<0);
    stroke_wb_L_bins_now_post = stroke_wb_L_bins_now(:,frame_nr_2resp_now>=0);
    
    stroke_wb_R_bins_now = stroke_wb_R_bins(:,seq_nr==seq_now);
    stroke_wb_R_bins_now_pre = stroke_wb_R_bins_now(:,frame_nr_2resp_now<0);
    stroke_wb_R_bins_now_post = stroke_wb_R_bins_now(:,frame_nr_2resp_now>=0);
    
    dev_wb_L_bins_now = dev_wb_L_bins(:,seq_nr==seq_now);
    dev_wb_L_bins_now_pre = dev_wb_L_bins_now(:,frame_nr_2resp_now<0);
    dev_wb_L_bins_now_post = dev_wb_L_bins_now(:,frame_nr_2resp_now>=0);
    
    dev_wb_R_bins_now = dev_wb_R_bins(:,seq_nr==seq_now);
    dev_wb_R_bins_now_pre = dev_wb_R_bins_now(:,frame_nr_2resp_now<0);
    dev_wb_R_bins_now_post = dev_wb_R_bins_now(:,frame_nr_2resp_now>=0);
    
    pitch_wb_L_bins_now = pitch_wb_L_bins(:,seq_nr==seq_now);
    pitch_wb_L_bins_now_pre = pitch_wb_L_bins_now(:,frame_nr_2resp_now<0);
    pitch_wb_L_bins_now_post = pitch_wb_L_bins_now(:,frame_nr_2resp_now>=0);
    
    pitch_wb_R_bins_now = pitch_wb_R_bins(:,seq_nr==seq_now);
    pitch_wb_R_bins_now_pre = pitch_wb_R_bins_now(:,frame_nr_2resp_now<0);
    pitch_wb_R_bins_now_post = pitch_wb_R_bins_now(:,frame_nr_2resp_now>=0);
    
    %% body kin
    roll_mean_wb_now = roll_mean_wb(seq_nr==seq_now);
    roll_mean_wb_now_pre = roll_mean_wb_now(frame_nr_2resp_now<0);
    roll_mean_wb_now_post = roll_mean_wb_now(frame_nr_2resp_now>=0);
    
    roll_dot_mean_wb_now = roll_dot_mean_wb(seq_nr==seq_now);
    roll_dot_mean_wb_now_pre = roll_dot_mean_wb_now(frame_nr_2resp_now<0);
    roll_dot_mean_wb_now_post = roll_dot_mean_wb_now(frame_nr_2resp_now>=0);
    
    roll_dot_dot_mean_wb_now = roll_dot_dot_mean_wb(seq_nr==seq_now);
    roll_dot_dot_mean_wb_now_pre = roll_dot_dot_mean_wb_now(frame_nr_2resp_now<0);
    roll_dot_dot_mean_wb_now_post = roll_dot_dot_mean_wb_now(frame_nr_2resp_now>=0);
    
    pitch_mean_wb_now = pitch_mean_wb(seq_nr==seq_now);
    pitch_mean_wb_now_pre = pitch_mean_wb_now(frame_nr_2resp_now<0);
    pitch_mean_wb_now_post = pitch_mean_wb_now(frame_nr_2resp_now>=0);
    
    pitch_dot_mean_wb_now = pitch_dot_mean_wb(seq_nr==seq_now);
    pitch_dot_mean_wb_now_pre = pitch_dot_mean_wb_now(frame_nr_2resp_now<0);
    pitch_dot_mean_wb_now_post = pitch_dot_mean_wb_now(frame_nr_2resp_now>=0);
    
    pitch_dot_dot_mean_wb_now = pitch_dot_dot_mean_wb(seq_nr==seq_now);
    pitch_dot_dot_mean_wb_now_pre = pitch_dot_dot_mean_wb_now(frame_nr_2resp_now<0);
    pitch_dot_dot_mean_wb_now_post = pitch_dot_dot_mean_wb_now(frame_nr_2resp_now>=0);
    
    yaw_mean_wb_now = yaw_mean_wb(seq_nr==seq_now);
    yaw_mean_wb_now_pre = yaw_mean_wb_now(frame_nr_2resp_now<0);
    yaw_mean_wb_now_post = yaw_mean_wb_now(frame_nr_2resp_now>=0);
    
    yaw_dot_mean_wb_now = yaw_dot_mean_wb(seq_nr==seq_now);
    yaw_dot_mean_wb_now_pre = yaw_dot_mean_wb_now(frame_nr_2resp_now<0);
    yaw_dot_mean_wb_now_post = yaw_dot_mean_wb_now(frame_nr_2resp_now>=0);
    
    yaw_dot_dot_mean_wb_now = yaw_dot_dot_mean_wb(seq_nr==seq_now);
    yaw_dot_dot_mean_wb_now_pre = yaw_dot_dot_mean_wb_now(frame_nr_2resp_now<0);
    yaw_dot_dot_mean_wb_now_post = yaw_dot_dot_mean_wb_now(frame_nr_2resp_now>=0);
    
%     drot_L_mean_wb_now = drot_L_mean_wb(seq_nr==seq_now);
%     drot_L_mean_wb_now_pre = drot_L_mean_wb_now(frame_nr_2resp_now<0);
%     drot_L_mean_wb_now_post = drot_L_mean_wb_now(frame_nr_2resp_now>=0);
%     
%     rot_dot_L_mean_wb_now = rot_dot_L_mean_wb(seq_nr==seq_now);
%     rot_dot_L_mean_wb_now_pre = rot_dot_L_mean_wb_now(frame_nr_2resp_now<0);
%     rot_dot_L_mean_wb_now_post = rot_dot_L_mean_wb_now(frame_nr_2resp_now>=0);
%     
%     rot_dot_dot_L_mean_wb_now = rot_dot_dot_L_mean_wb(seq_nr==seq_now);
%     rot_dot_dot_L_mean_wb_now_pre = rot_dot_dot_L_mean_wb_now(frame_nr_2resp_now<0);
%     rot_dot_dot_L_mean_wb_now_post = rot_dot_dot_L_mean_wb_now(frame_nr_2resp_now>=0);
%     
%     drot_R_mean_wb_now = drot_R_mean_wb(seq_nr==seq_now);
%     drot_R_mean_wb_now_pre = drot_R_mean_wb_now(frame_nr_2resp_now<0);
%     drot_R_mean_wb_now_post = drot_R_mean_wb_now(frame_nr_2resp_now>=0);
%     
%     rot_dot_R_mean_wb_now = rot_dot_R_mean_wb(seq_nr==seq_now);
%     rot_dot_R_mean_wb_now_pre = rot_dot_R_mean_wb_now(frame_nr_2resp_now<0);
%     rot_dot_R_mean_wb_now_post = rot_dot_R_mean_wb_now(frame_nr_2resp_now>=0);
%     
%     rot_dot_dot_R_mean_wb_now = rot_dot_dot_R_mean_wb(seq_nr==seq_now);
%     rot_dot_dot_R_mean_wb_now_pre = rot_dot_dot_R_mean_wb_now(frame_nr_2resp_now<0);
%     rot_dot_dot_R_mean_wb_now_post = rot_dot_dot_R_mean_wb_now(frame_nr_2resp_now>=0);
% 
%     %% torque
%     Mroll_mean_wb_accel_now = Mroll_mean_wb_accel(seq_nr==seq_now);
%     Mroll_mean_wb_accel_now_pre = Mroll_mean_wb_accel_now(frame_nr_2resp_now<0);
%     Mroll_mean_wb_accel_now_post = Mroll_mean_wb_accel_now(frame_nr_2resp_now>=0);
%     
%     Mroll_mean_wb_damp_now = Mroll_mean_wb_damp(seq_nr==seq_now);
%     Mroll_mean_wb_damp_now_pre = Mroll_mean_wb_damp_now(frame_nr_2resp_now<0);
%     Mroll_mean_wb_damp_now_post = Mroll_mean_wb_damp_now(frame_nr_2resp_now>=0);
%     
%     Mroll_mean_wb_now = Mroll_mean_wb(seq_nr==seq_now);
%     Mroll_mean_wb_now_pre = Mroll_mean_wb_now(frame_nr_2resp_now<0);
%     Mroll_mean_wb_now_post = Mroll_mean_wb_now(frame_nr_2resp_now>=0);
%     
%     Mpitch_mean_wb_accel_now = Mpitch_mean_wb_accel(seq_nr==seq_now);
%     Mpitch_mean_wb_accel_now_pre = Mpitch_mean_wb_accel_now(frame_nr_2resp_now<0);
%     Mpitch_mean_wb_accel_now_post = Mpitch_mean_wb_accel_now(frame_nr_2resp_now>=0);
%     
%     Mpitch_mean_wb_damp_now = Mpitch_mean_wb_damp(seq_nr==seq_now);
%     Mpitch_mean_wb_damp_now_pre = Mpitch_mean_wb_damp_now(frame_nr_2resp_now<0);
%     Mpitch_mean_wb_damp_now_post = Mpitch_mean_wb_damp_now(frame_nr_2resp_now>=0);
%     
%     Mpitch_mean_wb_now = Mpitch_mean_wb(seq_nr==seq_now);
%     Mpitch_mean_wb_now_pre = Mpitch_mean_wb_now(frame_nr_2resp_now<0);
%     Mpitch_mean_wb_now_post = Mpitch_mean_wb_now(frame_nr_2resp_now>=0);
%     
%     Myaw_mean_wb_accel_now = Myaw_mean_wb_accel(seq_nr==seq_now);
%     Myaw_mean_wb_accel_now_pre = Myaw_mean_wb_accel_now(frame_nr_2resp_now<0);
%     Myaw_mean_wb_accel_now_post = Myaw_mean_wb_accel_now(frame_nr_2resp_now>=0);
%     
%     Myaw_mean_wb_damp_now = Myaw_mean_wb_damp(seq_nr==seq_now);
%     Myaw_mean_wb_damp_now_pre = Myaw_mean_wb_damp_now(frame_nr_2resp_now<0);
%     Myaw_mean_wb_damp_now_post = Myaw_mean_wb_damp_now(frame_nr_2resp_now>=0);
%     
%     Myaw_mean_wb_now = Myaw_mean_wb(seq_nr==seq_now);
%     Myaw_mean_wb_now_pre = Myaw_mean_wb_now(frame_nr_2resp_now<0);
%     Myaw_mean_wb_now_post = Myaw_mean_wb_now(frame_nr_2resp_now>=0);
%     
%     M_L_mean_wb_accel_now = M_L_mean_wb_accel(seq_nr==seq_now);
%     M_L_mean_wb_accel_now_pre = M_L_mean_wb_accel_now(frame_nr_2resp_now<0);
%     M_L_mean_wb_accel_now_post = M_L_mean_wb_accel_now(frame_nr_2resp_now>=0);
%     
%     M_L_mean_wb_damp_now = M_L_mean_wb_damp(seq_nr==seq_now);
%     M_L_mean_wb_damp_now_pre = M_L_mean_wb_damp_now(frame_nr_2resp_now<0);
%     M_L_mean_wb_damp_now_post = M_L_mean_wb_damp_now(frame_nr_2resp_now>=0);
%     
%     M_L_mean_wb_now = M_L_mean_wb(seq_nr==seq_now);
%     M_L_mean_wb_now_pre = M_L_mean_wb_now(frame_nr_2resp_now<0);
%     M_L_mean_wb_now_post = M_L_mean_wb_now(frame_nr_2resp_now>=0);
%     
%     M_R_mean_wb_accel_now = M_R_mean_wb_accel(seq_nr==seq_now);
%     M_R_mean_wb_accel_now_pre = M_R_mean_wb_accel_now(frame_nr_2resp_now<0);
%     M_R_mean_wb_accel_now_post = M_R_mean_wb_accel_now(frame_nr_2resp_now>=0);
%     
%     M_R_mean_wb_damp_now = M_R_mean_wb_damp(seq_nr==seq_now);
%     M_R_mean_wb_damp_now_pre = M_R_mean_wb_damp_now(frame_nr_2resp_now<0);
%     M_R_mean_wb_damp_now_post = M_R_mean_wb_damp_now(frame_nr_2resp_now>=0);
%     
%     M_R_mean_wb_now = M_R_mean_wb(seq_nr==seq_now);
%     M_R_mean_wb_now_pre = M_R_mean_wb_now(frame_nr_2resp_now<0);
%     M_R_mean_wb_now_post = M_R_mean_wb_now(frame_nr_2resp_now>=0);
    
    %% pre response wbs
    wb_pre_max = length(frame_nr_2resp_now_pre);
    for wb_now = 1:wb_pre_max
        
        wb_pre = wb_max + wb_now - wb_pre_max;
        
        if wb_pre > 0
        % wb kin
        f_wb_seq_pre(seq_now,wb_pre) = f_wb_now_pre(wb_now);
        
        stroke_wb_L_seq_bins_pre(:,seq_now,wb_pre) = stroke_wb_L_bins_now_pre(:,wb_now);
        stroke_wb_R_seq_bins_pre(:,seq_now,wb_pre) = stroke_wb_R_bins_now_pre(:,wb_now);
        
        dev_wb_L_seq_bins_pre(:,seq_now,wb_pre) = dev_wb_L_bins_now_pre(:,wb_now);
        dev_wb_R_seq_bins_pre(:,seq_now,wb_pre) = dev_wb_R_bins_now_pre(:,wb_now);
        
        pitch_wb_L_seq_bins_pre(:,seq_now,wb_pre) = pitch_wb_L_bins_now_pre(:,wb_now);
        pitch_wb_R_seq_bins_pre(:,seq_now,wb_pre) = pitch_wb_R_bins_now_pre(:,wb_now);
        
        % body kin
        roll_mean_wb_seq_pre(seq_now,wb_pre) = roll_mean_wb_now_pre(wb_now);
        roll_dot_mean_wb_seq_pre(seq_now,wb_pre) = roll_dot_mean_wb_now_pre(wb_now);
        roll_dot_dot_mean_wb_seq_pre(seq_now,wb_pre) = roll_dot_dot_mean_wb_now_pre(wb_now);
        
        pitch_mean_wb_seq_pre(seq_now,wb_pre) = pitch_mean_wb_now_pre(wb_now);
        pitch_dot_mean_wb_seq_pre(seq_now,wb_pre) = pitch_dot_mean_wb_now_pre(wb_now);
        pitch_dot_dot_mean_wb_seq_pre(seq_now,wb_pre) = pitch_dot_dot_mean_wb_now_pre(wb_now);
        
        yaw_mean_wb_seq_pre(seq_now,wb_pre) = yaw_mean_wb_now_pre(wb_now);
        yaw_dot_mean_wb_seq_pre(seq_now,wb_pre) = yaw_dot_mean_wb_now_pre(wb_now);
        yaw_dot_dot_mean_wb_seq_pre(seq_now,wb_pre) = yaw_dot_dot_mean_wb_now_pre(wb_now);
        
%         drot_L_mean_wb_seq_pre(seq_now,wb_pre) = drot_L_mean_wb_now_pre(wb_now);
%         rot_dot_L_mean_wb_seq_pre(seq_now,wb_pre) = rot_dot_L_mean_wb_now_pre(wb_now);
%         rot_dot_dot_L_mean_wb_seq_pre(seq_now,wb_pre) = rot_dot_dot_L_mean_wb_now_pre(wb_now);
%         
%         drot_R_mean_wb_seq_pre(seq_now,wb_pre) = drot_R_mean_wb_now_pre(wb_now);
%         rot_dot_R_mean_wb_seq_pre(seq_now,wb_pre) = rot_dot_R_mean_wb_now_pre(wb_now);
%         rot_dot_dot_R_mean_wb_seq_pre(seq_now,wb_pre) = rot_dot_dot_R_mean_wb_now_pre(wb_now);
%         
%         % torque
%         Mroll_mean_wb_accel_seq_pre(seq_now,wb_pre) = Mroll_mean_wb_accel_now_pre(wb_now);
%         Mroll_mean_wb_damp_seq_pre(seq_now,wb_pre) = Mroll_mean_wb_damp_now_pre(wb_now);
%         Mroll_mean_wb_seq_pre(seq_now,wb_pre) = Mroll_mean_wb_now_pre(wb_now);
%         
%         Mpitch_mean_wb_accel_seq_pre(seq_now,wb_pre) = Mpitch_mean_wb_accel_now_pre(wb_now);
%         Mpitch_mean_wb_damp_seq_pre(seq_now,wb_pre) = Mpitch_mean_wb_damp_now_pre(wb_now);
%         Mpitch_mean_wb_seq_pre(seq_now,wb_pre) = Mpitch_mean_wb_now_pre(wb_now);
% 
%         Myaw_mean_wb_accel_seq_pre(seq_now,wb_pre) = Myaw_mean_wb_accel_now_pre(wb_now);
%         Myaw_mean_wb_damp_seq_pre(seq_now,wb_pre) = Myaw_mean_wb_damp_now_pre(wb_now);
%         Myaw_mean_wb_seq_pre(seq_now,wb_pre) = Myaw_mean_wb_now_pre(wb_now);
%         
%         M_L_mean_wb_accel_seq_pre(seq_now,wb_pre) = M_L_mean_wb_accel_now_pre(wb_now);
%         M_L_mean_wb_damp_seq_pre(seq_now,wb_pre) = M_L_mean_wb_damp_now_pre(wb_now);
%         M_L_mean_wb_seq_pre(seq_now,wb_pre) = M_L_mean_wb_now_pre(wb_now);
%         
%         M_R_mean_wb_accel_seq_pre(seq_now,wb_pre) = M_R_mean_wb_accel_now_pre(wb_now);
%         M_R_mean_wb_damp_seq_pre(seq_now,wb_pre) = M_R_mean_wb_damp_now_pre(wb_now);
%         M_R_mean_wb_seq_pre(seq_now,wb_pre) = M_R_mean_wb_now_pre(wb_now);
 
        %% clusters
        if sum(seq_now == n1) > 0
            % wb kin
            f_wb_seq_pre_c1(seq_now,wb_pre) = f_wb_now_pre(wb_now);

            stroke_wb_L_seq_bins_pre_c1(:,seq_now,wb_pre) = stroke_wb_L_bins_now_pre(:,wb_now);
            stroke_wb_R_seq_bins_pre_c1(:,seq_now,wb_pre) = stroke_wb_R_bins_now_pre(:,wb_now);

            dev_wb_L_seq_bins_pre_c1(:,seq_now,wb_pre) = dev_wb_L_bins_now_pre(:,wb_now);
            dev_wb_R_seq_bins_pre_c1(:,seq_now,wb_pre) = dev_wb_R_bins_now_pre(:,wb_now);

            pitch_wb_L_seq_bins_pre_c1(:,seq_now,wb_pre) = pitch_wb_L_bins_now_pre(:,wb_now);
            pitch_wb_R_seq_bins_pre_c1(:,seq_now,wb_pre) = pitch_wb_R_bins_now_pre(:,wb_now);

            % body kin
            roll_mean_wb_seq_pre_c1(seq_now,wb_pre) = roll_mean_wb_now_pre(wb_now);
            roll_dot_mean_wb_seq_pre_c1(seq_now,wb_pre) = roll_dot_mean_wb_now_pre(wb_now);
            roll_dot_dot_mean_wb_seq_pre_c1(seq_now,wb_pre) = roll_dot_dot_mean_wb_now_pre(wb_now);

            pitch_mean_wb_seq_pre_c1(seq_now,wb_pre) = pitch_mean_wb_now_pre(wb_now);
            pitch_dot_mean_wb_seq_pre_c1(seq_now,wb_pre) = pitch_dot_mean_wb_now_pre(wb_now);
            pitch_dot_dot_mean_wb_seq_pre_c1(seq_now,wb_pre) = pitch_dot_dot_mean_wb_now_pre(wb_now);

            yaw_mean_wb_seq_pre_c1(seq_now,wb_pre) = yaw_mean_wb_now_pre(wb_now);
            yaw_dot_mean_wb_seq_pre_c1(seq_now,wb_pre) = yaw_dot_mean_wb_now_pre(wb_now);
            yaw_dot_dot_mean_wb_seq_pre_c1(seq_now,wb_pre) = yaw_dot_dot_mean_wb_now_pre(wb_now);

%             drot_L_mean_wb_seq_pre_c1(seq_now,wb_pre) = drot_L_mean_wb_now_pre(wb_now);
%             rot_dot_L_mean_wb_seq_pre_c1(seq_now,wb_pre) = rot_dot_L_mean_wb_now_pre(wb_now);
%             rot_dot_dot_L_mean_wb_seq_pre_c1(seq_now,wb_pre) = rot_dot_dot_L_mean_wb_now_pre(wb_now);
% 
%             drot_R_mean_wb_seq_pre_c1(seq_now,wb_pre) = drot_R_mean_wb_now_pre(wb_now);
%             rot_dot_R_mean_wb_seq_pre_c1(seq_now,wb_pre) = rot_dot_R_mean_wb_now_pre(wb_now);
%             rot_dot_dot_R_mean_wb_seq_pre_c1(seq_now,wb_pre) = rot_dot_dot_R_mean_wb_now_pre(wb_now);
% 
%             % torque
%             Mroll_mean_wb_accel_seq_pre_c1(seq_now,wb_pre) = Mroll_mean_wb_accel_now_pre(wb_now);
%             Mroll_mean_wb_damp_seq_pre_c1(seq_now,wb_pre) = Mroll_mean_wb_damp_now_pre(wb_now);
%             Mroll_mean_wb_seq_pre_c1(seq_now,wb_pre) = Mroll_mean_wb_now_pre(wb_now);
% 
%             Mpitch_mean_wb_accel_seq_pre_c1(seq_now,wb_pre) = Mpitch_mean_wb_accel_now_pre(wb_now);
%             Mpitch_mean_wb_damp_seq_pre_c1(seq_now,wb_pre) = Mpitch_mean_wb_damp_now_pre(wb_now);
%             Mpitch_mean_wb_seq_pre_c1(seq_now,wb_pre) = Mpitch_mean_wb_now_pre(wb_now);
% 
%             Myaw_mean_wb_accel_seq_pre_c1(seq_now,wb_pre) = Myaw_mean_wb_accel_now_pre(wb_now);
%             Myaw_mean_wb_damp_seq_pre_c1(seq_now,wb_pre) = Myaw_mean_wb_damp_now_pre(wb_now);
%             Myaw_mean_wb_seq_pre_c1(seq_now,wb_pre) = Myaw_mean_wb_now_pre(wb_now);
% 
%             M_L_mean_wb_accel_seq_pre_c1(seq_now,wb_pre) = M_L_mean_wb_accel_now_pre(wb_now);
%             M_L_mean_wb_damp_seq_pre_c1(seq_now,wb_pre) = M_L_mean_wb_damp_now_pre(wb_now);
%             M_L_mean_wb_seq_pre_c1(seq_now,wb_pre) = M_L_mean_wb_now_pre(wb_now);
% 
%             M_R_mean_wb_accel_seq_pre_c1(seq_now,wb_pre) = M_R_mean_wb_accel_now_pre(wb_now);
%             M_R_mean_wb_damp_seq_pre_c1(seq_now,wb_pre) = M_R_mean_wb_damp_now_pre(wb_now);
%             M_R_mean_wb_seq_pre_c1(seq_now,wb_pre) = M_R_mean_wb_now_pre(wb_now);
            
        elseif sum(seq_now == n2) > 0
            % wb kin
            f_wb_seq_pre_c2(seq_now,wb_pre) = f_wb_now_pre(wb_now);

            stroke_wb_L_seq_bins_pre_c2(:,seq_now,wb_pre) = stroke_wb_L_bins_now_pre(:,wb_now);
            stroke_wb_R_seq_bins_pre_c2(:,seq_now,wb_pre) = stroke_wb_R_bins_now_pre(:,wb_now);

            dev_wb_L_seq_bins_pre_c2(:,seq_now,wb_pre) = dev_wb_L_bins_now_pre(:,wb_now);
            dev_wb_R_seq_bins_pre_c2(:,seq_now,wb_pre) = dev_wb_R_bins_now_pre(:,wb_now);

            pitch_wb_L_seq_bins_pre_c2(:,seq_now,wb_pre) = pitch_wb_L_bins_now_pre(:,wb_now);
            pitch_wb_R_seq_bins_pre_c2(:,seq_now,wb_pre) = pitch_wb_R_bins_now_pre(:,wb_now);

            % body kin
            roll_mean_wb_seq_pre_c2(seq_now,wb_pre) = roll_mean_wb_now_pre(wb_now);
            roll_dot_mean_wb_seq_pre_c2(seq_now,wb_pre) = roll_dot_mean_wb_now_pre(wb_now);
            roll_dot_dot_mean_wb_seq_pre_c2(seq_now,wb_pre) = roll_dot_dot_mean_wb_now_pre(wb_now);

            pitch_mean_wb_seq_pre_c2(seq_now,wb_pre) = pitch_mean_wb_now_pre(wb_now);
            pitch_dot_mean_wb_seq_pre_c2(seq_now,wb_pre) = pitch_dot_mean_wb_now_pre(wb_now);
            pitch_dot_dot_mean_wb_seq_pre_c2(seq_now,wb_pre) = pitch_dot_dot_mean_wb_now_pre(wb_now);

            yaw_mean_wb_seq_pre_c2(seq_now,wb_pre) = yaw_mean_wb_now_pre(wb_now);
            yaw_dot_mean_wb_seq_pre_c2(seq_now,wb_pre) = yaw_dot_mean_wb_now_pre(wb_now);
            yaw_dot_dot_mean_wb_seq_pre_c2(seq_now,wb_pre) = yaw_dot_dot_mean_wb_now_pre(wb_now);

%             drot_L_mean_wb_seq_pre_c2(seq_now,wb_pre) = drot_L_mean_wb_now_pre(wb_now);
%             rot_dot_L_mean_wb_seq_pre_c2(seq_now,wb_pre) = rot_dot_L_mean_wb_now_pre(wb_now);
%             rot_dot_dot_L_mean_wb_seq_pre_c2(seq_now,wb_pre) = rot_dot_dot_L_mean_wb_now_pre(wb_now);
% 
%             drot_R_mean_wb_seq_pre_c2(seq_now,wb_pre) = drot_R_mean_wb_now_pre(wb_now);
%             rot_dot_R_mean_wb_seq_pre_c2(seq_now,wb_pre) = rot_dot_R_mean_wb_now_pre(wb_now);
%             rot_dot_dot_R_mean_wb_seq_pre_c2(seq_now,wb_pre) = rot_dot_dot_R_mean_wb_now_pre(wb_now);
% 
%             % torque
%             Mroll_mean_wb_accel_seq_pre_c2(seq_now,wb_pre) = Mroll_mean_wb_accel_now_pre(wb_now);
%             Mroll_mean_wb_damp_seq_pre_c2(seq_now,wb_pre) = Mroll_mean_wb_damp_now_pre(wb_now);
%             Mroll_mean_wb_seq_pre_c2(seq_now,wb_pre) = Mroll_mean_wb_now_pre(wb_now);
% 
%             Mpitch_mean_wb_accel_seq_pre_c2(seq_now,wb_pre) = Mpitch_mean_wb_accel_now_pre(wb_now);
%             Mpitch_mean_wb_damp_seq_pre_c2(seq_now,wb_pre) = Mpitch_mean_wb_damp_now_pre(wb_now);
%             Mpitch_mean_wb_seq_pre_c2(seq_now,wb_pre) = Mpitch_mean_wb_now_pre(wb_now);
% 
%             Myaw_mean_wb_accel_seq_pre_c2(seq_now,wb_pre) = Myaw_mean_wb_accel_now_pre(wb_now);
%             Myaw_mean_wb_damp_seq_pre_c2(seq_now,wb_pre) = Myaw_mean_wb_damp_now_pre(wb_now);
%             Myaw_mean_wb_seq_pre_c2(seq_now,wb_pre) = Myaw_mean_wb_now_pre(wb_now);
% 
%             M_L_mean_wb_accel_seq_pre_c2(seq_now,wb_pre) = M_L_mean_wb_accel_now_pre(wb_now);
%             M_L_mean_wb_damp_seq_pre_c2(seq_now,wb_pre) = M_L_mean_wb_damp_now_pre(wb_now);
%             M_L_mean_wb_seq_pre_c2(seq_now,wb_pre) = M_L_mean_wb_now_pre(wb_now);
% 
%             M_R_mean_wb_accel_seq_pre_c2(seq_now,wb_pre) = M_R_mean_wb_accel_now_pre(wb_now);
%             M_R_mean_wb_damp_seq_pre_c2(seq_now,wb_pre) = M_R_mean_wb_damp_now_pre(wb_now);
%             M_R_mean_wb_seq_pre_c2(seq_now,wb_pre) = M_R_mean_wb_now_pre(wb_now);
            
        elseif sum(seq_now == n3) > 0
            % wb kin
            f_wb_seq_pre_c3(seq_now,wb_pre) = f_wb_now_pre(wb_now);

            stroke_wb_L_seq_bins_pre_c3(:,seq_now,wb_pre) = stroke_wb_L_bins_now_pre(:,wb_now);
            stroke_wb_R_seq_bins_pre_c3(:,seq_now,wb_pre) = stroke_wb_R_bins_now_pre(:,wb_now);

            dev_wb_L_seq_bins_pre_c3(:,seq_now,wb_pre) = dev_wb_L_bins_now_pre(:,wb_now);
            dev_wb_R_seq_bins_pre_c3(:,seq_now,wb_pre) = dev_wb_R_bins_now_pre(:,wb_now);

            pitch_wb_L_seq_bins_pre_c3(:,seq_now,wb_pre) = pitch_wb_L_bins_now_pre(:,wb_now);
            pitch_wb_R_seq_bins_pre_c3(:,seq_now,wb_pre) = pitch_wb_R_bins_now_pre(:,wb_now);

            % body kin
            roll_mean_wb_seq_pre_c3(seq_now,wb_pre) = roll_mean_wb_now_pre(wb_now);
            roll_dot_mean_wb_seq_pre_c3(seq_now,wb_pre) = roll_dot_mean_wb_now_pre(wb_now);
            roll_dot_dot_mean_wb_seq_pre_c3(seq_now,wb_pre) = roll_dot_dot_mean_wb_now_pre(wb_now);

            pitch_mean_wb_seq_pre_c3(seq_now,wb_pre) = pitch_mean_wb_now_pre(wb_now);
            pitch_dot_mean_wb_seq_pre_c3(seq_now,wb_pre) = pitch_dot_mean_wb_now_pre(wb_now);
            pitch_dot_dot_mean_wb_seq_pre_c3(seq_now,wb_pre) = pitch_dot_dot_mean_wb_now_pre(wb_now);

            yaw_mean_wb_seq_pre_c3(seq_now,wb_pre) = yaw_mean_wb_now_pre(wb_now);
            yaw_dot_mean_wb_seq_pre_c3(seq_now,wb_pre) = yaw_dot_mean_wb_now_pre(wb_now);
            yaw_dot_dot_mean_wb_seq_pre_c3(seq_now,wb_pre) = yaw_dot_dot_mean_wb_now_pre(wb_now);

%             drot_L_mean_wb_seq_pre_c3(seq_now,wb_pre) = drot_L_mean_wb_now_pre(wb_now);
%             rot_dot_L_mean_wb_seq_pre_c3(seq_now,wb_pre) = rot_dot_L_mean_wb_now_pre(wb_now);
%             rot_dot_dot_L_mean_wb_seq_pre_c3(seq_now,wb_pre) = rot_dot_dot_L_mean_wb_now_pre(wb_now);
% 
%             drot_R_mean_wb_seq_pre_c3(seq_now,wb_pre) = drot_R_mean_wb_now_pre(wb_now);
%             rot_dot_R_mean_wb_seq_pre_c3(seq_now,wb_pre) = rot_dot_R_mean_wb_now_pre(wb_now);
%             rot_dot_dot_R_mean_wb_seq_pre_c3(seq_now,wb_pre) = rot_dot_dot_R_mean_wb_now_pre(wb_now);
% 
%             % torque
%             Mroll_mean_wb_accel_seq_pre_c3(seq_now,wb_pre) = Mroll_mean_wb_accel_now_pre(wb_now);
%             Mroll_mean_wb_damp_seq_pre_c3(seq_now,wb_pre) = Mroll_mean_wb_damp_now_pre(wb_now);
%             Mroll_mean_wb_seq_pre_c3(seq_now,wb_pre) = Mroll_mean_wb_now_pre(wb_now);
% 
%             Mpitch_mean_wb_accel_seq_pre_c3(seq_now,wb_pre) = Mpitch_mean_wb_accel_now_pre(wb_now);
%             Mpitch_mean_wb_damp_seq_pre_c3(seq_now,wb_pre) = Mpitch_mean_wb_damp_now_pre(wb_now);
%             Mpitch_mean_wb_seq_pre_c3(seq_now,wb_pre) = Mpitch_mean_wb_now_pre(wb_now);
% 
%             Myaw_mean_wb_accel_seq_pre_c3(seq_now,wb_pre) = Myaw_mean_wb_accel_now_pre(wb_now);
%             Myaw_mean_wb_damp_seq_pre_c3(seq_now,wb_pre) = Myaw_mean_wb_damp_now_pre(wb_now);
%             Myaw_mean_wb_seq_pre_c3(seq_now,wb_pre) = Myaw_mean_wb_now_pre(wb_now);
% 
%             M_L_mean_wb_accel_seq_pre_c3(seq_now,wb_pre) = M_L_mean_wb_accel_now_pre(wb_now);
%             M_L_mean_wb_damp_seq_pre_c3(seq_now,wb_pre) = M_L_mean_wb_damp_now_pre(wb_now);
%             M_L_mean_wb_seq_pre_c3(seq_now,wb_pre) = M_L_mean_wb_now_pre(wb_now);
% 
%             M_R_mean_wb_accel_seq_pre_c3(seq_now,wb_pre) = M_R_mean_wb_accel_now_pre(wb_now);
%             M_R_mean_wb_damp_seq_pre_c3(seq_now,wb_pre) = M_R_mean_wb_damp_now_pre(wb_now);
%             M_R_mean_wb_seq_pre_c3(seq_now,wb_pre) = M_R_mean_wb_now_pre(wb_now);
            
        elseif sum(seq_now == n4) > 0
            % wb kin
            f_wb_seq_pre_c4(seq_now,wb_pre) = f_wb_now_pre(wb_now);

            stroke_wb_L_seq_bins_pre_c4(:,seq_now,wb_pre) = stroke_wb_L_bins_now_pre(:,wb_now);
            stroke_wb_R_seq_bins_pre_c4(:,seq_now,wb_pre) = stroke_wb_R_bins_now_pre(:,wb_now);

            dev_wb_L_seq_bins_pre_c4(:,seq_now,wb_pre) = dev_wb_L_bins_now_pre(:,wb_now);
            dev_wb_R_seq_bins_pre_c4(:,seq_now,wb_pre) = dev_wb_R_bins_now_pre(:,wb_now);

            pitch_wb_L_seq_bins_pre_c4(:,seq_now,wb_pre) = pitch_wb_L_bins_now_pre(:,wb_now);
            pitch_wb_R_seq_bins_pre_c4(:,seq_now,wb_pre) = pitch_wb_R_bins_now_pre(:,wb_now);

            % body kin
            roll_mean_wb_seq_pre_c4(seq_now,wb_pre) = roll_mean_wb_now_pre(wb_now);
            roll_dot_mean_wb_seq_pre_c4(seq_now,wb_pre) = roll_dot_mean_wb_now_pre(wb_now);
            roll_dot_dot_mean_wb_seq_pre_c4(seq_now,wb_pre) = roll_dot_dot_mean_wb_now_pre(wb_now);

            pitch_mean_wb_seq_pre_c4(seq_now,wb_pre) = pitch_mean_wb_now_pre(wb_now);
            pitch_dot_mean_wb_seq_pre_c4(seq_now,wb_pre) = pitch_dot_mean_wb_now_pre(wb_now);
            pitch_dot_dot_mean_wb_seq_pre_c4(seq_now,wb_pre) = pitch_dot_dot_mean_wb_now_pre(wb_now);

            yaw_mean_wb_seq_pre_c4(seq_now,wb_pre) = yaw_mean_wb_now_pre(wb_now);
            yaw_dot_mean_wb_seq_pre_c4(seq_now,wb_pre) = yaw_dot_mean_wb_now_pre(wb_now);
            yaw_dot_dot_mean_wb_seq_pre_c4(seq_now,wb_pre) = yaw_dot_dot_mean_wb_now_pre(wb_now);

%             drot_L_mean_wb_seq_pre_c4(seq_now,wb_pre) = drot_L_mean_wb_now_pre(wb_now);
%             rot_dot_L_mean_wb_seq_pre_c4(seq_now,wb_pre) = rot_dot_L_mean_wb_now_pre(wb_now);
%             rot_dot_dot_L_mean_wb_seq_pre_c4(seq_now,wb_pre) = rot_dot_dot_L_mean_wb_now_pre(wb_now);
% 
%             drot_R_mean_wb_seq_pre_c4(seq_now,wb_pre) = drot_R_mean_wb_now_pre(wb_now);
%             rot_dot_R_mean_wb_seq_pre_c4(seq_now,wb_pre) = rot_dot_R_mean_wb_now_pre(wb_now);
%             rot_dot_dot_R_mean_wb_seq_pre_c4(seq_now,wb_pre) = rot_dot_dot_R_mean_wb_now_pre(wb_now);
% 
%             % torque
%             Mroll_mean_wb_accel_seq_pre_c4(seq_now,wb_pre) = Mroll_mean_wb_accel_now_pre(wb_now);
%             Mroll_mean_wb_damp_seq_pre_c4(seq_now,wb_pre) = Mroll_mean_wb_damp_now_pre(wb_now);
%             Mroll_mean_wb_seq_pre_c4(seq_now,wb_pre) = Mroll_mean_wb_now_pre(wb_now);
% 
%             Mpitch_mean_wb_accel_seq_pre_c4(seq_now,wb_pre) = Mpitch_mean_wb_accel_now_pre(wb_now);
%             Mpitch_mean_wb_damp_seq_pre_c4(seq_now,wb_pre) = Mpitch_mean_wb_damp_now_pre(wb_now);
%             Mpitch_mean_wb_seq_pre_c4(seq_now,wb_pre) = Mpitch_mean_wb_now_pre(wb_now);
% 
%             Myaw_mean_wb_accel_seq_pre_c4(seq_now,wb_pre) = Myaw_mean_wb_accel_now_pre(wb_now);
%             Myaw_mean_wb_damp_seq_pre_c4(seq_now,wb_pre) = Myaw_mean_wb_damp_now_pre(wb_now);
%             Myaw_mean_wb_seq_pre_c4(seq_now,wb_pre) = Myaw_mean_wb_now_pre(wb_now);
% 
%             M_L_mean_wb_accel_seq_pre_c4(seq_now,wb_pre) = M_L_mean_wb_accel_now_pre(wb_now);
%             M_L_mean_wb_damp_seq_pre_c4(seq_now,wb_pre) = M_L_mean_wb_damp_now_pre(wb_now);
%             M_L_mean_wb_seq_pre_c4(seq_now,wb_pre) = M_L_mean_wb_now_pre(wb_now);
% 
%             M_R_mean_wb_accel_seq_pre_c4(seq_now,wb_pre) = M_R_mean_wb_accel_now_pre(wb_now);
%             M_R_mean_wb_damp_seq_pre_c4(seq_now,wb_pre) = M_R_mean_wb_damp_now_pre(wb_now);
%             M_R_mean_wb_seq_pre_c4(seq_now,wb_pre) = M_R_mean_wb_now_pre(wb_now);
        end
        end
      end
    
    %% post response wbs
    wb_post_max = length(frame_nr_2resp_now_post);
    for wb_post = 1:wb_post_max
        
        if wb_post <= wb_max
        
            % wb kin
            f_wb_seq_post(seq_now,wb_post) = f_wb_now_post(wb_post);

            stroke_wb_L_seq_bins_post(:,seq_now,wb_post) = stroke_wb_L_bins_now_post(:,wb_post);
            stroke_wb_R_seq_bins_post(:,seq_now,wb_post) = stroke_wb_R_bins_now_post(:,wb_post);

            dev_wb_L_seq_bins_post(:,seq_now,wb_post) = dev_wb_L_bins_now_post(:,wb_post);
            dev_wb_R_seq_bins_post(:,seq_now,wb_post) = dev_wb_R_bins_now_post(:,wb_post);

            pitch_wb_L_seq_bins_post(:,seq_now,wb_post) = pitch_wb_L_bins_now_post(:,wb_post);
            pitch_wb_R_seq_bins_post(:,seq_now,wb_post) = pitch_wb_R_bins_now_post(:,wb_post);

            % body kin
            roll_mean_wb_seq_post(seq_now,wb_post) = roll_mean_wb_now_post(wb_post);
            roll_dot_mean_wb_seq_post(seq_now,wb_post) = roll_dot_mean_wb_now_post(wb_post);
            roll_dot_dot_mean_wb_seq_post(seq_now,wb_post) = roll_dot_dot_mean_wb_now_post(wb_post);

            pitch_mean_wb_seq_post(seq_now,wb_post) = pitch_mean_wb_now_post(wb_post);
            pitch_dot_mean_wb_seq_post(seq_now,wb_post) = pitch_dot_mean_wb_now_post(wb_post);
            pitch_dot_dot_mean_wb_seq_post(seq_now,wb_post) = pitch_dot_dot_mean_wb_now_post(wb_post);

            yaw_mean_wb_seq_post(seq_now,wb_post) = yaw_mean_wb_now_post(wb_post);
            yaw_dot_mean_wb_seq_post(seq_now,wb_post) = yaw_dot_mean_wb_now_post(wb_post);
            yaw_dot_dot_mean_wb_seq_post(seq_now,wb_post) = yaw_dot_dot_mean_wb_now_post(wb_post);

%             drot_L_mean_wb_seq_post(seq_now,wb_post) = drot_L_mean_wb_now_post(wb_post);
%             rot_dot_L_mean_wb_seq_post(seq_now,wb_post) = rot_dot_L_mean_wb_now_post(wb_post);
%             rot_dot_dot_L_mean_wb_seq_post(seq_now,wb_post) = rot_dot_dot_L_mean_wb_now_post(wb_post);
% 
%             drot_R_mean_wb_seq_post(seq_now,wb_post) = drot_R_mean_wb_now_post(wb_post);
%             rot_dot_R_mean_wb_seq_post(seq_now,wb_post) = rot_dot_R_mean_wb_now_post(wb_post);
%             rot_dot_dot_R_mean_wb_seq_post(seq_now,wb_post) = rot_dot_dot_R_mean_wb_now_post(wb_post);
% 
%             % torque
%             Mroll_mean_wb_accel_seq_post(seq_now,wb_post) = Mroll_mean_wb_accel_now_post(wb_post);
%             Mroll_mean_wb_damp_seq_post(seq_now,wb_post) = Mroll_mean_wb_damp_now_post(wb_post);
%             Mroll_mean_wb_seq_post(seq_now,wb_post) = Mroll_mean_wb_now_post(wb_post);
% 
%             Mpitch_mean_wb_accel_seq_post(seq_now,wb_post) = Mpitch_mean_wb_accel_now_post(wb_post);
%             Mpitch_mean_wb_damp_seq_post(seq_now,wb_post) = Mpitch_mean_wb_damp_now_post(wb_post);
%             Mpitch_mean_wb_seq_post(seq_now,wb_post) = Mpitch_mean_wb_now_post(wb_post);
% 
%             Myaw_mean_wb_accel_seq_post(seq_now,wb_post) = Myaw_mean_wb_accel_now_post(wb_post);
%             Myaw_mean_wb_damp_seq_post(seq_now,wb_post) = Myaw_mean_wb_damp_now_post(wb_post);
%             Myaw_mean_wb_seq_post(seq_now,wb_post) = Myaw_mean_wb_now_post(wb_post);
% 
%             M_L_mean_wb_accel_seq_post(seq_now,wb_post) = M_L_mean_wb_accel_now_post(wb_post);
%             M_L_mean_wb_damp_seq_post(seq_now,wb_post) = M_L_mean_wb_damp_now_post(wb_post);
%             M_L_mean_wb_seq_post(seq_now,wb_post) = M_L_mean_wb_now_post(wb_post);
% 
%             M_R_mean_wb_accel_seq_post(seq_now,wb_post) = M_R_mean_wb_accel_now_post(wb_post);
%             M_R_mean_wb_damp_seq_post(seq_now,wb_post) = M_R_mean_wb_damp_now_post(wb_post);
%             M_R_mean_wb_seq_post(seq_now,wb_post) = M_R_mean_wb_now_post(wb_post);

            %% clusters
            if sum(seq_now == n1) > 0
                % wb kin
                f_wb_seq_post_c1(seq_now,wb_post) = f_wb_now_post(wb_post);

                stroke_wb_L_seq_bins_post_c1(:,seq_now,wb_post) = stroke_wb_L_bins_now_post(:,wb_post);
                stroke_wb_R_seq_bins_post_c1(:,seq_now,wb_post) = stroke_wb_R_bins_now_post(:,wb_post);

                dev_wb_L_seq_bins_post_c1(:,seq_now,wb_post) = dev_wb_L_bins_now_post(:,wb_post);
                dev_wb_R_seq_bins_post_c1(:,seq_now,wb_post) = dev_wb_R_bins_now_post(:,wb_post);

                pitch_wb_L_seq_bins_post_c1(:,seq_now,wb_post) = pitch_wb_L_bins_now_post(:,wb_post);
                pitch_wb_R_seq_bins_post_c1(:,seq_now,wb_post) = pitch_wb_R_bins_now_post(:,wb_post);

                % body kin
                roll_mean_wb_seq_post_c1(seq_now,wb_post) = roll_mean_wb_now_post(wb_post);
                roll_dot_mean_wb_seq_post_c1(seq_now,wb_post) = roll_dot_mean_wb_now_post(wb_post);
                roll_dot_dot_mean_wb_seq_post_c1(seq_now,wb_post) = roll_dot_dot_mean_wb_now_post(wb_post);

                pitch_mean_wb_seq_post_c1(seq_now,wb_post) = pitch_mean_wb_now_post(wb_post);
                pitch_dot_mean_wb_seq_post_c1(seq_now,wb_post) = pitch_dot_mean_wb_now_post(wb_post);
                pitch_dot_dot_mean_wb_seq_post_c1(seq_now,wb_post) = pitch_dot_dot_mean_wb_now_post(wb_post);

                yaw_mean_wb_seq_post_c1(seq_now,wb_post) = yaw_mean_wb_now_post(wb_post);
                yaw_dot_mean_wb_seq_post_c1(seq_now,wb_post) = yaw_dot_mean_wb_now_post(wb_post);
                yaw_dot_dot_mean_wb_seq_post_c1(seq_now,wb_post) = yaw_dot_dot_mean_wb_now_post(wb_post);

%                 drot_L_mean_wb_seq_post_c1(seq_now,wb_post) = drot_L_mean_wb_now_post(wb_post);
%                 rot_dot_L_mean_wb_seq_post_c1(seq_now,wb_post) = rot_dot_L_mean_wb_now_post(wb_post);
%                 rot_dot_dot_L_mean_wb_seq_post_c1(seq_now,wb_post) = rot_dot_dot_L_mean_wb_now_post(wb_post);
% 
%                 drot_R_mean_wb_seq_post_c1(seq_now,wb_post) = drot_R_mean_wb_now_post(wb_post);
%                 rot_dot_R_mean_wb_seq_post_c1(seq_now,wb_post) = rot_dot_R_mean_wb_now_post(wb_post);
%                 rot_dot_dot_R_mean_wb_seq_post_c1(seq_now,wb_post) = rot_dot_dot_R_mean_wb_now_post(wb_post);
% 
%                 % torque
%                 Mroll_mean_wb_accel_seq_post_c1(seq_now,wb_post) = Mroll_mean_wb_accel_now_post(wb_post);
%                 Mroll_mean_wb_damp_seq_post_c1(seq_now,wb_post) = Mroll_mean_wb_damp_now_post(wb_post);
%                 Mroll_mean_wb_seq_post_c1(seq_now,wb_post) = Mroll_mean_wb_now_post(wb_post);
% 
%                 Mpitch_mean_wb_accel_seq_post_c1(seq_now,wb_post) = Mpitch_mean_wb_accel_now_post(wb_post);
%                 Mpitch_mean_wb_damp_seq_post_c1(seq_now,wb_post) = Mpitch_mean_wb_damp_now_post(wb_post);
%                 Mpitch_mean_wb_seq_post_c1(seq_now,wb_post) = Mpitch_mean_wb_now_post(wb_post);
% 
%                 Myaw_mean_wb_accel_seq_post_c1(seq_now,wb_post) = Myaw_mean_wb_accel_now_post(wb_post);
%                 Myaw_mean_wb_damp_seq_post_c1(seq_now,wb_post) = Myaw_mean_wb_damp_now_post(wb_post);
%                 Myaw_mean_wb_seq_post_c1(seq_now,wb_post) = Myaw_mean_wb_now_post(wb_post);
% 
%                 M_L_mean_wb_accel_seq_post_c1(seq_now,wb_post) = M_L_mean_wb_accel_now_post(wb_post);
%                 M_L_mean_wb_damp_seq_post_c1(seq_now,wb_post) = M_L_mean_wb_damp_now_post(wb_post);
%                 M_L_mean_wb_seq_post_c1(seq_now,wb_post) = M_L_mean_wb_now_post(wb_post);
% 
%                 M_R_mean_wb_accel_seq_post_c1(seq_now,wb_post) = M_R_mean_wb_accel_now_post(wb_post);
%                 M_R_mean_wb_damp_seq_post_c1(seq_now,wb_post) = M_R_mean_wb_damp_now_post(wb_post);
%                 M_R_mean_wb_seq_post_c1(seq_now,wb_post) = M_R_mean_wb_now_post(wb_post);

            elseif sum(seq_now == n2) > 0
                % wb kin
                f_wb_seq_post_c2(seq_now,wb_post) = f_wb_now_post(wb_post);

                stroke_wb_L_seq_bins_post_c2(:,seq_now,wb_post) = stroke_wb_L_bins_now_post(:,wb_post);
                stroke_wb_R_seq_bins_post_c2(:,seq_now,wb_post) = stroke_wb_R_bins_now_post(:,wb_post);

                dev_wb_L_seq_bins_post_c2(:,seq_now,wb_post) = dev_wb_L_bins_now_post(:,wb_post);
                dev_wb_R_seq_bins_post_c2(:,seq_now,wb_post) = dev_wb_R_bins_now_post(:,wb_post);

                pitch_wb_L_seq_bins_post_c2(:,seq_now,wb_post) = pitch_wb_L_bins_now_post(:,wb_post);
                pitch_wb_R_seq_bins_post_c2(:,seq_now,wb_post) = pitch_wb_R_bins_now_post(:,wb_post);

                % body kin
                roll_mean_wb_seq_post_c2(seq_now,wb_post) = roll_mean_wb_now_post(wb_post);
                roll_dot_mean_wb_seq_post_c2(seq_now,wb_post) = roll_dot_mean_wb_now_post(wb_post);
                roll_dot_dot_mean_wb_seq_post_c2(seq_now,wb_post) = roll_dot_dot_mean_wb_now_post(wb_post);

                pitch_mean_wb_seq_post_c2(seq_now,wb_post) = pitch_mean_wb_now_post(wb_post);
                pitch_dot_mean_wb_seq_post_c2(seq_now,wb_post) = pitch_dot_mean_wb_now_post(wb_post);
                pitch_dot_dot_mean_wb_seq_post_c2(seq_now,wb_post) = pitch_dot_dot_mean_wb_now_post(wb_post);

                yaw_mean_wb_seq_post_c2(seq_now,wb_post) = yaw_mean_wb_now_post(wb_post);
                yaw_dot_mean_wb_seq_post_c2(seq_now,wb_post) = yaw_dot_mean_wb_now_post(wb_post);
                yaw_dot_dot_mean_wb_seq_post_c2(seq_now,wb_post) = yaw_dot_dot_mean_wb_now_post(wb_post);

%                 drot_L_mean_wb_seq_post_c2(seq_now,wb_post) = drot_L_mean_wb_now_post(wb_post);
%                 rot_dot_L_mean_wb_seq_post_c2(seq_now,wb_post) = rot_dot_L_mean_wb_now_post(wb_post);
%                 rot_dot_dot_L_mean_wb_seq_post_c2(seq_now,wb_post) = rot_dot_dot_L_mean_wb_now_post(wb_post);
% 
%                 drot_R_mean_wb_seq_post_c2(seq_now,wb_post) = drot_R_mean_wb_now_post(wb_post);
%                 rot_dot_R_mean_wb_seq_post_c2(seq_now,wb_post) = rot_dot_R_mean_wb_now_post(wb_post);
%                 rot_dot_dot_R_mean_wb_seq_post_c2(seq_now,wb_post) = rot_dot_dot_R_mean_wb_now_post(wb_post);
% 
%                 % torque
%                 Mroll_mean_wb_accel_seq_post_c2(seq_now,wb_post) = Mroll_mean_wb_accel_now_post(wb_post);
%                 Mroll_mean_wb_damp_seq_post_c2(seq_now,wb_post) = Mroll_mean_wb_damp_now_post(wb_post);
%                 Mroll_mean_wb_seq_post_c2(seq_now,wb_post) = Mroll_mean_wb_now_post(wb_post);
% 
%                 Mpitch_mean_wb_accel_seq_post_c2(seq_now,wb_post) = Mpitch_mean_wb_accel_now_post(wb_post);
%                 Mpitch_mean_wb_damp_seq_post_c2(seq_now,wb_post) = Mpitch_mean_wb_damp_now_post(wb_post);
%                 Mpitch_mean_wb_seq_post_c2(seq_now,wb_post) = Mpitch_mean_wb_now_post(wb_post);
% 
%                 Myaw_mean_wb_accel_seq_post_c2(seq_now,wb_post) = Myaw_mean_wb_accel_now_post(wb_post);
%                 Myaw_mean_wb_damp_seq_post_c2(seq_now,wb_post) = Myaw_mean_wb_damp_now_post(wb_post);
%                 Myaw_mean_wb_seq_post_c2(seq_now,wb_post) = Myaw_mean_wb_now_post(wb_post);
% 
%                 M_L_mean_wb_accel_seq_post_c2(seq_now,wb_post) = M_L_mean_wb_accel_now_post(wb_post);
%                 M_L_mean_wb_damp_seq_post_c2(seq_now,wb_post) = M_L_mean_wb_damp_now_post(wb_post);
%                 M_L_mean_wb_seq_post_c2(seq_now,wb_post) = M_L_mean_wb_now_post(wb_post);
% 
%                 M_R_mean_wb_accel_seq_post_c2(seq_now,wb_post) = M_R_mean_wb_accel_now_post(wb_post);
%                 M_R_mean_wb_damp_seq_post_c2(seq_now,wb_post) = M_R_mean_wb_damp_now_post(wb_post);
%                 M_R_mean_wb_seq_post_c2(seq_now,wb_post) = M_R_mean_wb_now_post(wb_post);

            elseif sum(seq_now == n3) > 0
                % wb kin
                f_wb_seq_post_c3(seq_now,wb_post) = f_wb_now_post(wb_post);

                stroke_wb_L_seq_bins_post_c3(:,seq_now,wb_post) = stroke_wb_L_bins_now_post(:,wb_post);
                stroke_wb_R_seq_bins_post_c3(:,seq_now,wb_post) = stroke_wb_R_bins_now_post(:,wb_post);

                dev_wb_L_seq_bins_post_c3(:,seq_now,wb_post) = dev_wb_L_bins_now_post(:,wb_post);
                dev_wb_R_seq_bins_post_c3(:,seq_now,wb_post) = dev_wb_R_bins_now_post(:,wb_post);

                pitch_wb_L_seq_bins_post_c3(:,seq_now,wb_post) = pitch_wb_L_bins_now_post(:,wb_post);
                pitch_wb_R_seq_bins_post_c3(:,seq_now,wb_post) = pitch_wb_R_bins_now_post(:,wb_post);

                % body kin
                roll_mean_wb_seq_post_c3(seq_now,wb_post) = roll_mean_wb_now_post(wb_post);
                roll_dot_mean_wb_seq_post_c3(seq_now,wb_post) = roll_dot_mean_wb_now_post(wb_post);
                roll_dot_dot_mean_wb_seq_post_c3(seq_now,wb_post) = roll_dot_dot_mean_wb_now_post(wb_post);

                pitch_mean_wb_seq_post_c3(seq_now,wb_post) = pitch_mean_wb_now_post(wb_post);
                pitch_dot_mean_wb_seq_post_c3(seq_now,wb_post) = pitch_dot_mean_wb_now_post(wb_post);
                pitch_dot_dot_mean_wb_seq_post_c3(seq_now,wb_post) = pitch_dot_dot_mean_wb_now_post(wb_post);

                yaw_mean_wb_seq_post_c3(seq_now,wb_post) = yaw_mean_wb_now_post(wb_post);
                yaw_dot_mean_wb_seq_post_c3(seq_now,wb_post) = yaw_dot_mean_wb_now_post(wb_post);
                yaw_dot_dot_mean_wb_seq_post_c3(seq_now,wb_post) = yaw_dot_dot_mean_wb_now_post(wb_post);

%                 drot_L_mean_wb_seq_post_c3(seq_now,wb_post) = drot_L_mean_wb_now_post(wb_post);
%                 rot_dot_L_mean_wb_seq_post_c3(seq_now,wb_post) = rot_dot_L_mean_wb_now_post(wb_post);
%                 rot_dot_dot_L_mean_wb_seq_post_c3(seq_now,wb_post) = rot_dot_dot_L_mean_wb_now_post(wb_post);
% 
%                 drot_R_mean_wb_seq_post_c3(seq_now,wb_post) = drot_R_mean_wb_now_post(wb_post);
%                 rot_dot_R_mean_wb_seq_post_c3(seq_now,wb_post) = rot_dot_R_mean_wb_now_post(wb_post);
%                 rot_dot_dot_R_mean_wb_seq_post_c3(seq_now,wb_post) = rot_dot_dot_R_mean_wb_now_post(wb_post);
% 
%                 % torque
%                 Mroll_mean_wb_accel_seq_post_c3(seq_now,wb_post) = Mroll_mean_wb_accel_now_post(wb_post);
%                 Mroll_mean_wb_damp_seq_post_c3(seq_now,wb_post) = Mroll_mean_wb_damp_now_post(wb_post);
%                 Mroll_mean_wb_seq_post_c3(seq_now,wb_post) = Mroll_mean_wb_now_post(wb_post);
% 
%                 Mpitch_mean_wb_accel_seq_post_c3(seq_now,wb_post) = Mpitch_mean_wb_accel_now_post(wb_post);
%                 Mpitch_mean_wb_damp_seq_post_c3(seq_now,wb_post) = Mpitch_mean_wb_damp_now_post(wb_post);
%                 Mpitch_mean_wb_seq_post_c3(seq_now,wb_post) = Mpitch_mean_wb_now_post(wb_post);
% 
%                 Myaw_mean_wb_accel_seq_post_c3(seq_now,wb_post) = Myaw_mean_wb_accel_now_post(wb_post);
%                 Myaw_mean_wb_damp_seq_post_c3(seq_now,wb_post) = Myaw_mean_wb_damp_now_post(wb_post);
%                 Myaw_mean_wb_seq_post_c3(seq_now,wb_post) = Myaw_mean_wb_now_post(wb_post);
% 
%                 M_L_mean_wb_accel_seq_post_c3(seq_now,wb_post) = M_L_mean_wb_accel_now_post(wb_post);
%                 M_L_mean_wb_damp_seq_post_c3(seq_now,wb_post) = M_L_mean_wb_damp_now_post(wb_post);
%                 M_L_mean_wb_seq_post_c3(seq_now,wb_post) = M_L_mean_wb_now_post(wb_post);
% 
%                 M_R_mean_wb_accel_seq_post_c3(seq_now,wb_post) = M_R_mean_wb_accel_now_post(wb_post);
%                 M_R_mean_wb_damp_seq_post_c3(seq_now,wb_post) = M_R_mean_wb_damp_now_post(wb_post);
%                 M_R_mean_wb_seq_post_c3(seq_now,wb_post) = M_R_mean_wb_now_post(wb_post);

            elseif sum(seq_now == n4) > 0
                % wb kin
                f_wb_seq_post_c4(seq_now,wb_post) = f_wb_now_post(wb_post);

                stroke_wb_L_seq_bins_post_c4(:,seq_now,wb_post) = stroke_wb_L_bins_now_post(:,wb_post);
                stroke_wb_R_seq_bins_post_c4(:,seq_now,wb_post) = stroke_wb_R_bins_now_post(:,wb_post);

                dev_wb_L_seq_bins_post_c4(:,seq_now,wb_post) = dev_wb_L_bins_now_post(:,wb_post);
                dev_wb_R_seq_bins_post_c4(:,seq_now,wb_post) = dev_wb_R_bins_now_post(:,wb_post);

                pitch_wb_L_seq_bins_post_c4(:,seq_now,wb_post) = pitch_wb_L_bins_now_post(:,wb_post);
                pitch_wb_R_seq_bins_post_c4(:,seq_now,wb_post) = pitch_wb_R_bins_now_post(:,wb_post);

                % body kin
                roll_mean_wb_seq_post_c4(seq_now,wb_post) = roll_mean_wb_now_post(wb_post);
                roll_dot_mean_wb_seq_post_c4(seq_now,wb_post) = roll_dot_mean_wb_now_post(wb_post);
                roll_dot_dot_mean_wb_seq_post_c4(seq_now,wb_post) = roll_dot_dot_mean_wb_now_post(wb_post);

                pitch_mean_wb_seq_post_c4(seq_now,wb_post) = pitch_mean_wb_now_post(wb_post);
                pitch_dot_mean_wb_seq_post_c4(seq_now,wb_post) = pitch_dot_mean_wb_now_post(wb_post);
                pitch_dot_dot_mean_wb_seq_post_c4(seq_now,wb_post) = pitch_dot_dot_mean_wb_now_post(wb_post);

                yaw_mean_wb_seq_post_c4(seq_now,wb_post) = yaw_mean_wb_now_post(wb_post);
                yaw_dot_mean_wb_seq_post_c4(seq_now,wb_post) = yaw_dot_mean_wb_now_post(wb_post);
                yaw_dot_dot_mean_wb_seq_post_c4(seq_now,wb_post) = yaw_dot_dot_mean_wb_now_post(wb_post);

%                 drot_L_mean_wb_seq_post_c4(seq_now,wb_post) = drot_L_mean_wb_now_post(wb_post);
%                 rot_dot_L_mean_wb_seq_post_c4(seq_now,wb_post) = rot_dot_L_mean_wb_now_post(wb_post);
%                 rot_dot_dot_L_mean_wb_seq_post_c4(seq_now,wb_post) = rot_dot_dot_L_mean_wb_now_post(wb_post);
% 
%                 drot_R_mean_wb_seq_post_c4(seq_now,wb_post) = drot_R_mean_wb_now_post(wb_post);
%                 rot_dot_R_mean_wb_seq_post_c4(seq_now,wb_post) = rot_dot_R_mean_wb_now_post(wb_post);
%                 rot_dot_dot_R_mean_wb_seq_post_c4(seq_now,wb_post) = rot_dot_dot_R_mean_wb_now_post(wb_post);
% 
%                 % torque
%                 Mroll_mean_wb_accel_seq_post_c4(seq_now,wb_post) = Mroll_mean_wb_accel_now_post(wb_post);
%                 Mroll_mean_wb_damp_seq_post_c4(seq_now,wb_post) = Mroll_mean_wb_damp_now_post(wb_post);
%                 Mroll_mean_wb_seq_post_c4(seq_now,wb_post) = Mroll_mean_wb_now_post(wb_post);
% 
%                 Mpitch_mean_wb_accel_seq_post_c4(seq_now,wb_post) = Mpitch_mean_wb_accel_now_post(wb_post);
%                 Mpitch_mean_wb_damp_seq_post_c4(seq_now,wb_post) = Mpitch_mean_wb_damp_now_post(wb_post);
%                 Mpitch_mean_wb_seq_post_c4(seq_now,wb_post) = Mpitch_mean_wb_now_post(wb_post);
% 
%                 Myaw_mean_wb_accel_seq_post_c4(seq_now,wb_post) = Myaw_mean_wb_accel_now_post(wb_post);
%                 Myaw_mean_wb_damp_seq_post_c4(seq_now,wb_post) = Myaw_mean_wb_damp_now_post(wb_post);
%                 Myaw_mean_wb_seq_post_c4(seq_now,wb_post) = Myaw_mean_wb_now_post(wb_post);
% 
%                 M_L_mean_wb_accel_seq_post_c4(seq_now,wb_post) = M_L_mean_wb_accel_now_post(wb_post);
%                 M_L_mean_wb_damp_seq_post_c4(seq_now,wb_post) = M_L_mean_wb_damp_now_post(wb_post);
%                 M_L_mean_wb_seq_post_c4(seq_now,wb_post) = M_L_mean_wb_now_post(wb_post);
% 
%                 M_R_mean_wb_accel_seq_post_c4(seq_now,wb_post) = M_R_mean_wb_accel_now_post(wb_post);
%                 M_R_mean_wb_damp_seq_post_c4(seq_now,wb_post) = M_R_mean_wb_damp_now_post(wb_post);
%                 M_R_mean_wb_seq_post_c4(seq_now,wb_post) = M_R_mean_wb_now_post(wb_post);
            end
        end
    end
end

%% separate mean wingbeat seq throughout maneuver
define_vars_timeline_means
define_vars_timeline_means_c1
define_vars_timeline_means_c2
define_vars_timeline_means_c3
define_vars_timeline_means_c4

for wb_now = 1:wb_max
    calc_vars_timeline_means_pre
    calc_vars_timeline_means_post
end
calc_vars_timeline_means_preNpost

for wb_now = 1:wb_max
    calc_vars_timeline_means_pre_c1
    calc_vars_timeline_means_post_c1
end
calc_vars_timeline_means_preNpost_c1


for wb_now = 1:wb_max
    calc_vars_timeline_means_pre_c2
    calc_vars_timeline_means_post_c2
end
calc_vars_timeline_means_preNpost_c2


for wb_now = 1:wb_max
    calc_vars_timeline_means_pre_c3
    calc_vars_timeline_means_post_c3
end
calc_vars_timeline_means_preNpost_c3


for wb_now = 1:wb_max
    calc_vars_timeline_means_pre_c4
    calc_vars_timeline_means_post_c4
end
calc_vars_timeline_means_preNpost_c4

%% normalized time
% bin time
t_norm_wb_seq_bins_pre_mean_all = [0:length(stroke_wb_L_seq_bins_pre_mean_all)-1]'/(n_bins-1);
t_norm_wb_seq_bins_pre_mean_all = -t_norm_wb_seq_bins_pre_mean_all;
t_norm_wb_seq_bins_pre_mean_all = flipud(t_norm_wb_seq_bins_pre_mean_all);

t_norm_wb_seq_bins_post_mean_all = [0:length(stroke_wb_L_seq_bins_post_mean_all)-1]'/(n_bins-1);

t_norm_wb_seq_bins_mean_all = [t_norm_wb_seq_bins_pre_mean_all;t_norm_wb_seq_bins_post_mean_all];

% wb average time
t_norm_wb_seq_pre_mean_all = [1:wb_max]';
t_norm_wb_seq_pre_mean_all = -t_norm_wb_seq_pre_mean_all;
t_norm_wb_seq_pre_mean_all = flipud(t_norm_wb_seq_pre_mean_all);

t_norm_wb_seq_post_mean_all = [0:wb_max-1]';

t_norm_wb_seq_mean_all = [t_norm_wb_seq_pre_mean_all;t_norm_wb_seq_post_mean_all];

%% real time
t_norm = [1:n_bins]'/(n_bins-1);
t_norm = t_norm(1:end-1);

% pre
t_wb_seq_bins_pre_mean_all = [];
t_wb_seq_pre_mean_all = [];
for i = 1:length(f_wb_seq_pre_mean_all)
    if isnan(f_wb_seq_pre_mean_all(i)) == 0
        
        f_now = f_wb_seq_pre_mean_all(i);
        t_now = t_norm / f_now;
        
        if ~isempty(t_wb_seq_bins_pre_mean_all) && isnan(t_wb_seq_bins_pre_mean_all(end))==0
            t_now = t_now + t_wb_seq_bins_pre_mean_all(end);
        end
    else
        t_now = nan(length(t_norm),1);
    end
    t_wb_seq_bins_pre_mean_all = [t_wb_seq_bins_pre_mean_all;t_now];
    t_wb_seq_pre_mean_all = [t_wb_seq_pre_mean_all;t_now(1)];
end
t_wb_seq_pre_mean_all = t_wb_seq_pre_mean_all - max(t_wb_seq_bins_pre_mean_all);
t_wb_seq_bins_pre_mean_all = t_wb_seq_bins_pre_mean_all - max(t_wb_seq_bins_pre_mean_all);


% post
t_wb_seq_bins_post_mean_all = [];
t_wb_seq_post_mean_all = [];
for i = 1:length(f_wb_seq_post_mean_all)
    if isnan(f_wb_seq_post_mean_all(i)) == 0
        
        f_now = f_wb_seq_post_mean_all(i);
        t_now = t_norm / f_now;
        
        if ~isempty(t_wb_seq_bins_post_mean_all) && isnan(t_wb_seq_bins_post_mean_all(end))==0
            t_now = t_now + t_wb_seq_bins_post_mean_all(end);
        end
    else
        t_now = nan(length(t_norm),1);
    end
    t_wb_seq_bins_post_mean_all = [t_wb_seq_bins_post_mean_all;t_now];
    t_wb_seq_post_mean_all = [t_wb_seq_post_mean_all;t_now(1)];
end

t_wb_seq_bins_mean_all = [t_wb_seq_bins_pre_mean_all;t_wb_seq_bins_post_mean_all];
t_wb_seq_mean_all = [t_wb_seq_pre_mean_all;t_wb_seq_post_mean_all];

%% wingbeat characteristics
calc_WBparams_timeline
calc_WBparams_timeline_c1
calc_WBparams_timeline_c2
calc_WBparams_timeline_c3
calc_WBparams_timeline_c4

%% save data
save('WBdataset_temporal_dynamics_INCclusters.mat')

% save('WBdataset_temporal_dynamics_TorqueNorm_subset.mat',...
%     'stroke_wb_L_seq_bins_mean_all','stroke_wb_L_seq_bins_mean',...
%     'dev_wb_L_seq_bins_mean_all','dev_wb_L_seq_bins_mean',...
%     'pitch_wb_L_seq_bins_mean_all','pitch_wb_L_seq_bins_mean',...
%     'f_wb_seq_mean_all','t_wb_seq_bins_mean_all','t_wb_seq_mean_all',...
%     't_norm_wb_seq_bins_mean_all','t_norm_wb_seq_mean_all')

%% PLOT ALL&MEAN DATA
% 
% clear
% load('WBdataset_temporal_dynamics_TorqueNorm_clusters.mat')
% 
% mkdir('MSfigs_WBnBodyKin_TempDynamics')
% cd('MSfigs_WBnBodyKin_TempDynamics')
% 
% plot_wb_seq_wbkin_allNmeanN95pCI_savefigs_LOOMING
% plot_wb_seq_bodydyn_allNmeanN95pCI_savefigs_LOOMING
% plot_wb_seq_torque_allNmeanN95pCI_savefigs_LOOMING
% 
% %% plot
% 
% %% stroke
% figure
% 
% % L+R
% subplot(3,1,1)
% hold on
% % plot(t_wb_seq_mean_all,stroke_wb_seq_bins_MAX,'-r.')
% % plot(t_wb_seq_mean_all,-stroke_wb_seq_bins_MIN,'-b.')
% % plot(t_wb_seq_mean_all,Astroke_wb_seq_bins,'-g.')
% plot(t_wb_seq_mean_all,Astroke_wb_R_seq_bins,'-c.')
% plot(t_wb_seq_mean_all,Astroke_wb_L_seq_bins,'-m.')
% % plot(t_wb_seq_mean_all,Astroke_wb_seq_bins,'-g.')
% 
% % xlabel('time [sec]')
% xlim([-.05 .05]) 
% ylabel('stroke angle amplitude [deg]')
% ylim([120 150])
% set(gca,'YTick',[-90:10:190])
% 
% % L-R
% subplot(3,1,2)
% hold on
% plot(t_wb_seq_mean_all,Dstroke_wb_seq_bins_MAX,'-r.')
% plot(t_wb_seq_mean_all,Dstroke_wb_seq_bins_MIN,'-b.')
% plot(t_wb_seq_mean_all,-dAstroke_wb_seq_bins,'-g.')
% 
% xlabel('time [sec]')
% xlim([-.05 .05]) 
% ylabel('stroke angle L-R [deg]')
% ylim([-5 5])
% set(gca,'YTick',[-90:5:90])
% 
% saveas(gca,['WBvsTimeSeqs_strokeParams.fig'])
% saveas(gca,['WBvsTimeSeqs_strokeParams.png'])
% plot2svg(['WBvsTimeSeqs_strokeParams.svg'])
% 
% %% pitch
% figure
% 
% % Downstroke
% subplot(3,1,1)
% hold on
% plot(t_wb_seq_mean_all,pitch_wb_R_seq_bins_MAXmidDS,'-c.')
% plot(t_wb_seq_mean_all,pitch_wb_L_seq_bins_MAXmidDS,'-m.')
% % plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MAXmidDS,'-g.')
% 
% % xlabel('time [sec]')
% xlim([-.05 .05]) 
% ylabel('rotation angle downstroke [deg]')
% ylim([50 60])
% set(gca,'YTick',[-90:5:90])
% 
% % Upstroke
% subplot(3,1,2)
% hold on
% plot(t_wb_seq_mean_all,pitch_wb_R_seq_bins_MINmidUS,'-c.')
% plot(t_wb_seq_mean_all,pitch_wb_L_seq_bins_MINmidUS,'-m.')
% % plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MINmidUS,'-g.')
% 
% % xlabel('time [sec]')
% xlim([-.05 .05]) 
% ylabel('rotation angle upstroke [deg]')
% ylim([-50 -40])
% set(gca,'YTick',[-90:5:90])
% 
% % L-R
% subplot(3,1,3)
% hold on
% plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MAXmidDS,'-r.')
% plot(t_wb_seq_mean_all,Dpitch_wb_seq_bins_MINmidUS,'-b.')
% % plot(t_wb_seq_mean_all,-dApitch_wb_seq_bins,'-g.')
% 
% xlabel('time [sec]')
% xlim([-.05 .05]) 
% ylabel('rotation angle L-R [deg]')
% ylim([-5 5])
% set(gca,'YTick',[-90:5:90])
% 
% saveas(gca,['WBvsTimeSeqs_rotationParams.fig'])
% saveas(gca,['WBvsTimeSeqs_rotationParams.png'])
% plot2svg(['WBvsTimeSeqs_rotationParams.svg'])
% 
% %% dev
% figure
% 
% % Downstroke
% subplot(3,1,1)
% hold on
% % plot(t_wb_seq_mean_all,dev_wb_R_seq_bins_MAXds,'-c.')
% % plot(t_wb_seq_mean_all,dev_wb_L_seq_bins_MAXds,'-m.')
% % % plot(t_wb_seq_mean_all,Ddev_wb_seq_bins_MAXds,'-g.')
% 
% plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_DS,'-c.')
% plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_DS,'-m.')
% % 
% % plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_US,'-c.')
% % plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_US,'-m.')
% 
% % xlabel('time [sec]')
% xlim([-.05 .05]) 
% ylabel('deviation amplitude downstroke [deg]')
% ylim([15 25])
% set(gca,'YTick',[-90:5:90])
% 
% % Upstroke
% subplot(3,1,2)
% hold on
% % plot(t_wb_seq_mean_all,dev_wb_R_seq_bins_MAXds,'-c.')
% % plot(t_wb_seq_mean_all,dev_wb_L_seq_bins_MAXds,'-m.')
% % % plot(t_wb_seq_mean_all,Ddev_wb_seq_bins_MAXds,'-g.')
% % 
% % plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_DS,'-r.')
% % plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_DS,'-b.')
% 
% plot(t_wb_seq_mean_all,Adev_wb_R_seq_bins_US,'-c.')
% plot(t_wb_seq_mean_all,Adev_wb_L_seq_bins_US,'-m.')
% 
% % xlabel('time [sec]')
% xlim([-.05 .05]) 
% ylabel('deviation amplitude upstroke [deg]')
% ylim([10 20])
% set(gca,'YTick',[-90:5:90])
% 
% % L-R
% subplot(3,1,3)
% hold on
% plot(t_wb_seq_mean_all,dAdev_wb_seq_bins_DS,'-r.')
% plot(t_wb_seq_mean_all,dAdev_wb_seq_bins_US,'-b.')
% % plot(t_wb_seq_mean_all,-dAdev_wb_seq_bins,'-g.')
% 
% xlabel('time [sec]')
% xlim([-.05 .05]) 
% ylabel('deviation amplitude L-R [deg]')
% ylim([-5 5])
% set(gca,'YTick',[-90:5:90])
% 
% saveas(gca,['WBvsTimeSeqs_deviationParams.fig'])
% saveas(gca,['WBvsTimeSeqs_deviationParams.png'])
% plot2svg(['WBvsTimeSeqs_deviationParams.svg'])
% 
% 
% cd ..
% 
% 




