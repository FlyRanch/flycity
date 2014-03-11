clc
clear
% close all

addpath('/home/florian/Dropbox/WORK/toolbox')
addpath('/home/florian/Dropbox/WORK/toolbox/plot2svg')
addpath('/home/florian/Dropbox/WORK/toolbox/CircStat')
addpath('/home/florian/Dropbox/WORK/toolbox/flytracker')

%% load data
loadname=dir('WBdataset_all_*')
loadname=loadname.name;
load(loadname)

steady_name=dir('WBdataset_steady_*')
steady_name=steady_name.name;
load(steady_name)

%% settings
n_bins = size(stroke_wb_L_bins,1);
wb_max = 30;
wb_min4mean = max(seq_nr)/2; % only mean if half of the sequences are providing data


%% wb kin vars
f_wb = (f_wb_L + f_wb_R)./2;
f_wb_seq_pre = nan(max(seq_nr),wb_max);
f_wb_seq_post = nan(max(seq_nr),wb_max);

stroke_wb_L_seq_bins_pre = nan(n_bins,max(seq_nr),wb_max);
dev_wb_L_seq_bins_pre = nan(n_bins,max(seq_nr),wb_max);
pitch_wb_L_seq_bins_pre = nan(n_bins,max(seq_nr),wb_max);

stroke_wb_R_seq_bins_pre = nan(n_bins,max(seq_nr),wb_max);
dev_wb_R_seq_bins_pre = nan(n_bins,max(seq_nr),wb_max);
pitch_wb_R_seq_bins_pre = nan(n_bins,max(seq_nr),wb_max);

stroke_wb_L_seq_bins_post = nan(n_bins,max(seq_nr),wb_max);
dev_wb_L_seq_bins_post = nan(n_bins,max(seq_nr),wb_max);
pitch_wb_L_seq_bins_post = nan(n_bins,max(seq_nr),wb_max);

stroke_wb_R_seq_bins_post = nan(n_bins,max(seq_nr),wb_max);
dev_wb_R_seq_bins_post = nan(n_bins,max(seq_nr),wb_max);
pitch_wb_R_seq_bins_post = nan(n_bins,max(seq_nr),wb_max);

%% body kin vars
% pre
roll_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
roll_dot_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
roll_dot_dot_mean_wb_seq_pre = nan(max(seq_nr),wb_max);

pitch_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
pitch_dot_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
pitch_dot_dot_mean_wb_seq_pre = nan(max(seq_nr),wb_max);

yaw_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
yaw_dot_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
yaw_dot_dot_mean_wb_seq_pre = nan(max(seq_nr),wb_max);

rot_L_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
rot_dot_L_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
rot_dot_dot_L_mean_wb_seq_pre = nan(max(seq_nr),wb_max);

rot_R_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
rot_dot_R_mean_wb_seq_pre = nan(max(seq_nr),wb_max);
rot_dot_dot_R_mean_wb_seq_pre = nan(max(seq_nr),wb_max);

% post
roll_mean_wb_seq_post = nan(max(seq_nr),wb_max);
roll_dot_mean_wb_seq_post = nan(max(seq_nr),wb_max);
roll_dot_dot_mean_wb_seq_post = nan(max(seq_nr),wb_max);

pitch_mean_wb_seq_post = nan(max(seq_nr),wb_max);
pitch_dot_mean_wb_seq_post = nan(max(seq_nr),wb_max);
pitch_dot_dot_mean_wb_seq_post = nan(max(seq_nr),wb_max);

yaw_mean_wb_seq_post = nan(max(seq_nr),wb_max);
yaw_dot_mean_wb_seq_post = nan(max(seq_nr),wb_max);
yaw_dot_dot_mean_wb_seq_post = nan(max(seq_nr),wb_max);

rot_L_mean_wb_seq_post = nan(max(seq_nr),wb_max);
rot_dot_L_mean_wb_seq_post = nan(max(seq_nr),wb_max);
rot_dot_dot_L_mean_wb_seq_post = nan(max(seq_nr),wb_max);

rot_R_mean_wb_seq_post = nan(max(seq_nr),wb_max);
rot_dot_R_mean_wb_seq_post = nan(max(seq_nr),wb_max);
rot_dot_dot_R_mean_wb_seq_post = nan(max(seq_nr),wb_max);

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
    
    drot_L_mean_wb_now = drot_L_mean_wb(seq_nr==seq_now);
    drot_L_mean_wb_now_pre = drot_L_mean_wb_now(frame_nr_2resp_now<0);
    drot_L_mean_wb_now_post = drot_L_mean_wb_now(frame_nr_2resp_now>=0);
    
    rot_dot_L_mean_wb_now = rot_dot_L_mean_wb(seq_nr==seq_now);
    rot_dot_L_mean_wb_now_pre = rot_dot_L_mean_wb_now(frame_nr_2resp_now<0);
    rot_dot_L_mean_wb_now_post = rot_dot_L_mean_wb_now(frame_nr_2resp_now>=0);
    
    rot_dot_dot_L_mean_wb_now = rot_dot_dot_L_mean_wb(seq_nr==seq_now);
    rot_dot_dot_L_mean_wb_now_pre = rot_dot_dot_L_mean_wb_now(frame_nr_2resp_now<0);
    rot_dot_dot_L_mean_wb_now_post = rot_dot_dot_L_mean_wb_now(frame_nr_2resp_now>=0);
    
    drot_R_mean_wb_now = drot_R_mean_wb(seq_nr==seq_now);
    drot_R_mean_wb_now_pre = drot_R_mean_wb_now(frame_nr_2resp_now<0);
    drot_R_mean_wb_now_post = drot_R_mean_wb_now(frame_nr_2resp_now>=0);
    
    rot_dot_R_mean_wb_now = rot_dot_R_mean_wb(seq_nr==seq_now);
    rot_dot_R_mean_wb_now_pre = rot_dot_R_mean_wb_now(frame_nr_2resp_now<0);
    rot_dot_R_mean_wb_now_post = rot_dot_R_mean_wb_now(frame_nr_2resp_now>=0);
    
    rot_dot_dot_R_mean_wb_now = rot_dot_dot_R_mean_wb(seq_nr==seq_now);
    rot_dot_dot_R_mean_wb_now_pre = rot_dot_dot_R_mean_wb_now(frame_nr_2resp_now<0);
    rot_dot_dot_R_mean_wb_now_post = rot_dot_dot_R_mean_wb_now(frame_nr_2resp_now>=0);
    
    %% pre response wbs
    wb_pre_max = length(frame_nr_2resp_now_pre);
    for wb_now = 1:wb_pre_max
        
        wb_pre = wb_max + wb_now - wb_pre_max;
        
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
        
        drot_L_mean_wb_seq_pre(seq_now,wb_pre) = drot_L_mean_wb_now_pre(wb_now);
        rot_dot_L_mean_wb_seq_pre(seq_now,wb_pre) = rot_dot_L_mean_wb_now_pre(wb_now);
        rot_dot_dot_L_mean_wb_seq_pre(seq_now,wb_pre) = rot_dot_dot_L_mean_wb_now_pre(wb_now);
        
        drot_R_mean_wb_seq_pre(seq_now,wb_pre) = drot_R_mean_wb_now_pre(wb_now);
        rot_dot_R_mean_wb_seq_pre(seq_now,wb_pre) = rot_dot_R_mean_wb_now_pre(wb_now);
        rot_dot_dot_R_mean_wb_seq_pre(seq_now,wb_pre) = rot_dot_dot_R_mean_wb_now_pre(wb_now);
    end
    
    %% post response wbs
    wb_post_max = length(frame_nr_2resp_now_post);
    for wb_post = 1:wb_post_max
        
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
        
        drot_L_mean_wb_seq_post(seq_now,wb_post) = drot_L_mean_wb_now_post(wb_post);
        rot_dot_L_mean_wb_seq_post(seq_now,wb_post) = rot_dot_L_mean_wb_now_post(wb_post);
        rot_dot_dot_L_mean_wb_seq_post(seq_now,wb_post) = rot_dot_dot_L_mean_wb_now_post(wb_post);
        
        drot_R_mean_wb_seq_post(seq_now,wb_post) = drot_R_mean_wb_now_post(wb_post);
        rot_dot_R_mean_wb_seq_post(seq_now,wb_post) = rot_dot_R_mean_wb_now_post(wb_post);
        rot_dot_dot_R_mean_wb_seq_post(seq_now,wb_post) = rot_dot_dot_R_mean_wb_now_post(wb_post);
    end
end

%% mean wingbeat seq throughout maneuver

% wb kin
f_wb_seq_pre_mean = [];
f_wb_seq_post_mean = [];

stroke_wb_L_seq_bins_pre_mean = [];
stroke_wb_R_seq_bins_pre_mean = [];
stroke_wb_L_seq_bins_post_mean = [];
stroke_wb_R_seq_bins_post_mean = [];

dev_wb_L_seq_bins_pre_mean = [];
dev_wb_R_seq_bins_pre_mean = [];
dev_wb_L_seq_bins_post_mean = [];
dev_wb_R_seq_bins_post_mean = [];

pitch_wb_L_seq_bins_pre_mean = [];
pitch_wb_R_seq_bins_pre_mean = [];
pitch_wb_L_seq_bins_post_mean = [];
pitch_wb_R_seq_bins_post_mean = [];

% body kin
roll_mean_wb_seq_pre_mean = [];
roll_dot_mean_wb_seq_pre_mean = [];
roll_dot_dot_mean_wb_seq_pre_mean = [];

roll_mean_wb_seq_post_mean = [];
roll_dot_mean_wb_seq_post_mean = [];
roll_dot_dot_mean_wb_seq_post_mean = [];

pitch_mean_wb_seq_pre_mean = [];
pitch_dot_mean_wb_seq_pre_mean = [];
pitch_dot_dot_mean_wb_seq_pre_mean = [];

pitch_mean_wb_seq_post_mean = [];
pitch_dot_mean_wb_seq_post_mean = [];
pitch_dot_dot_mean_wb_seq_post_mean = [];

yaw_mean_wb_seq_pre_mean = [];
yaw_dot_mean_wb_seq_pre_mean = [];
yaw_dot_dot_mean_wb_seq_pre_mean = [];

yaw_mean_wb_seq_post_mean = [];
yaw_dot_mean_wb_seq_post_mean = [];
yaw_dot_dot_mean_wb_seq_post_mean = [];

drot_L_mean_wb_seq_pre_mean = [];
rot_dot_L_mean_wb_seq_pre_mean = [];
rot_dot_dot_L_mean_wb_seq_pre_mean = [];

drot_L_mean_wb_seq_post_mean = [];
rot_dot_L_mean_wb_seq_post_mean = [];
rot_dot_dot_L_mean_wb_seq_post_mean = [];

drot_R_mean_wb_seq_pre_mean = [];
rot_dot_R_mean_wb_seq_pre_mean = [];
rot_dot_dot_R_mean_wb_seq_pre_mean = [];

drot_R_mean_wb_seq_post_mean = [];
rot_dot_R_mean_wb_seq_post_mean = [];
rot_dot_dot_R_mean_wb_seq_post_mean = [];

for wb_now = 1:wb_max

    %% pre
    if sum(isnan(stroke_wb_L_seq_bins_pre(1,:,wb_now))==0) >= wb_min4mean    % only mean if Nwb>wb_min4mean

        % wb kin
        f_wb_mean_now = nanmean(f_wb_seq_pre(:,wb_now)')';
        
        stroke_wb_L_mean_now = nanmean(stroke_wb_L_seq_bins_pre(:,:,wb_now)')';
        stroke_wb_R_mean_now = nanmean(stroke_wb_R_seq_bins_pre(:,:,wb_now)')';
        
        dev_wb_L_mean_now = nanmean(dev_wb_L_seq_bins_pre(:,:,wb_now)')';
        dev_wb_R_mean_now = nanmean(dev_wb_R_seq_bins_pre(:,:,wb_now)')';
        
        pitch_wb_L_mean_now = nanmean(pitch_wb_L_seq_bins_pre(:,:,wb_now)')';
        pitch_wb_R_mean_now = nanmean(pitch_wb_R_seq_bins_pre(:,:,wb_now)')';
        
        % body kin
        roll_mean_wb_mean_now = nanmean(roll_mean_wb_seq_pre(:,wb_now)')';
        roll_dot_mean_wb_mean_now = nanmean(roll_dot_mean_wb_seq_pre(:,wb_now)')';
        roll_dot_dot_mean_wb_mean_now = nanmean(roll_dot_dot_mean_wb_seq_pre(:,wb_now)')';
        
        pitch_mean_wb_mean_now = nanmean(pitch_mean_wb_seq_pre(:,wb_now)')';
        pitch_dot_mean_wb_mean_now = nanmean(pitch_dot_mean_wb_seq_pre(:,wb_now)')';
        pitch_dot_dot_mean_wb_mean_now = nanmean(pitch_dot_dot_mean_wb_seq_pre(:,wb_now)')';
        
        yaw_mean_wb_mean_now = nanmean(yaw_mean_wb_seq_pre(:,wb_now)')';
        yaw_dot_mean_wb_mean_now = nanmean(yaw_dot_mean_wb_seq_pre(:,wb_now)')';
        yaw_dot_dot_mean_wb_mean_now = nanmean(yaw_dot_dot_mean_wb_seq_pre(:,wb_now)')';
        
        drot_L_mean_wb_mean_now = nanmean(drot_L_mean_wb_seq_pre(:,wb_now)')';
        rot_dot_L_mean_wb_mean_now = nanmean(rot_dot_L_mean_wb_seq_pre(:,wb_now)')';
        rot_dot_dot_L_mean_wb_mean_now = nanmean(rot_dot_dot_L_mean_wb_seq_pre(:,wb_now)')';
        
        drot_R_mean_wb_mean_now = nanmean(drot_R_mean_wb_seq_pre(:,wb_now)')';
        rot_dot_R_mean_wb_mean_now = nanmean(rot_dot_R_mean_wb_seq_pre(:,wb_now)')';
        rot_dot_dot_R_mean_wb_mean_now = nanmean(rot_dot_dot_R_mean_wb_seq_pre(:,wb_now)')';
        
        %% add to list
        % wb kin
        f_wb_seq_pre_mean = [f_wb_seq_pre_mean; f_wb_mean_now];
        
        stroke_wb_L_seq_bins_pre_mean = [stroke_wb_L_seq_bins_pre_mean; stroke_wb_L_mean_now(1:end-1)];
        stroke_wb_R_seq_bins_pre_mean = [stroke_wb_R_seq_bins_pre_mean; stroke_wb_R_mean_now(1:end-1)];
        
        dev_wb_L_seq_bins_pre_mean = [dev_wb_L_seq_bins_pre_mean; dev_wb_L_mean_now(1:end-1)];
        dev_wb_R_seq_bins_pre_mean = [dev_wb_R_seq_bins_pre_mean; dev_wb_R_mean_now(1:end-1)];
        
        pitch_wb_L_seq_bins_pre_mean = [pitch_wb_L_seq_bins_pre_mean; pitch_wb_L_mean_now(1:end-1)];
        pitch_wb_R_seq_bins_pre_mean = [pitch_wb_R_seq_bins_pre_mean; pitch_wb_R_mean_now(1:end-1)];
        
        % body kin
        roll_mean_wb_seq_pre_mean = [roll_mean_wb_seq_pre_mean; roll_mean_wb_mean_now];
        roll_dot_mean_wb_seq_pre_mean = [roll_dot_mean_wb_seq_pre_mean; roll_dot_mean_wb_mean_now];
        roll_dot_dot_mean_wb_seq_pre_mean = [roll_dot_dot_mean_wb_seq_pre_mean; roll_dot_dot_mean_wb_mean_now];
        
        pitch_mean_wb_seq_pre_mean = [pitch_mean_wb_seq_pre_mean; pitch_mean_wb_mean_now];
        pitch_dot_mean_wb_seq_pre_mean = [pitch_dot_mean_wb_seq_pre_mean; pitch_dot_mean_wb_mean_now];
        pitch_dot_dot_mean_wb_seq_pre_mean = [pitch_dot_dot_mean_wb_seq_pre_mean; pitch_dot_dot_mean_wb_mean_now];
        
        yaw_mean_wb_seq_pre_mean = [yaw_mean_wb_seq_pre_mean; yaw_mean_wb_mean_now];
        yaw_dot_mean_wb_seq_pre_mean = [yaw_dot_mean_wb_seq_pre_mean; yaw_dot_mean_wb_mean_now];
        yaw_dot_dot_mean_wb_seq_pre_mean = [yaw_dot_dot_mean_wb_seq_pre_mean; yaw_dot_dot_mean_wb_mean_now];
        
        drot_L_mean_wb_seq_pre_mean = [drot_L_mean_wb_seq_pre_mean; drot_L_mean_wb_mean_now];
        rot_dot_L_mean_wb_seq_pre_mean = [rot_dot_L_mean_wb_seq_pre_mean; rot_dot_L_mean_wb_mean_now];
        rot_dot_dot_L_mean_wb_seq_pre_mean = [rot_dot_dot_L_mean_wb_seq_pre_mean; rot_dot_dot_L_mean_wb_mean_now];
        
        drot_R_mean_wb_seq_pre_mean = [drot_R_mean_wb_seq_pre_mean; drot_R_mean_wb_mean_now];
        rot_dot_R_mean_wb_seq_pre_mean = [rot_dot_R_mean_wb_seq_pre_mean; rot_dot_R_mean_wb_mean_now];
        rot_dot_dot_R_mean_wb_seq_pre_mean = [rot_dot_dot_R_mean_wb_seq_pre_mean; rot_dot_dot_R_mean_wb_mean_now];
        
    else
        % wb kin
        f_wb_seq_pre_mean = [f_wb_seq_pre_mean; nan];

        stroke_wb_L_seq_bins_pre_mean = [stroke_wb_L_seq_bins_pre_mean; nan(n_bins-1,1)];
        stroke_wb_R_seq_bins_pre_mean = [stroke_wb_R_seq_bins_pre_mean; nan(n_bins-1,1)];
        
        dev_wb_L_seq_bins_pre_mean = [dev_wb_L_seq_bins_pre_mean; nan(n_bins-1,1)];
        dev_wb_R_seq_bins_pre_mean = [dev_wb_R_seq_bins_pre_mean; nan(n_bins-1,1)];
        
        pitch_wb_L_seq_bins_pre_mean = [pitch_wb_L_seq_bins_pre_mean; nan(n_bins-1,1)];
        pitch_wb_R_seq_bins_pre_mean = [pitch_wb_R_seq_bins_pre_mean; nan(n_bins-1,1)];
        
        % body kin
        roll_mean_wb_seq_pre_mean = [roll_mean_wb_seq_pre_mean; nan];
        roll_dot_mean_wb_seq_pre_mean = [roll_dot_mean_wb_seq_pre_mean; nan];
        roll_dot_dot_mean_wb_seq_pre_mean = [roll_dot_dot_mean_wb_seq_pre_mean; nan];
        
        pitch_mean_wb_seq_pre_mean = [pitch_mean_wb_seq_pre_mean; nan];
        pitch_dot_mean_wb_seq_pre_mean = [pitch_dot_mean_wb_seq_pre_mean; nan];
        pitch_dot_dot_mean_wb_seq_pre_mean = [pitch_dot_dot_mean_wb_seq_pre_mean; nan];
        
        yaw_mean_wb_seq_pre_mean = [yaw_mean_wb_seq_pre_mean; nan];
        yaw_dot_mean_wb_seq_pre_mean = [yaw_dot_mean_wb_seq_pre_mean; nan];
        yaw_dot_dot_mean_wb_seq_pre_mean = [yaw_dot_dot_mean_wb_seq_pre_mean; nan];
        
        drot_L_mean_wb_seq_pre_mean = [drot_L_mean_wb_seq_pre_mean; nan];
        rot_dot_L_mean_wb_seq_pre_mean = [rot_dot_L_mean_wb_seq_pre_mean; nan];
        rot_dot_dot_L_mean_wb_seq_pre_mean = [rot_dot_dot_L_mean_wb_seq_pre_mean; nan];
        
        drot_R_mean_wb_seq_pre_mean = [drot_R_mean_wb_seq_pre_mean; nan];
        rot_dot_R_mean_wb_seq_pre_mean = [rot_dot_R_mean_wb_seq_pre_mean; nan];
        rot_dot_dot_R_mean_wb_seq_pre_mean = [rot_dot_dot_R_mean_wb_seq_pre_mean; nan];
    end
    
    %% post
    if sum(isnan(stroke_wb_L_seq_bins_post(1,:,wb_now))==0) >= wb_min4mean    % only mean if Nwb>wb_min4mean
        
        % wb kin
        f_wb_mean_now = nanmean(f_wb_seq_post(:,wb_now)')';
        
        stroke_wb_L_mean_now = nanmean(stroke_wb_L_seq_bins_post(:,:,wb_now)')';
        stroke_wb_R_mean_now = nanmean(stroke_wb_R_seq_bins_post(:,:,wb_now)')';
        
        dev_wb_L_mean_now = nanmean(dev_wb_L_seq_bins_post(:,:,wb_now)')';
        dev_wb_R_mean_now = nanmean(dev_wb_R_seq_bins_post(:,:,wb_now)')';
        
        pitch_wb_L_mean_now = nanmean(pitch_wb_L_seq_bins_post(:,:,wb_now)')';
        pitch_wb_R_mean_now = nanmean(pitch_wb_R_seq_bins_post(:,:,wb_now)')';
        
        % body kin
        roll_mean_wb_mean_now = nanmean(roll_mean_wb_seq_post(:,wb_now)')';
        roll_dot_mean_wb_mean_now = nanmean(roll_dot_mean_wb_seq_post(:,wb_now)')';
        roll_dot_dot_mean_wb_mean_now = nanmean(roll_dot_dot_mean_wb_seq_post(:,wb_now)')';
        
        pitch_mean_wb_mean_now = nanmean(pitch_mean_wb_seq_post(:,wb_now)')';
        pitch_dot_mean_wb_mean_now = nanmean(pitch_dot_mean_wb_seq_post(:,wb_now)')';
        pitch_dot_dot_mean_wb_mean_now = nanmean(pitch_dot_dot_mean_wb_seq_post(:,wb_now)')';
        
        yaw_mean_wb_mean_now = nanmean(yaw_mean_wb_seq_post(:,wb_now)')';
        yaw_dot_mean_wb_mean_now = nanmean(yaw_dot_mean_wb_seq_post(:,wb_now)')';
        yaw_dot_dot_mean_wb_mean_now = nanmean(yaw_dot_dot_mean_wb_seq_post(:,wb_now)')';
        
        drot_L_mean_wb_mean_now = nanmean(drot_L_mean_wb_seq_post(:,wb_now)')';
        rot_dot_L_mean_wb_mean_now = nanmean(rot_dot_L_mean_wb_seq_post(:,wb_now)')';
        rot_dot_dot_L_mean_wb_mean_now = nanmean(rot_dot_dot_L_mean_wb_seq_post(:,wb_now)')';
        
        drot_R_mean_wb_mean_now = nanmean(drot_R_mean_wb_seq_post(:,wb_now)')';
        rot_dot_R_mean_wb_mean_now = nanmean(rot_dot_R_mean_wb_seq_post(:,wb_now)')';
        rot_dot_dot_R_mean_wb_mean_now = nanmean(rot_dot_dot_R_mean_wb_seq_post(:,wb_now)')';
        
        %% add to list
        % wb kin
        f_wb_seq_post_mean = [f_wb_seq_post_mean; f_wb_mean_now];
        
        stroke_wb_L_seq_bins_post_mean = [stroke_wb_L_seq_bins_post_mean; stroke_wb_L_mean_now(1:end-1)];
        stroke_wb_R_seq_bins_post_mean = [stroke_wb_R_seq_bins_post_mean; stroke_wb_R_mean_now(1:end-1)];
        
        dev_wb_L_seq_bins_post_mean = [dev_wb_L_seq_bins_post_mean; dev_wb_L_mean_now(1:end-1)];
        dev_wb_R_seq_bins_post_mean = [dev_wb_R_seq_bins_post_mean; dev_wb_R_mean_now(1:end-1)];
        
        pitch_wb_L_seq_bins_post_mean = [pitch_wb_L_seq_bins_post_mean; pitch_wb_L_mean_now(1:end-1)];
        pitch_wb_R_seq_bins_post_mean = [pitch_wb_R_seq_bins_post_mean; pitch_wb_R_mean_now(1:end-1)];
        
        % body kin
        roll_mean_wb_seq_post_mean = [roll_mean_wb_seq_post_mean; roll_mean_wb_mean_now];
        roll_dot_mean_wb_seq_post_mean = [roll_dot_mean_wb_seq_post_mean; roll_dot_mean_wb_mean_now];
        roll_dot_dot_mean_wb_seq_post_mean = [roll_dot_dot_mean_wb_seq_post_mean; roll_dot_dot_mean_wb_mean_now];
        
        pitch_mean_wb_seq_post_mean = [pitch_mean_wb_seq_post_mean; pitch_mean_wb_mean_now];
        pitch_dot_mean_wb_seq_post_mean = [pitch_dot_mean_wb_seq_post_mean; pitch_dot_mean_wb_mean_now];
        pitch_dot_dot_mean_wb_seq_post_mean = [pitch_dot_dot_mean_wb_seq_post_mean; pitch_dot_dot_mean_wb_mean_now];
        
        yaw_mean_wb_seq_post_mean = [yaw_mean_wb_seq_post_mean; yaw_mean_wb_mean_now];
        yaw_dot_mean_wb_seq_post_mean = [yaw_dot_mean_wb_seq_post_mean; yaw_dot_mean_wb_mean_now];
        yaw_dot_dot_mean_wb_seq_post_mean = [yaw_dot_dot_mean_wb_seq_post_mean; yaw_dot_dot_mean_wb_mean_now];
        
        drot_L_mean_wb_seq_post_mean = [drot_L_mean_wb_seq_post_mean; drot_L_mean_wb_mean_now];
        rot_dot_L_mean_wb_seq_post_mean = [rot_dot_L_mean_wb_seq_post_mean; rot_dot_L_mean_wb_mean_now];
        rot_dot_dot_L_mean_wb_seq_post_mean = [rot_dot_dot_L_mean_wb_seq_post_mean; rot_dot_dot_L_mean_wb_mean_now];
        
        drot_R_mean_wb_seq_post_mean = [drot_R_mean_wb_seq_post_mean; drot_R_mean_wb_mean_now];
        rot_dot_R_mean_wb_seq_post_mean = [rot_dot_R_mean_wb_seq_post_mean; rot_dot_R_mean_wb_mean_now];
        rot_dot_dot_R_mean_wb_seq_post_mean = [rot_dot_dot_R_mean_wb_seq_post_mean; rot_dot_dot_R_mean_wb_mean_now];
        
    else
        % wb kin
        f_wb_seq_post_mean = [f_wb_seq_post_mean; nan];

        stroke_wb_L_seq_bins_post_mean = [stroke_wb_L_seq_bins_post_mean; nan(n_bins-1,1)];
        stroke_wb_R_seq_bins_post_mean = [stroke_wb_R_seq_bins_post_mean; nan(n_bins-1,1)];
        
        dev_wb_L_seq_bins_post_mean = [dev_wb_L_seq_bins_post_mean; nan(n_bins-1,1)];
        dev_wb_R_seq_bins_post_mean = [dev_wb_R_seq_bins_post_mean; nan(n_bins-1,1)];
        
        pitch_wb_L_seq_bins_post_mean = [pitch_wb_L_seq_bins_post_mean; nan(n_bins-1,1)];
        pitch_wb_R_seq_bins_post_mean = [pitch_wb_R_seq_bins_post_mean; nan(n_bins-1,1)];
        
        % body kin
        roll_mean_wb_seq_post_mean = [roll_mean_wb_seq_post_mean; nan];
        roll_dot_mean_wb_seq_post_mean = [roll_dot_mean_wb_seq_post_mean; nan];
        roll_dot_dot_mean_wb_seq_post_mean = [roll_dot_dot_mean_wb_seq_post_mean; nan];
        
        pitch_mean_wb_seq_post_mean = [pitch_mean_wb_seq_post_mean; nan];
        pitch_dot_mean_wb_seq_post_mean = [pitch_dot_mean_wb_seq_post_mean; nan];
        pitch_dot_dot_mean_wb_seq_post_mean = [pitch_dot_dot_mean_wb_seq_post_mean; nan];
        
        yaw_mean_wb_seq_post_mean = [yaw_mean_wb_seq_post_mean; nan];
        yaw_dot_mean_wb_seq_post_mean = [yaw_dot_mean_wb_seq_post_mean; nan];
        yaw_dot_dot_mean_wb_seq_post_mean = [yaw_dot_dot_mean_wb_seq_post_mean; nan];
        
        drot_L_mean_wb_seq_post_mean = [drot_L_mean_wb_seq_post_mean; nan];
        rot_dot_L_mean_wb_seq_post_mean = [rot_dot_L_mean_wb_seq_post_mean; nan];
        rot_dot_dot_L_mean_wb_seq_post_mean = [rot_dot_dot_L_mean_wb_seq_post_mean; nan];
        
        drot_R_mean_wb_seq_post_mean = [drot_R_mean_wb_seq_post_mean; nan];
        rot_dot_R_mean_wb_seq_post_mean = [rot_dot_R_mean_wb_seq_post_mean; nan];
        rot_dot_dot_R_mean_wb_seq_post_mean = [rot_dot_dot_R_mean_wb_seq_post_mean; nan];
    end
end

%% pre & post data
% wb kin
f_wb_seq_mean = [f_wb_seq_pre_mean;f_wb_seq_post_mean];

stroke_wb_L_seq_bins_mean = [stroke_wb_L_seq_bins_pre_mean;stroke_wb_L_seq_bins_post_mean];
stroke_wb_R_seq_bins_mean = [stroke_wb_R_seq_bins_pre_mean;stroke_wb_R_seq_bins_post_mean];

dev_wb_L_seq_bins_mean = [dev_wb_L_seq_bins_pre_mean;dev_wb_L_seq_bins_post_mean];
dev_wb_R_seq_bins_mean = [dev_wb_R_seq_bins_pre_mean;dev_wb_R_seq_bins_post_mean];

pitch_wb_L_seq_bins_mean = [pitch_wb_L_seq_bins_pre_mean;pitch_wb_L_seq_bins_post_mean];
pitch_wb_R_seq_bins_mean = [pitch_wb_R_seq_bins_pre_mean;pitch_wb_R_seq_bins_post_mean];

% body kin
roll_mean_wb_seq_mean = [roll_mean_wb_seq_pre_mean;roll_mean_wb_seq_post_mean];
roll_dot_mean_wb_seq_mean = [roll_dot_mean_wb_seq_pre_mean;roll_dot_mean_wb_seq_post_mean];
roll_dot_dot_mean_wb_seq_mean = [roll_dot_dot_mean_wb_seq_pre_mean;roll_dot_dot_mean_wb_seq_post_mean];

pitch_mean_wb_seq_mean = [pitch_mean_wb_seq_pre_mean;pitch_mean_wb_seq_post_mean];
pitch_dot_mean_wb_seq_mean = [pitch_dot_mean_wb_seq_pre_mean;pitch_dot_mean_wb_seq_post_mean];
pitch_dot_dot_mean_wb_seq_mean = [pitch_dot_dot_mean_wb_seq_pre_mean;pitch_dot_dot_mean_wb_seq_post_mean];

yaw_mean_wb_seq_mean = [yaw_mean_wb_seq_pre_mean;yaw_mean_wb_seq_post_mean];
yaw_dot_mean_wb_seq_mean = [yaw_dot_mean_wb_seq_pre_mean;yaw_dot_mean_wb_seq_post_mean];
yaw_dot_dot_mean_wb_seq_mean = [yaw_dot_dot_mean_wb_seq_pre_mean;yaw_dot_dot_mean_wb_seq_post_mean];

drot_L_mean_wb_seq_mean = [drot_L_mean_wb_seq_pre_mean;drot_L_mean_wb_seq_post_mean];
rot_dot_L_mean_wb_seq_mean = [rot_dot_L_mean_wb_seq_pre_mean;rot_dot_L_mean_wb_seq_post_mean];
rot_dot_dot_L_mean_wb_seq_mean = [rot_dot_dot_L_mean_wb_seq_pre_mean;rot_dot_dot_L_mean_wb_seq_post_mean];

drot_R_mean_wb_seq_mean = [drot_R_mean_wb_seq_pre_mean;drot_R_mean_wb_seq_post_mean];
rot_dot_R_mean_wb_seq_mean = [rot_dot_R_mean_wb_seq_pre_mean;rot_dot_R_mean_wb_seq_post_mean];
rot_dot_dot_R_mean_wb_seq_mean = [rot_dot_dot_R_mean_wb_seq_pre_mean;rot_dot_dot_R_mean_wb_seq_post_mean];

%% normalized time
% bin time
t_norm_wb_seq_bins_pre_mean = [0:length(stroke_wb_L_seq_bins_pre_mean)-1]'/(n_bins-1);
t_norm_wb_seq_bins_pre_mean = -t_norm_wb_seq_bins_pre_mean;
t_norm_wb_seq_bins_pre_mean = flipud(t_norm_wb_seq_bins_pre_mean);

t_norm_wb_seq_bins_post_mean = [0:length(stroke_wb_L_seq_bins_post_mean)-1]'/(n_bins-1);

t_norm_wb_seq_bins_mean = [t_norm_wb_seq_bins_pre_mean;t_norm_wb_seq_bins_post_mean];

% wb average time
t_norm_wb_seq_pre_mean = [1:wb_max]';
t_norm_wb_seq_pre_mean = -t_norm_wb_seq_pre_mean;
t_norm_wb_seq_pre_mean = flipud(t_norm_wb_seq_pre_mean);

t_norm_wb_seq_post_mean = [0:wb_max-1]';

t_norm_wb_seq_mean = [t_norm_wb_seq_pre_mean;t_norm_wb_seq_post_mean];

%% real time
t_norm = [1:n_bins]'/(n_bins-1);
t_norm = t_norm(1:end-1);

% pre
t_wb_seq_bins_pre_mean = [];
t_wb_seq_pre_mean = [];
for i = 1:length(f_wb_seq_pre_mean)
    if isnan(f_wb_seq_pre_mean(i)) == 0
        
        f_now = f_wb_seq_pre_mean(i);
        t_now = t_norm / f_now;
        
        if ~isempty(t_wb_seq_bins_pre_mean) && isnan(t_wb_seq_bins_pre_mean(end))==0
            t_now = t_now + t_wb_seq_bins_pre_mean(end);
        end
    else
        t_now = nan(length(t_norm),1);
    end
    t_wb_seq_bins_pre_mean = [t_wb_seq_bins_pre_mean;t_now];
    t_wb_seq_pre_mean = [t_wb_seq_pre_mean;t_now(1)];
end
t_wb_seq_pre_mean = t_wb_seq_pre_mean - max(t_wb_seq_bins_pre_mean);
t_wb_seq_bins_pre_mean = t_wb_seq_bins_pre_mean - max(t_wb_seq_bins_pre_mean);


% post
t_wb_seq_bins_post_mean = [];
t_wb_seq_post_mean = [];
for i = 1:length(f_wb_seq_post_mean)
    if isnan(f_wb_seq_post_mean(i)) == 0
        
        f_now = f_wb_seq_post_mean(i);
        t_now = t_norm / f_now;
        
        if ~isempty(t_wb_seq_bins_post_mean) && isnan(t_wb_seq_bins_post_mean(end))==0
            t_now = t_now + t_wb_seq_bins_post_mean(end);
        end
    else
        t_now = nan(length(t_norm),1);
    end
    t_wb_seq_bins_post_mean = [t_wb_seq_bins_post_mean;t_now];
    t_wb_seq_post_mean = [t_wb_seq_post_mean;t_now(1)];
end

t_wb_seq_bins_mean = [t_wb_seq_bins_pre_mean;t_wb_seq_bins_post_mean];
t_wb_seq_mean = [t_wb_seq_pre_mean;t_wb_seq_post_mean];

%% plot data
% wb kin
figure
hold on
plot(t_norm_wb_seq_pre_mean,f_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,f_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,f_wb_seq_mean,':k')

figure
hold on
plot(t_norm_wb_seq_bins_pre_mean,stroke_wb_L_seq_bins_pre_mean,'b')
plot(t_norm_wb_seq_bins_post_mean,stroke_wb_L_seq_bins_post_mean,'c')

plot(t_norm_wb_seq_bins_pre_mean,stroke_wb_R_seq_bins_pre_mean,'r')
plot(t_norm_wb_seq_bins_post_mean,stroke_wb_R_seq_bins_post_mean,'m')

figure
hold on
plot(t_norm_wb_seq_bins_mean,stroke_wb_L_seq_bins_mean,'b')
plot(t_norm_wb_seq_bins_mean,stroke_wb_R_seq_bins_mean,'r')

figure
hold on
plot(t_norm_wb_seq_bins_pre_mean,dev_wb_L_seq_bins_pre_mean,'b')
plot(t_norm_wb_seq_bins_post_mean,dev_wb_L_seq_bins_post_mean,'c')

plot(t_norm_wb_seq_bins_pre_mean,dev_wb_R_seq_bins_pre_mean,'r')
plot(t_norm_wb_seq_bins_post_mean,dev_wb_R_seq_bins_post_mean,'m')

figure
hold on
plot(t_norm_wb_seq_bins_mean,dev_wb_L_seq_bins_mean,'b')
plot(t_norm_wb_seq_bins_mean,dev_wb_R_seq_bins_mean,'r')

figure
hold on
plot(t_norm_wb_seq_bins_pre_mean,pitch_wb_L_seq_bins_pre_mean,'b')
plot(t_norm_wb_seq_bins_post_mean,pitch_wb_L_seq_bins_post_mean,'c')

plot(t_norm_wb_seq_bins_pre_mean,pitch_wb_R_seq_bins_pre_mean,'r')
plot(t_norm_wb_seq_bins_post_mean,pitch_wb_R_seq_bins_post_mean,'m')

figure
hold on
plot(t_norm_wb_seq_bins_mean,pitch_wb_L_seq_bins_mean,'b')
plot(t_norm_wb_seq_bins_mean,pitch_wb_R_seq_bins_mean,'r')

% body kin
figure
subplot(3,3,1)
hold on
plot(t_norm_wb_seq_pre_mean,roll_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,roll_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,roll_mean_wb_seq_mean,':k')

subplot(3,3,2)
hold on
plot(t_norm_wb_seq_pre_mean,roll_dot_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,roll_dot_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,roll_dot_mean_wb_seq_mean,':k')

subplot(3,3,3)
hold on
plot(t_norm_wb_seq_pre_mean,roll_dot_dot_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,roll_dot_dot_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,roll_dot_dot_mean_wb_seq_mean,':k')


subplot(3,3,4)
hold on
plot(t_norm_wb_seq_pre_mean,pitch_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,pitch_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,pitch_mean_wb_seq_mean,':k')

subplot(3,3,5)
hold on
plot(t_norm_wb_seq_pre_mean,pitch_dot_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,pitch_dot_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,pitch_dot_mean_wb_seq_mean,':k')

subplot(3,3,6)
hold on
plot(t_norm_wb_seq_pre_mean,pitch_dot_dot_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,pitch_dot_dot_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,pitch_dot_dot_mean_wb_seq_mean,':k')


subplot(3,3,7)
hold on
plot(t_norm_wb_seq_pre_mean,yaw_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,yaw_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,yaw_mean_wb_seq_mean,':k')

subplot(3,3,8)
hold on
plot(t_norm_wb_seq_pre_mean,yaw_dot_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,yaw_dot_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,yaw_dot_mean_wb_seq_mean,':k')

subplot(3,3,9)
hold on
plot(t_norm_wb_seq_pre_mean,yaw_dot_dot_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,yaw_dot_dot_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,yaw_dot_dot_mean_wb_seq_mean,':k')

figure
subplot(3,3,1)
hold on
plot(t_norm_wb_seq_pre_mean,drot_R_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,drot_R_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,drot_R_mean_wb_seq_mean,':k')

subplot(3,3,2)
hold on
plot(t_norm_wb_seq_pre_mean,rot_dot_R_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,rot_dot_R_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,rot_dot_R_mean_wb_seq_mean,':k')

subplot(3,3,3)
hold on
plot(t_norm_wb_seq_pre_mean,rot_dot_dot_R_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,rot_dot_dot_R_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,rot_dot_dot_R_mean_wb_seq_mean,':k')

subplot(3,3,4)
hold on
plot(t_norm_wb_seq_pre_mean,drot_L_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,drot_L_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,drot_L_mean_wb_seq_mean,':k')

subplot(3,3,5)
hold on
plot(t_norm_wb_seq_pre_mean,rot_dot_L_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,rot_dot_L_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,rot_dot_L_mean_wb_seq_mean,':k')

subplot(3,3,6)
hold on
plot(t_norm_wb_seq_pre_mean,rot_dot_dot_L_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,rot_dot_dot_L_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,rot_dot_dot_L_mean_wb_seq_mean,':k')

subplot(3,3,7)
hold on
plot(t_norm_wb_seq_pre_mean,yaw_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,yaw_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,yaw_mean_wb_seq_mean,':k')

subplot(3,3,8)
hold on
plot(t_norm_wb_seq_pre_mean,yaw_dot_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,yaw_dot_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,yaw_dot_mean_wb_seq_mean,':k')

subplot(3,3,9)
hold on
plot(t_norm_wb_seq_pre_mean,yaw_dot_dot_mean_wb_seq_pre_mean,'r')
plot(t_norm_wb_seq_post_mean,yaw_dot_dot_mean_wb_seq_post_mean,'c')
plot(t_norm_wb_seq_mean,yaw_dot_dot_mean_wb_seq_mean,':k')

%% plot data
% wb kin
figure
hold on
plot(t_wb_seq_pre_mean,f_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,f_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,f_wb_seq_mean,':k')

figure
hold on
plot(t_wb_seq_bins_pre_mean,stroke_wb_L_seq_bins_pre_mean,'b')
plot(t_wb_seq_bins_post_mean,stroke_wb_L_seq_bins_post_mean,'c')

plot(t_wb_seq_bins_pre_mean,stroke_wb_R_seq_bins_pre_mean,'r')
plot(t_wb_seq_bins_post_mean,stroke_wb_R_seq_bins_post_mean,'m')

figure
hold on
plot(t_wb_seq_bins_mean,stroke_wb_L_seq_bins_mean,'b')
plot(t_wb_seq_bins_mean,stroke_wb_R_seq_bins_mean,'r')

figure
hold on
plot(t_wb_seq_bins_pre_mean,dev_wb_L_seq_bins_pre_mean,'b')
plot(t_wb_seq_bins_post_mean,dev_wb_L_seq_bins_post_mean,'c')

plot(t_wb_seq_bins_pre_mean,dev_wb_R_seq_bins_pre_mean,'r')
plot(t_wb_seq_bins_post_mean,dev_wb_R_seq_bins_post_mean,'m')

figure
hold on
plot(t_wb_seq_bins_mean,dev_wb_L_seq_bins_mean,'b')
plot(t_wb_seq_bins_mean,dev_wb_R_seq_bins_mean,'r')

figure
hold on
plot(t_wb_seq_bins_pre_mean,pitch_wb_L_seq_bins_pre_mean,'b')
plot(t_wb_seq_bins_post_mean,pitch_wb_L_seq_bins_post_mean,'c')

plot(t_wb_seq_bins_pre_mean,pitch_wb_R_seq_bins_pre_mean,'r')
plot(t_wb_seq_bins_post_mean,pitch_wb_R_seq_bins_post_mean,'m')

figure
hold on
plot(t_wb_seq_bins_mean,pitch_wb_L_seq_bins_mean,'b')
plot(t_wb_seq_bins_mean,pitch_wb_R_seq_bins_mean,'r')

% body kin
figure
subplot(3,3,1)
hold on
plot(t_wb_seq_pre_mean,roll_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,roll_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,roll_mean_wb_seq_mean,':k')

subplot(3,3,2)
hold on
plot(t_wb_seq_pre_mean,roll_dot_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,roll_dot_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,roll_dot_mean_wb_seq_mean,':k')

subplot(3,3,3)
hold on
plot(t_wb_seq_pre_mean,roll_dot_dot_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,roll_dot_dot_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,roll_dot_dot_mean_wb_seq_mean,':k')


subplot(3,3,4)
hold on
plot(t_wb_seq_pre_mean,pitch_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,pitch_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,pitch_mean_wb_seq_mean,':k')

subplot(3,3,5)
hold on
plot(t_wb_seq_pre_mean,pitch_dot_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,pitch_dot_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,pitch_dot_mean_wb_seq_mean,':k')

subplot(3,3,6)
hold on
plot(t_wb_seq_pre_mean,pitch_dot_dot_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,pitch_dot_dot_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,pitch_dot_dot_mean_wb_seq_mean,':k')


subplot(3,3,7)
hold on
plot(t_wb_seq_pre_mean,yaw_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,yaw_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,yaw_mean_wb_seq_mean,':k')

subplot(3,3,8)
hold on
plot(t_wb_seq_pre_mean,yaw_dot_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,yaw_dot_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,yaw_dot_mean_wb_seq_mean,':k')

subplot(3,3,9)
hold on
plot(t_wb_seq_pre_mean,yaw_dot_dot_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,yaw_dot_dot_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,yaw_dot_dot_mean_wb_seq_mean,':k')

figure
subplot(3,3,1)
hold on
plot(t_wb_seq_pre_mean,drot_R_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,drot_R_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,drot_R_mean_wb_seq_mean,':k')

subplot(3,3,2)
hold on
plot(t_wb_seq_pre_mean,rot_dot_R_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,rot_dot_R_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,rot_dot_R_mean_wb_seq_mean,':k')

subplot(3,3,3)
hold on
plot(t_wb_seq_pre_mean,rot_dot_dot_R_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,rot_dot_dot_R_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,rot_dot_dot_R_mean_wb_seq_mean,':k')

subplot(3,3,4)
hold on
plot(t_wb_seq_pre_mean,drot_L_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,drot_L_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,drot_L_mean_wb_seq_mean,':k')

subplot(3,3,5)
hold on
plot(t_wb_seq_pre_mean,rot_dot_L_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,rot_dot_L_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,rot_dot_L_mean_wb_seq_mean,':k')

subplot(3,3,6)
hold on
plot(t_wb_seq_pre_mean,rot_dot_dot_L_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,rot_dot_dot_L_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,rot_dot_dot_L_mean_wb_seq_mean,':k')

subplot(3,3,7)
hold on
plot(t_wb_seq_pre_mean,yaw_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,yaw_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,yaw_mean_wb_seq_mean,':k')

subplot(3,3,8)
hold on
plot(t_wb_seq_pre_mean,yaw_dot_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,yaw_dot_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,yaw_dot_mean_wb_seq_mean,':k')

subplot(3,3,9)
hold on
plot(t_wb_seq_pre_mean,yaw_dot_dot_mean_wb_seq_pre_mean,'r')
plot(t_wb_seq_post_mean,yaw_dot_dot_mean_wb_seq_post_mean,'c')
plot(t_wb_seq_mean,yaw_dot_dot_mean_wb_seq_mean,':k')

