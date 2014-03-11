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
wb_min4mean = 10;


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
    
    % pre response wbs
    wb_pre_max = length(frame_nr_2resp_now_pre);
    for wb_now = 1:wb_pre_max
        
        wb_pre = wb_pre_max - wb_now +1;
        
        f_wb_seq_pre(seq_now,wb_pre) = f_wb_now_pre(wb_now);
        
        stroke_wb_L_seq_bins_pre(:,seq_now,wb_pre) = stroke_wb_L_bins_now_pre(:,wb_now);
        stroke_wb_R_seq_bins_pre(:,seq_now,wb_pre) = stroke_wb_R_bins_now_pre(:,wb_now);
        
        dev_wb_L_seq_bins_pre(:,seq_now,wb_pre) = dev_wb_L_bins_now_pre(:,wb_now);
        dev_wb_R_seq_bins_pre(:,seq_now,wb_pre) = dev_wb_R_bins_now_pre(:,wb_now);
        
        pitch_wb_L_seq_bins_pre(:,seq_now,wb_pre) = pitch_wb_L_bins_now_pre(:,wb_now);
        pitch_wb_R_seq_bins_pre(:,seq_now,wb_pre) = pitch_wb_R_bins_now_pre(:,wb_now);
    end
    
    % post response wbs
    wb_post_max = length(frame_nr_2resp_now_post);
    for wb_post = 1:wb_post_max
        
        f_wb_seq_post(seq_now,wb_post) = f_wb_now_post(wb_post);
        
        stroke_wb_L_seq_bins_post(:,seq_now,wb_post) = stroke_wb_L_bins_now_post(:,wb_post);
        stroke_wb_R_seq_bins_post(:,seq_now,wb_post) = stroke_wb_R_bins_now_post(:,wb_post);
        
        dev_wb_L_seq_bins_post(:,seq_now,wb_post) = dev_wb_L_bins_now_post(:,wb_post);
        dev_wb_R_seq_bins_post(:,seq_now,wb_post) = dev_wb_R_bins_now_post(:,wb_post);
        
        pitch_wb_L_seq_bins_post(:,seq_now,wb_post) = pitch_wb_L_bins_now_post(:,wb_post);
        pitch_wb_R_seq_bins_post(:,seq_now,wb_post) = pitch_wb_R_bins_now_post(:,wb_post);
    end
end

%% mean wingbeat seq throughout maneuver
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

for wb_now = 1:wb_max
    wb_pre = wb_max - wb_now +1;
    wb_post = wb_now;
    
    %% only mean if Nwb>wb_min4mean
    % pre
    if sum(isnan(stroke_wb_L_seq_bins_pre(1,:,wb_pre))==0) >= wb_min4mean
        stroke_wb_L_seq_bins_pre_mean = [stroke_wb_L_seq_bins_pre_mean; nanmean(stroke_wb_L_seq_bins_pre(:,:,wb_pre)')'];
        stroke_wb_R_seq_bins_pre_mean = [stroke_wb_R_seq_bins_pre_mean; nanmean(stroke_wb_R_seq_bins_pre(:,:,wb_pre)')'];
        
        dev_wb_L_seq_bins_pre_mean = [dev_wb_L_seq_bins_pre_mean; nanmean(dev_wb_L_seq_bins_pre(:,:,wb_pre)')'];
        dev_wb_R_seq_bins_pre_mean = [dev_wb_R_seq_bins_pre_mean; nanmean(dev_wb_R_seq_bins_pre(:,:,wb_pre)')'];
        
        pitch_wb_L_seq_bins_pre_mean = [pitch_wb_L_seq_bins_pre_mean; nanmean(pitch_wb_L_seq_bins_pre(:,:,wb_pre)')'];
        pitch_wb_R_seq_bins_pre_mean = [pitch_wb_R_seq_bins_pre_mean; nanmean(pitch_wb_R_seq_bins_pre(:,:,wb_pre)')'];
        
    else
        stroke_wb_L_seq_bins_pre_mean = [stroke_wb_L_seq_bins_pre_mean; nan(n_bins,1)];
        stroke_wb_R_seq_bins_pre_mean = [stroke_wb_R_seq_bins_pre_mean; nan(n_bins,1)];
        
        dev_wb_L_seq_bins_pre_mean = [dev_wb_L_seq_bins_pre_mean; nan(n_bins,1)];
        dev_wb_R_seq_bins_pre_mean = [dev_wb_R_seq_bins_pre_mean; nan(n_bins,1)];
        
        pitch_wb_L_seq_bins_pre_mean = [pitch_wb_L_seq_bins_pre_mean; nan(n_bins,1)];
        pitch_wb_R_seq_bins_pre_mean = [pitch_wb_R_seq_bins_pre_mean; nan(n_bins,1)];
        
    end
    
    % post
    if sum(isnan(stroke_wb_L_seq_bins_post(1,:,wb_post))==0) >= wb_min4mean
        stroke_wb_L_seq_bins_post_mean = [stroke_wb_L_seq_bins_post_mean; nanmean(stroke_wb_L_seq_bins_post(:,:,wb_post)')'];
        stroke_wb_R_seq_bins_post_mean = [stroke_wb_R_seq_bins_post_mean; nanmean(stroke_wb_R_seq_bins_post(:,:,wb_post)')'];
        
        dev_wb_L_seq_bins_post_mean = [dev_wb_L_seq_bins_post_mean; nanmean(dev_wb_L_seq_bins_post(:,:,wb_post)')'];
        dev_wb_R_seq_bins_post_mean = [dev_wb_R_seq_bins_post_mean; nanmean(dev_wb_R_seq_bins_post(:,:,wb_post)')'];
        
        pitch_wb_L_seq_bins_post_mean = [pitch_wb_L_seq_bins_post_mean; nanmean(pitch_wb_L_seq_bins_post(:,:,wb_post)')'];
        pitch_wb_R_seq_bins_post_mean = [pitch_wb_R_seq_bins_post_mean; nanmean(pitch_wb_R_seq_bins_post(:,:,wb_post)')'];
        
    else
        stroke_wb_L_seq_bins_post_mean = [stroke_wb_L_seq_bins_post_mean; nan(n_bins,1)];
        stroke_wb_R_seq_bins_post_mean = [stroke_wb_R_seq_bins_post_mean; nan(n_bins,1)];
        
        dev_wb_L_seq_bins_post_mean = [dev_wb_L_seq_bins_post_mean; nan(n_bins,1)];
        dev_wb_R_seq_bins_post_mean = [dev_wb_R_seq_bins_post_mean; nan(n_bins,1)];
        
        pitch_wb_L_seq_bins_post_mean = [pitch_wb_L_seq_bins_post_mean; nan(n_bins,1)];
        pitch_wb_R_seq_bins_post_mean = [pitch_wb_R_seq_bins_post_mean; nan(n_bins,1)];
        
    end
end

%% pre & post data
stroke_wb_L_seq_bins_mean = [stroke_wb_L_seq_bins_pre_mean;stroke_wb_L_seq_bins_post_mean(2:end)];
stroke_wb_R_seq_bins_mean = [stroke_wb_R_seq_bins_pre_mean;stroke_wb_R_seq_bins_post_mean(2:end)];

dev_wb_L_seq_bins_mean = [dev_wb_L_seq_bins_pre_mean;dev_wb_L_seq_bins_post_mean(2:end)];
dev_wb_R_seq_bins_mean = [dev_wb_R_seq_bins_pre_mean;dev_wb_R_seq_bins_post_mean(2:end)];

pitch_wb_L_seq_bins_mean = [pitch_wb_L_seq_bins_pre_mean;pitch_wb_L_seq_bins_post_mean(2:end)];
pitch_wb_R_seq_bins_mean = [pitch_wb_R_seq_bins_pre_mean;pitch_wb_R_seq_bins_post_mean(2:end)];

%% normalized time

t_wb_seq_bins_pre_mean = [0:length(stroke_wb_L_seq_bins_pre_mean)-1]'/n_bins;
t_wb_seq_bins_pre_mean = -t_wb_seq_bins_pre_mean;
t_wb_seq_bins_pre_mean = flipud(t_wb_seq_bins_pre_mean);

t_wb_seq_bins_post_mean = [0:length(stroke_wb_L_seq_bins_post_mean)-1]'/n_bins;

t_wb_seq_bins_mean = [t_wb_seq_bins_pre_mean;t_wb_seq_bins_post_mean(2:end)];


%% plot data
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















