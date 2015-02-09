clc
clear
close all

% % LOOMING
%     loadname = 'kinflightpathDB_pos_qbodyEKF_INCroll_2clusters_Ahor2.75mps2_strokeplane47.5deg_startframe2945.mat';

% % SACCADE (startframe 100, kalman filter issues at start of seq)
% loadname='kinflightpathDB_pos_qbodyEKF_INCroll_9clusters_1.4n-1.7n1.9_strokeplane47.5deg_startframe100.mat'

% CLIP
loadname = 'kinflightpathDB_pos_qbodyEKF_INCroll_strokeplane47.5deg_rotaxis0deg_startframe1.mat';

load(loadname)

const_name=dir('flyVars*')
const_name=const_name.name;
load(const_name)

%% settings
n=200; % bins
fr_max = 50;

fps = settings.fps;
cmap_180 = jet(180);

%% data import
steady_frames = pathDB.steady_frames;
strokeplane_angle = strokeplane_WBkin;

cmap_k = settings.cmap_k;
cmap_360 = settings.cmap_360;

t = pathDB.t;
dt = t(2)-t(1);
t_all = t;
for i=1:size(steady_frames,2)-1
    t_all = [t_all t];
end

%% clip data
clip_side = clipDB.clip_side;
clip_type = clipDB.clip_type;

SecondMomentClipped = clipDB.SecondMomentCC;
SecondMomentIntact = clipDB.SecondMomentCI;
SecondMomentRatio = clipDB.SecondMomentRatioCentroid;

AreaClipped = clipDB.WingAreaC;
AreaIntact = clipDB.WingAreaI;
AreaRatio = clipDB.WingAreaRatio;

LengthClipped = clipDB.WingLengthC;
LengthIntact = clipDB.WingLengthI;
LengthRatio = clipDB.WingLengthRatio;

%% wingkin data
n_down_start_L = kinDB.n_down_start_L;
n_up_start_L = kinDB.n_up_start_L;
wingtip_path_L = kinDB.wingtip_path_L;
joint_pos_L = kinDB.joint_pos_L;
stroke_L = kinDB.stroke_L;
pitch_L = kinDB.pitch_L;
dev_L = kinDB.dev_L;

n_down_start_R = kinDB.n_down_start_R;
n_up_start_R = kinDB.n_up_start_R;
wingtip_path_R = kinDB.wingtip_path_R;
joint_pos_R = kinDB.joint_pos_R;
stroke_R = kinDB.stroke_R;
pitch_R = kinDB.pitch_R;
dev_R = kinDB.dev_R;

wing_length = kinDB.wing_length;

% U&aoa@75% wingspan
U_L(:,:) = kinDB.Uwing_L(:,3,:);
U_R(:,:) = kinDB.Uwing_R(:,3,:);
aoa_L(:,:) = kinDB.aoa_L(:,3,:);
aoa_R(:,:) = kinDB.aoa_R(:,3,:);

%% reverse wingkin data with left wing clipped: LEFT WING = INTACT; RIGHT WING = CLIPPED
stroke_L_temp = stroke_L;
stroke_R_temp = stroke_R;
pitch_L_temp = pitch_L;
pitch_R_temp = pitch_R;
dev_L_temp = dev_L;
dev_R_temp = dev_R;
U_L_temp = U_L;
U_R_temp = U_R;
aoa_L_temp = aoa_L;
aoa_R_temp = aoa_R;

for i = 1:length(clip_side)
    if clip_side(i) == 1
        stroke_L(:,i) = stroke_R_temp(:,i);
        stroke_R(:,i) = stroke_L_temp(:,i);
        pitch_L(:,i) = pitch_R_temp(:,i);
        pitch_R(:,i) = pitch_L_temp(:,i);
        dev_L(:,i) = dev_R_temp(:,i);
        dev_R(:,i) = dev_L_temp(:,i);
        U_L(:,i) = U_R_temp(:,i);
        U_R(:,i) = U_L_temp(:,i);
        aoa_L(:,i) = aoa_R_temp(:,i);
        aoa_R(:,i) = aoa_L_temp(:,i);
    end
end

%% calc & save wb dataset
% calc_wb_data_subsets_INCrotAxes_InCnorm
calc_wb_data_subsets_INCrotAxes1n2_InCnoNorm

