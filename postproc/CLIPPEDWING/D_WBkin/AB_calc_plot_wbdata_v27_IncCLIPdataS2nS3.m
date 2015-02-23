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
clip_side_all = clipDB.clip_side;
clip_type_all = clipDB.clip_type;

SecondMomentClipped_all = clipDB.SecondMomentCC;
SecondMomentIntact_all = clipDB.SecondMomentCI;
SecondMomentRatio_all = clipDB.SecondMomentRatioCentroid;

ThirdMomentClipped_all = clipDB.ThirdMomentCC;
ThirdMomentIntact_all = clipDB.ThirdMomentCI;
ThirdMomentRatio_all = clipDB.ThirdMomentRatioCentroid;

AreaClipped_all = clipDB.WingAreaC;
AreaIntact_all = clipDB.WingAreaI;
AreaRatio_all = clipDB.WingAreaRatio;

LengthClipped_all = clipDB.WingLengthC;
LengthIntact_all = clipDB.WingLengthI;
LengthRatio_all = clipDB.WingLengthRatio;

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

%% BODY KIN DATA
V = pathDB.V;

slip = -pathDB.slip;                        % !!! negative slip !!!
pitch = pathDB.pitch;
roll = pathDB.roll;

slip_body = -pathDB.slip_body;              % !!! negative slip !!!
pitch_body = pathDB.pitch_body;
roll_body = pathDB.roll;

dyaw = pathDB.dyaw_sp;
dpitch = pathDB.dpitch_sp;
droll = pathDB.droll_sp;

dyaw_body = pathDB.dyaw_body;
dpitch_body = pathDB.dpitch_body;
droll_body = pathDB.droll_body;

slip_aero = pathDB.slip_aero;
pitch_aero = pathDB.pitch_aero;
slip_aero_hor = pathDB.slip_aero_hor;
pitch_aero_hor = pathDB.pitch_aero_hor;

Fsp_pitch = pathDB.Fsp_pitch;
Fsp_roll = pathDB.Fsp_roll;

Fb_pitch = pathDB.Fb_pitch;
Fb_roll = pathDB.Fb_roll;

%% MIRROR wing&body kin data with left wing clipped: LEFT WING = INTACT; RIGHT WING = CLIPPED
% wingbeat kin
stroke_L_mirror = stroke_L;
stroke_R_mirror = stroke_R;
pitch_L_mirror = pitch_L;
pitch_R_mirror = pitch_R;
dev_L_mirror = dev_L;
dev_R_mirror = dev_R;
U_L_mirror = U_L;
U_R_mirror = U_R;
aoa_L_mirror = aoa_L;
aoa_R_mirror = aoa_R;

% body kin
slip_mirror = slip;
roll_mirror = roll;
slip_body_mirror = slip_body;
roll_body_mirror = roll_body;
dyaw_mirror = dyaw;
droll_mirror = droll;
dyaw_body_mirror = dyaw_body;
droll_body_mirror = droll_body;
slip_aero_mirror = slip_aero;
slip_aero_hor_mirror = slip_aero_hor;
Fsp_roll_mirror = Fsp_roll;
Fb_roll_mirror = Fb_roll;


for i=1:length(clip_side_all)
    if clip_side_all(i) == 1
        
        % wb kin
        stroke_L_mirror(:,i) = stroke_R(:,i);
        stroke_R_mirror(:,i) = stroke_L(:,i);
        pitch_L_mirror(:,i) = pitch_R(:,i);
        pitch_R_mirror(:,i) = pitch_L(:,i);
        dev_L_mirror(:,i) = dev_R(:,i);
        dev_R_mirror(:,i) = dev_L(:,i);
        U_L_mirror(:,i) = U_R(:,i);
        U_R_mirror(:,i) = U_L(:,i);
        aoa_L_mirror(:,i) = aoa_R(:,i);
        aoa_R_mirror(:,i) = aoa_L(:,i);

        % body kin
        slip_mirror(:,i) = -slip_mirror(:,i);
        roll_mirror(:,i) = -roll_mirror(:,i);
        slip_body_mirror(:,i) = -slip_body_mirror(:,i);
        roll_body_mirror(:,i) = -roll_body_mirror(:,i);
        dyaw_mirror(:,i) = -dyaw_mirror(:,i);
        droll_mirror(:,i) = -droll_mirror(:,i);
        dyaw_body_mirror(:,i) = -dyaw_body_mirror(:,i);
        droll_body_mirror(:,i) = -droll_body_mirror(:,i);
        slip_aero_mirror(:,i) = -slip_aero_mirror(:,i);
        slip_aero_hor_mirror(:,i) = -slip_aero_hor_mirror(:,i);
        Fsp_roll_mirror(:,i) = -Fsp_roll_mirror(:,i);
        Fb_roll_mirror(:,i) = -Fb_roll_mirror(:,i);
    end
end

%% calc & save wb dataset
% calc_wb_data_subsets_INCrotAxes_InCnorm
% calc_wb_data_subsets_CLIP
% calc_wb_data_subsets_incCLIPdata
calc_wb_data_subsets_incCLIPdataS2nS3

