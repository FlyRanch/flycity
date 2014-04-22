%% Run variables
date_str = '20140410'
seq_number = 1
startframe = 9;
BGframe = 2;
endframe = 1365;

%% Deployment variables
software_root_dir = '/home/matt/flycity/flytracker_tethered'
data_root_dir = '/media/thad01b/'

%% Kine settings
kine_dir = fullfile(software_root_dir,'kine');
setup_path = fullfile(software_root_dir,'kine','setup');
kine_mfile_path = kine_dir;
setup_file = 'flightmuscles.mat';

%% Flytracker settings
% PAR is a static structure that defines various parameters of the tracking
% sequence
global PAR
% String formating conventions
seq_format_str = 'C%03dH%03dS%04d' %cam number, head number, sequence number
sol_format_str = '%s_S%04d'%date_str, sol number
% File path formating
MoveModelPath = fullfile(software_root_dir,'Models',filesep);
sol_root_path = fullfile(data_root_dir,'solutions');
seq_sol_path = fullfile(sol_root_path,sprintf(sol_format_str,date_str,seq_number),filesep)
FileFromKine = fullfile(seq_sol_path,['kine_',sprintf(sol_format_str,date_str,seq_number),'.mat']);
cam1_seq_path = fullfile(data_root_dir,'Photron','SEQS',date_str,sprintf(seq_format_str,1,1,1),filesep)
cam1_files = dir(fullfile(cam1_seq_path,'*tif'));
cam1_file1 = cam1_files(1).name
cih_file = dir(fullfile(cam1_seq_path,'*.cih'));
cih_file = cih_file.name;

getIC = true;
getBodyShape = false;

% to plot or not to plot
plot_data = true;

% This is the number of sample points used within the Fly model.
%It changes how fine the mesh is. 
L1 = 15; %# of steps for body along length
L2 = 6; %# of steps for head along length
L3 = 30; %# of steps for wing around the boundary
T1 = 13; %# of theta steps for head and body
T2 = 2; %# of steps towards center of wing
% -- These dimensions have to be defined by the user.
%The dimensions are for a single fly.
%The subsample scale of the image is 1/2^(pwr)
etamax = 0;
paramdim = 1;
%number of flies
flynum = [1 ];  %[1 2];
numfly = length(flynum);
OccludeShape = {[],[],[]};
%spline order
c = 4; 
statedim = 15*ones(1,numfly); 
pNoisedim = statedim;
% This is the number of frames to scan back when performing the pattern
% matching in the motion prediction.
NumPrevSol = 5;
% Dont know what this is?
streams = 2;
% Tracker search distance (moved from Closest_PtsNrmlGate1, FTM 20120605)
WingSearchDist = [120 120];
BodySearchDist = [5 5];
numfeatpts = 30;
% Movement Model scale (VARIABLE SCALE DOES NOT WORK YET)
MoveModelScale = 1
BWfilterRatio = 0.3; % reduce min threshold
WingBodyRatio = 0.7; % increase wing-body treshold
% angular variance for (initial fit?) 2*standard deviation in degrees;
angvar = 0.002;
%Variance for the wing joint locations
JointVar = [0.001 0.0004 0.0003];
%Variance for the body linear acceleration
LinVar = [0.183 0.639 1.33];
%Variance for the body angular acceleration
AngVar = [.168 .181 2.20];
%Settings for the movement models
move_model_groups = containers.Map();
move_model_groups('computer_track_tethered_hydei_7500fps') = {'20130412_S0011_successful_flytrack_MoveModel.mat',
                                                              '20130412_S0017_successful_flytrack_MoveModel.mat',
                                                              '20130412_S0022_successful_flytrack_MoveModel.mat',
                                                              '20130412_S0023_successful_flytrack_MoveModel.mat',
                                                              '20130412_S0026_successful_flytrack_MoveModel.mat',
                                                              '20130412_S0028_successful_flytrack_MoveModel.mat',
                                                              '20130412_S0030_successful_flytrack_MoveModel.mat',
                                                              '20130412_S0034_successful_flytrack_MoveModel.mat',
                                                              '20130412_S0036_successful_flytrack_MoveModel.mat',
                                                              '20130412_S0051_successful_flytrack_MoveModel.mat',
                                                              '20130412_S0052_successful_flytrack_MoveModel.mat',
                                                              '20130412_S0054_successful_flytrack_MoveModel.mat',};
move_model_groups('manual_track_tethered_hydei_7500fps') = {'MoveModel_IncSearch_20130412_S0002_frame0001-0036_NOstimulus_onewingbeat.mat';
                                                              'MoveModel_IncSearch_20130412_S0017_frame0026-0037_NOstimulus_strokerev.mat';
                                                              'MoveModel_IncSearch_20130412_S0027_frame0760-0774_90deg_strokerev.mat';
                                                              'MoveModel_NOSearch_20130412_S0002_frame0001-0036_NOstimulus_onewingbeat.mat';
                                                              'MoveModel_NOSearch_20130412_S0017_frame0026-0037_NOstimulus_strokerev.mat';
                                                              'MoveModel_NOSearch_20130412_S0027_frame0760-0774_90deg_strokerev.mat';}
move_model_groups('manual_track_tethered_melanogaster_6000fps') = {'MoveModel_IncSearch_20140310_S0001.mat',
                                                                   'MoveModel_IncSearch_20140327_S0001.mat',
                                                                   'MoveModel_IncSearch_20140327_S0005.mat',
                                                                   'MoveModel_IncSearch_20140401_S0001.mat'};
move_model_groups('computer_track_tethered_melanogaster_6000fps') = {'20140310_S0001_successful_flytrack_MoveModel_frm599_615.mat',
                                                                     '20140310_S0001_successful_flytrack_MoveModel_frm411_448.mat',
                                                                     '20140310_S0001_successful_flytrack_MoveModel_frm246_325.mat',
                                                                     '20140310_S0001_successful_flytrack_MoveModel_frm141_203.mat',
                                                                     '20140327_S0001_successful_flytrack_MoveModel_frm18_682.mat'};

MoveModelNameS = move_model_groups('manual_track_tethered_melanogaster_6000fps');



% addpath('/home/flycity/Dropbox/WORK/flytracker_thad');
% addpath('/home/flycity/Dropbox/WORK/flytracker_thad/mex/');
% addpath('/home/flycity/Dropbox/WORK/flytracker_thad/core/');
% addpath('/home/flycity/Dropbox/WORK/flytracker_thad/results/');
% addpath('/home/flycity/Dropbox/WORK/postproc/AA_check_wbkin/')
