% DEMSE_MODAL  Demonstrate state estimation of 3D drosophila
% kinematic model 
%
%
%   See also in this script:
%   auto_init - initializes the tracker by calculating calibration,
%   background, and model shape
%   
%   gssm_flyOcc - defines the motion and measurement model of the fly
%
%
%   srcdkf - performs the estimation
% 
%   More details are available in each file
%   -------------------------------------------------------------------
%
%   When 'getIC' is set to 'true', this program will calculate and save a
%   structure named 'ManualFit' that contains all the tracking parameters
%   associated with a particular video sequence.  It saves it in the same
%   directory that the image data is located.
%
%   When 'getIC' is set to 'false', you have already calculated the initial
%   parameters, so the program just loads them from disk.
%=============================================================================================

close all
clear
clc

start = clock;
%startframe = 681;
%ICframe = startframe;

BGframe = 2;


% % Determine endframe as last flytrack
% cd flytracks
% frame_number = dir('fly*.mat');
% A = frame_number.name;
% startframe = str2num(A(8:9));
% endframe = length(frame_number);
% cd ..



%% default settings
% start sequence for first time (fit model to KINEfit and do bodyscaling)
getIC = true;
getBodyShape = true;

% continue after tracker has lost fly (fit model to new KINEfit but no bodyscaling)
%  getIC = true;
%  getBodyShape = false;

% % continue after tracker has crashed (reload previous frame and model)
% getIC = false;
% getBodyShape = false;

%% to plot or not to plot
plot_data = false;
%plot_data = true;


fprintf('\nDEMSE_MODAL : This demonstration performs state estimation\n');
fprintf('              for 3D Drosophila body and wing trajectories\n'); 
fprintf('              The pixel observations are corrupted by additive\n');
fprintf('              white Gaussian noise.\n\n');

% PAR is a static structure that defines various parameters of the tracking
% sequence
global PAR

%% Add the subdirectories that contain important / necessary files
% flytracker program files
warning off

% cd '/home/matt/Dropbox/WORK/flytracker/flytracker'
%addpath(pwd);
% addpath('/home/matt/Dropbox/WORK/flytracker_thad');
% addpath('/home/matt/Dropbox/WORK/flytracker_thad/mex/');
% addpath('/home/matt/Dropbox/WORK/flytracker_thad/core/');
% addpath('/home/matt/Dropbox/WORK/flytracker_thad/results/');
% addpath('/home/matt/Dropbox/WORK/postproc/AA_check_wbkin/');
% addpath('/home/matt/Dropbox/WORK/flytracker_thad/SpinConv');
% 
% addpath('/home/florian/Dropbox/FlyCity/WORK/flytracker_thad');
% addpath('/home/florian/Dropbox/FlyCity/WORK/flytracker_thad/mex/');
% addpath('/home/florian/Dropbox/FlyCity/WORK/flytracker_thad/core/');
% addpath('/home/florian/Dropbox/FlyCity/WORK/flytracker_thad/results/');
% addpath('/home/florian/Dropbox/FlyCity/WORK/postproc/AA_check_wbkin/')
% 
% addpath('/home/flycity/Dropbox/WORK/flytracker_thad');
% addpath('/home/flycity/Dropbox/WORK/flytracker_thad/mex/');
% addpath('/home/flycity/Dropbox/WORK/flytracker_thad/core/');
% addpath('/home/flycity/Dropbox/WORK/flytracker_thad/results/');
% addpath('/home/flycity/Dropbox/WORK/postproc/AA_check_wbkin/')
% 
% addpath('/home/nicole/Dropbox/WORK/flytracker_thad');
% addpath('/home/nicole/Dropbox/WORK/flytracker_thad/mex/');
% addpath('/home/nicole/Dropbox/WORK/flytracker_thad/core/');
% addpath('/home/nicole/Dropbox/WORK/flytracker_thad/results/');
% addpath('/home/nicole/Dropbox/WORK/postproc/AA_check_wbkin/');

%% KINE file

%KinePath = '/media/DATA1/Video_pilot_1/solutions/';
KinePath = '/media/Matt01/solutions/';

% manual select
KineName = [];

%% Movement Model data

% MoveModelPath = '/home/flycity/Dropbox/Models/matt_5.2013/';
if isdir('/home/matt/Dropbox/WORK/flytracker_thad') == 1
    MoveModelPath = '/home/matt/Dropbox/Models/';
end
if isdir('/home/flycity/Dropbox/WORK/flytracker_thad') == 1;
    MoveModelPath = '/home/flycity/Dropbox/Models/';
end
if isdir('/home/nicole/Dropbox/WORK/flytracker_thad') == 1;
    MoveModelPath = '/home/nicole/Dropbox/Models/';
end
if isdir('/home/florian/Dropbox/FlyCity/WORK/flytracker_thad') == 1
    MoveModelPath = '/home/florian/Dropbox/FlyCity/Models/';
end

%% florians models
% MoveModelName1 = '20121010_S0001_successful_flytrack_MoveModel.mat';
% MoveModelName2 = '20121009_S0003_successful_flytrack_MoveModel.mat';
% MoveModelName3 = '20121002_S0003_successful_flytrack_MoveModel.mat';
% MoveModelNameS = {MoveModelName1 MoveModelName2 MoveModelName3}';

% Successfully tracked tethered flight models - Hydei 7500 fps
MoveModelName1 = '20130412_S0011_successful_flytrack_MoveModel.mat';
MoveModelName2 = '20130412_S0017_successful_flytrack_MoveModel.mat';
MoveModelName3 = '20130412_S0022_successful_flytrack_MoveModel.mat';
MoveModelName4 = '20130412_S0023_successful_flytrack_MoveModel.mat';
MoveModelName5 = '20130412_S0026_successful_flytrack_MoveModel.mat';
MoveModelName6 = '20130412_S0028_successful_flytrack_MoveModel.mat';
MoveModelName7 = '20130412_S0030_successful_flytrack_MoveModel.mat';
MoveModelName8 = '20130412_S0034_successful_flytrack_MoveModel.mat';
MoveModelName9 = '20130412_S0036_successful_flytrack_MoveModel.mat';
MoveModelName10 = '20130412_S0051_successful_flytrack_MoveModel.mat';
MoveModelName11 = '20130412_S0052_successful_flytrack_MoveModel.mat';
MoveModelName12 = '20130412_S0054_successful_flytrack_MoveModel.mat';

% Manually tracked tethered flight models - Hydei 7500 fps
MoveModelName13 = 'MoveModel_IncSearch_20130412_S0002_frame0001-0036_NOstimulus_onewingbeat.mat';
MoveModelName14 = 'MoveModel_IncSearch_20130412_S0017_frame0026-0037_NOstimulus_strokerev.mat';
MoveModelName15 = 'MoveModel_IncSearch_20130412_S0027_frame0760-0774_90deg_strokerev.mat';
MoveModelName16 = 'MoveModel_NOSearch_20130412_S0002_frame0001-0036_NOstimulus_onewingbeat.mat';
MoveModelName17 = 'MoveModel_NOSearch_20130412_S0017_frame0026-0037_NOstimulus_strokerev.mat';
MoveModelName18 = 'MoveModel_NOSearch_20130412_S0027_frame0760-0774_90deg_strokerev.mat';

%Manually tracked tethered melanogaster 6000 fps
MoveModelName19 = 'MoveModel_NOsearch_20140310_S0001_melanogaster.mat'
MoveModelName20 = 'MoveModel_IncSearch_20140310_S0001.mat'
MoveModelName21 = 'MoveModel_IncSearch_20140327_S0001.mat'
MoveModelName22 = 'MoveModel_IncSearch_20140327_S0005.mat'
MoveModelName23 = 'MoveModel_IncSearch_20140401_S0001.mat'

%Successfully tracked tethered melanogaster 6000 fps
MoveModelName24 = '20140310_S0001_successful_flytrack_MoveModel_frm599_615.mat'
MoveModelName25 = '20140310_S0001_successful_flytrack_MoveModel_frm411_448.mat'
MoveModelName26 = '20140310_S0001_successful_flytrack_MoveModel_frm246_325.mat'
MoveModelName27 = '20140310_S0001_successful_flytrack_MoveModel_frm141_203.mat'
MoveModelName28 = '20140327_S0001_successful_flytrack_MoveModel_frm18_682.mat'

MoveModelNameS = {
    MoveModelName20
    MoveModelName21
    MoveModelName22
    MoveModelName23
    MoveModelName28};

% MoveModelNameS = {
%     MoveModelName19
%     MoveModelName20
%     MoveModelName21
%     MoveModelName22
%     MoveModelName23};


% MoveModelNameS = {
%     MoveModelName19};

% MoveModelNameS = {
%     MoveModelName1
%     MoveModelName2
%     MoveModelName3
%     MoveModelName4
%     MoveModelName5
%     MoveModelName6
%     MoveModelName7
%     MoveModelName8
%     MoveModelName9
%     MoveModelName10
%     MoveModelName11
%     MoveModelName12
%     MoveModelName13
%     MoveModelName14
%     MoveModelName15
%     MoveModelName16
%     MoveModelName17
%     MoveModelName18};

%% KINE file
% %manual fit file
% cd(KinePath)
% if exist('KINE')==7
%     cd('KINE')
% end

FileFromKine = [KinePath,KineName];
if exist(FileFromKine)~=2
    [KineName,KinePath] = uigetfile({'*.mat'},'Select file with initial position from Kine');
    FileFromKine = [KinePath,KineName];
end
load(FileFromKine);

%endframe = data.images.frames;

%% Movie image files from kine
FileName = data.images.file;
PathName = data.images.path;

if exist(PathName) ~= 7
    [FileName,PathName] = uigetfile('*.tif','Click on any image in video sequence');
end

%% solutionpath from kine
solutionpath = [data.save.pathname,'/'];

if exist(solutionpath)~=7
    mkdir(solutionpath)
end

%% calibration file (from kine)
CalibFile = [data.cal.pathname,data.cal.filename];

if exist(CalibFile)~=2
    [CalibFileName,CalibPathName] = uigetfile({'*.mat'},'Select calibration file with DLT coefficients','DLT_coeff.mat');
    CalibFile = [CalibPathName,CalibFileName];
end

%% Movement Model data
% MoveModelfile = [MoveModelPath,MoveModelName];
% 
% if exist([MoveModelPath,MoveModelNameS{1}])~=2
%     if exist(MoveModelPath)==7
%         cd(MoveModelPath)
%     end
%     [MoveModelNameS,MoveModelPath] = uigetfile({'*.mat'},'Select movement model file');
% %     MoveModelfile = [MoveModelPath,MoveModelName];
% end

%% const
cd(PathName)
allfiles = dir(['*',FileName(end-3:end)]);

% ICframe = find(data.kine.body.data.length~=0, 1 )
% BGframe = ICframe + 8;
ICframe = find(data.kine.body.data.length~=0, 1, 'last' )

%% Endframe as offset of IC frame
% endframe = ICframe+399;
% endframe = ICframe+950;
%endframe = ICframe + 200;
endframe = 1365;

if exist('startframe') == 0
    startframe = ICframe;
end
% 
% % startframe = str2num(allfiles(1).name(end-6:end-4));
% startframe = ICframe;

% endframe = str2num(allfiles(end).name(end-9:end-4));
% endframe = ICframe;
% endframe = 2000;

% BGframe = 1;
% BGframe = endframe;
% BGframe = ICframe-50;
% BGframe = str2num(allfiles(end).name(end-9:end-4));

numframes = endframe-startframe+1;
frames2skip = [];
framesample = 1;

% - Camera Info
numcam = 3;

% frame rate
cd(PathName)
cam_info = importdata([FileName(1:13),'.cih'],' ',12);
fps = cam_info.data(1);
dt = 1/fps;
% dt = 1/7500;  %Framerate of the camera
% dt = 1/6000;  %Framerate of the camera

% resolution
im = imread(FileName);
imgres = [size(im,1) size(im,2)];
% imgres = [512 512];

% Set to true if you need to calculate state of fly at start frame 'PAR.ICframe'.
% getIC = true;
% getBodyShape = true;

useIC = getIC; % FTMmod: dont load frame<startframe, but use ICframe


%Change this depending on how many flies you are tracking
%More than 1 Fly Is Not Working Yet!!!!!
flynum = [1 ];  %[1 2];
numfly = length(flynum);

OccludeShape = {[],[],[]};
%spline order
c = 4; 

%% This is the number of sample points used within the Fly model.
%It changes how fine the mesh is. 
L1 = 15; %# of steps for body along length
L2 = 6; %# of steps for head along length
L3 = 30; %# of steps for wing around the boundary
T1 = 13; %# of theta steps for head and body
T2 = 2; %# of steps towards center of wing

%% -- These dimensions have to be defined by the user.
%The dimensions are for a single fly.

%The subsample scale of the image is 1/2^(pwr)
etamax = 0;
paramdim = 1;
statedim = 15*ones(1,numfly); 
pNoisedim = statedim;

% This is the number of frames to scan back when performing the pattern
% matching in the motion prediction.
NumPrevSol = 5;

streams = 2;

%% tracker search distance (moved from Closest_PtsNrmlGate1, FTM 20120605)
% 
% % ebraheems
% WingSearchDist = [15 20];
% BodySearchDist = [15 15];

% FTMmod 20120627 BEST SETTINGS Dm
%WingSearchDist = [10 15];
%BodySearchDist = [10 10];

% Dh (?)
%WingSearchDist = [25 30];
%BodySearchDist = [20 20];

WingSearchDist = [120 120];
BodySearchDist = [5 5];

numfeatpts = 30;

%% Movement Model scale (VARIABLE SCALE DOES NOT WORK YET)
% MoveModelScale = [8/12 8/10 1 1.2 4/3];
% MoveModelScale = .95:.01:1.05;
MoveModelScale = 1;

BWfilterRatio = 0.1; % reduce min threshold
WingBodyRatio = 0.7; % increase wing-body treshold
% BWfilterRatio = 1; % reduce min threshold
% WingBodyRatio = 1; % increase wing-body treshold
%% VARS

% 2*standard deviation in degrees;
%angvar = 0.6;  % FTMmod
%angvar = 0.002;
angvar = 0.002;
%angvar = 2.0;

%Variance for the wing joint locations
%JointVar = [0.001 0.0004 0.0003];

JointVar = [0.001 0.0004 0.0003];
%Variance for the body linear acceleration
LinVar = [0.183 0.639 1.33];
%Variance for the body angular acceleration
AngVar = [.168 .181 2.20];


%% END ADVANCED PARAMETERS
%=======================================================

%% FTMmod: load previous ManualFit if it exists
PAR = LoadVideo(PathName,FileName,solutionpath);

if exist([PAR.solutionpath PAR.solutiondirname '/ManualFit_' PAR.solutiondirname '.mat'])==2
    
    load([PAR.solutionpath PAR.solutiondirname '/ManualFit_' PAR.solutiondirname '.mat'])
    PAR = ManualFit.ImageData;
else
    % Make solution directories
    mkdir(PAR.solutionpath,[ PAR.solutiondirname]);
    mkdir(PAR.solutionpath,['Features']);
    
    if plot_data == true
        mkdir([PAR.solutionpath,'Images/cam1']);
        mkdir([PAR.solutionpath,'Images/cam2']);
        mkdir([PAR.solutionpath,'Images/cam3']);
        mkdir([PAR.solutionpath,'Projections/cam1']);
        mkdir([PAR.solutionpath,'Projections/cam2']);
        mkdir([PAR.solutionpath,'Projections/cam3']);
    end
    
    % save new ManualFit
    ManualFit.ImageData = PAR;
    save([PAR.solutionpath '/' PAR.solutiondirname '/ManualFit_' PAR.solutiondirname],'ManualFit');
end

%% store data in PAR
PAR.plot_data = plot_data;
PAR.digits = 4;

PAR.CalibFile = CalibFile;
PAR.FileFromKine = FileFromKine;
% PAR.MoveModelfile = MoveModelfile;
PAR.MoveModelNameS = MoveModelNameS;
PAR.MoveModelPath = MoveModelPath;

PAR.startframe = startframe;
PAR.endframe = endframe;
PAR.ICframe = ICframe;
PAR.BGframe = BGframe;
PAR.frames2skip = frames2skip;
PAR.framesample = framesample;

% - Camera Info
PAR.dt = dt;  %Framerate of the camera
PAR.numcam = numcam;
PAR.imgres = imgres;

% number of flies (More than 1 Fly Is Not Working Yet!!!!!)
PAR.flynum = flynum;
PAR.numfly = numfly;

PAR.OccludeShape = OccludeShape;
%spline order
PAR.c = c; 

% number of sample points used within the Fly model. 
PAR.L1 = L1; %# of steps for body along length
PAR.L2 = L2; %# of steps for head along length
PAR.L3 = L3; %# of steps for wing around the boundary
PAR.T1 = T1; %# of theta steps for head and body
PAR.T2 = T2; %# of steps towards center of wing

%% tracker search distance (moved from Closest_PtsNrmlGate1, FTM 20120605)
% ebraheems
PAR.WingSearchDist = WingSearchDist;
PAR.BodySearchDist = BodySearchDist;
PAR.numfeatpts = numfeatpts;

%The subsample scale of the image is 1/2^(PAR.pwr)
PAR.etamax = etamax;
PAR.paramdim = paramdim;
PAR.statedim = statedim; 
PAR.pNoisedim = pNoisedim;

% number of frames to scan back when performing the pattern
% matching in the motion prediction.
PAR.NumPrevSol = NumPrevSol;

PAR.streams = streams;

% MoveModelScale
PAR.MoveModelScale = MoveModelScale;

% BW filter ratio
PAR.BWfilterRatio = BWfilterRatio;
PAR.WingBodyRatio = WingBodyRatio;
%% Get Initial Conditions
if getIC == true
    %Calculate the initial condition
    ManualFit = auto_init(ManualFit,PAR.ICframe,'BG');
    ManualFit = auto_init(ManualFit,PAR.ICframe,'IC');
    
    save([PAR.solutionpath  PAR.solutiondirname '/ManualFit_' PAR.solutiondirname],'ManualFit');

else
    %Load initial condition
    load([PAR.solutionpath  PAR.solutiondirname '/ManualFit_' PAR.solutiondirname,'.mat']);
end

% Assign model parameters & image BG
PAR.params = ManualFit.params;
PAR.DLT = ManualFit.DLT;
PAR.cam = ManualFit.cam;
% Assign image BG
PAR.BG = ManualFit.IMBG;

% ================================================
% Modify Body Shape
if getBodyShape == 1
    ManualFit = auto_init(ManualFit,PAR.ICframe,'autoradius');
    save([PAR.solutionpath  PAR.solutiondirname '/ManualFit_' PAR.solutiondirname],'ManualFit');
end

PAR.params = ManualFit.params;

%Reassign values again because they are changed in 'auotradius' mode of
%'auto_init.m'
PAR.statedim = 15*ones(1,PAR.numfly); 
PAR.pNoisedim = PAR.statedim;

%--- Initialise GSSM model from external system description script.
model = gssm_flyOcc('init');

% Define start and end frames of calculation
frames = [startframe endframe];

Arg.type = 'state';                             
Arg.tag = ['State estimation for ' PAR.stub ' data.']; 
Arg.model = model;                      

% Create inference data structure and
InfDS = geninfds(Arg);                               

% generate process and observation noise sources
[pNoise, oNoise, InfDS] = gensysnoiseds(InfDS, 'srcdkf');       

%Initialize occlusion index
InfDS.model.Occ = cell(1,length(PAR.numfly));

%--- initial estimate of state E[X(0)]
p0 = reshape(ManualFit(1).solnQ',[],1);
Xh(:,1) = p0;

% initial state covariance
Px_ = [0.001.*ones(1,3) 0.00001.*ones(1,4) angvar.*ones(1,8)];

%Create a diagonal covariance matrix by replicating Px_ # of fish
%times and then placing it on the diagonal of a matrix. 
Px = diag(repmat(Px_,1,PAR.numfly));

%--- Call inference algorithm / estimator
% Square Root Central Difference Kalman Filter
%---------------

InfDS.spkfParams  = sqrt(3);    % scale factor (CDKF parameter h)
%InfDS.spkfParams  = sqrt(10);
%THL mod
%InfDS.spkfParams  = sqrt(0.2);


Sx = chol(Px)';
srcdkf_const(Xh(:,1),Sx,pNoise,oNoise,InfDS,frames,useIC);
% srcdkf_const_variablestep(Xh(:,1),Sx,pNoise,oNoise,InfDS,frames,useIC);

stop = clock
time_used = stop-start

secs_per_frame = ((time_used(3)*24*60*60 + time_used(4)*60*60 + time_used(5)*60 + time_used(6)) / (endframe-startframe)*framesample);

% make projection movie
cd(solutionpath);
paste_top_projection_makemovie
%A_Wing_kinematics_v2_checkFlyData