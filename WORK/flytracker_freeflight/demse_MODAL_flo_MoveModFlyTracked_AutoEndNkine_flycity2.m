% DEMSE_MODAL  Demonstrate state estimation of 3D drosophila
% kinematic model 
%
%
%   See also in this script:
%   auto_init - initializes the tracker by calculating calibration,
%   background, and model shape
%   
%   gssm_flyOcc - defines the motion and measurement model of the fish
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
function demse_MODAL_flo_MoveModFlyTracked_AutoEndNkine_flycity2
close all
% clear
clc

start = clock;
%ICframe = 2518
%startframe = ICframe
% startframe =1839

endframe =2700
% cd flytracks
% frame_dir = dir;
% endframe = frame_dir(end).name;
% endframe = str2num(endframe(4:end-4));
% cd ..

BGframe =1

% default settings
% start sequence for first time (fit model to KINEfit and do bodyscaling)
getIC = true;
getBodyShape = true;

% % continue after tracker has lost fly (fit model to new KINEfit but no bodyscaling)
% getIC = true;
% getBodyShape = false;
% % 
% % continue after tracker has crashed (reload previous frame and model)
% getIC = false;
% getBodyShape = false;

%% to plot or not to plot
% plot_data = true
plot_data = false


fprintf('\nDEMSE_MODAL : This demonstration performs state estimation\n');
fprintf('              for 3D Drosophila body and wing trajectories\n'); 
fprintf('              The pixel observations are corrupted by additive\n');
fprintf('              white Gaussian noise.\n\n');

% PAR is a static structure that defines various parameters of the tracking
% sequence
global PAR

%% Add the subdirectories that contain important / necessary files
% flytracker program files


% cd '/home/florian/Dropbox/WORK/flytracker/flytracker'
% addpath(pwd);
addpath('/home/matt/Dropbox/WORK/flytracker');
addpath('/home/matt/Dropbox/WORK/flytracker/mex/');
addpath('/home/matt/Dropbox/WORK/flytracker/core/');
addpath('/home/matt/Dropbox/WORK/flytracker/results/');
addpath('/home/matt/Dropbox/WORK/postproc/AA_check_wbkin/')

addpath('/home/florian/Dropbox/FlyCity/WORK/flytracker');
addpath('/home/florian/Dropbox/FlyCity/WORK/flytracker/mex/');
addpath('/home/florian/Dropbox/FlyCity/WORK/flytracker/core/');
addpath('/home/florian/Dropbox/FlyCity/WORK/flytracker/results/');
addpath('/home/florian/Dropbox/FlyCity/WORK/postproc/AA_check_wbkin/')

addpath('/home/flycity/Dropbox/WORK/flytracker');
addpath('/home/flycity/Dropbox/WORK/flytracker/mex/');
addpath('/home/flycity/Dropbox/WORK/flytracker/core/');
addpath('/home/flycity/Dropbox/WORK/flytracker/results/');
addpath('/home/flycity/Dropbox/WORK/postproc/AA_check_wbkin/')

%% KINE file
% manual select
KinePath = [];
KineName = [];

% AUTOLOAD
KinePath = [cd,'/'];
KineName_dir = dir('kine*');
KineName = KineName_dir(1).name;

%% Movement Model data

MoveModelPath = '/home/flycity/Dropbox/Models/';

MoveModelName1 = '20121002_S0003_successful_flytrack_MoveModel.mat'; % florians
MoveModelName2 = '20121009_S0003_successful_flytrack_MoveModel.mat';
MoveModelName3 = '20121010_S0001_successful_flytrack_MoveModel.mat';
MoveModelNameS = {MoveModelName1 MoveModelName2 MoveModelName3}';

% % Free flight models - Hydei 7500 fps
% MoveModelName1 = 'MoveModel_IncSearch_20120724_S0002_frame1654-1672_Dh_strokerev.mat';
% MoveModelName2 = 'MoveModel_IncSearch_20120724_S0002_frame1806-1932_movemodel_MDLTzup.mat';
% MoveModelName3 = 'MoveModel_NOsearch_20120724_S0002_frame1654-1672_Dh_strokerev.mat';
% MoveModelName4 = 'MoveModel_NOsearch_20120724_S0002_frame1806-1932_movemodel_MDLTzup.mat';
% 
% % Tethered flight models - Hydei 7500 fps
% MoveModelName5 = 'MoveModel_IncSearch_20130412_S0002_frame0001-0036_NOstimulus_onewingbeat.mat'; % one wing beat - no stimulus
% MoveModelName6 = 'MoveModel_NOSearch_20130412_S0002_frame0001-0036_NOstimulus_onewingbeat.mat';
% MoveModelName7 = 'MoveModel_IncSearch_20130412_S0017_frame0026-0037_NOstimulus_strokerev.mat'; % stroke reversal - no stimulus
% MoveModelName8 = 'MoveModel_NOSearch_20130412_S0017_frame0026-0037_NOstimulus_strokerev.mat';
% MoveModelName9 = 'MoveModel_IncSearch_20130412_S0027_frame0760-0774_90deg_strokerev.mat'; % stroke reversal - stimulus from right
% MoveModelName10 = 'MoveModel_NOSearch_20130412_S0027_frame0760-0774_90deg_strokerev.mat';
% 
% % Successful flytrack models - Hydei 7500 fps
% MoveModelName11 = '20130401_S0003_MoveModel.mat';
% MoveModelName12 = '20130401_S0005_MoveModel.mat';
% MoveModelName13 = '20130401_S0006_MoveModel.mat';
% MoveModelName14 = '20130401_S0007_MoveModel.mat';
% MoveModelName15 = '20130401_S0009_MoveModel.mat';
% MoveModelName16 = '20130401_S0011_MoveModel.mat';
% 
% MoveModelName17 = '20121002_S0003_MoveModel.mat'; % florians
% MoveModelName18 = '20121009_S0003_MoveModel.mat';
% MoveModelName19 = '20121010_S0001_MoveModel.mat';
% 
% % Use a single model
% % MoveModelNameS = {MoveModelName1};
% 
% % Use multiple models
% MoveModelNameS = {MoveModelName1 MoveModelName2 MoveModelName3 MoveModelName4 MoveModelName5 MoveModelName6 MoveModelName7 MoveModelName8 MoveModelName9 MoveModelName10 MoveModelName11 MoveModelName12 MoveModelName13 MoveModelName14 MoveModelName15 MoveModelName16 MoveModelName17 MoveModelName18 MoveModelName19}';


%% KINE file
% % manual fit file
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

if exist([MoveModelPath,MoveModelNameS{1}])~=2
    if exist(MoveModelPath)==7
        cd(MoveModelPath)
    end
    [MoveModelNameS,MoveModelPath] = uigetfile({'*.mat'},'Select movement model file');
%     MoveModelfile = [MoveModelPath,MoveModelName];
end

%% const
cd(PathName)
allfiles = dir(['*',FileName(end-3:end)]);

ICframe = find(data.kine.body.data.length~=0, 1 )
% ICframe = find(data.kine.body.data.length~=0, 1, 'last' )

% startframe = str2num(allfiles(1).name(end-6:end-4));

if exist('startframe') == 0
    startframe = ICframe;
end

% endframe = str2num(allfiles(end).name(end-9:end-4));
% endframe = ICframe;
% endframe = 4384;

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
 WingSearchDist = [10 15];
 BodySearchDist = [10 10];

% Dh (?)
WingSearchDist = [15 20];
BodySearchDist = [15 15];

numfeatpts = 30;

%% Movement Model scale (VARIABLE SCALE DOES NOT WORK YET)
% MoveModelScale = [8/12 8/10 1 1.2 4/3];
% MoveModelScale = .95:.01:1.05;
MoveModelScale = 1;

BWfilterRatio = 0.5; % reduce min threshold
WingBodyRatio = 0.5; % increase wing-body treshold
% BWfilterRatio = 1; % reduce min threshold
% WingBodyRatio = 1; % increase wing-body treshold
%% VARS

% 2*standard deviation in degrees;
angvar = 0.02;  % FTMmod
% angvar = 0.002;

%Variance for the wing joint locations
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


% getImageData = true;
% 
% if getImageData
%     PAR = LoadVideo(PathName,FileName,solutionpath);
%         
%     ManualFit.ImageData = PAR;
%     
%     % Make two directories to save the estimated state and the features into
%     % if they don't already exist
%     if exist([PAR.solutionpath,PAR.solutiondirname]) ~= 7
%         mkdir(PAR.solutionpath,[ PAR.solutiondirname]);
%         mkdir(PAR.solutionpath,['Features']);
%         mkdir([PAR.solutionpath,'Images/cam1']);
%         mkdir([PAR.solutionpath,'Images/cam2']);
%         mkdir([PAR.solutionpath,'Images/cam3']);
%         mkdir([PAR.solutionpath,'Projections/cam1']);
%         mkdir([PAR.solutionpath,'Projections/cam2']);
%         mkdir([PAR.solutionpath,'Projections/cam3']);
%     end
% 
%     save([PAR.solutionpath PAR.solutiondirname '/' 'ManualFit_' PAR.solutiondirname],'ManualFit');
% else
%     %Load the previously stored data
%     [FileName,PathName] = uigetfile({'*.mat'},'Select "ManualFit" data file for the video sequence',solutionpath);
%     load([PathName FileName]);
%     
%     % Check that the paths stored in 'ManualFit' match the location that you
%     % just selected the file from
%     % If 75% of the paths match , then we'll assume everything's okay.
%     endd = round(.75*length(ManualFit.ImageData.solutionpath));
%     
%     if strcmp(PathName(1:endd),ManualFit.ImageData.solutionpath(1:endd)) == 0
%         %if different, run Loadvideo routine
%         warning('The directories stored in ManualFit structure do not match its current location.\n You are prompted to relocate the directories',[]);
%         
%         PAR = LoadVideo;
%         ManualFit.ImageData = PAR;
%         save([PAR.solutionpath  PAR.solutiondirname '/' 'ManualFit_' PAR.solutiondirname],'ManualFit');
%     else
%         PAR = ManualFit.ImageData;
%     end
% end


%% store data in PAR
PAR.plot_data = plot_data;

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
    load([PAR.solutionpath  PAR.solutiondirname '/ManualFit_' PAR.solutiondirname]);
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
Sx = chol(Px)';
srcdkf_const(Xh(:,1),Sx,pNoise,oNoise,InfDS,frames,useIC);
% srcdkf_const_variablestep(Xh(:,1),Sx,pNoise,oNoise,InfDS,frames,useIC);

stop = clock
time_used = stop-start

secs_per_frame = (time_used(3)*24*60*60 + time_used(4)*60*60 + time_used(5)*60 + time_used(6)) / (endframe-startframe)*framesample

% make projection movie
% cd(solutionpath);
% paste_top_projection_makemovie_skip3

cd(solutionpath);
paste_top_projection_makemovie
A_Wing_kinematics_v2_checkFlyData