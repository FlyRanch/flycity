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

close all
clear
clc

ICframe = 1716


fprintf('\nDEMSE_MODAL : This demonstration performs state estimation\n');
fprintf('              for 3D Drosophila body and wing trajectories\n'); 
fprintf('              The pixel observations are corrupted by additive\n');
fprintf('              white Gaussian noise.\n\n');

% PAR is a static structure that defines various parameters of the tracking
% sequence
global PAR

%% Add the subdirectories that contain important / necessary files
% flytracker program files


cd '/home/matt/Dropbox/WORK/flytracker/flytracker'
% addpath(pwd);
addpath('/home/matt/Dropbox/WORK/flytracker/flytracker');
addpath('/home/matt/Dropbox/WORK/flytracker/flytracker/mex/');
addpath('/home/matt/Dropbox/WORK/flytracker/flytracker/core/');
addpath('/home/matt/Dropbox/WORK/flytracker/flytracker/results/');


%% KINE file
% % melano @ 7500fps & Zdown
% KineName = 'kine_20120702_S0002_1852up_MDLT.mat';
% KinePath = '/media/FLYDATA1/solutions/20120702_S0002/';

% hydei @ 7500fps & Zup
% KineName = 'kine_20120724_S0002_fr1806-1932_movemodel_MDLTzup.mat';
% KineName = 'kine_20120724_S0002_frame1365_MDLTzup.mat';
KineName = 'kine_20130412_S0017.mat';

% KinePath = '/media/FLYDATA1/solutions/20120724_S0002/';
KinePath = '/media/Matt01/solutions/20130412_S0017';

% % manual select
KineName = [];

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

%% const
cd(PathName)
allfiles = dir(['*',FileName(end-3:end)]);

% ICframe = find(data.kine.body.data.length~=0, 1 );
% ICframe = find(data.kine.body.data.length~=0, 1, 'last' );

% startframe = str2num(allfiles(1).name(end-6:end-4));
startframe = ICframe;

endframe = str2num(allfiles(end).name(end-9:end-4));
% endframe = ICframe;

% BGframe = 1;
BGframe = endframe;
% BGframe = ICframe-50;

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
getIC = 1;
getBodyShape = 1;

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
WingSearchDist = [10 20];
BodySearchDist = [15 15];

numfeatpts = 30;

BWfilterRatio = 0.5; % reduce min threshold
WingBodyRatio = 0.5; % increase wing-body treshold
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

getImageData = true;

if getImageData
    PAR = LoadVideo(PathName,FileName,solutionpath);
        
    ManualFit.ImageData = PAR;
    
    % Make two directories to save the estimated state and the features into
    % if they don't already exist
    if exist([PAR.solutionpath,PAR.solutiondirname]) ~= 7
        mkdir(PAR.solutionpath,[ PAR.solutiondirname]);
        mkdir(PAR.solutionpath,['Features']);
        mkdir([PAR.solutionpath,'Images/cam1']);
        mkdir([PAR.solutionpath,'Images/cam2']);
        mkdir([PAR.solutionpath,'Images/cam3']);
        mkdir([PAR.solutionpath,'Projections/cam1']);
        mkdir([PAR.solutionpath,'Projections/cam2']);
        mkdir([PAR.solutionpath,'Projections/cam3']);
    end

    save([PAR.solutionpath '/' PAR.solutiondirname '/' 'ManualFit_' PAR.solutiondirname],'ManualFit');
else
    %Load the previously stored data
    [FileName,PathName] = uigetfile({'*.mat'},'Select "ManualFit" data file for the video sequence',solutionpath);
    load([PathName FileName]);
    
    % Check that the paths stored in 'ManualFit' match the location that you
    % just selected the file from
    % If 75% of the paths match , then we'll assume everything's okay.
    endd = round(.75*length(ManualFit.ImageData.solutionpath));
    
    if strcmp(PathName(1:endd),ManualFit.ImageData.solutionpath(1:endd)) == 0
        %if different, run Loadvideo routine
        warning('The directories stored in ManualFit structure do not match its current location.\n You are prompted to relocate the directories',[]);
        
        PAR = LoadVideo;
        ManualFit.ImageData = PAR;
        save([PAR.solutionpath  PAR.solutiondirname '/' 'ManualFit_' PAR.solutiondirname],'ManualFit');
    else
        PAR = ManualFit.ImageData;
    end
end

%% store data in PAR

PAR.CalibFile = CalibFile;
PAR.FileFromKine = FileFromKine;

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

% BW filter ratio
PAR.BWfilterRatio = BWfilterRatio;
PAR.WingBodyRatio = WingBodyRatio;

%% Get Initial Conditions
if getIC == true
%     %Calculate the initial condition
%     ManualFit = auto_init(ManualFit,PAR.ICframe,'BG');
%     ManualFit = auto_init(ManualFit,PAR.ICframe,'IC');
%     
%     save([PAR.solutionpath  PAR.solutiondirname '/ManualFit_' PAR.solutiondirname],'ManualFit');
%     
    %% make movement model
%     MoveModel = auto_init_multi(ManualFit,PAR.ICframe,'IC_mult');
%     save([PAR.solutionpath,'MoveModel_IncSearch' KineName(5:end)],'MoveModel');

    MoveModel = auto_init_multi_nosearch(ManualFit,PAR.ICframe,'IC_mult');
    save([PAR.solutionpath,'MoveModel_NOsearch' KineName(5:end)],'MoveModel');
else
    %Load initial condition
    load([PAR.solutionpath  PAR.solutiondirname '/ManualFit_' PAR.solutiondirname]);
end

% % Assign model parameters & image BG
% PAR.params = ManualFit.params;
% PAR.DLT = ManualFit.DLT;
% PAR.cam = ManualFit.cam;
% % Assign image BG
% PAR.BG = ManualFit.IMBG;
% 
% % ================================================
% % Modify Body Shape
% if getBodyShape == 1
%     ManualFit = auto_init(ManualFit,PAR.ICframe,'autoradius');
%     save([PAR.solutionpath  PAR.solutiondirname '/ManualFit_' PAR.solutiondirname],'ManualFit');
% end
% 
% PAR.params = ManualFit.params;
% 
% %Reassign values again because they are changed in 'auotradius' mode of
% %'auto_init.m'
% PAR.statedim = 15*ones(1,PAR.numfly); 
% PAR.pNoisedim = PAR.statedim;
% 
% %--- Initialise GSSM model from external system description script.
% model = gssm_flyOcc('init');
% 
% % Define start and end frames of calculation
% frames = [startframe endframe];
% 
% Arg.type = 'state';                             
% Arg.tag = ['State estimation for ' PAR.stub ' data.']; 
% Arg.model = model;                      
% 
% % Create inference data structure and
% InfDS = geninfds(Arg);                               
% 
% % generate process and observation noise sources
% [pNoise, oNoise, InfDS] = gensysnoiseds(InfDS, 'srcdkf');       
% 
% %Initialize occlusion index
% InfDS.model.Occ = cell(1,length(PAR.numfly));
% 
% %--- initial estimate of state E[X(0)]
% p0 = reshape(ManualFit(1).solnQ',[],1);
% Xh(:,1) = p0;
% 
% % initial state covariance
% Px_ = [0.001.*ones(1,3) 0.00001.*ones(1,4) angvar.*ones(1,8)];
% 
% %Create a diagonal covariance matrix by replicating Px_ # of fish
% %times and then placing it on the diagonal of a matrix. 
% Px = diag(repmat(Px_,1,PAR.numfly));
% 
% %--- Call inference algorithm / estimator
% % Square Root Central Difference Kalman Filter
% %---------------
% InfDS.spkfParams  = sqrt(3);    % scale factor (CDKF parameter h)
% Sx = chol(Px)';
% srcdkf_const(Xh(:,1),Sx,pNoise,oNoise,InfDS,frames);