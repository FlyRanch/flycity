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

% software paths
addpath(fullfile(software_root_dir,'mex');
addpath(fullfile(software_root_dir,'core');
addpath(fullfile(software_root_dir,'results');
addpath(fullfile(software_root_dir,'postproc','AA_check_wbkin');

%% Add the subdirectories that contain important / necessary files
% flytracker program files
% ----get this from the deployment vars
warning off

%% Load kine and calibration data
% --- calculate from deplyment var strings
%Check to make sure the kine data file exists if so load
if exist(FileFromKine)~=2
    [KineName,KinePath] = uigetfile({'*.mat'},'Select file with initial position from Kine');
    FileFromKine = [KinePath,KineName];
end
load(FileFromKine);
%check to make sure that we have a good first file to work with
if exist(PathName) ~= 7
    [cam1_file1,cam1_seq_path] = uigetfile('*.tif','Click on any image in video sequence');
end
%check to make sure the solution path for the sequence exists
if exist(seq_sol_path)~=7
    mkdir(seq_sol_path)
end
% get the calibration file from the kine file - probably best to keep it 
%this way, but there is a potential platform dependency created here
CalibFile = [data.cal.pathname,data.cal.filename];
if exist(CalibFile)~=2
    [CalibFileName,CalibPathName] = uigetfile({'*.mat'},'Select calibration file with DLT coefficients','DLT_coeff.mat');
    CalibFile = [CalibPathName,CalibFileName];
end

%% const
cd(PathName)
ICframe = find(data.kine.body.data.length~=0, 1, 'last' )

%% Endframe as offset of IC frame

if exist('startframe') == 0
    startframe = ICframe;
end

numframes = endframe-startframe+1;
frames2skip = [];
framesample = 1;

% - Camera Info
numcam = 3;

% frame rate
cd(PathName)
cam_info = importdata(cih_file,' ',12);
fps = cam_info.data(1);
dt = 1/fps;

% resolution
im = imread(cam1_file1);
imgres = [size(im,1) size(im,2)];

% Set to true if you need to calculate state of fly at start frame 'PAR.ICframe'.
% getIC = true;
% getBodyShape = true;

useIC = getIC; % FTMmod: dont load frame<startframe, but use ICframe

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