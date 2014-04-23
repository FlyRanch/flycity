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
start = clock;
load_deployment_vars
getIC = true;
% software paths
addpath(fullfile(software_root_dir,'flytracker_library'));
addpath(fullfile(software_root_dir,'mex'));
addpath(fullfile(software_root_dir,'core'));
addpath(fullfile(software_root_dir,'results'));
addpath(fullfile(software_root_dir,'SpinConv'));
addpath(fullfile(software_root_dir,'postproc','AA_check_wbkin'));
addpath(fullfile(software_root_dir,'postproc','AA_check_wbkin','WBKinProg_Johan'));
addpath(fullfile(software_root_dir,'postproc','AA_check_wbkin','WBKinProg_Johan','Kalman_Filter'));
addpath(fullfile(software_root_dir,'postproc','AA_check_wbkin','WBKinProg_Johan','Plot'));
addpath(fullfile(software_root_dir,'postproc','AA_check_wbkin','WBKinProg_Johan','Plot','Smoothing_plots'));
addpath(fullfile(software_root_dir,'postproc','AA_check_wbkin','WBKinProg_Johan','Body_model'));
%addpath(fullfile(software_root_dir,'postproc','AA_check_wbkin','WBKinProg_Johan','Body_model','mex'));
addpath(fullfile(software_root_dir,'postproc','AA_check_wbkin','WBKinProg_Johan','Type2Regression'));

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
if exist(cam1_seq_path) ~= 7
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
cd(cam1_seq_path)
ICframe = find(data.kine.body.data.length~=0, 1, 'first' )

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
cd(cam1_seq_path)
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
PAR = LoadVideo(cam1_seq_path,cam1_file1,seq_sol_path);

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
    %% make movement model
    PAR.digits = 4;
    MoveModel = auto_init_multi(ManualFit,PAR.ICframe,'IC_mult');
    save([PAR.solutionpath,'MoveModel_IncSearch_' sprintf(sol_format_str,date_str,seq_number)],'MoveModel');

    % MoveModel = auto_init_multi_nosearch(ManualFit,PAR.ICframe,'IC_mult');
    % save([PAR.solutionpath,'MoveModel_NOsearch' KineName(5:end)],'MoveModel');
else
    %Load initial condition
    load([PAR.solutionpath  PAR.solutiondirname '/ManualFit_' PAR.solutiondirname]);
end
