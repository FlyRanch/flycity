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

fprintf('\nDEMSE_MODAL : This demonstration performs state estimation\n');
fprintf('              for 3D Drosophila body and wing trajectories\n'); 
fprintf('              The pixel observations are corrupted by additive\n');
fprintf('              white Gaussian noise.\n\n');

% PAR is a static structure that defines various parameters of the tracking
% sequence
global PAR

%% Add the subdirectories that contain important / necessary files
% flytracker program files


% cd '/home/florian/WORK/flytracker'
% addpath(pwd);
addpath('/home/florian/WORK/flytracker');
addpath('/home/florian/WORK/flytracker/mex/');
addpath('/home/florian/WORK/flytracker/core/');
addpath('/home/florian/WORK/flytracker/results/');


% select movie image files

[FileName,PathName] = uigetfile({'*.bmp';'*.tif'},'Click on any image in video sequence');

% tether test
% PathName = '/home/florian/DATA/flighttracker/flighttracks/Takeoff/exp095/cam001/';
% FileName = 'exp095000001.bmp';
% PathName = '/home/florian/DATA/flighttracker/tether_test/flighttracks/120328_6000fps_54000shutter_5.6aperture/cam001/';
% FileName = 'frame001.bmp';
% 
cd(PathName)
allfiles = dir(['*',FileName(end-3:end)]);


% solutions path
cd ..
startpath = [cd,'/solutions/'];

% calibration file

cd(PathName)
% cd ..
% cd ..
% cd ..
% cd calibration
[CalibFileName,CalibPathName] = uigetfile({'*.mat'},'Select calibration file with DLT coefficients','cali_points_rotated.mat');
CalibFile = [CalibPathName,CalibFileName];

% tether test
% CalibFile = '/home/florian/DATA/flighttracker/tether_test/calibration/120328_CALI_BACK2FRONT/cali_points_rotated.mat';
% CalibFile = '/home/florian/DATA/flighttracker/tether_test/calibration/120328_CALI_BACK2FRONT/cali_points_rotated_zup.mat';
% CalibFile = '/home/florian/DATA/flighttracker/tether_test/calibration/120328_CALI_BACK2FRONT/cali_points_rotated_switchyz_zup.mat';

% manual fit file

cd(PathName)
cd ..
[KineName,KinePath] = uigetfile({'*.mat'},'Select file with initial position from Kine');
FileFromKine = [KinePath,KineName];

% FileFromKine = '/home/florian/WORK/flytracker/ICFromKine/KineManualFit_exp095f430';
% FileFromKine = '/home/florian/DATA/flighttracker/tether_test/flighttracks/120328_6000fps_54000shutter_5.6aperture/kin_init.mat';
% FileFromKine = '/home/florian/DATA/flighttracker/tether_test/flighttracks/120328_6000fps_54000shutter_5.6aperture/kin_init_zup.mat';
% FileFromKine = '/home/florian/DATA/flighttracker/tether_test/flighttracks/120328_6000fps_54000shutter_5.6aperture/kin_init_switchyz_zup.mat';


%% const

% Kine frame
load(FileFromKine);

% ICframe = 20;
ICframe = min(find(data.kine.body.data.length~=0));
ICframe = max(find(data.kine.body.data.length~=0));

% startframe = str2num(allfiles(1).name(end-6:end-4));
startframe = ICframe;

endframe = str2num(allfiles(end).name(end-6:end-4));
% ICframe = 23;
% startframe = 23;
% endframe = 100;



numframes = endframe-startframe+1;
frames2skip = [];
framesample = 1;

% - Camera Info
dt = 1/6000;  %Framerate of the camera
numcam = 3;

cd(PathName)
im = imread(allfiles(1).name);
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

%% END ADVANCED PARAMETERS
%=======================================================

getImageData = true;

if getImageData
    PAR = LoadVideo(PathName,FileName,startpath);
        
    ManualFit.ImageData = PAR;
    
    % Make two directories to save the estimated state and the features into
    % if they don't already exist
    if exist([PAR.solutionpath ['fly_' PAR.stub]],'dir') ~= 7
        mkdir(PAR.solutionpath,['fly_' PAR.stub]);
        mkdir(PAR.solutionpath,['Features_' PAR.stub]);
    end

    save([PAR.solutionpath 'fly_' PAR.stub '/' 'ManualFit_' PAR.stub],'ManualFit');
else
    %Load the previously stored data
    [FileName,PathName] = uigetfile({'*.mat'},'Select "ManualFit" data file for the video sequence',startpath);
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
        save([PAR.solutionpath 'fly_' PAR.stub '/' 'ManualFit_' PAR.stub],'ManualFit');
    else
        PAR = ManualFit.ImageData;
    end
end

%% store data in PAR

PAR.CalibFile = CalibFile;
PAR.FileFromKine = FileFromKine;

PAR.ICframe = ICframe;
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

%The subsample scale of the image is 1/2^(PAR.pwr)
PAR.etamax = etamax;
PAR.paramdim = paramdim;
PAR.statedim = statedim; 
PAR.pNoisedim = pNoisedim;

% number of frames to scan back when performing the pattern
% matching in the motion prediction.
PAR.NumPrevSol = NumPrevSol;

PAR.streams = streams;

%% Get Initial Conditions
if getIC == true
    %Calculate the initial condition
    ManualFit = auto_init(ManualFit,PAR.ICframe,'BG');
    ManualFit = auto_init(ManualFit,PAR.ICframe,'IC');
    
    save([PAR.solutionpath 'fly_' PAR.stub '/ManualFit_' PAR.stub],'ManualFit');
    
else
    %Load initial condition
    load([PAR.solutionpath 'fly_' PAR.stub '/' 'ManualFit_' PAR.stub]);
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
    save([PAR.solutionpath 'fly_' PAR.stub '/ManualFit_' PAR.stub],'ManualFit');
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

% 2*standard deviation in degrees;
angvar = 0.002;
%Variance for the wing joint locations
JointVar = [0.001 0.0004 0.0003];

%Variance for the body linear acceleration
LinVar = [0.183 0.639 1.33];
%Variance for the body angular acceleration
AngVar = [.168 .181 2.20];
  
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
srcdkf_const(Xh(:,1),Sx,pNoise,oNoise,InfDS,frames);