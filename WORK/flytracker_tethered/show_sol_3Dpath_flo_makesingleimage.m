%show_sol_3Dpath.m
%
%show_sol_auto makes a movie of the tracking result for a particular video
%sequence.  It works similar to 'paste_image_auto' but iterates over each
%frame
%
%After running this program the user should call 'movie2avi(M2,...)' to
%save the movie to whatever filename they choose.
% addpath(pwd);
% addpath('/home/florian/WORK/flytracker');
% addpath('/home/florian/WORK/flytracker/mex/');
% addpath('/home/florian/WORK/flytracker/core/');
% addpath('/home/florian/WORK/flytracker/results/');

addpath('/home/florian/Dropbox/WORK/flytracker/flytracker');
addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/mex/');
addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/core/');
addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/results/');

clear;close all;clc

datadir = cd;
cd ..
% mkdir('3Drender')
% cd('3Drender')
soldir = cd;

% movdir = cd;
% 
% movfilename = '3Drender_startattriggerpoint_2500fps.avi';
% mov = avifile(movfilename, 'compression', 'none');

cd(datadir)

% PAR.videopath = '../FlyVideo/';
% PAR.filetag = 'exp101000';
% 
% PAR.solutionpath = [PAR.videopath '/solutions/'];
% PAR.stub = [PAR.filetag];
fignum = 1;

%changes the scaling of the final figure, as percentage of original image
%resolution to make fit on screen.
plotscale = .8;
plotscale = 1;

% --------------------------------
% Load the ManualFit Parameters
% load([PAR.solutionpath 'fly_' PAR.stub '/' 'ManualFit_' PAR.stub]);
% load('ManualFit_frame');
load('ManualFit_flytracks');
% load('ManualFit_exp095000');

%----------------------------------------------------
%load video data into buffer
%----------------------------------------------------
allfiles=dir('*.mat')
for n=2:length(allfiles)
    allnums(n-1) = str2num(allfiles(n).name(4:end-4));
end
startframe = min(allnums);
% startframe = 2577;
endframe = max(allnums);
samplerate = 126;

% startframe = 3866
% endframe = 5000
% samplerate = 1;
% samplerate = 3;

frames = [startframe endframe]; % number of frames in movie

movidx = frames(1):samplerate:frames(2);
numframes = length(movidx);


% %% standard colormaps
% cmap = colormap(jet)
% 
% cmap = colormap(autumn)
% cmap = flipdim(cmap,1)

%% custom colormaps
% white to dark
cmap1 = [1:-(1-.75)/100:0.75];
cmap2 = [1:-(1-0)/100:0];
cmap3 = [1:-(1-0)/100:0];

% % grey to dark
% cmap1 = [.5:-(.5-.75)/100:0.75];
% cmap2 = [.5:-(.5-0)/100:0];
% cmap3 = [.5:-(.5-0)/100:0];

% % light
% cmap1 = [1:-(1-.99)/100:0.99];
% cmap2 = [1:-(1-0)/100:0];
% cmap3 = [1:-(1-0)/100:0];

% grey to red
cmap = [cmap1; cmap2; cmap3]';

% % grey to green
% cmap = [cmap2; cmap1; cmap3]';
% 
% % grey to blue
% cmap = [cmap3; cmap2; cmap1]';

for i = 1:numframes
    n = round(((length(cmap)-1) / (numframes-1))*(i-1) +1);
    colors(i,:) = cmap(n,:);
end

% colormap(cmap)
% colormap(colors)

cptsTag = 0;     %Set to 1 to see the corresponding points
nrmlTag = 0;     %Set to 1 to see the normal vectors
KF_tag = 0;      %Set to 1 to see the Predicted Model and Gating Ellipses
interptag = 1;
%--------------------------------
% Define the Tracking Parameters
PAR.pixpermm = 1;
PAR.numfly = 1;
%Number of parameters of the model (i.e. 8 control points)
PAR.mdlpar = 15*ones(1,PAR.numfly);
PAR.statedim = PAR.mdlpar;
PAR.modelfun_H = @modcurvesplineP;
PAR.etamax = 0;

%spline order
PAR.c = 4;
PAR.L1 = 15; %# of steps for body along length
PAR.L2 = 6; %# of steps for head along length
PAR.L3 = 25; %# of steps for wing around the boundary
PAR.T1 = 13; %# of theta steps for head and body
PAR.T2 = 2; %# of steps towards center of wing

% - Camera Info
% PAR.dt = 1/6000;  %Framerate of the camera
PAR.dt = 1/7500;  %Framerate of the camera
PAR.numcam = 3;
solidx = [1 length(movidx)];

% Assign model parameters
PAR.params = ManualFit.params;
PAR.DLT = ManualFit.DLT;
PAR.cam = ManualFit.cam;


soln1 = zeros(length(movidx),PAR.numfly*PAR.statedim);
SOLN = zeros(length(frames(1):frames(2)),PAR.numfly*PAR.statedim);
if interptag == 1
    for i=1:length(movidx)
%         load(['fly' num2str(movidx(i)) '.mat']);
        load(['fly'  sprintf(['%0',num2str(6),'d'], movidx(i)) '.mat']);

        soln1(i,:) = xh(1:PAR.numfly*PAR.statedim)';
        clear xh InternalVariablesDS
    end
end

% %Now, perform interpolation to calculate the state at in between frames
% for k = 1:size(soln1,2)
%     SOLN(:,k) = interp1(movidx,soln1(:,k),frames(1):frames(2),'spline')';
% end
% 
% no inperp
SOLN = soln1;

%% Calculate the region of fly motion
buffer = 3;
minx = min(SOLN(:,1)) - buffer;
maxx = max(SOLN(:,1)) + buffer;
miny = min(SOLN(:,2)) - buffer;
maxy = max(SOLN(:,2)) + buffer;
minz = min(SOLN(:,3)) - buffer;
maxz = max(SOLN(:,3)) + buffer;

%--------------------------------------------
%calculate theta-dot
theta_vec = sqrt(sum(SOLN(:,4:6).^2,2));
theta_dot = gradient(theta_vec,PAR.dt);
dt = PAR.dt;
tsteps = 1;
dummy = [];
torder = 1;

M1dum = [];
M2dum = [];

T = 0;%linspace(0,12,60);

% M2 = moviein(size(SOLN));
for i= 1:size(SOLN,1)
%     kk = frames(1)+(i-1);

    clear flymodQ
    [x,y,z] = flymodQ(SOLN(i,:),PAR.params,PAR);
    for j = 1:length(x);
        PAR.modsample(j) = size(x{j},1);
    end

       
    
%     close all
    figure(fignum);
 set(fignum,'units','pixels','position',[0 0 1024 1024].*plotscale);
%  figure('units','pixels','position',[0 0 1024 1024].*plotscale);
view(3)
%  clf;
az = -127.5;
el = 45;
view(az, el);


    for k = 1:length(x);
%         surf(x{k},y{k},z{k},'facecolor','r','edgecolor','k','facelighting','phong');
        surf(x{k},y{k},z{k},'facecolor',colors(i,:),'edgecolor','k','facelighting','phong');
        hold on;
    end

    set(gca,'xlim',[minx maxx],'ylim',[miny maxy],'zlim',[minz maxz]);
%     set(gca,'xlim',[minx maxx],'ylim',[miny maxy],'zlim',[minz maxz],...
%         'zdir','reverse','ydir','reverse');
%     axis vis3d


    set(gca,'XTickLabel',[])
    set(gca,'YTickLabel',[])
    set(gca,'ZTickLabel',[])
    
    figure(fignum);
    set(fignum,'color','w');
    %title([filetag ', Frame ' num2str(kk) ]);%', x = ' num2str(soln1(i,:))])
    

%     M2(i) = getframe(fignum);
%     axis tight
%     axis equal
%     M2(kk) = getframe(fignum);
% cd(soldir)
% tic
%     saveas(gca,[soldir '/3DrenderIm' sprintf(['%0',num2str(6),'d'], movidx(i)) '.tif'])
% toc
% cd(datadir)

% cd(movdir)
%     mov = addframe(mov,gcf);
% cd(datadir)
end
camlight('headlight')
% colorbar
% colorbar('YTickLabel',[]);
% 
% cd(movdir)
% mov = close(mov);
% cd(datadir)
cd(soldir)
tic
    saveas(gca,[soldir '/3DrenderImage.fig'])
%     saveas(gca,[soldir '/3DrenderImage.ai'])
%     saveas(gca,[soldir '/3DrenderImage.pdf'])
%     saveas(gca,[soldir '/3DrenderImage.eps'])
%     saveas(gca,[soldir '/3DrenderImage.png'])
%     saveas(gca,[soldir '/3DrenderImage.tif'])
toc
cd(datadir)
