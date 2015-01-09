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
% 
% addpath('/home/florian/Dropbox/WORK/flytracker/flytracker');
% addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/mex/');
% addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/core/');
% addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/results/');

addpath('/home/florian/GitHub/flycity/flytracker_freeflight');
addpath('/home/florian/GitHub/flycity/flytracker_freeflight/mex/');
addpath('/home/florian/GitHub/flycity/flytracker_freeflight/core/');
addpath('/home/florian/GitHub/flycity/flytracker_freeflight/results/');

clear;close all;clc

db = dir('kinflight*mat')
db = db.name;
load(db)

%% SETTINGS
samplerate = 1;

plotscale = .8;
plotscale = 1;

parallaxAngle = 1 ;       % You can change this to suit yourself.
fignum = 1;

body_zoom = 1
% body_zoom = 2;

buffer = 3;

cutx = 600
cuty = 500
rez=300; %resolution (dpi) of final graphic

%% model
load('ManualFit_4bodymodel.mat');

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

% Assign model parameters
PAR.params = ManualFit.params;
PAR.DLT = ManualFit.DLT;
PAR.cam = ManualFit.cam;

%% SEQ NR LOOP
for seq_nr = 1:size(pathDB.pos,2)

    SEQ = size(pathDB.pos,2) - seq_nr
    % 
% seqs = [1 3 5 9 11 12];
% % rot1=[-155 -15 -190 -80 110 -10];
% rot1=[0 0 0 0 0 0];
% rot2=[30 30 30 30 30 30];
% 
% for n=1:length(seqs)
% for n=4
% 
%     seq_nr = seqs(n)
% 
%     SEQ = size(seqs,2) - n
%     

    
mkdir(['3Drender_phong_7500fps_seq',num2str(seq_nr),'_zoom',num2str(body_zoom),'x'])
cd(['3Drender_phong_7500fps_seq',num2str(seq_nr),'_zoom',num2str(body_zoom),'x'])

pos(:,:) = 1/body_zoom*1000*pathDB.pos(:,seq_nr,:);
qbody(:,:) = pathDB.qbody(:,seq_nr,:);
qwingL(:,:) = pathDB.qwingL_obs(:,seq_nr,:);
qwingR(:,:) = pathDB.qwingR_obs(:,seq_nr,:);
xh = [pos qbody qwingL qwingR];

allnums = pos(:,1)*0 + 1;
startframe = find((isnan(allnums)==0), 1 );
endframe = find((isnan(allnums)==0), 1, 'last' );

% startframe = 2795 % start at trigger

frames = [startframe endframe]; % number of frames in movie

movidx = frames(1):samplerate:frames(2);
numframes = length(movidx);
solidx = [1 length(movidx)];

soldir = cd;
soln1 = zeros(length(movidx),PAR.numfly*PAR.statedim);
SOLN = zeros(length(frames(1):frames(2)),PAR.numfly*PAR.statedim);
if interptag == 1
    for i=1:length(movidx)
        soln1(i,:) = xh(movidx(i),:);
    end
end

% % remove translation
% soln1(:,1:3)=0;
% 
% %Now, perform interpolation to calculate the state at in between frames
% for k = 1:size(soln1,2)
%     SOLN(:,k) = interp1(movidx,soln1(:,k),frames(1):frames(2),'spline')';
% end
% 
% no inperp
SOLN = soln1;

%% Calculate the region of fly motion
minx = min(SOLN(:,1)) - buffer;
maxx = max(SOLN(:,1)) + buffer;
miny = min(SOLN(:,2)) - buffer;
maxy = max(SOLN(:,2)) + buffer;
minz = min(SOLN(:,3)) - buffer;
maxz = max(SOLN(:,3)) + buffer;

cameraAnchor = [  mean([minx maxx])  mean([miny maxy])  mean([minz maxz]) ] ;

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


%% setup fig
% close all
figure(fignum);
set(fignum,'units','pixels','position',[0 0 1024 1024].*plotscale);

i=1;
clear flymodQ
[x,y,z] = flymodQ(SOLN(i,:),PAR.params,PAR);
for j = 1:length(x);
    PAR.modsample(j) = size(x{j},1);
end

for k = 1:length(x);
%         surf(x{k},y{k},z{k},'facecolor','r','edgecolor','k','facelighting','phong');
%         surf(x{k},y{k},z{k},'facecolor',[.95 .95 .95],'edgecolor','k','facelighting','phong');
        surf(x{k},y{k},z{k},'facecolor',[.95 .95 .95],'edgecolor',[.25 .25 .25],'facelighting','phong');
        hold on;
end
    
axis tight
axis equal
grid off
axis off
camlight headlight
lighting phong


    set(gca,'xlim',[minx maxx],'ylim',[miny maxy],'zlim',[minz maxz]);
%     set(gca,'xlim',[minx maxx],'ylim',[miny maxy],'zlim',[minz maxz],...
%         'zdir','reverse','ydir','reverse');
%     axis vis3d
    set(fignum,'color','w');
    
view(3)

    set(gca,'XTickLabel',[])
    set(gca,'YTickLabel',[])
    set(gca,'ZTickLabel',[])

axis vis3d ;                   % Lock the aspect ration for rotation.
for i= 1:size(SOLN,1)
%     kk = frames(1)+(i-1);

    clear flymodQ
    [x,y,z] = flymodQ(SOLN(i,:),PAR.params,PAR);
    for j = 1:length(x);
        PAR.modsample(j) = size(x{j},1);
    end

    
%     close all
%  figure('units','pixels','position',[0 0 1024 1024].*plotscale);

figure(fignum);
f=gcf;

% clf;
cla
% hold off
for k = 1:length(x);
%         surf(x{k},y{k},z{k},'facecolor','r','edgecolor','k','facelighting','phong');
%         surf(x{k},y{k},z{k},'facecolor',[.95 .95 .95],'edgecolor','k','facelighting','phong');
        surf(x{k},y{k},z{k},'facecolor',[.95 .95 .95],'edgecolor',[.25 .25 .25],'facelighting','phong');
        hold on;
    end
    
 camtarget(gca,[cameraAnchor])  % Point the camera at the anchor point.
%  view((rot1(n) - parallaxAngle),rot2(n)) ; % Rotate the figre to the correct angle.

    view(0-parallaxAngle,0)
 camlight headlight
% camlight right
lighting phong
%     pause(.001)
% tic

    view(0-parallaxAngle,0)
%     saveas(gca,['3DrenderLeft.tif'])
    print(f,'3DrenderLeft.tif','-dpng',['-r',num2str(rez)],'-opengl') %save file 
    saveas(gca,[soldir '/3DrenderLeft' sprintf(['%0',num2str(6),'d'], movidx(i)) '.fig'])
    
    view(0+parallaxAngle,0)
%     saveas(gca,['3DrenderRight.tif'])
    print(f,'3DrenderRight.tif','-dpng',['-r',num2str(rez)],'-opengl') %save file 


 leftEyeImage  = imread(['3DrenderLeft.tif'])  ; % Load the left eye image.
 rightEyeImage = imread(['3DrenderRight.tif']) ; % Load the right eye image.
 
 leftEyeImageCrop = leftEyeImage(cuty+1:end-cuty,cutx+1:end-cutx,:);
 imwrite(leftEyeImageCrop,['3DrenderLeft_crop_0deg_fr' sprintf(['%0',num2str(6),'d'], movidx(i)) '.tif'])
 
 leftEyeImage(:,:,2:3) = 0 ;               % Removes green and blue from the left eye image.
 rightEyeImage(:,:,1)  = 0 ;               % Removes red from the right eye image.
 anaglyph = leftEyeImage + rightEyeImage ; % Combines the two to produce the finished anaglyph.
 
 anaglyph = anaglyph(cuty+1:end-cuty,cutx+1:end-cutx,:);
 
 imwrite(anaglyph,['3DrenderAnaglyph_crop_0deg_fr' sprintf(['%0',num2str(6),'d'], movidx(i)) '.tif'])

end
% 
% cd(movdir)
% mov = close(mov);
% cd(datadir)
cd ..
end
