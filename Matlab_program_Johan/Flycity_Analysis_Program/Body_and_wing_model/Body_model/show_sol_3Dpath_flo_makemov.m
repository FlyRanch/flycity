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

% addpath('/home/florian/Dropbox/WORK/flytracker/flytracker');
% addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/mex/');
% addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/core/');
% addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/results/');

addpath('D:/Flytracker_Ebraheem/Fly_tracker/');
addpath('D:/Flytracker_Ebraheem/Fly_tracker/flytracker/');
addpath('D:/Flytracker_Ebraheem/Fly_tracker/flytracker/mex/');
addpath('D:/Flytracker_Ebraheem/Fly_tracker/flytracker/core/');
addpath('D:/Flytracker_Ebraheem/Fly_tracker/flytracker/results/');
addpath('D:/Flytracker_Ebraheem/Fly_tracker/flytracker/flytracks/');

clear;close all;clc

datadir = cd;
cd ..
movdir = cd;

movfilename = '3Drender_startattriggerpoint_2500fps.avi';
mov = avifile(movfilename, 'compression', 'none');
%mov = VideoWriter(movfilename);

cd(datadir)

PAR.videopath = '../FlyVideo/';
PAR.filetag = 'exp101000';

PAR.solutionpath = [PAR.videopath '/solutions/'];
PAR.stub = [PAR.filetag];
fignum = 1;

%changes the scaling of the final figure, as percentage of original image
%resolution to make fit on screen.
plotscale = .8;


% --------------------------------
% Load the ManualFit Parameters
% load([PAR.solutionpath 'fly_' PAR.stub '/' 'ManualFit_' PAR.stub]);
% load('ManualFit_frame');

% load('D:/Flytracker_Ebraheem/flytracker/flytracks/ManualFit_flytracks.mat');
load('ManualFit_flytracks');
% 
% % load('ManualFit_exp095000');
% %load('C:/Users/Johan/Desktop/Flytracker/Fly_tracker/flytracks/ManualFit_flytracks.mat')
% % 
% clmexfiles = dir(['D:/Flytracker_Ebraheem/Fly_tracker/flytracker/mex/' '*.c'])
% % mex D:/Flytracker_Ebraheem/Fly_tracker/flytracker/mex/*.c
% % % mex D:/Flytracker_Ebraheem/Fly_tracker/flytracker/mex/dbasis.c
% % % mex D:/Flytracker_Ebraheem/Fly_tracker/flytracker/mex/xformq_surf.c
% % % mex D:/Flytracker_Ebraheem/Fly_tracker/flytracker/mex/*.c
% 
% mex D:/Flytracker_Ebraheem/Fly_tracker/flytracker/mex/*.c

%----------------------------------------------------
%load video data into buffer
%----------------------------------------------------
allfiles=dir('*.mat')
for n=2:length(allfiles)
    allnums(n-1) = str2num(allfiles(n).name(4:end-4));
end
startframe = min(allnums);
endframe = max(allnums);
samplerate = 1;

startframe = 2794
% endframe = 2880
samplerate = 1;

frames = [startframe endframe]; % number of frames in movie

movidx = frames(1):samplerate:frames(2);
numframes = length(movidx);



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


%Calculate the postions of the left and right wing hinges w.r.t. the bodyframe

BL = ManualFit.params.bodyscale*(ManualFit.params.bodylen+ManualFit.params.headlen);
RJTrans = BL.*([0.2021 0.1055 -0.1477]);
LJTrans = BL.*([0.2021 -0.1055 -0.1477]);




%% Calculate the region of fly motion
buffer = 1;
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

M2 = moviein(size(SOLN));
for i= 1:size(SOLN,1)
    kk = frames(1)+(i-1);

    clear flymodQ
    [x,y,z] = flymodQ(SOLN(i,:),PAR.params,PAR);
    for j = 1:length(x);
        PAR.modsample(j) = size(x{j},1);
    end

    %added code to plot position and bodyframe-----------------------------
     
    position = SOLN(i,1:3);
    
    temp1 = quat2matNEW(SOLN(i,4:7));   %Direct cosine matrix body
    temp2 = quat2matNEW(SOLN(i,8:11));  %Direct cosine matrix left wing
    temp3 = quat2matNEW(SOLN(i,12:15)); %Direct cosine matrix right wing
    
    attitude_body = [temp1(:,1)', temp1(:,2)', temp1(:,3)']; % Attitude body

    temp4 = temp1 * temp2;  %
    temp5 = temp1 * temp3;
    
    attitude_Lwing = [temp4(:,1)', temp4(:,2)', temp4(:,3)']; % Attitude left wing
    attitude_Rwing = [temp5(:,1)', temp5(:,2)', temp5(:,3)']; % Attitude right wing
    
    Lhinge = temp1 * LJTrans'; % Position left hinge
    Rhinge = temp1 * RJTrans'; % Position right hinge
    
    Rmidwing = Rhinge + temp1*temp3*[0; 2; 0];
    
%     [body_yaw_Mat(i,1) body_pitch_Mat(i,1) body_roll_Mat(i,1)] = quat2angle([SOLN(i,7) SOLN(i,4) SOLN(i,5) SOLN(i,6)]);
    
    clear temp1 temp2 temp3 temp4 temp5
    
%     %----------------------------------------------------------------------
%     
%     %calculate body yaw, pitch and roll -----------------------------------
%     
%     body_yaw(i,1)=atan2(attitude_body(2),attitude_body(1));
%     body_pitch(i,1)=atan2(attitude_body(3),sqrt(attitude_body(1)^2+attitude_body(2)^2));
%     body_roll(i,1)=acos(sqrt(attitude_body(4)^2+attitude_body(5)^2)/sqrt(attitude_body(4)^2+attitude_body(5)^2+attitude_body(6)^2));
%     body_roll2(i,1)=atan2(sqrt(attitude_body(6)^2),1);
%     %body_roll(i,1)=atan2(attitude_body
    
    
    
    
    %----------------------------------------------------------------------
    
    figure(fignum);
    clf;
    
    for k = 1:length(x);
        surf(x{k},y{k},z{k},'facecolor','b','edgecolor','k','facelighting','phong');
        hold on;
    end

    set(gca,'xlim',[minx maxx],'ylim',[miny maxy],'zlim',[minz maxz]);
%     set(gca,'xlim',[minx maxx],'ylim',[miny maxy],'zlim',[minz maxz],...
%         'zdir','reverse','ydir','reverse');
%     axis vis3d
    
    figure(fignum);
    set(fignum,'color','w');
    %title([filetag ', Frame ' num2str(kk) ]);%', x = ' num2str(soln1(i,:))])
    
    %added code to plot attitude-------------------------------------------
    
    % Body attitude (located at body center)
    figure(fignum);
    % x-axis body
    quiver3(position(1),position(2),position(3),attitude_body(1),attitude_body(2),attitude_body(3),5,'r')
    % y-axis body
    quiver3(position(1),position(2),position(3),attitude_body(4),attitude_body(5),attitude_body(6),5,'r') 
    % z-axis body
    quiver3(position(1),position(2),position(3),attitude_body(7),attitude_body(8),attitude_body(9),5,'r')
    hold on;
    
    
    % Left wing attitude (located at left wing hinge)
    figure(fignum);
    % x-axis left wing
    quiver3(position(1)+Lhinge(1,1),position(2)+Lhinge(2,1),position(3)+Lhinge(3,1),attitude_Lwing(1),attitude_Lwing(2),attitude_Lwing(3),5,'y')
    % y-axis left wing
    quiver3(position(1)+Lhinge(1,1),position(2)+Lhinge(2,1),position(3)+Lhinge(3,1),attitude_Lwing(4),attitude_Lwing(5),attitude_Lwing(6),5,'y')
    % z-axis left wing
    quiver3(position(1)+Lhinge(1,1),position(2)+Lhinge(2,1),position(3)+Lhinge(3,1),attitude_Lwing(7),attitude_Lwing(8),attitude_Lwing(9),5,'y')
    hold on;
   
   % Right wing attitude (located at right wing hinge)
   figure(fignum);
   % x-axis right wing
   quiver3(position(1)+Rhinge(1,1),position(2)+Rhinge(2,1),position(3)+Rhinge(3,1),attitude_Rwing(1),attitude_Rwing(2),attitude_Rwing(3),5,'g')
   % y-axis right wing
   quiver3(position(1)+Rhinge(1,1),position(2)+Rhinge(2,1),position(3)+Rhinge(3,1),attitude_Rwing(4),attitude_Rwing(5),attitude_Rwing(6),5,'g')
   % z-axis right wing
   quiver3(position(1)+Rhinge(1,1),position(2)+Rhinge(2,1),position(3)+Rhinge(3,1),attitude_Rwing(7),attitude_Rwing(8),attitude_Rwing(9),5,'g')
   hold on;
    
   fly_position(i,:)= [position(1) position(2) position(3)];
    
    clear position attitude_body attitude_Lwing attitude_Rwing
    %----------------------------------------------------------------------

    M2(i) = getframe(fignum);
%     axis equal
%     M2(kk) = getframe(fignum);
cd(movdir)
    mov = addframe(mov,gcf);
cd(datadir)
end

% figure()
% plot(1:1:length(body_yaw),radtodeg(body_yaw),1:1:length(body_yaw),radtodeg(body_yaw_Mat))
% title('plot yaw angle of body')
% xlabel('iterations')
% ylabel('yaw angle (deg)')
% legend('own definition','Matlab definition')
% 
% figure()
% plot(1:1:length(body_pitch),radtodeg(body_pitch),1:1:length(body_pitch),radtodeg(body_pitch_Mat))
% title('plot pitch angle of body')
% xlabel('iterations')
% ylabel('pitch angle (deg)')
% legend('own definition','Matlab definition')
% 
% figure()
% plot(1:1:length(body_roll),radtodeg(body_roll),1:1:length(body_roll),radtodeg(body_roll_Mat))
% title('plot roll angle of body')
% xlabel('iterations')
% ylabel('roll angle (deg)')
% legend('own definition','Matlab definition')
% 
% figure()
% plot(1:1:length(body_roll),radtodeg(body_roll2))
% title('plot roll angle of body, own definition 2')
% xlabel('iterations')
% ylabel('roll angle (deg)')
% 
% figure()
% plot3(fly_position(:,1), fly_position(:,2), fly_position(:,3),SOLN(:,1),SOLN(:,2),SOLN(:,3))
% title('position body center')
% xlabel('x')
% ylabel('y')
% zlabel('z')

% pitch_FFT = fft(body_pitch,length(body_pitch)); %Fast Fourrier Transform of body pitch angle
% Power_dens = pitch_FFT.*conj(pitch_FFT)/length(body_pitch); %compute power density matrix with conjugate
% freq = (2500/length(body_pitch))*(1:length(body_pitch)); %Frequencies measurable in the domain
% figure()
% plot(freq,Power_dens)


cd(movdir)
mov = close(mov);
cd(datadir)
