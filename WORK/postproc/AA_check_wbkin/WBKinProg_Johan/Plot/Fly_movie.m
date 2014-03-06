function Fly_movie(settings,pathDB, i)


    % Make a movie of the body orientation and the wingkinematics per
    % wingbeat. Can only create movies when save_on_off == 1:
    
    if settings.save_on_off == 1
    


    close all;clc
    
    addpath(char(settings.path_names(7)))
    addpath(char(settings.path_names(8)))
        
    cd(char(settings.path_names(1)))
    
    datadir = cd;

    cd([char(settings.plot_folders(4)) '/' char(settings.sequence_names(i))]);
    
    movdir = cd;

%     movfilename = ['Trajectory_movie_' char(settings.sequence_names(i)) '.avi'];
%     mov = avifile(movfilename, 'compression', 'none');
%     mov.fps = 30;
    
    movfilename = ['Trajectory_movie_' char(settings.sequence_names(i)) '.mj2'];
    mov = VideoWriter(movfilename, 'Motion JPEG AVI');
    mov.FrameRate = 30;
    mov.Quality = 100;

    cd(datadir)

    PAR.videopath = '../FlyVideo/';
    PAR.filetag = 'exp101000';

    PAR.solutionpath = [PAR.videopath '/solutions/'];
    PAR.stub = [PAR.filetag];
    
    fignum = i;

    %changes the scaling of the final figure, as percentage of original image
    %resolution to make fit on screen.
    plotscale = .8;
    
    samplerate = 1;
    
    % start and stop point for the measurements
    startframe = find(isnan(pathDB.x(:,i))==0, 1 );
    endframe = find(isnan(pathDB.x(:,i))==0, 1, 'last' );

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
    
    
    % Load Manualfit
    
    cd(char(settings.sequence_names(i)))
    
    cd('flytracks')
    
    cd
    
    load('ManualFit_flytracks');
    


    % Assign model parameters
    PAR.params = ManualFit.params;
    PAR.DLT = ManualFit.DLT;
    PAR.cam = ManualFit.cam;
    
    
    
    soln1 = zeros(length(movidx),PAR.numfly*PAR.statedim);
    SOLN = zeros(length(frames(1):frames(2)),PAR.numfly*PAR.statedim);
    if interptag == 1
        for j=1:length(movidx)
            soln1(j,:) = [pathDB.x_filt(startframe-1+j,i); pathDB.y_filt(startframe-1+j,i); pathDB.z_filt(startframe-1+j,i); ...
                          pathDB.qb1_filt(startframe-1+j,i); pathDB.qb2_filt(startframe-1+j,i); pathDB.qb3_filt(startframe-1+j,i); pathDB.qb4_filt(startframe-1+j,i); ...
                          pathDB.qL1_filt2(startframe-1+j,i); pathDB.qL2_filt2(startframe-1+j,i); pathDB.qL3_filt2(startframe-1+j,i); pathDB.qL4_filt2(startframe-1+j,i); ...
                          pathDB.qR1_filt2(startframe-1+j,i); pathDB.qR2_filt2(startframe-1+j,i); pathDB.qR3_filt2(startframe-1+j,i); pathDB.qR4_filt2(startframe-1+j,i)];
                      
            % Noisy data
%             soln1(j,:) = [pathDB.x(startframe-1+j,i); pathDB.y(startframe-1+j,i); pathDB.z(startframe-1+j,i); ...
%                           pathDB.qb1(startframe-1+j,i); pathDB.qb2(startframe-1+j,i); pathDB.qb3(startframe-1+j,i); pathDB.qb4(startframe-1+j,i); ...
%                           pathDB.qL1(startframe-1+j,i); pathDB.qL2(startframe-1+j,i); pathDB.qL3(startframe-1+j,i); pathDB.qL4(startframe-1+j,i); ...
%                           pathDB.qR1(startframe-1+j,i); pathDB.qR2(startframe-1+j,i); pathDB.qR3(startframe-1+j,i); pathDB.qR4(startframe-1+j,i)];
        end
    end
    


    SOLN = soln1;

    
    %% Calculate the region of fly motion
    buffer = 1;
    minx = min(SOLN(:,1)) - buffer;
    maxx = max(SOLN(:,1)) + buffer;
    miny = min(SOLN(:,2)) - buffer;
    maxy = max(SOLN(:,2)) + buffer;
    minz = min(SOLN(:,3)) - buffer;
    maxz = max(SOLN(:,3)) + buffer;

%     buffer = 25;
%     minx = mean(SOLN(:,1))-buffer;
%     maxx = mean(SOLN(:,1))+buffer;
%     miny = mean(SOLN(:,2))-buffer;
%     maxy = mean(SOLN(:,2))+buffer;
%     minz = mean(SOLN(:,3))-buffer;
%     maxz = mean(SOLN(:,3))+buffer;
    
    
    %%
    
    tsteps = 1;
    dummy = [];
    torder = 1;

    M1dum = [];
    M2dum = [];
    
%     T = 0;%linspace(0,12,60);
% 
%     M2 = moviein(size(SOLN));

    open(mov);
    
    set(gca,'nextplot','replacechildren');
    set(gcf,'Renderer','zbuffer');
    
    for m= 1:size(SOLN,1)
        
        kk = frames(1)+(m-1);

        clear flymodQ
        [x,y,z] = flymodQ(SOLN(m,:),PAR.params,PAR);
        for j = 1:length(x);
            PAR.modsample(j) = size(x{j},1);
        end

        hfig = figure(fignum);
        set(hfig, 'position', [100 100 1800 900])
        clf;

%         for k = 1:length(x);
%             subplot(2,2,1); surf(x{k},y{k},z{k},'facecolor','b','edgecolor','k','facelighting','phong');
%             hold on;
%             plot3(SOLN(1:m,1),SOLN(1:m,2),SOLN(1:m,3),'r');
%             axis equal
%             axis([minx maxx miny maxy minz maxz])
%             xlabel('x')
%             ylabel('y')
%             zlabel('z')
%             hold on;
%             subplot(2,2,2); surf(x{k},y{k},z{k},'facecolor','b','edgecolor','k','facelighting','phong');
%             hold on;
%             plot3(SOLN(1:m,1),SOLN(1:m,2),SOLN(1:m,3),'r');
%             axis equal
%             axis([minx maxx miny maxy minz maxz])
%             xlabel('x')
%             ylabel('y')
%             zlabel('z')
%             view(0,90)
%             hold on;
%             subplot(2,2,3); surf(x{k},y{k},z{k},'facecolor','b','edgecolor','k','facelighting','phong');
%             hold on;
%             plot3(SOLN(1:m,1),SOLN(1:m,2),SOLN(1:m,3),'r');
%             axis equal
%             axis([minx maxx miny maxy minz maxz])
%             xlabel('x')
%             ylabel('y')
%             zlabel('z')
%             view(0,0)
%             hold on;
%             subplot(2,2,4); surf(x{k},y{k},z{k},'facecolor','b','edgecolor','k','facelighting','phong');
%             hold on;
%             plot3(SOLN(1:m,1),SOLN(1:m,2),SOLN(1:m,3),'r');
%             axis equal
%             axis([minx maxx miny maxy minz maxz])
%             xlabel('x')
%             ylabel('y')
%             zlabel('z')
%             view(90,0)
%             hold on;
%         end

        for k = 1:length(x);
            surf(x{k},y{k},z{k},'facecolor','b','edgecolor','k','facelighting','phong');
            hold on;
            plot3(SOLN(1:m,1),SOLN(1:m,2),SOLN(1:m,3),'r');
            axis equal
            axis([minx maxx miny maxy minz maxz])
            xlabel('x')
            ylabel('y')
            zlabel('z')
            hold on;
        end


        set(gca,'xlim',[minx maxx],'ylim',[miny maxy],'zlim',[minz maxz]);
    %     set(gca,'xlim',[minx maxx],'ylim',[miny maxy],'zlim',[minz maxz],...
    %         'zdir','reverse','ydir','reverse');
    %     axis vis3d

        figure(fignum);
        set(fignum,'color','w');
        %view(0,90)
        %title([filetag ', Frame ' num2str(kk) ]);%', x = ' num2str(soln1(i,:))])

% 
%         M2(m) = getframe(fignum);
        %     axis equal
        %     M2(kk) = getframe(fignum);
        
        frame = getframe(fignum);
        
        cd(movdir)
        writeVideo(mov,frame);
        cd(datadir)
        
%         
%         
%         cd(movdir)
%         mov = addframe(mov,gcf);
%         cd(datadir)

        
        end
    

    
%     %end
%     cd(movdir)
%     mov = close(mov);
%     cd(datadir)

    cd(movdir)
    close(mov);
    cd(datadir)



    
    clear ManualFit
    
    
    end
    
    clear all
    


end

