function Fly_movie3(settings,pathDB, i)


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

    movfilename = ['Trajectory_movie_' char(settings.sequence_names(i)) '.avi'];
    mov = avifile(movfilename, 'compression', 'none', 'fps', 5);

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
    
    %%
    
    if pathDB.L_wingbeat_loc(1,1,i) > pathDB.L_wingbeat_loc(1,2,i)
        start_L = 0;
        end_L = find(isnan(pathDB.L_wingbeat_loc(:,1,i))==0, 1 ,'last')-1;
    else
        start_L = 1;
        end_L = find(isnan(pathDB.L_wingbeat_loc(:,1,i))==0, 1 ,'last')-2;
    end
    
    
    startframe = 1+start_L;
    endframe = end_L;

    frames = [startframe endframe]; % number of frames in movie

    movidx = frames(1):samplerate:frames(2);
    numframes = length(movidx);

    %%

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
    
    theta = zeros(length(movidx),1);
    
    if interptag == 1
        for j=1:length(movidx)
            
            theta(j) = pathDB.theta_body_mean(startframe-1+j,1,i);
            
            theta_body = theta(j);
            
            soln1(j,:) = [5; 5; 5; ...
                          0; sin((theta_body)/2); 0; cos((theta_body)/2); ...
                          0; 0; 0; 1; ...
                          0; 0; 0; 1];
% 
%             soln1(j,:) = [5; 5; 5; ...
%                           0; 0; 0; 1; ...
%                           0; 0; 0; 1; ...
%                           0; 0; 0; 1];
            
        end
    end
    
    SOLN = soln1;

    
    %% Calculate the region of fly motion
    buffer = 5;
    minx = min(SOLN(:,1)) - buffer;
    maxx = max(SOLN(:,1)) + buffer;
    miny = min(SOLN(:,2)) - buffer;
    maxy = max(SOLN(:,2)) + buffer;
    minz = min(SOLN(:,3)) - buffer;
    maxz = max(SOLN(:,3)) + buffer;
    
    
    %%
    
    tsteps = 1;
    dummy = [];
    torder = 1;

    M1dum = [];
    M2dum = [];
    
    T = 0;%linspace(0,12,60);

    M2 = moviein(size(SOLN));
    
    start_meas = find(isnan(pathDB.qb1(:,i))==0, 1 );
    stop_meas = find(isnan(pathDB.qb1(:,i))==0, 1, 'last' );
    
    jlocL = pathDB.joint_pos_L(:,i);
    jlocR = pathDB.joint_pos_R(:,i);
    
    for m= 1:size(SOLN,1)
        
        kk = frames(1)+(m-1);

        clear flymodQ
        [x,y,z] = flymodQ(SOLN(m,:),PAR.params,PAR);
        for j = 1:length(x);
            PAR.modsample(j) = size(x{j},1);
        end


        a = (start_meas-1+pathDB.L_wingbeat_loc(m,2,i)):1:(start_meas-1+pathDB.L_wingbeat_loc(start_L+m,1,i));
        b = (start_meas-1+pathDB.L_wingbeat_loc(start_L+m,1,i)):1:(start_meas-1+pathDB.L_wingbeat_loc(m+1,2,i));
        c = (start_meas-1+pathDB.R_wingbeat_loc(m,2,i)):1:(start_meas-1+pathDB.R_wingbeat_loc(start_L+m,1,i));
        d = (start_meas-1+pathDB.R_wingbeat_loc(start_L+m,1,i)):1:(start_meas-1+pathDB.R_wingbeat_loc(m+1,2,i));

        figure(fignum);
        clf;
        
        k = 1;

        
        L_wingtip_x_down = cos(theta(m)).*(jlocL(1)+pathDB.Lwingtip(a,1,i))-sin(theta(m)).*(jlocL(3)+pathDB.Lwingtip(a,3,i));
        L_wingtip_z_down = sin(theta(m)).*(jlocL(1)+pathDB.Lwingtip(a,1,i))+cos(theta(m)).*(jlocL(3)+pathDB.Lwingtip(a,3,i));

        
        L_wingtip_x_up = cos(theta(m)).*(jlocL(1)+pathDB.Lwingtip(b,1,i))-sin(theta(m)).*(jlocL(3)+pathDB.Lwingtip(b,3,i));
        L_wingtip_z_up = sin(theta(m)).*(jlocL(1)+pathDB.Lwingtip(b,1,i))+cos(theta(m)).*(jlocL(3)+pathDB.Lwingtip(b,3,i));
        
        
        R_wingtip_x_down = cos(theta(m)).*(jlocR(1)+pathDB.Rwingtip(c,1,i))-sin(theta(m)).*(jlocR(3)+pathDB.Rwingtip(c,3,i));
        R_wingtip_z_down = sin(theta(m)).*(jlocR(1)+pathDB.Rwingtip(c,1,i))+cos(theta(m)).*(jlocR(3)+pathDB.Rwingtip(c,3,i));

        
        R_wingtip_x_up = cos(theta(m)).*(jlocR(1)+pathDB.Rwingtip(d,1,i))-sin(theta(m)).*(jlocR(3)+pathDB.Rwingtip(d,3,i));
        R_wingtip_z_up = sin(theta(m)).*(jlocR(1)+pathDB.Rwingtip(d,1,i))+cos(theta(m)).*(jlocR(3)+pathDB.Rwingtip(d,3,i));

        
        

            subplot(2,2,1); surf(x{k},y{k},z{k},'facecolor','b','edgecolor','k','facelighting','phong');
            hold on
            plot3(5+L_wingtip_x_down,10.*ones(length(a),1),5-L_wingtip_z_down,'--r')
            plot3(5+L_wingtip_x_up,10.*ones(length(b),1),5-L_wingtip_z_up,'--r')
            axis([minx maxx miny maxy minz maxz])
            xlabel('x')
            ylabel('y')
            zlabel('z')
            view(180,0)
            hold on;
            subplot(2,2,2); surf(x{k},y{k},z{k},'facecolor','b','edgecolor','k','facelighting','phong');
            hold on;
            plot3(SOLN(1:m,1),SOLN(1:m,2),SOLN(1:m,3),'r');
            axis([minx maxx miny maxy minz maxz])
            xlabel('x')
            ylabel('y')
            zlabel('z')
            view(0,90)
            hold on;
            subplot(2,2,3); surf(x{k},y{k},z{k},'facecolor','b','edgecolor','k','facelighting','phong');
            hold on;
            plot3(5+R_wingtip_x_down,zeros(length(c),1),5-R_wingtip_z_down,'--r')
            plot3(5+R_wingtip_x_up,zeros(length(d),1),5-R_wingtip_z_up,'--r')
            axis([minx maxx miny maxy minz maxz])
            xlabel('x')
            ylabel('y')
            zlabel('z')
            view(0,0)
            hold on;
            subplot(2,2,4); surf(x{k},y{k},z{k},'facecolor','b','edgecolor','k','facelighting','phong');
            hold on;
            plot3(SOLN(1:m,1),SOLN(1:m,2),SOLN(1:m,3),'r');
            axis([minx maxx miny maxy minz maxz])
            xlabel('x')
            ylabel('y')
            zlabel('z')
            view(90,0)
            hold on;

      
        clear a b Lwt_x_down Lwt_z_down Lwt_x_up Lwt_z_up L_wingtip_x_down L_wingtip_z_down L_wingtip_x_up L_wingtip_z_up

        set(gca,'xlim',[minx maxx],'ylim',[miny maxy],'zlim',[minz maxz]);
    %     set(gca,'xlim',[minx maxx],'ylim',[miny maxy],'zlim',[minz maxz],...
    %         'zdir','reverse','ydir','reverse');
    %     axis vis3d

        figure(fignum);
        set(fignum,'color','w');
        %view(0,90)
        %title([filetag ', Frame ' num2str(kk) ]);%', x = ' num2str(soln1(i,:))])


        M2(m) = getframe(fignum);
        %     axis equal
        %     M2(kk) = getframe(fignum);
        cd(movdir)
        mov = addframe(mov,gcf);
        mov = addframe(mov,gcf);
        mov = addframe(mov,gcf);
        mov = addframe(mov,gcf);
        mov = addframe(mov,gcf);
        cd(datadir)

        
        end
    

    
    %end
    cd(movdir)
    mov = close(mov);
    cd(datadir)
    
    clear ManualFit
    
    
    end
    
    clear all
    


end

