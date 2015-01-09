function Trajectory_plus_orientation_movie(settings,pathDB,seq_nr,save_on_off)


    % Program that plots the trajectory of the fruit fly + a close-up of
    % the 3D orientation of the fruitfly
    
    if save_on_off == 1
        
    close all;clc
    
    i = seq_nr;
    
        
    cd(char(settings.path_names(1)))
    
    datadir = cd;

    cd([char(settings.plot_folders(4)) '/' char(settings.sequence_names(i))]);
    
    movdir = cd;

    movfilename = ['Trajectory_moving_body_frame_' char(settings.sequence_names(i)) '.mj2'];
    mov = VideoWriter(movfilename, 'Motion JPEG AVI');
    mov.FrameRate = 30;
    mov.Quality = 100;

    
    cd(datadir)

    fignum = i;

    plotscale = .8;
    
    samplerate = 1;
    
%     if pathDB.L_wingbeat_loc(1,1,i) > pathDB.L_wingbeat_loc(1,2,i)
%         end_L = find(isnan(pathDB.L_wingbeat_loc(:,1,i))==0, 1 ,'last')-1;
%     else
%         end_L = find(isnan(pathDB.L_wingbeat_loc(:,1,i))==0, 1 ,'last')-2;
%     end
%     
%     startframe = 1;
%     endframe = end_L;
% 
%     frames = [startframe endframe]; % number of frames in movie
% 
%     movidx = frames(1):samplerate:frames(2);
%     numframes = length(movidx);
    
    
    
%     %%
%     
%     tsteps = 1;
%     dummy = [];
%     torder = 1;
% 
%     M1dum = [];
%     M2dum = [];
%     
%     T = 0;
  
    
    open(mov);
    
    set(gca,'nextplot','replacechildren');
    set(gcf,'Renderer','zbuffer');
    
    % start and stop point for the measurements
    start = find(isnan(pathDB.qb1(:,i))==0, 1 );
    stop = find(isnan(pathDB.qb1(:,i))==0, 1, 'last' );

    
    for m= start:stop
        
%         kk = frames(1)+(m-1);


        hfig = figure(fignum);
        set(hfig, 'position', [100 100 1800 900])
        clf;
        
        phi = pathDB.b_roll(m,i);
        theta = pathDB.b_pitch(m,i);
        xsi = pathDB.b_yaw(m,i);

            subplot(1,2,1); Trajectory_movie_plot2(settings,pathDB,i, m)
            axis([-10 30 0 30 10 40])
            xlabel('x [mm]','FontSize',14)
            set(gca, 'XTick', [-10,0,10,20,30],'FontSize',14);
            ylabel('y [mm]','FontSize',14)
            set(gca, 'YTick', [0,10,20,30],'FontSize',14);
            zlabel('z [mm]','FontSize',14)
            set(gca, 'ZTick', [10, 20 , 30 , 40],'FontSize',14);
            title('Trajectory','FontSize',16)
            view(-45,30)
            camlight(45,70)
            hold on;
            subplot(1,2,2); Body_orientation(settings, pathDB,i, phi, theta, xsi, m);
            axis([-4 4 -4 4 -4 4 -4 4])
            xlabel('x [mm]','FontSize',14)
            set(gca, 'XTick', [-4,0,4],'FontSize',14);
            ylabel('y [mm]','FontSize',14)
            set(gca, 'YTick', [-4,0,4],'FontSize',14);
            zlabel('z [mm]','FontSize',14)
            set(gca, 'ZTick', [-4,0,4],'FontSize',14);
            view(-45,30)
            title('Orientation','FontSize',16)
            camlight(45,70)
            hold on;

      
        clear phi theta xsi

        %set(gca,'xlim',[minx maxx],'ylim',[miny maxy],'zlim',[minz maxz]);


        figure(fignum);
        set(fignum,'color','w');
        
        frame = getframe(fignum);
        
        cd(movdir)
        writeVideo(mov,frame);
        cd(datadir)

%         % Save plot
% 
% 
%         saveas(fignum, [char(settings.plot_folders(3)) '/' char(settings.sequence_names(seq_nr)) '/Wing_kinematics_3D_moving/' int2str(m) ], 'fig')
% 
%         saveas(fignum, [char(settings.plot_folders(3)) '/' char(settings.sequence_names(seq_nr)) '/Wing_kinematics_3D_moving/' int2str(m) ], 'png')
        
    end
    
    cd(movdir)
    close(mov);
    cd(datadir)
    
    clear ManualFit
    
    
    end
    
    clear all
    


end
