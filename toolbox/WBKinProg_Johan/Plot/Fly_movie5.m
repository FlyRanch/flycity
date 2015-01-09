function Fly_movie5(settings,pathDB, i)


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

    movfilename = ['Moving_body_frame_movie_' char(settings.sequence_names(i)) '.mj2'];
    mov = VideoWriter(movfilename, 'Motion JPEG AVI');
    mov.FrameRate = 5;
    mov.Quality = 100;

    
    cd(datadir)

    fignum = i;

    %changes the scaling of the final figure, as percentage of original image
    %resolution to make fit on screen.
    plotscale = .8;
    
    samplerate = 1;
    
    if pathDB.L_wingbeat_loc(1,1,i) > pathDB.L_wingbeat_loc(1,2,i)
        end_L = find(isnan(pathDB.L_wingbeat_loc(:,1,i))==0, 1 ,'last')-1;
    else
        end_L = find(isnan(pathDB.L_wingbeat_loc(:,1,i))==0, 1 ,'last')-2;
    end
    
    startframe = 1;
    endframe = end_L;

    frames = [startframe endframe]; % number of frames in movie

    movidx = frames(1):samplerate:frames(2);
    numframes = length(movidx);
    
    


    
    %% Calculate the region of fly motion

    minx = -5;
    maxx = 5;
    miny = -5;
    maxy = 5;
    minz = -5;
    maxz = 5;
    
    
    %%
    
    tsteps = 1;
    dummy = [];
    torder = 1;

    M1dum = [];
    M2dum = [];
    
    T = 0;
  

    
    
    open(mov);
    
    set(gca,'nextplot','replacechildren');
    set(gcf,'Renderer','zbuffer');

    
    for m= 1:end_L
        
        kk = frames(1)+(m-1);


        hfig = figure(fignum);
        set(hfig, 'position', [100 100 1800 900])
        clf;
        
        phi = pathDB.phi_body_mean(m,1,i);
        theta = pathDB.theta_body_mean(m,1,i);
        xsi = pathDB.xsi_body_mean(m,1,i);

            
            subplot(2,3,1); Trajectory_movie_plot(settings,pathDB,i, m);
            axis([-10 30 0 30 10 40])
            xlabel('x')
            ylabel('y')
            zlabel('z')
            hold on;
            subplot(2,3,2); Sphere_plot_wing_kinematics(settings, pathDB, i, m, pi, theta, 0,1,0);
            axis([minx maxx miny maxy minz maxz])
            xlabel('x')
            ylabel('y')
            zlabel('z')
            view(180,0)
            hold on;
            subplot(2,3,3); Sphere_plot_wing_kinematics(settings, pathDB, i, m, pi, 0, xsi,1,1);
            axis([minx maxx miny maxy minz maxz])
            xlabel('x')
            ylabel('y')
            zlabel('z')
            view(0,90)
            hold on;
            subplot(2,3,4); Sphere_plot_wing_kinematics(settings, pathDB, i, m, phi, theta, xsi,1,1);
            axis([minx maxx miny maxy minz maxz])
            xlabel('x')
            ylabel('y')
            zlabel('z')
            view(45,30)
            hold on;
            subplot(2,3,5); Sphere_plot_wing_kinematics(settings, pathDB, i, m, pi, theta, 0,0,1);
            axis([minx maxx miny maxy minz maxz])
            xlabel('x')
            ylabel('y')
            zlabel('z')
            view(0,0)
            hold on;
            subplot(2,3,6); Sphere_plot_wing_kinematics(settings, pathDB, i, m, phi, 0, 0,1,1);
            axis([minx maxx miny maxy minz maxz])
            xlabel('x')
            ylabel('y')
            zlabel('z')
            view(90,0)
            hold on;

      
        clear a b Lwt_x_down Lwt_z_down Lwt_x_up Lwt_z_up L_wingtip_x_down L_wingtip_z_down L_wingtip_x_up L_wingtip_z_up

        set(gca,'xlim',[minx maxx],'ylim',[miny maxy],'zlim',[minz maxz]);


        figure(fignum);
        set(fignum,'color','w');
        
        frame = getframe(fignum);
        
        cd(movdir)
        writeVideo(mov,frame);
        writeVideo(mov,frame);
        writeVideo(mov,frame);
        writeVideo(mov,frame);
        writeVideo(mov,frame);
        cd(datadir)

        
    end
    
    cd(movdir)
    close(mov);
    cd(datadir)
    
    clear ManualFit
    
    
    end
    
    clear all
    


end