function Force_3D_stationary(settings,pathDB,seq_nr,save_on_off)


    % Make a movie from the average force during the up and downstroke and
    % the arms on which they act:
    
    if save_on_off == 1
        
    close all;clc
    
    i = seq_nr;
    
        
    cd(char(settings.path_names(1)))
    
    datadir = cd;

    cd([char(settings.plot_folders(4)) '/' char(settings.sequence_names(i))]);
    
    movdir = cd;

    movfilename = ['Force_movie_3D_stationary_' char(settings.sequence_names(i)) '.mj2'];
    mov = VideoWriter(movfilename, 'Motion JPEG AVI');
    mov.FrameRate = 5;
    mov.Quality = 100;

    
    cd(datadir)

    fignum = i;

    %changes the scaling of the final figure, as percentage of original image
    %resolution to make fit on screen.
    plotscale = .8;
    
    samplerate = 1;
    
    end_L = find(isnan(pathDB.wingbeat_time(:,1,i))==0,1,'last');
    
    startframe = 1;
    endframe = end_L;

    frames = [startframe endframe]; % number of frames in movie

    movidx = frames(1):samplerate:frames(2);
    numframes = length(movidx);
    
    
    
    
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
    
    
    phi = pi;
    theta = 0;
    xsi = 0;
    
    
    for m= 1:end_L
        
        kk = frames(1)+(m-1);


        hfig = figure(fignum);
        set(hfig, 'position', [100 100 1800 900])
        clf;

        
            subplot(2,3,1); Trajectory_movie_plot(settings,pathDB,i, m);
            axis([-10 30 0 30 10 40])
            xlabel('x')
            ylabel('y')
            zlabel('z')
            hold on;
            subplot(2,3,2); Force_plot4(settings, pathDB, i, m, phi, theta, xsi,1,0);
            xlabel('x')
            ylabel('y')
            zlabel('z')
            view(180,0)
            hold on;
            subplot(2,3,3); Force_plot4(settings, pathDB, i, m, phi, theta, xsi,1,1);
            xlabel('x')
            ylabel('y')
            zlabel('z')
            view(0,90)
            hold on;
            subplot(2,3,4); Force_plot4(settings, pathDB, i, m, phi, theta, xsi,1,1);
            xlabel('x')
            ylabel('y')
            zlabel('z')
            view(135,30)
            hold on;
            subplot(2,3,5); Force_plot4(settings, pathDB, i, m, phi, theta, xsi,0,1);
            xlabel('x')
            ylabel('y')
            zlabel('z')
            view(0,0)
            hold on;
            subplot(2,3,6); Force_plot4(settings, pathDB, i, m, phi, theta, xsi,1,1);
            xlabel('x')
            ylabel('y')
            zlabel('z')
            view(90,0)
            hold on;


        figure(fignum);
        set(fignum,'color','w');
        
        frame = getframe(fignum);
        
        cd(movdir)
        writeVideo(mov,frame);
        cd(datadir)

        % Save plot


        saveas(fignum, [char(settings.plot_folders(3)) '/' char(settings.sequence_names(seq_nr)) '/Force_3D_stationary/' int2str(m) ], 'fig')

        saveas(fignum, [char(settings.plot_folders(3)) '/' char(settings.sequence_names(seq_nr)) '/Force_3D_stationary/' int2str(m) ], 'png')
        
    end
    
    cd(movdir)
    close(mov);
    cd(datadir)
    
    clear ManualFit
    
    
    end
    
    clear all


end
