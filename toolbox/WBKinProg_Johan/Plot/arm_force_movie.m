function arm_force_movie(settings,pathDB,seq_nr,save_on_off)


    % Make a movie from the average force during the up and downstroke and
    % the arms on which they act:
    
    if save_on_off == 1
        
    close all;clc
    
    i = seq_nr;
    
        
    cd(char(settings.path_names(1)))
    
    datadir = cd;

    cd([char(settings.plot_folders(4)) '/' char(settings.sequence_names(i))]);
    
    movdir = cd;

    movfilename = ['Arm_force_movie_' char(settings.sequence_names(i)) '.mj2'];
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

    
    for m= 1:end_L
        
        kk = frames(1)+(m-1);


        hfig = figure(fignum);
        set(hfig, 'position', [100 100 1800 900])
        clf;

        subplot(2,2,1); Arm_force(settings,pathDB,i,m)
        view(60,45)
        subplot(2,2,2); Arm_force(settings,pathDB,i,m)
        view(0,0)
        subplot(2,2,3); Arm_force(settings,pathDB,i,m)
        view(0,90)
        subplot(2,2,4); Arm_force(settings,pathDB,i,m)
        view(90,0)


        figure(fignum);
        set(fignum,'color','w');
        
        frame = getframe(fignum);
        
        cd(movdir)
        writeVideo(mov,frame);
        cd(datadir)

        % Save plot


        saveas(fignum, [char(settings.plot_folders(3)) '/' char(settings.sequence_names(seq_nr)) '/Arm_force/' int2str(m) ], 'fig')

        saveas(fignum, [char(settings.plot_folders(3)) '/' char(settings.sequence_names(seq_nr)) '/Arm_force/' int2str(m) ], 'png')
        
    end
    
    cd(movdir)
    close(mov);
    cd(datadir)
    
    clear ManualFit
    
    
    end
    
    clear all


end

