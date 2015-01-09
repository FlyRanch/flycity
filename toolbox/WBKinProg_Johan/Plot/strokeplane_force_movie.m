function strokeplane_force_movie(settings,pathDB,seq_nr,save_on_off)


    % Make a movie of the strokeplane forces:
    
    if save_on_off == 1
        
    close all;clc
    
    i = seq_nr;
    
        
    cd(char(settings.path_names(1)))
    
    datadir = cd;

    cd([char(settings.plot_folders(4)) '/' char(settings.sequence_names(i))]);
    
    movdir = cd;

    movfilename = ['Wing_kinematics_2D_' char(settings.sequence_names(i)) '.mj2'];
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
    
    


    
%     %% Calculate the region of fly motion
% 
%     minx = -5;
%     maxx = 5;
%     miny = -5;
%     maxy = 5;
%     minz = -5;
%     maxz = 5;
    
    
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

        strokeplane_forces(settings,pathDB,i,m)

%         set(gca,'xlim',[minx maxx],'ylim',[miny maxy],'zlim',[minz maxz]);


        figure(fignum);
        set(fignum,'color','w');
        
        frame = getframe(fignum);
        
        cd(movdir)
        writeVideo(mov,frame);
        cd(datadir)

        % Save plot


        saveas(fignum, [char(settings.plot_folders(3)) '/' char(settings.sequence_names(seq_nr)) '/Wing_kinematics_2D/' int2str(m) ], 'fig')

        saveas(fignum, [char(settings.plot_folders(3)) '/' char(settings.sequence_names(seq_nr)) '/Wing_kinematics_2D/' int2str(m) ], 'png')
        
    end
    
    cd(movdir)
    close(mov);
    cd(datadir)
    
    clear ManualFit
    
    
    end
    
    clear all
    



end

