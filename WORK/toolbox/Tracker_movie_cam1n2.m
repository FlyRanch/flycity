% function Tracker_movie(settings, pathDB, seq_nr)


    % Create a plot of the body orientation and save as jpeg files:
    
    close all;clc
    
    samplerate = 3;
    
    % start and stop point for the measurements
    start = 2500;
    stop = 3600;

    xmin1=170;
    xmax1=620;
    dx1 = xmax1-xmin1
    ymin1=580;
    ymax1=930;
    dy1 = ymax1-ymin1
    
    xmin2=425;
    xmax2=725;
    dx2 = xmax2-xmin2
    ymin2=530;
    ymax2=880;
    dy2 = ymax2-ymin2
    
    xmin3=425;
    xmax3=725;
    dx3 = xmax3-xmin3
    ymin3=530;
    ymax3=880;
    dy3 = ymax3-ymin3
    
    
    datadir = cd;
    datadir1 = [datadir(1:end-10),'1',datadir(end-8:end)];
    datadir2 = [datadir(1:end-10),'2',datadir(end-8:end)];
    datadir3 = [datadir(1:end-10),'3',datadir(end-8:end)];
    
    cd(datadir1)
    figs1 = dir;
    figs1(1:3)=[];
       
    cd(datadir2)
    figs2 = dir;
    figs2(1:3)=[];
       
    cd(datadir3)
    figs3 = dir;
    figs3(1:3)=[];
       
    movdir = cd;
%     movfilename = [char(movdir) '.mj2'];
%     mov = VideoWriter(movfilename, 'Motion JPEG AVI');
%     mov.FrameRate = 30;
%     mov.Quality = 100;
%     open(mov);
%     
    
    
    for m = start:samplerate:stop
    
    A1 = imread([datadir1,'/',figs1(m).name]);
    
    A2 = imread([datadir2,'/',figs2(m).name]);
    
    A3 = imread([datadir3,'/',figs3(m).name]);


    A1n = im2double(imadjust(A1(ymin1:ymax1,xmin1:xmax1)));
    
    A2n = im2double(imadjust(A2(ymin2:ymax2,xmin2:xmax2)));
    
    A3n = im2double(imadjust(A3(ymin3:ymax3,xmin3:xmax3)));
   
       
%     Aall = [A1n A2n A3n];
    Aall = [A1n A2n];

    % Save the concatenated matrices as a movie
    
    imwrite(Aall, [ char(movdir) '_'  figs1(m).name]);

%     cd(movdir)
%     writeVideo(mov,A5);
%     cd(datadir)
    
    clear A1 A2 A3
    
    end
    
%     cd(movdir)
%     close(mov);
%     cd(datadir)
    
  
    
    



