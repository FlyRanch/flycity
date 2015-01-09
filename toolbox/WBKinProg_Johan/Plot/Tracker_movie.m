function Tracker_movie(settings, pathDB, seq_nr)


    % Create a plot of the body orientation and save as jpeg files:
    
    close all;clc
    
    i = seq_nr;
   
    fignum = i;
    
    % start and stop point for the measurements
    start = 3200;
    stop = 3600;

    
%     for m= start:stop
%         
%         hfig = figure(fignum);
%         set(hfig, 'position', [100 100 1200 900])
%         clf;
%         
%         phi = pathDB.b_roll(m,i);
%         theta = pathDB.b_pitch(m,i);
%         xsi = pathDB.b_yaw(m,i);
% 
% 
%         Body_orientation(settings,pathDB,i,phi,theta,xsi,m);
%         axis([-4 4 -4 4 -4 4 -4 4])
% %         xlabel('x [mm]','FontSize',14)
% %         set(gca, 'XTick', [-4,0,4],'FontSize',14);
% %         ylabel('y [mm]','FontSize',14)
% %         set(gca, 'YTick', [-4,0,4],'FontSize',14);
% %         zlabel('z [mm]','FontSize',14)
% %         set(gca, 'ZTick', [-4,0,4],'FontSize',14);
%         view(45,30)
% %         title('Orientation','FontSize',16)
%         camlight(-45,70)
%         axis off
%         hold on;
%       
%         clear phi theta xsi
%         
%         figure(fignum);
%         set(fignum,'color','w');
%         %set(fignum,'Units','pixels','Position',[500 400 600 500]);   
% 
%         f = getframe(gcf);
%         pause
%         imwrite(f.cdata, [ char(settings.path_names(9)) '/3D_model/'  int2str(m) '.tif'])
%         
%         clear f
%         
%         %saveas(gcf, [ char(settings.path_names(9)) '/3D_model/'  int2str(m) ], 'jpg')
%         
%     end
% 
%     close all

    % Open all jpeg files, load them into matrices and concatenate these
    % matrices:
    
    cd(char(settings.path_names(1)));
    
    datadir = cd;
       
    cd(char(settings.path_names(9)))
    
    movdir = cd;

    movfilename = ['Tracking_movie_' char(settings.sequence_names(i)) '.mj2'];
    mov = VideoWriter(movfilename, 'Motion JPEG AVI');
    mov.FrameRate = 30;
    mov.Quality = 100;

    
    %cd(datadir)

    fignum = i;

    plotscale = .8;
    
    samplerate = 1;
    
        
    open(mov);
    
    set(gca,'nextplot','replacechildren');
    set(gcf,'Renderer','zbuffer');
    
    
    
    for m = start:stop
    
    A1 = imread([char(settings.path_names(9)) '/3D_model/'  int2str(m) '.tif' ]);
    
%     A2 = imread([char(settings.path_names(9)) '/cam1/flyimage_00' int2str(m) '.tif' ]);
%     
%     A3 = imread([char(settings.path_names(9)) '/cam2/flyimage_00' int2str(m) '.tif' ]);
%     
%     A4 = imread([char(settings.path_names(9)) '/cam3/flyimage_00' int2str(m) '.tif' ]);

    A2 = imread([char(settings.path_names(9)) '/cam1/C001H001S000100' int2str(m) '.tif' ]);
    
    A3 = imread([char(settings.path_names(9)) '/cam2/C002H001S000100' int2str(m) '.tif' ]);
    
    A4 = imread([char(settings.path_names(9)) '/cam3/C003H001S000100' int2str(m) '.tif' ]);

    A1n = im2double(A1);
    
%     A2n = imadjust(A2(312:971,158:937),[0 1], [0.1 0.6]);
%     
%     A3n = imadjust(A3(298:957,1:780),[0 1], [0.1 0.6]);
%     
%     A4n = imadjust(A4(341:1000,200:979),[0 1], [0.1 0.6]);

    A2n = im2double(imadjust(A2(312:971,158:937)));
    
    A3n = im2double(imadjust(A3(298:957,1:780)));
    
    A4n = im2double(imadjust(A4(341:1000,200:979)));
        
%     figure, imshow(A1n), figure, imshow(A2n), figure, imshow(A3n), figure, imshow(A4n);
%     
%     pause
    
%     
%     pause

%     size(A1)
%     
%     size(A2)
%     
%     size(A3)
%     
%     size(A4)
    
    % Concatenate the matrices
    
%     A5 = [ A1(121:780,211:990,:) A2(41:700,421:1200,:); ...
%            A3(40:699,98:877,:) A4(136:795,421:1200,:)];

%     size(A1n(121:780,211:990,:))
%     
%     size(A2n(312:971,158:937,:))
%     
%     size(A3n(298:957,1:780,:))
%     
%     size(A4n(341:1000,200:979,:))
   
       
    A5 = [A1n(121:780,211:990) A2n; ...
          A3n A4n];

%     hfig = figure(fignum);
%     set(hfig, 'position', [100 100 1200 900])
%     clf;   
%        
%     figure, imshow(A5);
%     
%     pause
       
%     A5(:,:,1) = flipud(A5(:,:,1));
%     A5(:,:,2) = flipud(A5(:,:,2));
%     A5(:,:,3) = flipud(A5(:,:,3));
    

%     hfig = figure(fignum);
%     set(hfig, 'position', [100 100 1200 900])
%     clf;
%     
%     
%     image(A5)
%     axis off
%     axis image
%     
%     
    % Save the concatenated matrices as a movie
    
%     frame = getframe(imgcf);
    imwrite(A5, [ char(settings.path_names(9)) '/4views/'  int2str(m) '.tif']);

    cd(movdir)
    writeVideo(mov,A5);
    cd(datadir)
    
    clear A1 A2 A3 A4 A5
    
    end
    
    cd(movdir)
    close(mov);
    cd(datadir)
    
  
    
    

end

