function Fly_movie2(settings,pathDB, i)


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

    movfilename = ['Body_orientation_movie_' char(settings.sequence_names(i)) '.avi'];
    mov = avifile(movfilename, 'compression', 'none');

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
    startframe = find(isnan(pathDB.x(:,i))==0, 1 )
    endframe = find(isnan(pathDB.x(:,i))==0, 1, 'last' )

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
    
    for n = 1:1:find(isnan(pathDB.L_wingbeat_loc(:,1,i))==0, 1 ,'last')
        Lwingbeat_loc = pathDB.L_wingbeat_loc(n,:,i);
        Rwingbeat_loc = pathDB.R_wingbeat_loc(n,:,i);
    end
    
    for n = 1:1:find(isnan(pathDB.L_wingbeat_loc(:,1,i))==0, 1 ,'last')
        Lwingtip = pathDB.Lwingtip(n,:,i);
        Rwingtip = pathDB.Rwingtip(n,:,i);
    end
    
    

    if Lwingbeat_loc(1,1) > Lwingbeat_loc(1,2)
        start_L = 1;
        end_L = size(Lwingbeat_loc,1)-2;
    else
        start_L = 0;
        end_L = size(Lwingbeat_loc,1)-1;
    end

    if Rwingbeat_loc(1,1) > Rwingbeat_loc(1,2)
        start_R = 1;
        end_R = size(Rwingbeat_loc,1)-2;
    else
        start_R = 0;
        end_R = size(Rwingbeat_loc,1)-1;
    end
    
    
    
    soln1 = zeros(length(movidx),PAR.numfly*PAR.statedim);
    SOLN = zeros(length(frames(1):frames(2)),PAR.numfly*PAR.statedim);
    if interptag == 1
        for j=1:end_L
            
            theta = pathDB.b_pitch(startframe-1+j,i)-pi/2;
            
            soln1(j,:) = [5; 5; 5; ...
                          0; sin(theta/2); 0; cos(theta/2); ...
                          0; 0; 0; 1; ...
                          0; 0; 0; 1];
                      
            % Noisy data
%             soln1(j,:) = [pathDB.x(startframe-1+j,i); pathDB.y(startframe-1+j,i); pathDB.z(startframe-1+j,i); ...
%                           pathDB.qb1(startframe-1+j,i); pathDB.qb2(startframe-1+j,i); pathDB.qb3(startframe-1+j,i); pathDB.qb4(startframe-1+j,i); ...
%                           pathDB.qL1(startframe-1+j,i); pathDB.qL2(startframe-1+j,i); pathDB.qL3(startframe-1+j,i); pathDB.qL4(startframe-1+j,i); ...
%                           pathDB.qR1(startframe-1+j,i); pathDB.qR2(startframe-1+j,i); pathDB.qR3(startframe-1+j,i); pathDB.qR4(startframe-1+j,i)];

           clear theta
        end
    end
    



    SOLN = soln1;

    
    %% Calculate the region of fly motion
    buffer = 1;
%     minx = min(SOLN(:,1)) - buffer;
%     maxx = max(SOLN(:,1)) + buffer;
%     miny = min(SOLN(:,2)) - buffer;
%     maxy = max(SOLN(:,2)) + buffer;
%     minz = min(SOLN(:,3)) - buffer;
%     maxz = max(SOLN(:,3)) + buffer;
    
    minx = 0;
    maxx = 10;
    miny = 0;
    maxy = 10;
    minz = 0;
    maxz = 10;
    
    %%
    
    tsteps = 1;
    dummy = [];
    torder = 1;

    M1dum = [];
    M2dum = [];
    
    T = 0;%linspace(0,12,60);

    M2 = moviein(size(SOLN));
    
    
   
        
        
    
    
    for m= 1:end_L
        
        kk = frames(1)+(m-1);

        clear flymodQ
        [x,y,z] = flymodQ(SOLN(m,:),PAR.params,PAR);
        
%         for j = 1:length(x);
%             PAR.modsample(j) = size(x{j},1);
%         end
        
        
        PAR.modsample(1) = size(x{1},1);
        

        figure(fignum);
        clf;
        
%         wing_l = pathDB.wing_l(i);
%         jpL = pathDB.joint_pos_L(:,i)
%         jpR = pathDB.joint_pos_R(:,i);
%         
%         n_sphere = 2^5-1;
%         [x_sphere,y_sphere,z_sphere] = sphere(n_sphere);

%             hold on;
%             surf(-jpL(1)+5+(wing_l.*x_sphere),-jpL(2)+5+(wing_l.*y_sphere),-jpL(3)+5+(wing_l.*z_sphere),'FaceColor','black','EdgeColor','none');
%             alpha(0.2)
%             hold off;


%         for k = 1:length(x);
        
        a = Lwingbeat_loc(m,1):Lwingbeat_loc(start_L+m,2);
        b = Lwingbeat_loc(start_L+m,2):Lwingbeat_loc(m+1,1);
        
        Lwt_x_down = Lwingtip(a,1);
        Lwt_z_down = Lwingtip(a,3);
        
        Lwt_x_up = Lwingtip(b,1);
        Lwt_z_up = Lwingtip(b,3)
        
        pause
               
        xsi_L_down = real(atan2(2*(qL4(a).*qL1(a)+qL2(a).*qL3(a)),1-2*(qL1(a).^2+qL2(a).^2)));
        theta_L_down = real(asin(2*(qL4(a).*qL2(a)-qL3(a).*qL1(a))));
        phi_L_down = real(atan2(2*(qL4(a).*qL3(a)+qL1(a).*qL2(a)),1-2*(qL2(a).^2+qL3(a).^2)));
        
        xsi_L_up = real(atan2(2*(qL4(b).*qL1(b)+qL2(b).*qL3(b)),1-2*(qL1(b).^2+qL2(b).^2)));
        theta_L_up = real(asin(2*(qL4(b).*qL2(b)-qL3(b).*qL1(b))));
        phi_L_up = real(atan2(2*(qL4(b).*qL3(b)+qL1(b).*qL2(b)),1-2*(qL2(b).^2+qL3(b).^2)));
        
        
        for n = 1:length(a)
            
            T_phi_down = [1 0 0; 0 cos(phi_L_down(n)) sin(phi_L_down(n)); 0 -sin(phi_L_down(n)) cos(phi_L_down(n))];
            T_theta_down = [cos(theta_L_down(n)) 0 -sin(theta_L_down(n)); 0 1 0; sin(theta_L_down(n)) 0 cos(theta_L_down(n))];
            T_xsi_down = [cos(xsi_L_down(n)) -sin(xsi_L_down(n)) 0; sin(xsi_L_down(n)) cos(xsi_L_down(n)) 0; 0 0 1];
            
            
            
            quiv_x_down(n) = cos(theta_L_down(n));
            quiv_z_down(n) = -sin(theta_L_down(n));
            
            clear Vect_down
            
        end
            
        for n = 1:length(b)
            
            T_phi_up = [1 0 0; 0 cos(phi_L_up(n)) sin(phi_L_up(n)); 0 -sin(phi_L_up(n)) cos(phi_L_up(n))];
            T_theta_up = [cos(theta_L_up(n)) 0 -sin(theta_L_up(n)); 0 1 0; sin(theta_L_up(n)) 0 cos(theta_L_up(n))];
            T_xsi_up = [cos(xsi_L_up(n)) -sin(xsi_L_up(n)) 0; sin(xsi_L_up(n)) cos(xsi_L_up(n)) 0; 0 0 1];
            
                     
            quiv_x_up(n) = cos(theta_L_up(n));
            quiv_z_up(n) = -sin(theta_L_up(n));
            
            clear Vect_up
            
        end
        


        




        for k = 1:1
            subplot(2,2,1); surf(x{k},y{k},z{k},'facecolor',[0.6 0.8 1.0],'edgecolor','k','facelighting','phong'); %plot(-Lwt_x_down,-Lwt_z_down,'--r',-Lwt_x_up,-Lwt_z_up,'--r')
%             hold on
%             plot(-Lwt_x_down-0.1.*quiv_x_down', -Lwt_z_down-0.1.*quiv_z_down','.k')
%             plot(-Lwt_x_up-0.1.*quiv_x_up', -Lwt_z_up-0.1.*quiv_z_up','.k')
%             quiver(-Lwt_x_down,-Lwt_z_down,-quiv_x_down',-quiv_z_down',0.1,'k', 'ShowArrowHead','off')
%             quiver(-Lwt_x_up,-Lwt_z_up,-quiv_x_up',-quiv_z_up',0.1,'k', 'ShowArrowHead','off')
%             quiver(-Lwt_x_down,-Lwt_z_down,quiv_x_down',quiv_z_down',0.2,'k', 'ShowArrowHead','off')
%             quiver(-Lwt_x_up,-Lwt_z_up,quiv_x_up',quiv_z_up',0.2,'k','ShowArrowHead','off')
%             %plot(x_strkpln,z_strkpln,':')
%             quiver(0,0,-1,0,3,'--b')
%             quiver(0,0,0,-1,3,'--b')
%             axis equal
            %axis([minx maxx miny maxy minz maxz])
%             xlabel('x')
%             zlabel('z')
%             view(0,0)
            hold on;
            subplot(2,2,2); surf(x{k},y{k},z{k},'facecolor',[0.6 0.8 1.0],'edgecolor','k','facelighting','phong');
            axis([minx maxx miny maxy minz maxz])
            xlabel('x')
            ylabel('y')
            zlabel('z')
            view(90,radtodeg(pathDB.b_pitch(startframe-1+m,i))+90) %pathDB.b_pitch(startframe-1+m,i)
            hold on;
            subplot(2,2,3); surf(x{k},y{k},z{k},'facecolor',[0.6 0.8 1.0],'edgecolor','k','facelighting','phong');
            axis([minx maxx miny maxy minz maxz])
            xlabel('x')
            zlabel('z')
            view(180,0)
            hold on;
            subplot(2,2,4); surf(x{k},y{k},z{k},'facecolor',[0.6 0.8 1.0],'edgecolor','k','facelighting','phong');
            axis([minx maxx miny maxy minz maxz])
            xlabel('x')
            ylabel('y')
            zlabel('z')
            view(90,radtodeg(pathDB.b_pitch(startframe-1+m,i))-90)
            hold on;
        end
        
        clear Lwt_x_down Lwt_z_down Lwt_x_up Lwt_z_up quiv_x_down quiv_z_down quiv_x_up quiv_z_up

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


