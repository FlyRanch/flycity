
function Wing_kinematics_3D_moving(settings,pathDB,seq_nr,save_on_off)


    % Make a movie of the body orientation and the wingkinematics per
    % wingbeat. Can only create movies when save_on_off == 1:
    
    if save_on_off == 1
        
    close all;clc
    
    i = seq_nr;
    
        
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

    minx = -4;
    maxx = 4;
    miny = -4;
    maxy = 4;
    minz = -4;
    maxz = 4;
    
    
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
%         set(hfig, 'position', [100 100 1800 900])
        set(hfig, 'position', [100 100 900 900])
        clf;
        
%         phi = pathDB.phi_body_mean(m,1,i);
%         theta = pathDB.theta_body_mean(m,1,i);
%         xsi = pathDB.xsi_body_mean(m,1,i);

    % Get roll angle, angle of attack and angle of side-slip
    
%     k_end = find(isnan(pathDB.wingbeat_time(m,:,i))==0, 1 ,'last');
%     
%     k = pathDB.wingbeat_time(m,1:k_end,i);
    
    alfa = -pathDB.alfa_body_mean(m,1,i);
    
    beta = pi/2+pathDB.beta_body_mean(m,1,i);
    
    phi = pathDB.phi_body_mean(m,1,i);


%     qb = [pathDB.qb1_filt(m,i); pathDB.qb2_filt(m,i); pathDB.qb3_filt(m,i); pathDB.qb4_filt(m,i)];
% 
%     % R_test = quat2matNEW(q_test)
%     
%     % mean!!!
%     
%     find(isnan(pathDB.L_wingbeat_loc(:,1,i))==0, 1 ,'last')
%     
%     k = pathDB.wingbeat_time(m,1,i):pathDB.wingbeat_time(m,1,i)
% 
%     theta = asin(-2*qb(1)*qb(2)+2*qb(3)*qb(4));
% 
%     xsi = asin((2*qb(1)*qb(3)+2*qb(2)*qb(4))/cos(theta));
% 
%     phi = acos((-qb(1)^2+qb(2)^2-qb(3)^2+qb(4)^2)/cos(theta));

    % R_a = [ cos(theta_a)*cos(xsi_a) -sin(theta_a) cos(theta_a)*sin(xsi_a); ...
    %         sin(phi_a)*sin(xsi_a)+cos(phi_a)*cos(xsi_a)*sin(theta_a) cos(phi_a)*cos(theta_a) cos(phi_a)*sin(theta_a)*sin(xsi_a)-cos(xsi_a)*sin(phi_a); ...
    %         cos(xsi_a)*sin(phi_a)*sin(theta_a)-cos(phi_a)*sin(xsi_a) cos(theta_a)*sin(phi_a) cos(phi_a)*cos(xsi_a)+sin(phi_a)*sin(theta_a)*sin(xsi_a)]
    % 
    % phi_e = atan2(2*(q_test(4)*q_test(1)+q_test(2)*q_test(3)), 1-2*(q_test(1)^2+q_test(2)^2))
    % 
    % theta_e = asin(2*(q_test(4)*q_test(2)-q_test(3)*q_test(1)))
    % 
    % xsi_e = atan2(2*(q_test(4)*q_test(3)+q_test(1)*q_test(2)),1-2*(q_test(2)^2+q_test(3)^2))


            
%             subplot(2,3,1); Trajectory_movie_plot(settings,pathDB,i, m);
%             axis([-10 30 0 30 10 40])
%             xlabel('x')
%             ylabel('y')
%             zlabel('z')
%             hold on;
%             subplot(2,3,2); Sphere_plot_wing_kinematics(settings, pathDB, i, m, pi, theta, 0,1,0);
%             axis([minx maxx miny maxy minz maxz])
%             xlabel('x')
%             ylabel('y')
%             zlabel('z')
%             view(180,0)
%             hold on;
%             subplot(2,3,3); Sphere_plot_wing_kinematics(settings, pathDB, i, m, pi, 0, xsi,1,1);
%             axis([minx maxx miny maxy minz maxz])
%             xlabel('x')
%             ylabel('y')
%             zlabel('z')
%             view(0,90)
%             hold on;
%             subplot(2,3,4); Sphere_plot_wing_kinematics(settings, pathDB, i, m, phi, theta, xsi,1,1);
%             axis([minx maxx miny maxy minz maxz])
%             xlabel('x')
%             ylabel('y')
%             zlabel('z')
%             view(45,30)
%             hold on;
%             subplot(2,3,5); Sphere_plot_wing_kinematics(settings, pathDB, i, m, pi, theta, 0,0,1);
%             axis([minx maxx miny maxy minz maxz])
%             xlabel('x')
%             ylabel('y')
%             zlabel('z')
%             view(0,0)
%             hold on;
%             subplot(2,3,6); Sphere_plot_wing_kinematics(settings, pathDB, i, m, phi, 0, 0,1,1);
%             axis([minx maxx miny maxy minz maxz])
%             xlabel('x')
%             ylabel('y')
%             zlabel('z')
%             view(90,0)
%             hold on;


%             subplot(2,2,1); Sphere_plot_wing_kinematics(settings, pathDB, i, m, pi, theta, 0,1,0);
%             axis([minx maxx miny maxy minz maxz])
%             xlabel('x [mm]','FontSize',14)
%             set(gca, 'XTick', [-4,0,4],'FontSize',14);
%             ylabel('y [mm]','FontSize',14)
%             set(gca, 'YTick', [-4,0,4],'FontSize',14);
%             zlabel('z [mm]','FontSize',14)
%             set(gca, 'ZTick', [-4,0,4],'FontSize',14);
%             view(180,0)
%             title('Left view','FontSize',16)
%             hold on;
%             subplot(2,2,2); Sphere_plot_wing_kinematics(settings, pathDB, i, m, pi, theta, 0,0,1);
%             axis([minx maxx miny maxy minz maxz])
%             xlabel('x [mm]','FontSize',14)
%             set(gca, 'XTick', [-4,0,4],'FontSize',14);
%             ylabel('y [mm]','FontSize',14)
%             set(gca, 'YTick', [-4,0,4],'FontSize',14);
%             zlabel('z [mm]','FontSize',14)
%             set(gca, 'ZTick', [-4,0,4],'FontSize',14);
%             view(0,0)
%             title('Right view','FontSize',16)
%             hold on;
%             subplot(2,2,3); Sphere_plot_wing_kinematics(settings, pathDB, i, m, pi, 0, xsi,1,1);
%             axis([minx maxx miny maxy minz maxz])
%             xlabel('x [mm]','FontSize',14)
%             set(gca, 'XTick', [-4,0,4],'FontSize',14);
%             ylabel('y [mm]','FontSize',14)
%             set(gca, 'YTick', [-4,0,4],'FontSize',14);
%             zlabel('z [mm]','FontSize',14)
%             set(gca, 'ZTick', [-4,0,4],'FontSize',14);
%             view(0,90)
%             title('Top view','FontSize',16)
%             hold on;
%             subplot(2,2,4); Sphere_plot_wing_kinematics(settings, pathDB, i, m, phi, 0, 0,1,1);
%             axis([minx maxx miny maxy minz maxz])
%             xlabel('x [mm]','FontSize',14)
%             set(gca, 'XTick', [-4,0,4],'FontSize',14);
%             ylabel('y [mm]','FontSize',14)
%             set(gca, 'YTick', [-4,0,4],'FontSize',14);
%             zlabel('z [mm]','FontSize',14)
%             set(gca, 'ZTick', [-4,0,4],'FontSize',14);
%             view(90,0)
%             title('Front view','FontSize',16)
%             hold on;

    
            subplot(2,2,1); Sphere_plot_wing_kinematics(settings, pathDB, i, m, pi, alfa, 0,1,0);
            axis([minx maxx miny maxy minz maxz])
            xlabel('x [mm]','FontSize',14)
            set(gca, 'XTick', [-4,0,4],'FontSize',14);
            ylabel('y [mm]','FontSize',14)
            set(gca, 'YTick', [-4,0,4],'FontSize',14);
            zlabel('z [mm]','FontSize',14)
            set(gca, 'ZTick', [-4,0,4],'FontSize',14);
            view(180,0)
            title('Left view','FontSize',16)
            hold on;
            subplot(2,2,2); Sphere_plot_wing_kinematics(settings, pathDB, i, m, pi, alfa, 0,0,1);
            axis([minx maxx miny maxy minz maxz])
            xlabel('x [mm]','FontSize',14)
            set(gca, 'XTick', [-4,0,4],'FontSize',14);
            ylabel('y [mm]','FontSize',14)
            set(gca, 'YTick', [-4,0,4],'FontSize',14);
            zlabel('z [mm]','FontSize',14)
            set(gca, 'ZTick', [-4,0,4],'FontSize',14);
            view(0,0)
            title('Right view','FontSize',16)
            hold on;
            subplot(2,2,3); Sphere_plot_wing_kinematics(settings, pathDB, i, m, pi, 0, beta,1,1);
            axis([minx maxx miny maxy minz maxz])
            xlabel('x [mm]','FontSize',14)
            set(gca, 'XTick', [-4,0,4],'FontSize',14);
            ylabel('y [mm]','FontSize',14)
            set(gca, 'YTick', [-4,0,4],'FontSize',14);
            zlabel('z [mm]','FontSize',14)
            set(gca, 'ZTick', [-4,0,4],'FontSize',14);
            view(0,90)
            title('Top view','FontSize',16)
            hold on;
            subplot(2,2,4); Sphere_plot_wing_kinematics(settings, pathDB, i, m, phi, 0, 0,1,1);
            axis([minx maxx miny maxy minz maxz])
            xlabel('x [mm]','FontSize',14)
            set(gca, 'XTick', [-4,0,4],'FontSize',14);
            ylabel('y [mm]','FontSize',14)
            set(gca, 'YTick', [-4,0,4],'FontSize',14);
            zlabel('z [mm]','FontSize',14)
            set(gca, 'ZTick', [-4,0,4],'FontSize',14);
            view(90,0)
            title('Front view','FontSize',16)
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

        % Save plot


        saveas(fignum, [char(settings.plot_folders(3)) '/' char(settings.sequence_names(seq_nr)) '/Wing_kinematics_3D_moving/' int2str(m) ], 'fig')

        saveas(fignum, [char(settings.plot_folders(3)) '/' char(settings.sequence_names(seq_nr)) '/Wing_kinematics_3D_moving/' int2str(m) ], 'png')
        
    end
    
    cd(movdir)
    close(mov);
    cd(datadir)
    
    clear ManualFit
    
    
    end
    
    clear all
    


end