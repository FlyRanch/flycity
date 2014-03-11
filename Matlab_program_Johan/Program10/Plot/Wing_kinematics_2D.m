function Wing_kinematics_2D(settings,pathDB,seq_nr,save_on_off)


    % Plot the wing kinematics in a 2D plane and plot the forces normal and
    % tangential to the stroke plane in this 2D plane:
    
    
    if save_on_off == 1
        
    close all;clc
    
    i = seq_nr;

        
% %     cd(char(settings.path_names(1)))
% %     
% %     datadir = cd;
% % 
% %     cd([char(settings.plot_folders(4)) '/' char(settings.sequence_names(i))]);
% %     
% %     movdir = cd;
% % 
% %     movfilename = ['2D_wing_kinematics_movie_' char(settings.sequence_names(i)) '.mj2'];
% %     mov = VideoWriter(movfilename, 'Motion JPEG AVI');
% %     mov.FrameRate = 5;
% %     mov.Quality = 100;
% % 
% %     
% %     cd(datadir)
% % 
% %     fignum = i;
% % 
% %     %changes the scaling of the final figure, as percentage of original image
% %     %resolution to make fit on screen.
% %     plotscale = .8;
% %     
% %     samplerate = 1;
% %     
% %     if pathDB.L_wingbeat_loc(1,1,i) > pathDB.L_wingbeat_loc(1,2,i)
% %         end_L = find(isnan(pathDB.L_wingbeat_loc(:,1,i))==0, 1 ,'last')-1;
% %     else
% %         end_L = find(isnan(pathDB.L_wingbeat_loc(:,1,i))==0, 1 ,'last')-2;
% %     end
% %     
% %     startframe = 1;
% %     endframe = end_L;
% % 
% %     frames = [startframe endframe]; % number of frames in movie
% % 
% %     movidx = frames(1):samplerate:frames(2);
% %     numframes = length(movidx);
    
    
    %%
    
    % Calculate the normal and tangential strokeplane forces:
    
    % start and stop point for the measurements
    start = find(isnan(pathDB.x(:,i))==0, 1 );
    stop = find(isnan(pathDB.x(:,i))==0, 1, 'last' );
    
    
    Fr_L_down_b = zeros(length(pathDB.down_time_L(:,1,i)),2);
    Fr_L_up_b = zeros(length(pathDB.up_time_L(:,1,i)),2);
    Fr_R_down_b = zeros(length(pathDB.down_time_R(:,1,i)),2);
    Fr_R_up_b = zeros(length(pathDB.up_time_R(:,1,i)),2);
    
    Fx__down_st = zeros(length(pathDB.down_time_L(:,1,i)),1);
    Fx_L_up_st = zeros(length(pathDB.up_time_L(:,1,i)),1);
    Fx_R_down_st = zeros(length(pathDB.down_time_R(:,1,i)),1);
    Fx_R_up_st = zeros(length(pathDB.up_time_R(:,1,i)),1);
    
    Fy_L_down_st = zeros(length(pathDB.down_time_L(:,1,i)),1);
    Fy_L_up_st = zeros(length(pathDB.up_time_L(:,1,i)),1);
    Fy_R_down_st = zeros(length(pathDB.down_time_R(:,1,i)),1);
    Fy_R_up_st = zeros(length(pathDB.up_time_R(:,1,i)),1);
    
    Fz_L_down_st = zeros(length(pathDB.down_time_L(:,1,i)),1);
    Fz_L_up_st = zeros(length(pathDB.up_time_L(:,1,i)),1);
    Fz_R_down_st = zeros(length(pathDB.down_time_R(:,1,i)),1);
    Fz_R_up_st = zeros(length(pathDB.up_time_R(:,1,i)),1);    
    
    % Transform the forces first to the body frame and than to the
    % wing-stroke plane
    
    beta = -(55/180)*pi;
    
    Rot_mat_roll = [1 0 0; ...
                    0 -1 0; ...
                    0 0 -1];
    
    Rot_mat_pitch = [cos(beta) 0 -sin(beta); ...
                     0 1 0; ...
                     sin(beta) 0 cos(beta)];

    Rot_mat = Rot_mat_pitch;
    
    % Left wing:           
               
    for j = 1:find(isnan(pathDB.down_time_L(:,1,i))==0,1,'last')
        
        a_L = pathDB.down_time_L(j,1:find(isnan(pathDB.down_time_L(j,:,i))==0,1,'last'));
        b_L = pathDB.up_time_L(j,1:find(isnan(pathDB.up_time_L(j,:,i))==0,1,'last'));

        
        for k = 1:length(a_L)
            
        DCM = quat2matNEW([pathDB.qL1_filt2(start-1+a_L(k),i) pathDB.qL2_filt2(start-1+a_L(k),i) pathDB.qL3_filt2(start-1+a_L(k),i) pathDB.qL4_filt2(start-1+a_L(k),i)]);
        
        tempF1 = Rot_mat*DCM*[-pathDB.Ft_L_down(j,k,i); 0; -pathDB.Fn_L_down(j,k,i)]; 
        
        tempF2 = Rot_mat*pathDB.F_joint_L_down(start-1+a_L(k),:,i)';
        
        Fr_L_down_b(j,:) = 2e7.*[sign(tempF1(1))*sqrt(tempF1(1)^2+tempF2(2)^2) tempF1(3)];
        
        Fx__down_st(j) = tempF2(1);
        
        Fy_L_down_st(j) = tempF2(2);
        
        Fz_L_down_st(j) = tempF2(3);
            
        clear DCM tempF1 tempF2
        
        end
        
        for k = 1:length(b_L)
        
        DCM = quat2matNEW([pathDB.qL1_filt2(start-1+b_L(k),i) pathDB.qL2_filt2(start-1+b_L(k),i) pathDB.qL3_filt2(start-1+b_L(k),i) pathDB.qL4_filt2(start-1+b_L(k),i)]);
        
        tempF1 = Rot_mat*DCM*[-pathDB.Ft_L_up(j,k,i); 0; -pathDB.Fn_L_up(j,k,i)]; 
        
        tempF2 = Rot_mat*pathDB.F_joint_L_up(start-1+b_L(k),:,i)';
        
        Fr_L_up_b(j,:) = 2e7.*[sign(tempF1(1))*sqrt(tempF1(1)^2+tempF2(2)^2) tempF1(3)];
        
        Fx_L_up_st(j) = tempF2(1);
        
        Fy_L_up_st(j) = tempF2(2);
        
        Fz_L_up_st(j) = tempF2(3);
            
        clear DCM tempF1 tempF2
        
        end

    end
    
    % Right wing:           
               
    for j = 1:find(isnan(pathDB.down_time_R(:,1,i))==0,1,'last')
        
        a_R = pathDB.down_time_R(j,1:find(isnan(pathDB.down_time_R(j,:,i))==0,1,'last'));
        b_R = pathDB.up_time_R(j,1:find(isnan(pathDB.up_time_R(j,:,i))==0,1,'last'));

        
        for k = 1:length(a_R)
        
        DCM = quat2matNEW([pathDB.qR1_filt2(start-1+a_R(k),i) pathDB.qR2_filt2(start-1+a_R(k),i) pathDB.qR3_filt2(start-1+a_R(k),i) pathDB.qR4_filt2(start-1+a_R(k),i)]);
        
        tempF1 = Rot_mat*DCM*[-pathDB.Ft_R_down(j,k,i); 0; -pathDB.Fn_R_down(j,k,i)]; 
        
        tempF2 = Rot_mat*pathDB.F_joint_R_down(start-1+a_R(k),:,i)';
        
        Fr_R_down_b(j,:) = 2e7.*[sign(tempF1(1))*sqrt(tempF1(1)^2+tempF2(2)^2) tempF1(3)];
        
        Fx_R_down_st(j) = tempF2(1);
        
        Fy_R_down_st(j) = tempF2(2);
        
        Fz_R_down_st(j) = tempF2(3);
            
        clear DCM tempF1 tempF2
        
        end
        
        for k = 1:length(b_R)
        
        DCM = quat2matNEW([pathDB.qR1_filt2(start-1+b_R(k),i) pathDB.qR2_filt2(start-1+b_R(k),i) pathDB.qR3_filt2(start-1+b_R(k),i) pathDB.qR4_filt2(start-1+b_R(k),i)]);
        
        tempF1 = Rot_mat*DCM*[-pathDB.Ft_R_up(j,k,i); 0; -pathDB.Fn_R_up(j,k,i)]; 
        
        tempF2 = Rot_mat*pathDB.F_joint_R_up(start-1+b_R(k),:,i)';
        
        Fr_R_up_b(j,:) = 2e7.*[sign(tempF1(1))*sqrt(tempF1(1)^2+tempF2(2)^2) tempF1(3)];
        
        Fx_R_up_st(j) = tempF2(1);
        
        Fy_R_up_st(j) = tempF2(2);
        
        Fz_R_up_st(j) = tempF2(3);
            
        clear DCM tempF1 tempF2
        
        end
        

    end
    

    
    

    
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
  

    
%     
%     open(mov);
%     
%     set(gca,'nextplot','replacechildren');
%     set(gcf,'Renderer','zbuffer');

    
    for m= 1:find(isnan(pathDB.down_time_L(:,1,i))==0,1,'last')
        
%         kk = frames(1)+(m-1);
% 
% 
%         hfig = figure(fignum);
%         set(hfig, 'position', [100 100 1800 900])
%         clf;

            % Get strokeplane angles:

            %Left wing


            a_L = pathDB.down_time_L(m,1:find(isnan(pathDB.down_time_L(m,:,i))==0,1,'last'));
            b_L = pathDB.up_time_L(m,1:find(isnan(pathDB.up_time_L(m,:,i))==0,1,'last'));

            phi_L_down = pathDB.phi_L(start-1+a_L,i);
            theta_L_down = pathDB.theta_L(start-1+a_L,i);
            eta_L_down = pathDB.eta_L(start-1+a_L,i);

            phi_L_up = pathDB.phi_L(start-1+b_L,i);
            theta_L_up = pathDB.theta_L(start-1+b_L,i);
            eta_L_up = pathDB.eta_L(start-1+b_L,i);

            
            %Right wing


             a_R = pathDB.down_time_R(m,1:find(isnan(pathDB.down_time_R(m,:,i))==0,1,'last'));
             b_R = pathDB.up_time_R(m,1:find(isnan(pathDB.up_time_R(m,:,i))==0,1,'last'));

             phi_R_down = pathDB.phi_R(start-1+a_R,i);
             theta_R_down = pathDB.theta_R(start-1+a_R,i);
             eta_R_down = pathDB.eta_R(start-1+a_R,i);

             phi_R_up = pathDB.phi_R(start-1+b_R,i);
             theta_R_up = pathDB.theta_R(start-1+b_R,i);
             eta_R_up = pathDB.eta_R(start-1+b_R,i);
%              
%              scale_c = 5;
%              
%              figure()
%              
%              subplot(2,2,1); plot(radtodeg(phi_L_down),radtodeg(theta_L_down))
%              hold on
%              for k = 1:length(a_L)
%                  plot([radtodeg(phi_L_down(k))-0.75*scale_c*cos(eta_L_down(k)) radtodeg(phi_L_down(k))+0.25*scale_c*cos(eta_L_down(k))],[radtodeg(theta_L_down(k))-0.75*scale_c*sin(eta_L_down(k)) radtodeg(theta_L_down(k))+0.25*scale_c*sin(eta_L_down(k))],'k')
%                  quiver(radtodeg(phi_L_down(k)),radtodeg(theta_L_down(k)),Fr_L_down_b(k,1),Fr_L_down_b(k,2),'r')
%              end
%              axis([-70 100 -45 45])
%              axis equal
%              hold on
%              subplot(2,2,2); plot(radtodeg(phi_L_up),radtodeg(theta_L_up))
%              hold on
%              for k = 1:length(b_L)
%                  plot([radtodeg(phi_L_up(k))-0.75*scale_c*cos(eta_L_up(k)) radtodeg(phi_L_up(k))+0.25*scale_c*cos(eta_L_up(k))],[radtodeg(theta_L_up(k))-0.75*scale_c*sin(eta_L_up(k)) radtodeg(theta_L_up(k))+0.25*scale_c*sin(eta_L_up(k))],'k')
%                  quiver(radtodeg(phi_L_up(k)),radtodeg(theta_L_up(k)),Fr_L_up_b(k,1),Fr_L_up_b(k,2),'r')
%              end
%              axis([-70 100 -45 45])
%              axis equal
%              hold on
%              subplot(2,2,3); plot(radtodeg(phi_R_down),radtodeg(theta_R_down))
%              hold on
%              for k = 1:length(a_R)
%                  plot([radtodeg(phi_R_down(k))-0.75*scale_c*cos(eta_R_down(k)) radtodeg(phi_R_down(k))+0.25*scale_c*cos(eta_R_down(k))],[radtodeg(theta_R_down(k))-0.75*scale_c*sin(eta_R_down(k)) radtodeg(theta_R_down(k))+0.25*scale_c*sin(eta_R_down(k))],'k')
%                  quiver(radtodeg(phi_R_down(k)),radtodeg(theta_R_down(k)),Fr_R_down_b(k,1),Fr_R_down_b(k,2),'r')
%              end
%              axis([-70 100 -45 45])
%              axis equal
%              hold on
%              subplot(2,2,4); plot(radtodeg(phi_R_up),radtodeg(theta_R_up))
%              hold on
%              for k = 1:length(b_R)
%                  plot([radtodeg(phi_R_up(k))-0.75*scale_c*cos(eta_R_up(k)) radtodeg(phi_R_up(k))+0.25*scale_c*cos(eta_R_up(k))],[radtodeg(theta_R_up(k))-0.75*scale_c*sin(eta_R_up(k)) radtodeg(theta_R_up(k))+0.25*scale_c*sin(eta_R_up(k))],'k')
%                  quiver(radtodeg(phi_R_up(k)),radtodeg(theta_R_up(k)),Fr_R_up_b(k,1),Fr_R_up_b(k,2),'r')
%              end
%              axis([-70 100 -45 45])
%              axis equal
%              hold off
%                 
%              clear a_L b_L a_R b_R phi_L_down theta_L_down eta_L_down phi_L_up theta_L_up eta_L_up phi_R_down theta_R_down eta_R_down phi_R_up theta_R_up eta_R_up
             
figure()
plot(-Fz_R_down_st)

pause

     

%         %set(gca,'xlim',[minx maxx],'ylim',[miny maxy],'zlim',[minz maxz]);
% 
% 
%         figure(fignum);
%         set(fignum,'color','w');
%         
%         frame = getframe(fignum);
%         
%         cd(movdir)
%         writeVideo(mov,frame);
%         writeVideo(mov,frame);
%         writeVideo(mov,frame);
%         writeVideo(mov,frame);
%         writeVideo(mov,frame);
%         cd(datadir)
% 
%         % Save plot
% 
% 
%         saveas(fignum, [char(settings.plot_folders(3)) '/' char(settings.sequence_names(seq_nr)) '/Wing_kinematics_3D_moving/' int2str(m) ], 'fig')
% 
%         saveas(fignum, [char(settings.plot_folders(3)) '/' char(settings.sequence_names(seq_nr)) '/Wing_kinematics_3D_moving/' int2str(m) ], 'png')
        
    end
    
%     cd(movdir)
%     close(mov);
%     cd(datadir)
%     
%     clear ManualFit
    
    
    end
%     
%     clear all


end

