function Force_plot4(settings,pathDB,seq_nr,wingbeat_nr,phi,theta,xsi,left_on_off,right_on_off)


    % Generate trajectory plot with average forces and moment arms for the
    % upstroke and downstroke, left and right:
    
    clear i j start stop N M_down_L M_up_L M_down_R M_up_R F_joint_L_down F_joint_L_up F_joint_R_down F_joint_R_up a_L b_L a_R b_R
    
    i = seq_nr;
    
    j = wingbeat_nr;
    
    % start and stop point for the measurements
    start = find(isnan(pathDB.x(:,i))==0, 1 );
    stop = find(isnan(pathDB.x(:,i))==0, 1, 'last' );
    
    N = find(isnan(pathDB.wingbeat_time(:,1,i))==0,1,'last');
    
    
        % Get the variables for the left and right wings, upstroke and
        % downstroke:
        
        %Left wing---------------------------------------------------------
        
        
        M_down_L = find(isnan(pathDB.down_time_L(j,:,i))==0,1,'last');
        M_up_L = find(isnan(pathDB.up_time_L(j,:,i))==0,1,'last');
        
        
        a_L = pathDB.down_time_L(j,1:M_down_L,i);
        b_L = pathDB.up_time_L(j,1:M_up_L,i);
               
        
        % Forces at the joints:
        
        F_joint_L_down = pathDB.F_joint_L_down(start-1+a_L,:,i);
        F_joint_L_up   = pathDB.F_joint_L_up(start-1+b_L,:,i);
               
        
        % Calculate mean forces normal and tangential to the strokeplane:
        
        F_mean_L_down = mean(F_joint_L_down);
        
        F_mean_L_up = mean(F_joint_L_up);
        
        % Calculate the average arm on which the mean forces act:

        r_body_L_down = mean(pathDB.r_body_L_down(1:M_down_L,:,i).*norm(F_joint_L_down)/norm(F_mean_L_down))./(norm(F_joint_L_down)/norm(F_mean_L_down));
        r_body_L_up = mean(pathDB.r_body_L_up(1:M_up_L,:,i).*norm(F_joint_L_up)/norm(F_mean_L_up))./(norm(F_joint_L_up)/norm(F_mean_L_up));

%         r_body_L_down = mean(pathDB.r_body_L_down(1:M_down_L,:,i));
%         r_body_L_up = mean(pathDB.r_body_L_up(1:M_up_L,:,i));

        % Get wingtip path:
        
        Lwingtip_down = pathDB.Lwingtip(start-1+a_L,:,i);
        Lwingtip_up = pathDB.Lwingtip(start-1+b_L,:,i);

        % Get the joint position:
        
        joint_L = pathDB.joint_pos_L(:,i);
        
        %------------------------------------------------------------------
    
        
        %Right wing--------------------------------------------------------
        
        
        M_down_R = find(isnan(pathDB.down_time_R(j,:,i))==0,1,'last');
        M_up_R = find(isnan(pathDB.up_time_R(j,:,i))==0,1,'last');
        
        a_R = pathDB.down_time_R(j,1:M_down_R,i);
        b_R = pathDB.up_time_R(j,1:M_up_R,i);
        
        
        % Forces at the joints:
        
        F_joint_R_down = pathDB.F_joint_R_down(start-1+a_R,:,i);
        F_joint_R_up   = pathDB.F_joint_R_up(start-1+b_R,:,i);


        % Calculate mean forces normal and tangential to the strokeplane:
        
        F_mean_R_down = mean(F_joint_R_down);
        
        F_mean_R_up = mean(F_joint_R_up);
        
        % Calculate the average arm on which the mean forces act:
        
        r_body_R_down = mean(pathDB.r_body_R_down(1:M_down_R,:,i).*norm(F_joint_R_down)/norm(F_mean_R_down))./(norm(F_joint_R_down)/norm(F_mean_R_down));
        r_body_R_up = mean(pathDB.r_body_R_up(1:M_up_R,:,i).*norm(F_joint_R_up)/norm(F_mean_R_up))./(norm(F_joint_R_up)/norm(F_mean_R_up));
        
%         r_body_R_down = mean(pathDB.r_body_R_down(1:M_down_R,:,i));
%         r_body_R_up = mean(pathDB.r_body_R_up(1:M_up_R,:,i));

        % Get wingtip path:
        
        Rwingtip_down = pathDB.Rwingtip(start-1+a_R,:,i);
        Rwingtip_up = pathDB.Rwingtip(start-1+b_R,:,i);
        
        % Get the joint position:
        
        joint_R = pathDB.joint_pos_R(:,i);

        
        %------------------------------------------------------------------
        
        
        % Average postion and orientation of the body during the wingbeat:
        
        M_body = find(isnan(pathDB.wingbeat_time(j,:,i))==0,1,'last');
        
        t_wingbeat = pathDB.wingbeat_time(j,1:M_body,i);
        
        avg_pos = [mean(pathDB.x_filt(start-1+t_wingbeat,i)); mean(pathDB.y_filt(start-1+t_wingbeat,i)); mean(pathDB.z_filt(start-1+t_wingbeat,i))];
        
        avg_q_body = q_avg(pathDB.qb1_filt(start-1+t_wingbeat,i), pathDB.qb2_filt(start-1+t_wingbeat,i), pathDB.qb3_filt(start-1+t_wingbeat,i), pathDB.qb4_filt(start-1+t_wingbeat,i));
        
        
        %------------------------------------------------------------------
        
        % Get the force at the center of gravity:
        
        F_cg = pathDB.F_cg(start-1+t_wingbeat,:,i);
        
        % Calculate mean force over current wingbeat
        
        F_cg_mean = mean(F_cg);
        
        
        % Load body model:
        
        q_body = [sin(phi/2)*cos(theta/2)*cos(xsi/2)-cos(phi/2)*sin(theta/2)*sin(xsi/2); ...
                  cos(phi/2)*sin(theta/2)*cos(xsi/2)+sin(phi/2)*cos(theta/2)*sin(xsi/2); ...
                  cos(phi/2)*cos(theta/2)*sin(xsi/2)-sin(phi/2)*sin(theta/2)*cos(xsi/2); ...
                  cos(phi/2)*cos(theta/2)*cos(xsi/2)+sin(phi/2)*sin(theta/2)*sin(xsi/2)];
          
        q_body = q_body./norm(q_body);
        
        xh = [0 0 0 q_body(1) q_body(2) q_body(3) q_body(4) 0 0 0 1 0 0 0 1];
        
        [ x, y, z ] = load_body_model(settings, i, xh );
        
        x_path = pathDB.x_filt((start):(start-1+t_wingbeat(1)),i);
        
        y_path = pathDB.y_filt((start):(start-1+t_wingbeat(1)),i);
        
        z_path = pathDB.z_filt((start):(start-1+t_wingbeat(1)),i);
        
        %------------------------------------------------------------------
        
        % Transform all paramaters to global reference frame:
        
        DCM = quat2matNEW([q_body(1) q_body(2) q_body(3) q_body(4)]);
        
        r_L_down = DCM*r_body_L_down';
        r_L_up = DCM*r_body_L_up';
        
        r_R_down = DCM*r_body_R_down';
        r_R_up = DCM*r_body_R_up';
        
        F_L_down = DCM*F_mean_L_down';
        F_L_up = DCM*F_mean_L_up';
        
        F_R_down = DCM*F_mean_R_down';
        F_R_up = DCM*F_mean_R_up';
        
        j_L = DCM*joint_L;
        j_R = DCM*joint_R;

        F_cg_global = DCM*F_cg_mean';
        
        Lwt_down = DCM*Lwingtip_down';
        
        Lwt_up = DCM*Lwingtip_up';
        
        Rwt_down = DCM*Rwingtip_down';
        
        Rwt_up = DCM*Rwingtip_up';
        
        
        quiv_scale = 1.5e5;
        


        %------------------------------------------------------------------
        
        % Determine axis dimensions:
    
        % Calculate the region of fly motion


        buffer = 5;
        minx = -buffer;
        maxx = buffer;
        miny = -buffer;
        maxy = buffer;
        minz = -buffer;
        maxz = buffer;
        
        
        
        
        
        %------------------------------------------------------------------
        
        % Plot trajectory and arms and forces:
        
        hold on
        
        surf(x{1},y{1},z{1},'facecolor',[0.9 0.9 0.9],'edgecolor','k','facelighting','phong');
        
        hold on
        
        if left_on_off == 1
        
        plot3([j_L(1) j_L(1)+r_L_down(1)],[j_L(2) j_L(2)+r_L_down(2)],[j_L(3) j_L(3)+r_L_down(3)],'k')
        
        plot3([j_L(1) j_L(1)+r_L_up(1)],[j_L(2) j_L(2)+r_L_up(2)],[j_L(3) j_L(3)+r_L_up(3)],'k')
        
        quiver3(j_L(1)+r_L_down(1),j_L(2)+r_L_down(2),j_L(3)+r_L_down(3),quiv_scale*F_L_down(1),quiv_scale*F_L_down(2),quiv_scale*F_L_down(3),'r')
        
        quiver3(j_L(1)+r_L_up(1),j_L(2)+r_L_up(2),j_L(3)+r_L_up(3),quiv_scale*F_L_up(1),quiv_scale*F_L_up(2),quiv_scale*F_L_up(3),'b')
        
        plot3(j_L(1)+Lwt_down(1,:), j_L(2)+Lwt_down(2,:) , j_L(3)+Lwt_down(3,:),'r');
        
        plot3(j_L(1)+Lwt_up(1,:), j_L(2)+Lwt_up(2,:) , j_L(3)+Lwt_up(3,:),'r');
        
        plot3([j_L(1) j_L(1)+Lwt_down(1,1)], [j_L(2) j_L(2)+Lwt_down(2,1)], [j_L(3) j_L(3)+Lwt_down(3,1)], 'k')
        
        plot3([j_L(1) j_L(1)+Lwt_up(1,1)], [j_L(2) j_L(2)+Lwt_up(2,1)], [j_L(3) j_L(3)+Lwt_up(3,1)], 'k')
         
        end
        
        if right_on_off == 1
        
        plot3([j_R(1) j_R(1)+r_R_down(1)],[j_R(2) j_R(2)+r_R_down(2)],[j_R(3) j_R(3)+r_R_down(3)],'k')
        
        plot3([j_R(1) j_R(1)+r_R_up(1)],[j_R(2) j_R(2)+r_R_up(2)],[j_R(3) j_R(3)+r_R_up(3)],'k')
        
        quiver3(j_R(1)+r_R_down(1),j_R(2)+r_R_down(2),j_R(3)+r_R_down(3),quiv_scale*F_R_down(1),quiv_scale*F_R_down(2),quiv_scale*F_R_down(3),'r')
        
        quiver3(j_R(1)+r_R_up(1),j_R(2)+r_R_up(2),j_R(3)+r_R_up(3),quiv_scale*F_R_up(1),quiv_scale*F_R_up(2),quiv_scale*F_R_up(3),'b')
        
        plot3(j_R(1)+Rwt_down(1,:), j_R(2)+Rwt_down(2,:) , j_R(3)+Rwt_down(3,:),'g');
        
        plot3(j_R(1)+Rwt_up(1,:), j_R(2)+Rwt_up(2,:) , j_R(3)+Rwt_up(3,:),'g');
        
        plot3([j_R(1) j_R(1)+Rwt_down(1,1)], [j_R(2) j_R(2)+Rwt_down(2,1)], [j_R(3) j_R(3)+Rwt_down(3,1)], 'k')
        
        plot3([j_R(1) j_R(1)+Rwt_up(1,1)], [j_R(2) j_R(2)+Rwt_up(2,1)], [j_R(3) j_R(3)+Rwt_up(3,1)], 'k')
        
        end
                
        quiver3(0,0,0,quiv_scale*F_cg_global(1),quiv_scale*F_cg_global(2),quiv_scale*F_cg_global(3),'Color',[0 0.5 0])
        
        axis equal
        axis([minx maxx miny maxy minz maxz])
        hold off
        
        hold off
        
        
        
        
end


