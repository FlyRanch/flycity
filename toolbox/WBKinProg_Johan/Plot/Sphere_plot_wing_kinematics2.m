function Sphere_plot_wing_kinematics2(settings, pathDB, seq_nr, wingbeat_nr, phi, theta, xsi, left_on_off, right_on_off)


    % Program that generates 3D plots of the wingtip_path, postion of the
    % fruit fly body and the orientation of the wing for a single given
    % wingbeat
    
    start_meas = find(isnan(pathDB.qb1(:,seq_nr))==0, 1 );

    
    a = (start_meas-1+pathDB.L_wingbeat_loc(wingbeat_nr,2,seq_nr)):1:(start_meas-1+pathDB.L_wingbeat_loc(wingbeat_nr+1,2,seq_nr)-1);
    b = (start_meas-1+pathDB.R_wingbeat_loc(wingbeat_nr,2,seq_nr)):1:(start_meas-1+pathDB.R_wingbeat_loc(wingbeat_nr+1,2,seq_nr)-1);    
    
    % Calculate orientation and position cross-section left wing:
    
    LE_L = pathDB.wing_l(seq_nr).*[0.03; -1; 0];
    TE_L = pathDB.wing_l(seq_nr).*[-0.05; -1; 0];
    
    LE_R = pathDB.wing_l(seq_nr).*[0.03; 1; 0];
    TE_R = pathDB.wing_l(seq_nr).*[-0.05; 1; 0];
    
    L_cross_LE = zeros(3,length(a));
    L_cross_TE = zeros(3,length(a));
    
    R_cross_LE = zeros(3,length(b));
    R_cross_TE = zeros(3,length(b));
    
    if left_on_off == 1
    
    for i = 1:length(a)
        
    DCM_L = quat2matNEW([pathDB.qL1_filt2(a(i),seq_nr) pathDB.qL2_filt2(a(i),seq_nr) pathDB.qL3_filt2(a(i),seq_nr) pathDB.qL4_filt2(a(i),seq_nr)]);

    L_cross_LE(:,i) = DCM_L*LE_L;
    
    L_cross_TE(:,i) = DCM_L*TE_L;
    
    end
    
    end
    
    if right_on_off == 1
    
    for i = 1:length(b)
        
    DCM_R = quat2matNEW([pathDB.qR1_filt2(b(i),seq_nr) pathDB.qR2_filt2(b(i),seq_nr) pathDB.qR3_filt2(b(i),seq_nr) pathDB.qR4_filt2(b(i),seq_nr)]);

    R_cross_LE(:,i) = DCM_R*LE_R;
    
    R_cross_TE(:,i) = DCM_R*TE_R;
    
    end
    
    end

    
    % Retrieve the location of the left and right wing joints:
    
    jLocL = pathDB.joint_pos_L(:,seq_nr);
    
    jLocR = pathDB.joint_pos_R(:,seq_nr);
    
    % Retrieve wingtip path:
    
    Lwingtip = pathDB.Lwingtip(a,:,seq_nr);
    Rwingtip = pathDB.Rwingtip(b,:,seq_nr);
    
    % Calculate identity begin upstroke:
    
    if left_on_off == 1
    
    if pathDB.L_wingbeat_loc(wingbeat_nr,1,seq_nr) > pathDB.L_wingbeat_loc(wingbeat_nr,2,seq_nr)
    
        i_upst_L = pathDB.L_wingbeat_loc(wingbeat_nr,1,seq_nr)-pathDB.L_wingbeat_loc(wingbeat_nr,2,seq_nr)+1;
        
    else
        
        i_upst_L = pathDB.L_wingbeat_loc(wingbeat_nr+1,1,seq_nr)-pathDB.L_wingbeat_loc(wingbeat_nr,2,seq_nr)+1;
        
    end
    
    end
    
    if right_on_off == 1
    
    if pathDB.R_wingbeat_loc(wingbeat_nr,1,seq_nr) > pathDB.R_wingbeat_loc(wingbeat_nr,2,seq_nr)
    
        i_upst_R = pathDB.R_wingbeat_loc(wingbeat_nr,1,seq_nr)-pathDB.R_wingbeat_loc(wingbeat_nr,2,seq_nr)+1;
        
    else
        
        i_upst_R = pathDB.R_wingbeat_loc(wingbeat_nr+1,1,seq_nr)-pathDB.R_wingbeat_loc(wingbeat_nr,2,seq_nr)+1;
        
    end
    
    end
    
    % Retrieve velocity vector:
    
    u = pathDB.u_body_mean(wingbeat_nr,1,seq_nr); %/pathDB.U_body_mean(wingbeat_nr,1,seq_nr);
    v = pathDB.v_body_mean(wingbeat_nr,1,seq_nr); %/pathDB.U_body_mean(wingbeat_nr,1,seq_nr);
    w = pathDB.w_body_mean(wingbeat_nr,1,seq_nr); %/pathDB.U_body_mean(wingbeat_nr,1,seq_nr);
    
   
    % Rotate body and wingtip path accoring to given body roll, pitch and
    % yaw:

  
    R_phi   = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
    R_theta = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    R_xsi   = [cos(xsi) -sin(xsi) 0; sin(xsi) cos(xsi) 0; 0 0 1];
    
    R_matrix = R_phi*R_theta*R_xsi;
        
    % Create body:
    
    % Calculate body quaternion from Euler angles:
    
    q_body = [sin(phi/2)*cos(theta/2)*cos(xsi/2)-cos(phi/2)*sin(theta/2)*sin(xsi/2); ...
              cos(phi/2)*sin(theta/2)*cos(xsi/2)+sin(phi/2)*cos(theta/2)*sin(xsi/2); ...
              cos(phi/2)*cos(theta/2)*sin(xsi/2)-sin(phi/2)*sin(theta/2)*cos(xsi/2); ...
              cos(phi/2)*cos(theta/2)*cos(xsi/2)+sin(phi/2)*sin(theta/2)*sin(xsi/2)];
          
    q_body = q_body./norm(q_body);
    
    SOLN = [0 0 0 q_body(1) q_body(2) q_body(3) q_body(4) 0 0 0 1 0 0 0 1];
    
    % Load body model:
    
    [ x, y, z ] = load_body_model(settings, seq_nr, SOLN );
    
    DCM_body = quat2matNEW(q_body);
    
    R_matrix = DCM_body;
    
    
    jLocL_rot = R_matrix*jLocL;
    
    if left_on_off == 1
    
    for i = 1:length(a)
    
    Lwingtip_rot(:,i)   = R_matrix*[jLocL(1)+Lwingtip(i,1); jLocL(2)+Lwingtip(i,2); jLocL(3)+Lwingtip(i,3)];
    L_cross_LE_rot(:,i) = R_matrix*[jLocL(1)+L_cross_LE(1,i); jLocL(2)+L_cross_LE(2,i); jLocL(3)+L_cross_LE(3,i)];
    L_cross_TE_rot(:,i) = R_matrix*[jLocL(1)+L_cross_TE(1,i); jLocL(2)+L_cross_TE(2,i); jLocL(3)+L_cross_TE(3,i)];
    
    end
    
    
    end
    
    if right_on_off == 1
    
    jLocR_rot = R_matrix*jLocR;
    
    for i = 1:length(b)
        
    Rwingtip_rot(:,i)   = R_matrix*[jLocR(1)+Rwingtip(i,1); jLocR(2)+Rwingtip(i,2); jLocR(3)+Rwingtip(i,3)];
    R_cross_LE_rot(:,i) = R_matrix*[jLocR(1)+R_cross_LE(1,i); jLocR(2)+R_cross_LE(2,i); jLocR(3)+R_cross_LE(3,i)];
    R_cross_TE_rot(:,i) = R_matrix*[jLocR(1)+R_cross_TE(1,i); jLocR(2)+R_cross_TE(2,i); jLocR(3)+R_cross_TE(3,i)];   
    
    end
    
    end
    
    quiver_pos = R_matrix*[3+0.005*u; 0.005*v;0.005*w];
    
    quiver_vel = R_matrix*[-u; -v; -w];
    
    % Plot wingtip path and the wing orientation in one 3D plot

    hold on
    surf(x{1},y{1},z{1},'facecolor','b','edgecolor','k','facelighting','phong');
    
    if left_on_off == 1
    plot3(Lwingtip_rot(1,:),Lwingtip_rot(2,:),Lwingtip_rot(3,:),'r');
    plot3([jLocL_rot(1) Lwingtip_rot(1,1)], [jLocL_rot(2) Lwingtip_rot(2,1)], [jLocL_rot(3) Lwingtip_rot(3,1)],'k');
    plot3([jLocL_rot(1) Lwingtip_rot(1,i_upst_L)], [jLocL_rot(2) Lwingtip_rot(2,i_upst_L)], [jLocL_rot(3) Lwingtip_rot(3,i_upst_L)],'k');
    end
    
    if right_on_off == 1
    plot3(Rwingtip_rot(1,:),Rwingtip_rot(2,:),Rwingtip_rot(3,:),'g');
    plot3([jLocR_rot(1) Rwingtip_rot(1,1)], [jLocR_rot(2) Rwingtip_rot(2,1)], [jLocR_rot(3) Rwingtip_rot(3,1)],'k');
    plot3([jLocR_rot(1) Rwingtip_rot(1,i_upst_R)], [jLocR_rot(2) Rwingtip_rot(2,i_upst_R)], [jLocR_rot(3) Rwingtip_rot(3,i_upst_R)],'k');
    end
    
    quiver3(quiver_pos(1),quiver_pos(2),quiver_pos(3),quiver_vel(1),quiver_vel(2),quiver_vel(3),0.005,'k');
    
    if left_on_off == 1
        
    for k = 1:length(a)
        %plot3([jLocL_rot(1) Lwingtip_rot(1,k)], [jLocL_rot(2) Lwingtip_rot(2,k)], [jLocL_rot(3) Lwingtip_rot(3,k)],'k');
        plot3([L_cross_LE_rot(1,k) L_cross_TE_rot(1,k)], [L_cross_LE_rot(2,k) L_cross_TE_rot(2,k)], [L_cross_LE_rot(3,k) L_cross_TE_rot(3,k)],'k');
    end
    
    end
    
    if right_on_off == 1
        
    for k = 1:length(b)
        %plot3([jLocR_rot(1) Rwingtip_rot(1,k)], [jLocR_rot(2) Rwingtip_rot(2,k)], [jLocR_rot(3) Rwingtip_rot(3,k)],'k');
        plot3([R_cross_LE_rot(1,k) R_cross_TE_rot(1,k)], [R_cross_LE_rot(2,k) R_cross_TE_rot(2,k)], [R_cross_LE_rot(3,k) R_cross_TE_rot(3,k)],'k');
    end
    
    end
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold off



end

