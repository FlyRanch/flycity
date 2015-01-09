function strokeplane_forces(settings,pathDB,seq_nr,wingbeat_nr)

    clear i j start stop N M_down_L M_up_L M_down_R M_up_R F_joint_L_down F_joint_L_up F_joint_R_down F_joint_R_up a_L b_L a_R b_R


    % Plot 2D wing kinematics
    
    i = seq_nr;
    
    j = wingbeat_nr;
    
    % start and stop point for the measurements
    start = find(isnan(pathDB.x(:,i))==0, 1 );
    stop = find(isnan(pathDB.x(:,i))==0, 1, 'last' );
    
    N = find(isnan(pathDB.wingbeat_time(:,1,i))==0,1,'last');
    
    % Strokeplane rotation matrix:
    
    beta = -(55/180)*pi;
    
    Rot_mat_roll = [1 0 0; ...
                    0 -1 0; ...
                    0 0 -1];
    
    Rot_mat_pitch = [cos(beta) 0 -sin(beta); ...
                     0 1 0; ...
                     sin(beta) 0 cos(beta)];

    Rot_mat = Rot_mat_pitch;

    
    % Wing_length
    
    wing_l = pathDB.wing_l(i);
        
              
        % Get the variables for the left and right wings, upstroke and
        % downstroke:
        
        %Left wing---------------------------------------------------------
        
        
        M_down_L = find(isnan(pathDB.down_time_L(j,:,i))==0,1,'last');
        M_up_L = find(isnan(pathDB.up_time_L(j,:,i))==0,1,'last');
        
        
        a_L = pathDB.down_time_L(j,1:M_down_L,i);
        b_L = pathDB.up_time_L(j,1:M_up_L,i);
        
        % Strokeplane angles:
        
        phi_L_down = pathDB.phi_L(start-1+a_L,i);
        theta_L_down = pathDB.theta_L(start-1+a_L,i);
        eta_L_down = pathDB.eta_L(start-1+a_L,i);

        phi_L_up = pathDB.phi_L(start-1+b_L,i);
        theta_L_up = pathDB.theta_L(start-1+b_L,i);
        eta_L_up = pathDB.eta_L(start-1+b_L,i);
        
        
        % Forces at the joints:
        
        F_joint_L_down = pathDB.F_joint_L_down(start-1+a_L,:,i);
        F_joint_L_up   = pathDB.F_joint_L_up(start-1+b_L,:,i);
        
        
        % Transfer forces to the strokeplane:
        
        F_st_L_down(:,1:M_down_L) = Rot_mat*F_joint_L_down(1:M_down_L,:)';
        F_st_L_up(:,1:M_up_L) = Rot_mat*F_joint_L_up(1:M_up_L,:)';
        
        
        % Calculate mean forces normal and tangential to the strokeplane:
        
        F_mean_L_down = mean(F_st_L_down,2);
        
        F_mean_L_up = mean(F_st_L_up,2);
        
        % Calculate the average arm on which the mean forces act:
        
        r_body_L_down = mean(pathDB.r_body_L_down(1:M_down_L,:,i));
        r_body_L_up = mean(pathDB.r_body_L_up(1:M_up_L,:,i));
        
        % convert to strokeplane reference frame:
        
        r_st_L_down = Rot_mat*r_body_L_down';
        r_st_L_up = Rot_mat*r_body_L_up';
        
        % conver to strokeplane angles phi and theta:
        
        r_phi_L_down = real(atan2(-r_st_L_down(1),-r_st_L_down(2)));
        r_theta_L_down = real(atan2(-r_st_L_down(3),sqrt(r_st_L_down(1)^2+r_st_L_down(2)^2)));
        
        r_phi_L_up = real(atan2(-r_st_L_up(1),-r_st_L_up(2)));
        r_theta_L_up = real(atan2(-r_st_L_up(3),sqrt(r_st_L_up(1)^2+r_st_L_up(2)^2)));
        
        %------------------------------------------------------------------
    
        
        %Right wing--------------------------------------------------------
        
        
        M_down_R = find(isnan(pathDB.down_time_R(j,:,i))==0,1,'last');
        M_up_R = find(isnan(pathDB.up_time_R(j,:,i))==0,1,'last');
        
        a_R = pathDB.down_time_R(j,1:M_down_R,i);
        b_R = pathDB.up_time_R(j,1:M_up_R,i);
        
        % Strokeplane angles:
        
        phi_R_down = pathDB.phi_R(start-1+a_R,i);
        theta_R_down = pathDB.theta_R(start-1+a_R,i);
        eta_R_down = pathDB.eta_R(start-1+a_R,i);
        
        phi_R_up = pathDB.phi_R(start-1+b_R,i);
        theta_R_up = pathDB.theta_R(start-1+b_R,i);
        eta_R_up = pathDB.eta_R(start-1+b_R,i);
        
        % Forces at the joints:
        
        F_joint_R_down = pathDB.F_joint_R_down(start-1+a_R,:,i);
        F_joint_R_up   = pathDB.F_joint_R_up(start-1+b_R,:,i);

        
        
        % Transfer forces to the strokeplane:
        
        F_st_R_down(:,1:M_down_R) = Rot_mat*F_joint_R_down(1:M_down_R,:)';
        F_st_R_up(:,1:M_up_R) = Rot_mat*F_joint_R_up(1:M_up_R,:)';
       

        % Calculate mean forces normal and tangential to the strokeplane:
        
        F_mean_R_down = mean(F_st_R_down,2);
        
        F_mean_R_up = mean(F_st_R_up,2);
        
        % Calculate the average arm on which the mean forces act:
        
        r_body_R_down = mean(pathDB.r_body_R_down(1:M_down_R,:,i));
        r_body_R_up = mean(pathDB.r_body_R_up(1:M_up_R,:,i));
        
        % convert to strokeplane reference frame:
        
        r_st_R_down = Rot_mat*r_body_R_down';
        r_st_R_up = Rot_mat*r_body_R_up';
        
        % conver to strokeplane angles phi and theta:
        
        r_phi_R_down = real(atan2(-r_st_R_down(1),r_st_R_down(2)));
        r_theta_R_down = real(atan2(-r_st_R_down(3),sqrt(r_st_R_down(1)^2+r_st_R_down(2)^2)));
        
        r_phi_R_up = real(atan2(-r_st_R_up(1),r_st_R_up(2)));
        r_theta_R_up = real(atan2(-r_st_R_up(3),sqrt(r_st_R_up(1)^2+r_st_R_up(2)^2)));
        
        %------------------------------------------------------------------
        
        % Plot the normal and tangential strokeplane forces:

        
        scale_c = 8;

        scale_quiv = 1e6;
        
        scale_quiv2 = 5e4;

        hold on
        subplot(3,2,1); plot(radtodeg(phi_L_down),radtodeg(theta_L_down),'r','LineWidth',1.5)
        hold on   
        for k = 1:M_down_L
        plot([radtodeg(phi_L_down(k))-0.75*scale_c*cos(eta_L_down(k)) radtodeg(phi_L_down(k))+0.25*scale_c*cos(eta_L_down(k))],[radtodeg(theta_L_down(k))-0.75*scale_c*sin(eta_L_down(k)) radtodeg(theta_L_down(k))+0.25*scale_c*sin(eta_L_down(k))],'k')
        end            
        quiv_id = quiver(radtodeg(phi_L_down(1:M_down_L))',radtodeg(theta_L_down(1:M_down_L))',scale_quiv*-F_st_L_down(1,1:M_down_L),scale_quiv*-F_st_L_down(3,1:M_down_L),'k');
        adjust_quiver_arrowhead_size(quiv_id,0.5)
        clear quiv_id
        ylabel('\theta [deg]','FontSize',12)
        set(gca, 'YTick', [-20,0,20],'FontSize',12);
        title('Downstroke','FontSize',12)
        axis([-80 100 -30 60])
%         xlabel('\phi [deg]','FontSize',12)
%         set(gca, 'XTick', [-60,0,60],'FontSize',12);
        set(gca, 'XTick',[]);
        hold off
        
        
        subplot(3,2,2); plot(radtodeg(phi_L_up),radtodeg(theta_L_up),'r','LineWidth',1.5)
        hold on
        for k = 1:M_up_L
        plot([radtodeg(phi_L_up(k))-0.75*scale_c*cos(eta_L_up(k)) radtodeg(phi_L_up(k))+0.25*scale_c*cos(eta_L_up(k))],[radtodeg(theta_L_up(k))-0.75*scale_c*sin(eta_L_up(k)) radtodeg(theta_L_up(k))+0.25*scale_c*sin(eta_L_up(k))],'k')
        end
        quiv_id = quiver(radtodeg(phi_L_up(1:M_up_L)'),radtodeg(theta_L_up(1:M_up_L)'),scale_quiv*-F_st_L_up(1,1:M_up_L),scale_quiv*-F_st_L_up(3,1:M_up_L),'k');
        adjust_quiver_arrowhead_size(quiv_id,0.5)
        clear quiv_id
        title('Upstroke','FontSize',12)
        axis([-80 100 -30 60])
%         xlabel('\phi [deg]','FontSize',12)
%         set(gca, 'XTick', [-60,0,60],'FontSize',12);
        set(gca, 'XTick',[]);
        set(gca, 'YTick', []);
        hold off
        
        
        subplot(3,2,3); plot(radtodeg(phi_R_down),radtodeg(theta_R_down),'Color',[0 0.5 0],'LineWidth',1.5)
        hold on
        for k = 1:M_down_R
        plot([radtodeg(phi_R_down(k))-0.75*scale_c*cos(eta_R_down(k)) radtodeg(phi_R_down(k))+0.25*scale_c*cos(eta_R_down(k))],[radtodeg(theta_R_down(k))-0.75*scale_c*sin(eta_R_down(k)) radtodeg(theta_R_down(k))+0.25*scale_c*sin(eta_R_down(k))],'k')
        end
        quiv_id = quiver(radtodeg(phi_R_down(1:M_down_R)'),radtodeg(theta_R_down(1:M_down_R)'),scale_quiv*-F_st_R_down(1,1:M_down_R),scale_quiv*-F_st_R_down(3,1:M_down_R),'k');
        adjust_quiver_arrowhead_size(quiv_id,0.5)
        clear quiv_id
        set(gca, 'XTick',[]);
        ylabel('\theta [deg]','FontSize',12)
        set(gca, 'YTick', [-20,0,20],'FontSize',12);
        axis([-80 100 -30 60])
        hold off
              
        
        subplot(3,2,4); plot(radtodeg(phi_R_up),radtodeg(theta_R_up),'Color',[0 0.5 0],'LineWidth',1.5)
        hold on
        for k = 1:M_up_R
        plot([radtodeg(phi_R_up(k))-0.75*scale_c*cos(eta_R_up(k)) radtodeg(phi_R_up(k))+0.25*scale_c*cos(eta_R_up(k))],[radtodeg(theta_R_up(k))-0.75*scale_c*sin(eta_R_up(k)) radtodeg(theta_R_up(k))+0.25*scale_c*sin(eta_R_up(k))],'k')
        end
        quiv_id = quiver(radtodeg(phi_R_up(1:M_up_R)'),radtodeg(theta_R_up(1:M_up_R)'),scale_quiv*-F_st_R_up(1,1:M_up_R),scale_quiv*-F_st_R_up(3,1:M_up_R),'k');
        adjust_quiver_arrowhead_size(quiv_id,0.5)
        clear quiv_id
        axis([-80 100 -30 60])
        set(gca, 'XTick',[]);
        set(gca, 'YTick',[]);
        hold off      
        
        
        subplot(3,2,5); plot(radtodeg(phi_L_down),-1e6.*F_st_L_down(3,:),'r','LineWidth',1.5)
        hold on
        plot(radtodeg(phi_R_down),-1e6.*F_st_R_down(3,:),'Color',[0 0.5 0],'LineWidth',1.5)
        xlabel('\phi [deg]','FontSize',12) %,'Position',[0 -35 0]
        set(gca, 'XTick', [-60,0,60],'FontSize',12);
        ylabel('Fz \cdot 10^{-6} [N]','FontSize',12)
        set(gca, 'YTick', [0,40,80],'FontSize',12);
        axis([-80 100 -10 80])
        hold off
        
        subplot(3,2,6);plot(radtodeg(phi_L_up),-1e6.*F_st_L_up(3,:),'r','LineWidth',1.5)
        hold on
        plot(radtodeg(phi_R_up),-1e6.*F_st_R_up(3,:),'Color',[0 0.5 0],'LineWidth',1.5)
        xlabel('\phi [deg]','FontSize',12)
        set(gca, 'XTick', [-60,0,60],'FontSize',12);
        set(gca, 'YTick',[]);
        axis([-80 100 -10 80])
        hold off       
        
                
%         subplot(1,2,1); plot(cos(phi_L_down),sin(phi_L_down),'r','LineWidth',1.5)
%         hold on
%         plot(-cos(phi_R_down),sin(phi_R_down),'Color',[0 0.5 0],'LineWidth',1.5)
%         quiv_id1 = quiver(cos(phi_L_down),sin(phi_L_down),2*-F_st_L_down(2,:)',2*-F_st_L_down(1,:)','k');
%         quiv_id2 = quiver(-cos(phi_R_down),sin(phi_R_down),2*-F_st_R_down(2,:)',2*-F_st_R_down(1,:)','k');
%         quiv_id3 = quiver(-r_st_L_down(2)/(wing_l),r_st_L_down(1)/(wing_l),scale_quiv2*F_mean_L_down(2),-scale_quiv2*F_mean_L_down(1),'r');
%         quiv_id4 = quiver(-r_st_R_down(2)/(wing_l),r_st_R_down(1)/(wing_l),scale_quiv2*F_mean_R_down(2),-scale_quiv2*F_mean_R_down(1),'Color',[0 0.5 0]);
%         adjust_quiver_arrowhead_size(quiv_id1,0.5)
%         adjust_quiver_arrowhead_size(quiv_id2,0.5)
%         adjust_quiver_arrowhead_size(quiv_id3,0.5)
%         adjust_quiver_arrowhead_size(quiv_id4,0.5)
%         clear quiv_id1 quiv_id2 quiv_id3 quiv_id4
%         plot([0 -r_st_L_down(2)/(wing_l)],[0 r_st_L_down(1)/(wing_l)],'k')
%         plot([0 -r_st_R_down(2)/(wing_l)],[0 r_st_R_down(1)/(wing_l)],'k')
%         axis([-1.5 1.5 -1.5 1.5])
%         title('Downstroke','FontSize',12)
%         axis equal
%         set(gca, 'XTick',[]);
%         set(gca, 'YTick',[]);
%         %axis off
%         hold off
%         
%         subplot(1,2,2); plot(cos(phi_L_up),sin(phi_L_up),'r','LineWidth',1.5)
%         hold on
%         plot(-cos(phi_R_up),sin(phi_R_up),'Color',[0 0.5 0],'LineWidth',1.5)
%         quiv_id1 = quiver(cos(phi_L_up),sin(phi_L_up),2*-F_st_L_up(2,:)',2*-F_st_L_up(1,:)','k');
%         quiv_id2 = quiver(-cos(phi_R_up),sin(phi_R_up),2*-F_st_R_up(2,:)',2*-F_st_R_up(1,:)','k');
%         quiv_id3 = quiver(-r_st_L_up(2)/(wing_l),r_st_L_up(1)/(wing_l),-scale_quiv2*F_mean_L_up(2),-scale_quiv2*F_mean_L_up(1),'r');
%         quiv_id4 = quiver(-r_st_R_up(2)/(wing_l),r_st_R_up(1)/(wing_l),-scale_quiv2*F_mean_R_up(2),-scale_quiv2*F_mean_R_up(1),'Color',[0 0.5 0]);
%         adjust_quiver_arrowhead_size(quiv_id1,0.5)
%         adjust_quiver_arrowhead_size(quiv_id2,0.5)
%         adjust_quiver_arrowhead_size(quiv_id3,0.5)
%         adjust_quiver_arrowhead_size(quiv_id4,0.5)
%         clear quiv_id1 quiv_id2 quiv_id3 quiv_id4
%         plot([0 -r_st_L_up(2)/(wing_l)],[0 r_st_L_up(1)/(wing_l)],'k')
%         plot([0 -r_st_R_up(2)/(wing_l)],[0 r_st_R_up(1)/(wing_l)],'k')
%         axis([-1.5 1.5 -1.5 1.5])
%         title('Upstroke','FontSize',12)
%         axis equal
%         set(gca, 'XTick',[]);
%         set(gca, 'YTick',[]);
%         %axis off
%         hold off
        
        
        

        
        hold off

        
        clear i j a_L b_L a_R b_R M_down_L M_up_L M_down_R M_up_R phi_L_down theta_L_down eta_L_down phi_L_up theta_L_up eta_L_up phi_R_down theta_R_down eta_R_down phi_R_up theta_R_up eta_R_up F_joint_L_down F_joint_L_up F_st_L_down F_st_L_up F_joint_R_down F_joint_R_up F_st_R_down F_st_R_up
        

        
        
        
        
end
