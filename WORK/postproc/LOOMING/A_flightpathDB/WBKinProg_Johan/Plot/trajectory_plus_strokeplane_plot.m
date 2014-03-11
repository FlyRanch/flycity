function trajectory_plus_strokeplane_plot(settings,pathDB,seq_nr,save_on_off,fig_nr)

    % Plot the strokeplane for every wingbeat (represented as a circle) in
    % 3 D along the trajectory:
    
    
    addpath(char(settings.path_names(9)))
    
    
    
    % start and stop point for the measurements:
    
    start = find(isnan(pathDB.x(:,seq_nr))==0, 1 );
    stop = find(isnan(pathDB.x(:,seq_nr))==0, 1, 'last' );
    
    
    % create a circle:
    
    r_circle = 0.4;
    theta_circle = 0:(0.1*pi):2*pi;
    circle_global = [ r_circle*cos(theta_circle); r_circle*sin(theta_circle); zeros(1,21)];
    
    theta_circle_L = pi:(0.1*pi):2*pi;
    circle_global_L = [ r_circle*cos(theta_circle_L); r_circle*sin(theta_circle_L); zeros(1,11)];
    
    theta_circle_R = 0:(0.1*pi):pi;
    circle_global_R = [ r_circle*cos(theta_circle_R); r_circle*sin(theta_circle_R); zeros(1,11)];
    
    end_wingbeat = find(isnan(pathDB.wingbeat_time(:,1,seq_nr))==0, 1, 'last' );
    
    circle_strkpln = zeros(3,21,end_wingbeat);
    
    circle_strkpln_L = zeros(3,11,end_wingbeat);
    
    circle_strkpln_R = zeros(3,11,end_wingbeat);
    
    pos_strkpln = zeros(3,end_wingbeat);
    
    q_body = zeros(4,end_wingbeat);
    
    pos_joint = zeros(3,end_wingbeat);
    
    jpL = [pathDB.joint_pos_L(1,seq_nr); 0; pathDB.joint_pos_L(3,seq_nr)];
    
    for m = 1:end_wingbeat
        
        end_temp = find(isnan(pathDB.wingbeat_time(m,:,seq_nr))==0, 1, 'last' );
        
        frame_wb = round((pathDB.wingbeat_time(m,end_temp,seq_nr)+pathDB.wingbeat_time(m,1,seq_nr))/2);
        
        % Get average body quaternion during the wingbeat:
        
        n = (start-1+pathDB.wingbeat_time(m,1,seq_nr)):1:(start-1+pathDB.wingbeat_time(m,end_temp,seq_nr));
        
        q_body(:,m) = q_avg(pathDB.qb1_filt(n,seq_nr),pathDB.qb2_filt(n,seq_nr),pathDB.qb3_filt(n,seq_nr),pathDB.qb4_filt(n,seq_nr));
        
        R_body = quat2matNEW([q_body(1,m); q_body(2,m); q_body(3,m); q_body(4,m)]);
        
        beta = (55/180)*pi;
        
        R_stroke = [cos(beta) 0 -sin(beta); ...
                    0 1 0; ...
                    sin(beta) 0 cos(beta)];
                
        q_left = q_avg(pathDB.qL1_filt2(n,seq_nr),pathDB.qL2_filt2(n,seq_nr),pathDB.qL3_filt2(n,seq_nr),pathDB.qL4_filt2(n,seq_nr));
        
        q_right = q_avg(pathDB.qR1_filt2(n,seq_nr),pathDB.qR2_filt2(n,seq_nr),pathDB.qR3_filt2(n,seq_nr),pathDB.qR4_filt2(n,seq_nr));
        
% %         phi_L = atan2(2*(q_left(4).*q_left(1)+q_left(2).*q_left(3)),(1-2*(q_left(1).^2+q_left(2)^2)));
%         phi_L = 0;
%         theta_L = asin(2*(q_left(4).*q_left(2)-q_left(3).*q_left(1)));
% %         xsi_L = atan2(2*(q_left(4).*q_left(3)+q_left(1).*q_left(2)),1-2*(q_left(2).^2+q_left(3).^2));
%         xsi_L = 0;
%         
% %         phi_R = atan2(2*(q_right(4).*q_right(1)+q_right(2).*q_right(3)),(1-2*(q_right(1).^2+q_right(2)^2)));
%         phi_R = 0;
%         theta_R = asin(2*(q_right(4).*q_right(2)-q_right(3).*q_right(1)));
% %         xsi_R = atan2(2*(q_right(4).*q_right(3)+q_right(1).*q_right(2)),1-2*(q_right(2).^2+q_right(3).^2));   
%         xsi_R = 0;
%         
%         q_l = [sin(phi_L/2)*cos(theta_L/2)*cos(xsi_L/2)-cos(phi_L/2)*sin(theta_L/2)*sin(xsi_L/2); ...
%                cos(phi_L/2)*sin(theta_L/2)*cos(xsi_L/2)+sin(phi_L/2)*cos(theta_L/2)*sin(xsi_L/2); ...
%                cos(phi_L/2)*cos(theta_L/2)*sin(xsi_L/2)-sin(phi_L/2)*sin(theta_L/2)*cos(xsi_L/2); ...
%                cos(phi_L/2)*cos(theta_L/2)*cos(xsi_L/2)+sin(phi_L/2)*sin(theta_L/2)*sin(xsi_L/2)];
%           
%         q_r = [sin(phi_R/2)*cos(theta_R/2)*cos(xsi_R/2)-cos(phi_R/2)*sin(theta_R/2)*sin(xsi_R/2); ...
%                cos(phi_R/2)*sin(theta_R/2)*cos(xsi_R/2)+sin(phi_R/2)*cos(theta_R/2)*sin(xsi_R/2); ...
%                cos(phi_R/2)*cos(theta_R/2)*sin(xsi_R/2)-sin(phi_R/2)*sin(theta_R/2)*cos(xsi_R/2); ...
%                cos(phi_R/2)*cos(theta_R/2)*cos(xsi_R/2)+sin(phi_R/2)*sin(theta_R/2)*sin(xsi_R/2)];


        phi_L = pathDB.phi_L(n,seq_nr);
        
        theta_L = pathDB.theta_L(n,seq_nr);
        
        phi_R = pathDB.phi_R(n,seq_nr);
        
        theta_R = pathDB.theta_R(n,seq_nr);
        
        [eta_L,theta_mean_L,smL,sbL]=lsqfitgm(phi_L,theta_L);
        
        [eta_R,theta_mean_R,smR,sbR]=lsqfitgm(phi_R,theta_R);
        
        R_eta_L = [cos(eta_L) 0 -sin(eta_L); ...
                    0 1 0; ...
                    sin(eta_L) 0 cos(eta_L)];
        
        R_theta_L = [1 0 0; ...
                     0 cos(theta_mean_L) -sin(theta_mean_L); ...
                     0 sin(theta_mean_L) cos(theta_mean_L)];
        
        R_eta_R = [cos(eta_R) 0 -sin(eta_R); ...
                    0 1 0; ...
                    sin(eta_R) 0 cos(eta_R)];
        
        R_theta_R = [1 0 0; ...
                     0 cos(theta_mean_R) -sin(theta_mean_R); ...
                     0 sin(theta_mean_R) cos(theta_mean_R)];

        
%         R_L = quat2matNEW([q_l(1); q_l(2); q_l(3); q_l(4)]);
%         
%         R_R = quat2matNEW([q_r(1); q_r(2); q_r(3); q_r(4)]);
        
        circle_strkpln_L(:,1:11,m) = R_body*R_theta_L*R_eta_L*R_stroke*circle_global_L(:,1:11);
        
        circle_strkpln_R(:,1:11,m) = R_body*R_theta_R*R_eta_R*R_stroke*circle_global_R(:,1:11);
        
        circle_strkpln(:,1:21,m) = R_body*R_stroke*circle_global(:,1:21);

%         for k = 1:21
%             
%             circle_strkpln(:,k,m) = R_body*R_stroke*circle_global(:,k);
%             
%         end
        
        pos_joint(:,m) = R_body*jpL;
        
        pos_strkpln(:,m) = [pathDB.x_filt(start-1+frame_wb,seq_nr); pathDB.y_filt(start-1+frame_wb,seq_nr); pathDB.z_filt(start-1+frame_wb,seq_nr)];
        
    end
    
    
       
    figure(fig_nr)
    plot3(pathDB.x_filt(start:stop,seq_nr),pathDB.y_filt(start:stop,seq_nr),pathDB.z_filt(start:stop,seq_nr),'g')
    hold on
    for j = 1:end_wingbeat
    scale_body = 0.2;
    plot3(pos_strkpln(1,j)+scale_body*pos_joint(1,j)+circle_strkpln(1,:,j),pos_strkpln(2,j)+scale_body*pos_joint(2,j)+circle_strkpln(2,:,j),pos_strkpln(3,j)+scale_body*pos_joint(3,j)+circle_strkpln(3,:,j),'b')
    [ x, y, z ] = load_body_model(settings, seq_nr, [pos_strkpln(1,j)/scale_body pos_strkpln(2,j)/scale_body pos_strkpln(3,j)/scale_body q_body(1,j) q_body(2,j) q_body(3,j) q_body(4,j) 0 0 0 1 0 0 0 1]);
    surf(scale_body.*x{1},scale_body.*y{1},scale_body.*z{1},'facecolor',[0.5 0.5 0.5],'edgecolor','k','facelighting','phong');
    plot3(pos_strkpln(1,j)+scale_body*pos_joint(1,j)+circle_strkpln_L(1,:,j),pos_strkpln(2,j)+scale_body*pos_joint(2,j)+circle_strkpln_L(2,:,j),pos_strkpln(3,j)+scale_body*pos_joint(3,j)+circle_strkpln_L(3,:,j),'r')
    plot3(pos_strkpln(1,j)+scale_body*pos_joint(1,j)+circle_strkpln_R(1,:,j),pos_strkpln(2,j)+scale_body*pos_joint(2,j)+circle_strkpln_R(2,:,j),pos_strkpln(3,j)+scale_body*pos_joint(3,j)+circle_strkpln_R(3,:,j),'r')
    end
    axis equal
    hold off
    
    
    % Save plots
    
    if save_on_off == 1
    
    saveas(fig_nr, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/stkpln_plus_path'], 'fig')
    
    saveas(fig_nr, [char(settings.plot_folders(2)) '/flight_path/stkpln_plus_path_' int2str(seq_nr)], 'fig')
    
    end

    
end

