function Plot_wingbeats(Lwingtip, Rwingtip, Lwingbeat_loc, Rwingbeat_loc, qL1, qL2, qL3, qL4, qR1, qR2, qR3, qR4, t)

    % Program that plots the wingtip path for a single wingbeat, first the
    % downstroke then the upstroke. Besides the path the measurement
    % locations and the orientation of the wing will be plotted.
    
    
    % Insert a fly picture in matlab
    
    %img = im
    
    
    
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
        end_R = size(wingbeat_loc,1)-1;
    end
    
    
    % Define strokeplane
    
    x_strkpln = -2:0.1:2;
    z_strkpln = x_strkpln*tan((55/180)*pi);
    
    
    
    for i = 1:end_L
        
        a = Lwingbeat_loc(i,1):Lwingbeat_loc(start_L+i,2);
        b = Lwingbeat_loc(start_L+i,2):Lwingbeat_loc(i+1,1);
        
        Lwt_x_down = Lwingtip(a,1);
        Lwt_z_down = Lwingtip(a,3);
        
        Lwt_x_up = Lwingtip(b,1);
        Lwt_z_up = Lwingtip(b,3);
               
        xsi_L_down = real(atan2(2*(qL4(a).*qL1(a)+qL2(a).*qL3(a)),1-2*(qL1(a).^2+qL2(a).^2)));
        theta_L_down = real(asin(2*(qL4(a).*qL2(a)-qL3(a).*qL1(a))));
        phi_L_down = real(atan2(2*(qL4(a).*qL3(a)+qL1(a).*qL2(a)),1-2*(qL2(a).^2+qL3(a).^2)));
        
        xsi_L_up = real(atan2(2*(qL4(b).*qL1(b)+qL2(b).*qL3(b)),1-2*(qL1(b).^2+qL2(b).^2)));
        theta_L_up = real(asin(2*(qL4(b).*qL2(b)-qL3(b).*qL1(b))));
        phi_L_up = real(atan2(2*(qL4(b).*qL3(b)+qL1(b).*qL2(b)),1-2*(qL2(b).^2+qL3(b).^2)));
        
        
        for j = 1:length(a)
            
            T_phi_down = [1 0 0; 0 cos(phi_L_down(j)) sin(phi_L_down(j)); 0 -sin(phi_L_down(j)) cos(phi_L_down(j))];
            T_theta_down = [cos(theta_L_down(j)) 0 -sin(theta_L_down(j)); 0 1 0; sin(theta_L_down(j)) 0 cos(theta_L_down(j))];
            T_xsi_down = [cos(xsi_L_down(j)) -sin(xsi_L_down(j)) 0; sin(xsi_L_down(j)) cos(xsi_L_down(j)) 0; 0 0 1];
            
            
            
            quiv_x_down(j) = cos(theta_L_down(j));
            quiv_z_down(j) = -sin(theta_L_down(j));
            
            clear Vect_down
            
        end
            
        for j = 1:length(b)
            
            T_phi_up = [1 0 0; 0 cos(phi_L_up(j)) sin(phi_L_up(j)); 0 -sin(phi_L_up(j)) cos(phi_L_up(j))];
            T_theta_up = [cos(theta_L_up(j)) 0 -sin(theta_L_up(j)); 0 1 0; sin(theta_L_up(j)) 0 cos(theta_L_up(j))];
            T_xsi_up = [cos(xsi_L_up(j)) -sin(xsi_L_up(j)) 0; sin(xsi_L_up(j)) cos(xsi_L_up(j)) 0; 0 0 1];
            
                     
            quiv_x_up(j) = cos(theta_L_up(j));
            quiv_z_up(j) = -sin(theta_L_up(j));
            
            clear Vect_up
            
        end
        
        
        figure()
        plot(-Lwt_x_down,-Lwt_z_down,'--r',-Lwt_x_up,-Lwt_z_up,'--r')
        hold on
        plot(-Lwt_x_down-0.1.*quiv_x_down', -Lwt_z_down-0.1.*quiv_z_down','.k')
        plot(-Lwt_x_up-0.1.*quiv_x_up', -Lwt_z_up-0.1.*quiv_z_up','.k')
        quiver(-Lwt_x_down,-Lwt_z_down,-quiv_x_down',-quiv_z_down',0.1,'k', 'ShowArrowHead','off')
        quiver(-Lwt_x_up,-Lwt_z_up,-quiv_x_up',-quiv_z_up',0.1,'k', 'ShowArrowHead','off')
        quiver(-Lwt_x_down,-Lwt_z_down,quiv_x_down',quiv_z_down',0.2,'k', 'ShowArrowHead','off')
        quiver(-Lwt_x_up,-Lwt_z_up,quiv_x_up',quiv_z_up',0.2,'k','ShowArrowHead','off')
        plot(x_strkpln,z_strkpln,':')
        quiver(0,0,-1,0,3,'--b')
        quiver(0,0,0,-1,3,'--b')
        axis equal
        
        clear Lwt_x_down Lwt_z_down Lwt_x_up Lwt_z_up quiv_x_down quiv_z_down quiv_x_up quiv_z_up
        
    end
    
    clear a b quiv_x_down quiv_x_up quiv_z_down quiv_z_up
    
    for i = 1:end_R
        
        a = Rwingbeat_loc(i,1):Rwingbeat_loc(start_R+i,2);
        b = Rwingbeat_loc(start_R+i,2):Rwingbeat_loc(i+1,1);
        
        Rwt_x_down = Rwingtip(a,1);
        Rwt_z_down = Rwingtip(a,3);
        
        Rwt_x_up = Rwingtip(b,1);
        Rwt_z_up = Rwingtip(b,3);
               
        xsi_R_down = real(atan2(2*(qR4(a).*qR1(a)+qR2(a).*qR3(a)),1-2*(qR1(a).^2+qR2(a).^2)));
        theta_R_down = real(asin(2*(qR4(a).*qR2(a)-qR3(a).*qR1(a))));
        phi_R_down = real(atan2(2*(qR4(a).*qR3(a)+qR1(a).*qR2(a)),1-2*(qR2(a).^2+qR3(a).^2)));
        
        xsi_R_up = real(atan2(2*(qR4(b).*qR1(b)+qR2(b).*qR3(b)),1-2*(qR1(b).^2+qR2(b).^2)));
        theta_R_up = real(asin(2*(qR4(b).*qR2(b)-qR3(b).*qR1(b))));
        phi_R_up = real(atan2(2*(qR4(b).*qR3(b)+qR1(b).*qR2(b)),1-2*(qR2(b).^2+qR3(b).^2)));
        
        
        for j = 1:length(a)
            
            T_phi_down = [1 0 0; 0 cos(phi_R_down(j)) sin(phi_R_down(j)); 0 -sin(phi_R_down(j)) cos(phi_R_down(j))];
            T_theta_down = [cos(theta_R_down(j)) 0 -sin(theta_R_down(j)); 0 1 0; sin(theta_R_down(j)) 0 cos(theta_R_down(j))];
            T_xsi_down = [cos(xsi_R_down(j)) -sin(xsi_R_down(j)) 0; sin(xsi_R_down(j)) cos(xsi_R_down(j)) 0; 0 0 1];
            
            
            
            quiv_x_down(j) = cos(theta_R_down(j));
            quiv_z_down(j) = -sin(theta_R_down(j));

            
            clear Vect_down
            
        end
            
        for j = 1:length(b)
            
            T_phi_up = [1 0 0; 0 cos(phi_R_up(j)) sin(phi_R_up(j)); 0 -sin(phi_R_up(j)) cos(phi_R_up(j))];
            T_theta_up = [cos(theta_R_up(j)) 0 -sin(theta_R_up(j)); 0 1 0; sin(theta_R_up(j)) 0 cos(theta_R_up(j))];
            T_xsi_up = [cos(xsi_R_up(j)) -sin(xsi_R_up(j)) 0; sin(xsi_R_up(j)) cos(xsi_R_up(j)) 0; 0 0 1];
            
                     
            quiv_x_up(j) = cos(theta_R_up(j));
            quiv_z_up(j) = -sin(theta_R_up(j));

            
            clear Vect_up
            
        end
        
        
        figure()
        plot(-Rwt_x_down,-Rwt_z_down,'--r',-Rwt_x_up,-Rwt_z_up,'--r')
        hold on
        plot(-Rwt_x_down-0.1.*quiv_x_down', -Rwt_z_down-0.1.*quiv_z_down','.k')
        plot(-Rwt_x_up-0.1.*quiv_x_up', -Rwt_z_up-0.1.*quiv_z_up','.k')
        quiver(-Rwt_x_down,-Rwt_z_down,-quiv_x_down',-quiv_z_down',0.1, 'k', 'ShowArrowHead','off')
        quiver(-Rwt_x_up,-Rwt_z_up,-quiv_x_up',-quiv_z_up',0.1, 'k', 'ShowArrowHead','off')
        quiver(-Rwt_x_down,-Rwt_z_down,quiv_x_down',quiv_z_down',0.2, 'k', 'ShowArrowHead','off')
        quiver(-Rwt_x_up,-Rwt_z_up,quiv_x_up',quiv_z_up',0.2, 'k', 'ShowArrowHead','off')
        plot(x_strkpln,z_strkpln,':')
        quiver(0,0,-1,0,3,'--b')
        quiver(0,0,0,-1,3,'--b')
        axis equal
        
        
        clear Rwt_x_down Rwt_z_down Rwt_x_up Rwt_z_up quiv_x_down quiv_z_down quiv_x_up quiv_z_up
        
    end
    
    
    

end

