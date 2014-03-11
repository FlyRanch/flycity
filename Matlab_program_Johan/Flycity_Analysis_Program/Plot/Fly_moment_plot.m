function Fly_moment_plot( Rb, R_strk, xyz, wb, wb_raw, MA_b, MI_vel_b, MI_acc_b, body_model, wing_model )


    % Plot fly wings and body in 3D:
    
    x_mod_b = body_model.x_mod;
    y_mod_b = body_model.y_mod;
    z_mod_b = body_model.z_mod;
       
    n_b = size(x_mod_b,1);
    m_b = size(x_mod_b,2);
    
    x_mod_b_star = zeros(n_b,m_b);
    y_mod_b_star = zeros(n_b,m_b);
    z_mod_b_star = zeros(n_b,m_b);
    
    wing_length = wing_model.length;
    
    % Transform the body to current orientation

    for j = 1:n_b
        for k = 1:m_b

            mod_t = Rb'*[ x_mod_b(j,k); y_mod_b(j,k); z_mod_b(j,k)];
                
            x_mod_b_star(j,k) = mod_t(1);
            y_mod_b_star(j,k) = mod_t(2);
            z_mod_b_star(j,k) = mod_t(3);
                
            clear mod_t

        end
    end 

    % Transform the left wing to the current orientation
        
    JL       = body_model.Joint_left;
    JR       = body_model.Joint_right;
    
    JL_t     = xyz+Rb'*JL;
    JR_t     = xyz+Rb'*JR;
    
    % Create left and right strokeplane representation:
    
    x_mod_strk_L = zeros(20,2);
    y_mod_strk_L = zeros(20,2);
    z_mod_strk_L = zeros(20,2);
    
    x_mod_strk_R = zeros(20,2);
    y_mod_strk_R = zeros(20,2);
    z_mod_strk_R = zeros(20,2);
    
    x_circ = wing_length*cos((pi/4):((0.75*pi)/19):(pi));
    y_circ = wing_length*sin((pi/4):((0.75*pi)/19):(pi));
    
    
    for i = 1:20
            

            L_plane      = Rb'*R_strk'*[x_circ(i); -y_circ(i); 0];
            R_plane      = Rb'*R_strk'*[x_circ(i); y_circ(i); 0];
                    

            xyz_mod_L       = [L_plane(1); L_plane(2); L_plane(3)];
            xyz_mod_R       = [R_plane(1); R_plane(2); R_plane(3)];

            
            x_mod_strk_L(i,:) = [xyz_mod_L(1) 0];
            y_mod_strk_L(i,:) = [xyz_mod_L(2) 0];
            z_mod_strk_L(i,:) = [xyz_mod_L(3) 0];
    
            x_mod_strk_R(i,:) = [xyz_mod_R(1) 0];
            y_mod_strk_R(i,:) = [xyz_mod_R(2) 0];
            z_mod_strk_R(i,:) = [xyz_mod_R(3) 0];

        
    end
    
    
    
    % Convert the forces to the inertial reference frame:
    
    w_scale     = 1e-2;
    M_scale     = 5e5;
    
    MA          = Rb'*MA_b*M_scale;
    MI_vel      = Rb'*MI_vel_b*M_scale;
    MI_acc      = Rb'*MI_acc_b*M_scale;

    w_xyz       = Rb'*wb*w_scale;
    
%     a_xyz_raw   = Rb*ab_raw*a_scale;
    w_xyz_raw   = Rb'*wb_raw*w_scale;
        
    avg_color   = [0.5 0.5 0.5];
    
    % Plot the body and the average wingbeat:
    
    hold on
    
    surf(xyz(1)+x_mod_b_star,xyz(2)+y_mod_b_star,xyz(3)+z_mod_b_star,'facecolor','k','edgecolor',avg_color,'facelighting','phong');
    alpha(0.1)
    surf(JL_t(1)+x_mod_strk_L,JL_t(2)+y_mod_strk_L,JL_t(3)+z_mod_strk_L,'facecolor','k','edgecolor',avg_color)
    alpha(0.1)
    surf(JR_t(1)+x_mod_strk_R,JR_t(2)+y_mod_strk_R,JR_t(3)+z_mod_strk_R,'facecolor','k','edgecolor',avg_color)
    alpha(0.1)
    plot3(JL_t(1),JL_t(2),JL_t(3),'o','Color',avg_color);
    plot3(JR_t(1),JR_t(2),JR_t(3),'o','Color',avg_color);
    quiver3(xyz(1),xyz(2),xyz(3),MA(1),MA(2),MA(3),'Color','b')
    quiver3(xyz(1),xyz(2),xyz(3),MI_vel(1),MI_vel(2),MI_vel(3),'Color','r')
    quiver3(xyz(1),xyz(2),xyz(3),MI_acc(1),MI_acc(2),MI_acc(3),'Color','y')
    quiver3(xyz(1),xyz(2),xyz(3),0,0,-2,'Color',avg_color)
    quiver3(xyz(1),xyz(2),xyz(3),w_xyz(1),w_xyz(2),w_xyz(3),'Color','m')
    quiver3(xyz(1),xyz(2),xyz(3),w_xyz_raw(1),w_xyz_raw(2),w_xyz_raw(3),'Color','c')
    set(gca,'Xcolor',avg_color);
    set(gca,'Ycolor',avg_color);
    set(gca,'Zcolor',avg_color);
    axis equal
    hold off
    
    set(gca,'Color','k')

end


