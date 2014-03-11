function Fly_helicopter_model( Rb, R_strk, xyz, a_xyz, body_model, wing_model, body_scale )


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

            mod_t = Rb'*[ x_mod_b(j,k); y_mod_b(j,k); z_mod_b(j,k)]*body_scale;
                
            x_mod_b_star(j,k) = mod_t(1);
            y_mod_b_star(j,k) = mod_t(2);
            z_mod_b_star(j,k) = mod_t(3);
                
            clear mod_t

        end
    end 

    % Transform the left wing to the current orientation
        
    JL       = body_model.Joint_left;
    JR       = body_model.Joint_right;
    
    JL_t     = xyz+Rb'*JL*body_scale;
    JR_t     = xyz+Rb'*JR*body_scale;
    
    % Create left and right strokeplane representation:
    
    x_mod_strk_L = zeros(20,2);
    y_mod_strk_L = zeros(20,2);
    z_mod_strk_L = zeros(20,2);
    
    x_mod_strk_R = zeros(20,2);
    y_mod_strk_R = zeros(20,2);
    z_mod_strk_R = zeros(20,2);
    
    x_circ = wing_length*cos((pi/4):((0.75*pi)/19):(pi))*body_scale;
    y_circ = wing_length*sin((pi/4):((0.75*pi)/19):(pi))*body_scale;
    
    
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
   
    a_scale = 2e-4;
    
    a_xyz = a_xyz * a_scale;
        
    avg_color   = [0.3 0.3 0.3];
    
    % Plot the body and the average wingbeat:
    
    hold on
    a_wL = surf(JL_t(1)+x_mod_strk_L,JL_t(2)+y_mod_strk_L,JL_t(3)+z_mod_strk_L,'facecolor','b','edgecolor','none','facelighting','phong');
    a_wR = surf(JR_t(1)+x_mod_strk_R,JR_t(2)+y_mod_strk_R,JR_t(3)+z_mod_strk_R,'facecolor','b','edgecolor','none','facelighting','phong');
    alpha(a_wL,0.2)
    alpha(a_wR,0.2)
    surf(xyz(1)+x_mod_b_star,xyz(2)+y_mod_b_star,xyz(3)+z_mod_b_star,'facecolor',avg_color,'edgecolor','k','facelighting','phong');
    quiver3(xyz(1),xyz(2),xyz(3),a_xyz(1),a_xyz(2),a_xyz(3),'Color','b')
    axis equal
    hold off
    
%     set(gca,'Color','k')


end

