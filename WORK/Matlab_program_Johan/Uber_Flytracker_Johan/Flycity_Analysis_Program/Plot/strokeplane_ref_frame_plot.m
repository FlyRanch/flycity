function strokeplane_ref_frame_plot( Rb, RL, RR, body_model, wing_model )

    % Plot fly wings and body in 3D:
    
    x_mod_b = body_model.x_mod;
    y_mod_b = body_model.y_mod;
    z_mod_b = body_model.z_mod;
    
    R_strk = body_model.R_strk;
    
    L_w = wing_model.length;
       
    n_b = size(x_mod_b,1);
    m_b = size(x_mod_b,2);
    
    x_mod_b_star = zeros(n_b,m_b);
    y_mod_b_star = zeros(n_b,m_b);
    z_mod_b_star = zeros(n_b,m_b);
    
    % Transform the body to current orientation

    for j = 1:n_b
        for k = 1:m_b

            mod_t = Rb*[ x_mod_b(j,k); y_mod_b(j,k); z_mod_b(j,k)];
                
            x_mod_b_star(j,k) = mod_t(1);
            y_mod_b_star(j,k) = mod_t(2);
            z_mod_b_star(j,k) = mod_t(3);
                
            clear mod_t

        end
    end 
    
    % Transform the left wing to the current orientation
    
    x_mod_wL = wing_model.x_mod_L;
    y_mod_wL = wing_model.y_mod_L;
    z_mod_wL = wing_model.z_mod_L;
    
    x_mod_wR = wing_model.x_mod_R;
    y_mod_wR = wing_model.y_mod_R;
    z_mod_wR = wing_model.z_mod_R;
    
    JL       = body_model.Joint_left;
    JR       = body_model.Joint_right;
    
    JL_t     = Rb*JL;
    JR_t     = Rb*JR;
    
    n_w = size(x_mod_wL,1);
    m_w = size(x_mod_wL,2);
        
    x_mod_wL_star = zeros(n_w,m_w);
    y_mod_wL_star = zeros(n_w,m_w);
    z_mod_wL_star = zeros(n_w,m_w);
    
    x_mod_wR_star = zeros(n_w,m_w);
    y_mod_wR_star = zeros(n_w,m_w);
    z_mod_wR_star = zeros(n_w,m_w);
    
    for j = 1:n_w
        for k = 1:m_w

            mod_t_L = Rb*(JL+RL'*[ x_mod_wL(j,k)-JL(1); y_mod_wL(j,k)-JL(2); z_mod_wL(j,k)-JL(3)]);

            x_mod_wL_star(j,k) = mod_t_L(1);
            y_mod_wL_star(j,k) = mod_t_L(2);
            z_mod_wL_star(j,k) = mod_t_L(3);
            
            mod_t_R = Rb*(JR+RR'*[ x_mod_wR(j,k)-JR(1); y_mod_wR(j,k)-JR(2); z_mod_wR(j,k)-JR(3)]);
                
            x_mod_wR_star(j,k) = mod_t_R(1);
            y_mod_wR_star(j,k) = mod_t_R(2);
            z_mod_wR_star(j,k) = mod_t_R(3);
                
            clear mod_t_L mod_t_R

        end
    end 
    
    % Create left and right strokeplane representation:
    
    x_mod_strk_L = zeros(100,2);
    y_mod_strk_L = zeros(100,2);
    z_mod_strk_L = zeros(100,2);
    
    x_mod_strk_R = zeros(100,2);
    y_mod_strk_R = zeros(100,2);
    z_mod_strk_R = zeros(100,2);
    
    x_circ = L_w*cos((pi/4):((0.75*pi)/99):(pi));
    y_circ = L_w*sin((pi/4):((0.75*pi)/99):(pi));
    
    
    for i = 1:100
            

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
    
        
    x_axis_strk = R_strk*Rb*[2; 0; 0];
    y_axis_strk = R_strk*Rb*[0; 2; 0];
    z_axis_strk = R_strk*Rb*[0; 0; 2];
    
    x_axis_b = Rb*[4.7; 0; 0];
    
    wing_line_L = Rb'*RL'*[0; -L_w; 0];
    wing_line_R = Rb'*RR'*[0; L_w; 0];

    
    e_color = 'k';
    f_color = [0.3 0.3 0.3];

    hold on
    surf(JL_t(1)+x_mod_strk_L,JL_t(2)+y_mod_strk_L,JL_t(3)+z_mod_strk_L,'facecolor','k','edgecolor','none','facelighting','phong')
    surf(JR_t(1)+x_mod_strk_R,JR_t(2)+y_mod_strk_R,JR_t(3)+z_mod_strk_R,'facecolor','k','edgecolor','none','facelighting','phong')
    alpha(0.3)
    surf(x_mod_b_star,y_mod_b_star,z_mod_b_star,'facecolor',f_color,'edgecolor',e_color,'facelighting','phong');
    surf(x_mod_wL_star,y_mod_wL_star,z_mod_wL_star,'facecolor',f_color,'edgecolor','k','facelighting','phong');
    surf(x_mod_wR_star,y_mod_wR_star,z_mod_wR_star,'facecolor',f_color,'edgecolor','k','facelighting','phong');
    plot3([JL_t(1) JL_t(1)+wing_line_L(1)],[JL_t(2) JL_t(2)+wing_line_L(2)],[JL_t(3) JL_t(3)+wing_line_L(3)],'r')
    plot3([JR_t(1) JR_t(1)+wing_line_R(1)],[JR_t(2) JR_t(2)+wing_line_R(2)],[JR_t(3) JR_t(3)+wing_line_R(3)],'r')
%     alpha(0.5)
    plot3(JL_t(1),JL_t(2),JL_t(3),'o','Color','k');
    plot3(JR_t(1),JR_t(2),JR_t(3),'o','Color','k');
    quiver3(JL_t(1),0,JL_t(3),x_axis_strk(1),x_axis_strk(2),x_axis_strk(3),'Color','b','LineWidth',1.0)
    quiver3(JL_t(1),0,0,x_axis_b(1),x_axis_b(2),x_axis_b(3),'Color','k','LineWidth',1.0)
    quiver3(JL_t(1),0,JL_t(3),y_axis_strk(1),y_axis_strk(2),y_axis_strk(3),'Color','b','LineWidth',1.0)
    quiver3(JL_t(1),0,JL_t(3),z_axis_strk(1),z_axis_strk(2),z_axis_strk(3),'Color','b','LineWidth',1.0)
    xlim([-3 4])
    ylim([-3.5 3.5])
    zlim([-3 4])
    axis equal
    hold off
end