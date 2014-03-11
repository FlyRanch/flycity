function Blade_element_plot( Rb, RL, RR, body_model, wing_model )

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
    
   
    wing_line_L = Rb'*RL'*[0; -L_w; 0];
    wing_line_R = Rb'*RR'*[0; L_w; 0];

    
    e_color = 'k';
    f_color = [0.3 0.3 0.3];

    hold on
    surf(x_mod_b_star,y_mod_b_star,z_mod_b_star,'facecolor',f_color,'edgecolor',e_color,'facelighting','phong');
    surf(x_mod_wL_star,y_mod_wL_star,z_mod_wL_star,'facecolor','k','edgecolor','none','facelighting','phong');
    surf(x_mod_wR_star,y_mod_wR_star,z_mod_wR_star,'facecolor','k','edgecolor','none','facelighting','phong');
    alpha(0.5)
    plot3([JL_t(1) JL_t(1)+wing_line_L(1)],[JL_t(2) JL_t(2)+wing_line_L(2)],[JL_t(3) JL_t(3)+wing_line_L(3)],'k')
    plot3([JR_t(1) JR_t(1)+wing_line_R(1)],[JR_t(2) JR_t(2)+wing_line_R(2)],[JR_t(3) JR_t(3)+wing_line_R(3)],'k')
    plot3(JL_t(1),JL_t(2),JL_t(3),'o','Color','k');
    plot3(JR_t(1),JR_t(2),JR_t(3),'o','Color','k');
    xlim([-3 2])
    ylim([-2.5 2.5])
    zlim([-2 3])
    axis equal
    hold off
end
