function Fly_plot_3D( body_pos, Rb, RL, RR, body_model, wing_model )

    % Plot fly wings and body in 3D:
    
    x_mod_b = body_model.x_mod;
    y_mod_b = body_model.y_mod;
    z_mod_b = body_model.z_mod;
    
       
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
    
    JL_t     = body_pos+Rb*JL;
    JR_t     = body_pos+Rb*JR;
    
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

            x_mod_wL_star(j,k) = body_pos(1)+mod_t_L(1);
            y_mod_wL_star(j,k) = body_pos(2)+mod_t_L(2);
            z_mod_wL_star(j,k) = body_pos(3)+mod_t_L(3);
            
            mod_t_R = Rb*(JR+RR'*[ x_mod_wR(j,k)-JR(1); y_mod_wR(j,k)-JR(2); z_mod_wR(j,k)-JR(3)]);
                
            x_mod_wR_star(j,k) = body_pos(1)+mod_t_R(1);
            y_mod_wR_star(j,k) = body_pos(2)+mod_t_R(2);
            z_mod_wR_star(j,k) = body_pos(3)+mod_t_R(3);
                
            clear mod_t_L mod_t_R

        end
    end 
        
    x_axis = Rb*[3; 0; 0];
    y_axis = Rb*[0; 3; 0];
    z_axis = Rb*[0; 0; 3];

    figure()
    hold on
    surf(body_pos(1)+x_mod_b_star,body_pos(2)+y_mod_b_star,body_pos(3)+z_mod_b_star,'facecolor',[0.3 0.3 0.3],'edgecolor','k','facelighting','phong');
    surf(x_mod_wL_star,y_mod_wL_star,z_mod_wL_star,'facecolor',[0.3 0.3 0.3],'edgecolor','k','facelighting','phong');
    surf(x_mod_wR_star,y_mod_wR_star,z_mod_wR_star,'facecolor',[0.3 0.3 0.3],'edgecolor','k','facelighting','phong');
%     plot3(JL_t(1),JL_t(2),JL_t(3),'o');
%     plot3(JR_t(1),JR_t(2),JR_t(3),'o');
%     plot3([0 x_axis(1)],[0 x_axis(2)],[0 x_axis(3)],'r')
%     plot3([0 y_axis(1)],[0 y_axis(2)],[0 y_axis(3)],'r')
%     plot3([0 z_axis(1)],[0 z_axis(2)],[0 z_axis(3)],'r')
    axis equal
    hold off
    
end

