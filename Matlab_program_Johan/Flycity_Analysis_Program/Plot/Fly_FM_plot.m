function Fly_FM_plot( Rb, R_strk,body_model, wing_model, F_sim, kine_sim, kine_raw )

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
    
    JL_t     = Rb'*JL;
    JR_t     = Rb'*JR;
    
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
    
    a_scale     = 2e-4;
    F_scale     = 2e5;
    
    FA          = F_sim.FA*F_scale;
    FI_vel      = F_sim.FI_vel*F_scale;
    FI_acc      = F_sim.FI_acc*F_scale;
    Fg          = F_sim.Fg*F_scale;
    
    w_dot_scale = 4e-4;
    M_scale     = 1e6;
    
    MA          = F_sim.MA*M_scale;
    MI_vel      = F_sim.MI_vel*M_scale;
    MI_acc      = F_sim.MI_acc*M_scale;

    a_sim       = kine_sim.a_xyz*a_scale;
    w_dot_sim   = kine_sim.w_dot*w_dot_scale;
    
    a_raw       = kine_raw.a_xyz*a_scale;
    w_dot_raw   = kine_raw.w_dot*w_dot_scale;
        
    avg_color   = [0.7 0.7 0.7];
    
    % Plot the body and the average wingbeat:
    
%     figure(n)
%     hold on
%     surf(x_mod_b_star,y_mod_b_star,z_mod_b_star,'facecolor','k','edgecolor',avg_color,'facelighting','phong');
%     alpha(0.1)
%     surf(JL_t(1)+x_mod_strk_L,JL_t(2)+y_mod_strk_L,JL_t(3)+z_mod_strk_L,'facecolor','k','edgecolor',avg_color)
%     alpha(0.1)
%     surf(JR_t(1)+x_mod_strk_R,JR_t(2)+y_mod_strk_R,JR_t(3)+z_mod_strk_R,'facecolor','k','edgecolor',avg_color)
%     alpha(0.1)
%     plot3(JL_t(1),JL_t(2),JL_t(3),'o','Color',avg_color);
%     plot3(JR_t(1),JR_t(2),JR_t(3),'o','Color',avg_color);
%     quiver3(0,0,0,FA(1),FA(2),FA(3),'Color','b')
%     quiver3(0,0,0,FI_vel(1),FI_vel(2),FI_vel(3),'Color','r')
%     quiver3(0,0,0,FI_acc(1),FI_acc(2),FI_acc(3),'Color','y')
%     quiver3(0,0,0,Fg(1),Fg(2),Fg(3),'Color','g')
%     quiver3(0,0,0,0,0,-2,'Color',avg_color)
%     quiver3(0,0,0,a_sim(1),a_sim(2),a_sim(3),'Color','m')
%     quiver3(0,0,0,a_raw(1),a_raw(2),a_raw(3),'Color','c')
%     set(gca,'Xcolor',avg_color);
%     set(gca,'Ycolor',avg_color);
%     set(gca,'Zcolor',avg_color);
%     axis equal
%     hold off  
%     set(gca,'Color','k')
%     set(n,'Color','k')
%     
%     n = n+1;
    

    hold on
    subplot(1,2,1); hold on
    surf(x_mod_b_star,y_mod_b_star,z_mod_b_star,'facecolor','k','edgecolor',avg_color,'facelighting','phong');
    alpha(0)
    surf(JL_t(1)+x_mod_strk_L,JL_t(2)+y_mod_strk_L,JL_t(3)+z_mod_strk_L,'facecolor','k','edgecolor',avg_color)
    alpha(0)
    surf(JR_t(1)+x_mod_strk_R,JR_t(2)+y_mod_strk_R,JR_t(3)+z_mod_strk_R,'facecolor','k','edgecolor',avg_color)
    alpha(0)
    plot3(JL_t(1),JL_t(2),JL_t(3),'o','Color',avg_color);
    plot3(JR_t(1),JR_t(2),JR_t(3),'o','Color',avg_color);
    quiver3(0,0,0,FA(1),FA(2),FA(3),'Color','b')
    quiver3(0,0,0,FI_vel(1),FI_vel(2),FI_vel(3),'Color','r')
    quiver3(0,0,0,FI_acc(1),FI_acc(2),FI_acc(3),'Color','y')
    quiver3(0,0,0,Fg(1),Fg(2),Fg(3),'Color','g')
%     quiver3(0,0,0,0,0,-2,'Color',avg_color)
    quiver3(0,0,0,a_sim(1),a_sim(2),a_sim(3),'Color','m')
    quiver3(0,0,0,a_raw(1),a_raw(2),a_raw(3),'Color','c')
    set(gca,'Xcolor',avg_color);
    set(gca,'Ycolor',avg_color);
    set(gca,'Zcolor',avg_color);
    axis equal
    hold off  
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    title('Forces Dynamic Simulation')
    set(get(gca,'Title'),'Color',avg_color)
    xlim([-3.5 3.5])
    ylim([-3.5 3.5])
    zlim([-3.5 3.5])
    set(gca,'Color','k')
    subplot(1,2,2); hold on
    surf(x_mod_b_star,y_mod_b_star,z_mod_b_star,'facecolor','k','edgecolor',avg_color,'facelighting','phong');
    alpha(0)
    surf(JL_t(1)+x_mod_strk_L,JL_t(2)+y_mod_strk_L,JL_t(3)+z_mod_strk_L,'facecolor','k','edgecolor',avg_color)
    alpha(0)
    surf(JR_t(1)+x_mod_strk_R,JR_t(2)+y_mod_strk_R,JR_t(3)+z_mod_strk_R,'facecolor','k','edgecolor',avg_color)
    alpha(0)
    plot3(JL_t(1),JL_t(2),JL_t(3),'o','Color',avg_color);
    plot3(JR_t(1),JR_t(2),JR_t(3),'o','Color',avg_color);
    quiver3(0,0,0,MA(1),MA(2),MA(3),'Color','b')
    quiver3(0,0,0,MI_vel(1),MI_vel(2),MI_vel(3),'Color','r')
    quiver3(0,0,0,MI_acc(1),MI_acc(2),MI_acc(3),'Color','y')
%     quiver3(0,0,0,0,0,-2,'Color',avg_color)
    quiver3(0,0,0,w_dot_sim(1),w_dot_sim(2),w_dot_sim(3),'Color','m')
    quiver3(0,0,0,w_dot_raw(1),w_dot_raw(2),w_dot_raw(3),'Color','c')
    set(gca,'Xcolor',avg_color);
    set(gca,'Ycolor',avg_color);
    set(gca,'Zcolor',avg_color);
    axis equal
    hold off
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    title('Moments Dynamic Simulation')
    set(get(gca,'Title'),'Color',avg_color)
    xlim([-3.5 3.5])
    ylim([-3.5 3.5])
    zlim([-3.5 3.5])
    set(gca,'Color','k')
    hold off
    
 
    
end

