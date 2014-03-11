function Aero_force_plot( body_pos, Rb, RL, RR, down_up, FA_L, FA_R, body_model, wing_model, case_nr)

    % if case_nr = 0 (both left and right wing)
    % if case_nr = 1 (left wing only)
    % if case_nr = 2 (right wing only)


    % Plot fly wings and body in 3D:
    
    x_mod_b = body_model.x_mod;
    y_mod_b = body_model.y_mod;
    z_mod_b = body_model.z_mod;
  
    N = size(RL,3);
       
    n_b = size(x_mod_b,1);
    m_b = size(x_mod_b,2);
    
    x_mod_b_star = zeros(n_b,m_b);
    y_mod_b_star = zeros(n_b,m_b);
    z_mod_b_star = zeros(n_b,m_b);
    
   
    % Transform the body to current orientation

    for j = 1:n_b
        for k = 1:m_b

            mod_t = Rb'*[1 0 0; 0 -1 0; 0 0 -1]*[ x_mod_b(j,k); y_mod_b(j,k); z_mod_b(j,k)];
                
            x_mod_b_star(j,k) = mod_t(1);
            y_mod_b_star(j,k) = mod_t(2);
            z_mod_b_star(j,k) = mod_t(3);
                
            clear mod_t

        end
    end 

    % Transform the left wing to the current orientation
        
    JL       = body_model.Joint_left;
    JR       = body_model.Joint_right;
    
    JL_t     = body_pos+Rb'*JL;
    JR_t     = body_pos+Rb'*JR;

    wing_length = wing_model.length;

    wtL  = [0; -1; 0]*wing_length;
    wtR  = [0; 1; 0]*wing_length;
    wtL1 = [0.05; -1; 0]*wing_length;
    wtR1 = [0.05; 1; 0]*wing_length;
    wtL2 = [-0.05; -1; 0]*wing_length;
    wtR2 = [-0.05; 1; 0]*wing_length;
    
    xyz_L  = zeros(3,N);
    xyz_R  = zeros(3,N);
    xyz1_L = zeros(3,N);
    xyz2_L = zeros(3,N);
    xyz1_R = zeros(3,N);
    xyz2_R = zeros(3,N);
    
    for i = 1:N
                
        xyz_L(:,i)  = JL_t+Rb'*RL(:,:,i)'*wtL;
        xyz_R(:,i)  = JR_t+Rb'*RR(:,:,i)'*wtR;
        xyz1_L(:,i) = JL_t+Rb'*RL(:,:,i)'*wtL1;
        xyz2_L(:,i) = JL_t+Rb'*RL(:,:,i)'*wtL2;
        xyz1_R(:,i) = JR_t+Rb'*RR(:,:,i)'*wtR1;
        xyz2_R(:,i) = JR_t+Rb'*RR(:,:,i)'*wtR2;
        
    end
    
    
    % find the downstroke and upstroke point
    
    xyz_down_L  = xyz_L(:,1);
    xyz_down_R  = xyz_R(:,1);
    xyz_up_L    = xyz_L(:,round(N*down_up));
    xyz_up_R    = xyz_R(:,round(N*down_up));
    
    % Transfer aerodynamic forces to inertial reference frame:
    
    FA_L_xyz = zeros(3,N);
    FA_R_xyz = zeros(3,N);
    
    quiv_scale = 2e4;
    
    for i = 1:N
        
        FA_L_xyz(:,i) = Rb'*RL(:,:,i)'*FA_L(:,i)*quiv_scale;
        FA_R_xyz(:,i) = Rb'*RR(:,:,i)'*FA_R(:,i)*quiv_scale;
        
    end
    

    hold on
    
    avg_color = [0.5 0.5 0.5];
    
    % Plot the body and the average wingbeat:
    
    if case_nr == 0
    
        surf(body_pos(1)+x_mod_b_star,body_pos(2)+y_mod_b_star,body_pos(3)+z_mod_b_star,'facecolor',avg_color,'edgecolor','k','facelighting','phong');
        plot3(xyz_L(1,:),xyz_L(2,:),xyz_L(3,:),'Color','k')
        plot3(xyz_R(1,:),xyz_R(2,:),xyz_R(3,:),'Color','k')
        for i = 1:N
            plot3([xyz1_L(1,i) xyz2_L(1,i)],[xyz1_L(2,i) xyz2_L(2,i)],[xyz1_L(3,i) xyz2_L(3,i)],'Color','r')
            plot3([xyz1_R(1,i) xyz2_R(1,i)],[xyz1_R(2,i) xyz2_R(2,i)],[xyz1_R(3,i) xyz2_R(3,i)],'Color','g')
            quiver3(xyz_L(1,i),xyz_L(2,i),xyz_L(3,i),FA_L_xyz(1,i),FA_L_xyz(2,i),FA_L_xyz(3,i),'Color','r')
            quiver3(xyz_R(1,i),xyz_R(2,i),xyz_R(3,i),FA_R_xyz(1,i),FA_R_xyz(2,i),FA_R_xyz(3,i),'Color','g')
        end
        plot3([JL_t(1) xyz_down_L(1)],[JL_t(2) xyz_down_L(2)],[JL_t(3) xyz_down_L(3)],'Color','k');
        plot3([JL_t(1) xyz_up_L(1)],[JL_t(2) xyz_up_L(2)],[JL_t(3) xyz_up_L(3)],'Color','k');
        plot3([JR_t(1) xyz_down_R(1)],[JR_t(2) xyz_down_R(2)],[JR_t(3) xyz_down_R(3)],'Color','k');
        plot3([JR_t(1) xyz_up_R(1)],[JR_t(2) xyz_up_R(2)],[JR_t(3) xyz_up_R(3)],'Color','k');
        
    
    elseif case_nr == 1
        
        surf(body_pos(1)+x_mod_b_star,body_pos(2)+y_mod_b_star,body_pos(3)+z_mod_b_star,'facecolor',avg_color,'edgecolor','k','facelighting','phong');
        plot3(xyz_L(1,:),xyz_L(2,:),xyz_L(3,:),'Color','k')
        for i = 1:N
            plot3([xyz1_L(1,i) xyz2_L(1,i)],[xyz1_L(2,i) xyz2_L(2,i)],[xyz1_L(3,i) xyz2_L(3,i)],'Color','r')
            quiver3(xyz_L(1,i),xyz_L(2,i),xyz_L(3,i),FA_L_xyz(1,i),FA_L_xyz(2,i),FA_L_xyz(3,i),'Color','r')
        end
        plot3([JL_t(1) xyz_down_L(1)],[JL_t(2) xyz_down_L(2)],[JL_t(3) xyz_down_L(3)],'Color','k');
        plot3([JL_t(1) xyz_up_L(1)],[JL_t(2) xyz_up_L(2)],[JL_t(3) xyz_up_L(3)],'Color','k');

    elseif case_nr == 2

        surf(body_pos(1)+x_mod_b_star,body_pos(2)+y_mod_b_star,body_pos(3)+z_mod_b_star,'facecolor',avg_color,'edgecolor','k','facelighting','phong');
        plot3(xyz_R(1,:),xyz_R(2,:),xyz_R(3,:),'Color','k')
        for i = 1:N
            plot3([xyz1_R(1,i) xyz2_R(1,i)],[xyz1_R(2,i) xyz2_R(2,i)],[xyz1_R(3,i) xyz2_R(3,i)],'Color','g')
            quiver3(xyz_R(1,i),xyz_R(2,i),xyz_R(3,i),FA_R_xyz(1,i),FA_R_xyz(2,i),FA_R_xyz(3,i),'Color','g')
        end
        plot3([JR_t(1) xyz_down_R(1)],[JR_t(2) xyz_down_R(2)],[JR_t(3) xyz_down_R(3)],'Color','k');
        plot3([JR_t(1) xyz_up_R(1)],[JR_t(2) xyz_up_R(2)],[JR_t(3) xyz_up_R(3)],'Color','k');
    
    end
    axis equal
    hold off
    
    
end