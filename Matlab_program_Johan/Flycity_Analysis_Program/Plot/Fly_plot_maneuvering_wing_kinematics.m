function Fly_plot_maneuvering_wing_kinematics( body_pos, Rb, RL_avg, RR_avg, down_up_avg, RL_man, RR_man, down_up_man, body_model, wing_model, case_nr)

    % if case_nr = 0 (both left and right wing)
    % if case_nr = 1 (left wing only)
    % if case_nr = 2 (right wing only)


    % Plot fly wings and body in 3D:
    
    x_mod_b = body_model.x_mod;
    y_mod_b = body_model.y_mod;
    z_mod_b = body_model.z_mod;
  
    N = size(RL_avg,3);
       
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
        
    JL       = body_model.Joint_left;
    JR       = body_model.Joint_right;
    
    JL_t     = body_pos+Rb*JL;
    JR_t     = body_pos+Rb*JR;

    wing_length = wing_model.length;

    wtL  = [0; -1; 0]*wing_length;
    wtR  = [0; 1; 0]*wing_length;
    wtL1 = [0.05; -1; 0]*wing_length;
    wtR1 = [0.05; 1; 0]*wing_length;
    wtL2 = [-0.05; -1; 0]*wing_length;
    wtR2 = [-0.05; 1; 0]*wing_length;
    
    xyz_L_avg  = zeros(3,N);
    xyz_R_avg  = zeros(3,N);
    xyz1_L_avg = zeros(3,N);
    xyz2_L_avg = zeros(3,N);
    xyz1_R_avg = zeros(3,N);
    xyz2_R_avg = zeros(3,N);
    
    xyz_L_man  = zeros(3,N);
    xyz_R_man  = zeros(3,N);
    xyz1_L_man = zeros(3,N);
    xyz2_L_man = zeros(3,N);
    xyz1_R_man = zeros(3,N);
    xyz2_R_man = zeros(3,N);
    
    for i = 1:N
                
        xyz_L_avg(:,i)  = JL_t+Rb'*RL_avg(:,:,i)'*wtL;
        xyz_R_avg(:,i)  = JR_t+Rb'*RR_avg(:,:,i)'*wtR;
        xyz1_L_avg(:,i) = JL_t+Rb'*RL_avg(:,:,i)'*wtL1;
        xyz2_L_avg(:,i) = JL_t+Rb'*RL_avg(:,:,i)'*wtL2;
        xyz1_R_avg(:,i) = JR_t+Rb'*RR_avg(:,:,i)'*wtR1;
        xyz2_R_avg(:,i) = JR_t+Rb'*RR_avg(:,:,i)'*wtR2;
        
        xyz_L_man(:,i)  = JL_t+Rb'*RL_man(:,:,i)'*wtL;
        xyz_R_man(:,i)  = JR_t+Rb'*RR_man(:,:,i)'*wtR;
        xyz1_L_man(:,i) = JL_t+Rb'*RL_man(:,:,i)'*wtL1;
        xyz2_L_man(:,i) = JL_t+Rb'*RL_man(:,:,i)'*wtL2;
        xyz1_R_man(:,i) = JR_t+Rb'*RR_man(:,:,i)'*wtR1;
        xyz2_R_man(:,i) = JR_t+Rb'*RR_man(:,:,i)'*wtR2;
        
    end
    
    
    % find the downstroke and upstroke point
    
    xyz_down_L_avg  = xyz_L_avg(:,1);
    xyz_down_R_avg  = xyz_R_avg(:,1);
    xyz_up_L_avg    = xyz_L_avg(:,round(N*down_up_avg));
    xyz_up_R_avg    = xyz_R_avg(:,round(N*down_up_avg));
    
    xyz_down_L_man  = xyz_L_man(:,1);
    xyz_down_R_man  = xyz_R_man(:,1);
    xyz_up_L_man    = xyz_L_man(:,round(N*down_up_man));
    xyz_up_R_man    = xyz_R_man(:,round(N*down_up_man));

    hold on
    
    avg_color = [0.7 0.7 0.7];
    
    % Plot the body and the average wingbeat:
    
    if case_nr == 0
    
        surf(body_pos(1)+x_mod_b_star,body_pos(2)+y_mod_b_star,body_pos(3)+z_mod_b_star,'facecolor','k','edgecolor',avg_color,'facelighting','phong');
        plot3(xyz_L_avg(1,:),xyz_L_avg(2,:),xyz_L_avg(3,:),'Color',avg_color)
        plot3(xyz_R_avg(1,:),xyz_R_avg(2,:),xyz_R_avg(3,:),'Color',avg_color)
        for i = 1:N
            plot3([xyz1_L_avg(1,i) xyz2_L_avg(1,i)],[xyz1_L_avg(2,i) xyz2_L_avg(2,i)],[xyz1_L_avg(3,i) xyz2_L_avg(3,i)],'Color',avg_color)
            plot3([xyz1_R_avg(1,i) xyz2_R_avg(1,i)],[xyz1_R_avg(2,i) xyz2_R_avg(2,i)],[xyz1_R_avg(3,i) xyz2_R_avg(3,i)],'Color',avg_color)
        end
        plot3(JL_t(1),JL_t(2),JL_t(3),'o','Color',avg_color);
        plot3(JR_t(1),JR_t(2),JR_t(3),'o','Color',avg_color);
        plot3([JL_t(1) xyz_down_L_avg(1)],[JL_t(2) xyz_down_L_avg(2)],[JL_t(3) xyz_down_L_avg(3)],'Color',avg_color);
        plot3([JL_t(1) xyz_up_L_avg(1)],[JL_t(2) xyz_up_L_avg(2)],[JL_t(3) xyz_up_L_avg(3)],'Color',avg_color);
        plot3([JR_t(1) xyz_down_R_avg(1)],[JR_t(2) xyz_down_R_avg(2)],[JR_t(3) xyz_down_R_avg(3)],'Color',avg_color);
        plot3([JR_t(1) xyz_up_R_avg(1)],[JR_t(2) xyz_up_R_avg(2)],[JR_t(3) xyz_up_R_avg(3)],'Color',avg_color);
    
    elseif case_nr == 1
        
        surf(body_pos(1)+x_mod_b_star,body_pos(2)+y_mod_b_star,body_pos(3)+z_mod_b_star,'facecolor','k','edgecolor',avg_color,'facelighting','phong');
        plot3(xyz_L_avg(1,:),xyz_L_avg(2,:),xyz_L_avg(3,:),'Color',avg_color)
        for i = 1:N
            plot3([xyz1_L_avg(1,i) xyz2_L_avg(1,i)],[xyz1_L_avg(2,i) xyz2_L_avg(2,i)],[xyz1_L_avg(3,i) xyz2_L_avg(3,i)],'Color',avg_color)
        end
        plot3(JL_t(1),JL_t(2),JL_t(3),'o','Color',avg_color);
        plot3([JL_t(1) xyz_down_L_avg(1)],[JL_t(2) xyz_down_L_avg(2)],[JL_t(3) xyz_down_L_avg(3)],'Color',avg_color);
        plot3([JL_t(1) xyz_up_L_avg(1)],[JL_t(2) xyz_up_L_avg(2)],[JL_t(3) xyz_up_L_avg(3)],'Color',avg_color);

    elseif case_nr == 2

        surf(body_pos(1)+x_mod_b_star,body_pos(2)+y_mod_b_star,body_pos(3)+z_mod_b_star,'facecolor','k','edgecolor',avg_color,'facelighting','phong');
        plot3(xyz_R_avg(1,:),xyz_R_avg(2,:),xyz_R_avg(3,:),'Color',avg_color)
        for i = 1:N
            plot3([xyz1_R_avg(1,i) xyz2_R_avg(1,i)],[xyz1_R_avg(2,i) xyz2_R_avg(2,i)],[xyz1_R_avg(3,i) xyz2_R_avg(3,i)],'Color',avg_color)
        end
        plot3(JR_t(1),JR_t(2),JR_t(3),'o','Color',avg_color);
        plot3([JR_t(1) xyz_down_R_avg(1)],[JR_t(2) xyz_down_R_avg(2)],[JR_t(3) xyz_down_R_avg(3)],'Color',avg_color);
        plot3([JR_t(1) xyz_up_R_avg(1)],[JR_t(2) xyz_up_R_avg(2)],[JR_t(3) xyz_up_R_avg(3)],'Color',avg_color);
    
    end
    
    % Plot the maneuvering wingbeat:
    
    if case_nr == 0
    
        plot3(xyz_L_man(1,:),xyz_L_man(2,:),xyz_L_man(3,:),'r')
        plot3(xyz_R_man(1,:),xyz_R_man(2,:),xyz_R_man(3,:),'g')
        for i = 1:N
            plot3([xyz1_L_man(1,i) xyz2_L_man(1,i)],[xyz1_L_man(2,i) xyz2_L_man(2,i)],[xyz1_L_man(3,i) xyz2_L_man(3,i)],'r')
            plot3([xyz1_R_man(1,i) xyz2_R_man(1,i)],[xyz1_R_man(2,i) xyz2_R_man(2,i)],[xyz1_R_man(3,i) xyz2_R_man(3,i)],'g')
        end
        plot3([JL_t(1) xyz_down_L_man(1)],[JL_t(2) xyz_down_L_man(2)],[JL_t(3) xyz_down_L_man(3)],'Color','r');
        plot3([JL_t(1) xyz_up_L_man(1)],[JL_t(2) xyz_up_L_man(2)],[JL_t(3) xyz_up_L_man(3)],'Color','r');
        plot3([JR_t(1) xyz_down_R_man(1)],[JR_t(2) xyz_down_R_man(2)],[JR_t(3) xyz_down_R_man(3)],'Color','g');
        plot3([JR_t(1) xyz_up_R_man(1)],[JR_t(2) xyz_up_R_man(2)],[JR_t(3) xyz_up_R_man(3)],'Color','g');
    
    elseif case_nr == 1
        
        plot3(xyz_L_man(1,:),xyz_L_man(2,:),xyz_L_man(3,:),'r')
        for i = 1:N
            plot3([xyz1_L_man(1,i) xyz2_L_man(1,i)],[xyz1_L_man(2,i) xyz2_L_man(2,i)],[xyz1_L_man(3,i) xyz2_L_man(3,i)],'r')
        end
        plot3([JL_t(1) xyz_down_L_man(1)],[JL_t(2) xyz_down_L_man(2)],[JL_t(3) xyz_down_L_man(3)],'Color','r');
        plot3([JL_t(1) xyz_up_L_man(1)],[JL_t(2) xyz_up_L_man(2)],[JL_t(3) xyz_up_L_man(3)],'Color','r');
        
    elseif case_nr == 2
        
        plot3(xyz_R_man(1,:),xyz_R_man(2,:),xyz_R_man(3,:),'g')
        for i = 1:N
            plot3([xyz1_R_man(1,i) xyz2_R_man(1,i)],[xyz1_R_man(2,i) xyz2_R_man(2,i)],[xyz1_R_man(3,i) xyz2_R_man(3,i)],'g')
        end
        plot3([JR_t(1) xyz_down_R_man(1)],[JR_t(2) xyz_down_R_man(2)],[JR_t(3) xyz_down_R_man(3)],'Color','g');
        plot3([JR_t(1) xyz_up_R_man(1)],[JR_t(2) xyz_up_R_man(2)],[JR_t(3) xyz_up_R_man(3)],'Color','g');
    
    end
    
    set(gca,'Xcolor',avg_color);
    set(gca,'Ycolor',avg_color);
    set(gca,'Zcolor',avg_color);
    axis equal
    hold off
    
    set(gca,'Color','k')
    
%     set(fignum,'Color','k')
    
end

