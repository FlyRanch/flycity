function Fly_plot_3D_wingtip( body_pos, Rb, RL, RR, down_up, body_model, wing_model, fignum )


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
    




    figure(fignum)
    hold on
    surf(body_pos(1)+x_mod_b_star,body_pos(2)+y_mod_b_star,body_pos(3)+z_mod_b_star,'facecolor','k','edgecolor',[0.3 0.3 0.3],'facelighting','phong');
    plot3(xyz_L(1,:),xyz_L(2,:),xyz_L(3,:),'r')
    plot3(xyz_R(1,:),xyz_R(2,:),xyz_R(3,:),'g')
    for i = 1:N
        plot3([xyz1_L(1,i) xyz2_L(1,i)],[xyz1_L(2,i) xyz2_L(2,i)],[xyz1_L(3,i) xyz2_L(3,i)],'r')
        plot3([xyz1_R(1,i) xyz2_R(1,i)],[xyz1_R(2,i) xyz2_R(2,i)],[xyz1_R(3,i) xyz2_R(3,i)],'g')
    end
    plot3(JL_t(1),JL_t(2),JL_t(3),'o','Color',[0.3 0.3 0.3]);
    plot3(JR_t(1),JR_t(2),JR_t(3),'o','Color',[0.3 0.3 0.3]);
    plot3([JL_t(1) xyz_down_L(1)],[JL_t(2) xyz_down_L(2)],[JL_t(3) xyz_down_L(3)],'Color',[0.3 0.3 0.3]);
    plot3([JL_t(1) xyz_up_L(1)],[JL_t(2) xyz_up_L(2)],[JL_t(3) xyz_up_L(3)],'Color',[0.3 0.3 0.3]);
    plot3([JR_t(1) xyz_down_R(1)],[JR_t(2) xyz_down_R(2)],[JR_t(3) xyz_down_R(3)],'Color',[0.3 0.3 0.3]);
    plot3([JR_t(1) xyz_up_R(1)],[JR_t(2) xyz_up_R(2)],[JR_t(3) xyz_up_R(3)],'Color',[0.3 0.3 0.3]);
    set(gca,'Xcolor',[0.3 0.3 0.3]);
    set(gca,'Ycolor',[0.3 0.3 0.3]);
    set(gca,'Zcolor',[0.3 0.3 0.3]);
    axis equal
    hold off
    
    set(gca,'Color','k')
    
    set(fignum,'Color','k')
    
end

