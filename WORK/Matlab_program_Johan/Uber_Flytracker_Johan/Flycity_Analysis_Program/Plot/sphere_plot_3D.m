function sphere_plot_3D( settings, pathDB, Rb, RL_raw, RR_raw, RL_filt, RR_filt, seq_nr)

    % Plot fly wings and body in 3D:
    
    x_mod_b = pathDB.body_model.x_mod(:,:,seq_nr);
    y_mod_b = pathDB.body_model.y_mod(:,:,seq_nr);
    z_mod_b = pathDB.body_model.z_mod(:,:,seq_nr);
    
       
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
        
    JL       = pathDB.body_model.Joint_left(seq_nr,:);
    
    Lw       = pathDB.wing_model.length(seq_nr);
    
    J_LR     = Rb*[1 0 0; 0 -1 0; 0 0 -1]*[JL(1); 0; JL(3)];
    
    N = size(RL_raw,3);
    
    wtL_raw = zeros(3,N);
    wtR_raw = zeros(3,N);
    wtL_filt = zeros(3,N);
    wtR_filt = zeros(3,N);
    
    wtL = [0; -Lw; 0];
    wtR = [0; Lw; 0];
    
    for i = 1:N
        
        wtL_raw(:,i) = J_LR+[1 0 0; 0 -1 0; 0 0 -1]*Rb'*RL_raw(:,:,i)'*wtL;
        wtR_raw(:,i) = J_LR+[1 0 0; 0 -1 0; 0 0 -1]*Rb'*RR_raw(:,:,i)'*wtR;
        wtL_filt(:,i) = J_LR+[1 0 0; 0 -1 0; 0 0 -1]*Rb'*RL_filt(:,:,i)'*wtL;
        wtR_filt(:,i) = J_LR+[1 0 0; 0 -1 0; 0 0 -1]*Rb'*RR_filt(:,:,i)'*wtR;
        
    end
    
    
    k = 5;
    n = 2^k-1;
    [x,y,z] = sphere(n);

    hold on
    surf(J_LR(1)+x.*Lw,J_LR(2)+y.*Lw,J_LR(3)+z.*Lw,'FaceColor','black','EdgeColor','none');
    alpha(0.2)
    surf(x_mod_b_star,y_mod_b_star,z_mod_b_star,'facecolor',[0.3 0.3 0.3],'edgecolor','k','facelighting','phong');
    plot3(wtL_raw(1,:),wtL_raw(2,:),wtL_raw(3,:),'Color','g')
    plot3(wtR_raw(1,:),wtR_raw(2,:),wtR_raw(3,:),'Color','g')
    plot3(wtL_filt(1,:),wtL_filt(2,:),wtL_filt(3,:),'Color','r')
    plot3(wtR_filt(1,:),wtR_filt(2,:),wtR_filt(3,:),'Color','r')
    xlim([-1.05*Lw+J_LR(1) 1.05*Lw+J_LR(1)])
    ylim([-1.05*Lw+J_LR(2) 1.05*Lw+J_LR(2)])
    zlim([-1.05*Lw+J_LR(3) 1.05*Lw+J_LR(3)])
    hold off
    
end