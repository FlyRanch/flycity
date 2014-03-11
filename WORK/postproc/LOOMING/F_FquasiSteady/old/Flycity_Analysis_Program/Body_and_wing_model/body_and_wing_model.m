function body_and_wing_model( settings, pathDB )

    savefile = 'pathDB3.mat';

    % Rotation matrices of body and wing reference frame and strokeplane
    % reference frame will be given. Body and wing kinematics will be
    % defined in strokeplane reference frame. Body and wing mass and
    % inertia will be computed.
    
    N = settings.nr_of_seq;
    
    M = settings.frame_end;

    rho_air         = settings.rho_air;
    C_mfly          = settings.C_mfly;
    rho_cuticle     = settings.rho_cuticle;
    h_wing          = settings.h_wing;
    nr_sect         = settings.nr_chord_sect;
    g               = settings.g;
       
    % Compute Rotation Matrices:
    
    rot_mat.RB = nan(3,3,M,N);
    rot_mat.RL = nan(3,3,M,N);
    rot_mat.RR = nan(3,3,M,N);
    
    beta_strk = settings.beta_strk;
    
    rot_mat.Rstr = [ cos(beta_strk) 0 -sin(beta_strk); ...
                     0              1              0; ...
                     sin(beta_strk) 0 cos(beta_strk)];
    
    for i = 1:N
        
        start = settings.start_stop(i,1);
        stop = settings.start_stop(i,2);
        
        for j = start:stop
            
            rot_mat.RB(:,:,j,i) = quat2mat(pathDB.filt.qB(j,:,i));
            rot_mat.RL(:,:,j,i) = quat2mat(pathDB.filt.qL(j,:,i));
            rot_mat.RR(:,:,j,i) = quat2mat(pathDB.filt.qR(j,:,i));
            
        end
        
    end
    
    clear start stop
    
    
    % Compute dimensions, mass, cg and inertia of body and wing:
    
    body_model.mass_fly           = nan(N,1);
    body_model.mass_body          = nan(N,1);
    body_model.length             = nan(N,1);
    body_model.Inertia            = nan(3,3,N);
    body_model.Joint_left         = nan(N,3);
    body_model.Joint_right        = nan(N,3);
    body_model.cg                 = nan(N,3);
    body_model.x_mod              = nan(21,13,N);
    body_model.y_mod              = nan(21,13,N);
    body_model.z_mod              = nan(21,13,N);
    
    wing_model.length             = nan(N,1);
    wing_model.mass               = nan(N,1);
    wing_model.virtual_mass       = nan(N,1);
    wing_model.Inertia            = nan(3,3,N);
    wing_model.virtual_Inertia    = nan(3,3,N);
    wing_model.area               = nan(N,1);
    
    wing_model.wing_cg_L          = nan(N,3);
    wing_model.y_sect_L           = nan(nr_sect,3,N);
    wing_model.chords_L           = nan(N,nr_sect);
    wing_model.x_mod_L            = nan(25,2,N);
    wing_model.y_mod_L            = nan(25,2,N);
    wing_model.z_mod_L            = nan(25,2,N);
    
    wing_model.wing_cg_R          = nan(N,3);
    wing_model.y_sect_R           = nan(nr_sect,3,N);
    wing_model.chords_R           = nan(N,nr_sect);
    wing_model.x_mod_R            = nan(25,2,N);
    wing_model.y_mod_R            = nan(25,2,N);
    wing_model.z_mod_R            = nan(25,2,N);
    
    for i = 1:N
        
        % Load body model:        
        
        xh = [0 0 0 0 0 0 1 0 0 0 1 0 0 0 1];
        
        [x_mod,y_mod,z_mod,mod_fit] = load_body_model(settings,i,xh );
        
        body_model.x_mod(:,:,i) = x_mod{1};
        body_model.y_mod(:,:,i) = y_mod{1};
        body_model.z_mod(:,:,i) = z_mod{1};
               
        wing_model.x_mod_L(:,:,i) = x_mod{2};
        wing_model.y_mod_L(:,:,i) = y_mod{2};
        wing_model.z_mod_L(:,:,i) = z_mod{2};
        
        wing_model.x_mod_R(:,:,i) = x_mod{3};
        wing_model.y_mod_R(:,:,i) = y_mod{3};
        wing_model.z_mod_R(:,:,i) = z_mod{3};
        
        body_model.length(i)           = max(body_model.x_mod(:,1,i))-min(body_model.x_mod(:,1,i));
        wing_model.length(i)           = max(wing_model.y_mod_R(:,1,i))-min(wing_model.y_mod_R(:,1,i));
        body_model.Joint_left(i,:)     = body_model.length(i).*([0.2021 -0.1055 -0.1477]);
        body_model.Joint_right(i,:)    = body_model.length(i).*([0.2021 0.1055 -0.1477]);

        body_model.mass_fly(i)        = C_mfly*(wing_model.length(i)/3)^3; % [kg]
                
        body_model_temp.x_mod = body_model.x_mod(:,:,i);
        body_model_temp.y_mod = body_model.y_mod(:,:,i);
        body_model_temp.z_mod = body_model.z_mod(:,:,i);
        wing_model_temp.x_mod = wing_model.x_mod_R(:,:,i);
        wing_model_temp.y_mod = wing_model.y_mod_R(:,:,i);
        wing_model_temp.z_mod = wing_model.z_mod_R(:,:,i);
        
        [cg_body,cg_L,cg_R,I_body,I_wing,I_v_wing,m_w,m_v_w,y_sect,chords,area_w] = ...
            comp_Inertia(rho_air,rho_cuticle,h_wing,body_model.length(i),wing_model.length(i),body_model.mass_fly(i), ...
            body_model_temp,wing_model_temp,body_model.Joint_right(i,:),nr_sect);
        
%         figure()
%         hold on
%         surf(x_mod{1},y_mod{1},z_mod{1},'facecolor',[0.5 0.5 0.5],'edgecolor','k','facelighting','phong')
%         surf(x_mod{2},y_mod{2},z_mod{2},'facecolor',[0.5 0.5 0.5],'edgecolor','k','facelighting','phong')
%         surf(x_mod{3},y_mod{3},z_mod{3},'facecolor',[0.5 0.5 0.5],'edgecolor','k','facelighting','phong')
%         axis equal
%         hold off

        body_model.cg(i,:)                   = cg_body;
        body_model.Inertia(:,:,i)            = I_body;
        
        body_model.mass_body(i)              = body_model.mass_fly(i)-2*m_w;

        wing_model.mass(i)                   = m_w;
        wing_model.virtual_mass(i)           = m_v_w;
        wing_model.Inertia(:,:,i)            = I_wing;
        wing_model.virtual_Inertia(:,:,i)    = I_v_wing;
        wing_model.area(i)                   = area_w;
        wing_model.wing_cg_L(i,:)            = cg_L;
        wing_model.wing_cg_R(i,:)            = cg_R;
        wing_model.y_sect_L(:,:,i)           = -y_sect';
        wing_model.chords_L(i,:)             = chords;
        wing_model.y_sect_R(:,:,i)           = y_sect';
        wing_model.chords_R(i,:)             = chords;
        
    
    end
    
    
    % Compute the wing kinematics defined in the strokeplane reference
    % frame:
      

    strkpln_kin.alfa        = nan(N,M);
    strkpln_kin.beta        = nan(N,M);
    strkpln_kin.roll        = nan(N,M);
    strkpln_kin.pitch       = nan(N,M);
    strkpln_kin.yaw         = nan(N,M);
    strkpln_kin.uvw         = nan(M,3,N);
    strkpln_kin.w           = nan(M,3,N);
    strkpln_kin.a_xyz       = nan(M,3,N);
    strkpln_kin.Fg          = nan(M,3,N);
    
    wing_kin.theta_L        = nan(N,M);
    wing_kin.eta_L          = nan(N,M);
    wing_kin.phi_L          = nan(N,M);
    wing_kin.U_tip_L        = nan(M,3,N);
    wing_kin.alfa_tip_L     = nan(N,M);
    wing_kin.beta_tip_L     = nan(N,M);
    
    wing_kin.theta_R        = nan(N,M);
    wing_kin.eta_R          = nan(N,M);
    wing_kin.phi_R          = nan(N,M);
    wing_kin.U_tip_R        = nan(M,3,N);
    wing_kin.alfa_tip_R     = nan(N,M);
    wing_kin.beta_tip_R     = nan(N,M);
    
    Rstr = rot_mat.Rstr;
    
        
    for i = 1:N
        
        start = settings.start_stop(i,1);
        stop = settings.start_stop(i,2);
        
        for j = start:stop

        
            RB = rot_mat.RB(:,:,j,i);
            RL = rot_mat.RL(:,:,j,i);
            RR = rot_mat.RR(:,:,j,i);
            
            U_inf  = pathDB.filt.uvw(j,:,i)';
            a_inf  = pathDB.filt.a_xyz(j,:,i)';
            w_body = pathDB.filt.wB(j,:,i)';
            
            strkpln_kin.uvw(j,:,i)  = Rstr*RB*U_inf;
            strkpln_kin.a_xyz(j,:,i)= Rstr*RB*a_inf;
            strkpln_kin.w(j,:,i)    = Rstr*w_body;
            strkpln_kin.alfa(i,j)   = real(atan2(strkpln_kin.uvw(j,3,i),strkpln_kin.uvw(j,1,i)));
            strkpln_kin.beta(i,j)   = real(atan2(strkpln_kin.uvw(j,2,i),abs(strkpln_kin.uvw(j,1,i))));

            R_180 = [1 0 0; ...
                     0 -1 0; ...
                     0 0 -1];
                 
            q_b_str = quat2mat(Rstr*RB*R_180);
            
            [roll_st,pitch_st,yaw_st]  = roll_pitch_yaw(q_b_str);
            
            strkpln_kin.roll(i,j)   = roll_st;
            strkpln_kin.pitch(i,j)  = pitch_st;
            strkpln_kin.yaw(i,j)    = yaw_st;
            
            strkpln_kin.Fg(j,:,i)   = Rstr*RB*[0; 0; -g*body_model.mass_fly(i)];
            
            r_tip_L = [0; -wing_model.length(i); 0];
            
            wL = pathDB.filt.wL(j,:,i)';

            
            wing_kin.U_tip_L(j,:,i)  = RL*RB*U_inf + RL*cross(w_body,body_model.Joint_left(i,:)'+RL'*r_tip_L) + cross(wL,r_tip_L);

            
            wing_kin.alfa_tip_L(i,j) = real(atan2(wing_kin.U_tip_L(j,3,i),wing_kin.U_tip_L(j,1,i)));
            wing_kin.beta_tip_L(i,j) = real(atan2(wing_kin.U_tip_L(j,2,i),abs(wing_kin.U_tip_L(j,1,i))));
            
            r_tip_R = [0; wing_model.length(i); 0];
            
            wR = pathDB.filt.wR(j,:,i)';
            
            wing_kin.U_tip_R(j,:,i)  = RR*RB*U_inf + RR*cross(w_body,body_model.Joint_right(i,:)'+RR'*r_tip_R) + cross(wR,r_tip_R);
            
            wing_kin.alfa_tip_R(i,j) = real(atan2(wing_kin.U_tip_R(j,3,i),wing_kin.U_tip_R(j,1,i)));
            wing_kin.beta_tip_R(i,j) = real(atan2(wing_kin.U_tip_R(j,2,i),abs(wing_kin.U_tip_R(j,1,i))));

            
            % Convert the left and right wing rotation matrices into wing
            % kinematics:
            
            [ theta_L, eta_L, phi_L, theta_R, eta_R, phi_R ] = theta_eta_phi(RL,RR,Rstr);

            wing_kin.theta_L(i,j)   = theta_L;
            wing_kin.eta_L(i,j)     = eta_L;
            wing_kin.phi_L(i,j)     = phi_L;
            
            wing_kin.theta_R(i,j)   = theta_R;
            wing_kin.eta_R(i,j)     = eta_R;
            wing_kin.phi_R(i,j)     = phi_R;
            
            clear theta_L eta_L phi_L theta_R eta_R phi_R
            
        end
        
%         figure()
%         hold on
%         plot(radtodeg(strkpln_kin.roll(i,start:stop)),'r')
%         plot(radtodeg(strkpln_kin.pitch(i,start:stop)),'g')
%         plot(radtodeg(strkpln_kin.yaw(i,start:stop)),'b')
%         hold off
%         
%         figure()
%         hold on
%         plot(strkpln_kin.uvw(start:stop,1,i),'r')
%         plot(strkpln_kin.uvw(start:stop,2,i),'g')
%         plot(strkpln_kin.uvw(start:stop,3,i),'b')
%         hold off
%         
%         figure()
%         hold on
%         plot(strkpln_kin.Fg(start:stop,1,i),'r')
%         plot(strkpln_kin.Fg(start:stop,2,i),'g')
%         plot(strkpln_kin.Fg(start:stop,3,i),'b')
%         hold off
        
%         figure()
%         hold on
%         subplot(3,1,1); plot(wing_kin.theta_L(i,start:stop))
%         subplot(3,1,2); plot(wing_kin.eta_L(i,start:stop))
%         subplot(3,1,3); plot(wing_kin.phi_L(i,start:stop))
%         hold off
%         
%         figure()
%         hold on
%         subplot(3,1,1); plot(wing_kin.theta_R(i,start:stop))
%         subplot(3,1,2); plot(wing_kin.eta_R(i,start:stop))
%         subplot(3,1,3); plot(wing_kin.phi_R(i,start:stop))
%         hold off
        
%         figure()
%         plot(wing_kin.alfa_tip_L(i,start:stop))
%         
%         figure()
%         plot(wing_kin.beta_tip_L(i,start:stop))
%         
%         figure()
%         plot(wing_kin.alfa_tip_R(i,start:stop))
%         
%         figure()
%         plot(wing_kin.beta_tip_R(i,start:stop))
%         
%         figure()
%         hold on
%         plot(wing_kin.U_tip_L(start:stop,1,i),'r')
%         plot(wing_kin.U_tip_L(start:stop,2,i),'g')
%         plot(wing_kin.U_tip_L(start:stop,3,i),'b')
%         hold off
%         
%         figure()
%         hold on
%         plot(wing_kin.U_tip_R(start:stop,1,i),'r')
%         plot(wing_kin.U_tip_R(start:stop,2,i),'g')
%         plot(wing_kin.U_tip_R(start:stop,3,i),'b')
%         hold off

%         wingtip_traject_L = nan(3,M);
%         
%         for k = start:stop
%             wingtip_traject_L(:,k) = rot_mat.Rstr*rot_mat.RL(:,:,k,i)'*[0; -1; 0];
%         end
%         
%         wingtip_traject_R = nan(3,M);
%         
%         for k = start:stop
%             wingtip_traject_R(:,k) = rot_mat.Rstr*rot_mat.RR(:,:,k,i)'*[0; 1; 0];
%         end
% 
%         figure()
%         hold on
%         plot3(wingtip_traject_L(1,start:stop),wingtip_traject_L(2,start:stop),wingtip_traject_L(3,start:stop),'r')
%         plot3(wingtip_traject_R(1,start:stop),wingtip_traject_R(2,start:stop),wingtip_traject_R(3,start:stop),'b')
%         axis equal
%         hold off
%         
%         pause
        
    end
    
    
    
    % Save the structures:
    
    save(savefile,'rot_mat','body_model','wing_model','strkpln_kin','wing_kin')
     

end

